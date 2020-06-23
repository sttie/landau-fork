# Description
Landau is a **lan**guage for **d**ynamical systems with **au**tomatic differentiation.

Landau is a Turing incomplete statically typed domain-specific differentiable
language. The Turing incompleteness provides the ability of sophisticated
source code analysis and, as a result, a highly optimized compiled code. Among
other things, the language syntax supports functions, compile-time ranged for
loops, if/else branching constructions, real variables and arrays, and the
ability to manually discard calculation where the automatic derivatives values
are expected to be negligibly small. In spite of reasonable restrictions, the
language is rich enough to express and differentiate any cumbersome
paper-equation with practically no effort.

## Automatic differentiation
The automatic differentiation (forward mode) technique can be built on the
concept of dual numbers. Dual numbers can be thought as a pair $`\langle x, x'
\rangle`$ of a numerical value and its numerical derivative. Arithmetic applies
to the right part according to differentiation rules:

```math
\begin{aligned}
\langle u, u' \rangle \pm \langle v, v' \rangle =& \langle u \pm v, u' \pm v' \rangle\\
\langle u, u' \rangle \langle v, v' \rangle =& \langle uv, u'v + uv' \rangle\\
\langle u, u' \rangle / \langle v, v' \rangle =& \left\langle \frac{u}{v},
  \frac{u'v - uv'}{v^2} \right\rangle, v \neq 0\\ \end{aligned}
```

However, the manual usage of dual numbers makes the implementation of the equations more
tedious and harder to maintain.

## Problem formulation
The task is to calculate the (automatic) derivatives of functions
$`f_i(x_1,\ldots,x_j)`$ w.r.t. parameters $`\{p_k\}`$. Numerical values of
$`\frac{\mathrm{d} x_j}{\mathrm{d} p_k}`$ are given externally.

### Example
```math
    \begin{aligned}
      \frac{\mathrm{d} x_i}{\mathrm{d} x_0^{(j)}} = &\ J_{ij} \quad \text{(numerical jacobian)}\\
      f_i = &\ 2\sin(x_i) + x_i\\
      \frac{\mathrm{d} f_i}{\mathrm{d} x_0^{(j)}} = &\ ?
    \end{aligned}
```

Describing this problem in Landau:

```python
    parameter[k] x0 # Declare parameters
                    # w.r.t which derivatives should be taken
    real[k * k] func(real[k] x, real[k * k] dxdx0) { 
      x[:] ' x0[:] = dxdx0[:]          # Tell the compiler 
                                       # about dx/dx0 jacobian

      real[k] f = 2 * sin(x[:]) + x[:] # Compose the function's body

      func[:] = f[:] ' x0[:]           # Write the numerical value
                                       # of df/dx0 jacobian 
                                       # to the return array
    }
```

## Differentiation method
Landau uses a hybrid approach of derivatives calculation: symbolic
differentiation at the expression level and automatic (forward mode) at the
assignation level. Consider the following piece of Landau code:
```python
      real[k] g, f
      g[:] = 2 * sin(x[:]) + x[:]
      for i = [0 : k]
        f[i] = 2 * g[i] + g[i] * x[i]
```
The right-hand side of the assignation is differentiated symbolically. The
result of computation is stored in the corresponding array (`dgdx0`). Whenever
the result of $`\mathrm{d}g/\mathrm{d}x_0`$ is needed, it is taken from the
corresponding array.
```python
      dgdx0[..] = 2 * cos(x[..]) * dxdx0[..] + dxdx0[..]
      g[:] = 2 * sin(x[:]) + x[:]
      for i = [0 : k]
        dfdx0[..] = 2 * dgdx0[..] 
                      + dgdx0[..] * x[i] + g[i] * dxdx0[..]
        f[i] = 2 * g[i] + g[i] * x[i]
```
The `..` is not a valid syntax construction, but rather a placeholder for a
slice index expression not mentioned in the listing for the simplicity purpose.
Indexing mechanics is revealed later in the Implementation details section.

## Implementation
The Landau transpiler is built on top of the [Racket](https://racket-lang.org)
platform and consists of two main parts: the syntax analyzer and the code
generator. The job of the syntax analyzer is:

- unrolling all loops into an assignation sequence
- tracing dependencies of variables and their derivatives
- checking the program for correctness, e.g. catching out-of-range errors.

For example, this piece of Landau code will cause the compile-time error:
    
```python
            real[4] x
            for i = [0 : 5]
              x[i] = 1. # Compile error
            ~~~~^ index 4 out of range [0, 3]
```

The code generator is capable to emit Racket and ANSI-C code. It is possible
to use 80-bits floating-point numbers (extended precision) without changing
the Landau source.

## Landau advantages
Landau facilitates the development of dynamical models used in celestial
mechanics and possibly other areas. Landau:

- generates highly efficient C code (or not-so-efficient Racket syntax objects
  that are useful if you are programming in Racket);
- exempts user of C pointer managing trouble and low-level debugging;
- allows easy notation of mathematical equations, while also allowing loops and
  functions;
- and, of course, provides first-class automatic derivatives (while you can use
  Landau without them, too).

The Landau-generated C code is portable on different compilers on Windows,
Linux or MacOS.

## Development in Landau
The model development consists of three stages:

1. Manual model description in Landau language (`model.dau`)
2. Transpiling to ANSI-C (code generation `model.c`)
3. Compiling to a shared library for a specific platform
  (`model.dll`, `model.so`, `model.dylib`)

## Example: Spacecraft movement
As an example, consider the Landau program for modeling spacecraft movement
around a planet -- [spacecraft.dau](tests/spacecraft.dau). Spacecraft's
initial position and velocity, as well as the gravitational parameter of the
planet, are supposed to be determined by an external implementation of the
least-squares method. Derivatives of the right part of the movement equation
w.r.t. the initial state vector and the planet mass are handled by Landau.

## Publications
The work is published in papers:
- I. Dolgakov, D. Pavlov. "Landau: language for dynamical systems with
  automatic differentiation". Zap. Nauch. Sem. POMI 485 (2019)
  [(Full text)](ftp://ftp.pdmi.ras.ru/pub/publicat/znsl/v485/p078.pdf).
- I. Dolgakov, D. Pavlov. "Using capabilities of Racket to implement a
  domain-specific language". Computer Assisted Mathematics 1 (2019)
  [(Full text)](http://cte.eltech.ru/ojs/index.php/cam/article/view/1625).

# Installation
1. Download and install [Racket](https://download.racket-lang.org/).
2. Clone repository: `$ git clone https://gitlab.iaaras.ru/iaaras/landau.git`.
3. Link Landau as a package: `$ cd Landau/ && raco link landau`.

# Configuration
`config.json`:
```
{   
    // Language to compile .dau file to
    "target_language": "racket" | "ansi-c", 

    // Whether to use 80-bits floating point
    "use_extfloat": true | false,

    // Where to write output .c file, 
    // in case the `target_language` is "ansi-c"
    "output_directory": DIRECTORY_PATH
}
```

# Examples
`$ cd tests`
```
$ racket spacecraft.dau -c -extfl
Success.
The output is written to DIRECTORY_PATH/spacecraft.c
```

# Comparison with other AD tools 
There exist a large number of forward-mode AD software tools for
differentiating functions that are written in general-purpose programming
languages, like Fortran (ADIFOR) [1], C (ADIC) [2] or C++ (ADOL-C) [3]. Rich
features of the ``host'' languages, like arrays, loops, conditions, and
recursion, often make it difficult to implement a practically usable AD system
without imposing limitations on the language and/or extra technical work when
specifying the function, especially in presence of multi-dimensional functions
with many independent variables.

On the other hand, there exist a number of languages developed specially for AD
tasks, like Taylor [4] and VLAD [5, 6]. Taylor syntax, while very simple and
natural, is very limited (no conditionals, loops, arrays, or subprocedures).
VLAD, a functional Scheme-like language, has conditionals, loops, recursion,
and subprocedures, but does not have arrays or mutability.

Finally, there are tools for differentiating functions specified as
mathematical expressions in mathematical computing systems, like MATLAB (ADMAT)
[7] or Mathematica (TIDES) [8, 9]. Such tools often require a bigger effort
(as compared to a general-purpose languages) to input a practical dynamical
system of large dimension with a lot of free variables.

Like VLAD, Landau is a domain-specific language designed with automatic
differentiation in mind. Like TIDES and Taylor, Landau offers C code
generation. Like general-purpose languages, Landau has common control flow
constructs, arrays, and mutability; but unlike general-purpose languages,
Landau embraces Turing incompleteness to perform static source analysis and
generate efficient code.

# Implementation details
There are AD approaches where an initial program is unrolled to a sequence
of simple assignments (Wengert list) and AD problem is solved in linear algebra
formalism. The initial problem is reduced to a sparse linear equations
system. Even though there are known lots of sparse matrix algorithms it is not
obvious if it is a good idea to unroll all loops to a giant sparse matrix.

We decided to take the code generation approach preserving the initial
program's loops, e.g. Landau loop is transformed into the C or Racket loop
where the loop-body is augmented with a derivative calculation code if needed.

Unlike the other AD tools, Landau uses derivative annotations inside a 
program to be able to predict a number of nonzero derivatives and reduce a
machine resources used for jacobian computation and storage. Other AD tools
often do not use annotations and hence, make it impossible to perform radical
compile-time optimizations of array's derivatives computations.

To illustrate our jacobian storage approach let us take a look at the listing
from the Differentiation method section. One of the simplest possible storage
scheme is to keep all jacobian components in array row by row, e.g. to put
`df_i/dx0_j` in `dfdx0[i * length(x0) + j]`.

The `..` dots, hence, should be substituted as follows: 
```
dfdx0[i:] = 2 * dgdx0[i:]
              + dgdx0[i:] * x[i] + g[i] * dxdx0[i:]
```
This approach performs well in case of dense jacobians but causes memory and
computation time waste in sparse cases.

There is a widespread jacobian compression method based on grouping so-called
orthogonal jacobian columns. Two columns are called orthogonal if they don't
have nonzero value in the same row. Compression is made by multiplication of the
initial sparse jacobian by the specific permutation matrix (so-called seed matrix).

Example:
Sparse jacobian where columns 1 and 3 are orthogonal 
```math
  J = 
  \begin{bmatrix}
    a_{1,1} & a_{1,2} &  \\
    & & a_{2,3}  \\
    & & a_{3,3}  \\
  \end{bmatrix}.
```

Seed matrix:
```math
  S = 
  \begin{bmatrix}
    1 & 0  \\
    0 & 1  \\
    1 & 0  \\
  \end{bmatrix}.
```

Multiplying J and S we get the compressed jacobian B:
```math
  B = J \cdot S = 
  \begin{bmatrix}
    a_{1,1} & a_{1,2} \\
    a_{2,3} &  \\
    a_{3,3} &  \\
  \end{bmatrix}.
```
The major shortcoming of this approach is a poor compression in the absence of
orthogonal columns which is the case if, say, jacobian has a single row and
a column full of nonzero elements.

We chose to use another compression technique storing only nonzero jacobian
values row by row. The sparsity is handled by packing useful jacobian elements
into smaller arrays and generating mappings from the packed derivative indexes to
the original ones and inverse mappings, which map the original indexes to the
packed ones. Here is a brief demonstration of discussed approach
applied to the `dfdx0[..]` assignation from the listing:

```
for dfdx0_compressed_idx in length(dfdx0_mappings[i]):
 dfdx0_true_idx = dfdx0_mappings[i][dfdx0_compressed_idx]
 dfdx0[i * length(dfdx0_mappings[i]) + dfdx0_compressed_idx] = 
   2 * dgdx0[i * length(dgdx0_mappings[i]) + dgdx0_inverse_mapping[i][dfdx0_true_idx]]
     + dgdx0[i * length(dgdx0_mappings[i]) + dgdx0_inverse_mapping[i][dfdx0_true_idx]] * x[i]
     + g[i] * dxdx0[i * length(dxdx0_mappings[i]) + dxdx0_inverse_mapping[i][dfdx0_true_idx]]
```

Where `dfdx0`, `dgdx0`, `dxdx0` are compressed jacobians, `dfdx0_mappings` is
the mappings array for `dfdx0` and `dgdx0_inverse_mapping`,
`dxdx0_inverse_mapping` are inverse mappings for `dgdx0` and `dxdx0`.

While being close to optimal in terms of computation time, the chosen AD
algorithm is not very memory efficient. The main problem is that inverse
mapping of a single array cell takes O(maxidx) in space, where maxidx is a
maximal index of nonzero derivative for the array element. To illustrate the
problem, say, the mapping is equal to:
```
value [0 5 6 7 99 100]
index  0 1 2 3 4  5
```
The inverse mapping is then
```
value [0 ... 1 2 3 ... 4  5]
index  0 ... 5 6 7 ... 99 100 
```
To store this 5 elements inverse mapping should be an array of 100 elements,
which is an absolute waste of memory. Since a set of inverse mapping values is
static (read-only) and known at the compile-time the ideal solution of inverse
mapping storage problem would be an open addressed hash table with a minimal 
perfect hash function used for a key construction. But as we know there is no
simple way to construct one in case of open addressed hash table. So a better
jacobian compression algorithm it is still under consideration. 

# Planned improvements
## Syntax
- FFI (calling C functions inside Landau program).
- Complex numbers.
- Syntactic sugar for array operations.
- Multidimensional arrays.

## Core
Better Jacobian compression algorithm.

## Static analysis
Unrolling all loops and evaluating all compile-time computable variables is a
simple but inefficient approach. We are looking for better ways of performing
the static syntax analysis.

## Code generation
- Generate AVX intrinsic functions for the C backend.
- Generate unsafe vector operations for the Racket backend.

# References
[1]	Ch. Bischof, A. Carle, G. Corliss, A. Griewank, P. Hovland,
“ADIFOR–generating derivative codes from Fortran programs”, Scientific
Programming, 1:1 (1992), 11–29

[2]	Ch. Bischof, L. Roh, A. J. Mauer-Oats, “ADIC: an extensible automatic
differentiation tool for ANSI-C”, Software: Practice and Experience, 27:12
(1997), 1427–1456  

[3]	A. Griewank, D. Juedes, J. Utke, “Algorithm 755: ADOL-C: a package for
the automatic differentiation of algorithms written in C/C++”, ACM
Transactions on Mathematical Software (TOMS), 22:2 (1996), 131–167  

[4]	À. Jorba, M. Zou, “A Software Package for the Numerical Integration
of ODEs by Means of High-Order Taylor Methods”, Experimental Mathematics,
14:1 (2005), 99–117  

[5]	J. M. Siskind, B. A. Pearlmutter, “Nesting forward-mode AD in a
functional framework”, Higher-Order and Symbolic Computation, 21:4 (2008),
361–376  

[6]	J. M. Siskind, B, A. Pearlmutter, “Efficient Implementation of a
Higher-Order Language with Built-In AD”, AD2016 – 7th International
Conference on Algorithmic Differentiation (Oxford, UK, 2016)

[7]	T. F. Coleman, A. Verma, “ADMAT: An automatic differentiation toolbox
for MATLAB”, Proceedings of the SIAM Workshop on Object Oriented Methods for
Inter-Operable Scientific and Engineering Computing, v. 2, SIAM, Philadelphia,
PA, 1998  

[8]	A. Abad, R. Barrio, M. Marco-Buzunariz, M. Rodríguez, “Automatic
implementation of the numerical Taylor series method”, Appl. Math. Comput.,
268 (2015), 227–245  

[9]	“Automatic implementation of the numerical Taylor series method: A
Mathematica and Sage approach”, Appl. Math. Comput., 268 (2015), 227–245

