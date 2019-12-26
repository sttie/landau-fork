# Description
Landau is a *lan*guage for *d*ynamical systems with *au*tomatic differentiation.

Landau is a Turing incomplete statically typed domain-specific differentiable language. 
The Turing incompleteness provides 
the ability of sophisticated source code analysis and, as a result, a highly optimized compiled code. 
Among other things, the language syntax supports functions, compile-time ranged for loops, if/else 
branching constructions, real variables and arrays, and the ability to manually discard calculation 
where the automatic derivatives values are expected to be negligibly small. In spite of reasonable 
restrictions, the language is rich enough to express and differentiate any cumbersome paper-equation 
with practically no effort.

## Automatic differentiation

The automatic differentiation (forward mode) technique can be built on the concept of dual numbers.
Dual numbers can be thought as a pair $`\langle x, x' \rangle`$ of a numerical value
and its numerical derivative. Arithmetic applies to the right part according to
differentiation rules:

```math
\begin{aligned}
\langle u, u' \rangle \pm \langle v, v' \rangle =& \langle u \pm v, u' \pm v' \rangle\\
\langle u, u' \rangle \langle v, v' \rangle =& \langle uv, u'v + uv' \rangle\\
\langle u, u' \rangle / \langle v, v' \rangle =& \left\langle \frac{u}{v}, \frac{u'v - uv'}{v^2} \right\rangle, v \neq 0\\
\end{aligned}
```

However, the manual usage of dual numbers makes the implementation of the equations more
tedious and harder to maintain.

## Problem formulation

The task is to calculate the (automatic) derivatives of functions $`f_i(x_1,\ldots,x_j)`$
w.r.t. parameters $`\{p_k\}`$. Numerical values of $`\frac{\mathrm{d} x_j}{\mathrm{d} p_k}`$
are given externally.

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

Landau uses a hybrid approach of derevatives calculation: symbolic differentiation 
at the expression level and automatic (forward mode) at the assignation level.
Consider the following peice of Landau code:
```python
      real[k] g, f
      g[:] = 2 * sin(x[:]) + x[:]
      for i = [0 : k]
        f[i] = 2 * g[i] + g[i] * x[i]
```
The left-hand side of the assignation is differentiated symbolically.
The result of computation is stored in the corresponding array (`dgdx0`).
Whenever the result of $`\mathrm{d}g/\mathrm{d}x_0`$ is needed, it is
taken from the corresponding array.
```python
      dgdx0[..] = 2 * cos(x[..]) * dxdx0[..] + dxdx0[..]
      g[:] = 2 * sin(x[:]) + x[:]
      for i = [0 : k]
        dfdx0[..] = 2 * dgdx0[..] 
                      + dgdx0[..] * x[i] + g[i] * dxdx0[..]
        f[i] = 2 * g[i] + g[i] * x[i]
```
## Implementation

The Landau transpiler is built on top of the [Racket](https://racket-lang.org) platform and consists of two main parts: the syntax analyzer 
and the code generator. The job of the syntax analyzer is:

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

- generates highly efficient C code (or not-so-efficient Racket syntax
objects that are useful if you are programming in Racket);
- exempts user of C pointer managing trouble and low-level debugging;
- allows easy notation of mathematical equations, while also allowing loops and functions;
- and, of course, provides first-class automatic derivatives (while you can use Landau without them, too).

The Landau-generated C code is portable on different compilers on Windows, Linux or MacOS.

## Development in Landau

The model development consists of three stages:

1. Manual model description in Landau language (`model.dau`)
2. Transpiling to ANSI-C (code generation `model.c`)
3. Compiling to a shared library for a specific platform
  (`model.dll`, `model.so`, `model.dylib`)

## Example: Spacecraft movement

As an example, consider the Landau program for modeling spacecraft movement around a planet -- [spacecraft.dau](tests/spacecraft.dau).
Spacecraft's initial position and velocity, as well as the gravitational parameter of the planet, are supposed
to be determined by an external implementation of the least-squares method. Derivatives of the right part of the movement equation
w.r.t. the initial state vector and the planet mass are handled by Landau.

## Publications

The work is published in papers:
 - I. Dolgakov, D. Pavlov. “Landau: language for dynamical systems with automatic differentiation”. Zap. Nauch. Sem. POMI **485** (2019), 78--89. [Preprint in English](https://arxiv.org/abs/1905.10206).
 - I. Dolgakov, D. Pavlov. “Using capabilities of Racket to implement a domain-specifc language” (accepted to a journal, coming soon).

# Installation

1. Download and install [the Racket](https://download.racket-lang.org/).
2. Clone repository: `$ git clone https://gitlab.iaaras.ru/dolgakov/Landau.git`.
3. Link Landau as a package: `$ cd autodiff/ && raco link landau`.

# Configuraiton

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
