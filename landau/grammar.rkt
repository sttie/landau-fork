#lang ragg
program : constant* parameter* (func | expr)*

constant : 'const' type IDENTIFIER '=' (expr | "{" parlist "}")

expr : term
     | expr '+' term
     | expr '-' term

term : factor
     | term '*' factor
     | term '/' factor

factor : primary '^' factor 
       | primary

primary : unop primary 
        | element

element : number 
        | '(' expr ')' 
        | get-value 
        | func-call

number : INTEGER | FLOAT
unop : '-' | '+'

if-expr : "if" bool-expr "{" expr-body "}" [ "else" "{" expr-body "}" ]
single-term : var-decl | assignation | if-expr | for-expr | der-annot | der-apply | discard | print
expr-body : (var-decl | assignation | if-expr | for-expr | der-annot | der-apply | discard | print)*
; NOTE: It is needed for allocation parameter's derivative arrays
func-body : (var-decl | assignation | if-expr | for-expr | der-annot | der-apply | discard | print)*

;; NOTE: Hack to return a sequence of syntax objects after expansion.
decl-block : (var-decl | assignation | func | expr-body)*

print : "print" IDENTIFIER (expr | get-derivative) | "print" (expr | get-derivative)

get-derivative : get-value "'" get-value

bool-expr : bool-term | bool-expr "or" bool-expr
bool-term : bool-factor | bool-term "and" bool-term
bool-factor : "(" bool-expr ")" | "(" expr comp-op expr ")" | "not" bool-factor

comp-op : ">" | "<" | ">=" | "<=" | "==" | "!="

for-expr : "for" loop-var-decl ("{" expr-body "}" | single-term)
loop-var-decl : IDENTIFIER "=" "[" expr ":" expr "]"

; TODO: make slice a regular syntax
der-annot : get-value "'" get-value "=" expr
der-apply : get-value "=" get-value "'" get-value
discard : "discard" get-value "'" get-value
parameter : "parameter" "[" expr "]" IDENTIFIER

func: type IDENTIFIER '(' arglist ')' '{' func-body '}'

arglist : [ argument other-argument* ]
other-argument : ',' type IDENTIFIER
argument : type IDENTIFIER

basic-type : 'real' | 'int'
array-type : basic-type '[' expr ']'

type : basic-type | array-type

get-value: IDENTIFIER
         | IDENTIFIER "[" expr "]"
         | IDENTIFIER "[" expr ':' expr "]"
         | IDENTIFIER "[" ':' expr "]"
         | IDENTIFIER "[" expr ':' "]"
         | IDENTIFIER "[" ':' "]"

func-call: IDENTIFIER "(" parlist ")"
parlist : [ par other-par* ]
other-par : ',' expr
par : expr

mut-assign: "+=" | "-=" | "*=" | "/="

var-decl: type IDENTIFIER
        | type IDENTIFIER ("=" | mut-assign) expr

assignation: IDENTIFIER ("=" | mut-assign) expr
           | IDENTIFIER "[" expr "]" ("=" | mut-assign) expr
           | IDENTIFIER "[" expr ":" expr "]" ("=" | mut-assign) expr
           | IDENTIFIER "[" ":" expr "]" ("=" | mut-assign) expr
           | IDENTIFIER "[" expr ":" "]" ("=" | mut-assign) expr
           | IDENTIFIER "[" ":" "]" ("=" | mut-assign) expr
