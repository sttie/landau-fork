1. 
```
real[k * k] func(real[k] x, real[k * k] dfdx0) {
    func = dfdx0
}

assignation: the right-hand expression is a slice, but the left-hand one is not. A slice must be assigned to a slice
in: (assignation func "=" (expr (term (factor (primary (element (get-value dfdx0)))))))
```

Тупая ошибка в том, что rhs - не slice. Да и в принципе, почему бы не разрешить array1 = array2?

2. 
```
real[k * k] func(real[k] x, real[k * k] dfdx0) {
    func[:] = dfdx0
}
```
Ошибок здесь нет, но такой кейс генерит плохой код вида double_array1[i] = double_array2

3. Такой код
```
const int k = 3
parameter[k] x0

real[k * k] func(real[k] x, real[k * k] dfdx0) {
    x[:] ' x0[:] = dfdx0[:]
    real[k] f = 2 * x[:]
    func[:] = f[:]
}
```

вызывает ошибку
```
t-quote-syntax: undefined;
 cannot reference an identifier before its definition
  in module: top-level
  context...:
   body of top-level
   /usr/local/Cellar/minimal-racket/8.2/share/racket/collects/syntax/parse/private/parse.rkt:1020:13: dots-loop
   /usr/local/Cellar/minimal-racket/8.2/share/racket/collects/racket/private/stxparam.rkt:61:2
   /Users/sttiemath/projects/landau-fork/landau/semantics.rkt:346:0
   /usr/local/Cellar/minimal-racket/8.2/share/racket/collects/syntax/wrap-modbeg.rkt:46:4
```

но если исправить последнюю строку на ```func[:] = f[:] ' x0[:]```, то все починится

4. В коде
```
#lang landau

const int k = 2
parameter[k] x0

real[k * k] func(real[k] x, real[k * k] dfdx0) {
    # x[:] ' x0[:] = dfdx0[:]
    real[k] f = 2 * x[:]
    func[:] = f[:] ' x0[:]
}
```

func[:] = f[:] ' x0[:] никак не компилируется в С без der-annot. ЭТО СТРАННО! Нужно поддержать хотя бы ошибку

5. something размера 4, x размера 4, выражения something[0:2] = x[0:3] и something[0:3] = x[0:2] выдают
   ошибку "cannot reference an identifier before its definition". something[0:2] = x[0:2] работает нормально

6. Стопроцентный баг: 
   x[1:5] парсится в ```(#%app c-array-ref (quote "x12") (#%app c+ (quote 1) (quote "slice_idx")))```
   x[0:4] парсится в ```(#%app c-array-ref (quote "x12") (#%app c+ (quote 0) (quote "slice_idx")))```
   Короче, берется только первый индекс... да еще и c-array-ref
   А еще тупость в том, что получается (c+ 0 "slice_idx"), то есть число плюс строка)))

7. При local-expand-memo не производится никаких проверок на легальность операций. Например, print x + something + 2 - нормальный код,
   который спокойно доходит до генерации Си-кода. Но вообще, проверки производиться **должны**. Наверно, раскрывать expr нужно по-другому.
   Непонятно, где происходит проверка типов в выражении. В backrun или semantics?