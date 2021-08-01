1. Зачем код в функции оборачивается в отдельный скоуп с помощью { ... }?
2. 
```
real[k * k] func(real[k] x, real[k * k] dfdx0) {
    x[:] ' x0[:] = dfdx0[:]
}
```
В скомпилированном коде на Си нет никакого кода в func. Это что, не инструкция, а декларация для компилятора?


2. У вас пример из-за sin(x[:]) не компилируется...
```
parameter[k] x0 # Declare parameters
                # w.r.t which derivatives should be taken
real[k * k] func(real[k] x, real[k * k] dxdx0) { 
    x[:] ' x0[:] = dxdx0[:] # Tell the compiler 
                            # about dx/dx0 jacobian

    real[k] f = 2 * sin(x[:]) + x[:] # Compose the function's body

    func[:] = f[:] ' x0[:] # Write the numerical value
                           # of df/dx0 jacobian 
                           # to the return array
}
```


3. Что записывается в переменную x0 в ```x[:] ' x0[:] = dfdx0[:]```, где x0 == parameter[k]?

4. Что за магия происходит в сишном коде при транспиляции чего-то типа ```func[:] = f[:] ' x0[:]```? ОЧЕНЬ много непонятного кода.
    (классная идея - попробовать разобраться самому и скинуть Павлову и Долгакову свои догадки, чтобы они знали, что ты работаешь над проектом)

5. По какому принципу выбираются названия массивов типа f14 и т.д?
   Ответ: 14 - номер столбца в строке.

6. В чем отличия производных слева и справа от знака равно?

7. Что такое parameter?

8. Почему parameter'ы нельзя использовать в выражениях?

9. Почему x[:] ' x0[:] = dfdx0[:] - просто декларация (не компилируется в Си) и для чего она нужна?

10. Что означает debug в коде? У der-apply даже есть ```(displayln "der-apply: debug is not implemented")```

11. Что за типы 'dual-b, 'dual-l, 'dual-r, 'dual-b-derivative? Что они означают?

12. Что за файл combinators.rkt? Почему так называется?

13. В der-apply есть expanded-value, в print есть expr-expanded-typecheck (и наверняка это есть во множестве других мест).
    Почему все эти pattern-variables равны такой фигне: #<syntax:landau/semantics.rkt:1801:36 (#%app format (quote \"bug: nonsense in runtime\"))>?

14. Основной вопрос: что делает и зачем нужен файл metalang.rkt?
    Побочные вопросы: Что такое "metalang syntax" в шапке semantics.rkt? Как они связаны с кодом на Си? Как они переводятся в текст?
    Почему #, заставляет трансформироваться синтаксический объект с этими metalang-syntax-сущностями в ТЕКСТ-код на Си?

15. Что по сути своей такое pattern variable? А attribute? В чем их отличия?

16. Что такое parametrization? Что несет в себе current-parameterization?

17. В metalang.rkt повсюду возвращаемые значения - это  (syntax/loc stx str). Зачем это надо, если вернется обычная строка?

18. Зачем нужны _0.0, _0.5 и т.д. в metalang.rkt?

19. semantics.rkt:1765. Что за get-value.name?


## ОТВЕТЫ

6. В грамматике есть два нетерминала: ```der-annot : get-value "'" get-value "=" expr``` и ```der-apply : get-value "=" get-value "'" get-value```.

7. parameter[k] - это k параметров, относительно которых берется производная

8. Потому что parameter[i] - не число, а некоторый служебный для комилятора объект

17. Видимо (!), это хак для того, чтобы раскрытие макросов дало в итоге не синтаксический объект, а строку.
    Если это правда, то Павлов и Долгаков - очень умные люди.