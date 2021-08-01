1. Вопрос: "When should I use syntax/loc instead of #' (aka syntax)?"
   Ответ: https://stackoverflow.com/questions/10862665/when-should-i-use-syntax-loc-instead-of-aka-syntax

2. Что значат "\`" (quasiquote), "," (unquote), #` (quasisyntax), #, (unsyntax)? В чем отличие от '?
   Ответ: 

3. **ЗАПОМНИ!**
   * При использовании pattern variables в шаблонах, вместо их имен подставляется их значение.
   * 

4. Что такое "hygiene macro"? Какие есть подводные камни во всех этих лексических контекстах?

5. Почему ```(type (format "~a" (syntax-property #'expr-expanded-typecheck 'landau-type)))``` в with-syntax становится pattern variable?
   Это же просто строка...
   Ответ: ...

6. Всегда ли можно заменить всякие квазиквотинги на обычные? А в коде Landau?
