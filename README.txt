This is TDL -- a Tiny Data Language for numerical processing.

Notes and Bugs
--------------

4-8-06 T2
- In order for a fcn__cmdout to recieve additional keyword arguments, 
two small changes are needed: 1) in Expression.py the call to __cmdout__
needs to pass **kws, and the function it self must have **kws defined.  
Is there a better way to do that (see the ls example).

- Unfortunatley still get funky output for path like strings on windows... 