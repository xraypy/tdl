==============================
Tiny Data Language: A Tutorial
==============================


Introduction: What is TDL?
--------------------------

The Tiny Data Language (Tdl for short) is a simple, easy to use computer
language for numerical processing of scientific data.  While Tdl is easy to
use, it is also a full-featured programming language and provides access to
fast numerical algorithms for number crunching.  In addition, it is easy to
add new or custom functionality to Tdl.

Tdl is aimed at scientists needing to process and analyze data.  By itself,
Tdl is not a complete 'canned' package for analyzing certain types of
scientific data, and it is not a general-purpose programming language.
Instead, it tries to makes it very easy to manipulate data, either
interactively or by writing programs.

Given the aim for scientific computing, this tutorial presents Tdl with
some assumptions about computer skills, such that the user understands how
to run a program.  Some familiarity with computer programming is probably
also assumed, so that the reader should be aware that programming languages 
have distinct data types (integer, floating point or real numbers, strings,
and also arrays of these) as well as control structures like 'if-then-else'
statements, for-loops and sub-programs or functions.

TDL and Python
--------------

Tdl is written in Python, an elegant and dynamic high-level object-oriented
general purpose programming language.  Tdl makes heavy use of Python's
dynamic capabilities for its implementation, and borrows much of its syntax
from Python.  Tdl also relies on a few fantastic extensions to Python for
numerical processing: NumPy (which adds fast numerical arrays to Python),
SciPy (which adds interfaces to many sophisticated and well-tested
numerical algorithms using NumPy arrays), and MatPlotlib (which adds
beautiful plotting capabilities to Python).  Extensions to Tdl are written
Python, and Tdl can be called from Python code.

Tdl does have several differences from Python, primarily reflecting the
fact that tdl aims at providing a simple approach to numerical processing
of data, while Python is a general purpose language.  Some of the major
differences are then:

#. Python is object-oriented, while Tdl is not.

#. Python code needs to explicit import functionality, while Tdl provides
   array processing, plotting, and many numerical processing routines
   built in.

#. Python uses namespaces to distinguish different functional
   capabilities, while Tdl uses namespaces to distinguish both functionality
   and groups of data.

#. Python provides much more introspection and sophisticated approaches
   for dynamic programming.

#. Tdl allows a 'defined variable', which is a variable that acts
   somewhat more like a function than a normal variable.

These and other other distinctions between Tdl and Python are further
discussed below and in companion document ExtendingTdl.



TDL Overview
------------

Tdl is an interpreted procedural language that should look pretty familiar
to most people with some programming experience.  Named variables can be
used to hold a variety of native data types, including both numeric and
non-numeric types of data. Simple expressions are used for calculations and
assignment, and conditionals and control-flow are supported.  This section
will discuss the basics of data types, variable naming, and program flow
and execution.

----------------------------------
using the interactive interpreter
----------------------------------

Tdl executes interactively in a *shell*, which can be a useful way to test
and learn about tdl, as well as a simple way to work with your data.  You
may also put a set or `program` of tdl statements into a file and execute
that.  Throughout this tutorial, we'll use the convention of showing a Tdl
session, in which you type at the `tdl>` prompt.  We encourage you to try
many of these examples as you read this tutorial.

Tdl supports simple evaluation of numerical and non-numerical expressions,
so that you can simply type::

   tdl> 1 + 5
   6.0 

One of the most important tdl functions is `print`, which simply prints the
result of evaluating each of its arguments::
   
   tdl> print 1,  15/3., 'hello'
   1.0 5.0 hello

Named variables can be used to hold values of expressions, so that
assignment statement like::

   tdl> a = 1
   tdl  b = a/3.
   tdl> print a, 2*a, b
   1.0 2.0 0.333333333333


-------------
getting help
-------------

The Tdl interpreter has a built-in help system.....


Numeric Data and Arrays
-----------------------

Tdl is designed for numeric processing, and the ability to manipulate and
do calculations with numerical data is central to how the language works.
In this section, we'll go over the basics of how Tdl thinks about numbers
and arrays (or vectors) of numbers.  Tdl borrows almost all the
functionality hear from Python and the wonderful NumPy extension which
provides a basic 'numerical array' type to Python.

------------------------
real and complex numbers
------------------------

For numeric processing, Tdl has built-in data types for real (or floating
point) numbers, and complex numbers (which have both real and imaginary
parts).  In Tdl, a real number is represented with 8-bytes (double
precision) .  In contrast to many programming environments, Tdl generally
does not assume that numbers are integers, but rather assumes that they are
real, and allows you to convert numbers to integers or does it interally.
We've seen simple numerical values above.  Complex numbers are represented
using 'j' as the base of imaginary numbers.  Thus numbers look like::

   tdl> imag = 0 + 1.j
   tdl> print  image * image
   (-1+0j)

Tdl includes many primitive algebraic and trigonometric functions,
including `sin`, `cos`, `tan` for sine, cosine, and tangent functions,
`arcsin`, `arccos`, `arctan` for the inverse of these, `log`, `log10`, and
`exp` for logarithm (base e), logarithm base 10, and exponential, and
`sqrt` for square-root.  There are actually several more built in
mathematical functions -- see the Reference Manual for a complete list.
For complex numbers, these functions will return the appropriate complex
result, but you must explicit use complex numbers to get such a result.
That is,::

   tdl> sqrt(-1)
   evaluation error: Math Error (undefined number) at line:
     'sqrt(-1)'

indicating an error as occured, whereas::

   tdl> sqrt(complex(-1))
   1j
   tdl> sqrt(-1 + 0.j)
   1j

where the builtin function 'complex' automatically coerces otherwise purely
real numbers to become complex values, essentially by adding '0.j' to
them.  The second version above is more explicit: -1+0.j is a complex number.

-----------------
numerical arrays
-----------------

In addition to supporting scalars (individual numbers), Tdl supports
multi-dimensional numerical arrays, which represent a sequence of numeric
values all with the same data type.  Arrays may be of any size and
dimension, and each element may be an integer, floating point number, or
complex number..  The dimension and "shape" of numerical arrays and
ordering of elements can be manipulated easily.  

A simple way to make an array is to use a comma-separated list inside
delimiters '[', and ']', for example::

   tdl> arr = [1,2,3,4,5,6,7,8,9,10]
   tdl> print arr
   [ 1.  2.  3.  4.  5.  6.  7.  8.  9.  10.]

In addition, there are a few builtin functions for creating arrays,
including 'range' which generates a sequence of numbers.  The above array
could have been created with::

   tdl> arr = range(10)+1

Note that 'range(10)' would have made the first element of the array 0,
which is the convention used in Tdl.   The 'range' function has many
options to specify the starting and stopping values of sequence, the step
interval to use, and the shape of the resulting array.

Arrays are treated as single objects, and expressions and built-in
functions act on each of their elements as a whole::
   
   tdl> print arr + 1
   [ 2.  3.  4.  5.  6.  7.  8.  9.  10.  11.]

For the many primitive algebraic and trigonometric functions  built-in
to Tdl, including `sin`, `cos`, `tan` for sine, cosine, and tangent
functions, `arcsin`, `arccos`, `arctan` for the inverse of these, `log`,
`log10`, and `exp` for logarithm (base e), logarithm base 10, and
exponential, and `sqrt` for square-root.  There are actually several more
built in mathematical functions -- see the Reference Manual for a complete
list.  For arrays, these functions also work on the entire array::


   tdl> print sqrt(arr)
   [ 1.          1.41421356  1.73205081  2.          2.23606798  2.44948974
     2.64575131  2.82842712  3.          3.16227766]

Following the convention of many programming languages, arrays can be
indexed by position , with index 0 being the first element in the array,
using a notation using square brackets::

   tdl> print arr[1]
   2.0

The index can be a full expression, of course,::

   tdl> i = 4
   tdl> print arr[i-2]
   3.0
 
If the index is negative, the indexing will count backward from the last
element in the array, with index=-1 meaning 'the last element'::

   tdl> print arr[-1]
   4.0

If the index is beyond the length of an array, an 'index out of bonds'
error will be generated.

A sub-array can be taken from a starting array using a *slice*.  The
indexing syntax above is the simple form of a slice.  Other slices extend
this syntax to include 1 or 2 colons to separate starting and stopping
indices, and an optional index step.  Ignoring the optional step for a
moment, a subarray can be taken like this::

   tdl> print arr[2:5]
   [ 3.  4.  5.]
  
Leaving out the first or second index implies that the slice should start
at the beginning, or stop at the end of the array::

   tdl> print arr[:5]
   [ 1.  2. 3.  4.  5.]

   tdl> print arr[8:]
   [ 9. 10.]

Note that a slice with starting and ending index differing by 1 will return
a one-element array::

   tdl> print arr[8:9]
   [ 9.]

We note a few special cases here.  First, indices in a slice that are out
of range for that array are re-set to be the boundaries of the array
(though remember that arrays can be indexed with negative numbers).
Second, an array with one element is a slightly different from a scalar
value::

   tdl> print arr[8]
   9.

and is something of a special case.  In many cases, you'll be able to
safely ignore this distinction, and one-element arrays often act just like
scalars. 


An array slice may also include a third argument for the index step or
*stride* of the slice, allowing you to take every nth element of an array::

   tdl> print arr[::2]
   [ 0.  2.  4.  6.  8.]
   tdl> print arr[1::2]
   [ 1.  3.  5.  7.  9.]
   tdl> print arr[1::3]
   [ 1.  4.  7.]


--------------------------
multi-dimensional arrays
--------------------------

Arrays in Tdl can have any number of elements, and can have multiple
dimensions.  The layout of an array (the number of dimensions and the
length of each dimension) is called its *shape*.  The shape of array 
can be specified during array creation (for example, with the range()
function) or modified with the reshape() function.  A literal
multidimensional array uses multiple, nested brackets::

   tdl> md = [[1,2,3,4,5],[6,5,4,3,2],[9,7,5,3,1]]


Slices for multi-dimensional arrays can be complicated.


Non-Numeric Data Types
----------------------

Although it is primarily a language for processing numerical data, Tdl has
a few data types for storing non-numerical data.  Compared to low level
languages like C or Fortran, these data types act as simple data
structures.  As such, these data types are extremly when handling real
data, and in writing programs.

-------------------
logical or booleans
-------------------

Tdl supports boolean values (which have either the value True or the value
False).  These are most useful in combination with logical expressions in
conditional expressions (if-then-else) and while loops, as we will see
below.  


---------
strings
---------


The simplest non-numeric data type is a *string*, which is a sequence of
text characters.  Strings in Tdl are nearly identical to their python
counterpart, and typically specified as literal strings::
    
   tdl> greeting = 'Hello, World'
 
Literal strings can be enclosed either in single quotes or in double
quotes.  There is a difference between single and double quotes, which is
the handling of "escape characters".

Triple-quote strings can span multiple lines, and the newline characters
are preserved.

-----------
lists
-----------

The second non-numeric data type in Tdl is a *list*.  A list is similar to
an array in that it contains a sequence of values.  Unlike an numerical
array, in which all elements must contain the same type of data, a list is
a sequence of values that may be of different types -- numbers, arrays,
strings, other lists, or dictionaries (introduced below). In tdl, a list
looks quite a bit like an array, and the same syntax is used to create
them.  When creating a list, Tdl automatically checks whether it might be
an array and, if so, automatically converts it.  For this reason, the
difference between lists and arrays should be well understood.  As a simple
example:: 

   tdl> mylist = [1,'H', 1.008, 'Hydrogen']

is a list, as it consists of both numbers and strings.  The following are
also lists::

   tdl> list1 = [1,2,3,'a']
   tdl> list2 = [[1,2,3],[5,6,'']]
   tdl> list3 = [[1,2,3,4],[5,6,7,8,9]]

The last one is a list becaues it would not have rows of equal size as an
array.

like arrays, the values contained in a list can be 'indexed' using the
integer position of the element.....
  
List functions (append, etc).


--------------
dictionaries
--------------

Tdl supports a dictionary data type. Like lists, dictionaries are
collections of 

------------------
other data types
------------------

Tdl does actually allow other data types, and essentially any 'Python
object' can be held in a named Tdl variable.  This can be quite useful when
writing extensions.  The drawback is that most of Tdl won't know what to do
with that variable -- you'll have to be careful to use it only in functions
that know what to do with it.  One very common, builtin data type that is
not one of the ones listed above is a 'file handle', ....

------------------


Python notes:  In the above discussion, Tdl shows a few differences with
Python.

#. Python does not distinguish single and double quote strings. Tdl does
   make a distinction between these for how escaped strings are processed.

#. Python has tuples (which are sequences like lists, but that have a
   fixed size.  This is not so much of an issue in non-object-oriented Tdl,
   and so to avoid a 3rd sequence type, Tdl does not have tuples.  In Python,
   tuples are used many places, including string formatting and return values
   from functions.  Tdl uses lists for these functionalities.

#. Python has a strong distinction between 'immutable' types, and allows 
   dictionary keys to be any immutable types.  For simple uses, strings are
   the most common and natural type for dictionary keys, so Tdl just enforces
   this. 



Variables, Name Spaces and Groups
----------------------------------


In all the above examples, we've assumed that all the named program
variables are global variables

Operators, Expressions, and Statements
----------------------------------------

Thugh we've already seen several examples of using Tdl and have talked
about data types and variable names, at this point we should take a step
back and discuss the basic structure of Tdl, including how values are
assigned to variables and how calculations are done.

A Tdl script or program consists of a list of 'statements'.  Usually a
statement corresponds to one line of text, but they can extend beyond one
linee.

There are several types of statements in Tdl.  The most common kind is the
'assignment statement', which assigns some value to a named variable. That
is, the statement::
    tdl>  x = 1
 
creates 
String Formatting
-----------------




Conditional Statements
----------------------

In complex scripts, it is ofetern desirable or even necessary to execute
some part of the script only under certain conditions, or to change the
calculation done.  This is done with conditional statements or
'if-then-else' statements.  These are so vital to scripts that Tdl supports
a few different variations on these.  The simplest form is a one-line 'if
statement', that looks like this:
    
    if <condition>:  <single-line statement>

The <condition> here is a boolen expression.  If it True, the statement
will be run, if the condition is False, it will not be run.  Note that the
colon ':' after the condition is important.  A simple example::
  
    if x == 0:  print 'x is zero!'

If you need to execute more than one line of a script, you can use the
'if-endif' form::
  
    if <condition>:
       <block>
    endif

Here <block> is an arbitrarily long set of Tdl statements (one statement
per line), ending with an 'endif'.  A simple example of this form::

    if x == 0:
       logx = 0.
       print 'warning: set logx to zero!'
    endif



Program Flow Control
--------------------

It is often necessary to repeat some calculation or nearly-the-same
calculation.  Tdl provides two mechanisms for doing this: For loops and
While loops.

----------
For loops
----------

A for loop executes a block of Tdl statements repeatedly.  One value is
changed for each iteration of the block using a pre-determined list of 
values to take.  A generic for loop looks like this::

    for <value> in <list>:
        <block>
    endfor 

Because a list of values is used, a for loop will execute a predictable
number of times -- the length of the list.  The most common use of a for
loop is to increment some index that indicates the number of times through
the block.  Here, we use the 'range()' function::

    for x in range(10):
        print x
    endfor

Note that the value x is put into the current default namespace and 
exists when the loop finishes. 

------------
While loops
------------

A while loop executes a block of Tdl statements as long as some logical
(boolean) condition is met.  The generic while loop looks like this::

    while <condition>:
        <block>
    endwhile

The <condition> here is a boolean, and as long as the condition is True,
the block will be evaluated.  The expectation is that the block will do
something to alter the condition.  A simple example::
    x = 0
    while x < 10:
       print 'x = ', x
       x = x + 1
    endwhile

This form could just as well be implemented with a for-loop::
    for x in range(10):
       print 'x = ', x
       x = x + 1
    endfor

but many while-loops would be more difficult to implement with for loops::

    x = 0
    y = 9
    while x < 10 and y > 0:
         x = x + 1
         y = 10 - exp(x/3.)
         print x, y 
    endwhile


--------------------
break, continue
--------------------

It is sometimes necessary to interrupt the execution of the block of a for-
or while-loop.  Tdl has two mechanisms to do this.  The first is the
'break' statement, which jumps completely out of the loop as soon as it is
encoutered.  With a break statement, the above while loop could be
re-implemented as::

    x = 0
    y = 9
    while x < 10:
         x = x + 1
         y = 10 - exp(x/3.)
	 if y < 0: break
         print x, y 
    endwhile

In contrast to 'break', the 'continue' statement stops the evaluation of
the current block, but continues on to the next step in the loop.


Exceptions and Errors
---------------------

try / except


Builtin Functions
-----------------

-------------------------------
Basic Mathematical Functions
-------------------------------

-----------------------------------------
Functions about the state of the program
-----------------------------------------

---------------------------
Functions for getting help
---------------------------

User-Defined Functions
----------------------

argument lists

isolated namespace

return values


Defined Variables
------------------






Extending with Python
---------------------

