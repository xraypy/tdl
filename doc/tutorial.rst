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

Tdl executes interactively for a *shell*, which can be a useful way to test
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


Numeric Data and Arrays
-----------------------

For numeric processing, Tdl supports both scalar and numerical arrays
(representing a sequence of numeric values all with the same data type) for
integers, floating point numbers, and complex numbers (with real and
imaginary parts).  Numerical arrays can be multi-dimensional and can have
their shape and ordering of elements manipulated easily.  One simple way to
make an array is to use a comma-separated list inside delimiters '[', and
']', for example::

   tdl> arr = [1,2,3,4]
   tdl> print arr
   [ 1.  2.  3.  4.]

Arrays are treated as single objects, and expressions and built-in
functions act on each of their elements as a whole::
   
   tdl> print arr + 1
   [ 2.  3.  4.  5.]


Many primitive algebraic and trigonometric functions are built-in to Tdl,
including `sin`, `cos`, `tan` for sine, cosine, and tangent functions,
`arcsin`, `arccos`, `arctan` for the inverse of these, `log`, `log10`, and
`exp` for logarithm (base e), logarithm base 10, and exponential, and
`sqrt` for square-root.  There are actually several more built in
mathematical functions -- see the Reference Manual for a complete list.
For arrays, these functions also work on the entire array::


   tdl> print sqrt(arr)
   [ 1. 1.41421356 1.73205081 2. ]

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


   tdl> a = [1,2,3,4,5,6,7,8,9,10]
   tdl> print a[2:5]
   [ 3.  4.  5.]
  
Leaving out the first or second index implies that the slice should start
at the beginning, or stop at the end of the array::

   tdl> print a[:5]
   [ 1.  2. 3.  4.  5.]

   tdl> print a[8:]
   [ 9. 10.]

Note that a slice with starting and ending index differing by 1 will return
a one-element array::

   tdl> print a[8:9]
   [ 9.]

We note a few special cases here.  First, indices in a slice that are out
of range for that array are re-set to be the boundaries of the array
(though remember that arrays can be indexed with negative numbers).
Second, an array with one element is a slightly different from a scalar
value::

   tdl> print a[8]
   9.

and is something of a special case.  In many cases, you'll be able to
safely ignore this distinction, and one-element arrays often act just like
scalars. 


An array slice may also include a third argument for the index step or
*stride* of the slice, allowing you to take every nth element of an array::

   tdl> print a[::2]
   [ 0.  2.  4.  6.  8.]
   tdl> print a[1::2]
   [ 1.  3.  5.  7.  9.]
   tdl> print a[1::3]
   [ 1.  4.  7.]


Non-Numeric Data Types
----------------------


Strings

Lists

Dictionaries


String Formatting
-----------------




Conditional Statements
----------------------



Program Flow Control
--------------------

