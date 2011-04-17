The Levenberg-Marquardt algorithm is an iterative technique that
finds a local minimum of a function that is expressed as the sum of
squares of nonlinear functions. It has become a standard technique
for nonlinear least-squares problems and can be thought of as a
combination of steepest descent and the Gauss-Newton method. When
the current solution is far from the correct one, the algorithm
behaves like a steepest descent method: slow, but guaranteed to
converge. When the current solution is close to the correct
solution, it becomes a Gauss-Newton method.

Optional box- and linear constraints can be given. Both single and
double precision floating point types are supported.

The actual algorithm is implemented in a [C library] which is bundled
with [bindings-levmar] which this package depends on.

License
=======

This library depends on [bindings-levmar] which is bundled together
with a [C library] which falls under the GPL. Please be aware of this
when distributing programs linked with this library. For details see
the description and license of [bindings-levmar].

[bindings-levmar]: http://hackage.haskell.org/package/bindings-levmar
[C library]:       http://www.ics.forth.gr/~lourakis/levmar
