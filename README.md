# Semi-analytical differentiation

This package is designed with the purpose of computing higher order analytical derivatives. The primary reasons for its
exsistence are
  * Hyperdual numbers become unfeasible somewhere around fourth order derivatives (both regarding implementation and runtime)
  * Symbolic differentiation is way too slow.

Thus, the point of this package is that we want to be able to define explicit functions like

$$
f(x) = (sin^2(x^2 + exp(x)) - sin(cos(x)) * log(x)^2 - (x^3 + x^4)/(2 + x - x^2))
$$

and compute what ever derivative we wish, efficiently and without the overhead associated with full symbolic differentiation.
We don't care what the analytical expressions for the derivatives actually look like, as long as we can evaluate them.

For mathematical details, see the memo.

# Structure

The core of the package is the C++ implementation (`polynomials.[h/cpp]`), this is the only version with full functionality.

Some of the funcitonality has been ported to the Fortran implementation (`polynomials.f90`), and some has been exposed to
python (`polynomial_pybind.cpp`/`libpyanalytical.pyi`).

# Building and running the demos

To build the package:
```bash
git clone https://github.com/vegardjervell/analytical_derivatives.git # Clone this package
cd analytical_derivatives
git submodule update --init --recursive # Fetches pybind11 
python3 -m venv venv # Set up a python virtual environment
source venv/bin/activate
mkdir build
cd build
cmake ..
make install
```
Modify the above as needed if you don't need the python wrapper (you will need to remove a couple lines from `CMakeLists.txt`).

Once the package is built:
```bash
# Still in build directory
./fpoly # Run Fortran demo
./cpoly # Run C++ demo
cd .. && python derivatives.py # Run Python demo
```

# Core capabilities

The functionality exposed in all three languages is
  * The `Polynomial` class - Representing Laurent polynomials.
  * The `PolyExp` class - Representing functions of the form $f(x) * exp[g(x)]$, where $f$ and $g$ are Laurent polynomials.
  * The `PolyFrac` class - Representing functions of the form $f(x) / g(x)$, where $f$ and $g$ are Laurent polynomials.

These are special cases that show up often, and have been optimized as compared to the generalized functionality found
only in the C++ module:
  * `Analytical` - Abstract parent class representing some arbitrary function for which the $n$'th derivative is known.
    * `sin`, `cos`, `pow`, `log`, and `exp` - Analytical representation of these functions.
  * `Sum` - A sum of `Analytical` objects
  * `Product` - A Product of two `Analytical` objects
  * `Composed` - A composition of two `Analytical` objects, for example:
    * $exp(x^{-2})$ is represented by `Composed(exp(), pow(-2))`
    * $sin^2(cos(x))$ is represented by `Composed(pow(2), Composed(sin(), cos()))`
  * Operator overloads for `+`, `-`, `*`, `/`, `+=`, and `*=`, that ensure "ordinary" math can be done on the objects in order to generate `Sum`'s and `Product`'s.

All the above is contained in the `analytical` namespace to avoid conflicts for the common function names. Be careful about
`using namespace analytical` if you're also `include`'ing `cmath`, because `cmath` puts stuff in the global namespace.

# Usage

Take a look at the demo files.