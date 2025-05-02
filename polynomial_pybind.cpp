#include "polynomials.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

namespace py = pybind11;
using namespace analytical;

PYBIND11_MODULE(libpyanalytical, handle){
    handle.doc() = "The best docs. The very best.";

    handle.def("get_partitions", py::overload_cast<int>(&get_partitions));
    handle.def("partition_multiplicity", &partition_multiplicity);

    py::class_<Analytical>(handle, "Analytical")
        .def("derivative", &Analytical::derivative)
        .def("numerical_derivative", &Analytical::numerical_derivative)
        ;

    py::class_<Polynomial, Analytical>(handle, "Polynomial")
        .def(py::init<int, int, std::vector<double>, int>())
        ;

    py::class_<PolyExp, Analytical>(handle, "PolyExp")
        .def(py::init<Polynomial, Polynomial>())
        ;

    py::class_<PolyFrac, Analytical>(handle, "PolyFrac")
        .def(py::init<Polynomial, Polynomial>())
        ;
}