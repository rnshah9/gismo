#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsXmlCollection.h>
#include <gsIO/gsXmlCollection.hpp>

#ifdef GISMO_BUILD_PYBIND11
#include <gsCore/gsMultiPatch.h>
#endif

namespace gismo
{

#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;

void pybind11_init_gsXmlCollection(py::module &m)
{

  using Class = gsXmlCollection;
  py::class_<Class>fd(m, "gsXmlCollection");

  fd.def(py::init<const std::string&>())

  // Member functions
  .def("save", (void (Class::*) ()) &Class::save)
  .def("addFile", static_cast<void (Class::*)(std::string const & ) > (&Class::addFile), "Add String filename to Xml Collection")
  .def("addFile", static_cast<void (Class::*)(std::string const &, std::string const & ) > (&Class::addFile), "Add String filename with label to Xml Collection")

  .def("getId", static_cast<void (Class::*)(const int &, gsMultiPatch<real_t> &) const > (&Class::getId), "Get from file with id the MultiPatch")
  .def("getLabel", static_cast<void (Class::*)(const std::string &, gsMultiPatch<real_t> &) const > (&Class::getLabel), "Get from file with label the MultiPatch")
  .def("getMatrix", &Class::getMatrix, "Get from file with label the Matrix")

  ;
}

#endif

}


