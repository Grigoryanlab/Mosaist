/*
 * This inclusion should be put at the beginning.  It will include <Python.h>.
 */
#include <boost/python.hpp>
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include "msttypes.h"

/*
 * This is a macro Boost.Python provides to signify a Python extension module.
 */
BOOST_PYTHON_MODULE(mstpython) {
    // An established convention for using boost.python.
    using namespace boost::python;

    // expose classes
    class_<MST::Chain>("Chain",
        init<>())
        .def("residueSize", &MST::Chain::residueSize)
        .def("atomSize", &MST::Chain::atomSize)
        .def("getID", &MST::Chain::getID)
        .def("setID", &MST::Chain::setID)
        .def("getSegID", &MST::Chain::getSegID)
        .def("setSegID", &MST::Chain::setSegID)
    ;

    class_<MST::Structure>("Structure",
        init<std::string>())
        .def("chainSize", &MST::Structure::chainSize)
        .def("residueSize", &MST::Structure::residueSize)
        .def("atomSize", &MST::Structure::atomSize)
        .def("getChain", &MST::Structure::getChain, return_value_policy<reference_existing_object>())
        // the static_cast is needed to disambiguate an overloaded function
        .def("writePDB", static_cast<void (MST::Structure::*) (std::string, std::string)>(&MST::Structure::writePDB))
        .add_property("name", &MST::Structure::getName, &MST::Structure::setName)
    ;
}
