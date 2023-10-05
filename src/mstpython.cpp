/*
 * This inclusion should be put at the beginning.  It will include <Python.h>.
 */
#include <boost/python.hpp>
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <iostream>
#include <list>
#include <boost/python/stl_iterator.hpp>
#include "mstrotlib.h"
#include "msttypes.h"
#include "mstfasst.h"
#include "mstcondeg.h"
#include "mstsequence.h"
#include "mstfuser.h"
#include "mstrotlib.h"

// Source: https://stackoverflow.com/questions/15842126/feeding-a-python-list-into-a-function-taking-in-a-vector-with-boost-python/15940413#15940413

/// @brief Type that allows for registration of conversions from
///                python iterable types.

//can implement function here - mstreal of
//win body, return the more complex one

struct iterable_converter
{
    /// @note Registers converter from a python interable type to the
    ///             provided type.
    template <typename Container>
    iterable_converter&
    from_python()
    {
        boost::python::converter::registry::push_back(
                                                      &iterable_converter::convertible,
                                                      &iterable_converter::construct<Container>,
                                                      boost::python::type_id<Container>());

        // Support chaining.
        return *this;
    }

    /// @brief Check if PyObject is iterable.
    static void* convertible(PyObject* object)
    {
        return PyObject_GetIter(object) ? object : NULL;
    }

    /// @brief Convert iterable PyObject to C++ container type.
    ///
    /// Container Concept requirements:
    ///
    ///     * Container::value_type is CopyConstructable.
    ///     * Container can be constructed and populated with two iterators.
    ///         I.e. Container(begin, end)
    template <typename Container>
    static void construct(
                          PyObject* object,
                          boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        namespace python = boost::python;
        // Object is a borrowed reference, so create a handle indicting it is
        // borrowed for proper reference counting.
        python::handle<> handle(python::borrowed(object));

        // Obtain a handle to the memory block that the converter has allocated
        // for the C++ type.
        typedef python::converter::rvalue_from_python_storage<Container>
        storage_type;
        void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;

        typedef python::stl_input_iterator<typename Container::value_type>
        iterator;

        // Allocate the C++ type into the converter's memory block, and assign
        // its handle to the converter's convertible variable.    The C++
        // container is populated by passing the begin and end iterators of
        // the python object to the container's constructor.
        new (storage) Container(
                                iterator(python::object(handle)), // begin
                                iterator());                                            // end
        data->convertible = storage;
    }
};
//
//bool operator<(const fasstSolution& si, const fasstSolution& sj);
//
class Py_Structure: public MST::Structure {
public:
    Py_Structure(string pdbFile): MST::Structure(pdbFile) {};
    string structureToString(const Structure& structure) {
        std::stringstream ostream;
        structure.writePDB(ostream);
        return ostream.str();
    }
};



// Overloads for functions with optional arguments
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(addTargetsOverloads, FASST::addTargets, 1, 2)

/*
 * This is a macro Boost.Python provides to signify a Python extension module.
 */
BOOST_PYTHON_MODULE(mstpython) {
    // An established convention for using boost.python.
    using namespace boost::python;

    iterable_converter()
    .from_python<vector<double>>()
    .from_python<vector<int>>()
    .from_python<vector<vector<int>>>()
    .from_python<vector<Atom*>>()
    .from_python<vector<vector<Atom*>>>()

    .from_python<vector<Residue*>>()
    .from_python<vector<vector<Residue*>>>()
    .from_python<vector<std::string>>()
    .from_python<vector<MST::Sequence>>()
    ;
    // Various translations of vectors to lists needed by the classes below
    class_<vector<Atom*>>("AtomList")
    .def(vector_indexing_suite<vector<Atom*>>());

    class_<vector<vector<Atom*>>>("AtomListList")
    .def(vector_indexing_suite<vector<vector<Atom*>>>());

    class_<vector<Residue*>>("ResidueList")
    .def(vector_indexing_suite<vector<Residue*>>());

    class_<vector<int>>("IntList")
    .def(vector_indexing_suite<vector<int>>());

    class_<vector<mstreal>>("PointList")
    .def(vector_indexing_suite<vector<double>>());

    class_<vector<std::string>>("StringList")
    .def(vector_indexing_suite<vector<std::string>>());

    class_<vector<Sequence>>("SequenceList")
    .def(vector_indexing_suite<vector<MST::Sequence>>());

    class_<AtomPointerVector>("AtomPointerVector", init<const vector<Atom *> &>())
    .def("__getitem__", +[](const AtomPointerVector &a, int i) { return a[i]; }, return_value_policy<reference_existing_object>())
    .def("__len__", &AtomPointerVector::size);

    // expose classes

    class_<MST::Atom>("Atom",
                      init<>())
    .add_property("x", &MST::Atom::getX)
    .add_property("y", &MST::Atom::getY)
    .add_property("z", &MST::Atom::getZ)
    .def("getCoor", &MST::Atom::getCoor)
    .add_property("name", &MST::Atom::getName)
    .def("getAtomParent", &MST::Atom::getParent, return_value_policy<reference_existing_object>())
    ;

    class_<MST::Residue>("Residue",
                         init<>())
    .def("atomSize", &MST::Residue::atomSize)
    .add_property("name", &MST::Residue::getName)
    .def("getChainID", &MST::Residue::getChainID)
    .add_property("num", &MST::Residue::getNum, &MST::Residue::setNum)
    .add_property("iCode", &MST::Residue::getIcode, &MST::Residue::setIcode)
    .def("getAtom", &MST::Residue::getAtom, return_value_policy<reference_existing_object>())
    .def("getAtoms", &MST::Residue::getAtoms)
    .def("findAtom", &MST::Residue::findAtom, return_value_policy<reference_existing_object>())
    .def("previousResidue", &MST::Residue::previousResidue, return_value_policy<reference_existing_object>())
    .def("nextResidue", &MST::Residue::nextResidue, return_value_policy<reference_existing_object>())
    .add_property("phi", &MST::Residue::getPhi)
    .add_property("psi", &MST::Residue::getPsi)
    .add_property("omega", &MST::Residue::getOmega)
    .def("isBadDihedral", &MST::Residue::isBadDihedral)
    .staticmethod("isBadDihedral")
    .def("areBonded", static_cast<bool (*) (const MST::Residue&, const MST::Residue&, MST::mstreal)>(&MST::Residue::areBonded))
    .staticmethod("areBonded")
    .def("getResidueIndex", &MST::Residue::getResidueIndex)
    .def("getResidueIndexInChain", &MST::Residue::getResidueIndexInChain)
    .def("getStructure", &MST::Residue::getStructure, return_value_policy<reference_existing_object>())
    .def("getParent", &MST::Residue::getParent, return_value_policy<reference_existing_object>())
    ;

    def("emptyChain", +[]() { return new MST::Chain(); }, return_value_policy<reference_existing_object>());

    class_<MST::Chain>("Chain",
                       init<>())
    .def("residueSize", &MST::Chain::residueSize)
    .def("atomSize", &MST::Chain::atomSize)
    .add_property("id", &MST::Chain::getID, &MST::Chain::setID)
    .add_property("getSegID", &MST::Chain::getSegID, &MST::Chain::setSegID)
    .def("getResidue", &MST::Chain::getResidue, return_value_policy<reference_existing_object>())
    .def("getResidues", &MST::Chain::getResidues)
    .def("getStructure", &MST::Chain::getStructure, return_value_policy<reference_existing_object>())
    .def("getStructure", &MST::Chain::getStructure, return_value_policy<reference_existing_object>())
    .def("appendResidue", &MST::Chain::appendResidue)
    .def("appendResidueCopies", &MST::Chain::appendResidueCopies)
    ;

    class_<MST::Sequence>("Sequence", init<>())
    .def(init<const Chain>())
    .add_property("name", &MST::Sequence::getName, &MST::Sequence::setName)
    .def("toString", &MST::Sequence::toString)
    .def("getResidue", &MST::Sequence::getResidue)
    .def("__len__", &MST::Sequence::length)
    .def("appendResidue", static_cast<void (MST::Sequence::*) (const string&)>(&MST::Sequence::appendResidue))
    ;

    class_<MST::SeqTools>("SeqTools")
    .def("tripleToSingle", &MST::SeqTools::tripleToSingle)
    .staticmethod("tripleToSingle")
    .def("singleToTriple", &MST::SeqTools::singleToTriple)
    .staticmethod("singleToTriple")
    .def("seqToIdx", &MST::SeqTools::seqToIdx)
    .staticmethod("seqToIdx")
    .def("aaToIdx", &MST::SeqTools::aaToIdx)
    .staticmethod("aaToIdx")
    .def("isUnknown", &MST::SeqTools::isUnknown)
    .staticmethod("isUnknown")
    .def("idxToTriple", &MST::SeqTools::idxToTriple)
    .staticmethod("idxToTriple")
    .def("idxToSingle", &MST::SeqTools::idxToSingle)
    .staticmethod("idxToSingle")
    .def("readFasta", static_cast<vector<Sequence> (*) (const string&)>(&MST::SeqTools::readFasta))
    .staticmethod("readFasta")
    ;


    class_<RotamerLibrary>("RotamerLibrary", init<string>())
    .def("getBackbone", static_cast<vector<Atom *> (*) (const MST::Residue&, bool)>(&RotamerLibrary::getBackbone))
    .def("getBackbone", static_cast<vector<Atom *> (*) (const MST::Structure&, bool)>(&RotamerLibrary::getBackbone))
    ;

    def("emptyStructure", +[]() { return new MST::Structure(); }, return_value_policy<reference_existing_object>());

    class_<Structure>("Structure",
                      init<string,string>())
    .def(init<const vector<Atom *> &>())
    .def(init<const vector<Residue *> &>())
    .def("chainSize", &MST::Structure::chainSize)
    .def("residueSize", &MST::Structure::residueSize)
    .def("atomSize", &MST::Structure::atomSize)
    .def("getChain", &MST::Structure::getChain, return_value_policy<reference_existing_object>())
    .def("getResidue", &MST::Structure::getResidue,
         return_value_policy<reference_existing_object>())
    .def("getAtoms", &MST::Structure::getAtoms)
    .def("getResidues", &MST::Structure::getResidues)
    .def("appendChain", +[](MST::Structure& structure, MST::Chain *chain) {
        structure.appendChain(new MST::Chain(*chain));
    })/*static_cast<bool (MST::Structure::*) (Chain*, bool)>(&MST::Structure::appendChain))*/
    .def("deleteChain", &MST::Structure::deleteChain)
    .def("getChainByID", &MST::Structure::getChainByID, return_value_policy<reference_existing_object>())
    .def("addAtom", static_cast<void (MST::Structure::*) (MST::Atom*)>(&MST::Structure::addAtom))
    .def("addAtoms", static_cast<void (MST::Structure::*) (std::vector<MST::Atom*>)>(&MST::Structure::addAtoms))
    .def("addResidue", &MST::Structure::addResidue, return_value_policy<manage_new_object>())
    .def("getResidueIndex", &MST::Structure::getResidueIndex)
    .def("__eq__", &MST::Structure::operator==)
    .def("__ne__", &MST::Structure::operator!=)
    // the static_cast is needed to disambiguate an overloaded function
    .def("writePDB", static_cast<void (MST::Structure::*) (const std::string&, std::string) const>(&MST::Structure::writePDB))
    .def("__str__", &Py_Structure::structureToString)
    .add_property("name", &MST::Structure::getName, &MST::Structure::setName)
    .def("reassignChainsByConnectivity", static_cast<MST::Structure (MST::Structure::*) (MST::mstreal)> (&MST::Structure::reassignChainsByConnectivity))
    ;

    class_<RMSDCalculator>("RMSDCalculator", init<>())
    .def("rmsd", static_cast<mstreal (*) (const Structure&, const Structure&)>(&RMSDCalculator::rmsd))
    .staticmethod("rmsd")
    .def("rmsdCutoff", static_cast<mstreal (*) (const Structure&, mstreal, mstreal)>(&RMSDCalculator::rmsdCutoff))
    .def("rmsdCutoff", static_cast<mstreal (*) (const vector<int>&, mstreal, mstreal)>(&RMSDCalculator::rmsdCutoff))
    .def("apvRMSD", +[](const AtomPointerVector &a1, const AtomPointerVector &a2) { return RMSDCalculator::rmsd(a1, a2); })
    .staticmethod("apvRMSD")
    .def("apvBestRMSD", +[](RMSDCalculator &calc, const AtomPointerVector &a1, const AtomPointerVector &a2) { return calc.bestRMSD(a1, a2); })
    .def("allToOneBestRMSD", static_cast<std::vector<mstreal> (RMSDCalculator::*) (const vector<vector<Atom*>>&, const vector<Atom*>&, int, int)>(&RMSDCalculator::bestRMSD))
    ;

//    class_<fasstSolution, boost::noncopyable>("fasstSolution",init<>())
    class_<fasstSolution>("fasstSolution",init<>())
    .add_property("rmsd", &fasstSolution::getRMSD, &fasstSolution::setRMSD)
    .add_property("targetIndex", &fasstSolution::getTargetIndex)
    .add_property("alignment", &fasstSolution::getAlignment)
    .def("getSegLengths", &fasstSolution::getSegLengths)
    .def("numSegments", &fasstSolution::numSegments)
    .def("__getitem__", &fasstSolution::operator[])
    .def(self < self)
    //  the line below was giving me a compilation error, and the boost docs suggest that what I added above does the same thing
    //  https://www.boost.org/doc/libs/1_70_0/libs/python/doc/html/tutorial/tutorial/exposing.html#tutorial.exposing.class_operators_special_function
    //    .def("__lt__", static_cast<bool (*) (const fasstSolution&, const fasstSolution&)> (&operator<))
    ;

//    class_<fasstSolutionSet, boost::noncopyable>("fasstSolutionSet", init<>())
    class_<fasstSolutionSet>("fasstSolutionSet", init<>())
    .def("__iter__", boost::python::range(&fasstSolutionSet::begin, &fasstSolutionSet::end))
    .def("__getitem__", &fasstSolutionSet::operator[], return_value_policy<reference_existing_object>())
    .def("__len__", &fasstSolutionSet::size)
    .def("clear", &fasstSolutionSet::clear)
    .def("worstRMSD", &fasstSolutionSet::worstRMSD)
    .def("bestRMSD", &fasstSolutionSet::bestRMSD)
    ;

    class_<fasstSeqConstSimple>("fasstSeqConstSimple", init<int>())
    .def("addConstraint", &fasstSeqConstSimple::addConstraint)
    ;

    class_<fasstSearchOptions>("fasstSearchOptions", init<>())
    .add_property("minNumMatches", &fasstSearchOptions::getMinNumMatches, &fasstSearchOptions::setMinNumMatches)
    .add_property("maxNumMatches", &fasstSearchOptions::getMaxNumMatches, &fasstSearchOptions::setMaxNumMatches)
    .add_property("sufficientNumMatches", &fasstSearchOptions::getSufficientNumMatches, &fasstSearchOptions::setSufficientNumMatches)
    .add_property("rmsdCutoff", &fasstSearchOptions::getRMSDCutoff, &fasstSearchOptions::setRMSDCutoff)
    .add_property("minGap", &fasstSearchOptions::getMinGap, &fasstSearchOptions::setMinGap)
    .add_property("maxGap", &fasstSearchOptions::getMaxGap, &fasstSearchOptions::setMaxGap)
    .add_property("contextLength", &fasstSearchOptions::getContextLength, &fasstSearchOptions::setContextLength)
    .add_property("redundancyCut", &fasstSearchOptions::getRedundancyCut, &fasstSearchOptions::setRedundancyCut)
    .add_property("redundancyProperty", &fasstSearchOptions::getRedundancyProperty, &fasstSearchOptions::setRedundancyProperty)
    .def("unsetMinNumMatches", &fasstSearchOptions::unsetMinNumMatches)
    .def("unsetMaxNumMatches", &fasstSearchOptions::unsetMaxNumMatches)
    .def("unsetSufficientNumMatches", &fasstSearchOptions::unsetSufficientNumMatches)
    .def("resetGapConstraints", &fasstSearchOptions::resetGapConstraints)
    .def("unsetRedundancyCut", &fasstSearchOptions::unsetRedundancyCut)
    .def("unsetRedundancyProperty", &fasstSearchOptions::unsetRedundancyProperty)
    .def("isMinNumMatchesSet", &fasstSearchOptions::isMinNumMatchesSet)
    .def("isMaxNumMatchesSet", &fasstSearchOptions::isMaxNumMatchesSet)
    .def("isSufficientNumMatchesSet", &fasstSearchOptions::isSufficientNumMatchesSet)
    .def("isRedundancyCutSet", &fasstSearchOptions::isRedundancyCutSet)
    .def("isRedundancyPropertySet", &fasstSearchOptions::isRedundancyPropertySet)
    .def("minGapConstrained", &fasstSearchOptions::minGapConstrained)
    .def("maxGapConstrained", &fasstSearchOptions::maxGapConstrained)
    .def("gapConstrained", &fasstSearchOptions::gapConstrained)
    .def("gapConstraintsExist", &fasstSearchOptions::gapConstraintsExist)
    .def("sequenceConstraintsSet", &fasstSearchOptions::sequenceConstraintsSet)
    .def("isVerbose", &fasstSearchOptions::isVerbose)
    .def("validateGapConstraints", &fasstSearchOptions::validateGapConstraints)
    .def("validateSearchRequest", &fasstSearchOptions::validateSearchRequest)
    .def("areNumMatchConstraintsConsistent", &fasstSearchOptions::areNumMatchConstraintsConsistent)
    .def("setChainsDiff", &fasstSearchOptions::setChainsDiff)
    .def("resetDiffChainConstraints", &fasstSearchOptions::resetDiffChainConstraints)
    .def("setSequenceConstraints", +[](fasstSearchOptions &opts, const fasstSeqConstSimple &seqCons) { return opts.setSequenceConstraints(seqCons); })
    ;

    class_<contactList>("ContactList", init<>())
    .def("addContact", &contactList::addContact)
    .def("__len__", &contactList::size)
    .def("residueA", &contactList::residueA, return_value_policy<reference_existing_object>())
    .def("residueB", &contactList::residueB, return_value_policy<reference_existing_object>())
    .def("degree", static_cast<mstreal (contactList::*) (int)>(&contactList::degree))
    .def("sortByDegree", &contactList::sortByDegree)
    .def("areInContact", &contactList::areInContact)
    ;

    class_<ConFind, boost::noncopyable>("ConFind", init<string, Structure&>())
    .def("cache", static_cast<void (ConFind::*) (const Structure&)>(&ConFind::cache))
    .def("getNeighbors", static_cast<std::vector<Residue *> (ConFind::*) (Residue *)>(&ConFind::getNeighbors))
    .def("contactDegree", &ConFind::contactDegree)
    .def("getContacts", static_cast<contactList (ConFind::*) (Structure&, mstreal, contactList *)>(&ConFind::getContacts))
    .def("getResidueContacts", static_cast<contactList (ConFind::*) (Residue *, mstreal, contactList *)>(&ConFind::getContacts))
    ;

    class_<FASST>("FASST", init<>())
    .add_property("query", &FASST::getQuery)
    .def("setRMSDCutoff", &fasstSearchOptions::setRMSDCutoff)
    .def("setRedundancyCut", &fasstSearchOptions::setRedundancyCut)
    .def("setQuery", static_cast<void (FASST::*) (const Structure&, bool)>(&FASST::setQuery))
    .add_property("numQuerySegments", &FASST::getNumQuerySegments)
    .def("addTargetStructure", static_cast<void (FASST::*) (const Structure&, short)>(&FASST::addTarget))
    .def("addTarget", static_cast<void (FASST::*) (const string&, short)>(&FASST::addTarget))
    .def("addTargets", &FASST::addTargets, addTargetsOverloads())
    .add_property("options", make_function(&FASST::options, return_value_policy<reference_existing_object>()), &FASST::setOptions)
    .add_property("numTargets", &FASST::numTargets)
    .def("search", &FASST::search)
    .add_property("numMatches", &FASST::numMatches)
    .def("getMatches", &FASST::getMatches)
    .def("getTargetCopy",&FASST::getTargetCopy)
    .def("getMatchStructure", static_cast<Structure (FASST::*) (const fasstSolution&, bool, FASST::matchType, bool)>(&FASST::getMatchStructure))
    .def("getMatchResidueIndices", &FASST::getMatchResidueIndices)
    .def("getMatchSequence", &FASST::getMatchSequence)
    .def("getMatchSequences", &FASST::getMatchSequences)
    .def("readDatabase", &FASST::readDatabase)
    ;

    boost::python::enum_<FASST::matchType>("matchType")
    .value("REGION", FASST::REGION)
    .value("FULL", FASST::FULL)
    .value("WITHGAPS", FASST::WITHGAPS)
    .export_values()
    ;

    class_<fusionTopology>("fusionTopology", init<int>()) // init<const vector<vector<Residue*>> &>())
    .def("addFragment", static_cast<void (fusionTopology::*) (vector<Residue*>&, const vector<int>&, mstreal)>(&fusionTopology::addFragment))
    .def("addFixedPositions", &fusionTopology::addFixedPositions)
    .def("addFixedPosition", &fusionTopology::addFixedPosition)
    .def("numAlignedFrags", &fusionTopology::numAlignedFrags)
    .def("numFrags", &fusionTopology::numFrags)
    .def("numUniqueOverlaps", &fusionTopology::numUniqueOverlaps)
    .def("numFragsOverlapping", &fusionTopology::numFragsOverlapping)
    .def("getFragOverlapping", &fusionTopology::getFragOverlapping)
    .def("numOverlappingResidues", &fusionTopology::numOverlappingResidues)
    .def("getFixedPositions", &fusionTopology::getFixedPositions)
    .def("numFixedPositions", &fusionTopology::numFixedPositions)
    .def("__len__", &fusionTopology::length)
    .def("numMobileAtoms", static_cast<int (fusionTopology::*) ()>(&fusionTopology::numMobileAtoms))
    .def("numChains", &fusionTopology::numChains)
    .def("getChainLengths", &fusionTopology::getChainLengths)
    ;

    class_<fusionParams>("FusionParams", init<>());

    class_<fusionOutput>("FusionOutput", no_init)
    .add_property("bondScore", &fusionOutput::getBondScore)
    .add_property("angleScore", &fusionOutput::getAngleScore)
    .add_property("dihedralScore", &fusionOutput::getDihedralScore)
    .add_property("rmsdScore", &fusionOutput::getRMSDScore)
    .add_property("totRMSDScore", &fusionOutput::getTotRMSDScore)
    .add_property("score", &fusionOutput::getScore)
    ;

    class_<Fuser>("Fuser")
    // .def("fuse", static_cast<Structure (*) (const fusionTopology&, const fusionParams&)>(&Fuser::fuse))
    .def("fuse", +[](const fusionTopology &topo) {
        fusionOutput output;
        fusionParams params;
        Structure result = Fuser::fuse(topo, output, params);
        return boost::python::make_tuple(result, output);
    })
    .staticmethod("fuse")
    ;

    class_<MST::CartesianPoint>("CartesianPoint",init<Atom&>())
    .def(self + self)
    .def(self - self)
    .def(self * mstreal())
    .def(self / mstreal())
    .def(self += self)
    .def(self -= self)

    .def("norm", &CartesianPoint::norm)
    .def("norm2", &CartesianPoint::norm2)
    .def("mean", &CartesianPoint::mean)
    .def("cross", &CartesianPoint::cross)
    .def("dot", &CartesianPoint::dot)
    .def("getUnit", &CartesianPoint::getUnit)

    .def("getX", &CartesianPoint::getX)
    .def("getY", &CartesianPoint::getY)
    .def("getZ", &CartesianPoint::getZ)

    .def("distance", static_cast<mstreal (CartesianPoint::*) (const CartesianPoint&) const> (&CartesianPoint::distance))
    .def("distance2", static_cast<mstreal (CartesianPoint::*) (const CartesianPoint&) const> (&CartesianPoint::distance2))
    ;

    class_<MST::ProximitySearch>("ProximitySearch",
    init<const AtomPointerVector, mstreal>())
    .def("pointsWithin", &ProximitySearch::pointsWithin)
    .def("getPointsWithin", &ProximitySearch::getPointsWithin)
    ;
}
