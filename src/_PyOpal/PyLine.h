#ifndef PyOpal_PyLine_h


#include "Lines/Sequence.h"
#include "Lines/Line.h"
#include "Beamlines/TBeamline.h"
#include "PyOpal/PyOpalObject.h"

namespace PyOpal {

// to get the inheritance from PyElement I need to make PyLine_ templated; but
// really I only ever want to associate a PyLine_ with a TBeamline (it is a list
// of elements in a beamline). So I typedef a PyLine class which is the
// specialisation to PyLine_<TBeamline>.
template <class C>
class PyLine_; // forwards declaration of the PyLine_
typedef PyLine_<TBeamline<FlaggedElmPtr> > PyLine;

/** PyLine_ is a list of elements in a beamline
 *  
 *  PyLine_ is designed to look like a python list; but internally it holds
 *  a TBeamline object (which is used by OPAL to do a line of elements). 
 *  TBeamline internally holds a linked list (std::list) of ElementBase. In
 *  order to keep the type information, which is valuable for users, we maintain
 *  in PyLine_ a std::vector of python objects of the appropriate Component 
 *  type (the child class). This is in addition to the ElementBases stored in
 *  TBeamline. PyLine_ has the job of keeping the TBeamline and 
 *  the std::vector up to date and consistent.
 * 
 *  Normally I expect PyLine_<TBeamline<FlaggedElmPtr> > to be instantiated,
 *  which is typedef'd to a PyLine. I have not targetted any other
 *  specialization - caveat emptor if you try to do anything else.
 */
template <class C>
class PyLine_ : public PyOpalObjectNS::PyOpalObject<C> {
public:
    /** Get a python object from the vector of stored objects 
     * 
     *  i - index of the element in the python object 
     * 
     *  Returns the stored python object.
     */
    boost::python::object getElement(int i);

    /** Set a python object in the line vector and update the TBeamline. 
     * 
     *  i - index of the element in the python object 
     *  element - the element to be stored. element must have a python method
     *            get_opal_element that returns python::object that is 
     *            convertible to a PyElement<OpalElement> using python::extract.
     *  Note that OpalElement::update is called at this time to update the 
     *  underlying element properties.
     */
    void setElement(int i, boost::python::object element);

    /** Extend the length of the vector and TBeamline. 
     *  
     *  element - the element to be appended. Must be suitable for use in 
     *            setElement (see setElement doc for details).
     */
    void append(boost::python::object element);

    /** Returns the number of elements in the line */
    int getLength() const {return line.size();}

    /** Make a python::class_ object for a PyLine.
     * 
     *  Should normally be called during module definition.
     */
    boost::python::class_<PyLine> make_class(const char* className);

    void registerObject();

private:
    /** The python objects stored in the line */
    std::vector<boost::python::object> line;
};

template <>
void PyLine_<TBeamline<FlaggedElmPtr> >::registerObject() {
    TBeamline<FlaggedElmPtr>* wrapped = getOpalShared().get();
    std::cerr << "RegisterObject Line 1 " << wrapped << std::endl;
    Line* line = new Line();
    line->setElement(wrapped);
    Object* objectPtr = dynamic_cast<Object*>(line);
    std::cerr << "RegisterObject Line 2 " << objectPtr << std::endl;
    if (objectPtr == NULL) {
        throw OpalException("PyLine_<TBeamline<FlaggedElmPtr> >::register",
                            "Trying to register something that was not a Opal Object");
    }
    OpalData::getInstance()->define(objectPtr);
}

template<>
boost::python::class_<PyLine> PyLine_<TBeamline<FlaggedElmPtr> >::make_class(const char* className) {
    boost::python::docstring_options docop(true, true, false); // user_def, py_sig, cpp_sig
    auto pyclass = boost::python::class_<PyLine>(className);
    return pyclass;
}

template<>
boost::python::object PyLine_<TBeamline<FlaggedElmPtr> >::getElement(int i) {
    try {
        return line.at(i);
    } catch (std::exception& exc) {
        throw OpalException("PyLine::getElement", "Out of range");
    }
}

template<>
void PyLine_<TBeamline<FlaggedElmPtr> >::setElement(int i, boost::python::object pyelement) {
    // TBeamline is implemented as a double linked list??
    std::cerr << "PYLINE setElement " << i << std::endl;
    typedef TBeamline<FlaggedElmPtr> BL;
    try {
        line.at(i) = pyelement;
    } catch (std::exception& exc) {
        throw OpalException("PyLine::setElement", "Failed to set element");
    }
    boost::python::object pyopalelement = pyelement.attr("get_opal_element")();
    PyOpal::PyOpalObjectNS::PyOpalObject<OpalElement>& cpyelement =
            boost::python::extract<PyOpal::PyOpalObjectNS::PyOpalObject<OpalElement>& >(pyopalelement);
    int index = 0;
    for (BL::iterator it = object_m->begin(); it != object_m->end(); ++it) {
        if (index == i) {
            // 4 layers of nested inheritance template wrapper darkness
            std::shared_ptr<OpalElement> opalElementShared = cpyelement.getOpalShared();
            OpalElement* opalElement = opalElementShared.get();
            if (!opalElement) {
                throw OpalException("PyLine::setElement", "Failed to extract element");
            } 
            opalElement->update();
            ElementBase* elmbase = opalElement->getElement()->removeWrappers();
            std::cerr << "    PYLINE setElement loop " << elmbase << std::endl;
            if (!elmbase) {
                throw OpalException("PyLine::setElement", "Failed to cast element");
            } 
            ElmPtr elmptr(elmbase);
            FlaggedElmPtr felmptr(elmptr);
            object_m->insert(it, felmptr);
            return;
        }
        ++index;
    }
    throw OpalException("PyLine::getElement", "Out of range");
}

template <>
void PyLine_<TBeamline<FlaggedElmPtr> >::append(boost::python::object pyelement) {
    int i = line.size();
    // first we extend the size of the python object list by 1 with dummy variables
    line.push_back(boost::python::object());
    FlaggedElmPtr felmptr(ElmPtr(NULL));
    object_m->push_back(felmptr);
    // then we set the python object using setElement method
    setElement(i, pyelement);
}

}
#endif // PyOpal_PyLine_h

/*
OPAL Beamline insanity:

There exists:
1. Beamline is a beamline
2. TBeamline which is another sort of beamline that inherits from Beamline (but has no reason to exist)
3. BeamSequence which does nothing
4. Line which has some bindings to the UI
5. FlaggedBeamline which has a flag to indicate the Beamline is reflected (WTF?)
6. Sequence which is a list of Sequence elements whatever they are

Probably some other stuff. It's insane! Should be two classes
1. Classic Beamline
2. UI bindings
Done.

*/

