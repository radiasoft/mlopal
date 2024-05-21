//
// Class OpalBeamline
//   :FIXME: add class description
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Elements/OpalBeamline.h"

#include "AbsBeamline/Bend2D.h"
#include "AbstractObjects/OpalData.h"
#include "Physics/Units.h"
#include "Structure/MeshGenerator.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <fstream>

OpalBeamline::OpalBeamline():
    elements_m(),
    prepared_m(false),
    containsSource_m(false)
{
}

OpalBeamline::OpalBeamline(const Vector_t& origin,
                           const Quaternion& rotation):
    elements_m(),
    prepared_m(false),
    containsSource_m(false),
    coordTransformationTo_m(origin, rotation)
{
}

OpalBeamline::~OpalBeamline() {
    elements_m.clear();
}

std::set<std::shared_ptr<Component>> OpalBeamline::getElements(const Vector_t &x) {
    std::set<std::shared_ptr<Component> > elementSet;
    FieldList::iterator it = elements_m.begin();
    const FieldList::iterator end = elements_m.end();
    for (; it != end; ++ it) {
        std::shared_ptr<Component> element = (*it).getElement();
        Vector_t r = element->getCSTrafoGlobal2Local().transformTo(x);

        if (element->isInside(r)) {
            elementSet.insert(element);
        }
    }

    return elementSet;
}

unsigned long OpalBeamline::getFieldAt(const unsigned int &/*index*/, const Vector_t &/*pos*/, const long &/*sindex*/, const double &/*t*/, Vector_t &/*E*/, Vector_t &/*B*/) {

    unsigned long rtv = 0x00;

    return rtv;
}

unsigned long OpalBeamline::getFieldAt(const Vector_t &position,
                                       const Vector_t &momentum,
                                       const double &t,
                                       Vector_t &Ef,
                                       Vector_t &Bf) {
    unsigned long rtv = 0x00;

    std::set<std::shared_ptr<Component>> elements = getElements(position);

    std::set<std::shared_ptr<Component>>::const_iterator it = elements.begin();
    const std::set<std::shared_ptr<Component>>::const_iterator end = elements.end();

    for (; it != end; ++ it) {
        ElementType type = (*it)->getType();
        if (type == ElementType::MONITOR ||
            type == ElementType::MARKER ||
            type == ElementType::CCOLLIMATOR) continue;

        Vector_t localR = transformToLocalCS(*it, position);
        Vector_t localP = rotateToLocalCS(*it, momentum);
        Vector_t localE(0.0), localB(0.0);

        (*it)->applyToReferenceParticle(localR, localP, t, localE, localB);

        Ef += rotateFromLocalCS(*it, localE);
        Bf += rotateFromLocalCS(*it, localB);
    }

    //         if(section.hasWake()) {
    //             rtv |= BEAMLINE_WAKE;
    //         }
    //         if(section.hasParticleMatterInteraction()) {
    //             rtv |= BEAMLINE_PARTICLEMATTERINTERACTION;
    //         }

    return rtv;
}

void OpalBeamline::switchElements(const double &min, const double &max, const double &kineticEnergy, const bool &nomonitors) {

    FieldList::iterator fprev;
    for(FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit) {
        // don't set online monitors if the centroid of the bunch is allready inside monitor
        // or if explicitly not desired (eg during auto phasing)
        if(flit->getElement()->getType() == ElementType::MONITOR) {
            double spos = (max + min) / 2.;
            if(!nomonitors && spos < (*flit).getStart()) {
                if(!(*flit).isOn() && max > (*flit).getStart()) {
                    (*flit).setOn(kineticEnergy);
                }
            }

        } else {
            if(!(*flit).isOn() && max > (*flit).getStart() && min < (*flit).getEnd()) {
                (*flit).setOn(kineticEnergy);
            }

            //check if multiple degraders follow one another with no other elements in between
            //if element is off and it is a degrader
            if (!(*flit).isOn() && flit->getElement()->getType() == ElementType::DEGRADER) {
                //check if previous element: is on, is a degrader, ends where new element starts
                if ((*fprev).isOn() && fprev->getElement()->getType() == ElementType::DEGRADER &&
                    ((*fprev).getEnd() + 0.01 > (*flit).getStart()) ) {
                    (*flit).setOn(kineticEnergy);
                }
            }
        }

        fprev = flit;
    }
}

void OpalBeamline::switchElementsOff() {
    for(FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit)
        (*flit).setOff();
}

void OpalBeamline::prepareSections() {
    if (elements_m.empty()) {
        prepared_m = true;
        return;
    }
    elements_m.sort(ClassicField::SortAsc);
    prepared_m = true;
}

void OpalBeamline::print(Inform &/*msg*/) const {
}

void OpalBeamline::swap(OpalBeamline & rhs) {
    std::swap(elements_m, rhs.elements_m);
    std::swap(prepared_m, rhs.prepared_m);
    std::swap(coordTransformationTo_m, rhs.coordTransformationTo_m);
}

void OpalBeamline::merge(OpalBeamline &rhs) {
    elements_m.insert(elements_m.end(),
                      rhs.elements_m.begin(),
                      rhs.elements_m.end());
    prepared_m = false;
    containsSource_m = containsSource_m || rhs.containsSource_m;
}


FieldList OpalBeamline::getElementByType(ElementType type) {
    if (type == ElementType::ANY) {
        return elements_m;
    }

    FieldList elements_of_requested_type;
    for(FieldList::iterator fit = elements_m.begin(); fit != elements_m.end(); ++ fit) {
        if((*fit).getElement()->getType() == type) {
            elements_of_requested_type.push_back((*fit));
        }
    }
    return elements_of_requested_type;
}

void OpalBeamline::positionElementRelative(std::shared_ptr<ElementBase> element) {
    if (!element->isPositioned()) {
        return;
    }

    element->releasePosition();
    CoordinateSystemTrafo toElement = element->getCSTrafoGlobal2Local();
    toElement *= coordTransformationTo_m;

    element->setCSTrafoGlobal2Local(toElement);
    element->fixPosition();
}

void OpalBeamline::compute3DLattice() {
    static unsigned int order = 0;
    const FieldList::iterator end = elements_m.end();

    unsigned int minOrder = order;
    {
        double endPriorPathLength = 0.0;
        CoordinateSystemTrafo currentCoordTrafo = coordTransformationTo_m;

        FieldList::iterator it = elements_m.begin();
        for (; it != end; ++ it) {
            std::shared_ptr<Component> element = (*it).getElement();
            if (element->isPositioned()) {
                continue;
            }
            (*it).order_m = minOrder;

            if (element->getType() != ElementType::SBEND &&
                element->getType() != ElementType::RBEND &&
                element->getType() != ElementType::RBEND3D) {
                continue;
            }

            double beginThisPathLength = element->getElementPosition();
            Vector_t beginThis3D(0, 0, beginThisPathLength - endPriorPathLength);
            BendBase * bendElement = static_cast<BendBase*>(element.get());
            double thisLength = bendElement->getChordLength();
            double bendAngle = bendElement->getBendAngle();
            double entranceAngle = bendElement->getEntranceAngle();
            double arcLength = (thisLength * std::abs(bendAngle) / (2 * sin(std::abs(bendAngle) / 2)));

            double rotationAngleAboutZ = bendElement->getRotationAboutZ();
            Quaternion_t rotationAboutZ(cos(0.5 * rotationAngleAboutZ),
                                        sin(-0.5 * rotationAngleAboutZ) * Vector_t(0, 0, 1));

            Vector_t effectiveRotationAxis = rotationAboutZ.rotate(Vector_t(0, -1, 0));
            effectiveRotationAxis /= euclidean_norm(effectiveRotationAxis);

            Quaternion_t rotationAboutAxis(cos(0.5 * bendAngle),
                                           sin(0.5 * bendAngle) * effectiveRotationAxis);
            Quaternion_t halfRotationAboutAxis(cos(0.25 * bendAngle),
                                               sin(0.25 * bendAngle) * effectiveRotationAxis);
            Quaternion_t entryFaceRotation(cos(0.5 * entranceAngle),
                                           sin(0.5 * entranceAngle) * effectiveRotationAxis);

            if (!Options::idealized) {
                std::vector<Vector_t> truePath = bendElement->getDesignPath();
                Quaternion_t directionExitHardEdge(cos(0.5 * (0.5 * bendAngle - entranceAngle)),
                                                   sin(0.5 * (0.5 * bendAngle - entranceAngle)) * effectiveRotationAxis);
                Vector_t exitHardEdge = thisLength * directionExitHardEdge.rotate(Vector_t(0, 0, 1));
                double distanceEntryHETruePath = euclidean_norm(truePath.front());
                double distanceExitHETruePath = euclidean_norm(rotationAboutZ.rotate(truePath.back()) - exitHardEdge);
                double pathLengthTruePath = (*it).getEnd() - (*it).getStart();
                arcLength = pathLengthTruePath - distanceEntryHETruePath - distanceExitHETruePath;
            }

            Vector_t chord = thisLength * halfRotationAboutAxis.rotate(Vector_t(0, 0, 1));
            Vector_t endThis3D = beginThis3D + chord;
            double endThisPathLength = beginThisPathLength + arcLength;

            CoordinateSystemTrafo fromEndLastToBeginThis(beginThis3D,
                                                         (entryFaceRotation * rotationAboutZ).conjugate());
            CoordinateSystemTrafo fromEndLastToEndThis(endThis3D,
                                                       rotationAboutAxis.conjugate());

            element->setCSTrafoGlobal2Local(fromEndLastToBeginThis * currentCoordTrafo);

            currentCoordTrafo = (fromEndLastToEndThis * currentCoordTrafo);

            endPriorPathLength = endThisPathLength;
        }
    }

    double endPriorPathLength = 0.0;
    CoordinateSystemTrafo currentCoordTrafo = coordTransformationTo_m;

    FieldList::iterator it = elements_m.begin();
    for (; it != end; ++ it) {
        std::shared_ptr<Component> element = (*it).getElement();
        if (element->isPositioned()) continue;

        (*it).order_m = order ++;

        double beginThisPathLength = element->getElementPosition();
        double thisLength = element->getElementLength();
        Vector_t beginThis3D(0, 0, beginThisPathLength - endPriorPathLength);

        if (element->getType() == ElementType::SOURCE) {
            beginThis3D(2) -= thisLength;
        }

        Vector_t endThis3D;
        if (element->getType() == ElementType::SBEND ||
            element->getType() == ElementType::RBEND ||
            element->getType() == ElementType::RBEND3D) {

            BendBase * bendElement = static_cast<BendBase*>(element.get());
            thisLength = bendElement->getChordLength();
            double bendAngle = bendElement->getBendAngle();

            double rotationAngleAboutZ = bendElement->getRotationAboutZ();
            Quaternion_t rotationAboutZ(cos(0.5 * rotationAngleAboutZ),
                                        sin(-0.5 * rotationAngleAboutZ) * Vector_t(0, 0, 1));

            Vector_t effectiveRotationAxis = rotationAboutZ.rotate(Vector_t(0, -1, 0));
            effectiveRotationAxis /= euclidean_norm(effectiveRotationAxis);

            Quaternion_t rotationAboutAxis(cos(0.5 * bendAngle),
                                           sin(0.5 * bendAngle) * effectiveRotationAxis);
            Quaternion halfRotationAboutAxis(cos(0.25 * bendAngle),
                                             sin(0.25 * bendAngle) * effectiveRotationAxis);

            double arcLength = (thisLength * std::abs(bendAngle) /
                                (2 * sin(bendAngle / 2)));
            if (!Options::idealized) {
                std::vector<Vector_t> truePath = bendElement->getDesignPath();
                double entranceAngle = bendElement->getEntranceAngle();
                Quaternion_t directionExitHardEdge(cos(0.5 * (0.5 * bendAngle - entranceAngle)),
                                                   sin(0.5 * (0.5 * bendAngle - entranceAngle)) * effectiveRotationAxis);
                Vector_t exitHardEdge = thisLength * directionExitHardEdge.rotate(Vector_t(0, 0, 1));
                double distanceEntryHETruePath = euclidean_norm(truePath.front());
                double distanceExitHETruePath = euclidean_norm(rotationAboutZ.rotate(truePath.back()) - exitHardEdge);
                double pathLengthTruePath = (*it).getEnd() - (*it).getStart();
                arcLength = pathLengthTruePath - distanceEntryHETruePath - distanceExitHETruePath;
            }

            endThis3D = (beginThis3D +
                         halfRotationAboutAxis.rotate(Vector_t(0, 0, thisLength)));
            CoordinateSystemTrafo fromEndLastToEndThis(endThis3D,
                                                       rotationAboutAxis.conjugate());
            currentCoordTrafo = fromEndLastToEndThis * currentCoordTrafo;

            endPriorPathLength = beginThisPathLength + arcLength;
        } else {
            double rotationAngleAboutZ = (*it).getElement()->getRotationAboutZ();
            Quaternion_t rotationAboutZ(cos(0.5 * rotationAngleAboutZ),
                                        sin(-0.5 * rotationAngleAboutZ) * Vector_t(0, 0, 1));

            CoordinateSystemTrafo fromLastToThis(beginThis3D, rotationAboutZ);

            element->setCSTrafoGlobal2Local(fromLastToThis * currentCoordTrafo);
        }

        element->fixPosition();
    }
}

void OpalBeamline::save3DLattice() {
    if (Ippl::myNode() != 0 || OpalData::getInstance()->isOptimizerRun()) return;

    elements_m.sort([](const ClassicField& a, const ClassicField& b) {
            return a.order_m < b.order_m;
        });

    FieldList::iterator it = elements_m.begin();
    FieldList::iterator end = elements_m.end();

    std::ofstream pos;
    std::string fileName = Util::combineFilePath({
        OpalData::getInstance()->getAuxiliaryOutputDirectory(),
        OpalData::getInstance()->getInputBasename() + "_ElementPositions.txt"
    });
    if (OpalData::getInstance()->getOpenMode() == OpalData::OpenMode::APPEND &&
        boost::filesystem::exists(fileName)) {
        pos.open(fileName, std::ios_base::app);
    } else {
        pos.open(fileName);
    }

    MeshGenerator mesh;
    for (; it != end; ++ it) {
        std::shared_ptr<Component> element = (*it).getElement();
        CoordinateSystemTrafo toBegin = element->getEdgeToBegin() * element->getCSTrafoGlobal2Local();
        CoordinateSystemTrafo toEnd = element->getEdgeToEnd() * element->getCSTrafoGlobal2Local();
        Vector_t entry3D = toBegin.getOrigin();
        Vector_t exit3D = toEnd.getOrigin();

        mesh.add(*(element.get()));

        if (element->getType() == ElementType::SBEND ||
            element->getType() == ElementType::RBEND) {

            Bend2D * bendElement = static_cast<Bend2D*>(element.get());
            std::vector<Vector_t> designPath = bendElement->getDesignPath();
            toEnd = bendElement->getBeginToEnd_local() * element->getCSTrafoGlobal2Local();
            exit3D = toEnd.getOrigin();

            unsigned int size = designPath.size();
            unsigned int minNumSteps = std::max(20.0,
                                                std::abs(bendElement->getBendAngle() * Units::rad2deg));
            unsigned int frequency = std::floor((double)size / minNumSteps);

            pos << std::setw(30) << std::left << std::string("\"ENTRY EDGE: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << entry3D(2)
                << std::setw(18) << std::setprecision(10) << entry3D(0)
                << std::setw(18) << std::setprecision(10) << entry3D(1)
                << "\n";

            Vector_t position = element->getCSTrafoGlobal2Local().transformFrom(designPath.front());
            pos << std::setw(30) << std::left << std::string("\"BEGIN: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << position(2)
                << std::setw(18) << std::setprecision(10) << position(0)
                << std::setw(18) << std::setprecision(10) << position(1)
                << std::endl;

            for (unsigned int i = frequency; i + 1 < size; i += frequency) {

                position = element->getCSTrafoGlobal2Local().transformFrom(designPath[i]);
                pos << std::setw(30) << std::left << std::string("\"MID: ") + element->getName() + std::string("\"")
                    << std::setw(18) << std::setprecision(10) << position(2)
                    << std::setw(18) << std::setprecision(10) << position(0)
                    << std::setw(18) << std::setprecision(10) << position(1)
                    << std::endl;
            }

            position = element->getCSTrafoGlobal2Local().transformFrom(designPath.back());
            pos << std::setw(30) << std::left << std::string("\"END: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << position(2)
                << std::setw(18) << std::setprecision(10) << position(0)
                << std::setw(18) << std::setprecision(10) << position(1)
                << std::endl;

            pos << std::setw(30) << std::left << std::string("\"EXIT EDGE: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << exit3D(2)
                << std::setw(18) << std::setprecision(10) << exit3D(0)
                << std::setw(18) << std::setprecision(10) << exit3D(1)
                << std::endl;
        } else {
            pos << std::setw(30) << std::left << std::string("\"BEGIN: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << entry3D(2)
                << std::setw(18) << std::setprecision(10) << entry3D(0)
                << std::setw(18) << std::setprecision(10) << entry3D(1)
                << "\n";

            pos << std::setw(30) << std::left << std::string("\"END: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << exit3D(2)
                << std::setw(18) << std::setprecision(10) << exit3D(0)
                << std::setw(18) << std::setprecision(10) << exit3D(1)
                << std::endl;
        }
    }
    elements_m.sort(ClassicField::SortAsc);
    mesh.write(OpalData::getInstance()->getInputBasename());
}

namespace {
    std::string parseInput() {

        std::ifstream in(OpalData::getInstance()->getInputFn());
        std::string source("");
        std::string str;
        char testBit;
        const std::string commentFormat("");
        const boost::regex empty("^[ \t]*$");
        const boost::regex lineEnd(";");
        const std::string lineEndFormat(";\n");
        const boost::regex cppCommentExpr("//.*");
        const boost::regex cCommentExpr("/\\*.*?\\*/"); // "/\\*(?>[^*/]+|\\*[^/]|/[^*])*(?>(?R)(?>[^*/]+|\\*[^/]|/[^*])*)*\\*/"
        bool priorEmpty = true;

        in.get(testBit);
        while (!in.eof()) {
            in.putback(testBit);

            std::getline(in, str);
            str = boost::regex_replace(str, cppCommentExpr, commentFormat);
            str = boost::regex_replace(str, empty, commentFormat);
            if (!str.empty()) {
                source += str;// + '\n';
                priorEmpty = false;
            } else if (!priorEmpty) {
                source += "##EMPTY_LINE##";
                priorEmpty = true;
            }

            in.get(testBit);
        }

        source = boost::regex_replace(source, cCommentExpr, commentFormat);
        source = boost::regex_replace(source, lineEnd, lineEndFormat, boost::match_default | boost::format_all);

        // Since the positions of the elements are absolute in the laboratory coordinate system we have to make
        // sure that the line command doesn't provide an origin and orientation. Everything after the sequence of
        // elements can be deleted and only "LINE = (...);", the first sub-expression (denoted by '\1'), should be kept.
        const boost::regex lineCommand("(LINE[ \t]*=[ \t]*\\([^\\)]*\\))[ \t]*,[^;]*;", boost::regex::icase);
        const std::string firstSubExpression("\\1;");
        source = boost::regex_replace(source, lineCommand, firstSubExpression);

        return source;
    }

    unsigned int getMinimalSignificantDigits(double num, const unsigned int maxDigits) {
        char buf[32];
        snprintf(buf, 32, "%.*f", maxDigits + 1, num);
        std::string numStr(buf);
        unsigned int length = numStr.length();

        unsigned int numDigits = maxDigits;
        unsigned int i = 2;
        while (i < maxDigits + 1 && numStr[length - i] == '0') {
            --numDigits;
            ++ i;
        }

        return numDigits;
    }

    std::string round2string(double num, const unsigned int maxDigits) {
        char buf[64];

        snprintf(buf, 64, "%.*f", getMinimalSignificantDigits(num, maxDigits), num);

        return std::string(buf);
    }
}

void OpalBeamline::save3DInput() {
    if (Ippl::myNode() != 0 || OpalData::getInstance()->isOptimizerRun()) return;

    FieldList::iterator it = elements_m.begin();
    FieldList::iterator end = elements_m.end();

    std::string input = parseInput();
    std::string fname = Util::combineFilePath({
        OpalData::getInstance()->getAuxiliaryOutputDirectory(),
        OpalData::getInstance()->getInputBasename() + "_3D.opal"
    });
    std::ofstream pos(fname);

    for (; it != end; ++ it) {
        std::shared_ptr<Component> element = (*it).getElement();
        std::string elementName = element->getName();
        const boost::regex replacePSI("(" + elementName + "\\s*:[^\\n]*)PSI\\s*=[^,;]*,?", boost::regex::icase);
        input = boost::regex_replace(input, replacePSI, "\\1\\2");

        const boost::regex replaceELEMEDGE("(" + elementName + "\\s*:[^\\n]*)ELEMEDGE\\s*=[^,;]*(.)", boost::regex::icase);

        CoordinateSystemTrafo cst = element->getCSTrafoGlobal2Local();
        Vector_t origin = cst.getOrigin();
        Vector_t orient = Util::getTaitBryantAngles(cst.getRotation().conjugate(), elementName);
        for (unsigned int d = 0; d < 3; ++ d)
            orient(d) *= Units::rad2deg;

        std::string x = (std::abs(origin(0)) > 1e-10? "X = " + round2string(origin(0), 10) + ", ": "");
        std::string y = (std::abs(origin(1)) > 1e-10? "Y = " + round2string(origin(1), 10) + ", ": "");
        std::string z = (std::abs(origin(2)) > 1e-10? "Z = " + round2string(origin(2), 10) + ", ": "");

        std::string theta = (orient(0) > 1e-10? "THETA = " + round2string(orient(0), 6) + " * PI / 180, ": "");
        std::string phi = (orient(1) > 1e-10? "PHI = " + round2string(orient(1), 6) + " * PI / 180, ": "");
        std::string psi = (orient(2) > 1e-10? "PSI = " + round2string(orient(2), 6) + " * PI / 180, ": "");
        std::string coordTrafo = x + y + z + theta + phi + psi;
        if (coordTrafo.length() > 2) {
            coordTrafo = coordTrafo.substr(0, coordTrafo.length() - 2); // remove last ', '
        }

        std::string position = ("\\1" + coordTrafo + "\\2");

        input = boost::regex_replace(input, replaceELEMEDGE, position);

        if (element->getType() == ElementType::RBEND ||
            element->getType() == ElementType::SBEND) {
            const Bend2D* dipole = static_cast<const Bend2D*>(element.get());
            double angle = dipole->getBendAngle();
            double E1 = dipole->getEntranceAngle();
            double E2 = dipole->getExitAngle();

            const boost::regex angleR("(" + elementName + "\\s*:[^\\n]*ANGLE\\s*=)[^,;]*(.)");
            const std::string angleF("\\1 " + round2string(angle * Units::rad2deg, 6) + " / 180 * PI\\2");
            const boost::regex E1R("(" + elementName + "\\s*:[^\\n]*E1\\s*=)[^,;]*(.)");
            const std::string E1F("\\1 " + round2string(E1 * Units::rad2deg, 6) + " / 180 * PI\\2");
            const boost::regex E2R("(" + elementName + "\\s*:[^\\n]*E2\\s*=)[^,;]*(.)");
            const std::string E2F("\\1 " + round2string(E2 * Units::rad2deg, 6) + " / 180 * PI\\2");
            const boost::regex noRotation("(" + elementName + "\\s*:[^\\n]*),\\s*ROTATION\\s*=[^,;]*(.)");
            const std::string noRotationFormat("\\1\\2  ");

            input = boost::regex_replace(input, angleR, angleF);
            input = boost::regex_replace(input, E1R, E1F);
            input = boost::regex_replace(input, E2R, E2F);
            input = boost::regex_replace(input, noRotation, noRotationFormat);
        }
    }

    const boost::regex empty("##EMPTY_LINE##");
    const std::string emptyFormat("\n");
    input = boost::regex_replace(input, empty, emptyFormat);

    pos << input << std::endl;
}

void OpalBeamline::activateElements() {
    auto it = elements_m.begin();
    const auto end = elements_m.end();

    double designEnergy = 0.0;
    for (; it != end; ++ it) {
        std::shared_ptr<Component> element = (*it).getElement();
        if (element->getType() == ElementType::SBEND ||
            element->getType() == ElementType::RBEND) {
            Bend2D * bendElement = static_cast<Bend2D*>(element.get());
            designEnergy = bendElement->getDesignEnergy() * Units::eV2MeV;
        }
        (*it).setOn(designEnergy);
        // element->goOnline(designEnergy);
    }
}