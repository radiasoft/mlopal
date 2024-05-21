//
// Class OpalData
//   The global OPAL structure.
//   The OPAL object holds all global data required for a OPAL execution.
//   In particular it contains the main Directory, which allows retrieval
//   of command objects by their name.  For other data refer to the
//   implementation file.
//
// Copyright (c) 200x - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_OpalData_HH
#define OPAL_OpalData_HH

#include "AbstractObjects/ObjectFunction.h"

#include <iosfwd>
#include <string>
#include <vector>
#include <map>
#include <stack>

template <class T = double, unsigned Dim = 3>
class PartBunchBase;

class AttributeBase;
class Object;
class Table;
class ValueDefinition;
class DataSink;
class BoundaryGeometry;

// store element name, max phase
typedef std::pair<std::string, double > MaxPhasesT;
typedef std::map<double, double> energyEvolution_t;


/// The global OPAL structure.
class OpalData {

public:

    static OpalData *getInstance();

    static void deleteInstance();

    static void stashInstance();

    static OpalData *popInstance();

    ~OpalData();

    /// Enum for writing to files
    enum class OpenMode: unsigned short {
        UNDEFINED,
        WRITE,
        APPEND
    };

    /// reset object for consecutive runs
    void reset();

    /// Apply a function to all objects.
    //  Loop over the directory and apply the given functor object to each
    //  object in turn.
    void apply(const ObjectFunction &);

    /// Create new object.
    //  No replacement is allowed; if an object with the same name exists,
    //  throw [b]OpalException[/b].
    void create(Object *newObject);

    /// Define a new object.
    //  Replacement is allowed; however [b]OpalException[/b] is thrown,
    //  if the replacement cannot be done.
    void define(Object *newObject);

    /// Delete existing entry.
    //  Identified by [b]name[/b].
    void erase(const std::string &name);

    /// Find entry.
    //  Identified by [b]name[/b].
    Object *find(const std::string &name);

    /// Return value of global reference momentum.
    double getP0() const;

    /// Invalidate expressions.
    //  Force re-evaluation of all expressions before next command is
    //  executed.
    //  Also set the [b]modified[/b] flag in [b]object[/b], if not nullptr.
    void makeDirty(Object *object);

    /// Print all objects.
    //  Loop over the directory and print each object whose name matches
    //  the regular expression [b]pattern[/b].
    void printNames(std::ostream &stream, const std::string &pattern);

    /// Get a list of all objects.
    //  Loop over the directory and return a list of names of all objects
    std::vector<std::string> getAllNames() const;

    /// Register table.
    //  Register the table [b]t[/b].
    //  Registered tables are invalidated to be refilled when an object
    //  on which they depend is changed or replaced.
    void registerTable(Table *t);

    /// Unregister table.
    void unregisterTable(Table *t);

    /// Register expression.
    //  Registered expressions are invalidated to be recomputed when
    //  any object in the directory is changed or replaced.
    void registerExpression(AttributeBase *);

    /// Unregister expression.
    void unregisterExpression(AttributeBase *);

    /// Set the global momentum.
    void setP0(ValueDefinition *p0);

    /// Store the page title.
    void storeTitle(const std::string &);

    /// Print the page title.
    void printTitle(std::ostream &);

    /// Get the title string
    std::string getTitle();

    /// Update all objects.
    //  Loop over the directory and notify all objects to update themselves.
    void update();

    /// Clear Reference.
    //  This functor is used to clear the reference count stored in an object.
    struct ClearReference: public ObjectFunction {
        virtual void operator()(Object *) const;
    };

    std::map<std::string, std::string> getVariableData();
    std::vector<std::string> getVariableNames();

    bool isInOPALCyclMode();
    bool isInOPALTMode();
    bool isOptimizerRun();

    void setInOPALCyclMode();
    void setInOPALTMode();
    void setOptimizerFlag();

    bool isInPrepState();
    void setInPrepState(bool state);

    /// true if in follow-up track
    bool hasPriorTrack();

    /// true if in follow-up track
    void setPriorTrack(const bool &value = true);

    /// true if we do a restart run
    bool inRestartRun();

    /// set OPAL in restart mode
    void setRestartRun(const bool &value = true);

    /// store the location where to restart
    void setRestartStep(int s);

    /// get the step where to restart
    int getRestartStep();

    /// get the name of the the additional data directory
    std::string getAuxiliaryOutputDirectory() const;

    /// get opals input filename
    std::string getInputFn();

    /// get input file name without extension
    std::string getInputBasename();

    /// store opals input filename
    void storeInputFn(const std::string &fn);

    /// checks the output file names of all items to avoid duplicates
    void checkAndAddOutputFileName(const std::string &outfn);

    /// get opals restart h5 format filename
    std::string getRestartFileName();

    /// store opals restart h5 format filename
    void setRestartFileName(std::string s);

    /// true if we do a restart from specified h5 file
    bool hasRestartFile();

    /// set the dump frequency as found in restart file
    void setRestartDumpFreq(const int &N);

    /// get the dump frequency as found in restart file
    int getRestartDumpFreq() const;

    void setOpenMode(OpenMode openMode);
    OpenMode getOpenMode() const;

    /// set the last step in a run for possible follow-up run
    void setLastStep(const int &step);

    /// get the last step from a possible previous run
    int getLastStep() const;

    /// true if we already allocated a ParticleBunch object
    bool hasBunchAllocated();

    void bunchIsAllocated();

    PartBunchBase<double, 3> *getPartBunch();

    void setPartBunch(PartBunchBase<double, 3> *p);

    /// true if we already allocated a DataSink object
    bool hasDataSinkAllocated();

    DataSink *getDataSink();

    void setDataSink(DataSink *s);

    /// units: (sec)
    void   setGlobalPhaseShift(double shift);
    /// units: (sec)
    double getGlobalPhaseShift();

    ///
    void setGlobalGeometry(BoundaryGeometry *bg);
    ///
    BoundaryGeometry *getGlobalGeometry();

    bool hasGlobalGeometry();

    void setMaxPhase(std::string elName, double phi);

    std::vector<MaxPhasesT>::iterator getFirstMaxPhases();
    std::vector<MaxPhasesT>::iterator getLastMaxPhases();
    int getNumberOfMaxPhases();

    void addEnergyData(double spos, double ekin);
    energyEvolution_t::iterator getFirstEnergyData();
    energyEvolution_t::iterator getLastEnergyData();

    unsigned long long getMaxTrackSteps();
    void setMaxTrackSteps(unsigned long long s);
    void incMaxTrackSteps(unsigned long long s);

    void addProblemCharacteristicValue(const std::string &name, unsigned int value);
    const std::map<std::string, unsigned int> &getProblemCharacteristicValues() const;

    void storeArguments(int argc, char *argv[]);
    std::vector<std::string> getArguments();

private:

    static bool isInstantiated;
    static OpalData *instance;
    static std::stack<OpalData*> stashedInstances;

    OpalData();

    // Not implemented.
    OpalData(const OpalData &);
    void operator=(const OpalData &);

    // The private implementation details.
    struct OpalDataImpl *p;
};

//extern OpalData OPAL;

#endif // OPAL_OpalData_HH
