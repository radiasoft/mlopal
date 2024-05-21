#ifndef OPAL_ELEMENT_POSITION_WRITER_H
#define OPAL_ELEMENT_POSITION_WRITER_H

#include "SDDSWriter.h"

class ElementPositionWriter : public SDDSWriter {

public:
    explicit ElementPositionWriter(const std::string& fname);

    void addRow(double spos,
                const std::vector<double>& row,
                const std::string& elements);

private:
    void fillHeader();
};


#endif