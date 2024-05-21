#include "gtest/gtest.h"

#include "Index/NDIndex.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Field/BCond.h"
#include "Field/GuardCellSizes.h"
#include "AppTypes/Vektor.h"
#include "Utility/FieldDebug.h"

#include <fstream>

namespace {
    void hardCodedOutput(std::string filename); // Prototype of function defined below.
    bool thediff(std::string filename1, std::string filename2);
}

TEST(Field, FieldDebug)
{
    bool docomm = true; // Should try the test both ways, really....

    const unsigned Dim3=3;

    int nCells[Dim3];       // Number of grid cells in each direction
    unsigned nVNodes[Dim3]; // Number of vnodes (subdomains) in each direction.

    // Hardwired values for automated test, as in regression testing:
    for (unsigned int d=0; d<Dim3; d++) nCells[d] = 4;
    for (unsigned int d=0; d<Dim3; d++) nVNodes[d] = 2;

    int nx, ny, nz;
    nx = nCells[0]; ny = nCells[1]; nz = nCells[2];

    Index I(nx); Index J(ny); Index K(nz);
    // Specify multipple vnodes (8) to make sure this works right:
    //  FieldLayout<Dim3> layout3(I,J,K,PARALLEL,PARALLEL,PARALLEL,8);
    NDIndex<Dim3> ndi; ndi[0] = I; ndi[1] = J; ndi[2] = K;
    e_dim_tag serialParallelSpec[Dim3];
    for (unsigned int d=0; d<Dim3; d++) serialParallelSpec[d] = PARALLEL;
    FieldLayout<Dim3> layout3(ndi, serialParallelSpec, nVNodes);

    // New Inform-based version (tjw):
    Inform* fdip =
        new Inform(nullptr,"text.test.TestFieldDebug",Inform::OVERWRITE,0);
    Inform& fdi = *fdip;
    setInform(fdi);

    // Put guard cells and red-flag (value = -999) boundary conditions on
    // Fields, to make sure nothing funny is happening:
    GuardCellSizes<Dim3> gc(2);
    BConds<double,Dim3> sbc;
    for (unsigned int face=0; face < 2*Dim3; face++) {
        sbc[face] = new ConstantFace<double,Dim3>(face,-999.0);
    }
    BConds<Vektor<double,Dim3>,Dim3> vbc;
    for (unsigned int face=0; face < 2*Dim3; face++) {
        vbc[face] = new ConstantFace<Vektor<double,Dim3>,Dim3>(face,-999.0);
    }

    // Scalar Field -------------------------------------------------------------
    Field<double,Dim3> A3(layout3,sbc,gc);
    assign(A3[I][J][K], I + J + K);

    fdi << endl << "--------setFormat(8,3)-------" << endl;
    setFormat(8,3);

    fdi << endl << "--------fp3(A3)-------" << endl;
    fp3(A3,docomm);

    fdi << endl << "--------sfp3(A3,nx-1,1,0,ny-1,1,0,nz-1,1)-------" << endl;
    sfp3(A3,0,nx-1,1,0,ny-1,1,0,nz-1,1,docomm);

    fdi << endl << "--------sfp3(A3,nx-1,1,0,ny-1,2,0,nz-1,2)-------" << endl;
    sfp3(A3,0,nx-1,1,0,ny-1,1,0,nz-1,1,docomm);


    // Vector Field--------------------------------------------------------------
    Field<Vektor<double,Dim3>,Dim3> B3(layout3,vbc,gc);
    Vektor<double, Dim3 > Vinit3(1.0,2.0,3.0);
    assign(B3,Vinit3);

    fdi << endl << "--------setFormat(1,8)-------" << endl;
    setFormat(1,8);

    fdi << endl << "--------fp3(B3)-------" << endl;
    fp3(B3,docomm);

    // Write out "by hand" into another file what the previous field-printing
    // functions should have produced; this will be compared with what they
    // actually did produce:
    hardCodedOutput("text.correct.TestFieldDebug");

    // Compare the two files by mocking up the Unix "diff" command:
    delete fdip;
    bool passed = thediff("text.test.TestFieldDebug","text.correct.TestFieldDebug");
    EXPECT_TRUE(passed);
}

namespace {
    //-----------------------------------------------------------------------------
    // Mock up the Unix "diff" utility to compare two files:
    //-----------------------------------------------------------------------------
    bool thediff(std::string filename1, std::string filename2)
    {
        bool same = true;
        char ch1, ch2;
        std::ifstream file1(filename1);
        std::ifstream file2(filename2);
        while (file1.get(ch1)) {          // Read file 1 char-by-char until eof
            if (file2.get(ch2)) {           // Read equivalent char from file 2
                if (ch1 != ch2) same = false; // If they're different,files are different
            } else {
                same = false;                 // If file 2 ends before file 1, different
            }
        }
        return(same);
    }

    //-----------------------------------------------------------------------------
    void hardCodedOutput(std::string filename)
    {
        std::ofstream of(filename);
        of << std::endl
           << "--------setFormat(8,3)-------" << std::endl
           << "" << std::endl
           << "--------fp3(A3)-------" << std::endl
           << "~~~~~~~~ field slice (0:3:1, 0:3:1, 0:3:1) ~~~~~~~~" << std::endl
           << "==================================================K = 0" << std::endl
           << "--------------------------------------------------J = 0" << std::endl
           << "0.000e+00 1.000e+00 2.000e+00 3.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 1" << std::endl
           << "1.000e+00 2.000e+00 3.000e+00 4.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 2" << std::endl
           << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 3" << std::endl
           << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << std::endl
           << "" << std::endl
           << "==================================================K = 1" << std::endl
           << "--------------------------------------------------J = 0" << std::endl
           << "1.000e+00 2.000e+00 3.000e+00 4.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 1" << std::endl
           << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 2" << std::endl
           << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 3" << std::endl
           << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << std::endl
           << "" << std::endl
           << "==================================================K = 2" << std::endl
           << "--------------------------------------------------J = 0" << std::endl
           << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 1" << std::endl
           << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 2" << std::endl
           << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 3" << std::endl
           << "5.000e+00 6.000e+00 7.000e+00 8.000e+00 " << std::endl
           << "" << std::endl
           << "==================================================K = 3" << std::endl
           << "--------------------------------------------------J = 0" << std::endl
           << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 1" << std::endl
           << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 2" << std::endl
           << "5.000e+00 6.000e+00 7.000e+00 8.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 3" << std::endl
           << "6.000e+00 7.000e+00 8.000e+00 9.000e+00 " << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------sfp3(A3,nx-1,1,0,ny-1,1,0,nz-1,1)-------" << std::endl
           << "~~~~~~~~ field slice (0:3:1, 0:3:1, 0:3:1) ~~~~~~~~" << std::endl
           << "==================================================K = 0" << std::endl
           << "--------------------------------------------------J = 0" << std::endl
           << "0.000e+00 1.000e+00 2.000e+00 3.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 1" << std::endl
           << "1.000e+00 2.000e+00 3.000e+00 4.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 2" << std::endl
           << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 3" << std::endl
           << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << std::endl
           << "" << std::endl
           << "==================================================K = 1" << std::endl
           << "--------------------------------------------------J = 0" << std::endl
           << "1.000e+00 2.000e+00 3.000e+00 4.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 1" << std::endl
           << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 2" << std::endl
           << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 3" << std::endl
           << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << std::endl
           << "" << std::endl
           << "==================================================K = 2" << std::endl
           << "--------------------------------------------------J = 0" << std::endl
           << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 1" << std::endl
           << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 2" << std::endl
           << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 3" << std::endl
           << "5.000e+00 6.000e+00 7.000e+00 8.000e+00 " << std::endl
           << "" << std::endl
           << "==================================================K = 3" << std::endl
           << "--------------------------------------------------J = 0" << std::endl
           << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 1" << std::endl
           << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 2" << std::endl
           << "5.000e+00 6.000e+00 7.000e+00 8.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 3" << std::endl
           << "6.000e+00 7.000e+00 8.000e+00 9.000e+00 " << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------sfp3(A3,nx-1,1,0,ny-1,2,0,nz-1,2)-------" << std::endl
           << "~~~~~~~~ field slice (0:3:1, 0:3:1, 0:3:1) ~~~~~~~~" << std::endl
           << "==================================================K = 0" << std::endl
           << "--------------------------------------------------J = 0" << std::endl
           << "0.000e+00 1.000e+00 2.000e+00 3.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 1" << std::endl
           << "1.000e+00 2.000e+00 3.000e+00 4.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 2" << std::endl
           << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 3" << std::endl
           << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << std::endl
           << "" << std::endl
           << "==================================================K = 1" << std::endl
           << "--------------------------------------------------J = 0" << std::endl
           << "1.000e+00 2.000e+00 3.000e+00 4.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 1" << std::endl
           << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 2" << std::endl
           << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 3" << std::endl
           << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << std::endl
           << "" << std::endl
           << "==================================================K = 2" << std::endl
           << "--------------------------------------------------J = 0" << std::endl
           << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 1" << std::endl
           << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 2" << std::endl
           << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 3" << std::endl
           << "5.000e+00 6.000e+00 7.000e+00 8.000e+00 " << std::endl
           << "" << std::endl
           << "==================================================K = 3" << std::endl
           << "--------------------------------------------------J = 0" << std::endl
           << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 1" << std::endl
           << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 2" << std::endl
           << "5.000e+00 6.000e+00 7.000e+00 8.000e+00 " << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 3" << std::endl
           << "6.000e+00 7.000e+00 8.000e+00 9.000e+00 " << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------setFormat(1,8)-------" << std::endl
           << "" << std::endl
           << "--------fp3(B3)-------" << std::endl
           << "~~~~~~~~ field slice (0:3:1, 0:3:1, 0:3:1) ~~~~~~~~" << std::endl
           << "==================================================K = 0" << std::endl
           << "--------------------------------------------------J = 0" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 1" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 2" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 3" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "" << std::endl
           << "" << std::endl
           << "==================================================K = 1" << std::endl
           << "--------------------------------------------------J = 0" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 1" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 2" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 3" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "" << std::endl
           << "" << std::endl
           << "==================================================K = 2" << std::endl
           << "--------------------------------------------------J = 0" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 1" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 2" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 3" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "" << std::endl
           << "" << std::endl
           << "==================================================K = 3" << std::endl
           << "--------------------------------------------------J = 0" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 1" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 2" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------------------------------------------------J = 3" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << std::endl
           << "" << std::endl
           << "" << std::endl;
        of.close();
    }
}
