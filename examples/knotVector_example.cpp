/** @file knotVector_example.cpp

    @brief Tutorial on gsKnotVector class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#include <iostream>
#include <string>
#include <algorithm>
#include <gismo.h>

// #include </home/dv/jliu/Localdisk/A_Process/01_custom_header/my_print.h>
// #include </home/dv/jliu/Localdisk/workspace/gismo/src/gsKnotVector.h>
// #include </home/dv/jliu/Localdisk/workspace/gismo/src/gsNurbs/test.h>


// #include </home/dv/jliu/Localdisk/workspace/gismo/build/include_custom/gsNurbs/gsKnotVector.h>



using namespace gismo;
using namespace std;

// forward declaration of some utility functions
void printKnotVector(const gsKnotVector<>& kv, const std::string& name);
void printKnotVector(const gsKnotVector<>& kv);
std::vector<real_t> makeVectorOfKnots();
void print(const real_t& el);


int main(int argc, char* argv[])
{
    gsCmdLine cmd("Tutorial on gsKnotVector class.");
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    // ======================================================================
    // different construction of a knot vector
    // ======================================================================

    gsInfo << "------------- Constructions -----------------------------\n";

//     gsKnotVector<> kv0;
//     printKnotVector(kv0, "kv0");
    

//     // only with degree
//     gsKnotVector<> kv1(3);
//     printKnotVector(kv1, "kv1");

    // interval [0,1], 3 interior knots, multiplicity 3 at the ends (clamped knots).
    gsKnotVector<> kv2(0,1,3,3);                                // according to the definition, the template name of class gsKnotVector is indicated from the type of the first argumant, see gsKnotVector.h apr 23, 2020
    
    printKnotVector(kv2, "kv2");

    real_t a = 0.0; // starting knot
    real_t b = 1.0; // ending knot
    unsigned interior = 4; // number of interior knots
    unsigned multEnd = 3; // multiplicity at the two end knots              // the size of multiplicity determines the degree, i.e. p=multEnd-1, Nov 20, 2019
    gsKnotVector<> kv3(a, b, interior, multEnd);
    printKnotVector(kv3, "kv3");

    std::vector<real_t> knots = makeVectorOfKnots();                        // definition below, allowing knots set manually by user
    
    cout << "knot vector defined by user: ";
    for (unsigned int i = 0; i<knots.size(); ++i)
    {
        cout << knots[i] << " "; 
    }
    gsInfo << "\n\n";
    
    
    
    gsKnotVector<> kv4(knots, 2, 1); // knots, degree, regularity                       // line 582 in gsKnotVector.hpp, for now we only consider the regular case
    printKnotVector(kv4, "kv4");

//     gsKnotVector<> kv5(knots, 2); // knots, degree
//     printKnotVector(kv5, "kv5");

    gsKnotVector<> kv6;
    kv6.initUniform(5, 3); // number of knots, multiple ends
    printKnotVector(kv6, "kv6");

    gsKnotVector<> kv7;
    kv7.initClamped(a, b, 3, 5); // start, end, degree, number of interior knots
    printKnotVector(kv7, "kv7");


    // ======================================================================
    // looping over knots
    // ======================================================================

    // looping over all knots

//     gsInfo << "\n\n"
//               << "------------- Looping over knots -----------------------\n"
//               << "kv7: \n";
//     for (gsKnotVector<>::iterator it = kv7.begin(); it != kv7.end(); it++)
//     {
//         gsInfo << *it << " ";
//     }
//     gsInfo << "\n\n";
// 
//     // looping over unique knots
//     for (gsKnotVector<>::uiterator it = kv7.ubegin(); it != kv7.uend(); it++)
//     {
//         gsInfo << *it << " ";
//     }
//     gsInfo << "\n\n\n";


    // ======================================================================
    // some properties
    // ======================================================================

#if 1

//     gsInfo << "------------- Some properties    -----------------------\n"
//               << "kv7: \n\n";
// 
//     printKnotVector(kv7);
// 
//     gsInfo << "kv7.size(): " << kv7.size() << "\n\n"
//               << "kv7.findspan(1.5): " << kv7.iFind(1.5) - kv7.begin() << "\n\n"
//               << "kv7.findspan(2): " << kv7.iFind(2) - kv7.begin() << "\n\n"
//               << "kv7.has(2): " << kv7.has(2) << "\n\n"
//               << "kv7.has(2.1): " << kv7.has(2.1) << "\n\n"
//               << "kv7.isUniform(): " << kv7.isUniform() << "\n\n"
//               << "kv7.numKnotSpans(): " << kv7.uSize() - 1 << "\n\n"
//               << "kv7.isOpen(): " << kv7.isOpen() << "\n\n"
//               << "kv7.multiplicity(2): " << kv7.multiplicity(2) << "\n\n"
//               << "kv7.multiplicity(1): " << kv7.multiplicity(1) << "\n\n\n";


    // ======================================================================
    // some operations
    // ======================================================================

    gsInfo << "------------- Some operations    -----------------------\n";
//     printKnotVector(kv6, "kv6");
// 
// 
//     std::vector<real_t> unique = kv6.unique();
//     gsInfo << "\nUnique knots: \n";
//     std::for_each(unique.begin(), unique.end(), print);

    gsMatrix<>* greville = kv6.greville();
    gsInfo << "\nGreville points: \n" << *greville << "\n\n";
    delete greville;


//     std::vector<index_t> mult = kv6.multiplicities();
//     gsInfo << "Multiplicities: ";
//     std::for_each(mult.begin(), mult.end(), print);
//     gsInfo << "\n\n";
// 
//     printKnotVector(kv6, "kv6");
// 
//     gsInfo << "kv6.uniformRefine()\n";
//     kv6.uniformRefine();
//     printKnotVector(kv6);
// 
//     gsInfo << "kv6.degreeElevate()\n";
//     kv6.degreeElevate();
//     printKnotVector(kv6);

//     gsInfo << "For other capabilites of gsKnotVector look at "
//         "src/gsNurbs/gsKnotVector.h\n" << "\n";

#endif
        
    return 0;
}

void print(const real_t& el)
{
    gsInfo << el << " ";
}


void printKnotVector(const gsKnotVector<>& kv,
                     const std::string& name)
{
    gsInfo << name << ":\n";
    kv.print(gsInfo);
    gsInfo << "\n" << "\n";
}


void printKnotVector(const gsKnotVector<>& kv)
{
    for (gsKnotVector<>::const_iterator it = kv.begin(); it != kv.end(); it++)
    {
        gsInfo << *it << " ";
    }
    gsInfo << "\n\n";
}



std::vector<real_t> makeVectorOfKnots()
{
    std::vector<real_t> knots;
    knots.push_back(0);
    knots.push_back(0.1);
    knots.push_back(0.5);
    knots.push_back(0.6);
    knots.push_back(0.9);
    knots.push_back(1);

    return knots;
}
