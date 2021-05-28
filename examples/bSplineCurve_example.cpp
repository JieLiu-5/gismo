/** @file bSplineCurve_example.cpp

    @brief Tutorial on gsBSpline class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot = false; // If set to true, paraview file is generated and launched on exit


    gsCmdLine cmd("Tutorial 01 shows the use of BSpline curves.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    // Make a BSpline curve
    gsKnotVector<> kv(0, 1, 1, 3);//start,end,interior knots, start/end multiplicites of knots1
    
    gsInfo << "\n";
    gsInfo << "The knot vector reads:\n";
    gsInfo << kv;
    gsInfo << "\n\n";
    
    
    
    gsMatrix<> coefs(4, 3);                                         // unclear how coefs affect the curve, may 5 2020
    coefs << 0, 0, 0,
             1, 1, 0,
             2, 0, 1,
             4, 1, 0;
    gsBSpline<> curve( kv, give(coefs));
    
//     gsInfo << "basis of the curve:\n";
//     gsInfo << curve.basis() << "\n";
    
    
/*    gsInfo << "output basis functions\n";
    gsWriteParaview(curve.basis(), "bsplinecurve_basis");   */ 

    // Print the Bspline curve
    gsInfo << "\n";
    gsInfo << "I am a " << curve << "\n";
    
    
    gsInfo << "coefs of curve:\n";
    gsInfo << curve.coefs() << "\n";
    
    gsInfo << "\n";
    gsInfo << "number of coefs: " 
           << curve.coefsSize() << "\n"
           << "dimension of the parameter domain: "
           << curve.domainDim() << "\n"
           << "dimension of the ambient physical space: "
           << curve.targetDim() << "\n"
           << "support(): "
           << curve.support() << "\n"
           << "\n";
           
           
    unsigned n_interpolation = 10;
    
    gsMatrix<> pts = gsPointGrid(curve.support()[0],curve.support()[1],n_interpolation) ;
           
    gsMatrix<>  eval_func = curve.eval  ( pts ) ;//pts 
           
           
    gsInfo << "coordinates and function values on the points defined by support points: \n";
    for ( index_t j=0; j<pts.cols(); ++j)
    {
        gsInfo << "col " << j << ": ";
        gsInfo << pts(j) << " ";
        for ( index_t i=0; i<eval_func.rows(); ++i)
        {
            gsInfo << eval_func(i,j) << " ";
        }
        gsInfo << "\n";
    }    
    gsInfo << "\n";
           

    if (plot)
    {
        // Output a paraview file
        gsWriteParaview( curve, "bsplinecurve", n_interpolation);
//         gsFileManager::open("bsplinecurve.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    return 0;
}
