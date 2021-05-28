#include <gismo.h>
using namespace gismo;
int main(int, char**)
{
    
#if 0
    gsInfo.precision(3);
    gsKnotVector<>   kv (0, 1, 1, 4, 1 );
    
    
    gsInfo << "information of kv: " << "\n";
    kv.print(gsInfo);
    gsInfo << "\n\n";
    
    
    gsBSplineBasis<> bsp(kv);
    
    gsInfo << "detail of bsp: " << "\n";
    gsInfo << bsp.detail() << "\n";
    
#endif    
    
    
    std::string input("domain1d/bspline1d_01.xml");
    
    gsInfo << "creating an object of gsFileData" << "\n";
    gsFileData<> fileData(input);
    
    fileData.print(gsInfo);

#if 1
    
    gsInfo << "creating a pointer of gsGeometry" << "\n";
    gsGeometry<>::uPtr patch;
    
    gsInfo << "initializing patch" << "\n";
    
    if (fileData.has< gsGeometry<> >())
    {
        patch = fileData.getFirst< gsGeometry<> >();
    }
    else
    {
        gsInfo << "Input file doesn't have a geometry inside." << "\n";
        return -1;
    }

    if (!patch)
    {
        gsInfo << "Didn't find any geometry." << "\n";
        return -1;
    }
    
#endif    
    
#if 0
    
    gsInfo << "creating a pointer of gsBspline" << "\n";
    gsBSpline<>::uPtr patch;    
    
    patch = fileData.getFirst< gsBSpline<> >();
    
#endif
              
              
#if 0
              
    gsNurbsCreator<>::TensorBSpline2Ptr patch;
              
    patch = gsNurbsCreator<>::BSplineSquare(2,1);
    
#endif
    
    gsInfo << "\npatch: \n" << *patch << "\n";
    
    gsMatrix<> support = patch->support();
    gsInfo << "Support: \n"
           << support << "\n\n";    
    
    gsMesh<> mesh;
    patch->controlNet(mesh);
    gsInfo << "\nMesh: \n" << mesh << "\n";
    
    const gsBasis<>& basis = patch->basis();
    gsInfo << "\nBasis: \n" << basis << "\n";
    
    gsInfo << "\ndegree: " << patch->degree(0) << "\n";

    const gsMatrix<>& coefs = patch->coefs();
    gsInfo << "\nCoefficients: \n" << coefs << "\n" << "\n";

    gsInfo << "Dimension of the parameter space: " << patch->parDim() << "\n"
              << "Dimension of the geo patch: " << patch->geoDim() << "\n";      
              
              
#if 1
    std::string output = "test_bsplines_";
    
    std::string out = output + "Geometry";
    gsInfo << "Writing the geometry to a paraview file: " << out
                << "\n\n";
    gsWriteParaview(*patch, out);   
    
    out = output + "Basis";
    gsInfo << "Writing the basis to a paraview file: " << out
                << "\n\n";
    gsWriteParaview(basis, out);
    
//     
//     out = output + "ContolNet";
//     gsInfo << "Writing the control net to a paraview file: " << out
//                 << "\n" << "\n";
//     gsWriteParaview(mesh, out); 

#endif
    
    return 0;
}
