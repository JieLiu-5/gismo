#include <gismo.h>
using namespace gismo;
int main(int, char**)
{
    
#if 0
    
//     std::string input("planar/two_squares.xml");
    std::string input("domain2d/square.xml");
        
    
    gsInfo << "creating an object of gsFileData" << "\n";
    gsFileData<> fileData(input);
    
    fileData.print(gsInfo);

    gsMultiPatch<>::uPtr patch;
    
    patch = fileData.getFirst<gsMultiPatch<>>();
    
    gsInfo << "\npatch: \n" << *patch << "\n";
    
#endif

    gsMultiPatch<> patch;
    
    unsigned int id_dimension = 0;
    
    
    std::string fn;
    
    std::string output = "test_multipatch_";
    
    if (id_dimension == 0)
    {
//         gsReadFile<>("domain1d/bspline1d_01.xml", patch);
        
        fn="domain1d/bspline1d_01.xml";
        output = output + "1d_";
    }else
    {
//         gsReadFile<>("domain2d/square.xml", patch);
//         gsReadFile<>("planar/test_two_squares.xml", patch);
//         gsReadFile<>("pde/poisson2d_bvp.xml", patch);
    
//         fn="pde/poisson2d_bvp.xml";
        fn="domain2d/square.xml";
//         fn="planar/test_one_square.xml";
        
        output = output + "2d_";
    }
    
    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n\n";
    
    gsGeometry<>::uPtr pGeom=fd.getFirst< gsGeometry<> >();
    gsInfo << "geometry info: \n" << *pGeom << "\n";
    
//     fd.getId(0, patch);    
    patch=*pGeom;
    
    gsInfo << "\n";
    gsInfo << "patch.detail() \n" << patch.detail() << "\n\n";
    
    
    
    for (unsigned int i=0; i<patch.nPatches(); ++i)
    {
        gsInfo << "Info of the " << i << "th patch:\n";
        gsInfo << "coefs(): \n" << patch.patch(i).coefs() << "\n";
        gsInfo << "support(): \n" << patch.patch(i).support() << "\n";
        gsInfo << "\n";
    }
    
    
    gsFunctionExpr<> f;
    fd.getId(1, f); // id=1: source function
    gsInfo << "\n";
    gsInfo<<"Source function "<< f << "\n";

    gsBoundaryConditions<> bc;
    fd.getId(2, bc); // id=2: boundary conditions
    gsInfo << "\n";
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";
    for(unsigned int i=0; i< bc.dirichletSides().size(); i++)
    {
        gsInfo << i << ": " << *bc.dirichletSides()[i].function() << "\n";
    }
    gsInfo << "\n";
    
    gsFunctionExpr<> ms;
    fd.getId(3, ms); // id=3: reference solution
    gsInfo<<"Exact solution: "<< ms << "\n";
    gsInfo << "\n";     
    
    gsOptionList Aopt;
    fd.getId(4, Aopt); // id=4: assembler options
    gsInfo << "\n";
    gsInfo << "Aopt.print(gsInfo)\n";
    Aopt.print(gsInfo);   
    
    
    gsInfo << "\n";
    std::string out = output + "Geometry";
    gsInfo << "Writing the geometry to a paraview file: " << out
                << "\n";
    gsWriteParaview(patch, out);  
    
    
    const gsBasis<>& basis = patch.basis(0);
    gsInfo << "\nBasis: \n" << basis << "\n";        
    
    gsInfo << "\n";
    out = output + "Basis";
    gsInfo << "Writing the basis to a paraview file: " << out
                << "\n";
    gsWriteParaview(basis, out);      
    
    return 0;
}
