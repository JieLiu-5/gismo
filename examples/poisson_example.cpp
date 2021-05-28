/** @file poisson_example.cpp

    @brief Tutorial on how to use G+Smo to solve the Poisson equation,
    see the \ref PoissonTutorial

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

//! [Include namespace]
# include <gismo.h>

#include <gsAssembler/gsAssembler.h>

using namespace gismo;
using namespace std;

//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    int is_2d=0;
    
    //! [Function data]
    
    // Define source function
    
    gsFunctionExpr<> f("-exp(-(x-0.5)^2)*(4*x^2-4*x-1)",1);
    gsFunctionExpr<> g("exp(-(x-0.5)^2)",1);

    // Print out source function and solution
    gsInfo<<"Source function "<< f << "\n";
    gsInfo<<"Exact solution "<< g <<"\n\n";
    //! [Function data]

    gsInfo << "\n";
    gsInfo << "############################################" << "\n";
    gsInfo << "reading the xml file\n";
    
    //! [Geometry data]
    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patches;
    
    if(is_2d==1)
    {
        
        gsReadFile<>("domain2d/square.xml", patches);                          //square       lake
    
    }else if(is_2d==0)
    {
        gsReadFile<>("domain1d/bspline1d_01.xml", patches);           // an object of class gsFiledata is created as a private member in gsReadFile                  lshape_p2 read_temp 
    }

    gsInfo << "\n";
    gsInfo << "############################################" << "\n";
    gsInfo << "patch information: \n";
        
    gsInfo << "The domain is a "<< patches <<"\n";                      // this syntax prints the number of boundaries and interfaces, Nov 20, 2019
    
    
    gsInfo << "  basis(0): \n" << patches.basis(0) << "\n";             // basis(i) returning the basis of the i-th patch, which is of type class gsBasis
    gsInfo << "\n";
    
    cout << "patches.detail(): \n" << patches.detail() << "\n\n";
    
    gsInfo << "  nBoundary(): " << patches.nBoundary() << "\n";
    gsInfo << "  nBoxes(): " << patches.nBoxes() << "\n";
    gsInfo << "  nInterfaces(): " << patches.nInterfaces() << "\n";
    gsInfo << "  nPatches(): " << patches.nPatches() << "\n";
    gsInfo << "  nPieces(): " << patches.nPieces() << "\n";
    
    
    gsInfo << "  coefs(): \n" << patches.patch(0).coefs() << "\n";                  // patch(i) returns an object of gsGeometry
    gsInfo << "  support(): " << patches.patch(0).support() << "\n";
    
    
//     gsInfo << "boundaries(): \n" << patches.boundaries() << "\n";
    

    
    //! [Geometry data]

    gsInfo << "\n";
    gsInfo << "############################################" << "\n";
    gsInfo << "boundary information: \n";
    
    //! [Boundary conditions]
    gsBoundaryConditions<> bcInfo;                                  // this class contains objects of class BoundaryConditions<T>
    
    
    // Every patch with a boundary need to be specified. In this
    // there are in total 8 sides (two for each patch)

    if(is_2d==1)
    {
    
            // Neumann Boundary conditions
            gsFunctionExpr<> hEast ("exp(-(x-0.5)^2)*(-2*(x-0.5))", "2*(x-0.5)",2);        // 1*pi*cos(pi*1)*sin(pi*2*y)   3*pi*cos(pi*3)*sin(pi*4*y)
            gsFunctionExpr<> hSouth("0","0",2);      // -pi*2*sin(pi*x*1)    -pi*4*sin(pi*x*3)

        #if 0
            // Dirichlet Boundary conditions
            // First argument is the patch number
            bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, &g);             // for 4 patches
            bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, &g);

            bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, &g);
            bcInfo.addCondition(3, boundary::north, condition_type::dirichlet, &g);

            bcInfo.addCondition(3, boundary::east,  condition_type::neumann, &hEast);
            bcInfo.addCondition(2, boundary::east,  condition_type::neumann, &hEast);

            bcInfo.addCondition(0, boundary::south, condition_type::neumann, &hSouth);
            bcInfo.addCondition(2, boundary::south, condition_type::neumann, &hSouth);
            //! [Boundary conditions]
        #endif

        #if 1

            bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, &g);             // for one patch
            bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &g);
            bcInfo.addCondition(0, boundary::east,  condition_type::neumann, &hEast);
            bcInfo.addCondition(0, boundary::south, condition_type::neumann, &hSouth);
            
        #endif

    }else if(is_2d==0)
    {

        bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, &g);
        bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, &g);    
    
    }

    /*
      //Alternatively: You can automatically create Dirichlet boundary
      //conditions using one function (the exact solution) for all
      //boundaries like this:

    for (gsMultiPatch<>::const_biterator
             bit = patches.bBegin(); bit != patches.bEnd(); ++bit)
    {
        bcInfo.addCondition( *bit, condition_type::dirichlet, &g );
    }
    */
    
    
    
    gsInfo << "\n";
    gsInfo << "bcInfo.print(gsInfo): \n";   
    bcInfo.print(gsInfo);                           // this prints the number of certain type of boundary conditions that are imposed
    
    gsInfo << "\n";
           
    gsInfo << "size of bcInfo.dirichletSides(): " << bcInfo.dirichletSides().size() << "\n";        
                                                    // dirichletSides() returning data of type bcContainer, that is, std::deque<boundary_condition<T> >
           
    gsInfo << "\n";
           
           
    gsInfo << "*bcInfo.dirichletSides()\n";
    
    gsInfo << "[0] type: " << bcInfo.dirichletSides()[0].ctype()                            // content of first Dirichlet boundary condition
           << ", function: " << *bcInfo.dirichletSides()[0].function()
           << ", side: " << bcInfo.dirichletSides()[0].side()
           << ", patch: " << bcInfo.dirichletSides()[0].patch()
           << ", isHomogeneous? " << bcInfo.dirichletSides()[0].isHomogeneous()
           << ", parametric? " << bcInfo.dirichletSides()[0].parametric()
           << "\n";
           
    gsInfo << "[1] type: " << bcInfo.dirichletSides()[1].ctype()                            // content of second Dirichlet boundary condition
           << ", function: " << *bcInfo.dirichletSides()[1].function()
           << ", side: " << bcInfo.dirichletSides()[1].side()
           << ", patch: " << bcInfo.dirichletSides()[1].patch()
           << ", isHomogeneous? " << bcInfo.dirichletSides()[1].isHomogeneous()
           << ", parametric? " << bcInfo.dirichletSides()[1].parametric()
           << "\n";
           
           
    gsInfo << "\n";
    
    gsInfo << "size of bcInfo.neumannSides(): " << bcInfo.neumannSides().size() << "\n";
           
           
           
#if 1
           
    
    gsInfo << "\n";
    gsInfo << "############################################" << "\n";
    gsInfo << "basis information, which is initialized by patches\n";
    
    // Number for h-refinement of the computational (trial/test) basis.
    const int numRefine  = 0;

    // Number for p-refinement of the computational (trial/test) basis.
    const int degree     = 1;    

    
    
    //! [Refinement]
    // Copy basis from the geometry
    gsMultiBasis<> refine_bases( patches );
    
    gsInfo << "  refine_bases.degree(0): " << refine_bases.degree(0) << "\n";
    
    
//     gsMatrix<real_t> u(1, 1);
//     gsInfo << "  u reads:\n";
//     gsInfo << u << "\n";
//     gsInfo << refine_bases.eval(u);
    
    cout << "  refine_bases.print(gsInfo): " << "\n";
    refine_bases.print(gsInfo);
    gsInfo << "\n\n";
    
    
    gsInfo << "  refine_bases.basis(0): \n" << refine_bases.basis(0) << "\n";                   // basis(i) returning the i-th basis block, which is of type class gsBasis 

//     gsInfo << "refine_bases.size(): " << refine_bases.size() << "\n";
//     gsInfo << "refine_bases.maxCwiseDegree(): " << refine_bases.maxCwiseDegree() << ", minCwiseDegree(): " << refine_bases.minCwiseDegree() << "\n";
//     
    
    gsInfo << "\n";
    cout << "numRefine: " << numRefine << ", degree: " << degree << endl;
    gsInfo << "\n";

    // h-refine each basis (4, one for each patch)
    for ( int i = 0; i < numRefine; ++i)
      refine_bases.uniformRefine();

    // k-refinement (set degree)
    for ( size_t i = 0; i < refine_bases.nBases(); ++ i )
        refine_bases[i].setDegreePreservingMultiplicity(degree);

    //! [Refinement]
    
    gsInfo << "  refine_bases.basis(0): \n" << refine_bases.basis(0) << "\n";
    
    ////////////// Setup solver and solve //////////////
    // Initialize Solver
    // Setup method for handling Dirichlet boundaries, options:
    //
    // * elimination: Eliminate the Dirichlet DoFs from the linear system.
    //
    // * nitsche: Keep the Dirichlet DoFs and enforce the boundary
    //
    // condition weakly by a penalty term.
    // Setup method for handling patch interfaces, options:
    //
    // * glue:Glue patches together by merging DoFs across an interface into one.
    //   This only works for conforming interfaces
    //
    // * dg: Use discontinuous Galerkin-like coupling between adjacent patches.
    //       (This option might not be available yet)
    //! [Assemble]
    
    
//     gsInfo << "  coefs(): \n" << patches.patch(0).coefs() << "\n"; 
//     gsInfo << "  support(): " << patches.patch(0).support() << "\n";
    
    gsInfo << "\n";
    gsInfo << "############################################" << "\n";
    gsPoissonPde<real_t> pde(patches,bcInfo,f);                                             // members belonging to class gsPde
    gsInfo << "\n";
    
    gsInfo << "pde information: initialized by patches (containing basis, number of boundaries, coefs, support), bcInfo and f\n";
    
    gsInfo << "\n";
    gsInfo << "  pde.domain(): " << pde.domain();                                           // domain() returns an object of class gsMultiPatch
    gsInfo << "\n";
    gsInfo << "  pde.boundaryConditions(): " << pde.boundaryConditions();                   // boundaryConditions() returns an object of class gsBoundaryConditions
    gsInfo << "\n";
    
    gsInfo << "  pde.dim(): " << pde.dim() << "\n";                                         // 
    gsInfo << "  pde.numRhs(): " << pde.numRhs() << "\n";
    gsInfo << "  pde.rhs()->print(gsInfo): ";                                               // rhs() returning a pointer to the first piece of the rhs, which is of type gsFunction
    pde.rhs()->print(gsInfo);
    gsInfo << "\n";
    
    gsInfo << "  pde.rhs()->basis(0): " << pde.rhs()->basis(0) << "\n";                     // basis(i) returning an object of gsBasis
    
    
    gsInfo << "\n";
    gsInfo << "  pde.unknownDim(): " << pde.unknownDim() << "\n";    
    
    gsInfo << "\n";
    
    gsInfo << "############################################" << "\n";
    
    gsInfo << "creating an object of gsPoissonAssembler<T>\n\n";
    
    gsInfo << "dirichlet::elimination: " << dirichlet::elimination << "\n";
    gsInfo << "dirichlet::nitsche: " << dirichlet::nitsche << "\n";
    gsInfo << "dirichlet::penalize: " << dirichlet::penalize << "\n";
    
    gsInfo << "\n";
    
    gsInfo << "dirichlet::homogeneous: " << dirichlet::homogeneous << "\n";
    gsInfo << "dirichlet::interpolation: " << dirichlet::interpolation << "\n";
    gsInfo << "dirichlet::l2Projection: " << dirichlet::l2Projection << "\n";
    gsInfo << "dirichlet::user: " << dirichlet::user << "\n";
    
    
    gsInfo << "\n";
    
    gsInfo << "type of iFace: " << typeid(iFace).name() << "\n";
    gsInfo << "iFace::glue: " << iFace::glue << "\n";
    
    gsInfo << "\n";
           
           
    gsPoissonAssembler<real_t> assembler(pde,refine_bases,                                 // invoking gsAssembler<T>::initialize()
                                         dirichlet::elimination    , iFace::glue);      //
//                                          dirichlet::nitsche    , iFace::glue);             // in gsAssembler.h, iFace is a reference to an object of class boundaryInterface, which is defined in gsBoundary.h
    
    
    gsInfo << "assembler.options(): \n"
           << "size(): " << assembler.options().size() << "\n"
           << "DirichletStrategy: " << assembler.options().getInt("DirichletStrategy") << "\n"
           << "DirichletValues: " << assembler.options().getInt("DirichletValues") << "\n"
           << "InterfaceStrategy: " << assembler.options().getInt("InterfaceStrategy") << "\n"
//            << "quRule: " << assembler.options().getInt("quRule") << "\n"
           << "bdA: " << assembler.options().getReal("bdA") << "\n"
           << "quA: " << assembler.options().getReal("quA") << "\n"; 
           
    gsInfo << "\n";
    assembler.options().print_custom();
    
    
    gsInfo << "\n";
    gsInfo << "assembler.options().print(gsInfo)\n";
    assembler.options().print(gsInfo);
    
    gsInfo << "\n";
    gsInfo << "############################################" << "\n";
    assembler.assemble();                               // a function defined in gsPoissonAssembler.hpp
    
    gsInfo << "Have assembled a system (matrix and load vector) with " << assembler.numDofs() << " dofs.\n";
    //! [Assemble]
           
    gsInfo << "\n";
    gsInfo << "############################################" << "\n";
    
    cout << "matrix info("
         << assembler.matrix().dim().first << "*" << assembler.matrix().dim().second
         << "): " << endl;
    
    for (int i = 0; i < assembler.numDofs(); ++i)
    {
        for (int j = 0; j < assembler.numDofs(); ++j)
        {
            cout <<  assembler.matrix().at(i,j) << " ";    
        }
        cout << endl;
    }
    cout << endl;
    
    
    gsInfo << "rhs info\n";
    gsInfo << assembler.rhs();
    
    
#if 1
    gsInfo << "\n\n";
    gsInfo << "############################################" << "\n";
        
    //! [Solve]
    // Initialize the conjugate gradient solver
    gsInfo << "Solving...\n";
    gsSparseSolver<>::CGDiagonal solver( assembler.matrix() );
    gsMatrix<> solVector = solver.solve( assembler.rhs() );
    gsInfo << "Solved the system with CG solver.\n";
    //! [Solve]
    
    
    gsInfo << "\n";
    cout << "solution: \n";
    for (int i = 0; i < assembler.numDofs(); ++i)
    {
        cout << "[" << i << "]: " << solVector.at(i) << "\n";    
    }

    //! [Construct solution]
    // Construct the solution as a scalar field
    gsMultiPatch<> mpsol;
    assembler.constructSolution(solVector, mpsol);
    gsField<> sol( assembler.patches(), mpsol);
    //! [Construct solution]
    
    gsInfo << "\n";
    gsInfo << "solution field:\n";
//     sol.print(cout);
//     gsInfo << "\n";
    cout << "sol.nPieces(): " << sol.nPieces() << endl;
    cout << "sol.dim(): " << sol.dim() << endl;
    gsInfo << "\n";
    
    
//     =========================== operations for computing the error
    
    gsInfo << "\n";
    gsInfo << "############################################" << "\n";
    
    gsInfo << "Preparing for computing the error\n";
    
    gsFileData<> fd("domain1d/bspline1d_01.xml");
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";
    
    
    gsExprAssembler<> A(1,1);
    A.setOptions(assembler.options());


    typedef gsExprAssembler<>::geometryMap geometryMap;                         
                                                                                
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;    
    typedef gsExprAssembler<>::solution    solution;    
    
    // Elements used for numerical integration
    A.setIntegrationElements(refine_bases);
    gsExprEvaluator<> ev(A);   
    
    // Set the geometry map
    geometryMap G = A.getMap(patches);                                          // 
    
    // Set the discretization space
    space u = A.getSpace(refine_bases);
    u.setInterfaceCont(0);
    u.addBc( bcInfo.get("Dirichlet") );

    // Set the source term
    variable ff = A.getCoeff(f, G);
    
    gsInfo << "ff: \n" << ff << "\n";
    
    

    gsInfo << "\n";
    
    // Display manufactured solution
    gsInfo<<"Exact solution: "<< g << "\n";
    variable u_ex = ev.getVariable(g, G);

    gsInfo << "u_ex: " << u_ex << "\n";
    
    // Solution variable
    solution u_sol = A.getSolution(u, solVector);    
    
    gsVector<> l2err(numRefine+1), h1err(numRefine+1);
    
    
    
    
    gsInfo << "\n";
    gsInfo << "G.source(): \n" << G.source() << "\n";
    
    gsInfo << "G.data().points: \n" << G.data().points << "\n";
//     gsInfo << "G.data().values[0]: \n" << G.data().values[0].col(0) << "\n";
    gsInfo << "meas(G): " << meas(G) << "\n";
    
//     l2err[0]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
    
    gsInfo << "l2err: \n" << l2err << "\n";
    

    gsInfo << "\n";
    gsInfo << "\n";
    gsInfo << "\n";
    
    
#if 0
    if (plot)
    {
        //! [Plot in Paraview]
        // Write approximate and exact solution to paraview files
        
        gsInfo<<"Plotting in Paraview...\n";
        
        if(is_2d==1)
        {
            
            gsWriteParaview<>(sol, "poisson2d", 1000);
            const gsField<> exact( assembler.patches(), g, false );
            gsWriteParaview<>( exact, "poisson2d_exact", 1000);
            
        }else if(is_2d==0)
        {
            
            gsWriteParaview<>(sol, "poisson1d", 1000);
            const gsField<> exact( assembler.patches(), g, false );
            gsWriteParaview<>( exact, "poisson1d_exact", 1000);
                        
        }

        // Run paraview
//         gsFileManager::open("poisson2d.pvd");
//         gsFileManager::open("poisson2d_exact.pvd");
        //! [Plot in Paraview]
    }
    else
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    }
    return EXIT_SUCCESS;

#endif
    
#endif
    
#endif
}// end main
