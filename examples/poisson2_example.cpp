/** @file poisson2_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
#include <time.h>
#include <fstream>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    
    gsInfo << "content of argv: " << "\n";
    for (int i =0; i<argc; i++)
    {
        gsInfo << "[" << i << "]: " << *argv[i] << "\n";
    }
    gsInfo << "\n";
    
    //! [Parse command line]
    bool plot = false;
    index_t id_dimension = 0;
    index_t numRefine  = 1;
    index_t numElevate = 0;
    bool last = false;
    
    index_t is_matrix_stored = 0;
    
    std::string fn;
    
    
    unsigned int is_single_patch = 1;
    

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");                        // function defined in gsIO/gsCmdLine.h
//     gsInfo << cmd.getMessage() << "\n";
    
    cmd.addInt( "d", "dimension", "space dimension",  id_dimension );
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try 
    { 
//         gsInfo << "succeed in try\n";
        cmd.getValues(argc,argv);                                                   // using external/tclap/ValueArg.h,   mar 15 2020
    } 
    catch (int rv) 
    { 
//         gsInfo << "rv: " << rv << "\n";
        return rv; 
    }
    //! [Parse command line]
    
    
    if (id_dimension==0)
    {
        fn="domain1d/bspline1d_01.xml";
    } else
    {
        if (is_single_patch==1)
        {
            fn="domain2d/square.xml";
        }else{
            fn="pde/poisson2d_bvp.xml";
//             fn="planar/test_two_squares.xml";
        }
    }
    

    //! [Read input file]

    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    
    gsInfo << "\n";
    gsMultiPatch<> multipatch;
    
    if (is_single_patch==1)
    {
        gsGeometry<>::uPtr pGeom=fd.getFirst< gsGeometry<> >();
        multipatch=*pGeom;
    } else
    {
        fd.getId(0, multipatch); // id=0: Multipatch domain
    }
    
    
    gsInfo << "\n";
    for (unsigned int i = 0; i<multipatch.nPatches(); i++)
    {
        gsInfo << "patch " << i << "\n";
        gsInfo << "\n";
        gsInfo << "basis(): \n" << multipatch.patch(i).basis() << "\n";             
        gsInfo << "coefs(): \n" << multipatch.patch(i).coefs() << "\n";
        gsInfo << "support(): \n" << multipatch.patch(i).support() << "\n";
        gsInfo << "\n";
    }
    
    gsInfo << "\n";
    gsFunctionExpr<> f;
    fd.getId(1, f); // id=1: source function
    gsBoundaryConditions<> bc;
    fd.getId(2, bc); // id=2: boundary conditions
    gsOptionList Aopt;
    fd.getId(4, Aopt); // id=4: assembler options
    
    
    gsInfo << "\n";
    gsInfo<<"Source function "<< f << "\n";
//     gsInfo<<"Boundary conditions:\n";
//     for(unsigned int i=0; i< bc.dirichletSides().size(); i++)
//     {
//         gsInfo << i << " " << *bc.dirichletSides()[i].function() << "\n";
//     }
    
/*    gsFunctionExpr<> bc_neumann ("exp(-(x-0.5)^2)*(-2*(x-0.5))", "2*(x-0.5)",2);
    
    bc.addCondition(0, boundary::north,  condition_type::neumann, &bc_neumann);
    bc.addCondition(0, boundary::south,  condition_type::neumann, &bc_neumann);  */  
    
    bc.print(gsInfo);
    
//     gsInfo << "OptionList info: \n";
//     Aopt.print(gsInfo);
    
    gsInfo << "number of quadrature points: ";
    gsInfo << Aopt.getReal("quA")*(numElevate+1)+Aopt.getInt("quB") << "\n";
    
    
//     gsInfo << "numRefine: " << numRefine << "\n"
//            << "numElevate: " << numElevate << "\n";
           
    gsInfo << "\n";
    gsMultiBasis<> dbasis(multipatch);                              // gsMultiBasis<T> inheriting from gsFunctionSet<T>
    
    
    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);

    // h-refine each basis
    if (last)
    {
        gsInfo << "h-refine each basis\n";
        for (int r =0; r < numRefine-1; ++r)
            dbasis.uniformRefine();
        numRefine = 0;
    }
    
    gsInfo << "\n";
    gsInfo << "dbasis after p-refine and h-refine: \n";
    for (unsigned int i = 0; i<dbasis.nBases(); i++)
    {
        gsInfo << "[" << i << "] \n" << dbasis.basis(i) << "\n";             
    }    

# if 1
    
    gsInfo << "\n";
    gsInfo << "############################################" << "\n";

    gsExprAssembler<> A(1,1);                                      // end up with mutVar.setData()
    
    gsInfo << "\n";
    A.setOptions(Aopt);                                            

    typedef gsExprAssembler<>::geometryMap geometryMap;            // gsExprAssembler<>::geometryMap is an alias of gsExprHelper< T >::geometryMap 
                                                                   // in /src/gsAssembler/gsExprHelper.h,
                                                                   // gsExprHelper< T >::geometryMap is an alias of expr::gsGeometryMap<T>, see /src/gsAssembler/gsExprHelper.h
                                                                   // expr is a namespace defined in src/gsAssembler/gsExpression.h and it contains class gsGeometryMap<T>
                                                                   
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;               // in summary, the above four terms represent the classes defined in gsExpression.h

    A.setIntegrationElements(dbasis);                              // end up with gsExprHelper<T>::setMultiBasis()
                                                                   // mesh_ptr in gsExprHelper.h pointing to dbasis
    
    
    gsInfo << "\n";
    gsInfo << "############################################" << "\n";
    gsExprEvaluator<> ev(A);                                       // m_exprdata in gsExprEvaluator.h initialized by m_exprdata in gsExprAssembler.h
    
    gsInfo << "\n";
    geometryMap G = A.getMap(multipatch);                          // invoking getmap() in gsExprHelper.h 
                                                                   // dealing with mapVar, an object of gsGeometryMap<T>, defined in gsExpression.h, created in gsExprHelper.h
                                                                   // m_fs of mapVar pointing to multipatch
                                                                   // m_fd pointing to mapData
    
//     gsInfo << "\n";
//     gsInfo << "G.source(): \n" << G.source() << "\n";
//     gsInfo << "G.data().points: \n" << G.data().points << "\n";                 // data() returning mapData created in an object of gsExprHelper<T>
//     gsInfo << "G.data().values[0]: \n" << G.data().values[0].col(0) << "\n";
//     gsInfo << "meas(G): " << meas(G) << "\n";
    
    
    gsInfo << "\n";
    space u = A.getSpace(dbasis);                                           // invoking getSpace() in gsExprHelper.h
                                                                            // initializing m_ptable by dbasis
                                                                            // u links to dbasis by m_fs, to the second item of m_ptable by m_fd
    
                                                                            // In an object of gsFeVariable<T>, created when creating an object of gsFeSpace<T>
                                                                            // m_fs, a pointer to gsFunctionSet<T>, pointing to dbasis
                                                                            // m_fd, a pointer to gsFuncData<T>, pointing to m_ptable[&dbasis], the second item of m_ptable
                                                                            // m_d = dim
                                                                            // m_md pointing to NULL, see line 682 of gsExpression.h
    u.setInterfaceCont(0);
    u.addBc( bc.get("Dirichlet") );

    
    gsInfo << "\n";
    variable ff = A.getCoeff(f, G);                                         // initializing m_itable by f
                                                                            // in ff
                                                                            // m_fs pointing to f
                                                                            // m_fd pointing to m_itable[&f]
                                                                            // m_md pointing to mapData

    gsInfo << "\n";
    gsFunctionExpr<> ms;
    fd.getId(3, ms); // id=3: reference solution
    
    
    gsInfo << "\n";
    variable u_ex = ev.getVariable(ms, G);                                  // initializing m_itable by f
                                                                            // in u_ex
                                                                            // m_fs pointing to ms
                                                                            // m_fd pointing to m_itable[&ms]
                                                                            // m_md pointing to mapData
    
    gsInfo << "\n";
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

//     gsInfo << "source of the variable rhs: " << ff.source() << "\n";                    // note source() is a function defined in class gsFeVariable<T>
//     gsInfo<<"Exact solution: "<< ms << "\n";
    gsInfo << "source of the variable exact solution: " << u_ex.source() << "\n";
    
    gsSparseSolver<>::CGDiagonal solver;                                    // found in \src\gsMatrix\gsSparseSolver.h
    
# endif    
    

#if 1  
    //! [Solver loop]
    
    clock_t t;
    double time_cpu_per_refine;
    
    
    gsVector<> n_refine(numRefine+1);
    gsVector<> ndofs(numRefine+1);
    gsVector<> l2err(numRefine+1), h1_semierr(numRefine+1);
    gsVector<> time_cpu(numRefine+1);
//     gsInfo<< "(dot1=assembled, dot2=solved, dot3=got_error)\n";
    
    std::ofstream file;
    file.open("data_error_iga.txt",std::ios::app);
    file << "refine " << "ndofs " << "l2err " << "h1_semierr " << "\n";
    file.close();    
    
    for (int r=0; r<=numRefine; ++r)
    {
        gsInfo << "\n";     
        gsInfo << r << "th refinement\n";
        
        t = clock();
        
        n_refine[r] = r;
        
        
        if (r>0)
        {
            dbasis.uniformRefine();
        }
        
        gsInfo << "\n";
        gsInfo << "dbasis: \n";
        for (unsigned int i = 0; i<dbasis.nBases(); i++)
        {
            gsInfo << "[" << i << "] \n" << dbasis.basis(i) << "\n";             
        }        

//         gsInfo << "\n";
        A.initSystem();

        gsInfo << "DoFs: " << A.numDofs() <<std::flush;
        ndofs[r]=A.numDofs();

/*        gsInfo << "\n\n";
        gsInfo << "A.m_matrix:\n";
        gsInfo << A.matrix();
        gsInfo << "A.m_rhs:\n";
        gsInfo << A.rhs(); */       
        
        gsInfo << "\n\n";
        gsInfo << "////////////////////////////\n"
               << "assembling started\n";
               
//         gsInfo << "G.source(): \n";
//         gsInfo << G.source() << "\n";
//         gsInfo << "G.source().support(): \n";
//         gsInfo << G.source().support() << "\n";
//         gsInfo << "G.source().basis(): \n";
//         gsInfo << G.source().basis(0) << "\n";
//         gsInfo << "G.source().coefs(): " << G.source().coefs() << "\n";
        
//                
               
        gsInfo << "\n";
        gsInfo << "u.source(): \n";
        gsInfo << u.source();
        
        gsInfo << "u.source().support()\n";
        gsInfo << u.source().support();
        
        gsInfo << "\n";
        gsInfo << "u.data().actives: \n";
        gsInfo << u.data().actives;

        gsInfo << "\n";
        gsInfo << "u.data().values: \n";
        for(unsigned int i=0; i< u.data().values.size(); i++)
        {
            gsInfo << "[" << i << "]: \n";
            gsInfo << u.data().values[i] << "\n";
        }            
        
               
        gsInfo << "\n";
        gsInfo << "left-hand side and right-hand side\n";
        
        
//         gismo::expr::mult_expr<> expr_argu_first(igrad(u, G) * igrad(u, G).tr() * meas(G));
        
//         gsInfo << "typenames: \n";
//         gsInfo << "LHS: " << typeid(igrad(u, G) * igrad(u, G).tr() * meas(G)).name() << "\n";
//         gsInfo << "RHS: " << typeid(u * ff * meas(G)).name() << "\n";
        
//         gsInfo << "first argument: ";
//         (igrad(u, G) * igrad(u, G).tr() * meas(G)).print(gsInfo);
//         gsInfo << "second argument: ";
//         (u * ff * meas(G)).print(gsInfo);
        
        
        gsInfo << "\n";
        A.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G) );
        
        gsInfo << "\n";
        gsInfo << "A.m_matrix:\n";
        gsInfo << A.matrix();
        gsInfo << "A.m_rhs:\n";
        gsInfo << A.rhs();
        
//         gsInfo << "\n\n";
//         gsInfo << "G.data().points: \n" << G.data().points << "\n";
//         gsInfo << "info of G.data().values\n";
//         for(unsigned i=0; i<G.data().values.size(); i++)
//         {
//             gsInfo << "["<< i << "]: \n" << G.data().values[i] << "\n";
//         }
//         gsInfo << "meas(G): " << meas(G) << "\n";        
//         

        // Enforce Neumann conditions to right-hand side
        
#if 0
        gsInfo << "\n\n";
        gsInfo << "right-hand side and boundary\n";
        variable g_N = A.getBdrFunction();
        A.assembleRhsBc(u * g_N.val() * nv(G).norm(), bc.neumannSides() );
#endif
        
        gsInfo << "\n";
        gsInfo << "assembling finished\n"
               << "////////////////////////////\n";
               
//         gsInfo << "\n";
//         gsInfo<<"Sparse Matrix:\n"<< A.matrix().toDense() <<"\n";
//         gsInfo<<"Rhs vector:\n"<< A.rhs().transpose() <<"\n";
        
//         gsInfo << "A.matrix().printSparsity():\n";
//         gsInfo << A.matrix().printSparsity();
        
        if (is_matrix_stored!=0 && r>0)
        {
            std::ofstream fid;
            fid.open("matrix_deg_"+std::to_string(numElevate+1)+"_refine_"+std::to_string(r)+".txt", std::ios::trunc);
            for (index_t i = 0; i!=A.matrix().rows(); ++i)              // print A.matrix() element-wisely
            {
                for (index_t j = 0; j!=A.matrix().cols(); ++j)
                {
//                     gsInfo << A.matrix().coeff(i,j) << " ";
                    fid << A.matrix().coeff(i,j) << " ";
                }
//                 gsInfo << "\n";
                fid << "\n";
            }
            fid.close();
        }

               
#if 0
        gsInfo << "\n";
        gsInfo<< "solving" <<std::flush; // Linear solving done
        gsInfo << "\n";
        solver.compute( A.matrix() );
        solVector = solver.solve(A.rhs());

//         gsInfo << "solution reads:\n";
//         gsInfo << solVector << "\n\n";
        
        gsInfo << "computing errors\n";
        l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
        h1_semierr[r]= math::sqrt(ev.integral( ( igrad(u_ex) - grad(u_sol)*jac(G).inv() ).sqNorm() * meas(G) ));
        
        
#endif

        t = clock() - t;
        time_cpu_per_refine = ((float)t)/CLOCKS_PER_SEC;
        time_cpu[r] = time_cpu_per_refine;
//         printf ("It took me %ld clicks (%f seconds).\n",t,time_cpu_per_refine); 
                
        
        if (r>0)
        {
            file.open("data_error_iga.txt",std::ios::app);
            file << n_refine[r] << " " << ndofs[r] << " " << l2err[r] << " " << h1_semierr[r] << " " << time_cpu[r] << "\n";
            file.close();
        }

    }

    
#if 1
    gsInfo<< "\n\nn_refine: " << n_refine.transpose()<<"\n";
    gsInfo<< "ndofs: " << ndofs.transpose()<<"\n";
    gsInfo<< "L2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<h1_semierr.transpose()<<"\n";
    gsInfo<< "CPU time: "<<std::scientific<<time_cpu.transpose()<<"\n";

    if (!last && numRefine>0)
    {
        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              << ( l2err.head(numRefine).array() /
                   l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              <<( h1_semierr.head(numRefine).array() /
                  h1_semierr.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }

#endif     
    
    std::string out;
     
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
//         ev.options().setSwitch("plot.elements", true);
        if (id_dimension==1)
        {
            ev.writeParaview( u_sol, G, "poisson2_solution_");
        }
//         ev.writeParaview( u_ex    , G, "solution_ex");
        //ev.writeParaview( u, G, "aa");
//         gsFileManager::open("solution.pvd");
        
        
        const gsBasis<>& basis = dbasis.basis(0);
        
        gsInfo << "\n";
        out = "poisson2_basis_";
        gsInfo << "Writing the basis to a paraview file: " << out << "\n";
        gsWriteParaview(basis, out); 
        
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;
    
#endif

}// end main
