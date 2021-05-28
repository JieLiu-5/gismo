/** @file gsPoissonAssembler.hpp

    @brief Provides assembler implementation for the Poisson equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss, A. Mantzaflaris, J. Sogn
*/


#include <gsAssembler/gsVisitorPoisson.h> // Stiffness volume integrals
#include <gsAssembler/gsVisitorNeumann.h> // Neumann boundary integrals
#include <gsAssembler/gsVisitorNitsche.h> // Nitsche boundary integrals
#include <gsAssembler/gsVisitorDg.h>      // DG interface integrals

namespace gismo
{

template<class T>
void gsPoissonAssembler<T>::refresh()
{
    // We use predefined helper which initializes the system matrix
    // rows and columns using the same test and trial space
    Base::scalarProblemGalerkinRefresh();
}

template<class T>
void gsPoissonAssembler<T>::assemble()
{
    
    std::cout << "gsPoissonAssembler::assemble()" << std::endl;
    gsInfo << "\n";
    
    GISMO_ASSERT(m_system.initialized(), 
                 "Sparse system is not initialized, call initialize() or refresh()");

    // Reserve sparse system
    m_system.reserve(m_bases[0], m_options, this->pde().numRhs());
    
    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    Base::computeDirichletDofs();
    
    gsInfo << "\n";
    gsInfo << "m_ddof:\n";
    for (unsigned int i =0; i<m_ddof.size(); ++i)
    {
        gsInfo << i << ": \n" << m_ddof[i] << "\n";               
    }    

    // Clean the sparse system
   // m_system.setZero(); //<< this call leads to a quite significant performance degrade!
    
//     gsInfo << m_system.matrix();
       
    gsInfo << "\n";
    gsInfo << "########################################################################################" << "\n";
    gsInfo << "assembling initial LHS and RHS\n";
    gsInfo << "\n";
    // Assemble volume integrals
    Base::template push<gsVisitorPoisson<T> >();                // push(), a function of class gsAssembler defined in gsAssembler.h, invoking apply()
                                                                // Base, an alias of gsAssembler<T>, see gsPoissonAssembler.h; 
                                                                // gsVisitorPoisson.h is included in this file, gsAssembler.h is included in gsPoissonAssembler.h
                                                                // gsVisitorPoisson<T> is used to create an object for each patch
                                                                
                                                                
    
    
    gsInfo << "\n";
    gsInfo << "########################################################################################" << "\n";     
    gsInfo << "enforcing boundary conditions\n";
    gsInfo << "\n";
    
#if 1
    
//     gsInfo << "enforcing Neumann boundary conditions\n";
//     gsInfo << "\n";
    
    
    // Enforce Neumann boundary conditions
    
    if (m_pde_ptr->bc().neumannSides().size()!=0)
    {
        Base::template push<gsVisitorNeumann<T> >(m_pde_ptr->bc().neumannSides() );         // m_pde_ptr, a member defined in gsAssembler.h, pointing to an object of gsPde<T>
    } else                                                       // bc(), a function defined in gsPde<T>, returning an object of gsBoundaryConditions<T>, i.e. the boundary conditions being imposed
    {
        gsInfo << "    no Neumann boundary found\n";
    }
                                                                 // neumannSides(), a function defined in gsBoundaryConditions<t>, returning a reference to data of bcContainer type
//     gsInfo << "\n";
    
    
//     gsInfo << "done for enforcing Neumann boundary conditions\n";
    
//     gsInfo << "\n";
//     gsInfo << "enforcing Dirichlet boundary conditions\n";
//     gsInfo << "\n";
    
    const int dirStr = m_options.getInt("DirichletStrategy");
    
//     gsInfo << "DirichletStrategy: " << dirStr << "\n";
    
    
    // If requested, enforce Dirichlet boundary conditions by Nitsche's method
    if ( dirStr == dirichlet::nitsche )
    {
//         gsInfo << "  dirStr == dirichlet::nitsche\n\n";
        gsInfo << "    enforcing Dirichlet boundary conditions by Nitsche's method\n\n";
        
        Base::template push<gsVisitorNitsche<T> >(m_pde_ptr->bc().dirichletSides());            // gsVisitorNitsche<T> is used to create an object for each boundary
    }

     // If requested, enforce Dirichlet boundary conditions by diagonal penalization
     else if ( dirStr == dirichlet::penalize )
     {
//          gsInfo << "  dirStr == dirichlet::penalize\n";
        gsInfo << "    enforcing Dirichlet boundary conditions by diagonal penalization\n\n";
        
         Base::penalizeDirichletDofs();
     }else
     {
         gsInfo << "    enforcing Dirichlet boundary conditions by elimination\n\n";
     }

    if ( m_options.getInt("InterfaceStrategy") == iFace::dg )
    {
        gsInfo << "    InterfaceStrategy==iFace::dg\n\n";
        
        Base::template pushInterface<gsVisitorDg<T> >();
    }
    
    // Assembly is done, compress the matrix
    Base::finalize();
#endif
}


}// namespace gismo
