 /** @file gsPoissonAssembler.h

    @brief Provides assembler for the Poisson equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss A. Mantzaflaris, J. Sogn
*/

#pragma once

#include <gsAssembler/gsAssembler.h>
#include <gsPde/gsPoissonPde.h>


namespace gismo
{

/** @brief
    Implementation of an (multiple right-hand side) Poisson assembler.

    The Poisson equation: \f$-\Delta\mathbf{u}=\mathbf{f} \f$

    It sets up an assembler and assembles the system patch wise and combines
    the patch-local stiffness matrices into a global system by various methods
    (see gismo::gsInterfaceStrategy). It can also enforce Dirichlet boundary
    conditions in various ways (see gismo::gsDirichletStrategy).
    
    \ingroup Assembler
*/
template <class T>
class gsPoissonAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:

    gsPoissonAssembler()
    {
        gsInfo << "gsPoissonAssembler constructor default" << "\n";
    }

    /** @brief Main Constructor of the assembler object.

    \param[in] pde A boundary value Poisson problem
    \param[in] bases a multi-basis that contains patch-wise bases
    */
    gsPoissonAssembler( const gsPoissonPde<T>          & pde,
                        const gsMultiBasis<T>          & bases)
    {
        Base::initialize(pde, bases, m_options);
    }

    
    /** @brief Main Constructor of the assembler object.

    \param[in] pde A boundary value Poisson problem
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] dirStrategy option for the treatment of Dirichlet boundary
    \param[in] intStrategy option for the treatment of patch interfaces
    */
    gsPoissonAssembler( const gsPoissonPde<T>          & pde,
                        const gsMultiBasis<T>          & bases,
                        dirichlet::strategy           dirStrategy,
                        iFace::strategy               intStrategy = iFace::glue)
    {
        gsInfo << "gsPoissonAssembler<T>::gsPoissonAssembler()\n"
               << "  using members upon creating an object of class gsAssembler<T>\n"
               << "  receiving arguments pde, bases, dirStrategy and intStrategy\n" ;
        
        gsInfo << "\n";
        
        
        gsInfo << "adjusting m_options by dirStrategy and intStrategy\n";
        gsInfo << "\n";
        
        m_options.setInt("DirichletStrategy", dirStrategy);     // m_options is a member in gsAssembler.h, Jan 23 2020      
                                                                // id of DirichletStrategy is 11 for dirichlet::elimination, and 12 for dirichlet::nitsche
        
        
        
        m_options.setInt("InterfaceStrategy", intStrategy);     // we initialize them here because arguments dirStrategy and intStrategy belong to gsPoissonAssembler<T>
        
//         gsInfo << "m_options reads:\n";
//         m_options.print(gsInfo);
//         gsInfo << "\n";
        
        
        Base::initialize(pde, bases, m_options);                // Base is an alias of gsAssembler<T>, this page upper          
                                                                // m_options is assigned to itselft
        
        gsInfo << "\n";
        gsInfo << "\n";
        
        gsInfo << "  m_system.matrix():\n"
               << m_system.matrix()
               << "\n";
        gsInfo << "  m_system.rhs():\n"
               << m_system.rhs()
               << "\n";          
        
//         gsInfo << "gsPoissonAssembler::gsPoissonAssembler() finished\n";
        
    }

    /** @brief
    Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] basis a multi-basis that contains patch-wise bases
    \param[in] bconditions is a gsBoundaryConditions object that holds all boundary conditions.
    \param[in] rhs is the right-hand side of the Poisson equation, \f$\mathbf{f}\f$.
    \param[in] dirStrategy option for the treatment of Dirichlet boundary
    \param[in] intStrategy option for the treatment of patch interfaces
    
    \ingroup Assembler
    */
    gsPoissonAssembler( gsMultiPatch<T> const         & patches,
                        gsMultiBasis<T> const         & basis,
                        gsBoundaryConditions<T> const & bconditions,
                        const gsFunction<T>           & rhs,
                        dirichlet::strategy           dirStrategy = dirichlet::elimination,
                        iFace::strategy               intStrategy = iFace::glue)
    {
        
//         gsInfo << "gsPoissonAssembler constructor 1" << "\n";
        
        m_options.setInt("DirichletStrategy", dirStrategy);
        m_options.setInt("InterfaceStrategy", intStrategy);

        typename gsPde<T>::Ptr pde( new gsPoissonPde<T>(patches,bconditions,rhs) );
        
        
        
        Base::initialize(pde, basis, m_options);
    }

    virtual gsAssembler<T>* clone() const
    {
        return new gsPoissonAssembler<T>(*this);
    }
    
    virtual gsAssembler<T>* create() const
    {
        return new gsPoissonAssembler<T>();
    }

    // Refresh routine
    virtual void refresh();

    // Main assembly routine
    virtual void assemble();
    
    /// Returns an expression of the "full" assembled sparse
    /// matrix. Note that matrix() might return a lower diagonal
    /// matrix, if we exploit possible symmetry during assembly
    /// (check: m_matrix.symmetry() == true )
    Eigen::SparseSelfAdjointView< typename gsSparseMatrix<T>::Base, Lower> fullMatrix()
    {
        return m_system.matrix().template selfadjointView<Lower>();
    }


protected:

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;
};


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPoissonAssembler.hpp)
#endif
