/** @file gsVisitorNitsche.h

    @brief Weak (Nitsche-type) BC imposition visitor for elliptic problems.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris , S. Moore
*/

#pragma once

namespace gismo
{
/** \brief Visitor for the weak imposition of the dirichlet boundary condition.
 *
 * Adds this term to the bilinear terms
 * \f[ (\nabla u, v)_{\partial \Omega} + (u, \nabla v )_{\partial \Omega} 
 *                                     + (\mu*u, v)_{\partial \Omega} \f]
 * 
 * The following term is also added to the linear form
 * \f[ (g_D, \mu*v + \nabla v)_{\partial \Omega} \f],
 * where the dirichlet term is given as \f[ g_D \f].
 */

template <class T>
class gsVisitorNitsche
{
public:

    gsVisitorNitsche(const gsPde<T> & , const boundary_condition<T> & s)
    : dirdata_ptr( s.function().get() ), side(s.side())
    { 
        gsInfo << "gsVisitorNitsche<T>::gsVisitorNitsche()\n"
               << "    creating dirdata_ptr and side which are initialized by the arguments\n"
               << "    creating members md, actives, basisData, dirdata, unormal, pGrads\n\n";
    }

/** @brief
    Constructor of the assembler object.

    \param[in] dirdata  is a gsBoundaryConditions object that holds boundary conditions of the form:
    \f[ \text{Dirichlet: } u = g_D \text{ on } \Gamma.\f]
    \f$ v \f$ is the test function and \f$ \Gamma\f$ is the boundary side.
    \param[in] _penalty for inputing the penalty choice
    \param[in] s
*/
    gsVisitorNitsche(const gsFunction<T> & dirdata, T _penalty, boxSide s) : 
    dirdata_ptr(&dirdata),penalty(_penalty), side(s)
    { }

    void initialize(const gsBasis<T> & basis, 
                    gsQuadRule<T> & rule)
    {
        const int dir = side.direction();
        gsVector<int> numQuadNodes ( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i)
            numQuadNodes[i] = 2* basis.degree(i) + 1;
        numQuadNodes[dir] = 1;
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE|NEED_JACOBIAN|NEED_GRAD_TRANSFORM;
    }

    void initialize(const gsBasis<T> & basis,
                    const index_t ,
                    const gsOptionList & options, 
                    gsQuadRule<T>    & rule)
    {
        
        gsInfo << "\n";
        gsInfo << "*******************\n";            
        std::cout << "gsVisitorNitsche<T>::initialize()\n"
                  << "    setting up quadrature\n\n";
                  
        gsInfo << "side.direction(): " << side.direction() << "\n";
                  
        // Setup Quadrature (harmless slicing occurs)
        rule = gsQuadrature::get(basis, options, side.direction());

        gsInfo << "\n";
        gsInfo << "setting flags for md\n";
        
        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;

        // Compute penalty parameter
        const int deg = basis.maxDegree();
        penalty = (deg + basis.dim()) * (deg + 1) * T(2.5);
        
        gsInfo << "penalty parameter: " << penalty << "\n";
        
    }

    // Evaluate on element.
    inline void evaluate(const gsBasis<T>       & basis, // to do: more unknowns
                         const gsGeometry<T>    & geo,
                         // todo: add element here for efficiency
                         const gsMatrix<T>      & quNodes)
    {
        
        gsInfo << "\n";
        gsInfo << "*******************\n";            
        std::cout << "gsVisitorNitsche<T>::evaluate()\n\n";
//                   << "    setting up quadrature\n\n";    
                  
        md.points = quNodes;
        
        gsInfo << "  md.points: " << "\n"
               << md.points << "\n";
               
               
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the current element
        basis.active_into(md.points.col(0) , actives);
                
        const index_t numActive = actives.rows();

        gsInfo << "actives(" << numActive << " in total): \n";         
        gsInfo << actives.asVector();  
        gsInfo << "\n";

        // Evaluate basis values and derivatives on element
        basis.evalAllDers_into( md.points, 1, basisData);
        
        
        gsInfo << "\n";
        gsInfo << "  size of basisData: " << basisData.size() << "\n";
        gsInfo << "  size of basisData[0]: " << basisData[0].size()
               << ", nrows: " << basisData[0].rows()
               << ", ncols: " << basisData[0].cols()
               << "\n";
               
        gsInfo << "\n";
               
        gsInfo << "  basisData[0](values): \n";  
        for (auto i=0; i<basisData[0].rows(); ++i)
        {
            gsInfo << "  " << i << "th basis: " << basisData[0].row(i) << "\n";
            
        }
        
        gsInfo << "\n";
        gsInfo << "  basisData[1](gradients): \n";
        for (auto i=0; i<basisData[1].rows(); ++i)
        {
            gsInfo << "  " << i << "th basis: " << basisData[1].row(i) << "\n";
            
        }
        
        gsInfo << "\n";
        

        // Compute geometry related values
        geo.computeMap(md);

        // Evaluate the Dirichlet data
        dirdata_ptr->eval_into(md.values[0], dirData);
        
        gsInfo << "\n";
        gsInfo << "  dirData: \n"
               << dirData
               << "\n";        

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive);
        localRhs.setZero(numActive, dirdata_ptr->targetDim() );
    }

    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {
        gsInfo << "\n";
        gsInfo << "*******************\n";            
        std::cout << "gsVisitorNitsche<T>::assemble()\n\n";
        
        gsMatrix<T> & bGrads = basisData[1];
        const index_t numActive = actives.rows();
        
        gsInfo << "bGrads: \n"
               << bGrads
               << "\n";        

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            
        gsInfo << "\n";
        gsInfo << k << "th quadrature" << "\n\n";              
        
        const typename gsMatrix<T>::Block bVals =
            basisData[0].block(0,k,numActive,1);
            
        gsInfo << "bVals: \n"
               << bVals
               << "\n";                

        // Compute the outer normal vector on the side
        outerNormal(md, k, side, unormal);
        
        gsInfo << "unormal: \n" << unormal << "\n"
               << "unormal.norm(): " << unormal.norm()
               << "\n\n";         

        // Multiply quadrature weight by the geometry measure
        const T weight = quWeights[k] *unormal.norm();   

        
        gsInfo << "weight: " << weight
               << "\n";  
               
        // Compute the unit normal vector 
        unormal.normalize();
        
        gsInfo << "unormal: \n" << unormal << "\n\n";           
        
        // Compute physical gradients at k as a Dim x NumActive matrix
        transformGradients(md, k, bGrads, pGrads);
        
        
        gsInfo << "\n";
        gsInfo << "physical gradients: \n"
               << pGrads
               << "\n";          
        
        // Get penalty parameter
        const T h = element.getCellSize();
        const T mu = penalty / (0!=h?h:1);
        
        gsInfo << "h: " << h << ", mu: " << mu << "\n";

        // Sum up quadrature point evaluations
        localRhs.noalias() -= weight * (( pGrads.transpose() * unormal - mu * bVals )
                                        * dirData.col(k).transpose() );

        localMat.noalias() -= weight * ( bVals * unormal.transpose() * pGrads
                           +  (bVals * unormal.transpose() * pGrads).transpose()
                           -  mu * bVals * bVals.transpose() );
        }
        
        gsInfo << "\n";
        gsInfo << "localRhs: \n"
               << localRhs
               << "\n";  
        
        gsInfo << "localMat: \n"                                          
               << localMat
               << "\n";  

    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> >    & ,
                              gsSparseSystem<T>     & system)
    {
        gsInfo << "\n";
        gsInfo << "*******************\n";            
        std::cout << "gsVisitorNitsche<T>::localToGlobal()\n\n";
        
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);

        // Add contributions to the system matrix and right-hand side
        system.pushAllFree(localMat, localRhs, actives, 0);
    }
    
    void localToGlobal(const gsDofMapper  & mapper,
                       const gsMatrix<T>     & eliminatedDofs,
                       const int patchIndex,
                       gsSparseMatrix<T>     & sysMatrix,
                       gsMatrix<T>           & rhsMatrix )
    {
        // Local DoFs to global DoFs
        mapper.localToGlobal(actives, patchIndex, actives);
        const index_t numActive = actives.rows();

        // Push element contribution to the global matrix and load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            const unsigned jj = actives(j);
            rhsMatrix.row(jj) += localRhs.row(j);
            for (index_t i=0; i!=numActive; ++i)
            {
                const unsigned ii = actives(i);
//                if ( jj <= ii ) // assuming symmetric problem
                    sysMatrix( ii, jj ) += localMat(i,j);
            }
        }

    }

private:
    // Dirichlet function
    const gsFunction<T> * dirdata_ptr;

    // Penalty constant
    T penalty;

    // Side
    boxSide side;

private:
    // Basis values
    std::vector<gsMatrix<T> > basisData;
    gsMatrix<T>      pGrads;
    gsMatrix<unsigned> actives;

    // Normal and Neumann values
    gsVector<T> unormal;
    gsMatrix<T> dirData;

    // Local  matrix and rhs
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

    gsMapData<T> md;
};


} // namespace gismo
