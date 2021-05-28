/** @file gsVisitorPoisson.h

    @brief Poisson equation element visitor.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsAssembler/gsQuadrature.h>

namespace gismo
{

/** \brief Visitor for the Poisson equation.
 *
 * Assembles the bilinear terms
 * \f[ (\nabla u,\nabla v)_\Omega \text{ and } (f,v)_\Omega \f]
 * For \f[ u = g \quad on \quad \partial \Omega \f],
 *
 */

template <class T, bool paramCoef = false>
class gsVisitorPoisson
{
public:

    /** \brief Constructor for gsVisitorPoisson.
     */
    gsVisitorPoisson(const gsPde<T> & pde)
    { 
        gsInfo << "\n";
//         gsInfo << "*******************";          
        std::cout << "gsVisitorPoisson::gsVisitorPoisson()\n"
                  << "    creating members pde_ptr, rhs_ptr, md, actives, basisData, rhsVals\n";            // md: gsMapData, see below
        
        pde_ptr = static_cast<const gsPoissonPde<T>*>(&pde);            // pde_ptr points to the argument received
    }
    
    void initialize(const gsBasis<T> & basis,                                                                           // gsVisitorPoisson::initialize() definition
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T>    & rule)
    {
        gsInfo << "\n";
        gsInfo << "*******************\n";            
        std::cout << "gsVisitorPoisson<T>::initialize()\n"
                  << "    setting up rhs and quadrature\n\n";
        
        // Grab right-hand side for current patch
        
        gsInfo << "initializing rhs_ptr by pde_ptr and patchIndex\n";
        rhs_ptr = &pde_ptr->rhs()->piece(patchIndex);                                   // rhs_ptr: a pointer to gsFunction<T>, defined below, used in evaluate()
        
//         gsInfo << "  rhs reads: ";
//         rhs_ptr->print(gsInfo);
//         gsInfo << "\n";
//         
        // Setup Quadrature
        
        
        gsInfo << "\n";
        gsInfo << "setting up quRule using basis and options\n";
        
        rule = gsQuadrature::get(basis, options);                                       // harmless slicing occurs here       
                                                                                        // gsQuadrature.h is included in this file
        
        
        gsInfo << "\n";
        gsInfo << "setting flags for md\n";
        
        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;                     // md, an object of gsMapData<T>, declared below    
        
    }

    // Evaluate on element.
    inline void evaluate(const gsBasis<T>       & basis,                                                                // gsVisitorPoisson::evaluate() definition
                         const gsGeometry<T>    & geo,
                         const gsMatrix<T>      & quNodes)
    {
        gsInfo << "\n";
        gsInfo << "*******************\n";            
        std::cout << "gsVisitorPoisson<T>::evaluate()" << std::endl;
        
        gsInfo << "\n";
        
        md.points = quNodes;
        
        gsInfo << "  md.points: " << "\n"
               << md.points << "\n";
        
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(md.points.col(0), actives);                           // actives is a member in this class
        
        numActive = actives.rows();

/*        gsInfo << "  active basis functions(" << numActive << " in total): ";
        for (int i = 0; i<numActive; ++i)
        {
            gsInfo << actives.at(i) << " ";
        }
        gsInfo << "\n";  */       
        
        gsInfo << "actives(" << numActive << " in total): \n";         
        gsInfo << actives.asVector();  
        gsInfo << "\n";
        
        
        // Evaluate basis functions on element
        
/*        gsInfo << "\n";             
        gsInfo << "evaluating values and gradients of all basis functions at quadrature points \n";
        gsInfo << "\n";   */          
//         gsInfo << "comprising multiple rows, each row for each basis function, each column for each quadrature point\n";
        
        basis.evalAllDers_into( md.points, 1, basisData);       // definitions of evalAllDers_into() can be found in gsFunctionSet.hpp, gsBasis.hpp, but neither of them are used here, Jan. 30, 2020
                                                                // basisData containing values, i.e. basisData[0], and gradients, i.e. basisData[1], of the related basis function on the quadrature points
        
        gsInfo << "  size of basisData: " << basisData.size() << "\n";
        gsInfo << "  size of basisData[0]: " << basisData[0].size()
               << ", n_rows: " << basisData[0].rows()
               << ", n_cols: " << basisData[0].cols()
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
        
        
//         gsInfo << "  md.values[0]: ";
//         gsInfo << md.values[0];
//         gsInfo << "\n";
        
//         gsInfo << "md.measures: " << md.measures << "\n";                             // values not available before computeMap()   
        
        
//         gsInfo << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";        
//         
        gsInfo << "\n";
        gsInfo << "md.maxDeriv(): " << md.maxDeriv() << "\n";                         // used in gsFunctionSet<T>::compute()
//         
        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geo.computeMap(md);
        
        
        
//         gsInfo << "\n";
//         gsInfo << "  md.flags: " << md.flags << "\n";
//         gsInfo << "  md.deriv2Size(): " << md.deriv2Size() << "\n";
//         gsInfo << "  md.derivSize(): " << md.derivSize() << "\n";
        
        
/*        gsInfo << "\n";
        gsInfo << "md.values(" << md.values.size() << " in total): \n"                // not clear how these data come, mar 12, 2020
               << "[0]: \n" << md.values[0] << "\n"                                  
               << "[1]: \n" << md.values[1] << "\n";   */        
               
               
        gsInfo << "\n";
        gsInfo << "md.measures: " << md.measures << "\n";             // measures, defined in gsFuncData.h, is determined by coefs, 
        gsInfo << "\n";                                               // and hence, measures for each cell are different
                                                                        // used in assemble() to determine weight
                                                                        // we set measure to be 1 for each quadrature point
        
        
//         gsInfo << "  md.dim.first: " << md.dim.first                    // used in transformGradients()
//                << ", second: " << md.dim.second
//                << "\n";
//         gsInfo << "  md.jacobians(): \n"                                // definition of jacobians not clear
//                << md.jacobians() << "\n";                               // used in transformGradients() in assemble()
//         
/*        gsInfo << "  md.jacobian(0): " << md.jacobian(0)                
               << ", (1): " << md.jacobian(1) << "\n";   */             
//         gsInfo << "\n";             
//         gsInfo << "  size of md.values: " << md.values.size() << "\n";
        
//         gsInfo << "\n";             
//         gsInfo << "  contribution of basis function and coefs at quadrature points\n";
//         
                  
        // Evaluate right-hand side at the geometry points paramCoef
        // specifies whether the right hand side function should be
        // evaluated in parametric(true) or physical (false)
        
//         
//         gsInfo << "  rhs reads: ";
//         rhs_ptr->print(gsInfo);
//         gsInfo << "\n";
        
/*        gsInfo << "  md.points: \n"
               << md.points
               << "\n";  */      

        gsInfo << "\n";
        gsInfo << "paramCoef: " << paramCoef << "\n";

        rhs_ptr->eval_into( (paramCoef ?  md.points :  md.values[0] ), rhsVals );                   
        
/*        gsInfo << "  rhsVals: \n";
        gsInfo << rhsVals;
        gsInfo << "\n";      */  
        
//         gsInfo << "  rhs_ptr->eval(md.points): \n"
//                << rhs_ptr->eval( md.points)
//                << "\n";
        
//         rhs_ptr->eval_into( md.points, rhsVals );                               // usage of values not clear, we use rhsVals directly
                                                                                // location of eval_into() not clear
                                                                                // (paramCoef ?  md.points :  md.values[0] )            basisData[0].row(0)
        
        gsInfo << "rhsVals: \n"
               << rhsVals
               << "\n";      
        
        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive      );
        localRhs.setZero(numActive, rhsVals.rows() );//multiple right-hand sides
    }
    
    inline void assemble(gsDomainIterator<T>    & ,                                                                         // gsVisitorPoisson::assemble() definition
                         gsVector<T> const      & quWeights)
    {
        gsInfo << "\n";
        gsInfo << "*******************\n";            
        std::cout << "gsVisitorPoisson<T>::assemble()" << std::endl;
        
        gsMatrix<T> & bVals  = basisData[0];                                // values of basis functions at quadrature points
        gsMatrix<T> & bGrads = basisData[1];                                // gradients of basis functions at quadrature points

//         gsInfo << "  bVals: \n"
//                << bVals
//                << "\n";               
//                
//         gsInfo << "  bGrads: \n"
//                << bGrads
//                << "\n";
        
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * md.measure(k);    
            
            
            gsInfo << "\n";   
            gsInfo << k << "th quadrature point\n";
            
            gsInfo << "    measure: " << md.measure(k)                        // measure is valid after computeMap()
                   << ", weight: " << weight
                   << "\n";
            
            // Compute physical gradients at k as a Dim x NumActive matrix
            transformGradients(md, k, bGrads, physGrad);                    // a function defined in gsAssembler.h
                                                                            // physGrad equals bGrads divided by md.measure(k), evaluated by each basis function
            
                   
            gsInfo << "    physical gradients: " << physGrad
                   << "\n";               
            
            localRhs.noalias() += weight * ( bVals.col(k) * rhsVals.col(k).transpose() ) ;                          // definition of noalias(), col() can be found in library Eigen
            localMat.noalias() += weight * (physGrad.transpose() * physGrad);
            
//             gsInfo << "\n";            
//             gsInfo << "  localRhs: \n";                                           
//             gsInfo << localRhs;
//             gsInfo << "\n";  
//             
//             gsInfo << "  localMat: \n";                                           
//             gsInfo << localMat;
//             gsInfo << "\n";  
            
        }
        
/*        gsInfo << "\n";
        gsInfo << "############################################" << "\n"; 
        gsInfo << "assemble() results in\n";*/   
        
        gsInfo << "\n";
        gsInfo << "localRhs: \n"
               << localRhs
               << "\n";  
        gsInfo << "localMat: \n"                                          
               << localMat
               << "\n";  
               
        
    }

    inline void localToGlobal(const int patchIndex,                                                                         // gsVisitorPoisson::localToGlobal() definition
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        gsInfo << "\n";
        gsInfo << "*******************\n";            
        std::cout << "gsVisitorPoisson<T>::localToGlobal()" << std::endl;
        
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);
        
        gsInfo << "\n";
        gsInfo << "actives globally(" << numActive << " in total): \n";         
        gsInfo << actives.asVector();  
        gsInfo << "\n";        

        gsInfo << "eliminatedDofs.front(): \n" << eliminatedDofs.front() << "\n";
        
        // Add contributions to the system matrix and right-hand side
        system.push(localMat, localRhs, actives, eliminatedDofs.front(), 0, 0);
        
        gsInfo << "\n";
        gsInfo << "system.matrix():\n"
               << system.matrix();
        gsInfo << "system.rhs():\n"
               << system.rhs()
               << "\n";               
        
    }

protected:
    // Pointer to the pde data
    const gsPoissonPde<T> * pde_ptr;
    
protected:
    // Basis values
    std::vector<gsMatrix<T> > basisData;
    gsMatrix<T>        physGrad;
    gsMatrix<unsigned> actives;
    index_t numActive;

protected:
    // Right hand side ptr for current patch
    const gsFunction<T> * rhs_ptr;

    // Local values of the right hand side
    gsMatrix<T> rhsVals;

protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

    gsMapData<T> md;                                                                        // a class defined in gsFuncData.h
};


} // namespace gismo

