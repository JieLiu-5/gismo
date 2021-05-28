/** @file gsQuadrature.h

    @brief Creates a variety of quadrature rules

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsIO/gsOptionList.h>
#include <gsAssembler/gsGaussRule.h>
#include <gsAssembler/gsLobattoRule.h>

namespace gismo
{

/// Helper class for obtaining a quadrature rule
struct gsQuadrature
{
    /// Quadrature rule types
    enum rule
    {
        GaussLegendre = 1, ///< Gauss-Legendre quadrature
        GaussLobatto  = 2  ///< Gauss-Lobatto quadrature
    };

    /// Constructs a quadrature rule based on input \a options
    template<class T>
    static gsQuadRule<T> get(const gsBasis<T> & basis,
                             const gsOptionList & options, short_t fixDir = -1)
    {
        gsInfo << "gsQuadrature<T>::get() 1\n";
        
        const int     qu  = options.askInt("quRule", GaussLegendre);            // askInt() returning the value of GaussLegendre
        const T       quA = options.getReal("quA");
        const index_t quB = options.getInt ("quB");
        
        gsInfo << "    qu: " << qu
               << ", quA: " << quA
               << ", quB: " << quB
               << "\n";
        
        const gsVector<index_t> nnodes = numNodes(basis,quA,quB,fixDir);        // nnodes denoting nr. of quadrature points in each direction
                                                                                // it is related to the nr. of basis functions
        
        gsInfo << "    nnodes: " << nnodes << "\n";
        
        return get<T>(qu, nnodes);                  // defined below
    }

    /// Constructs a quadrature rule based on input \a options
    template<class T>
    static inline gsQuadRule<T> get(int qu, gsVector<index_t> const & numNodes, unsigned digits = 0)
    {
        
        gsInfo << "gsQuadrature<T>::get() 2\n";
        
//         gsInfo << "qu: " << qu << "\n"
//                << "numNodes: \n" << numNodes << "\n"
//                << "digits: " << digits << "\n"
//                << "\n";
        
        switch (qu)
        {
        case GaussLegendre :
//             gsInfo << "GaussLegendre\n";
            return gsGaussRule<T>(numNodes, digits);
        case GaussLobatto :
            return gsLobattoRule<T>(numNodes, digits);
        default:
            GISMO_ERROR("Invalid Quadrature rule request ("<<qu<<")");
        };
    }

    /// Computes and integer quA*deg_i + quB where deg_i is the degree
    /// of \a basis
    template<class T>
    static gsVector<index_t> numNodes(const gsBasis<T> & basis,
                               const T quA, const int quB, short_t fixDir = -1)
    {
        
        gsInfo << "gsQuadrature<T>::numNodes()\n";
        
        const short_t d  = basis.dim();
        GISMO_ASSERT( fixDir < d && fixDir>-2, "Invalid input fixDir = "<<fixDir);
        gsVector<index_t> nnodes(d);

        if (-1==fixDir)
        {
            fixDir = d;
        }
        else
        {
            gsInfo << "fixDir: " << fixDir << "\n";
            nnodes[fixDir] = 1;
        }

        short_t i;
        for(i=0; i!=fixDir; ++i )
        {
            //note: +0.5 for rounding
//             gsInfo << "first for loop accessed\n";
            nnodes[i] = cast<T,index_t>(quA * basis.degree(i) + quB + 0.5);
        }
        for(++i; i<d; ++i )
        {
//             gsInfo << "second for loop accessed\n";
            nnodes[i] = cast<T,index_t>(quA * basis.degree(i) + quB + 0.5);
        }
        return nnodes;
    }
};

}// namespace gismo
