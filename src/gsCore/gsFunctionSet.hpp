/** @file gsFunctionSet.hpp

    @brief implementation of default functions of the gsFunctionSet

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#include <gsCore/gsFuncData.h>
#include <gsCore/gsFunction.h>
#include <gsCore/gsBasis.h>

namespace gismo                                             // not much compiling time
{

template <class T>
gsFunctionSet<T>::gsFunctionSet() {}

template <class T>
gsFunctionSet<T>::gsFunctionSet(const gsFunctionSet &) {}

template <class T>
gsFunctionSet<T>::~gsFunctionSet () {}

template <class T>
const gsFunction<T>& gsFunctionSet<T>::function(const index_t k) const
{
    gsInfo << "gsFunctionSet<T>::function()\n";
    
    GISMO_ASSERT(dynamic_cast<const gsFunction<T>*>(&piece(k)),
                 "No function found, instead: "<< piece(k));
    
    return static_cast<const gsFunction<T>&>(piece(k));                         // piece() is a function in class gsMultipatch<T>
}

template <class T>
const gsBasis<T>& gsFunctionSet<T>::basis(const index_t k) const
{
    GISMO_ASSERT(dynamic_cast<const gsBasis<T>*>(&piece(k)),
                 "No basis found, instead: "<< piece(k));
    return static_cast<const gsBasis<T>&>(piece(k));                            
}

// support (domain of definition)
template <class T>
gsMatrix<T> gsFunctionSet<T>::support() const
{
    return gsMatrix<T>();
}

// actives

template <typename T>
void gsFunctionSet<T>::active_into (const gsMatrix<T> &, gsMatrix<unsigned> &) const
{
    GISMO_NO_IMPLEMENTATION
    // Single function 0 globally active:
    // result.setConstant(1,u.cols(),0);
}

// evaluation

template <typename T>
void gsFunctionSet<T>::eval_into (const gsMatrix<T> &, gsMatrix<T> &) const
{
//     gsInfo << "gsFunctionSet<T>::eval_into()\n";
    
    GISMO_NO_IMPLEMENTATION
}

template <typename T>
void gsFunctionSet<T>::deriv_into (const gsMatrix<T> &, gsMatrix<T> &) const
{GISMO_NO_IMPLEMENTATION}

template <typename T>
void gsFunctionSet<T>::deriv2_into (const gsMatrix<T> &, gsMatrix<T> &) const
{GISMO_NO_IMPLEMENTATION}

template <typename T>
void gsFunctionSet<T>::evalAllDers_into(const gsMatrix<T> & u, const int n,
                                        std::vector<gsMatrix<T> > & result) const
{
    
    gsInfo << "gsFunctionSet<T>::evalAllDers_into()\n";
    
    result.resize(n+1);

    switch(n)
    {
    case 0:
//         gsInfo << "value computed\n";
        eval_into(u, result[0]);
        break;
    case 1:
//         gsInfo << "value and gradient computed\n";
        eval_into (u, result[0]);
        deriv_into(u, result[1]);
        break;
    case 2:
        eval_into  (u, result[0]);
        deriv_into (u, result[1]);
        deriv2_into(u, result[2]);
        break;
    default:
        GISMO_ERROR("evalAllDers implemented for order up to 2<"<<n ); //<< " for "<<*this);
        break;
    }
}

template <class T>
gsMatrix<T>
gsFunctionSet<T>::eval(const gsMatrix<T>& u) const
{
//     gsInfo << "gsFunctionSet<T>::eval()\n";
    
    gsMatrix<T> result;
    this->eval_into( u, result );
    return result;
}

template <class T>
gsMatrix<T>
gsFunctionSet<T>::deriv(const gsMatrix<T>& u) const
{
    gsInfo << "gsFunctionSet<T>::deriv()\n";
    
    gsMatrix<T> result;
    this->deriv_into( u, result );
    return result;
}

template <class T>
gsMatrix<T>
gsFunctionSet<T>::deriv2(const gsMatrix<T>& u) const
{
    gsInfo << "gsFunctionSet<T>::deriv2()\n";
    
    gsMatrix<T> result;
    this->deriv2_into( u, result );
    return result;
}

/*
template <typename T>
void gsFunctionSet<T>::div_into       (const gsMatrix<T> & u, gsMatrix<T> &result) const
{
    gsMatrix<T> tmp;
    deriv_into(u,tmp);
    convertValue<T>::derivToDiv(tmp, result, info());
}

template <typename T>
void gsFunctionSet<T>::curl_into      (const gsMatrix<T> & u, gsMatrix<T> &result) const
{
    gsMatrix<T> tmp;
    deriv_into(u,tmp);
    convertValue<T>::derivToCurl(tmp, result, info());
}

template <typename T>
void gsFunctionSet<T>::laplacian_into (const gsMatrix<T> & u, gsMatrix<T> &result) const
{
    gsMatrix<T> tmp;
    deriv2_into(u,tmp);
    convertValue<T>::deriv2ToLaplacian(tmp, result, info());
}
*/


// Returns quantities either on the target domain or on the parametric
// domain depending on the representation of the object
template <typename T>
void gsFunctionSet<T>::compute(const gsMatrix<T> & in,                      
                               gsFuncData<T> & out   ) const                // out is md(an object of gsMapData<T>) in gsVisitorPoisson.h
                                                                            // this function is not only called for mapData, but also for m_ptable and m_itable
{
    gsInfo << "gsFunctionSet<T>::compute()\n";
//     
//     gsInfo << "in: \n";
//     gsInfo << in;
    
    const unsigned flags = out.flags;
    
//     gsInfo << "\n";
//     gsInfo << "flags of out: " << flags << "\n";

    out.dim = this->dimensions();
    //gsDebugVar(&out);

    const int md = out.maxDeriv();                                          // maxDeriv(), a function defined in gsFuncData.h
    
    gsInfo << "md: " << md << "\n";
    
    gsInfo << "\n";
    if (md != -1)
    {   
        gsInfo << "md != 1\n";
        
        gsInfo << "before evalAllDers_into(), out.values: \n";                      // values are independent of the physical length
        for(unsigned int id_item=0; id_item< out.values.size(); id_item++)
        {
            gsInfo << "[" << id_item << "]: \n";

            for (int i = 0; i!=out.values[id_item].rows(); ++i)
            {
                gsInfo << "basis " << i << ": ";
                
                for (int j = 0; j!=out.values[id_item].cols(); ++j)
                {
                    gsInfo << out.values[id_item].coeff(i,j) << " ";
                }
                gsInfo << "\n";
            }
        }     
        
        
        evalAllDers_into(in, md, out.values);                               // function defined at line 87
                                                                            // since md is 1, not only the values but also the gradients are evaluated
                                                                            // this function changes out.values
        
        // used for evaluating rhs at quadrature points
        // when measures equal 1, which we consider now, its value equals coordinates of quadrature points
        // coefs in .xml file affects values here.
        
        gsInfo << "\n";
        gsInfo << "after evalAllDers_into(), out.values: \n";           // this gives the values and gradients of quadrature points of the physical cell on a reference cell
        for(unsigned int id_item=0; id_item< out.values.size(); id_item++)
        {
            gsInfo << "[" << id_item << "]: \n";

            for (int i = 0; i!=out.values[id_item].rows(); ++i)
            {
                gsInfo << "basis " << i << ": ";
                
                for (int j = 0; j!=out.values[id_item].cols(); ++j)
                {
                    gsInfo << out.values[id_item].coeff(i,j) << " ";
                }
                gsInfo << "\n";
            }
        }
    }

    if (flags & NEED_ACTIVE && flags & SAME_ELEMENT)
    {
        gsInfo << "\n";
        gsInfo << "flags & NEED_ACTIVE && flags & SAME_ELEMENT\n";
        
        gsInfo << "coordinates of quadrature points: \n";
        for (int i = 0; i<in.cols(); ++i)
        {
            gsInfo << "col(" << i << "): " << in.col(i) << "\n";
        }
        
        GISMO_ASSERT(0!=in.cols(), "The points are empty.");
        
        active_into(in.col(0), out.actives);                                // this gives the active basis functions at the first quadrature point
        
        gsInfo << "out.actives: \n"
               << out.actives << "\n";   
               
/*        gsInfo << "out.values: \n";
        for(unsigned int i=0; i< out.values.size(); i++)
        {
            gsInfo << "[" << i << "]: \n";
            gsInfo << out.values[i] << "\n";
        }*/               
        
    }
    else if (flags & NEED_ACTIVE)
    {
        gsInfo << "\n";
        gsInfo << "flags & NEED_ACTIVE\n";
//                << "    active_into()\n";
               
        active_into(in, out.actives);
        
//         gsInfo << "md.actives: \n"
//                << out.actives << "\n";
               
/*        gsInfo << "\n";
        gsInfo << "md.values(" << out.values.size() << " in total): \n"
               << "[0]: \n" << out.values[0] << "\n"                                  
               << "[1]: \n" << out.values[1] << "\n";    */             
        
    }
    
//     gsInfo << "flags of out: " << flags << "\n";
    
    // if ( flags & NEED_DIV )
    //     convertValue<T>::derivToDiv(out.values[1], out.divs, info());
    // if ( flags & NEED_CURL )
    //     convertValue<T>::derivToCurl(out.values[1], out.curls, info());
    if (flags & NEED_LAPLACIAN)
    {
        gsInfo << "  flags & NEED_LAPLACIAN\n";
        
        const index_t dsz    = out.deriv2Size();
        const index_t numact = out.values[2].rows() / dsz;
        out.laplacians.resize(numact, in.cols());
        for (index_t i=0; i!= numact; ++i)
            out.laplacians.row(i) =
                out.values[2].middleRows(dsz*i,out.dim.first).colwise().sum();
    }
    
    gsInfo << "\n";
    gsInfo << "finished gsFunctionSet<T>::compute(), resulting in \n";                 
        
    gsInfo << "\n";
    gsInfo << "out.actives: \n"
           << out.actives << "\n";      
    gsInfo << "out.values: \n";
    for(unsigned int i=0; i< out.values.size(); i++)
    {
        gsInfo << "[" << i << "]: \n";
        gsInfo << out.values[i] << "\n";
    } 
}

} // namespace gismo
