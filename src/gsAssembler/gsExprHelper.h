/** @file gsExprHelper.h

    @brief Generic expressions helper

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsAssembler/gsExpressions.h>

#include <fstream>

namespace gismo
{
/**
   Class holding an expression environment
 */
template<class T>
class gsExprHelper                                                                  // not much compiling time
{
private:
    gsExprHelper(const gsExprHelper &);

    gsExprHelper() : mesh_ptr(NULL)
    {
        gsInfo << "gsExprHelper<T>::gsExprHelper()\n\n";   
        mutVar.setData(mutData); 
    }

private:
    typedef std::map<const gsFunctionSet<T>*,gsFuncData<T> > FunctionTable;
    typedef typename FunctionTable::iterator ftIterator;
    typedef typename FunctionTable::const_iterator const_ftIterator;

    // variable/space list
    std::deque<expr::gsFeVariable<T> > m_vlist;
    std::deque<expr::gsFeSpace<T> >    m_slist;

    // background functions
    FunctionTable m_itable;
    FunctionTable m_ptable;
    //FunctionTable i_map;

    // geometry map
    expr::gsGeometryMap<T> mapVar;
public:
    gsMapData<T> mapData;                                       // a class defined in gsFuncData.h
private:

    // mutable pair of variable and data,
    // ie. not uniquely assigned to a gsFunctionSet
    expr::gsFeVariable<T> mutVar ;
    gsFuncData<T>         mutData;
    bool mutParametric;

    gsSortedVector<const gsFunctionSet<T>*> evList;

    const gsMultiBasis<T> * mesh_ptr;                           // declaring a pointer variable or assigning values to it does not call any constructor
                                                                // https://stackoverflow.com/questions/44612995/is-there-no-constructor-call-when-pointer-is-created
//     const gsMultiBasis<T> mesh_ptr_custom;

public:
    typedef const expr::gsGeometryMap<T> & geometryMap;
    typedef const expr::gsFeElement<T>   & element;
    typedef const expr::gsFeVariable<T>  & variable;
    typedef const expr::gsFeSpace<T>     & space;
    typedef const expr::gsNullExpr<T>      nullExpr;


    typedef expr::gsFeVariable<T>  & nonConstVariable;
    typedef expr::gsFeSpace<T>     & nonConstSpace;

    typedef memory::unique_ptr<gsExprHelper> uPtr;
    typedef memory::shared_ptr<gsExprHelper>  Ptr;
public:

    gsMatrix<T> & points() { return mapData.points; }

    static uPtr make() 
    { 
        gsInfo << "gsExprHelper<T>::make()\n\n";
        return uPtr(new gsExprHelper());                            // constructor in line 31
    }

    void reset()
    {
        m_ptable.clear();
        m_itable.clear();
        points().clear();
        m_vlist .clear();
        m_slist .clear();
        //mapVar.reset();
    }

    void print() const
    {
        gsInfo << "gsExprHelper<T>::print()\n"
               << "    print info on mapVar, mutVar, m_ptable, m_itable\n\n";
        
        //mapData.side
        if ( mapVar.isValid() ) // list ?
        {
            gsInfo << "mapVar: "<< &mapData <<"\n";
        }

        if ( mutVar.isValid() && 0!=mutData.flags)
        {
            gsInfo << "mutVar: "<< &mutVar <<"\n";
        }

        gsInfo << "ptable:\n";
        for (const_ftIterator it = m_ptable.begin(); it != m_ptable.end(); ++it)
        {
            gsInfo << " * "<< &it->first <<" --> "<< &it->second <<"\n";
        }

        gsInfo << "itable:\n";
        for (const_ftIterator it = m_itable.begin(); it != m_itable.end(); ++it)
        {
            gsInfo << " * "<< &it->first <<" --> "<< &it->second <<"\n";
        }
        
        gsInfo << "\n";
        
    }
    
    
    void print_flags() const
    {
        
        gsInfo << "\nprinting flags:\n";
        
        gsInfo << "mapData: ";
        gsInfo << mapData.flags << "\n";
        
        gsInfo << "mutData: ";
        gsInfo << mutData.flags << "\n";
        
        gsInfo << "ptable: ";
        for (const_ftIterator it = m_ptable.begin(); it != m_ptable.end(); ++it)
        {
            gsInfo << it->second.flags <<" ";
        }

        gsInfo << "\nitable: ";
        for (const_ftIterator it = m_itable.begin(); it != m_itable.end(); ++it)
        {
            gsInfo << it->second.flags <<" ";
        }        
        gsInfo << "\n";
    }

    void cleanUp()
    {
        mapData.clear();
        mutData.clear();
        for (ftIterator it = m_ptable.begin(); it != m_ptable.end(); ++it)
            it->second.clear();
        for (ftIterator it = m_itable.begin(); it != m_itable.end(); ++it)
            it->second.clear();
    }

    void setMultiBasis(const gsMultiBasis<T> & mesh) 
    {
        
        gsInfo << "gsExprHelper<T>::setMultiBasis()\n";
//                << "    mesh_ptr pointing to dbasis\n";
        
        mesh_ptr = &mesh; 
    }

    bool multiBasisSet() { return NULL!=mesh_ptr;}

    const gsMultiBasis<T> & multiBasis()
    {
        GISMO_ASSERT(multiBasisSet(), "Integration elements not set.");
        return *mesh_ptr;
    }


    geometryMap getMap(const gsFunction<T> & mp)
    {
//         gsInfo << "gsExprHelper<T>::getMap() 1\n";
        //mapData.clear();
        mapVar.registerData(mp, mapData);
        return mapVar;
    }

    geometryMap getMap(const gsMultiPatch<T> & mp)
    {
        gsInfo << "gsExprHelper<T>::getMap() 2\n";
               
        //mapData.clear();
        mapVar.registerData(mp, mapData);                                   // mapVar, an object of expr::gsGeometryMap<T>, declared on line 50
                                                                            // mapData, created on line 52, is an object of gsMapData<T> defined in gsFuncData.h
                
/*        gsInfo << "\n";
        print();     */   
        
        return mapVar;
    }

    geometryMap getMap() const
    {
        GISMO_ASSERT(mapVar.isValid(), "The Geometry map is not initialized)");
        return mapVar;
    }

    nonConstVariable getVar(const gsFunctionSet<T> & mp, index_t dim = 1)
    {
        m_vlist.push_back( expr::gsFeVariable<T>() );
        expr::gsFeVariable<T> & var = m_vlist.back();
        gsFuncData<T> & fd = m_ptable[&mp];
        //fd.dim = mp.dimensions();
        //gsDebugVar(&fd);
        var.registerData(mp, fd, dim);
        return var;
    }

    nonConstVariable getVar(const gsFunctionSet<T> & mp, geometryMap G)                        // gsInfo << "G not used\n\n";
    {
//         gsInfo << "gsExprHelper<T>::getVar()\n";
        
        GISMO_UNUSED(G);
        GISMO_ASSERT(&G==&mapVar, "geometry map not known");
        
        m_vlist.push_back( expr::gsFeVariable<T>() );
        
        expr::gsFeVariable<T> & var = m_vlist.back();
        
//         gsInfo << "\ncreating a reference to m_itable[j], where j is f for A.getCoeff(), and ms for ev.getVariable()\n";
        
//         gsInfo << "size of m_itable: " << m_itable.size() << "\n";
        gsFuncData<T> & fd = m_itable[&mp];
        
        
/*        gsInfo << "\n";
        print();     */   
        
        //fd.dim = mp.dimensions();
        //gsDebugVar(&fd);
        var.registerData(mp, fd, 1, mapData);
        return var;
    }

    nonConstSpace getSpace(const gsFunctionSet<T> & mp, index_t dim = 1)
    {
        gsInfo << "gsExprHelper<T>::getSpace()\n";
        
        m_slist.push_back( expr::gsFeSpace<T>() );                  // expr::gsFeSpace<T> inheriting from expr::gsFeVariable<T>
        
//         gsInfo << "\n";
        
        expr::gsFeSpace<T> & var = m_slist.back();
        
//         print();
        
//         gsInfo << "\ncreating a reference to m_ptable[dbasis]\n";
        
        gsFuncData<T> & fd = m_ptable[&mp];                         // first item in m_ptable referring to dbasis
                                                                    // fd referring to the second item in m_ptable
                                                                    // not clear how this assigns values to m_ptable
        
//         gsInfo << "\n";
//         print();
        
//         gsInfo << "\n";
        
//         gsInfo << "fd: " << fd << "\n";
        
        
        //fd.dim = mp.dimensions();
        var.registerData(mp, fd, dim);
        return var;
    }

    //void rmVar(

    bool exists(variable a)
    {
        typedef typename std::deque<expr::gsFeSpace<T> >::const_iterator siter;
        for (siter it = m_slist.begin(); it!=m_slist.end(); ++it)
            if ( &a == &(*it) ) return true;

        typedef typename std::deque<expr::gsFeVariable<T> >::const_iterator viter;
        for (viter it = m_vlist.begin(); it!=m_vlist.end(); ++it)
            if ( &a == &(*it) ) return true;

        return false;
    }

    variable getMutVar() const { return mutVar; }

    void setMutSource(const gsFunction<T> & func, bool param)
    {
        mutVar.setSource(func);
        mutParametric = param;
    }

    template<class E>
    void check(const expr::_expr<E> & testExpr) const
    {
        if ( testExpr.isVector() )
            GISMO_ENSURE(m_ptable.find(&testExpr.rowVar().source())!=m_ptable.end(), "Check failed");
        if ( testExpr.isMatrix() )
            GISMO_ENSURE(m_ptable.find(&testExpr.colVar().source())!=m_ptable.end(), "Check failed");

        // todo: varlist ?
    }

    void initFlags(const unsigned fflag = 0,
                   const unsigned mflag = 0)
    {
        gsInfo << "gsExprHelper<T>::initFlags()\n"; 
//         gsInfo << "    first argu: " << fflag << ", second argu: " << mflag << "\n";
        
        mapData.flags = mflag;
        
        mutData.flags = fflag;
        
        for (ftIterator it = m_ptable.begin(); it != m_ptable.end(); ++it)
            it->second.flags = fflag;
        
        for (ftIterator it = m_itable.begin(); it != m_itable.end(); ++it)
            it->second.flags = fflag;

    }

    template<class Expr> // to remove
    void setFlags(const Expr & testExpr,
                  const unsigned fflag = 0,
                  const unsigned mflag = 0)
    {
        // todo:
        //testExpr.variables_into(m_ptable);
        // plus auto-registration

        initFlags(fflag, mflag);
        testExpr.setFlag(); // protected ?
    }

    //void precompute(const gsMatrix<T> & points, const index_t patchIndex = 0)

    void precompute(const index_t patchIndex = 0, unsigned int id_domain = 0)
    {
        
        gsInfo << "gsExprHelper<T>::precompute()\n";
        
        GISMO_ASSERT(0!=points().size(), "No points");

//         gsInfo << "before dealing with mapData\n";
//         gsInfo << "flags of mapData: " << mapData.flags << "\n";
//         gsInfo << "Info of mapData.values: \n";
        
//         for(unsigned int i=0; i< mapData.values.size(); i++)
//         {
//             gsInfo << "[" << i << "]: \n";
//             gsInfo << mapData.values[i] << "\n";
//         };
        
//         gsInfo << "patchId of mapData: " << mapData.patchId << "\n\n";
        
        gsInfo << "\n";
        //mapData.side
        if ( mapVar.isValid() ) // list ?                               // mapVar, an object of expr::gsGeometryMap<T>
        {
            gsInfo << "#### mapVar\n";
            
            //gsDebugVar("MAPDATA-------***************");
            mapData.flags |= NEED_VALUE;                                // mapData, an object of gsMapData<T>
                                                                        // used in gsFunction<T>::computeMap() and gsFunctionSet<T>::compute()
//             gsInfo << "flags of mapData absorbing NEED_VALUE, which is " << NEED_VALUE << ", resulting in " << mapData.flags << "\n";
            
            
//             const gsMultiPatch<T> & multipatch_inter = mapVar.source();
            
            const gsFunction<T> & bspline_inter=mapVar.source().function(patchIndex);   
                                                                        // source() returning multipatch   
                                                                        // function() returning an object of gsBSpline<T> which inherits from gsGeometry<T>, which inherits from gsFunction<T>
            
            gsInfo << "typeid: " << typeid(bspline_inter).name() << "\n";
            
//             gsInfo << "\n";
//             gsInfo << "basis: \n";
//             gsInfo << &bspline_inter.basis() << "\n";
//             gsInfo << "support: \n";
//             gsInfo << bspline_inter.support() << "\n";
            
//             gsInfo << "nr. of control points: \n";
//             gsInfo << bspline_inter.coefDim() << "\n";            
//             gsInfo << "control points: \n";
//             gsInfo << bspline_inter.coefs() << "\n";            
            
//             gsInfo << "print as a string: \n";
//             bspline_inter.print(gsInfo);
            
            
            
            
            gsInfo << "\n";
            bspline_inter.computeMap(mapData);                          // computeMap(), a function in gsFunction.hpp
                                                                        
                                  
            gsInfo << "\n";
            gsInfo << "after gsFunction<T>::computeMap(), \n";
/*            gsInfo << "mapData.points: \n";            
            gsInfo << mapData.points << "\n"; */                          // same as before  
            gsInfo << "mapData.values: \n";            
            for(unsigned int i=0; i< mapData.values.size(); i++)
            {
                gsInfo << "[" << i << "]: \n";
                gsInfo << mapData.values[i] << "\n";                    // first row equal to the coordinates of the quadrature points
            }                                                           // second row equal to 1
            
            gsInfo << "mapData.measures: \n";                           // equal to 1
            gsInfo << mapData.measures << "\n";            
            
            mapData.patchId = patchIndex;
        }
        
//         gsInfo << "\n";
//         gsInfo << "after gsFunction<T>::computeMap()\n";
//         gsInfo << "flags of mapData: " << mapData.flags << "\n";
//         gsInfo << "mapData.values: \n";
//         for(unsigned int i=0; i< mapData.values.size(); i++)
//         {
//             gsInfo << "[" << i << "]: \n";
//             gsInfo << mapData.values[i] << "\n";
//         }


#if 1
        if ( mutVar.isValid() && 0!=mutData.flags)
        {
            gsInfo << "#### mutVar\n";
//             
//             gsInfo << "mutParametric: " << mutParametric << "\n";
//             gsInfo << "mapData.values.size(): " << mapData.values.size() << "\n";
//             
//             
//             gsInfo << "mapData.values[0]: " << mapData.values[0] << "\n";
            
            GISMO_ASSERT( mutParametric || 0!=mapData.values.size(), "Map values not computed");
            //mutVar.source().piece(patchIndex).compute(mapData.points, mutData);
            mutVar.source().piece(patchIndex)
                .compute( mutParametric ? mapData.points : mapData.values[0], mutData);
        }
        
        unsigned int id_m_ptable = 0;
        unsigned int id_m_itable = 0;

        for (ftIterator it = m_ptable.begin(); it != m_ptable.end(); ++it)                                 // declared on line 46
        {
            gsInfo << "\n";
            gsInfo << "####//////####//////#### m_ptable " << id_m_ptable++ << "\n";
            
//             gsInfo << "flags of it->second: " << it->second.flags << "\n";                                  
                           
            
            gsInfo << "\n";
            gsInfo << "before compute(), actives of the second item of m_ptable: \n";
            gsInfo << it->second.actives;
            
            gsInfo << "\n";
            gsInfo << "values of the second item of m_ptable:\n";
            for(unsigned int i=0; i< it->second.values.size(); i++)
            {
                gsInfo << "[" << i << "]: \n";
                gsInfo << it->second.values[i] << "\n";
            }            
            
            
            gsInfo << "\n";
            //gsDebugVar("-------");
            //gsDebugVar(&it->second);
            //gsDebugVar(it->second.dim.first);
            it->first->piece(patchIndex).compute(mapData.points, it->second); // ! piece(.) ?           // evaluating values and gradients of basis functions at mapData.points
                                                                                                        // horizontally quad-ordered, vertically basis-ordered
                                                                                                        // give values of basis functions at the quadrature points of the physical cell
                                                                                                            
            //gsDebugVar(&it->second);                                                                      
            //gsDebugVar(it->second.dim.first);
            //gsDebugVar("-------");
            it->second.patchId = patchIndex;
            
            std::ofstream fid_value;
            fid_value.open("bspline_fe_values_"+std::to_string(id_domain)+".txt", std::ios::trunc);    
            
            std::ofstream fid_grad;
            fid_grad.open("bspline_fe_gradients_"+std::to_string(id_domain)+".txt", std::ios::trunc);                        
            
            gsInfo << "\n";
            gsInfo << "after compute(), actives of the second item of m_ptable: \n";                    // not changed
            gsInfo << it->second.actives;
            
            gsInfo << "\n";
            gsInfo << "values of the second item of m_ptable:\n";
            for(unsigned int id_item=0; id_item< it->second.values.size(); id_item++)                                     // changed
            {
                gsInfo << "[" << id_item << "]: \n";
                
                for (int i = 0; i!=it->second.values[id_item].rows(); ++i)
                {
                    gsInfo << "basis " << i << ": ";
                    
                    for (int j = 0; j!=it->second.values[id_item].cols(); ++j)
                    {
                        gsInfo << it->second.values[id_item].coeff(i,j) << " ";
                    }
                    gsInfo << "\n";
                    
                    if (id_item==0)
                    {                    
                        for (int j = 0; j!=it->second.values[id_item].cols(); ++j)
                        {
                            fid_value << it->second.values[id_item].coeff(i,j) << " ";
                        }
                        fid_value << "\n";                    
                    }                      
                    else if (id_item==1)
                    {                    
                        for (int j = 0; j!=it->second.values[id_item].cols(); ++j)
                        {
                            fid_grad << it->second.values[id_item].coeff(i,j) << " ";
                        }
                        fid_grad << "\n";                    
                    }                   
                }
            }
            
            fid_value.close();
            fid_grad.close();
            
            
        }
        
        
//         gsInfo << "size of m_itable: " << m_itable.size() << "\n";
        
        GISMO_ASSERT( m_itable.empty() || 0!=mapData.values.size(), "Map values not computed");
        if ( 0!=mapData.values.size() && 0!= mapData.values[0].rows() ) // avoid left-over from previous expr.
        {
            for (ftIterator it = m_itable.begin(); it != m_itable.end(); ++it)
            {
                gsInfo << "\n";
                gsInfo << "#### m_itable " << id_m_itable++ << "\n";  
                
//                 gsInfo << "flags of it->second: " << it->second.flags << "\n";                            
                
                //gsDebugVar(&it->second);
                //gsDebugVar(it->second.dim.first);
                it->first->piece(patchIndex).compute(mapData.values[0], it->second);                        // dealing with the value of the rhs at the input points for the first m_itable
                                                                                                            // do nothing for the second m_itable
                //gsDebugVar(it->second.dim.first);
                it->second.patchId = patchIndex;
            
            }
            
        }
#endif        
        
    }

    template<class E>
    void parse(const expr::_expr<E> & expr)
    {
        //evList.reserve(m_ptable.size()+m_itable.size());
        evList.clear();
        expr.parse(evList);
    }

/*
    void precompute(const index_t patch1, const index_t patch2);

    void precompute(const index_t patch)
    {
        for (ftIterator it = evList.begin(); it != evList.end(); ++it)
        {
            (*it)->piece(patchIndex).compute(mapData.points, it->second); // ! piece(.) ?
            it->second.patchId = patchIndex;
        }
    }
//*/


};//class


} //namespace gismo
