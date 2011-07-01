/****************************************************************************************[IntMap.h]
Copyright (c) 2011, Niklas Sorensson
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef Minisat_IntMap_h
#define Minisat_IntMap_h

#include "minisat/mtl/Vec.h"

namespace Minisat {

    template<class T> struct MkIndexDefault {
        typename vec<T>::Size operator()(T t) const { return (typename vec<T>::Size)t; }
    };
    
    template<class K, class V, class MkIndex = MkIndexDefault<K> >
    class IntMap {
        vec<V>   map;
        MkIndex  index;
    public:
        explicit IntMap(MkIndex _index = MkIndex()) : index(_index){}
        
        const V& operator[](K k) const { assert(index(k) < map.size()); return map[index(k)]; }
        V&       operator[](K k)       { assert(index(k) < map.size()); return map[index(k)]; }

        void     reserve(K key, V pad)       { map.growTo(index(key)+1, pad); }
        void     reserve(K key)              { map.growTo(index(key)+1); }
        void     insert (K key, V val, V pad){ reserve(key, pad); operator[](key) = val; }
        void     insert (K key, V val)       { reserve(key); operator[](key) = val; }

        void     clear  (bool dispose = false) { map.clear(false); }
        void     moveTo (IntMap& to)           { map.moveTo(to.map); to.index = index; }
        void     copyTo (IntMap& to) const     { map.copyTo(to.map); to.index = index; }
    };


    #if 0
    template<class K, class V, V nil, class MkIndex = MkIndexDefault<K> >
    class IntMapNil {
        vec<V> map;
        V      nil;

    public:
        IntMap(){}
        
        void     reserve(K);
        V&       find   (K);
        const V& operator[](K k) const;

    };
    #endif

//=================================================================================================
} // namespace Minisat
#endif