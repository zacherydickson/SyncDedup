#ifndef LOCAL_SYNCMER_MAP_HEADER_GAURD
#define LOCAL_SYNCMER_MAP_HEADER_GAURD

#include <cstddef>
#include <cstdint>
#include "Sketcher.h"
#include <unordered_map>
#include <vector>

struct LocationElement {
    size_t pos;
    size_t id;
};

using Location = std::vector<LocationElement>;

class LocalSyncmerMap {
    struct IdentityFunctor {
        size_t operator()(const uint64_t & a) const {return a; }
    };
    using HashMapType = std::unordered_map<uint64_t,Location,IdentityFunctor>;
    public:
        void clear() { data.clear(); } 
        auto insert(const std::pair<uint64_t,Location> & value) { return data.insert(value); }
        auto insert(std::pair<uint64_t,Location> && value) { return data.insert(value); }
        auto insert(size_t id, const SketchElement & se);
        void insert(size_t id, const Sketch & sketch);
        size_t erase( const uint64_t & key);
        size_t erase( size_t id, const SketchElement & se);
        Location & at(const uint64_t & key ) { return data.at(key); }
        const Location & at(const uint64_t & key ) const { return data.at(key); }
        Location & at(const SketchElement & se ) { return data.at(se.hash); }
        const Location & at(const SketchElement & se ) const { return data.at(se.hash); }
        size_t count( const uint64_t key) const { return data.count(key); }
        size_t count( const SketchElement & se) const { return data.count(se.hash); }
        size_t hits() const { return this->nHit; }
        HashMapType::iterator find( const uint64_t  & key) { return data.find(key); }
        HashMapType::const_iterator find( const uint64_t  & key) const { return data.find(key); }
    protected:
        size_t nHit;
        HashMapType data; 
};


auto LocalSyncmerMap::insert(size_t id, const SketchElement & se) {
    auto it = data.find(se.hash);
    bool bInsert = false;
    if(it == data.end()) {
        bInsert = true;
        data.insert(std::make_pair(se.hash,Location{}));
        it = data.find(se.hash);
    }
    it->second.push_back({se.position,id});
    return std::make_pair(it,bInsert);
}

void LocalSyncmerMap::insert(size_t id, const Sketch & sketch) {
    for(const auto & se : sketch) { this->insert(id,se); }
}


size_t LocalSyncmerMap::erase( const uint64_t & key) {
    auto it = data.find(key);
    if(it != data.end()) {
        nHit -= it->second.size();
        data.erase(it);
        return 1;
    }
    return 0;
}

size_t LocalSyncmerMap::erase( size_t id, const SketchElement & se) {
    auto it = data.find(se.hash);
    size_t initNhit = this->nHit;
    if(it != data.end()){
        for(auto vit = it->second.begin(); vit != it->second.end(); ){
            if(vit->id == id){
                vit = it->second.erase(vit);
                this->nHit--;
            } else {
                vit++;
            }
        }
        //Erase the whole location if there are no hits left
        if(!it->second.size()){
            data.erase(it);
        }
    }
    return this->nHit - initNhit;
}

#endif //LOCAL_SYNCMER_MAP_HEADER_GAURD



