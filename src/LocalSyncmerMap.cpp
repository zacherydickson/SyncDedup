#include "LocalSyncmerMap.h"

void LocalSyncmerMap::clear() {
    nHit = 0;
    data.clear();
}


std::pair<LocalSyncmerMap::HashMapType::const_iterator,bool>
    LocalSyncmerMap::insert(const std::pair<uint64_t,Location> & value)
{
    nHit += value.second.size();
    return data.insert(value);
}

std::pair<LocalSyncmerMap::HashMapType::const_iterator,bool>
    LocalSyncmerMap::insert(std::pair<uint64_t,Location> && value)
{
    nHit += value.second.size();
    return data.insert(value);
}

std::pair<LocalSyncmerMap::HashMapType::const_iterator,bool>
    LocalSyncmerMap::insert(size_t id, const SketchElement & se)
{
    auto it = data.find(se.hash);
    bool bInsert = false;
    if(it == data.end()) {
        bInsert = true;
        data.insert(std::make_pair(se.hash,Location{}));
        it = data.find(se.hash);
    }
    it->second.push_back({se.position,id});
    nHit++;
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
