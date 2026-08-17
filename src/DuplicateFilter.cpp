#include <DuplicateFilter.h>

bool DuplicateFilter::operator()(   const Sketcher & sketcher, 
                                    size_t fqtIdx,
                                    const HashedFastqSet & hfqSet) const
{
#ifndef NDEBUG
    if(fqtIdx > hfqSet.templates.size()) {
        throw std::logic_error("Attempt to use a Duplicate filter on an unknown template index " + std::to_string(fqtIdx));
    }
#endif
    const FastqTemplate_t & fqt = hfqSet.templates[fqtIdx];
    bool bPaired = fqt.segVec.size() > 1;
    SketchPair sp = GeneratePairedSketch(sketcher,fqt);
    //Find candidate matches for the template
    HitCandidateMap candMap = FindCandidateHits(fqtIdx,bPaired,sp,hfqSet); 
    //TODO: Construct synteny blocks from candidates
    //TODO: Filter on excessive indels
    //TODO: Filter on excessive subs
    return true;
}


DuplicateFilter::HitCandidateMap DuplicateFilter::FindCandidateHits(
        size_t fqtIdx, bool bPaired, const SketchPair & sp,
        const HashedFastqSet & hfqSet )
{
    //TODO: Finish Implementation
    HitCandidateMap candMap;
    for(uint8_t parity = 0; parity <= (bPaired ? 1 : 0); parity++){
        const Sketch & sketch = (parity == 0) ? sp.first : sp.second;
        for(const SketchElement & se : sketch) {
            std::vector<size_t> locPosVec;
            for(const LocationElement & le : hfqSet.sketchMap.at(se)) {
                //Skip hits with different parity
                if( HashedFastqSet::skparity(le.id) != parity ) { continue; }
                //Skip the template itself, after noting its location
                size_t curFqtIdx = HashedFastqSet::sk2fqt(le.id);
                if( curFqtIdx == fqtIdx ) {
                    locPosVec.push_back(le.pos);
                    continue;
                }
                std::vector<size_t> & candVec = candMap[curFqtIdx]
            }
        }
    }
    return candMap;
}

SketchPair DuplicateFilter::GeneratePairedSketch(   const Sketcher & sketcher, 
                                                    const FastqTemplate_t & fqt) 
{
    SketchPair sp;
    sp.first = sketcher.generate_sketch(fqt.segVec[0].seq); 
    if(fqt.segVec.size() > 1){
        sp.second = sketcher.generate_sketch(fqt.segVec[1].seq); 
    }
    return sp;
}

