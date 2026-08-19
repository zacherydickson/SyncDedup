#include <DuplicateFilter.h>
#include <cmath>
#include <limits>

namespace DuplicateFilter {

// TEMPLATE SUMMARY

bool TemplateSummary_t::operator<(const TemplateSummary_t & other) const {
    if(std::abs(meanQual - other.meanQual) > TemplateSummary_t::QualTolerance) {
        return meanQual < other.meanQual; 
    }
    if(length != other.length) {
        return length < other.length; 
    }
    return idx > other.idx;
}

// CANDIDATE DUPLICATE FINDER

HitCandidateMap CandidateDuplicateFinder::operator()(
        const HashedFastqSet & hfqSet ) const
{
    HitCandidateMap candMap;
    HitCandidate defaultHitCand = { std::vector<int>(
                                        sp.first.size(),
                                        std::numeric_limits<int>::max() ),
                                    std::vector<int>(
                                        (bPaired ? sp.first.size() : 0),
                                        std::numeric_limits<int>::max() ) };
    for(uint8_t parity = 0; parity <= (bPaired ? 1 : 0); parity++){
        const Sketch & sketch = (parity == 0) ? sp.first : sp.second;
        //Sketches appear in their order within a sequence
        //we can track the last place we considered in order to determine
        //which position the current se of interest occured at in its sequence
        //if the syncmer occured multiple times
        //i.e. the smallest pos greater than the last pos
        size_t lastSEoIPos = 0;
        for(size_t sketchElementIndex = 0; sketchElementIndex < sketch.size(); sketchElementIndex++) {
            const SketchElement & se = sketch[sketchElementIndex];
            size_t curSEoIPos = -1ULL;
            std::vector<size_t> hitIdxVec;
            //Track all positions within a particular template the syncmer occured
            std::unordered_map<size_t,std::vector<size_t>> posMap;
            //Find the non-self hits and determine the current SE of interest position
            for(const LocationElement & le : hfqSet.sketchMap.at(se)) {
                //Skip hits with different parity
                if( HashedFastqSet::skparity(le.id) != parity ) { continue; }
                //Determine the location of sketch in the current template of interest
                size_t curFqtIdx = HashedFastqSet::sk2fqt(le.id);
                if( curFqtIdx == fqtIdx ) {
                    if(le.pos > lastSEoIPos && le.pos < curSEoIPos){
                        curSEoIPos = le.pos;
                    }
                }
                posMap[curFqtIdx].push_back(le.pos);
            }
#ifndef NDEBUG
            if(curSEoIPos == -1ULL){
                throw std::logic_error("Attempt to search for a sketchElement which is not in the HashedFastqSet");
            }
#endif
            //Determine the hit position which is absolutely closest to the curSEoIPos
            for(const auto & pair : posMap) {
                if(pair.first == fqtIdx) { continue; }
                int minOff = std::numeric_limits<int>::max();
                int globalMinOff = minOff;
                for(size_t i = 0; i < pair.second.size(); i++){
                    int curOff = pair.second[i] - curSEoIPos;
                    if(std::abs(curOff) < std::abs(minOff)){
                        minOff = curOff;
                    }
                    //Check against all sketch of interest positions
                    for(size_t seoIPos : posMap[fqtIdx]){
                        int off = pair.second[i] - seoIPos;
                        if(std::abs(off) < std::abs(globalMinOff)){
                            globalMinOff = off;
                        }
                    }
                }
                //If a syncmer is in multiple positions (query or subject) then
                // the syncmer pos should only be assigned to the closest match
                if(minOff == globalMinOff) { 
                    auto resPair = candMap.insert({pair.first,defaultHitCand});
                    auto it = resPair.first;
                    std::vector<int> & offSetVec =  (parity == 0) ?
                                                    it->second.first :
                                                    it->second.second;
                    offSetVec[sketchElementIndex] = minOff;
                }
            }
        }
    }
    return candMap;
}


// INDEL FILTER

bool IndelFilter::operator()(   const FastqTemplate_t & fqt,
                                const HitCandidate & hc) const
{
    bool bPaired = fqt.segVec.size() > 1;
    for(uint8_t parity = 0; parity < (bPaired ? 2 : 1); parity++) {
        const std::vector<int> & offVec = (parity == 0) ? hc.first : hc.second;
        size_t maxIndels = std::ceil(double(fqt.segVec[parity].seq.length()-1) * MaxIndelRate);
        size_t count = 0;
        size_t i0 = 0;
        while(i0 < offVec.size() && offVec[i0] == std::numeric_limits<int>::max()) { i0++; }
        if(i0 >= offVec.size()) { continue; }
        size_t lastOffset = offVec[i0];
        size_t indelCount = 0;
        size_t zeroOffCount = (lastOffset == 0) ? 1 : 0;
        for(size_t i = i0+1; i < offVec.size(); i++){
            size_t off = offVec[i];
            if(off == std::numeric_limits<int>::max()) { continue; } //Skip unmatched sketches
            count++;
            if(off != lastOffset) { indelCount++; }
            if(off == 0) { zeroOffCount++; }
            lastOffset = off;
        }
        //std::cerr << count << "\t" << zeroOffCount << "\t" << indelCount << " vs " << maxIndels << "\n";
        if( !count ) { continue; }
        if( !zeroOffCount ) { return false; } //There are no zero-offset sketches
        if( indelCount > maxIndels ) { return false; } //There are too many offset swaps
    }
    return true;
}

// SUBSTITUTION FILTER

bool SubstitutionFilter::operator()(const Sketcher & sketcher,
                                    const FastqTemplate_t & fqt,
                                    const SketchPair & sp,
                                    const HitCandidate & hc) const
{
    using MMInterval = std::pair<size_t,size_t>;
    size_t k = sketcher.k();
    bool bPaired = fqt.segVec.size() > 1;
    std::vector<MMInterval> mmIntervals;
    for(uint8_t parity = 0; parity < (bPaired ? 2 : 1); parity++) {
        const Sketch & sketch = (parity == 0) ? sp.first : sp.second;
        const std::vector<int> & offVec = (parity == 1) ? hc.first : hc.second;
        size_t maxSubs = std::ceil(double(fqt.segVec[parity].seq.length()) * MaxSubRate);
        //Collect the unique intervals which must contain at least one substitution
        for(size_t i = 0; i < sketch.size(); i++){
            size_t pos = sketch[i].position;
            if(offVec[i] == std::numeric_limits<int>::max()) {
                //This sketch element has no match in the candidate
                //There must be a mismatch within this kmer
                MMInterval mmI = {pos,pos+k-1};
                //Add this interval, or expand the last interval as appropriate
                if( mmIntervals.empty() ||  //No intervals to compare
                    mmI.first > mmIntervals.back().second || // starts after the last ends
                    mmI.first >= mmIntervals.back().first + k) // starts after the first kmer in the last ends
                {
                    mmIntervals.push_back(mmI);
                } else {
                    mmIntervals.back().second = mmI.first;
                }
            }
        }
        size_t nMinSub = mmIntervals.size();
        if(nMinSub > maxSubs) { return false; }
    }
    return true;
}

// DUPLICATEFILTER

bool DuplicateFilter::operator()(   const Sketcher & sketcher, 
                                    size_t fqtIdx,
                                    const HashedFastqSet & hfqSet) const
{
#ifndef NDEBUG
    if(fqtIdx > hfqSet.templates.size()) {
        throw std::logic_error("Attempt to use a Duplicate filter on an unknown template index " + std::to_string(fqtIdx));
    }
#endif
    const ExtendedFastqTemplate_t & fqt = hfqSet.templates[fqtIdx];
    TemplateSummary_t fqtSummary{ fqt.meanQual, fqt.length(), fqtIdx };
    //Find candidate matches for the template
    CandidateDuplicateFinder cdf{   fqtIdx, bool(fqt.segVec.size() > 1),
                                    SketchPair::GeneratePairedSketch(sketcher,fqt) };
    //Cannot determine dup status for an unsketchable sequence
    if(cdf.sp.first.size() + cdf.sp.second.size() <= 0) {
        return true;
    }
    HitCandidateMap candMap = cdf(hfqSet);
    //Filter candidates with excessive indels or substitutions
    //  Indels must be handled first, as the assumption for the sub filter
    //  is that all positions correspond
    for(const auto & pair : candMap) {
        if(indelFilter_(fqt,pair.second) && subsFilter_(sketcher,fqt,cdf.sp,pair.second)) {
            const ExtendedFastqTemplate_t & hitFqt = hfqSet.templates[pair.first];
            TemplateSummary_t hitSummary{hitFqt.meanQual,hitFqt.length(),pair.first};
            if(fqtSummary < hitSummary) { return false; }
        }
    }
    return true;
}


} // namespace
