#ifndef DUPLICATE_FILTER_HEADER_GAURD_
#define DUPLICATE_FILTER_HEADER_GAURD_

#include <LocalSyncmerMap.h>
#include <SketchHashFactory.h>
#include <unordered_map>

class DuplicateFilter {
    public:
        DuplicateFilter() = delete;
        DuplicateFilter(DuplicateFilter & other)
            : DuplicateFilter(other.substitutionRate_, other.indelRate_) {}
        DuplicateFilter(double subRate, double indelRate) 
            : substitutionRate_(subRate), indelRate_(indelRate) {}

        bool operator()(const Sketcher & sketcher,
                        size_t fqtIdx,
                        const HashedFastqSet & hfqSet) const;
    protected:
        const double substitutionRate_ = 0.01;
        const double indelRate_ = 0.0001;

        using HitCandidate = std::pair<std::vector<size_t>,std::vector<size_t>>;
        using HitCandidateMap = std::unordered_map<size_t,HitCandidate>;
        static HitCandidateMap FindCandidateHits(
                    size_t fqtIdx, bool bPaired, const SketchPair & sp,
                    const HashedFastqSet & hfqSet );
        static SketchPair GeneratePairedSketch( const Sketcher & sketcher, 
                                                const FastqTemplate_t & fqt);

};

#endif //DUPLICATE_FILTER_HEADER_GAURD_



