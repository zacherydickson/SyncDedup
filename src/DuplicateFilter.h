#ifndef DUPLICATE_FILTER_HEADER_GAURD_
#define DUPLICATE_FILTER_HEADER_GAURD_

#include <LocalSyncmerMap.h>
#include <SketchHashFactory.h>
#include <unordered_map>

namespace DuplicateFilterNS {

typedef std::pair<std::vector<int>,std::vector<int>> HitCandidate;
typedef std::unordered_map<size_t,HitCandidate> HitCandidateMap;

struct TemplateSummary_t {
    const double meanQual;
    const size_t length;
    const size_t idx;
    static constexpr double QualTolerance = 0.01;
    bool operator<(const TemplateSummary_t & other) const;
};

struct CandidateDuplicateFinder {
    size_t fqtIdx;
    bool bPaired;
    SketchPair sp;
    HitCandidateMap operator()(const LocalSyncmerMap & sketchMap ) const ;
};

struct IndelFilter {
    double MaxIndelRate;
    //True if the input passes the filter, false otherwise
    bool operator()(const FastqTemplate_t & fqt,
                    const HitCandidate & hc) const ;
};


struct SubstitutionFilter {
    double MaxSubRate;
    //True if the input passes the filter, false otherwise
    bool operator()(size_t k, const FastqTemplate_t & fqt,
                    const SketchPair & sp, const HitCandidate & hc) const ;
};

class DuplicateFilter {
    public:
        DuplicateFilter()
            : DuplicateFilter(DefaultIndelRate,DefaultSubstitutionRate) {};
        DuplicateFilter(DuplicateFilter & other)
            : DuplicateFilter(other.indelFilter_.MaxIndelRate,
                                other.subsFilter_.MaxSubRate) {}
        DuplicateFilter(double indelRate, double subRate )
            : indelFilter_{indelRate}, subsFilter_{subRate} {}

        bool operator()(const Sketcher & sketcher,
                        size_t fqtIdx,
                        const HashedFastqSet & hfqSet) const;

        static constexpr double DefaultIndelRate = 0.0001;
        static constexpr double DefaultSubstitutionRate = 0.01;
    protected:
        const IndelFilter indelFilter_;
        const SubstitutionFilter subsFilter_;

};

} // namespace

#endif //DUPLICATE_FILTER_HEADER_GAURD_



