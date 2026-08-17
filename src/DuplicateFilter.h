#ifndef DUPLICATE_FILTER_HEADER_GAURD_
#define DUPLICATE_FILTER_HEADER_GAURD_

#include <SketchHashFactory.h>
#include <LocalSyncmerMap.h>

class DuplicateFilter {
    public:
        DuplicateFilter() = delete;
        DuplicateFilter(DuplicateFilter & other)
            : DuplicateFilter(other.substitutionRate_, other.indelRate_) {}
        DuplicateFilter(double subRate, double indelRate) 
            : substitutionRate_(subRate), indelRate_(indelRate) {}

        bool operator()(const ExtendedFastqTemplate_t & fqt,
                        const LocalSyncmerMap & lsMap) const;
    protected:
        const double substitutionRate_ = 0.01;
        const double indelRate_ = 0.0001;

};

#endif //DUPLICATE_FILTER_HEADER_GAURD_



