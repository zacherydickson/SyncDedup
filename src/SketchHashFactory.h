#ifndef SKETCH_HASH_FACTORY_HEADER_GAURD_
#define SKETCH_HASH_FACTORY_HEADER_GAURD_

#include <ctpl_stl.h>
#include <Sketcher.h>
#include <LocalSyncmerMap.h>
#include <FastqIO.h>
#include <FastqTemplateSource.h>
#include <vector>

struct ExtendedFastqTemplate_t : public FastqTemplate_t {
    double meanQual = 0.0;
};

struct HashedFastqSet {
    LocalSyncmerMap sketchMap;
    std::vector<ExtendedFastqTemplate_t> templateVec;
};

struct SketchPair {
    Sketch first;
    Sketch second;
    double meanQuality;
};

class SketchHashFactory {
    public:
    SketchHashFactory() = delete;
    SketchHashFactory(  Sketcher && sketcher,size_t nProcThread = 1,
                        int phredOffset = SketchHashFactory::DefaultPhredOffset);
    
    size_t FillHashedFastqSet(FastqTemplateSource & src, HashedFastqSet & hfqSet);
    SketchPair GeneratePairedSketch(const FastqTemplate_t & fqt) const ;
    double CalculateMeanQuality(const FastqTemplate_t & fqt) const ;
    void InsertFqt( size_t idx, const SketchPair& sp, const FastqTemplate_t & fqt,
                    HashedFastqSet & hfqSet);

    static const int DefaultPhredOffset = 33;
    protected:
    //Members
    Sketcher sketcher_;
    ctpl::thread_pool threadPool_;
    const int phredOffset_;

    //Methods
    size_t NoWorkerFill(FastqTemplateSource & src, HashedFastqSet &);
    size_t ParallelFill(FastqTemplateSource & src, HashedFastqSet &);
    std::vector<std::future<std::vector<SketchPair>>> BatchLaunchSketching(
            FastqTemplateSource & src,size_t batchSize, HashedFastqSet & hfqSet);
};

#endif // SKETCH_HASH_FACTORY_HEADER_GAURD_



