#ifndef SKETCH_HASH_FACTORY_HEADER_GAURD_
#define SKETCH_HASH_FACTORY_HEADER_GAURD_

#include <ctpl_stl.h>
#include <Sketcher.h>
#include <LocalSyncmerMap.h>
#include <FastqIO.h>
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
};

class SketchHashFactory {
    public:
    SketchHashFactory() = delete;
    SketchHashFactory(  Sketcher && sketcher,size_t nProcThread = 1,
                        int phredOffset = 33);
    
    HashedFastqSet BuildHashedFastqSet(FastqIO & in);
    protected:
    //Members
    Sketcher sketcher_;
    ctpl::thread_pool threadPool_;
    const int phredOffset_;

    //Methods
    HashedFastqSet NoWorkerBuild(FastqIO & in);
    HashedFastqSet ParallelBuild(FastqIO & in);
    SketchPair GeneratePairedSketch(const FastqTemplate_t & fqt);
    double CalculateMeanQuality(const FastqTemplate_t & fqt);
};

#endif // SKETCH_HASH_FACTORY_HEADER_GAURD_



