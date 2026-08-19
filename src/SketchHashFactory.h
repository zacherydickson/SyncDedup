#ifndef SKETCH_HASH_FACTORY_HEADER_GAURD_
#define SKETCH_HASH_FACTORY_HEADER_GAURD_

#include <deque>
#include <ctpl_stl.h>
#include <FastqIO.h>
#include <FastqTemplateSource.h>
#include <LocalSyncmerMap.h>
#include <Sketcher.h>
#include <vector>

struct ExtendedFastqTemplate_t : public FastqTemplate_t {
    ExtendedFastqTemplate_t() : FastqTemplate_t() {}
    ExtendedFastqTemplate_t(const FastqTemplate_t & fqt) : FastqTemplate_t(fqt) {}
    //TODO: Idea Move Mean Quality Calc here?
    double meanQual = 0.0;
};

struct SketchPair {
    Sketch first;
    Sketch second;
    double meanQuality;
    static SketchPair GeneratePairedSketch( const Sketcher & sketcher,
                                            const FastqTemplate_t & fqt);
};

struct HashedFastqSet {
    LocalSyncmerMap sketchMap;
//TODO: Idea, have this be an "IndexibleFastqTemplateSource" to allow for
//easier flipping between high mem (load it all) and low mem (store in
//indexed file) modes
    std::deque<ExtendedFastqTemplate_t> templates;
    void insert(const SketchPair& sp, const FastqTemplate_t & fqt);
    void add_sketch(size_t idx, const SketchPair & sp);
    static bool skparity(size_t sketchId) { return sketchId & 1; } 
    static size_t sk2fqt(size_t sketchId) { return sketchId >> 1; } 
    static size_t fqt2sk(size_t fqtId, bool parity) { 
        return (fqtId << 1) | parity; } 
};

class SketchHashFactory {
    public:
    SketchHashFactory() = delete;
    SketchHashFactory(  Sketcher && sketcher,size_t nProcThread = 1,
                        int phredOffset = SketchHashFactory::DefaultPhredOffset);
    
    size_t FillHashedFastqSet(FastqTemplateSource & src, HashedFastqSet & hfqSet);
    SketchPair GeneratePairedSketchWithQual(const FastqTemplate_t & fqt) const ;
    double CalculateMeanQuality(const FastqTemplate_t & fqt) const ;
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



