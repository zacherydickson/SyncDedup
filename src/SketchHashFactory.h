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

struct FastqTemplateSource {
    virtual FastqTemplate_t operator()() = 0;
};

struct FastqIOAsSource : public FastqTemplateSource {
    FastqIOAsSource(FastqIO && handler) : in(std::move(handler)) {
        if(!handler.isReader()) { throw std::invalid_argument("Attempt to use a non reader FastqIO object as a FastQTemplateSource"); }
    }
    FastqTemplate_t operator()() override {
        //TODO: Clean up error reporting and move implementation
        FastqTemplate_t fqt;
        FastqIO::READ_RESULT res = in.next_template(fqt);
        if(res != FastqIO::READ_PASS && res != FastqIO::READ_EOF ) {
            throw std::invalid_argument("Malformed fastq entry");
        }
        return fqt;
    }
    FastqIO in;
};

class SketchHashFactory {
    public:
    SketchHashFactory() = delete;
    SketchHashFactory(  Sketcher && sketcher,size_t nProcThread = 1,
                        int phredOffset = 33);
    
    HashedFastqSet BuildHashedFastqSet(FastqIO & in);
    //TODO: Switch implementation over to generic fqtSources with a HashedFastqSet reference
    //  This lets one expand an existing hash from multiple sources
    //  Ultimately this allows Reading part of a temporary stream and starting the hash,
    //  and then finishing with the rest of the stream
    void BuildHashedFastqSet(FastqTemplateSource & src, HashedFastqSet &);
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



