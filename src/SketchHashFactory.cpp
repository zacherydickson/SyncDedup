#include "SketchHashFactory.h"
#include <future>

SketchHashFactory::SketchHashFactory(Sketcher && sketcher, size_t nThread,
        int phredOffset)
    :   sketcher_(std::move(sketcher)), threadPool_(nThread),
        phredOffset_(phredOffset)
{
}


SketchPair SketchHashFactory::GeneratePairedSketch(const FastqTemplate_t & fqt)
    const 
{
#ifndef NDEBUG
    if(fqt.segVec.size() == 0) {
        throw std::invalid_argument("Attempt to generate paired sketch with an incomplete fastq template");
    }
#endif
    SketchPair sp;
    sp.first = sketcher_.generate_sketch(fqt.segVec[0].seq);
    if(fqt.segVec.size() > 1 ){
        sp.second = sketcher_.generate_sketch(fqt.segVec[1].seq);
    }
    sp.meanQuality = CalculateMeanQuality(fqt);
    return sp;
}

double SketchHashFactory::CalculateMeanQuality(const FastqTemplate_t & fqt) const
{
    double mean = 0.0;
    size_t n = 0;
    for(const auto & seg : fqt.segVec) {
        n += seg.qual.size();
        for(char q : seg.qual) {
           mean += q - this->phredOffset_; 
        }
    }
    mean /= double(n);
    return mean;
}


size_t SketchHashFactory::FillHashedFastqSet( FastqTemplateSource & src,
                                            HashedFastqSet & hfqSet)
{
    if(threadPool_.size() == 0ULL) {
        return NoWorkerFill(src,hfqSet);
    } else {
        return ParallelFill(src,hfqSet);
    }
    return 0ULL;
}


void SketchHashFactory::InsertFqt(  size_t idx, const SketchPair & sp,
                                    const FastqTemplate_t & fqt,
                                    HashedFastqSet & hfqSet)
{
    ExtendedFastqTemplate_t efqt;
    efqt.meanQual = sp.meanQuality;
    hfqSet.sketchMap.insert(2*idx+0,sp.first);
    hfqSet.sketchMap.insert(2*idx+1,sp.second);
    hfqSet.templateVec.push_back(efqt);
}

size_t SketchHashFactory::NoWorkerFill( FastqTemplateSource & src,
                                       HashedFastqSet & hfqSet)
{
    FastqTemplate_t fqt;
    size_t nInserted = 0;
    while( src(fqt) ) {
        size_t fqtIdx = hfqSet.templateVec.size();
        const SketchPair & sp = GeneratePairedSketch(fqt);
        InsertFqt(fqtIdx, sp, fqt, hfqSet);
        nInserted++;
    }
    return nInserted;
}

std::vector<std::future<std::vector<SketchPair>>>
    SketchHashFactory::BatchLaunchSketching(
        FastqTemplateSource & src,size_t batchSize, HashedFastqSet & hfqSet)
{
    if(!batchSize) { batchSize = 1; }
    std::vector<std::future<std::vector<SketchPair>>> futureVec;
    bool bContinue = true;
    while(bContinue) {
        std::vector<FastqTemplate_t> batch;
        if(batchSize > 1){ //For Actual batches
            batch = src.get_block(batchSize);
        } else { //For singlet batches, can use the lower overhead functor
            batch.emplace_back();
            if( ! src(batch.front()) ) { batch.clear(); }
        }
        if(batch.size()) {
            for(auto & fqt : batch){
                ExtendedFastqTemplate_t efqt;
                efqt.name = fqt.name;
                efqt.segVec = fqt.segVec;
                hfqSet.templateVec.push_back(efqt);
            }
            std::future<std::vector<SketchPair>> future = threadPool_.push(
                        [&batch,this](int id) {
                            std::vector<SketchPair> spVec;
                            for(auto & fqt : batch) {
                                spVec.push_back(this->GeneratePairedSketch(fqt));
                            }
                            return spVec;
                        } );
            futureVec.push_back(std::move(future));
        }
        bContinue = (batch.size() == batchSize);
    }
    return futureVec;
}


size_t SketchHashFactory::ParallelFill( FastqTemplateSource & src,
                                        HashedFastqSet & hfqSet)
{
    //Determine batch size 
    size_t batchSize;
    if ( src.get_size(batchSize)) { 
        batchSize /= threadPool_.size();
        if(!batchSize) { batchSize = 1; }
    } else { //src doesn't have defined size
        batchSize = 1;
    }

    //Store where the end of the hfqSet template vec ended
    size_t fqtIdx = hfqSet.templateVec.size();

    std::vector<std::future<std::vector<SketchPair>>> sketchFutureVec =
        this->BatchLaunchSketching(src,batchSize,hfqSet);
    
    size_t nInserted = 0;
    for(auto & future : sketchFutureVec){
        const std::vector<SketchPair> & spVec = future.get();
        for(const auto & sp : spVec) {
            InsertFqt(fqtIdx, sp, hfqSet.templateVec[fqtIdx], hfqSet);
            fqtIdx++;
            nInserted++;
        }
    }
    return nInserted;
}
