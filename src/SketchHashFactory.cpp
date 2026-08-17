#include "SketchHashFactory.h"
#include <future>


// ################### HashedFastqSet #####################################

void HashedFastqSet::insert(const SketchPair & sp, const FastqTemplate_t & fqt)
{
    size_t idx = templates.size();
    templates.emplace_back(fqt);
    this->add_sketch(idx,sp);
}

void HashedFastqSet::add_sketch(size_t idx, const SketchPair & sp)
{
#ifndef NDEBUG
    if(idx >= templates.size()){
        throw std::logic_error("Attempt to add a sketch for a template that has not been added");
    }
#endif
    templates[idx].meanQual = sp.meanQuality;
    sketchMap.insert(2*idx+0,sp.first);
    sketchMap.insert(2*idx+1,sp.second);
}


// ################### SKETCH HASH FACTORY #####################################

SketchHashFactory::SketchHashFactory(Sketcher && sketcher, size_t nThread,
        int phredOffset)
    :   sketcher_(std::move(sketcher)), threadPool_(nThread),
        phredOffset_(phredOffset)
{
}

std::vector<std::future<std::vector<SketchPair>>>
    SketchHashFactory::BatchLaunchSketching(
        FastqTemplateSource & src,size_t batchSize, HashedFastqSet & hfqSet)
{
    if(!batchSize) { batchSize = 1; }
    std::vector<std::future<std::vector<SketchPair>>> futureVec;
    std::vector<FastqTemplate_t> fqtVec;
    bool bContinue = true;
    while(bContinue) {
        size_t batchStart = hfqSet.templates.size();
        size_t batchEnd = batchStart;
        if(batchSize > 1){ //For Actual batches
            std::vector<FastqTemplate_t> batch;
            batch = src.get_block(batchSize);
            for(auto & fqt : batch){
                hfqSet.templates.emplace_back(fqt);
            }
            batchEnd += batch.size();
        } else { //For singlet batches, can use the functor - lower overhead
            FastqTemplate_t fqt;
            if( src(fqt) ) {
                hfqSet.templates.emplace_back(fqt);
                batchEnd++;
            }
        }
        if(batchEnd > batchStart) {
            std::future<std::vector<SketchPair>> future = threadPool_.push(
                        [&hfqSet,batchStart,batchEnd,this](int id) {
                            std::vector<SketchPair> spVec;
                            for(size_t i = batchStart; i < batchEnd; i++) {
                                FastqTemplate_t & fqt = hfqSet.templates[i];
                                spVec.push_back( this->GeneratePairedSketch(fqt));
                            }
                            return spVec;
                        } );
            futureVec.push_back(std::move(future));
        }
        bContinue = ( (batchEnd - batchStart) >= batchSize);
    } 
    return futureVec;
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





size_t SketchHashFactory::NoWorkerFill( FastqTemplateSource & src,
                                       HashedFastqSet & hfqSet)
{
    FastqTemplate_t fqt;
    size_t nInserted = 0;
    while( src(fqt) ) {
        SketchPair sp = GeneratePairedSketch(fqt);
        hfqSet.insert(sp, fqt);
        nInserted++;
    }
    return nInserted;
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
    size_t fqtIdx = hfqSet.templates.size();

    std::vector<std::future<std::vector<SketchPair>>> sketchFutureVec =
        this->BatchLaunchSketching(src,batchSize,hfqSet);
    
    size_t nInserted = 0;
    for(auto & future : sketchFutureVec){
        const std::vector<SketchPair> & spVec = future.get();
        for(const auto & sp : spVec) {
            hfqSet.add_sketch(fqtIdx++, sp);
            nInserted++;
        }
    }

    return nInserted;
}
