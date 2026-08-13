#include "SketchHashFactory.h"
#include <future>

SketchHashFactory::SketchHashFactory(Sketcher && sketcher, size_t nThread,
        int phredOffset)
    :   sketcher_(std::move(sketcher)), threadPool_(nThread),
        phredOffset_(phredOffset)
{
}


SketchPair SketchHashFactory::GeneratePairedSketch(const FastqTemplate_t & fqt)
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
    return sp;
}

double SketchHashFactory::CalculateMeanQuality(const FastqTemplate_t & fqt) {
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


HashedFastqSet SketchHashFactory::BuildHashedFastqSet(FastqIO & in) {
    if(threadPool_.size() == 0) {
        return NoWorkerBuild(in);
    } else {
        return ParallelBuild(in);
    }
}


HashedFastqSet SketchHashFactory::NoWorkerBuild(FastqIO & in) {
    HashedFastqSet hfqSet;
    FastqIO::READ_RESULT res = (in.canRead() ? FastqIO::READ_PASS : FastqIO::READ_EOF);
    while( res == FastqIO::READ_PASS) {
        size_t fqtIdx = hfqSet.templateVec.size();
        hfqSet.templateVec.emplace_back();
        ExtendedFastqTemplate_t & fqt = hfqSet.templateVec.back();
        res = in.next_template(fqt);
        if(res == FastqIO::READ_PASS) {
            const SketchPair & sp = GeneratePairedSketch(fqt);
            fqt.meanQual = CalculateMeanQuality(fqt);
            hfqSet.sketchMap.insert(2*fqtIdx+0,sp.first);
            hfqSet.sketchMap.insert(2*fqtIdx+1,sp.second);
        } else if(res != FastqIO::READ_EOF){
            throw std::invalid_argument("Improperly formatted fastq input");
        }
    }
    return hfqSet;
}

HashedFastqSet SketchHashFactory::ParallelBuild(FastqIO & in) {
    HashedFastqSet hfqSet;
    std::vector<std::future<SketchPair>> sketchFutureVec;
    FastqIO::READ_RESULT res = (in.canRead() ? FastqIO::READ_PASS : FastqIO::READ_EOF);
    while( res == FastqIO::READ_PASS) {
        hfqSet.templateVec.emplace_back();
        ExtendedFastqTemplate_t & fqt = hfqSet.templateVec.back();
        res = in.next_template(fqt);
        if(res == FastqIO::READ_PASS) {
            std::future<SketchPair> future = threadPool_.push(
                    [&fqt,this](int id) {
                        SketchPair sp;
                        return this->GeneratePairedSketch(fqt);
                    } );
            sketchFutureVec.push_back(std::move(future));
        } else if(res != FastqIO::READ_EOF){
            throw std::invalid_argument("Improperly formatted fastq input");
        }
    }
    size_t fqtIdx = 0;
    for(auto & future : sketchFutureVec){
        const SketchPair & sp = future.get();
        ExtendedFastqTemplate_t & fqt = hfqSet.templateVec[fqtIdx];
        fqt.meanQual = CalculateMeanQuality(fqt);
        hfqSet.sketchMap.insert(2*fqtIdx+0,sp.first);
        hfqSet.sketchMap.insert(2*fqtIdx+1,sp.second);
    }
    return hfqSet;
}
