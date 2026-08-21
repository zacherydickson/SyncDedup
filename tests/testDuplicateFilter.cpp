#include <catch2/catch_all.hpp>
#include <Catch2Extensions.hpp>
#include <DuplicateFilter.h>
#include <limits>

using namespace DuplicateFilterNS;


struct OffsetConstructor {
    static const size_t lengthMult = 2;
    virtual double operator()(std::vector<int> & offVec, size_t length, double rate) const = 0;
};

//Constructs a Hit candidate where the offset vectors are in the form:
// 0, 1, 2 ..., nIndel, 0, 0, 0, ...
// The number of indels is the closest integer to the rate * the number of indel slots
// in a sequence (seqlength - 1)
//The length may be longer than a possible given denom
struct IndelOffsetConstructor : public OffsetConstructor
{
    static const size_t minSize = 10;
    double operator()(std::vector<int> & offVec, size_t length, double rate) const override
    {
        size_t denom = (length) ? length-1 : 0;
        size_t nIndel = std::round(denom * rate);
        size_t outlength = lengthMult * nIndel;
        if(outlength < minSize) { outlength = minSize; }
        offVec = std::vector<int>(outlength,int(0));
        for(size_t i = 0; i < nIndel; i++) {
            offVec[i] += i+1;
        }
        return nIndel / double(denom);
    }
};

//Constructs Offset vectors with length matching a sketch element
//Substitutions are indicated by maxInt values
//Substitutions are spaced such that the corresponding sketch elements positions
//are at least k_ bp apart
//If the input sketch's elements are all at least k_ bp apart and the sketch is denom length long
// then any proportion of subtitution rates is possible
//Otherwise the rate is capped at
struct SubsOffsetConstructor : public OffsetConstructor
{
    SubsOffsetConstructor() = delete;
    SubsOffsetConstructor(size_t k):
        k_(k) {}
    size_t k_;
    static size_t calcNSub(size_t denom,double rate) {
        return std::round(denom * rate);
    }
    double operator()(std::vector<int> & offVec, size_t length, double rate) const override
    {
        size_t denom = length;
        size_t nSub = calcNSub(denom,rate);
        //size_t length = nSub;
        offVec = std::vector<int>(nSub,std::numeric_limits<int>::max());
        //size_t sketchIdx = 0;
        //size_t actualSubs = 0;
        //for(size_t i = 0; i < nSub && sketchIdx < length; i++){
        //    offVec[sketchIdx] = std::numeric_limits<int>::max();
        //    actualSubs++;
        //    size_t pos = sketch_ptr_->at(sketchIdx).position;
        //    while(++sketchIdx < length && sketch_ptr_->at(sketchIdx).position < pos + k_){}
        //}
        return nSub / double(denom);
    }
    Sketch idealizedSketch(size_t denom, double rate) {
        size_t nSub = calcNSub(denom,rate);
        Sketch sketch;
        for(size_t i = 0; i < nSub; i++){
            SketchElement se; 
            se.hash = i;
            se.position = i*k_;
            sketch.push_back(se);
        }
        return sketch;
    }
    SketchPair idealizedSketchPair(std::pair<size_t,size_t> denomPair, std::pair<double,double> ratePair, bool bPaired) {
        SketchPair sp;
        sp.first = idealizedSketch(denomPair.first,ratePair.first);
        if(bPaired) {
            sp.second = idealizedSketch(denomPair.second,ratePair.second);
        }
        return sp;
    }
};

//Input - A reference to a hit candidate in which to store the result
//      - a const reference to an fqt template to inform pairing/ sequence lengths
//      - an OffsetConstructor object defining how to construct each offset vector
//      - a double target Indel rate for the first segment
//      - a double target Indel rate for the second segment (defaults to 0)
//Output - A pair of doubles for the actual indel rates in the constructed hit candidate
std::pair<double,double> ConstructHC(   HitCandidate & hc,
                                        const FastqTemplate_t & fqt,
                                        const OffsetConstructor & offConst, 
                                        double rate,
                                        double rate2=0)
{
    std::pair<double,double> actualRatePair;
    actualRatePair.first = offConst(hc.first,fqt.segVec[0].seq.length(),rate);
    
    //for(auto x : hc.first) {
    //    std::cerr << x << "\t";
    //}
    //std::cerr << "\n";
    if(fqt.segVec.size() > 1) {
        HitCandidate tmphc;
        FastqTemplate_t tmpFqt(fqt);
        FastqSegment_t seg= fqt.segVec[1];
        tmpFqt.segVec.clear();
        tmpFqt.segVec.push_back(seg);
        auto rp = ConstructHC(tmphc,tmpFqt,offConst,rate2);
        actualRatePair.second = rp.first;
        hc.second = tmphc.first;
    }
    return actualRatePair;
}

FastqTemplate_t ConstructFQT(std::string name, size_t length, size_t length2 = 0) {
    FastqTemplate_t fqt;
    fqt.name = name;
    fqt.segVec.push_back({std::string(length,'A'),"",std::string(length,'I')});
    if(length2) {
        fqt.segVec.push_back({std::string(length2,'C'),"",std::string(length2,'F')});
    }
    return fqt;
}


double AdjustRate(double rate, double logitEpsilon, double minAdj = 0.01) {
    if(rate == 0.0) { return (logitEpsilon < 0.0) ? 0.0 : minAdj; }
    if(rate == 1.0) { return (logitEpsilon > 0.0) ? 1.0 : 1-minAdj; }
    double logit = std::log(rate / (1.0- rate)) + logitEpsilon;
    return 1.0 / (1 + std::exp(-logit));
}

bool ExpectedPass(double maxRate, double rate, size_t denom) {
    size_t maxN = std::ceil(maxRate * denom);
    size_t obsN = std::round(rate * denom);
    return (obsN <= maxN);
}

bool PairedExpectedPass(double maxRate, std::pair<double,double> ratePair, bool bPaired, std::pair<size_t,size_t> denomPair) {
    bool expected = ExpectedPass(maxRate, ratePair.first, denomPair.first);
    if(bPaired) {
        expected = expected && ExpectedPass(maxRate,ratePair.second,denomPair.second);
    }
    return expected;
}


// ############## TESTS ####################3
// ### TemplateSummary_t

SCENARIO ("Template Summaries correctly compare less than",
        "[TemplateSummary_t][Unit]")
{
    TemplateSummary_t base{30.0,200,123};
    double qualAdj = GENERATE(-1.0, 0.0, 1.0);
    INFO("and Given: qualAdj is " << qualAdj);
    int lenAdj = GENERATE(-1, 0, 1);
    INFO("and Given: lenAdj is " << lenAdj); 
    int idxAdj = GENERATE(-1, 0, 1);
    INFO("and Given: idxAdj is " << idxAdj); 
    bool exp = (idxAdj == -1);
    if(lenAdj == 1) {exp = true;} else if(lenAdj == -1) {exp = false;}
    if(qualAdj < -TemplateSummary_t::QualTolerance) {
        exp = false;
    } else if (qualAdj > TemplateSummary_t::QualTolerance) {
        exp = true;
    }
    GIVEN("An adjusted template") {
        TemplateSummary_t adj{  base.meanQual+qualAdj,
                                base.length+lenAdj,
                                base.idx+idxAdj};
        THEN("Then base template is less than the adjusted when appropriate") {
            REQUIRE( (base < adj) == exp );
        }
    }
}


// ### IndelFilter

SCENARIO ("Indel filtration on known cases", "[IndelFilter][Unit]") {
    GIVEN("An indel filter") {
        double maxIndelRate = GENERATE(0, 0.5, 1.0,
                DuplicateFilter::DefaultIndelRate  );
        CAPTURE(maxIndelRate);
        IndelFilter indelFilter{maxIndelRate};
        AND_GIVEN("A template") {
            size_t seqLen1 = GENERATE(150);
            size_t seqLen2 = GENERATE(100);
            bool bDiffLen = GENERATE(true,false);
            if(!bDiffLen) { seqLen2 = seqLen1; }
            CAPTURE(seqLen1,seqLen2);
            bool bPaired = GENERATE(true,false);
            CAPTURE(bPaired);
            FastqTemplate_t fqt = ConstructFQT( "mytemplate",seqLen1,
                                                bPaired ? seqLen2 : 0);
            AND_GIVEN("A hit candidate"){
                double logitEpsilon1 = GENERATE(-1.0,0.0,1.0);
                double logitEpsilon2 = GENERATE(-1.0,0.0,1.0);
                CAPTURE(logitEpsilon1,logitEpsilon2);
                double targetIndelRate1 = AdjustRate(maxIndelRate,logitEpsilon1);
                double targetIndelRate2 = AdjustRate(maxIndelRate,logitEpsilon2);
                if(!bPaired) {targetIndelRate2 = 0; }
                INFO("and Given: targetIndelRates are " <<
                        targetIndelRate1 << ", " << targetIndelRate2);
                HitCandidate hc;
                IndelOffsetConstructor ioc;
                auto ratePair = ConstructHC(hc,fqt,ioc,targetIndelRate1,
                                                    targetIndelRate2);
                INFO("and Given: indelRates are " << ratePair.first <<
                        (bPaired ? ", " + std::to_string(ratePair.second) : "" ));
                bool expected = PairedExpectedPass( maxIndelRate, ratePair,
                                                    bPaired, {seqLen1-1,seqLen2-1});
                WHEN("The filter is called") {
                    bool bRes = indelFilter(fqt,hc);
                    THEN( "The result matches expectation" ) {
                            REQUIRE( bRes == expected );
                    }
                }
            }
        }
    }
}



TEST_CASE ("Specific Indel Case","[IndelFilter][.SpecialCase]") {
    double maxIndelRate = 0.43087373746039748;
    size_t seqLen1 = 182;
    size_t seqLen2 = 155;
    bool bPaired = true;
    double targetIndelRate1 = 0.17560484551515509;
    double targetIndelRate2 = 0.43799889436149331;
    CAPTURE( maxIndelRate, seqLen1, seqLen2, bPaired, targetIndelRate1, targetIndelRate2 );
    IndelFilter indelFilter{maxIndelRate};
    FastqTemplate_t fqt = ConstructFQT( "mytemplate",seqLen1,
                                        bPaired ? seqLen2 : 0);
    HitCandidate hc;
    IndelOffsetConstructor ioc;
    auto ratePair = ConstructHC(hc,fqt,ioc,targetIndelRate1,targetIndelRate2);
    INFO("and Given: indelRates are " << ratePair.first <<
            (bPaired ? ", " + std::to_string(ratePair.second) : "" ));
    bool expected = PairedExpectedPass( maxIndelRate, ratePair,
                                        bPaired, {seqLen1-1,seqLen2-1});
    bool bRes = indelFilter(fqt,hc);
    REQUIRE( bRes == expected );
}

SCENARIO ("Indel issue Discovery", "[IndelFilter][.Discovery]") {
    GIVEN("An indel filter") {
        double maxIndelRate = GENERATE(take(10,random(0.0,1.0)));
        
        CAPTURE(maxIndelRate);
        IndelFilter indelFilter{maxIndelRate};
        AND_GIVEN("A template") {
            size_t seqLen1 = GENERATE(take(10,random(1,1000)));
            size_t seqLen2 = GENERATE(take(10,random(1,1000)));
            bool bDiffLen = GENERATE(true, false);
            if(!bDiffLen) {
                seqLen2 = seqLen1;
            }
            CAPTURE(seqLen1,seqLen2);
            bool bPaired = GENERATE(true,false);
            CAPTURE(bPaired);
            FastqTemplate_t fqt = ConstructFQT( "mytemplate",seqLen1,
                                                bPaired ? seqLen2 : 0);
            AND_GIVEN("A hit candidate"){
                double targetIndelRate1 = GENERATE(take(10,random(0.0,1.0)));
                double targetIndelRate2 = GENERATE(take(10,random(0.0,1.0)));
                bool bDiffRate = GENERATE(true,false);
                if(!bPaired){
                    targetIndelRate2 = 0;
                } else if(!bDiffRate) {
                    targetIndelRate2 = targetIndelRate1;
                }
                CAPTURE(targetIndelRate1,targetIndelRate2);
                IndelOffsetConstructor ioc;
                HitCandidate hc;
                auto ratePair = ConstructHC(hc,fqt,ioc,targetIndelRate1,targetIndelRate2);
                INFO("and Given: indelRates are " << ratePair.first <<
                        (bPaired ? ", " + std::to_string(ratePair.second) : "" ));
                bool expected = PairedExpectedPass( maxIndelRate, ratePair,
                                                    bPaired, {seqLen1-1,seqLen2-1});
                WHEN("The filter is called") {
                    bool bRes = indelFilter(fqt,hc);
                    THEN( "The result matches expectation" ) {
                            REQUIRE( bRes == expected );
                    }
                }
            }
        }
    }
}


// ### SubstitutionFilter

SCENARIO("Substitution filtering on known cases","[SubstitutionFilter][Unit]") {
   GIVEN("An indel filter") {
        double maxSubsRate = GENERATE(0, 0.5, 1.0,
                DuplicateFilter::DefaultSubstitutionRate  );
        CAPTURE(maxSubsRate);
        SubstitutionFilter subsFilter{maxSubsRate};
        AND_GIVEN("A template") {
            size_t seqLen1 = GENERATE(150);
            size_t seqLen2 = GENERATE(100);
            bool bDiffLen = GENERATE(true,false);
            if(!bDiffLen) { seqLen2 = seqLen1; }
            CAPTURE(seqLen1,seqLen2);
            bool bPaired = GENERATE(true,false);
            CAPTURE(bPaired);
            FastqTemplate_t fqt = ConstructFQT( "mytemplate",seqLen1,
                                                bPaired ? seqLen2 : 0);
            AND_GIVEN("A hit candidate and idealized sketchPair"){
                double logitEpsilon1 = GENERATE(-1.0,0.0,1.0);
                double logitEpsilon2 = GENERATE(-1.0,0.0,1.0);
                CAPTURE(logitEpsilon1,logitEpsilon2);
                double targetSubsRate1 = AdjustRate(maxSubsRate,logitEpsilon1);
                double targetSubsRate2 = AdjustRate(maxSubsRate,logitEpsilon2);
                if(!bPaired) {targetSubsRate2 = 0; }
                INFO("and Given: targetSubsRates are " <<
                        targetSubsRate1 << ", " << targetSubsRate2);
                size_t k = GENERATE(17);
                CAPTURE(k);
                SubsOffsetConstructor soc(k);
                SketchPair sp = soc.idealizedSketchPair(
                        {seqLen1,seqLen2}, {targetSubsRate1,targetSubsRate2},
                        bPaired);
                HitCandidate hc;
                auto ratePair = ConstructHC(hc,fqt,soc,targetSubsRate1,
                                                    targetSubsRate2);
                INFO("and Given: indelRates are " << ratePair.first <<
                        (bPaired ? ", " + std::to_string(ratePair.second) : "" ));
                bool expected = PairedExpectedPass( maxSubsRate, ratePair,
                                                    bPaired, {seqLen1,seqLen2});
                WHEN("The filter is called") {
                    bool bRes = subsFilter(k,fqt,sp,hc);
                    THEN( "The result matches expectation" ) {
                            REQUIRE( bRes == expected );
                    }
                }
            }
        }
    }
}


TEST_CASE("Specific Substitution filtering Case","[SubstitutionFilter][.SpecialCase]") {
    size_t k = 16;
    double maxSubsRate = 0.35743931412539642;
    size_t seqLen1 = 711;
    size_t seqLen2 = 711;
    bool bPaired = false;
    double targetSubsRate1 = 0.35817472859355326;
    double targetSubsRate2 = 0.0;
    CAPTURE( maxSubsRate, seqLen1, seqLen2, bPaired, targetSubsRate1, targetSubsRate2 );
    SubstitutionFilter subsFilter{maxSubsRate};
    FastqTemplate_t fqt = ConstructFQT( "mytemplate",seqLen1,
                                        bPaired ? seqLen2 : 0);
    SubsOffsetConstructor soc(k);
    SketchPair sp = soc.idealizedSketchPair(
            {seqLen1,seqLen2}, {targetSubsRate1,targetSubsRate2},
            bPaired);
    HitCandidate hc;
    auto ratePair = ConstructHC(hc,fqt,soc,targetSubsRate1,targetSubsRate2);
    INFO("and Given: subsRates are " << ratePair.first <<
            (bPaired ? ", " + std::to_string(ratePair.second) : "" ));
    bool expected = PairedExpectedPass( maxSubsRate, ratePair,
                                        bPaired, {seqLen1,seqLen2});
    bool bRes = subsFilter(k,fqt,sp,hc);
    REQUIRE( bRes == expected );

}

SCENARIO("Substitution filtering issue Discovery","[SubstitutionFilter][.Discovery]") {
    size_t k = GENERATE(take(5,random(13,31)));
    CAPTURE(k);
    GIVEN("A substitution filter") {
        double maxSubsRate = GENERATE(take(5,random(0.0,1.0)));
        CAPTURE(maxSubsRate);
        SubstitutionFilter subsFilter{maxSubsRate};
        AND_GIVEN("A template") {
            size_t seqLen1 = GENERATE(take(5,random(1,1000)));
            size_t seqLen2 = GENERATE(take(5,random(1,1000)));
            bool bDiffLen = GENERATE(true, false);
            if(!bDiffLen) {
                seqLen2 = seqLen1;
            }
            CAPTURE(seqLen1,seqLen2);
            bool bPaired = GENERATE(true,false);
            CAPTURE(bPaired);
            FastqTemplate_t fqt = ConstructFQT( "mytemplate",seqLen1,
                                                bPaired ? seqLen2 : 0);
            AND_GIVEN("SketchPair and A hit candidate"){
                double targetSubsRate1 = GENERATE(take(5,random(0.0,1.0)));
                double targetSubsRate2 = GENERATE(take(5,random(0.0,1.0)));
                bool bDiffRate = GENERATE(true,false);
                if(!bPaired){
                    targetSubsRate2 = 0;
                } else if(!bDiffRate) {
                    targetSubsRate2 = targetSubsRate1;
                }
                CAPTURE(targetSubsRate1,targetSubsRate2);
                SubsOffsetConstructor soc(k);
                SketchPair sp = soc.idealizedSketchPair(
                        {seqLen1,seqLen2}, {targetSubsRate1,targetSubsRate2},
                        bPaired);
                HitCandidate hc;
                auto ratePair = ConstructHC(hc,fqt,soc,targetSubsRate1,targetSubsRate2);
                INFO("and Given: indelRates are " << ratePair.first <<
                        (bPaired ? ", " + std::to_string(ratePair.second) : "" ));
                bool expected = PairedExpectedPass( maxSubsRate, ratePair,
                                                    bPaired, {seqLen1,seqLen2});
                WHEN("The filter is called") {
                    bool bRes = subsFilter(k,fqt,sp,hc);
                    THEN( "The result matches expectation" ) {
                            REQUIRE( bRes == expected );
                    }
                }
            }
        }
    }
}


// ### CandidateDuplicateFinder

void InsertSketchPair(   LocalSyncmerMap & lsMap, size_t fqtIdx,
                                const Sketch & fwd,
                                const Sketch & rev,
                                bool bPaired, bool bFlipParity = false,
                                bool bIdxMatch = true) 
{
    bool fwdParity = bFlipParity;
    bool revParity = !bFlipParity;
    if(!bIdxMatch) { fqtIdx++; }
    lsMap.insert(HashedFastqSet::fqt2sk(fqtIdx,fwdParity),fwd);
    if(bPaired) {
        lsMap.insert(HashedFastqSet::fqt2sk(fqtIdx,revParity),rev);
    }
}

SCENARIO("CandidateDuplicateFinder finds hits only when it should",
            "[CandidateDuplicateFinder][Unit][EmptyResults]")
{
    size_t soI = 1;
    SketchElement fwdSe = {0x01,1};
    SketchElement revSe = {0x11,1};
    size_t hitI = 5;
    size_t nonSoI = 125;
    SketchElement nonFwdSe = {0x0F,1};
    SketchElement nonRevSe = {0x1F,1};
    GIVEN("A cdf functor"){
        bool bEmptySketch = GENERATE(true,false);
        CAPTURE(bEmptySketch);
        bool bPaired = GENERATE(true,false);
        CAPTURE(bPaired);
        CandidateDuplicateFinder cdf;
        cdf.fqtIdx = 1;
        cdf.bPaired = bPaired;
        if(!bEmptySketch) {
            cdf.sp.first.push_back(fwdSe);
            if(bPaired) { cdf.sp.second.push_back(revSe); }
        }
        AND_GIVEN("A local syncmer map") {
            LocalSyncmerMap lsMap;
            bool bHasSoI = GENERATE(false,true);
            bool bSoIParityFlipped = GENERATE(true,false);
            bool bSoIIdxMatch = GENERATE(true,false);
            if(bHasSoI){
                InsertSketchPair(lsMap, soI, {fwdSe}, {revSe}, bPaired,
                                        bSoIParityFlipped, bSoIIdxMatch);
            }
            bool bHasHit = GENERATE(false,true);
            bool bHitParityFlipped = GENERATE(true,false);
            //The Hit idx matching doesn't matter as long as it is different from the SoI
            bool bHitIdxMatch = GENERATE(true,false); 
            if(bHasHit){
                InsertSketchPair(lsMap, hitI, {fwdSe}, {revSe}, bPaired,
                                        bHitParityFlipped, bHitIdxMatch);
            }
            //These test that other elements in the dictionary don't change the result
            bool bHasNonSoI = GENERATE(false,true);
            bool bNonSoIParityFlipped = GENERATE(true,false);
            bool bNonSoIIdxMatch = GENERATE(true,false);
            if(bHasNonSoI){
                InsertSketchPair(lsMap, nonSoI, {nonFwdSe}, {nonRevSe}, bPaired,
                                        bNonSoIParityFlipped, bNonSoIIdxMatch);
            }
            CAPTURE( bHasSoI, bHasHit, bHasNonSoI,
                     bSoIParityFlipped, bHitParityFlipped, bNonSoIParityFlipped,
                     bSoIIdxMatch, bHitIdxMatch, bNonSoIIdxMatch );
            bool bExpectedNonEmpty = ( !bEmptySketch && 
                bHasSoI && !bSoIParityFlipped && !bHitParityFlipped &&
                bSoIIdxMatch && bHasHit );
            WHEN("The functor is called") {
                HitCandidateMap res = cdf(lsMap);
                THEN("The result is empty") {
                    REQUIRE( (res.size() != 0) == bExpectedNonEmpty );
                }
            }
        }
    }
}

SCENARIO("CandidateDuplicateFinder finds the correct number of hits","[CandidateDuplicateFinder][ResSize]") {
    size_t soI = 0;
    SketchElement fwdSe = {0x01,1};
    SketchElement revSe = {0x11,1};
    SketchElement fwdSeRpt = {0x01,5};
    SketchElement revSeRpt = {0x11,5};
    SketchElement fwdSeSep = {0x02,7};
    SketchElement revSeSep = {0x12,7};
    GIVEN("A cdf functor"){
        bool bPaired = GENERATE(true,false);
        //The number of hits should not depend on the sketch structure
        bool bQueryHasRpt = GENERATE(true,false);
        bool bQueryHasSep = GENERATE(true,false);
        CAPTURE(bPaired, bQueryHasRpt, bQueryHasSep);
        CandidateDuplicateFinder cdf;
        cdf.fqtIdx = soI;
        cdf.bPaired = bPaired;
        cdf.sp.first.push_back(fwdSe);
        if(bPaired) { cdf.sp.second.push_back(revSe); }
        if(bQueryHasRpt){
            cdf.sp.first.push_back(fwdSeRpt);
            if(bPaired) { cdf.sp.second.push_back(revSeRpt); }
        }
        if(bQueryHasSep) {
            cdf.sp.first.push_back(fwdSeSep);
            if(bPaired) { cdf.sp.second.push_back(revSeSep); }
        }
        size_t nHits = GENERATE(0,1,5,10);
        AND_GIVEN("A local syncmer Map with the SoI and " << nHits << " hits" ) {
            //The number of hits should not depend on the db structure
            bool bDbHasRpt = GENERATE(true,false);
            bool bDbHasSep = GENERATE(true,false);
            CAPTURE(bDbHasRpt,bDbHasSep);
            LocalSyncmerMap lsMap;
            InsertSketchPair( lsMap, soI, cdf.sp.first, cdf.sp.second, bPaired );
            for(size_t i = 0; i < nHits; i++){
                Sketch fwd{fwdSe};
                Sketch rev{revSe};
                if(bDbHasRpt) {
                    fwd.push_back(fwdSeRpt);
                    rev.push_back(revSeRpt);
                }
                if(bDbHasSep) {
                    fwd.push_back(fwdSeSep);
                    rev.push_back(revSeSep);
                }
                InsertSketchPair( lsMap, soI + i + 1, fwd, rev, bPaired );
            }
            WHEN("The functor is called") {
                HitCandidateMap res = cdf(lsMap);
                size_t nOrig = res.size();
                THEN("Then "<< nHits <<" hits are found") {
                    REQUIRE ( nOrig == nHits );
                    AND_THEN("The indexes of those candidates are correct") {
                        for(size_t i = 0; i < nHits; i++){
                            size_t idx = soI + i + 1;
                            CAPTURE(idx);
                            REQUIRE ( res.count(idx) > 0 );
                        }
                    }
                }
            }
        }
    }
}

//We assume that we will not find hits we shouldn't and the number of hits is correct
// we must only test that the ids and offsets are correct
SCENARIO("CandidateDuplicateFinder finds the correct number of hits","[CandidateDuplicateFinder][ResCorrect]") {
    SketchElement fwdSe = {0x01,1};
    SketchElement revSe = {0x11,2};
    GIVEN("A one element sketch") {
        bool bPaired = GENERATE(true,false);
        CAPTURE(bPaired);
        CandidateDuplicateFinder cdf;
        cdf.fqtIdx = 1;
        cdf.bPaired = bPaired;
        cdf.sp.first.push_back(fwdSe);
        if(bPaired) { cdf.sp.second.push_back(revSe); }
        int offset = GENERATE(-1,0,1);
        AND_GIVEN("A local syncmer map with a hit offset by " << offset) {
            LocalSyncmerMap lsMap;
            InsertSketchPair(lsMap, cdf.fqtIdx, cdf.sp.first, cdf.sp.second, bPaired);
            size_t hitIdx = 2;
            SketchElement fwd = fwdSe; fwd.position += offset;
            SketchElement rev = revSe; rev.position += offset;
            InsertSketchPair(lsMap, hitIdx, {fwd}, {rev}, bPaired);
            WHEN("The functor is called") {
                HitCandidateMap res = cdf(lsMap);
                REQUIRE(res.size() == 1);
                auto it = res.begin();
                REQUIRE( it->first == hitIdx );
                REQUIRE( it->second.first.size() == 1 );
                if(bPaired) {
                    REQUIRE( it->second.second.size() == 1 );
                }
                THEN("The hit has the correct offset") {
                    REQUIRE( it->second.first.front() == offset);
                    if(bPaired) {
                        REQUIRE( it->second.second.front() == offset );
                    }
                }
            }
        }
    }
    GIVEN("A two element sketch with a repeat") {
        SketchElement fwdSeRpt = {0x01,5};
        SketchElement revSeRpt = {0x11,6};
        bool bPaired = GENERATE(true,false);
        CAPTURE(bPaired);
        CandidateDuplicateFinder cdf;
        cdf.fqtIdx = 1;
        cdf.bPaired = bPaired;
        cdf.sp.first = {fwdSe, fwdSeRpt };
        if(bPaired) { cdf.sp.second = {revSe, revSeRpt}; }
        int offset = GENERATE(-1,0,1);
        AND_GIVEN("A local syncmer map with a hit offset by " << offset) {
            LocalSyncmerMap lsMap;
            InsertSketchPair(lsMap, cdf.fqtIdx, cdf.sp.first, cdf.sp.second, bPaired);
            size_t hitIdx = 2;
            Sketch fwd = cdf.sp.first;
            fwd[0].position += offset; fwd[1].position += offset+5;
            Sketch rev = cdf.sp.second;
            if(bPaired) {
                rev[0].position += offset; rev[1].position += offset+5;
            }
            InsertSketchPair(lsMap, hitIdx, fwd, rev, bPaired);
            WHEN("The functor is called") {
                HitCandidateMap res = cdf(lsMap);
                REQUIRE(res.size() == 1);
                auto it = res.begin();
                REQUIRE( it->first == hitIdx );
                REQUIRE( it->second.first.size() == 2 );
                if(bPaired) {
                    REQUIRE( it->second.second.size() == 2 );
                }
                THEN("The offsets for both elements of the hit are correct"){
                    REQUIRE( it->second.first[0] == offset);
                    REQUIRE( it->second.first[1] == offset + 5);
                    if(bPaired) {
                        REQUIRE( it->second.second[0] == offset );
                        REQUIRE( it->second.second[1] == offset + 5);
                    }
                }
            }
        }
    }
    GIVEN("A two distinct element sketch") {
        SketchElement fwdSeSep = {0x02,7};
        SketchElement revSeSep = {0x12,8};
        bool bPaired = GENERATE(true,false);
        CAPTURE(bPaired);
        CandidateDuplicateFinder cdf;
        cdf.fqtIdx = 1;
        cdf.bPaired = bPaired;
        cdf.sp.first = {fwdSe, fwdSeSep };
        if(bPaired) { cdf.sp.second = {revSe, revSeSep}; }
        AND_GIVEN("A local syncmer map with a hit for only one of the two elements ") {
            int offset = GENERATE(-1,0,1);
            size_t seIdx = GENERATE(0,1);
            CAPTURE (offset, seIdx);
            LocalSyncmerMap lsMap;
            InsertSketchPair(lsMap, cdf.fqtIdx, cdf.sp.first, cdf.sp.second, bPaired);
            size_t hitIdx = 2;
            SketchElement fwd = cdf.sp.first[seIdx];
            SketchElement rev;
            fwd.position += offset;
            if(bPaired) {
                rev = cdf.sp.second[seIdx];
                rev.position += offset;
            }
            InsertSketchPair(lsMap, hitIdx, {fwd}, {rev}, bPaired);
            WHEN("The functor is called") {
                HitCandidateMap res = cdf(lsMap);
                REQUIRE(res.size() == 1);
                auto it = res.begin();
                REQUIRE( it->first == hitIdx );
                THEN("There are still 2 offsets in the result") {
                    REQUIRE( it->second.first.size() == 2 );
                    if(bPaired) {
                        REQUIRE( it->second.second.size() == 2 );
                    }
                    AND_THEN("The offsets for the one elements of the hit are correct"){
                        for(size_t i = 0; i < 2; i++){
                            int expect = (i == seIdx) ? offset : std::numeric_limits<int>::max();
                            REQUIRE( it->second.first[i] == expect );
                        }
                        if(bPaired) {
                            for(size_t i = 0; i < 2; i++){
                                int expect = (i == seIdx) ? offset : std::numeric_limits<int>::max();
                                REQUIRE( it->second.second[i] == expect );
                            }
                        }
                    }
                }
            }
        }
    }
}

    
// ### DuplicateFilter


//TODO





