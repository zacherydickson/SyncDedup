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

TEST_CASE("CandidateDuplicateFinder correctly finds candidates","[CandidateDuplicateFinder][Unit]") {
    SketchPair sp1;
    sp1.first = {{0x01,1},{0x02,2},{0x03,3},{0x01,4},{0x04,5}};
    sp1.second = {{0x11,1},{0x12,2},{0x13,3},{0x11,4},{0x14,5}};
    SketchPair sp2;
    sp2.first = {{0x01,1},{0x02,2},{0x03,3},{0x01,4},{0x24,5}};
    sp2.second = {{0x21,1},{0x12,2},{0x13,3},{0x11,4},{0x14,5}};
    GIVEN("A Local SyncmerMap") {
        //TODO:
    }
}

//TODO

// ### DuplicateFilter


//TODO





