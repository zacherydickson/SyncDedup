#include <catch2/catch_all.hpp>
#include <Catch2Extensions.hpp>
#include <DuplicateFilter.h>

using namespace DuplicateFilterNS;


struct OffsetConstructor {
    static const size_t minSize = 10;
    virtual double operator()(std::vector<int> & offVec, size_t denom, double rate) const = 0;
};

//Constructs a Hit candidate where the offset vectors are in the form:
// 0, 1, 2 ..., nIndel, 0, 0, 0, ...
// The number of indels is the closest integer to the rate * the number of indel slots
// in a sequence (seqlength - 1)
//The length may be longer than a possible given denom
struct IndelOffsetConstructor : public OffsetConstructor
{
    static const size_t lengthMult = 2;
    double operator()(std::vector<int> & offVec, size_t denom, double rate) const override
    {
        size_t nIndel = std::round(denom * rate);
        size_t length = lengthMult * nIndel;
        if(length < minSize) { length = minSize; }
        offVec = std::vector<int>(length,int(0));
        for(size_t i = 0; i < nIndel; i++) {
            offVec[i] += i+1;
        }
        return nIndel / double(denom);
    }
};

struct SubsOffsetConstructor : public OffsetConstructor
{
    const Sketcher * sketcher_ptr;
    const SketchPair * sp_ptr;
    double operator()(std::vector<int> & offVec, size_t denom, double rate) const override
    {
        //TODO
        return 0.0;
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
    actualRatePair.first = offConst(hc.first,fqt.segVec[0].seq.length()-1,rate);
    
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



TEST_CASE ("Specific Indel Case","[.SpecIndelCase]") {
    //TODO: Debug why this case fails
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
        size_t seqLen1 = GENERATE(take(10,random(1,1000)));
        size_t seqLen2 = GENERATE(take(10,random(1,1000)));
        bool bDiffLen = GENERATE(true, false);
        if(!bDiffLen) {
            seqLen2 = seqLen1;
        }
        CAPTURE(seqLen1,seqLen2);
        CAPTURE(maxIndelRate);
        IndelFilter indelFilter{maxIndelRate};
        AND_GIVEN("A template") {
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

SCENARIO("Substitution filtering on known cases","[SubtitutionFilter][Unit]") {
    //TODO
}


TEST_CASE("Specific Substitution filtering Case","[SubtitutionFilter][.SpecSubsCase]") {
    //TODO
}

SCENARIO("Substitution filtering issue Discovery","[SubtitutionFilter][.Discovery]") {
    //TODO
}


// ### CandidateDuplicateFinder

//TODO

// ### DuplicateFilter


//TODO





