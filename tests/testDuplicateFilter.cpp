#include <catch2/catch_all.hpp>
#include <Catch2Extensions.hpp>
#include <DuplicateFilter.h>

using namespace DuplicateFilter;


////Constructs a Hit candidate where the offset vectors are in the form:
//// 0, 1, 2 ..., nIndel, 0, 0, 0, ...
//// The total length is 2x the number of Indels
//// The number of indels is the closest integer to the rate * the number of indel slots
//// in a sequence (seqlength - 1)
////The length may be longer than a possible sketch for the template
////Input - A reference to a hit candidate in which to store the result
////      - a const reference to an fqt template to inform pairing/ sequence lengths
////      - a double target Indel rate for the first segment
////      - a double target Indel rate for the second segment (defaults to 0)
////Output - A pair of doubles for the actual indel rates in the constructed hit candidate
std::pair<double,double> ConstructIndelHC(  HitCandidate & hc,
                                            const FastqTemplate_t fqt,
                                            double targetIndelRate,
                                            double targetIndelRate2=0)
{
    size_t minSize = 10;
    size_t nIndel = std::round((fqt.segVec[0].seq.length()-1) * targetIndelRate);
    size_t length = 2 * nIndel;
    if(length < minSize) { length = minSize; }
    hc.first = std::vector<int>(length,int(0));
    for(size_t i = 0; i < nIndel; i++) {
        hc.first[i] += i+1;
    }
    //for(auto x : hc.first) {
    //    std::cerr << x << "\t";
    //}
    //std::cerr << "\n";
    std::pair<double,double> actualRatePair;
    actualRatePair.first = nIndel / double((fqt.segVec[0].seq.length()-1));
    if(fqt.segVec.size() > 1) {
        HitCandidate tmphc;
        FastqTemplate_t tmpFqt(fqt);
        FastqSegment_t seg= fqt.segVec[1];
        tmpFqt.segVec.clear();
        tmpFqt.segVec.push_back(seg);
        auto rp = ConstructIndelHC(tmphc,tmpFqt,targetIndelRate2);
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
        fqt.segVec.push_back({std::string(length,'C'),"",std::string(length,'F')});
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

bool PairedExpectedPass(double maxRate, std::pair<double,double> ratePair, bool bPaired, size_t denom) {
    bool expected = ExpectedPass(maxRate, ratePair.first, denom);
    if(bPaired) {
        expected = expected && ExpectedPass(maxRate,ratePair.second,denom);
    }
    return expected;
}


// ############## TESTS ####################3

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

SCENARIO ("Indel filtration on known cases", "[IndelFilter][Unit]") {
    GIVEN("An indel filter") {
        double maxIndelRate = GENERATE(0, 0.5, 1.0,
                DuplicateFilter::DuplicateFilter::DefaultIndelRate  );
        CAPTURE(maxIndelRate);
        IndelFilter indelFilter{maxIndelRate};
        AND_GIVEN("A template") {
            size_t seqLen = GENERATE(150);
            CAPTURE(seqLen);
            bool bPaired = GENERATE(true,false);
            CAPTURE(bPaired);
            FastqTemplate_t fqt = ConstructFQT( "mytemplate",seqLen,
                                                bPaired ? seqLen : 0);
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
                auto ratePair = ConstructIndelHC(hc,fqt,targetIndelRate1,
                                                    targetIndelRate2);
                INFO("and Given: indelRates are " << ratePair.first <<
                        (bPaired ? ", " + std::to_string(ratePair.second) : "" ));
                bool expected = PairedExpectedPass( maxIndelRate, ratePair,
                                                    bPaired, seqLen-1);
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
    double maxIndelRate = 0.1970667069006693;
    size_t seqLen = 111;
    bool bPaired = true;
    double targetIndelRate = 0.20381634913327029;
    CAPTURE( maxIndelRate, seqLen, bPaired, targetIndelRate );
    IndelFilter indelFilter{maxIndelRate};
    FastqTemplate_t fqt = ConstructFQT( "mytemplate",seqLen,
                                        bPaired ? seqLen : 0);
    HitCandidate hc;
    auto ratePair = ConstructIndelHC(hc,fqt,targetIndelRate,targetIndelRate);
    INFO("and Given: indelRates are " << ratePair.first <<
            (bPaired ? ", " + std::to_string(ratePair.second) : "" ));
    bool expected = PairedExpectedPass( maxIndelRate, ratePair,
                                        bPaired, seqLen-1);
    bool bRes = indelFilter(fqt,hc);
    REQUIRE( bRes == expected );
}

SCENARIO ("Indel issue Discovery", "[IndelFilter][.Discovery]") {
    GIVEN("An indel filter") {
        double maxIndelRate = GENERATE(take(10,random(0.0,1.0)));
        CAPTURE(maxIndelRate);
        IndelFilter indelFilter{maxIndelRate};
        AND_GIVEN("A template") {
            size_t seqLen = GENERATE(take(10,random(1,1000)));
            CAPTURE(seqLen);
            bool bPaired = GENERATE(true,false);
            CAPTURE(bPaired);
            FastqTemplate_t fqt = ConstructFQT( "mytemplate",seqLen,
                                                bPaired ? seqLen : 0);
            AND_GIVEN("A hit candidate"){
                double targetIndelRate = GENERATE(take(10,random(0.0,1.0)));
                CAPTURE(targetIndelRate);
                HitCandidate hc;
                auto ratePair = ConstructIndelHC(hc,fqt,targetIndelRate,targetIndelRate);
                INFO("and Given: indelRates are " << ratePair.first <<
                        (bPaired ? ", " + std::to_string(ratePair.second) : "" ));
                bool expected = PairedExpectedPass( maxIndelRate, ratePair,
                                                    bPaired, seqLen-1);
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
