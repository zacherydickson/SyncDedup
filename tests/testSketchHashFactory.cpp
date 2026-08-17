#include <catch2/catch_all.hpp>
#include <Catch2Extensions.hpp>
#include <SketchHashFactory.h>

FastqTemplate_t initFqt1pair = { "Read1",
    {   { "GACTTTTCCAG", "1:N:0:ACTG", "ABCDEFGHIJK" },
        { "CTGGAAAAGTC", "2:N:0:ACTG", "KJIHGFEDCBA" } } };
FastqTemplate_t initFqt1fwd = { initFqt1pair.name,
    { initFqt1pair.segVec[0] } };
FastqTemplate_t initFqt1rev = { initFqt1pair.name,
    { initFqt1pair.segVec[1] } };
double initMeanQual1fwd = 37;
double initMeanQual1rev = 37;
double initMeanQual1pair = 37;

FastqTemplate_t initFqt2pair = { "Read2",
    {   { "CCCCCCCCC", "1:N:0:ACCC", "IIIIIIIII" },
        { "GGGGGGGGGG", "2:N:0:ACCC", "[[[[[[[[[[" } } };
FastqTemplate_t initFqt2fwd = { initFqt2pair.name,
    { initFqt2pair.segVec[0] } };
FastqTemplate_t initFqt2rev = { initFqt2pair.name,
    { initFqt2pair.segVec[1] } };
double initMeanQual2fwd = 40;
double initMeanQual2rev = 58;
double initMeanQual2pair =  (   initMeanQual2fwd *
                                    initFqt2fwd.segVec[0].seq.length() +
                                initMeanQual2rev *
                                    initFqt2rev.segVec[0].seq.length() ) / 
                            (
                                initFqt2fwd.segVec[0].seq.length() +
                                initFqt2rev.segVec[0].seq.length()
                            );

TEST_CASE("SketchHashFactory Construction","[SketchHashFactory][Construction]") {
    auto nWorker = GENERATE(0,1,2,10);
    INFO("and Given: nWorker is " << nWorker);
    auto phredOffset = GENERATE(0,33,64);
    INFO("and Given: phredOffset is " << phredOffset);
    REQUIRE_NOTHROW( SketchHashFactory(Sketcher(7,5,1)) );
}

TEST_CASE("Mean Quality is correctly calculated","[SketchHashFactory]") {
    SketchHashFactory shFactory(Sketcher(7,5,1));
    struct testData {
        FastqTemplate_t fqt;
        double qual;
    };
    std::vector <testData> fqtVec = {   { initFqt1fwd, initMeanQual1fwd },
                                        { initFqt1rev, initMeanQual1rev },
                                        { initFqt1pair, initMeanQual1pair },
                                        { initFqt2fwd, initMeanQual2fwd },
                                        { initFqt2rev, initMeanQual2rev },
                                        { initFqt2pair, initMeanQual2pair } };
    for( const testData & td : fqtVec) {
        REQUIRE_THAT(   shFactory.CalculateMeanQuality(td.fqt),
                        Catch::Matchers::WithinAbs(td.qual,0.001) );
    }
}


SCENARIO("Paired Sketches are generated correctly","[SketchHashFactory]") {
    Sketcher externalSketcher(7,5,1);
    Sketch fwdSk = externalSketcher.generate_sketch(initFqt1fwd.segVec[0].seq);
    Sketch revSk = externalSketcher.generate_sketch(initFqt1rev.segVec[0].seq);
    bool bPaired = GENERATE(true,false);
    INFO("and Given: bPaired is " << bPaired);
    FastqTemplate_t & input = (bPaired) ? initFqt1pair : initFqt1fwd;
    GIVEN("A valid Sketch Factory") {
        SketchHashFactory shFactory(Sketcher(7,5,1));
        WHEN("A Sketch Pair object is generated") {
            SketchPair sp = shFactory.GeneratePairedSketch(input);
            THEN("The sketches match what an external sketcher would generate") {
                REQUIRE( sp.first == fwdSk );
                if(bPaired) {
                    REQUIRE( sp.second == revSk );
                }
            }
            THEN("The meanQuality matches what would be directly calculated") {
                REQUIRE_THAT(   shFactory.CalculateMeanQuality(input),
                                Catch::Matchers::WithinAbs(sp.meanQuality,0.001) );
            }
        }
    }
}


SCENARIO("Fqt's are properly inserted into the HashedFastqSet ","[HashedFastqSet]") {
    bool bPaired = GENERATE(true,false);
    INFO("and Given: bPaired is " << bPaired);
    FastqTemplate_t & input = (bPaired) ? initFqt1pair : initFqt1fwd;
    SketchHashFactory shFactory(Sketcher(7,5,1));
    GIVEN("A fully sketched fastqTemplate and an hfqSet") {
        HashedFastqSet hfqSet;
        size_t startingNTemplate = hfqSet.templateVec.size();
        size_t startingNHit = hfqSet.sketchMap.hits();
        SketchPair sp = shFactory.GeneratePairedSketch(input);
        size_t nHit = sp.first.size();
        if(bPaired) { nHit += sp.second.size(); }
        WHEN("The template and sketch are inserted") {
            hfqSet.insert(sp,input);
            THEN("The number of stored templates and hits increases correctly") {
                REQUIRE( hfqSet.templateVec.size() == startingNTemplate + 1 );
                REQUIRE( hfqSet.sketchMap.hits() == startingNHit + nHit );
            }
        }
        WHEN("The hfqSet has the template inserted prior to sketching") {
            size_t fqtIdx = hfqSet.templateVec.size();
            hfqSet.templateVec.emplace_back(input);
            startingNTemplate = hfqSet.templateVec.size();
            REQUIRE ( startingNHit == hfqSet.sketchMap.hits() );
            AND_WHEN("The sketch is added") {
                hfqSet.add_sketch(fqtIdx,sp);
                THEN("The meanQuality for the template is updated") {
                    REQUIRE ( hfqSet.templateVec[fqtIdx].meanQual == sp.meanQuality );
                }
                THEN("The number of stored hits increases"){
                    REQUIRE ( hfqSet.sketchMap.hits() == startingNHit + nHit );
                }
                THEN("The number of templates is unchanged") {
                    REQUIRE( hfqSet.templateVec.size() == startingNTemplate );
                }
            }
        }
    }
}


SCENARIO("Filling from FastqIO source","[SketchHashFactory]") {
    //TODO:
}


SCENARIO("Filling from a vector source","[SketchHashFactory]") {
    //TODO:
}
