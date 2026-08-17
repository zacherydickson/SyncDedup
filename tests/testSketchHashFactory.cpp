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


TEST_CASE("Mean Quality is correctly calculated","[SketchHashFactory]") {
    auto phredOffsetOffset = GENERATE(0, 1, 5, -5);
    INFO("and Given: phred offset from default is " << phredOffsetOffset);
    int nThread = 1;
    SketchHashFactory shFactory(
            Sketcher(13,5,1), nThread,
            SketchHashFactory::DefaultPhredOffset + phredOffsetOffset);
    struct testData {
        FastqTemplate_t fqt;
        double expqual;
    };
    std::vector <testData> fqtVec = {
        { initFqt1fwd, initMeanQual1fwd   + phredOffsetOffset  },
        { initFqt1rev, initMeanQual1rev   + phredOffsetOffset  },
        { initFqt1pair, initMeanQual1pair + phredOffsetOffset  },
        { initFqt2fwd, initMeanQual2fwd   + phredOffsetOffset  },
        { initFqt2rev, initMeanQual2rev   + phredOffsetOffset  },
        { initFqt2pair, initMeanQual2pair + phredOffsetOffset  } };
    for( const testData & td : fqtVec) {
        REQUIRE_THAT(   shFactory.CalculateMeanQuality(td.fqt),
                        Catch::Matchers::WithinAbs(td.expqual,0.001)
                     );
    }
}


SCENARIO("Paired Sketches are generated correctly","[SketchHashFactory]") {
    //TODO:
}


SCENARIO("Fqt's are properly inserted into the HashedFastqSet ","[SketchHashFactory]") {
    //TODO:
}


SCENARIO("Filling from FastqIO source","[SketchHashFactory]") {
    //TODO:
}


SCENARIO("Filling from a vector source","[SketchHashFactory]") {
    //TODO:
}
