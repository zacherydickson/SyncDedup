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

std::unique_ptr<std::iostream > constructSS(
        std::vector<FastqTemplate_t> initVec = {}, bool postappend=false,
        FastqIO::READ_RESULT res = FastqIO::READ_PASS)
{
    std::unique_ptr<std::iostream> ss_ptr(new std::stringstream());
    for(const FastqTemplate_t & fqt : initVec ){
        for(const FastqSegment_t & seg : fqt.segVec) {
            if(res != FastqIO::READ_MISSING_LEADER1) {
                *ss_ptr <<  '@';
            }
            *ss_ptr << fqt.name << ' ' << seg.desc << "\n";
            if(res == FastqIO::READ_EOF) { goto endloop; }
            *ss_ptr <<seg.seq << "\n";
            if(res != FastqIO::READ_MISSING_LEADER2) {
                *ss_ptr << '+';
            }
            *ss_ptr << "\n" << seg.qual << "\n";
        }
    }
    //std::cerr << dynamic_cast<std::stringstream*>(ss_ptr.get())->str() << "\n";
endloop:
    if(!postappend){
        ss_ptr->seekg(0);
        ss_ptr->seekp(0);
    }
    return ss_ptr;
}

TEST_CASE("SketchHashFactory Construction","[SketchHashFactory][Construction]") {
    auto nWorker = GENERATE(0,1,2,10);
    INFO("and Given: nWorker is " << nWorker);
    auto phredOffset = GENERATE(0,33,64);
    INFO("and Given: phredOffset is " << phredOffset);
    REQUIRE_NOTHROW( SketchHashFactory(Sketcher(7,5,1),nWorker,phredOffset) );
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
    bool bPreInsert = GENERATE(true,false);
    INFO("and Given: bPreInsert is " << bPreInsert);
    FastqTemplate_t & input = (bPaired) ? initFqt1pair : initFqt1fwd;
    SketchHashFactory shFactory(Sketcher(7,5,1));
    GIVEN("A fully sketched fastqTemplate and an hfqSet") {
        HashedFastqSet hfqSet;
        size_t startingNTemplate = hfqSet.templateVec.size();
        size_t startingNHit = hfqSet.sketchMap.hits();
        size_t fqtIdx = 0;
        size_t nTemplateIncrease = (bPreInsert) ? 0 : 1;
        if(bPreInsert) {
            fqtIdx = hfqSet.templateVec.size();
            hfqSet.templateVec.emplace_back(input);
            startingNTemplate = hfqSet.templateVec.size();
            REQUIRE ( startingNHit == hfqSet.sketchMap.hits() );
        }
        SketchPair sp = shFactory.GeneratePairedSketch(input);
        size_t nHit = sp.first.size();
        if(bPaired) { nHit += sp.second.size(); }
        WHEN("The template and sketch are added/inserted") {
            if(bPreInsert){
                hfqSet.add_sketch(fqtIdx,sp);
            } else {
                hfqSet.insert(sp,input);
            }
            THEN("The meanQuality for the template is correct") {
                REQUIRE ( hfqSet.templateVec[fqtIdx].meanQual == sp.meanQuality );
            }
            THEN("The number of stored templates and hits increases correctly") {
                REQUIRE( hfqSet.templateVec.size() == startingNTemplate + nTemplateIncrease );
                REQUIRE( hfqSet.sketchMap.hits() == startingNHit + nHit );
            }
            THEN("The stored templates are correct") {
                REQUIRE ( hfqSet.templateVec.front() == input );
            }
            THEN("The input sketch elements can be found") {
                for(const SketchElement & elem : sp.first) {
                    REQUIRE ( hfqSet.sketchMap.count(elem) > 0 );
                    AND_THEN("The ids for the sketch elements are correct") {
                        for(const LocationElement & loc : hfqSet.sketchMap.at(elem)) {
                            REQUIRE( ( loc.id / 2 ) == fqtIdx );
                        }
                    }
                }
            }
        }
    }
}


SCENARIO("Filling from FastqIO source","[SketchHashFactory][FastqIOFill][Filling]") {
    auto nWorker = GENERATE(0,1,2,10);
    INFO("and Given: nWorker is " << nWorker);
    auto phredOffset = GENERATE(0,33,64);
    INFO("and Given: phredOffset is " << phredOffset);
    SketchHashFactory shFactory(Sketcher(7,5,1),nWorker,phredOffset);
    GIVEN("A FastqIO as a source and an empty result") {
        std::shared_ptr<FastqIO> fqIO_ptr;
        auto [nStream, bInterleaved] = GENERATE(table<int,bool>({
                    {1,false}, {1,true}, {2,false} } ));
        INFO("and Given: mode (nStream,bInterleaved) is (" << nStream << ", " << bInterleaved << ")" );
        fqIO_ptr = (nStream == 2) ? std::make_shared<FastqIO>(
                                        constructSS({initFqt1fwd,initFqt2fwd}),
                                        constructSS({initFqt1rev,initFqt2rev}),
                                        FastqIO::IO_IN) :
                                    std::make_shared<FastqIO>(
                                        constructSS({initFqt1pair,initFqt2pair}),
                                        FastqIO::IO_IN, bInterleaved);
        size_t nTemplates = (nStream == 1 && !bInterleaved) ? 4 : 2;
        FastqIOAsFQTSource src(fqIO_ptr);
        HashedFastqSet hfqSet;
        WHEN("Filling is requested") {
            size_t nFill;
            REQUIRE_NOTHROW( nFill = shFactory.FillHashedFastqSet(src, hfqSet) );
            THEN( "The number of templates inserted is correct" ) {
                REQUIRE( nFill == nTemplates );
            }
        }
    }
}


SCENARIO("Filling from a vector source","[SketchHashFactory][VectorFill][Filling]") {
    auto nWorker = GENERATE(0,1,2,10);
    INFO("and Given: nWorker is " << nWorker);
    auto phredOffset = GENERATE(0,33,64);
    INFO("and Given: phredOffset is " << phredOffset);
    SketchHashFactory shFactory(Sketcher(7,5,1),nWorker,phredOffset);
    GIVEN("A vector as a fastq template source and an empty result") {
        std::vector<FastqTemplate_t> fqtVec = {initFqt1pair, initFqt2pair};
        VectorAsFQTSource src(std::make_shared<std::vector<FastqTemplate_t>>(fqtVec));
        HashedFastqSet hfqSet;
        WHEN("Filling is requested") {
            size_t nFill = shFactory.FillHashedFastqSet(src, hfqSet);
            THEN( "The number of templates inserted is correct" ) {
                REQUIRE( nFill == fqtVec.size() );
            }
        }
    }
}
