#include <catch2/catch_all.hpp>
#include <Catch2Extensions.hpp>
#include <FastqTemplateSource.h>

#include <memory>

enum SOURCE_TYPES {
    FASTQIO_SRC,
    VECTOR_SRC
};

FastqTemplate_t initFqt1pair = { "Read1",
    {   { "GACTTTTCCAG", "1:N:0:ACTG", "ABCDEFGHIJK" },
        { "CTGGAAAAGTC", "2:N:0:ACTG", "KJIHGFEDCBA" } } };
FastqTemplate_t initFqt1fwd = { initFqt1pair.name,
    { initFqt1pair.segVec[0] } };
FastqTemplate_t initFqt1rev = { initFqt1pair.name,
    { initFqt1pair.segVec[1] } };

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

SCENARIO ("FastqIO Construction","[FastqTemplateSource][FastqIOAsFQTSource][Construction]")
{
    GIVEN("A valid input handler") {
        std::shared_ptr<FastqIO> handler_ptr(
                new FastqIO(constructSS({initFqt1pair}),FastqIO::IO_IN));
        REQUIRE ( handler_ptr->canRead() );
        THEN("Construction doesn't throw"){
            REQUIRE_NOTHROW( FastqIOAsFQTSource(handler_ptr) );
        }
    }
    GIVEN("An output handler") {
        std::shared_ptr<FastqIO> handler_ptr(
                new FastqIO(constructSS({initFqt1pair}),FastqIO::IO_OUT));
        THEN("Construction throws") {
            REQUIRE_THROWS_AS(  FastqIOAsFQTSource(handler_ptr) , 
                                std::invalid_argument );
        }
    }
}


SCENARIO("Implemented Source Types have appropriate zie definition states","[FastqTemplateSource]") {
    auto nCopy = GENERATE(0,1,2,3);
    INFO("and Given: nCopy is " << nCopy);
    GIVEN("A FastqIO as a source") {
        std::shared_ptr<FastqIO> handler_ptr(
                new FastqIO(constructSS(std::vector<FastqTemplate_t>(nCopy,initFqt1pair)),FastqIO::IO_IN));
        FastqIOAsFQTSource source(handler_ptr);
        THEN("The source does not have defined size") {
            size_t srcSize;
            REQUIRE( !source.get_size(srcSize) );
        }
    }
    GIVEN("A Vector as a source") {
        VectorAsFQTSource source(std::make_shared<std::vector<FastqTemplate_t>>(nCopy,initFqt1pair));
        THEN("The source has a defined size"){
            size_t srcSize;
            REQUIRE( source.get_size(srcSize) );
            AND_THEN("The size is correct") {
                REQUIRE( srcSize == nCopy );
            }
        }
    }
}

SCENARIO("Fastq Templates can be retreived from sources","[FastqTemplateSource]") {
    auto nCopy = GENERATE(0,1,2,3);
    INFO("and Given: nCopy is " << nCopy);
    std::vector<FastqTemplate_t> dataVec(nCopy,initFqt1fwd);
    GIVEN("A FastqTemplateSource") {
        SOURCE_TYPES testType = GENERATE(FASTQIO_SRC,VECTOR_SRC);
        std::unique_ptr<FastqTemplateSource> src_ptr;
        std::string srcType;
        switch(testType) {
            case FASTQIO_SRC:
                srcType = "FastqIO";
                {
                    std::shared_ptr<FastqIO> handler_ptr(
                        new FastqIO(constructSS(dataVec),FastqIO::IO_IN));
                    src_ptr.reset(new FastqIOAsFQTSource(handler_ptr) );
                }
                break;
            case VECTOR_SRC:
                srcType = "Vector";
                src_ptr.reset(new VectorAsFQTSource(
                    std::make_shared<std::vector<FastqTemplate_t>>(dataVec)));
                break;
        }
        INFO("and Given: Source type is " << srcType);
        WHEN("Used as a functor") {
            FastqTemplate_t fqt;
            THEN("A template can be retreived if available") {
                REQUIRE( (*src_ptr)(fqt) == nCopy > 0 );
            }
        }
        WHEN("A block is retreived") {
            auto blockSize = GENERATE(2,5);
            INFO("and Given: blockSize is " << blockSize);
            std::vector<FastqTemplate_t> block = src_ptr->get_block(blockSize);
            THEN("The block size is appropriate") {
                REQUIRE ( block.size() == std::min(nCopy,blockSize) );
            }
        }
    }
    GIVEN("A Vector as a source") {
    }
}
