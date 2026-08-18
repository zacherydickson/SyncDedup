#include <catch2/catch_all.hpp>
#include <Catch2Extensions.hpp>
#include <FastqIO.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unistd.h>


FastqTemplate_t initFqt1pair = { "Read1",
    {   { "GACTTTTCCAG", "1:N:0:ACTG", "ABCDEFGHIJK" },
        { "CTGGAAAAGTC", "2:N:0:ACTG", "KJIHGFEDCBA" } } };
FastqTemplate_t initFqt1fwd = { initFqt1pair.name,
    { initFqt1pair.segVec[0] } };
FastqTemplate_t initFqt1rev = { initFqt1pair.name,
    { initFqt1pair.segVec[1] } };

FastqTemplate_t initFqt2pair = { "Read2",
    {   { "CCCCCCCCC", "1:N:0:ACCC", "IIIIIIIII" },
        { "GGGGGGGGGG", "2:N:0:ACCC", "[[[[[[[[[[" } } };
FastqTemplate_t initFqt2fwd = { initFqt2pair.name,
    { initFqt2pair.segVec[0] } };
FastqTemplate_t initFqt2rev = { initFqt2pair.name,
    { initFqt2pair.segVec[1] } };

FastqTemplate_t initFqtBadLenpair = { "Read3",
    {   { "ACTGACTG", "1:N:0:ACCC", "III" },
        { "GACTG", "2:N:0:ACCC", "[[[[[[[[[[" } } };
FastqTemplate_t initFqtBadLenfwd = { initFqtBadLenpair.name,
    { initFqtBadLenpair.segVec[0] } };
FastqTemplate_t initFqtBadLenrev = { initFqtBadLenpair.name,
    { initFqtBadLenpair.segVec[1] } };

std::string FakeFileStem1 = "NoT_A_rEaL_F1Le";
std::string FakeFileStem2 = "NoT_A_rEaL_F2Le";

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

class TempFile {
public:
    explicit TempFile(std::string_view contents, std::string ext)
    {
        auto path =
            std::filesystem::temp_directory_path() / ("FileIO-test-XXXXXX" + ext);

        std::string filename = path.string();

        int fd = mkstemps(filename.data(),ext.length());
        if (fd == -1)
            throw std::runtime_error("mkstemps failed");

        ::close(fd);

        path_ = std::move(filename);

        std::ofstream out(path_);
        if (!out)
            throw std::runtime_error("failed to open temp file");

        out << contents;
    }

    ~TempFile()
    {
        std::error_code ec;
        std::filesystem::remove(path_, ec);
    }

    const std::filesystem::path& path() const { return path_; }

private:
    std::filesystem::path path_;
};

TempFile constructTF(std::string ext, std::vector<FastqTemplate_t> initVec = {}, bool postappend=false,
        FastqIO::READ_RESULT res = FastqIO::READ_PASS) {
    auto ss = constructSS(initVec,postappend,res);
    return TempFile(dynamic_cast<std::stringstream*>(ss.get())->str(),ext);
}


TEST_CASE ("FastqSegment equality works","[FastqSegment_t]") {
    for(size_t i = 0; i < initFqt1pair.segVec.size(); i++) {
        for(size_t j = 0; j < initFqt1pair.segVec.size(); j++) {
            INFO("and Given: segment_i ("<< initFqt1pair.segVec[i].to_string()
                    << ") vs segment_j ("<< initFqt1pair.segVec[j].to_string()
                    <<")");
            REQUIRE( (initFqt1pair.segVec[i] == initFqt1pair.segVec[j]) ==
                     (i == j) );
        }
    }
}

TEST_CASE("FastqTemplate equality works","[FastqTemplate_t]") {
    std::vector<FastqTemplate_t> vec{initFqt1pair,initFqt1fwd,initFqt2pair};
    for(size_t i = 0; i < vec.size(); i++) {
        for(size_t j = 0; j < vec.size(); j++) {
            INFO("and Given: template_i (" << vec[i].to_string() << 
                    ") vs template_j (" << vec[j].to_string() << ")");
            REQUIRE( (vec[i] == vec[j]) == (i == j) );
        }
    }
}

TEST_CASE("FastqTemplate length works","[FastqTemplate_t]") {
    std::vector<FastqTemplate_t> vec{initFqt1pair,initFqt1fwd,initFqt1rev,
                                     initFqt2pair,initFqt2fwd,initFqt2rev};
    std::vector<size_t> exp{22,11,11,19,9,10};
    auto idx = GENERATE(0,1,2,3,4,5);
    INFO("and Given: template is " << vec[idx].to_string());
    REQUIRE( vec[idx].length() == exp[idx] );
}

SCENARIO ("Object Construction from Injection", "[FastqIO][Construction]") {
    auto mode = GENERATE(FastqIO::IO_IN, FastqIO::IO_OUT);
    INFO( "Given: mode is " << ( (mode == FastqIO::IO_IN) ? "IO_IN" : "IO_OUT"));
    AND_GIVEN("A single Stream") {
        auto bInterleaved = GENERATE(false, true);
        INFO( "and Given: Interleaved is " << bInterleaved );
        FastqIO handler(constructSS(),mode,bInterleaved);
        REQUIRE( handler.isGood() );
        REQUIRE( !handler.isBad() );
        REQUIRE( handler.isWriter() == (mode == FastqIO::IO_OUT) );
        REQUIRE( handler.isReader() == (mode == FastqIO::IO_IN) );
        REQUIRE( handler.canWrite() == (mode == FastqIO::IO_OUT) );
        REQUIRE( handler.canRead() == (mode == FastqIO::IO_IN) );
        REQUIRE( handler.fromInjection() );
        REQUIRE( handler.isInterleaved() == bInterleaved );
        REQUIRE( handler.isPaired() == bInterleaved );
    }
    AND_GIVEN("A Pair of streams") {
        FastqIO handler(constructSS(),constructSS(),mode);
        THEN("The state is correct"){
            REQUIRE( handler.isGood() );
            REQUIRE( !handler.isBad() );
            REQUIRE( handler.isWriter() == (mode == FastqIO::IO_OUT) );
            REQUIRE( handler.isReader() == (mode == FastqIO::IO_IN) );
            REQUIRE( handler.canWrite() == (mode == FastqIO::IO_OUT) );
            REQUIRE( handler.canRead() == (mode == FastqIO::IO_IN) );
            REQUIRE( handler.fromInjection() );
            REQUIRE( !handler.isInterleaved() );
            REQUIRE( handler.isPaired() );
        }
    }
}

SCENARIO ("Input handler construction from non-existent files", "[FastqIO][Construction]") {
    auto ext = GENERATE(".fq",".fq.gz");
    INFO( "and Given: ext is " << ext );
    TempFile tf = constructTF(ext);
    GIVEN("Construction with a single input file") {
        bool bInterleaved = GENERATE(true,false);
        INFO("And Given: bInterleaved is " << bInterleaved);
        THEN("Object construction throws invalid argument") {
            REQUIRE_THROWS_AS(
                    FastqIO(FakeFileStem1 + ext,FastqIO::IO_IN, bInterleaved) ,
                std::invalid_argument );
        }
    }
    GIVEN("Construction with two input files"){
        auto missingFlag = GENERATE(0b01,0b10,0b11);
        INFO( "and Given: missingFlag is " << missingFlag );
        std::string path1 = (missingFlag & 0b10) ? FakeFileStem1 + ext : tf.path().string();
        std::string path2 = (missingFlag & 0b01) ? FakeFileStem2 +ext : tf.path().string();
        THEN("Object construction throws invalid argument") {
            REQUIRE_THROWS_AS(  FastqIO(path1,path2,FastqIO::IO_IN) ,
                                std::invalid_argument );
        }
    }
}

SCENARIO("Handler Construction with valid files", "[FastqIO][Construction]" ) {
    auto ext = GENERATE(".fq",".fq.gz");
    INFO( "and Given: ext is " << ext );
    auto mode = GENERATE(FastqIO::IO_IN, FastqIO::IO_OUT);
    INFO("and Given: mode is " << ( (mode == FastqIO::IO_IN) ? "IO_IN" : "IO_OUT"));
    GIVEN("Construction with a single file") {
        bool bInterleaved = GENERATE(true,false);
        INFO("And Given: bInterleaved is " << bInterleaved);
        TempFile tf = constructTF(ext,{initFqt1pair});
        std::string path = tf.path().string();
        THEN("The object can be constructed") {
            std::unique_ptr<FastqIO> handler_ptr;
            REQUIRE_NOTHROW( handler_ptr = std::make_unique<FastqIO>(path,mode,bInterleaved) );
            AND_THEN("The Handler is in the correct state") {
                REQUIRE( handler_ptr->isGood() );
                REQUIRE( !handler_ptr->isBad() );
                REQUIRE( handler_ptr->isWriter() == (mode == FastqIO::IO_OUT) );
                REQUIRE( handler_ptr->isReader() == (mode == FastqIO::IO_IN) );
                REQUIRE( handler_ptr->canWrite() == (mode == FastqIO::IO_OUT) );
                REQUIRE( handler_ptr->canRead() == (mode == FastqIO::IO_IN) );
                REQUIRE( !handler_ptr->fromInjection() );
                REQUIRE( handler_ptr->isInterleaved() == bInterleaved );
                REQUIRE( handler_ptr->isPaired() == bInterleaved );
            }
        }
    }
    GIVEN("Construction with a pair of files") {
        TempFile tf1 = constructTF(ext,{initFqt1fwd});
        TempFile tf2 = constructTF(ext,{initFqt1rev});
        std::string path1 = tf1.path().string();
        std::string path2 = tf2.path().string();
        THEN("The object can be constructed") {
            std::unique_ptr<FastqIO> handler_ptr;
            REQUIRE_NOTHROW( handler_ptr = std::make_unique<FastqIO>(path1,path2,mode) );
            AND_THEN("The Handler is in the correct state") {
                REQUIRE( handler_ptr->isGood() );
                REQUIRE( !handler_ptr->isBad() );
                REQUIRE( handler_ptr->isWriter() == (mode == FastqIO::IO_OUT) );
                REQUIRE( handler_ptr->isReader() == (mode == FastqIO::IO_IN) );
                REQUIRE( handler_ptr->canWrite() == (mode == FastqIO::IO_OUT) );
                REQUIRE( handler_ptr->canRead() == (mode == FastqIO::IO_IN) );
                REQUIRE( !handler_ptr->fromInjection() );
                REQUIRE( !handler_ptr->isInterleaved() );
                REQUIRE( handler_ptr->isPaired() );
            }
        }
    }
}

//NOTE: It is assumed that, regardless of whether a handler was created from a file or
// from injection, its internal operations are consistent (and its easier to test injection)

SCENARIO ("Move Construction", "[MoveConstruction][FastqIO][Construction]") {
    auto mode = GENERATE(FastqIO::IO_IN, FastqIO::IO_OUT);
    INFO("Given: mode is " << ( (mode == FastqIO::IO_IN) ? "IO_IN" : "IO_OUT"));
    GIVEN("A valid handler") {
        FastqIO handler(constructSS({initFqt1fwd},mode == FastqIO::IO_OUT),mode);
        REQUIRE( handler.canWrite() == (mode == FastqIO::IO_OUT) );
        REQUIRE( handler.canRead() == (mode == FastqIO::IO_IN));
        bool fromInjection = handler.fromInjection();
        bool isInterleaved = handler.isInterleaved();
        bool isPaired = handler.isPaired();
        WHEN("Another handler is move constructed from the first") {
            FastqIO movedToHandler(std::move(handler));
            THEN("The original handler is in a bad state") {
                REQUIRE( handler.isBad() );
            }
            THEN("The new handler is in an appropriate valid state") {
                REQUIRE( movedToHandler.isGood() );
                REQUIRE( !movedToHandler.isBad() );
                REQUIRE( movedToHandler.canWrite() == (mode == FastqIO::IO_OUT) );
                REQUIRE( movedToHandler.canRead() == (mode == FastqIO::IO_IN));
                REQUIRE( movedToHandler.fromInjection()  == fromInjection );
                REQUIRE( movedToHandler.isInterleaved() == isInterleaved );
                REQUIRE( movedToHandler.isPaired() == isPaired );
            }
        }
    }
}

SCENARIO ("Inappropriate IO fails", "[FastqIO][Reading][Writing]" ) {
    auto mode = GENERATE(FastqIO::IO_IN, FastqIO::IO_OUT);
    INFO("Given: mode is " << ( (mode == FastqIO::IO_IN) ? "IO_IN" : "IO_OUT"));
    FastqIO handler(constructSS(),mode);
    REQUIRE( handler.canWrite() == (mode == FastqIO::IO_OUT) );
    REQUIRE( handler.canRead() == (mode == FastqIO::IO_IN));
    FastqTemplate_t fqt;
    THEN("Requesting the wrong operation throws") {
        if(mode == FastqIO::IO_IN) {
            REQUIRE_THROWS_AS( handler.write(fqt), std::logic_error );
        } else if(mode == FastqIO::IO_OUT) {
            REQUIRE_THROWS_AS( handler.next_template(fqt), std::logic_error );
        }
    }
    THEN("An IO handler in a bad state") {
        FastqIO fineHandler(std::move(handler));
        if(mode == FastqIO::IO_IN) {
            REQUIRE_THROWS_AS( handler.next_template(fqt), std::logic_error );
        } else if(mode == FastqIO::IO_OUT) {
            REQUIRE_THROWS_AS( handler.write(fqt), std::logic_error );
        }
    }
}

SCENARIO ("Writing to injected streams", "[FastqIO][Writing]") {
    bool bEmpty = GENERATE(true, false);
    INFO("Given: Input streams are " << (bEmpty ? "Empty" : "Non-Empty") );
    GIVEN("An unpaired output handler"){
        std::unique_ptr<std::iostream> ss = (bEmpty) ?
            constructSS() :
            constructSS({initFqt1fwd},true);
        FastqIO out(std::move(ss), FastqIO::IO_OUT);
        REQUIRE( out.canWrite() );
        WHEN("Single or multi-segmented templates are written") {
            for(const FastqTemplate_t & fqt : {initFqt1fwd, initFqt1pair} ) {
                bool res = out.write(fqt);
                THEN("The write succeeds") {
                    REQUIRE( res );
                    AND_THEN("The handler is still in a good state") {
                        REQUIRE( out.isGood() );
                    }
                }
            }
        }
    }
    GIVEN("A paired output handler") {
        bool bInterleaved  = GENERATE(true, false);
        INFO("Given: bInterleaved is " << bInterleaved);
        std::unique_ptr<std::iostream> ss1 = (bEmpty) ?
            constructSS() :
            constructSS({bInterleaved ? initFqt1pair : initFqt1fwd},true);
        std::unique_ptr<std::iostream> ss2 = (bEmpty) ?
            constructSS() :
            constructSS({initFqt1rev},true);
        std::unique_ptr<FastqIO> out_ptr( (bInterleaved) ?
                new FastqIO(std::move(ss1), FastqIO::IO_OUT, true) :
                new FastqIO(std::move(ss1),std::move(ss2), FastqIO::IO_OUT) );
        REQUIRE( out_ptr->canWrite() );
        WHEN("A paired template is written to the handler") {
            bool res = out_ptr->write(initFqt1pair);
            THEN("The write succeeds") {
                REQUIRE( res );
                AND_THEN("The handler is still in a good state") {
                    REQUIRE( out_ptr->isGood() );
                }
            }
        }
        THEN("Attempting to write an unpaired template throws invalid::argument")
        {
            REQUIRE_THROWS_AS( out_ptr->write(initFqt1fwd), std::invalid_argument );
            AND_THEN("The handler is still in a good state") {
                    REQUIRE( out_ptr->isGood() );
            }
        }
    }
}

SCENARIO ("Reading from empty injected streams", "[FastqIO][Reading]") {
    auto [nStream, bInterleaved] = GENERATE(table<int,bool>({
                {1,false}, {1,true}, {2,false} } ));
    INFO("and Given: mode (nStream,bInterleaved) is (" << nStream << ", " << bInterleaved << ")" );
    std::unique_ptr<FastqIO> in_ptr ( (nStream == 1) ?
        new FastqIO(constructSS(),FastqIO::IO_IN,bInterleaved) :
        new FastqIO(constructSS(),constructSS(),FastqIO::IO_IN) );
    REQUIRE( in_ptr->canRead() );
    THEN("Attempting to read fails") {
        FastqTemplate_t fqt;
        REQUIRE( in_ptr->next_template(fqt) != FastqIO::READ_PASS );
        AND_THEN ("The hander is in a bad state") {
            REQUIRE ( in_ptr->isBad() );
        }
    }
}

SCENARIO ("Reading from malformed input", "[FastqIO][Reading]") {
    GIVEN("A paired stream with mismatched input") {
        auto bInterleaved = GENERATE(true, false);
        INFO("And Given: bInterleaved is " << bInterleaved);
        std::unique_ptr<FastqIO> in_ptr( (bInterleaved) ?
                    new FastqIO(constructSS({initFqt1fwd,initFqt2rev}),FastqIO::IO_IN,true) :
                    new FastqIO(constructSS({initFqt1fwd}),constructSS({initFqt2rev}),FastqIO::IO_IN) );
        WHEN("Reading is attempted") {
            FastqTemplate_t fqt;
            FastqIO::READ_RESULT res = in_ptr->next_template(fqt);
            THEN("Reading fails with MISPAIRED") {
                REQUIRE ( res != FastqIO::READ_PASS );
                REQUIRE ( res == FastqIO::READ_MISPAIRED );
                AND_THEN("The handler is in a bad state") {
                    REQUIRE ( in_ptr->isBad() );
                }
            }
        }
    }
    GIVEN("A stream with mismatched sequence lengths") {
        auto [nStream, bInterleaved] = GENERATE(table<int,bool>({
                    {1,false}, {1,true}, {2,false} } ));
        INFO("and Given: mode (nStream,bInterleaved) is (" << nStream << ", " << bInterleaved << ")" );
        std::unique_ptr<FastqIO> in_ptr ( (nStream == 1) ?
            new FastqIO(constructSS({initFqtBadLenpair}),FastqIO::IO_IN,bInterleaved) :
            new FastqIO(constructSS({initFqtBadLenfwd}),constructSS({initFqtBadLenrev}),FastqIO::IO_IN) );
        REQUIRE( in_ptr->canRead() );
        WHEN("Reading is attempted") {
            FastqTemplate_t fqt;
            FastqIO::READ_RESULT res = in_ptr->next_template(fqt);
            THEN("Reading fails with mismatch seq len") {
                REQUIRE( res != FastqIO::READ_PASS );
                REQUIRE( res == FastqIO::READ_SEQ_QUAL_LEN );
                AND_THEN ("The hander is in a bad state") {
                    REQUIRE ( in_ptr->isBad() );
                }
            }
        }
    }
    auto expectedFailure = GENERATE(FastqIO::READ_MISSING_LEADER1,
                                        FastqIO::READ_MISSING_LEADER2,
                                        FastqIO::READ_EOF );
    std::string failureStr;
    switch(expectedFailure) {
        case FastqIO::READ_MISSING_LEADER1:
            failureStr = "MISSING '@'"; break;
        case FastqIO::READ_MISSING_LEADER2:
            failureStr = "MISSING '+'"; break;
        case FastqIO::READ_EOF:
            failureStr = "End of File"; break;
        default:
            failureStr = "No Error"; break;
    }
    INFO("Given: Expected failure is " << failureStr);
    GIVEN("An unpaired stream with partial input") {
        FastqIO in(constructSS({initFqt1fwd},false,expectedFailure),FastqIO::IO_IN);
        REQUIRE( in.canRead() );
        WHEN("Reading is attempted") {
            FastqTemplate_t fqt;
            FastqIO::READ_RESULT res = in.next_template(fqt);
            THEN("Reading fails") {
                REQUIRE ( res != FastqIO::READ_PASS );
                REQUIRE ( res == expectedFailure );
                AND_THEN("The handler is in a bad state") {
                    REQUIRE ( in.isBad() );
                }
            }
        }
    }
    GIVEN("A paired stream with partial input") {
        auto failingReadFlag = GENERATE(0b01,0b01,0b11);
        INFO("and Given: Failing read(s) are " << ((failingReadFlag == 0b11) ? "both" : std::to_string(failingReadFlag)));
        auto bInterleaved = GENERATE(true, false);
        INFO("And Given: bInterleaved is " << bInterleaved);
        std::unique_ptr<std::iostream> ss1 =
            constructSS({initFqt1fwd},false,
                        (failingReadFlag & 0b10) ?  expectedFailure :
                                                    FastqIO::READ_PASS);
        std::unique_ptr<std::iostream> ss2 =
            constructSS({initFqt1rev},false,
                        (failingReadFlag & 0b01) ?  expectedFailure :
                                                    FastqIO::READ_PASS);
        if(bInterleaved && !((failingReadFlag == 0b11) && (expectedFailure == FastqIO::READ_EOF))) {
            *ss1 << ss2->rdbuf();
        }
        std::unique_ptr<FastqIO> in_ptr( (bInterleaved) ?
                    new FastqIO(std::move(ss1),FastqIO::IO_IN,true) :
                    new FastqIO(std::move(ss1),std::move(ss2),FastqIO::IO_IN) );
        REQUIRE( in_ptr->canRead() );
        WHEN("Reading is attempted") {
            FastqTemplate_t fqt;
            FastqIO::READ_RESULT res = in_ptr->next_template(fqt);
            THEN("Reading fails") {
                REQUIRE ( res != FastqIO::READ_PASS );
                REQUIRE ( res == expectedFailure );
                AND_THEN("The handler is in a bad state") {
                    REQUIRE ( in_ptr->isBad() );
                }
            }
        }
    }
}

SCENARIO ("Reading from valid input", "[FastqIO][Reading]") {
    auto [nStream, bInterleaved] = GENERATE(table<int,bool>({
                {1,false}, {1,true}, {2,false} } ));
    INFO("Given: mode (nStream,bInterleaved) is (" << nStream << ", " << bInterleaved << ")" );
    GIVEN("(A) stream(s) which contain valid templates") {
        std::unique_ptr<FastqIO> in_ptr ( (nStream == 1) ?
            new FastqIO(constructSS({initFqt1pair}),FastqIO::IO_IN,bInterleaved) :
            new FastqIO(constructSS({initFqt1fwd}),constructSS({initFqt1rev}),FastqIO::IO_IN) );
        REQUIRE( in_ptr->canRead() );
        WHEN("A template is read from the handler") {
            FastqTemplate_t fqt;
            FastqIO::READ_RESULT res = in_ptr->next_template(fqt);
            INFO("the read is " << fqt.to_string() << "\n");
            THEN("The read succeeds") {
                REQUIRE( res == FastqIO::READ_PASS );
                AND_THEN("The handler is still in a good state") {
                    REQUIRE( in_ptr->isGood() );
                }
                AND_THEN("The read is correct") {
                    REQUIRE ( fqt == ((in_ptr->isPaired()) ? initFqt1pair : initFqt1fwd) );
                }
            }
        }
    }
}

SCENARIO ( "Releasing Injected Streams","[FastqIO][Injection]") {
    auto mode = GENERATE(FastqIO::IO_OUT, FastqIO::IO_IN);
    INFO("Given: mode is " << mode);
    GIVEN("A single stream handler") {
        bool bInterleaved = GENERATE(true, false);
        INFO("And Given: bInterleaved is " << bInterleaved);
        FastqIO handler(constructSS({}), mode, bInterleaved);
        WHEN("The stream is released") {
            FastqStreamPair_t sp = std::move(handler).releaseStreams();
            THEN("Only the fwd stream is defined") {
                REQUIRE ( bool(sp.first) );
                REQUIRE ( ! bool(sp.second) );
            }
            THEN("The handler is in a bad state") {
                REQUIRE ( handler.isBad() );
            }
        }
    }
    GIVEN("A two stream handler") {
        FastqIO handler(constructSS({}),constructSS({}), mode);
        WHEN("The stream is released") {
            FastqStreamPair_t sp = std::move(handler).releaseStreams();
            THEN("Both streams are defined") {
                REQUIRE ( bool(sp.first) );
                REQUIRE ( bool(sp.second) );
            }
            THEN("The handler is in a bad state") {
                REQUIRE ( handler.isBad() );
            }
        }
    }
}

SCENARIO ("Round Trip with injected streams", "[FastqIO][Reading][Writing]") {
    auto [nStream, bInterleaved] = GENERATE(table<int,bool>({
                {1,false}, {1,true}, {2,false} } ));
    INFO("Given: mode (nStream,bInterleaved) is (" << nStream << ", " << bInterleaved << ")" );
    GIVEN("A valid output handler") {
        std::unique_ptr<FastqIO> out_ptr ( (nStream == 1) ?
            new FastqIO(constructSS(),FastqIO::IO_OUT,bInterleaved) :
            new FastqIO(constructSS(),constructSS(),FastqIO::IO_OUT) );
        REQUIRE( out_ptr->canWrite() );
        WHEN("The template is successfully written to an output handler"){
            REQUIRE ( out_ptr->write(initFqt1pair) );
            AND_WHEN("An input handler is constructed from the released stream"){
                FastqStreamPair_t sp = std::move(*out_ptr).releaseStreams();
                REQUIRE( sp.first );
                REQUIRE( bool(sp.second) == (nStream > 1) );
                sp.first->seekg(0);
                INFO("Written stream1 is: '\n" << dynamic_cast<std::stringstream*>(sp.first.get())->str() << "\n'");
                if(sp.second){ 
                    sp.second->seekg(0);
                    INFO("Written stream2 is: '\n" << dynamic_cast<std::stringstream*>(sp.second.get())->str() << "\n'");
                }
                std::unique_ptr<FastqIO> in_ptr ( (nStream == 1) ?
                    new FastqIO(std::move(sp.first),FastqIO::IO_IN,bInterleaved) :
                    new FastqIO(std::move(sp.first),std::move(sp.second),FastqIO::IO_IN));
                REQUIRE(in_ptr->canRead());
                AND_WHEN("A template is successfully read back") {
                    FastqTemplate_t fqt;
                    REQUIRE( in_ptr->next_template(fqt) == FastqIO::READ_PASS );
                    INFO("the read is " << fqt.to_string() << "\n");
                    THEN("The template matches the original input") {
                        REQUIRE ( fqt == ( (in_ptr->isPaired()) ? initFqt1pair : initFqt1fwd) );
                    }
                }
            }
        }
    }
}

SCENARIO ("Round Trip with temporary files", "[FastqIO][Reading][Writing]") {
    auto ext = GENERATE(".fq",".fq.gz");
    INFO( "and Given: ext is " << ext );
    auto [nStream, bInterleaved] = GENERATE(table<int,bool>({
                {1,false}, {1,true}, {2,false} } ));
    INFO("and Given: mode (nStream,bInterleaved) is (" << 
            nStream << ", " << bInterleaved << ")" );
    TempFile tf1 = constructTF(ext);
    TempFile tf2 = constructTF(ext);
    GIVEN("A valid output handler") {
        std::unique_ptr<FastqIO> out_ptr ( (nStream == 1) ?
            new FastqIO(tf1.path().string(),FastqIO::IO_OUT,bInterleaved) :
            new FastqIO(tf1.path().string(),tf2.path().string(),FastqIO::IO_OUT) );
        REQUIRE( out_ptr->canWrite() );
        WHEN("The template is successfully written to an output handler"){
            REQUIRE ( out_ptr->write(initFqt1pair) );
            std::move(*out_ptr).close(); //Explicitly close the output handler prior to reading
            AND_WHEN("An input handler is constructed from the files"){
                std::unique_ptr<FastqIO> in_ptr ( (nStream == 1) ?
                    new FastqIO(tf1.path().string(),FastqIO::IO_IN,bInterleaved) :
                    new FastqIO(tf1.path().string(),tf2.path().string(),FastqIO::IO_IN));
                REQUIRE(in_ptr->canRead());
                AND_WHEN("A template is successfully read back") {
                    FastqTemplate_t fqt;
                    REQUIRE( in_ptr->next_template(fqt) == FastqIO::READ_PASS );
                    INFO("the read is " << fqt.to_string() << "\n");
                    THEN("The template matches the original input") {
                        REQUIRE ( fqt == ( (in_ptr->isPaired()) ? initFqt1pair : initFqt1fwd) );
                    }
                }
            }
        }
    }
}

SCENARIO ("Closing a FastqIO handler","[FastqIO]") {
    auto mode = GENERATE(FastqIO::IO_OUT, FastqIO::IO_IN);
    INFO("Given: mode is " << mode);
    GIVEN("A valid handler") {
        FastqIO handler(constructSS(),mode);
        REQUIRE( handler.canWrite() == (mode == FastqIO::IO_OUT) );
        REQUIRE( handler.canRead() == (mode == FastqIO::IO_IN) );
        WHEN("close() is called") {
            std::move(handler).close();
            THEN("The handler is in a closed state") {
                REQUIRE( !handler.isGood() );
                REQUIRE(  handler.isBad() );
                REQUIRE( !handler.canWrite() );
                REQUIRE( !handler.canRead() );
                REQUIRE( !handler.fromInjection() );
            }
        }
    }
}

SCENARIO ("Compressed vs Uncompressed output","[FastqIO][Writing]") {
    auto [ext,bCompressed] = GENERATE(table<std::string,bool>({
                {".fq",false},{".fq.gz",true},{".gz.fq",false} } ) );
    INFO("and Given: ext is " << ext);
    GIVEN("An output handler") {
        TempFile tf = constructTF(ext);
        FastqIO out(tf.path().string(),FastqIO::IO_OUT);
        REQUIRE( out.canWrite() );
        WHEN("A template has been successfully written to the file") {
            REQUIRE ( out.write(initFqt1fwd) );
            std::move(out).close(); //Explicitly close the output handler prior to reading
            std::ifstream in(tf.path().string());
            THEN("The magic bytes for the output file are correct") {
                REQUIRE( in.get() == (bCompressed ? 0x1F : '@') );
                REQUIRE( in.get() == (bCompressed ? 0x8B : initFqt1fwd.name[0] ) );
            }
        }
    }
}

SCENARIO ("Skipping templates", "[FastqIO][Reading]" ) {
    GIVEN("A valid output handler") {
        FastqIO out(constructSS(),FastqIO::IO_OUT);
        REQUIRE ( out.canWrite() );
        THEN("Attempting to skip templates throws std::logic_error"){
            REQUIRE_THROWS_AS( out.skip_templates(1) , std::logic_error );
        }
    }
    GIVEN("A valid input handler with four templates") {
        std::vector<FastqTemplate_t> input = {initFqt1fwd,initFqt1rev,initFqt2fwd,initFqt2rev};
        auto nSkip = GENERATE(0, 1, 2, 3);
        INFO("and Given: nSkip is" << nSkip);
        FastqIO in(constructSS(input), FastqIO::IO_IN);
        REQUIRE ( in.canRead() );
        WHEN("skip_templates is succesfully called") {
            FastqIO::READ_RESULT res;
            REQUIRE_NOTHROW( res = in.skip_templates(nSkip) );
            THEN("The result is a pass") {
                REQUIRE( res == FastqIO::READ_PASS );
                AND_WHEN("The next template is successfully read") {
                    FastqTemplate_t fqt;
                    res = in.next_template(fqt);
                    REQUIRE (res == FastqIO::READ_PASS );
                    THEN("The template matches the correct template") {
                        REQUIRE( fqt == input[nSkip] );
                    }
                }
            }
        }
    }
    GIVEN("A valid input handler with no templates") {
        FastqIO in(constructSS(), FastqIO::IO_IN);
        REQUIRE ( in.canRead() );
        WHEN("skip_templates is called") {
            FastqIO::READ_RESULT res = in.skip_templates(1);
            THEN("The result is eof"){
                REQUIRE( res == FastqIO::READ_EOF );
            }
        }
    }
    GIVEN("A valid input handler with mismatched qual and seq lengths"){
        FastqIO in(constructSS({initFqtBadLenpair}), FastqIO::IO_IN);
        REQUIRE ( in.canRead() );
        WHEN("skip_templates is called") {
            FastqIO::READ_RESULT res = in.skip_templates(1);
            THEN("The result is Seq len mismatch"){
                REQUIRE( res == FastqIO::READ_SEQ_QUAL_LEN );
            }
        }
    }
    GIVEN("A valid paired handler with mismatched read names"){
        FastqIO in(constructSS({initFqt1fwd,initFqt2rev}), FastqIO::IO_IN, true);
        REQUIRE ( in.canRead() );
        WHEN("skip_templates is called") {
            FastqIO::READ_RESULT res = in.skip_templates(1);
            THEN("The result is Pair mismatch"){
                REQUIRE( res == FastqIO::READ_MISPAIRED );
            }
        }
    }
    auto expectedFailure = GENERATE(FastqIO::READ_MISSING_LEADER1,
                                        FastqIO::READ_MISSING_LEADER2,
                                        FastqIO::READ_EOF );
    std::string failureStr;
    switch(expectedFailure) {
        case FastqIO::READ_MISSING_LEADER1:
            failureStr = "MISSING '@'"; break;
        case FastqIO::READ_MISSING_LEADER2:
            failureStr = "MISSING '+'"; break;
        case FastqIO::READ_EOF:
            failureStr = "End of File"; break;
        default:
            failureStr = "No Error"; break;
    }
    INFO("Given: Expected failure is " << failureStr);
    GIVEN("An unpaired stream with partial input") {
        FastqIO in( constructSS({initFqt1fwd},false,expectedFailure),
                    FastqIO::IO_IN);
        REQUIRE ( in.canRead() );
        WHEN("skip_templates is called") {
            FastqIO::READ_RESULT res = in.skip_templates(1);
            THEN("The result is the expected failure"){
                REQUIRE( res == expectedFailure );
            }
        }
    }
}

SCENARIO("Telling handler positions","[FastqIO][Telling]") {
    auto mode = GENERATE(FastqIO::IO_IN, FastqIO::IO_OUT);
    INFO("and Given: mode is " << ( (mode == FastqIO::IO_IN) ? "IO_IN" : "IO_OUT"));
    auto postappend = GENERATE(true, false);
    INFO("and Given: postappend is " << postappend);
    GIVEN("A handler in a bad state") {
        FastqIO handler(constructSS(),mode);
        FastqIO(std::move(handler));
        REQUIRE( handler.isBad() );
        WHEN("tell() is called") {
            auto posPair = handler.tell();
            THEN("both positions are -1") {
                REQUIRE( posPair.first == -1ULL);
                REQUIRE( posPair.second == -1ULL);
            }
        }
    }
    GIVEN("A valid single stream handler with a non-empty stream") {
        bool bInterleaved = GENERATE(true, false);
        INFO("And Given: bInterleaved is " << bInterleaved);
        auto ss_ptr = constructSS({initFqt1pair},postappend);
        size_t pos = (mode == FastqIO::IO_IN) ? ss_ptr->tellg() : ss_ptr->tellp(); 
        FastqIO handler(std::move(ss_ptr),mode,bInterleaved);
        REQUIRE( handler.canWrite() == (mode == FastqIO::IO_OUT) );
        REQUIRE( handler.canRead() == (mode == FastqIO::IO_IN) );
        WHEN("tell() is called") {
            auto posPair = handler.tell();
            THEN("The first position is correct and the second is -1") {
                REQUIRE( posPair.first == pos);
                REQUIRE( posPair.second == -1ULL);
            }
        }
    }
    GIVEN("A two stream handler") {
        auto ss1_ptr = constructSS({initFqt1fwd},postappend);
        auto ss2_ptr = constructSS({initFqt1rev},postappend);
        size_t pos1 = (mode == FastqIO::IO_IN) ? ss1_ptr->tellg() : ss1_ptr->tellp(); 
        size_t pos2 = (mode == FastqIO::IO_IN) ? ss2_ptr->tellg() : ss2_ptr->tellp(); 
        FastqIO handler(std::move(ss1_ptr),std::move(ss2_ptr),mode);
        REQUIRE( handler.canWrite() == (mode == FastqIO::IO_OUT) );
        REQUIRE( handler.canRead() == (mode == FastqIO::IO_IN) );
        WHEN("tell() is called") {
            auto posPair = handler.tell();
            THEN("Both positions are correct") {
                REQUIRE( posPair.first == pos1);
                REQUIRE( posPair.second == pos2);
            }
        }
    }
}

SCENARIO("Seeking to a valid position changes the results of tell","[FastqIO][Seeking][Telling]") {
    size_t seekPos = constructSS({initFqt1fwd},true)->tellp() / 2;
    if(!seekPos){ seekPos = 1; }
    INFO("and Given: seekPos is " << seekPos);
    auto [nStream, bInterleaved] = GENERATE(table<int,bool>({
                {1,false}, {1,true}, {2,false} } ));
    INFO("and Given: mode (nStream,bInterleaved) is (" << nStream << ", " << bInterleaved << ")" );
    auto mode = GENERATE(FastqIO::IO_IN, FastqIO::IO_OUT);
    INFO("and Given: mode is " << ( (mode == FastqIO::IO_IN) ? "IO_IN" : "IO_OUT"));
    GIVEN("A valid handler with at least a complete read in the stream(s)") {
        std::unique_ptr<FastqIO> handler_ptr ( (nStream == 1) ?
           new FastqIO(constructSS({initFqt1pair}),mode,bInterleaved) :
           new FastqIO(constructSS({initFqt1fwd}),constructSS({initFqt1rev}),mode) );
        REQUIRE( handler_ptr->canWrite() == (mode == FastqIO::IO_OUT) );
        REQUIRE( handler_ptr->canRead() == (mode == FastqIO::IO_IN) );
        auto posPair = handler_ptr->tell();
        REQUIRE ( posPair.first != -1ULL );
        REQUIRE ( (posPair.second == -1ULL) == (nStream == 1) );
        WHEN("seek() is successfully called into a position inside the read") {
            REQUIRE ( handler_ptr->seek({seekPos,seekPos}) );
            auto posPairAfter = handler_ptr->tell();
            THEN("The stream positions change") {
                REQUIRE ( posPair.first != posPairAfter.first );
                if(nStream == 2) {
                    REQUIRE ( posPair.second != posPairAfter.second );
                }
            }
        }
    }
}

SCENARIO("Seeking beyond the end of an input handler","[FastqIO][Seeking]") {
    size_t seekPos = 1ULL + constructSS({initFqt1pair},true)->tellp();
    INFO("and Given: seekPos is " << seekPos);
    auto [nStream, bInterleaved] = GENERATE(table<int,bool>({
                {1,false}, {1,true}, {2,false} } ));
    INFO("and Given: format (nStream,bInterleaved) is (" << nStream << ", " << bInterleaved << ")" );
    GIVEN("A valid input handler with at least one complete segment") {
        std::unique_ptr<FastqIO> in_ptr ( (nStream == 1) ?
           new FastqIO(constructSS({initFqt1pair}),FastqIO::IO_IN,bInterleaved) :
           new FastqIO(constructSS({initFqt1fwd}),constructSS({initFqt1rev}),FastqIO::IO_IN) );
        REQUIRE( in_ptr->canRead() );
            INFO("isReaderPRe is" << in_ptr->isReader() );
        WHEN("seek() is called beyond the stream's end"){
            bool res = in_ptr->seek({seekPos,seekPos});
            INFO("isReader is" << in_ptr->isReader() );
            if( !res) {
                THEN("Failure buts the reader in a bad state") {
                    REQUIRE ( in_ptr->isBad() );
                }
            } else {
                THEN("Attempts to read fail") {
                    FastqTemplate_t fqt;
                    REQUIRE( in_ptr->next_template(fqt) != FastqIO::READ_PASS );
                    REQUIRE( in_ptr->isBad() );
                }
            }
        }
    }
}

SCENARIO("Seeking on input Handlers","[FastqIO][Seeking]") {
    GIVEN("A single streamed input with two paired templates") {
        size_t endOfPair1 = constructSS({initFqt1pair},true)->tellp();
        INFO("and Given: eop is " << endOfPair1 );
        bool bInterleaved = GENERATE(true, false);
        const FastqTemplate_t & targetFqt = (bInterleaved) ? initFqt2pair : initFqt2fwd;
        INFO("and Given: bInterleaved is " << bInterleaved);
        FastqIO in(constructSS({initFqt1pair,initFqt2pair}),FastqIO::IO_IN,
                bInterleaved);
        REQUIRE( in.canRead() );
        bool bConvenience = GENERATE(true, false);
        INFO("and Given: bConvenience is " << bConvenience);
        WHEN("seek() is called to the end of the first pair") {
            bool res = bConvenience ?   in.seek(endOfPair1) :
                                        in.seek({endOfPair1,-1});
            THEN("The seek is successful") {
                REQUIRE ( res );
                REQUIRE ( in.canRead() );
                AND_WHEN("A template is successfully read") {
                    FastqTemplate_t fqt;
                    REQUIRE( in.next_template(fqt) == FastqIO::READ_PASS );
                    INFO("the read is " << fqt.to_string() << "\n");
                    THEN ( "The read is correct" ) {
                        REQUIRE ( fqt == targetFqt );
                    }
                }
            }
        }
    }
    GIVEN("A dual streamed input with two paired templates") {
        size_t endOfFwd2 = constructSS({initFqt2fwd},true)->tellp();
        size_t endOfRev2 = constructSS({initFqt2rev},true)->tellp();
        FastqIO in( constructSS({initFqt2fwd,initFqt1fwd}),
                    constructSS({initFqt2rev,initFqt1rev}),FastqIO::IO_IN );
        REQUIRE( in.canRead() );
        THEN("Convenience seek methods throw") {
            REQUIRE_THROWS_AS( in.seek(0) , std::logic_error );
        }
        WHEN("seek() is called in a proper paired manner") {
            bool res = in.seek({endOfFwd2,endOfRev2});
            THEN("The seek is successful") {
                REQUIRE ( res );
                REQUIRE ( in.canRead() );
                AND_WHEN ("A template is successfully read") {
                    FastqTemplate_t fqt;
                    REQUIRE( in.next_template(fqt) == FastqIO::READ_PASS );
                    INFO("the read is " << fqt.to_string() << "\n");
                    THEN("The template matches the second template") {
                        REQUIRE ( fqt == initFqt1pair );
                    }
                }
            }
        }
        WHEN("seek() is successfully called to mismatched read positions") {
            REQUIRE ( in.seek({endOfFwd2,0}) );
            REQUIRE ( in.canRead() );
            AND_WHEN("An attempt to read is made") {
                FastqTemplate_t fqt;
                FastqIO::READ_RESULT res = in.next_template(fqt);
                THEN("The attempt fails with mismatched pair") {
                    REQUIRE ( res != FastqIO::READ_PASS );
                    REQUIRE ( res == FastqIO::READ_MISPAIRED );
                }
            }
        }
    }
}

