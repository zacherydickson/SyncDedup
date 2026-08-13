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

SCENARIO ("Object Construction from Injection", "[FastqIO][Construction]") {
    auto mode = GENERATE(FastqIO::IO_IN, FastqIO::IO_OUT);
    INFO( "Given: mode is " << ( (mode == FastqIO::IO_IN) ? "IO_IN" : "IO_OUT"));
    AND_GIVEN("A single Stream") {
        auto bInterleaved = GENERATE(false, true);
        INFO( "and Given: Interleaved is " << bInterleaved );
        FastqIO handler(constructSS(),mode,bInterleaved);
        REQUIRE( handler.isGood() );
        REQUIRE( !handler.isBad() );
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
            AND_THEN("The Handler supports the correct operations") {
                REQUIRE( handler_ptr->canWrite() == (mode == FastqIO::IO_OUT) );
                REQUIRE( handler_ptr->canRead() == (mode == FastqIO::IO_IN));
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
            AND_THEN("The Handler supports the correct operations") {
                REQUIRE( handler_ptr->canWrite() == (mode == FastqIO::IO_OUT) );
                REQUIRE( handler_ptr->canRead() == (mode == FastqIO::IO_IN));
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
            out_ptr.reset(NULL);
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
