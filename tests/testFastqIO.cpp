#include <FastqIO.h>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <catch2/catch_all.hpp>
#include <Catch2Extensions.hpp>


FastqTemplate_t initFqt1 = { "Read1", 
    { { "GACTTTTCCAG", "1:N:0:ACTG", "ABCDEFGHIJK" } } };
FastqTemplate_t initFqt2 = { "Read2", 
    {   { "GACTTTTCCAG", "1:N:0:ACTG", "ABCDEFGHIJK" }, 
        { "CTGGAAAAGTC", "2:N:0:ACTG", "KJIHGFEDCBA" } } };

std::unique_ptr<std::iostream > constructSS(
        std::vector<FastqTemplate_t> initVec = {}, bool postappend=false)
{
    std::unique_ptr< std::iostream> ss_ptr(new std::stringstream());
    for(const FastqTemplate_t & fqt : initVec ){
        for(const FastqSegment_t & seg : fqt.segVec) {
            *ss_ptr <<  '@' << fqt.name << ' ' << seg.desc << "\n" <<
                        seg.seq << "\n" << "+\n" <<
                        seg.qual << "\n";
        }
    }
    if(!postappend){
        ss_ptr->seekg(0);
    }
    return ss_ptr;
}

SCENARIO ("Object Construction from Injection", "[FastqIO]") {
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

SCENARIO ("Writing fastq entries", "[FastqIO]") {
    GIVEN("A FastqTemplate with one segment and an empty output handler") {
        bool res = true;
        FastqIO out(constructSS({}), FastqIO::IO_OUT);
        REQUIRE( out.canWrite() );
        WHEN("A template is written to the handler") {
            res = out.write(initFqt1);
            THEN("The write succeeds") {
                REQUIRE( res );
                AND_THEN("The handler is still in a good state") {
                    REQUIRE( out.isGood() );
                }
            }
        }
    }
    GIVEN("A FastqTemplate with two segments and an empty paired output handler") {
        bool res = true;
        FastqIO out(constructSS({}),constructSS({}), FastqIO::IO_OUT);
        REQUIRE( out.canWrite() );
        WHEN("A template is written to the handler") {
            res = out.write(initFqt2);
            THEN("The write succeeds") {
                REQUIRE( res );
                AND_THEN("The handler is still in a good state") {
                    REQUIRE( out.isGood() );
                }
            }
        }
    }
    //TODO: Writing in other cases
}

SCENARIO ("Move Construction", "[MoveConstruction][FastqIO]") {
    auto mode = GENERATE(FastqIO::IO_IN, FastqIO::IO_OUT);
    INFO("Given: mode is " << ( (mode == FastqIO::IO_IN) ? "IO_IN" : "IO_OUT"));
    GIVEN("A valid handler") {
        FastqIO handler(constructSS({initFqt1},mode == FastqIO::IO_OUT),mode);
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

SCENARIO ("Inappropriate IO fails", "[FastqIO]" ) {
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

SCENARIO ("Reading from injected Fastq", "[FastqIO]") {
    GIVEN("An empty, unpaired IO handler"){
        FastqIO in(constructSS({}), FastqIO::IO_IN);
        WHEN("A template is read from the handler") {
            FastqTemplate_t fqt; 
            bool res = in.next_template(fqt);
            THEN("Reading fails") {
                REQUIRE ( !res );
                AND_THEN ("The hander is in a bad state") {
                    REQUIRE ( in.isBad() );
                }
            }
        }
    }
    GIVEN("An input handler on a string stream with one unpaired entry") {
        FastqIO in(constructSS({initFqt1}), FastqIO::IO_IN);
        REQUIRE( in.canRead() );
        WHEN("A template is read from the handler") {
            FastqTemplate_t fqt;
            bool res = in.next_template(fqt);
            INFO("the read is " << fqt.to_string() << "\n");
            THEN("The read succeeds") {
                REQUIRE( res );
                AND_THEN("The handler is still in a good state") {
                    REQUIRE( in.isGood() );
                }
                AND_THEN( "The read is correct") {
                    REQUIRE ( fqt == initFqt1 );
                }
            }
        }
    }
    //TODO: Writing in other cases
}

SCENARIO ( "Releasing Injected Streams ") {
    GIVEN("An unpaired output handler") {
        FastqIO out(constructSS({}), FastqIO::IO_OUT);
        WHEN("The stream is released") {
            FastqStreamPair_t sp = std::move(out).releaseStreams();
            THEN("Only the fwd stream is defined") {
                REQUIRE ( bool(sp.first) );
                REQUIRE ( ! bool(sp.second) );
            }
            THEN("The output handler is in a bad state") {
                REQUIRE ( out.isBad() );
            }
        }
    }
    //TODO: Release in other cases
}

SCENARIO ("Writing and Reading one unpaired fastq entry", "[FastqIO]") {
    GIVEN("A FastqTemplate with one segment") {
        WHEN("The template is written to an input handler"){
            FastqIO out(constructSS({}), FastqIO::IO_OUT);
            REQUIRE( out.canWrite() );
            REQUIRE ( out.write(initFqt1) );
            AND_WHEN("An input handler is constructed from the released stream"){
                FastqStreamPair_t sp = std::move(out).releaseStreams();
                REQUIRE( sp.first );
                sp.first->seekg(0);
                FastqIO in(std::move(sp.first),FastqIO::IO_IN);
                REQUIRE(in.canRead());
                AND_WHEN("A template is read back") {
                    FastqTemplate_t fqt;
                    REQUIRE( in.next_template(fqt) );
                    INFO("the read is " << fqt.to_string() << "\n");
                    THEN("The template matches the original input") {
                        REQUIRE ( fqt == initFqt1 );
                    }
                }
            }
        }
    }
}
