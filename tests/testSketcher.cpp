#include <iostream>
#include <Sketcher.h>
#include <string>
#include <catch2/catch_test_macros.hpp>
#include <Catch2Extensions.hpp>


Sketcher::Alphabet alpha2 = {{'A',0},{'C',1}};
Sketcher::Alphabet alpha3 = {{'A',0},{'C',1},{'G',2}};
Sketcher::Alphabet alpha4 = {{'A',0},{'C',1},{'G',2},{'T',3}};


TEST_CASE( "LexicographicCoding" , "[Sketch]") {
    SECTION ("Throws invalid argument for strings which are too long") {
        REQUIRE_THROWS_AS (
                Sketcher::LexicographicCoder(alpha2,std::string(65,'A')),
                std::invalid_argument);
        REQUIRE_THROWS_AS (
                Sketcher::LexicographicCoder(alpha3,std::string(41,'A')),
                std::invalid_argument);
        REQUIRE_THROWS_AS (
                Sketcher::LexicographicCoder(alpha4,std::string(33,'A')),
                std::invalid_argument);
    }
    SECTION ("Generates correct codes") {
        REQUIRE(Sketcher::LexicographicCoder(alpha4,"ACTG") == 30ULL);
        REQUIRE(Sketcher::LexicographicCoder(alpha4,"TTTT") == 255ULL);
        REQUIRE(Sketcher::LexicographicCoder(alpha4,"ACACGTGT") == 4539ULL);
        REQUIRE(Sketcher::LexicographicCoder(alpha4,"AAAAAAAA") == 0ULL);
    }
}

TEST_CASE( "Sketch object construction", "[Sketch]" ) {
    SECTION ( "No exceptions with valid input" ) {
        REQUIRE_NOTHROW( Sketcher(15,5,1.0) );
        REQUIRE_NOTHROW( Sketcher(15,5,2.5) );
        REQUIRE_NOTHROW( Sketcher(15,5,1.0,1) );
        REQUIRE_NOTHROW( Sketcher(15,5,2.5,4) );
        REQUIRE_NOTHROW( Sketcher(15,5,2.5,4,
                    Sketcher::BindLexicographicCoder(Sketcher::DNA_Alphabet)) );
    }
    SECTION ( "Alphabet Shortcut Constructor" ) {
        REQUIRE_NOTHROW( Sketcher(15,5,2.5,Sketcher::DNA_Alphabet) );
    }
    SECTION ( "Throws invalid arguments for zero k" ) {
        REQUIRE_THROWS_AS( Sketcher(0,5,1.0), std::invalid_argument);
    }
    SECTION ( "Throws invalid arguments for zero s" ) {
        REQUIRE_THROWS_AS( Sketcher(15,0,1.0), std::invalid_argument);
    }
    SECTION ( "Throws invalid arguments for s >= k " ) {
        REQUIRE_THROWS_AS( Sketcher(15,15,1.0), std::invalid_argument);
        REQUIRE_THROWS_AS( Sketcher(15,16,2.5), std::invalid_argument);
    }
    SECTION ( "Throws invalid arguments for d < 1 " ) {
        REQUIRE_THROWS_AS( Sketcher(15,5,0.0), std::invalid_argument);
        REQUIRE_THROWS_AS( Sketcher(15,5,-1.0), std::invalid_argument);
    }
    SECTION ( "Throws invalid arguments for |Σ| < 1 " ) {
        REQUIRE_THROWS_AS( Sketcher(15,5,1.0,0), std::invalid_argument);
    }
    SECTION ( "Throws invalid argument for bad hasher-k combos" ) {
        REQUIRE_THROWS_AS(
                Sketcher(33,5,2.5,Sketcher::DNA_Alphabet),
                std::invalid_argument);
    }
}


TEST_CASE( "Sketch generation", "[Sketch]") {
    Sketcher sketcher(15,5,1.0);
    std::string seq = "ACTGACTGGATCAGAACAGGG";
    Sketch sketch = sketcher.generate_sketch(seq);
    REQUIRE(sketch.size() > 0);
}

//
//TEST_CASE( "Hash function initialization works", "[Sketch]" ) {
//    Sketcher defaultSketcher(15,5,1.0);
//    Sketcher::HashFunction F([](std::string_view a ) {
//            return std::hash<std::string_view>{}(a);});
//    Sketcher stdhashSketcher(15,5,1.0,4,F);
//    std::string seq = "ACTGACTGGATCAGAACAGGG";
//    Sketch defaultSketch = defaultSketcher.generate_sketch(seq);
//    Sketch stdhashSketch = defaultSketcher.generate_sketch(seq);
//
//    for(auto & elem : defaultSketch){
//        std::cout << elem.position << "\t";
//    }
//    std::cerr << "\n";
//    for(auto & elem : stdhashSketch){
//        std::cout << elem.position << "\t";
//    }
//    std::cerr << "\n";
//}
