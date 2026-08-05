#include <iostream>
#include <Sketcher.h>
#include <string>
#include <catch2/catch_all.hpp>
#include <Catch2Extensions.hpp>
#include <algorithm>


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

TEST_CASE( "FNV Hash", "[Sketch]") {
    SECTION ("Generates correct codes") {
        REQUIRE(Sketcher::FNVHash("ACGT")       == 0x9a90178ba8feda4e);
        REQUIRE(Sketcher::FNVHash("ACGA")       == 0x9a900c8ba8fec79d);
        REQUIRE(Sketcher::FNVHash("ACGTTGCA")   == 0xcba7310ee49c735b);
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


TEST_CASE( "Correct Sketches can be generated", "[Sketch]") {
    struct testParam_t {
        std::string seq;
        Sketcher sketcher;
        Sketch sketch;
    };
    auto [seq,sketcher,expected] = GENERATE(values<testParam_t>({
       {   "ACTGACTGGA",
           Sketcher(6,4,1.0),
           {   {0xC0DEB678B0D9E8AD,2},
               {0xE61F35FF38889D0D,3},
               {0xAD5B20031CE9F118,4} } },
       {    "ACGATGCTACTTGACGT",
            Sketcher(8,5,1.0,4,[](std::string_view sv){
                    return (size_t) std::hash<std::string_view>{}(sv); }),
            {   {0x3676e1ac56f7b2eb,1},
                {0xc6a6dedf91006e09,2},
                {0xc23740410cd80b2b,5},
                {0x867f1fe436983e24,6},
                {0xc242b6e4fd5ba2a4,9} } },
        {
            "GTCAGTCGTAGCTAGCTGACTGCAT",
            Sketcher(17,5,1.0,Sketcher::DNA_Alphabet),
            {   {0x00000001B2727879,7}
            }
        }
    }) );

    Sketch sketch = sketcher.generate_sketch(seq);
    REQUIRE( std::equal(sketch.begin(),sketch.end(),expected.begin()));
 //   Sketcher sketcher();
 //   std::equal()
 //   static const int nSeq = 2;
 //   std::string seqs[nSeq] = {  "ACTGACTGGATCAGAACAGGGTGAGAT",
 //                               "GTACCGATACGTACGACGTACGTCAGT" };
    //auto i 
    //GIVEN( "An FNV Hash Sketcher") {
    //    Sketcher(15,5,1.0),
    //    WHEN( "Sketches are generated" ) {
    //        Sketch sketches[nSeq];
    //        for(int i = 0; i < nSeq; i++){
    //            sketches[i] = Sketcher.generate_
    //    }
    //}
    //static const int nSketcher = 3;
    //Sketcher sketchers[nSketcher] = {
    //                            Sketcher(15,6,1.0,Sketcher::DNA_Alphabet),
    //                            Sketcher(21,5,1.0,4,stdHasher) };
    //Sketch sketches[nSeq * nSketcher];
    //for(int i = 0; i < nSeq; i++){
    //    for (int j = 0; j < nSketcher; j++){
    //        sketches[i*nSeq+j] = sketchers[j].generate_sketch(seqs[i]);
    //    }
    //}
    //SECTION("Generates correct Sketches") {
    //    REQUIRE(
    //}
    //Sketch sketch = sketcher.generate_sketch(seq);
    //REQUIRE(sketch.size() > 0);
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
