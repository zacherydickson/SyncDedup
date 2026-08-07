#include <iostream>
#include <Sketcher.h>
#include <string>
#include <catch2/catch_all.hpp>
#include <Catch2Extensions.hpp>
#include <algorithm>
#include <random>
#include <unordered_set>


Sketcher::Alphabet alpha2 = {{'A',0},{'C',1}};
Sketcher::Alphabet alpha3 = {{'A',0},{'C',1},{'G',2}};
Sketcher::Alphabet alpha4 = {{'A',0},{'C',1},{'G',2},{'T',3}};
std::map<uint8_t,char> bit2nuc4 = {{0,'A'},{1,'C'},{2,'G'},{3,'T'}};


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
        {
            "GTCAGTCGTAGCTAGCTGACTGCAT",
            Sketcher(17,5,1.0,Sketcher::DNA_Alphabet),
            {   {0x00000001B2727879,7}
            }
        },
        { // Can't generate a sketch for a equence smaller than k
            "ACTG",
            Sketcher(10,4,1.0),
            {}
        }
    }) );

    INFO("The input sequence is " << seq);
    Sketch sketch = sketcher.generate_sketch(seq);
    REQUIRE (sketch.size() == expected.size());
    REQUIRE( std::equal(sketch.begin(),sketch.end(),expected.begin()));
}

TEST_CASE("Downsampling occurs correctly") {
    std::string seq = "ACTGACTGGA";
    const Sketch Expected =
           {   {0xC0DEB678B0D9E8AD,2},
               {0xE61F35FF38889D0D,3},
               {0xAD5B20031CE9F118,4} };
    auto [dsfactor, expIdxVec] = GENERATE(table<double, std::vector<int>>({
                {1.0, {0,1,2}},
                {2.0, {2}},
                {7.0, {}},
                {19.0, {1}} } ) );
    Sketcher sketcher(6,4,dsfactor);
    Sketch sketch = sketcher.generate_sketch(seq);
    REQUIRE(sketch.size() == expIdxVec.size());
    for(size_t i = 0; i < sketch.size(); i++){
        REQUIRE(sketch[i] == Expected[expIdxVec[i]]);
    }
}

//Assumes that the question is: how much of a does b cover?
double SketchSimilarity(const Sketch & a, const Sketch & b) {
    size_t intersect = 0;
    auto it_a = a.begin();
    auto it_b = b.begin();
    while(it_a != a.end() && it_b != b.end()){
        if(*it_a == *it_b){
            intersect++;
            it_a++;
            it_b++;
        } else if(it_a->position > it_b->position){
            it_b++;
        } else {
            it_a++;
        }
    }
    double ji = double(intersect) / double(a.size());
    return ji;
}

double SketchSimilarity( const Sketch & a, const std::string bStr,
                                const Sketcher & sketcher)
{
    Sketch b = sketcher.generate_sketch(bStr);
    return SketchSimilarity(a,b);
}

std::string RandomDNA(size_t length) {
    std::string seq;
    while(seq.length() < length){
        int64_t byte = rand();
        for(int i =0; i < 16; i++){
            seq += bit2nuc4[byte & 3];
            byte >>= 2;
        }
    }
    return  seq.substr(0,length);
}

std::string CompositionMatchedRandomDNA(const std::string & other) {
    std::string seq(other);
    std::random_device rd;
    std::mt19937 g(rd());
    while(seq == other){
        std::shuffle(seq.begin(),seq.end(),g);
    }
    return seq;
}

//nVar caps out at half the sequence length (rounded up)
std::string MutatedRandomDNA(   const std::string & other, size_t nVar,
                                std::unordered_set<int> * seenPos_pt = NULL )
{
    std::string seq(other);
    std::random_device rd;
    std::mt19937 g(rd());
    std::uniform_int_distribution<> rpos(0,seq.length()-1);
    std::uniform_int_distribution<> rmut(1,3);
    size_t maxVar = std::ceil(other.length() * 0.5);
    if(nVar > maxVar){
        nVar = maxVar;
    }
    //Local default
    std::unordered_set<int> seenPos;
    if(!seenPos_pt){
        seenPos_pt = & seenPos;
    }
    for(size_t i = 0; i < nVar; i++){
        int pos;
        //Select a position which has not yet been changed
        do {
            pos = rpos(g);
        } while(seenPos_pt->count(pos));
        seenPos_pt->insert(pos);
        uint8_t bits = Sketcher::DNA_Alphabet.at(other[pos]);
        char nuc = (bits + rmut(g)) % 4;
        seq[pos] = nuc;
    }
    return seq;
}

double factorial(double n) {
    double res = 1;
    for(size_t i = 2; i <= n; i++){
        res *= i;
    }
    return res;
}

//This probability assumes all characters in l are distinct
//Which operates as a conservative lower bound for the the probability
//with repeated characters
//l is the length of the permuted sequence 
//k is the number of positions which do not change
double ProbPermutedPass(size_t l, size_t k) {
    double p = 0.0;
    for(size_t m = 0; m <= k; m++){
        double derange = 0;
        for(size_t i = 0; i <= m; i++){
            derange += ((i % 2 == 1) ? -1 : 1) / factorial(i);
        }
        p += derange / factorial(l - m);
    }
    return p;
}


TEST_CASE ("Sketches from mutated duplicates are usually more similar than shuffled seqs", "[Sketcher]") {
    auto seqLen = GENERATE(25,take(3, random(50, 150)));
    auto sketcher = GENERATE(Sketcher(13,7,2.0),Sketcher(17,8,1.0));
    std::string templateSeq = RandomDNA(seqLen);
    Sketch templateSketch = sketcher.generate_sketch(templateSeq);

    size_t nScenario = 10;
    double expectedFalseFail = 0;
    size_t nFail = 0;
    for(size_t i = 0; i < nScenario; i++){
        std::string mutatedSeq = templateSeq;
        std::unordered_set<int> seenPos;
        double lastSim = 1.0;
        size_t lastNVar = 0;
        for(double varProp = 0; varProp < 0.5; varProp += 0.05){
            size_t nVar = std::ceil(varProp * seqLen);
            if(nVar == 0) { nVar++; }
            if(nVar <= lastNVar) { continue; }
            double probFalseFail = ProbPermutedPass(seqLen,nVar);
            expectedFalseFail += probFalseFail;
            mutatedSeq = MutatedRandomDNA(mutatedSeq,nVar-lastNVar,&seenPos);
            double mutatedSim = SketchSimilarity( templateSketch,
                                                  mutatedSeq,
                                                  sketcher );
            lastNVar = nVar;
            CHECK(mutatedSim <= lastSim);
            lastSim = mutatedSim;
            double matchedSim = SketchSimilarity(
                                templateSketch,
                                CompositionMatchedRandomDNA(templateSeq),
                                sketcher );
            if(mutatedSim < matchedSim){
                nFail++;
            }
        }
    }
    INFO("The number of failures was " << nFail << " compared to " << expectedFalseFail);
    REQUIRE(nFail < std::ceil(expectedFalseFail)  );
}
