//#include "Sketcher.h"
#include <iostream>
#include <Sketcher.h>

int main (int argc, char ** argv){

    std::cout << "Hello World\n"; 

    Sketcher defaultSketcher(15,5,1.0);
    Sketcher::HashFunction F([](std::string_view a ) {
            return std::hash<std::string_view>{}(a);});
    Sketcher stdhashSketcher(15,5,1.0,4,F);
    std::string seq = "ACTGACTGGATCAGAACAGGG";
    Sketch defaultSketch = defaultSketcher.generate_sketch(seq);
    Sketch stdhashSketch = defaultSketcher.generate_sketch(seq);
    std::cerr << defaultSketch.size() << "\n";

    for(auto & elem : defaultSketch){
        std::cerr << elem.position << "\t";
    }
    std::cerr << "\n";
    for(auto & elem : stdhashSketch){
        std::cerr << elem.position << "\t";
    }
    std::cerr << "\n";


    return 1;
}
