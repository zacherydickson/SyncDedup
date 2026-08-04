#include <iostream>
#include <Sketcher.h>
#include <string>

int main (int argc, char ** argv){

    //std::cerr << std::string('A',15);
    Sketcher defaultSketcher(15,5,1.0);
    //Sketcher lexicoSketcher(15,5,1.0,4,
    //        Sketcher::BindLexicographicCoder(Sketcher::DNA_Alphabet));
    //std::cout << "Hello World\n"; 
    std::string seq = "ACTGACTGGATCAGAACAGGG";
    Sketch sketch = defaultSketcher.generate_sketch(seq);
    if(sketch.size() > 0) {
        std::cerr << "HUHN\n";
    }
    //std::cerr << sketch.size() << "\n";
    //for(auto & elem : sketch){
    //    std::cerr << elem.hash << ":" << elem.position << "\t";
    //}
    //std::cerr << "\n";
    //
    //Sketcher::HashFunction F([](std::string_view a ) {
    //        return std::hash<std::string_view>{}(a);});
    //Sketcher stdhashSketcher(15,5,1.0,0,F);
    //std::string seq = "ACTGACTGGATCAGAACAGGG";
    //Sketch defaultSketch = defaultSketcher.generate_sketch(seq);
    //Sketch stdhashSketch = defaultSketcher.generate_sketch(seq);
    //std::cerr << defaultSketch.size() << "\n";

    //for(auto & elem : defaultSketch){
    //    std::cerr << elem.position << "\t";
    //}
    //std::cerr << "\n";
    //for(auto & elem : stdhashSketch){
    //    std::cerr << elem.position << "\t";
    //}
    //std::cerr << "\n";


    return 1;
}
