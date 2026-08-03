#include <iostream>
#include <Sketcher.h>
#include <string>

int main (int argc, char ** argv) {

    std::string op = argv[0];

    //std::cout << "Hello Test World\n";

    Sketcher sketcher(15,5,1);
    Sketch sketch = sketcher.generate_sketch("ACTGACTGGTACGTAACTAGAC");

    if(op == "Syncmer Length"){
        return 0;
    }

    return 0;
}
