#ifndef SKETCHER_HEADER_GAURD_
#define SKETCHER_HEADER_GAURD_

#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>

// Represents a single element in the sketch
struct SketchElement {
    uint64_t hash;
    size_t position;
    std::string to_string() const;
    bool operator==(const SketchElement & other) const {
        return (this->hash == other.hash && this->position == other.position);
    }
};


//The behaviour for Sketches is that the order of SketchElements in the vector
// matches the order of the syncmers in the sequence
using Sketch = std::vector<SketchElement>;

class Sketcher {
public:
    using HashFunction = std::function<uint64_t(const std::string &)>;
    using Alphabet = std::map<char,size_t>;

    // Constructor: Initializes parameters. 
    // alphabet_size defines the valid symbol range if needed, or used for
    // coding. custom_hasher allows injecting a hash function
    // (defaults to std::hash equivalent).
    Sketcher(   size_t k, size_t s, double d,
                size_t alphabet_size = 4,
                HashFunction hasher = Sketcher::FNVHash);
    Sketcher(   size_t k, size_t s, double d, const Alphabet & alpha) :
        Sketcher(k,s,d,alpha.size(),BindLexicographicCoder(alpha)) {};
    Sketcher( const Sketcher & other);

    // Generates a sketch from a string
    //  (accepts rvalue/lvalue, uses move semantics)
    [[nodiscard]] Sketch generate_sketch(const std::string& sequence) const {
        return generate_sketch_impl(sequence);
    }

    [[nodiscard]] Sketch generate_sketch(std::string&& sequence) const {
        // Since the sketch only stores hashes and positions,
        // the underlying string data isn't moved into the sketch,
        // but rvalue overload supports string temporaries.
        return generate_sketch_impl(sequence);
    }

    double ExpectedSyncmerOverlap() const;

    static const HashFunction FNVHash;
    static HashFunction BindLexicographicCoder(const Alphabet & alpha)
    {
        return std::bind(   Sketcher::LexicographicCoder, std::cref(alpha),
                            std::placeholders::_1);
    }
    static uint64_t LexicographicCoder( const Alphabet & alpha,
                                        const std::string & sv);

    static const Alphabet DNA_Alphabet;
    static const Alphabet RNA_Alphabet;

    size_t k() const { return k_; }
    size_t s() const { return s_; }
    double d() const { return double(c_); }

private:
    const size_t k_;
    const size_t s_;
    const size_t w_;
    const size_t c_;
    const size_t alphabet_size_;
    HashFunction hasher_;

    size_t fill_shashVec( const std::string & seq, size_t start, size_t count,
                        std::vector<uint64_t> & vec) const;

    Sketch generate_sketch_impl(const std::string & seq) const;
};


#endif //SKETCHER_HEADER_GAURD_
