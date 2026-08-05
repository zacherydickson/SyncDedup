#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <string_view>
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

using Sketch = std::vector<SketchElement>;

class Sketcher {
public:
    using HashFunction = std::function<uint64_t(std::string_view)>;
    using Alphabet = std::map<char,size_t>;

    // Constructor: Initializes parameters. 
    // alphabet_size defines the valid symbol range if needed, or used for
    // coding. custom_hasher allows injecting a hash function
    // (defaults to std::hash equivalent).
    Sketcher(   size_t k, size_t s, double d,
                size_t alphabet_size = 4,
                HashFunction hasher = Sketcher::FNVHash);
    Sketcher(   size_t k, size_t s, double d, Alphabet alpha) :
        Sketcher(k,s,d,alpha.size(),BindLexicographicCoder(alpha)) {};

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

    static const HashFunction FNVHash;
    static HashFunction BindLexicographicCoder(const Alphabet & alpha)
    {
        return std::bind(   Sketcher::LexicographicCoder, std::cref(alpha),
                            std::placeholders::_1);
    }
    static uint64_t LexicographicCoder( const Alphabet & alpha,
                                        std::string_view sv);

    static Alphabet DNA_Alphabet;
    static Alphabet RNA_Alphabet;

private:
    const size_t k_;
    const size_t s_;
    const size_t w_;
    const size_t c_;
    const size_t alphabet_size_;
    HashFunction hasher_;

    void fill_shashVec( std::string_view seq, size_t start, size_t end,
                        std::vector<uint64_t> & vec) const;

    Sketch generate_sketch_impl(std::string_view seq) const;
};
