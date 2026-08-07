#include <vector>
#include <string>
#include <cstdint>
#include <stdexcept>
#include <functional>
#include <algorithm>
#include <Sketcher.h>
#include <cmath>

std::string SketchElement:: to_string() const {
    return std::to_string(hash) + ":" + std::to_string(position);
}

Sketcher::Sketcher( size_t k, size_t s, double d,
                    size_t alphabet_size, HashFunction hasher)
    :   k_(k), s_(s), w_(k-s), c_(int(d)),
        alphabet_size_(alphabet_size), hasher_(std::move(hasher))
{
    if (s_ <= 0 || k_ <= 0) {
        throw std::invalid_argument("k- and s-mer length must be greater than zero.");
    }
    if (s_ >= k_) {
        throw std::invalid_argument("s-mer length (s) must be smaller than k-mer length (k).");
    }
    if(d < 1.0) {
        throw std::invalid_argument("downsampling factor (d) must at least 1");
    }
    if(alphabet_size_ < 1) {
        throw std::invalid_argument("Alphabet size must be at least 1.");
    }
    try {
        hasher_(std::string(k_,'A'));
    } catch (std::invalid_argument & e) {
        throw std::invalid_argument("Provided hasher is incompatible with requested kmer length");
    }
}

//Returns the first position in the sequence which is still in the vector
//Assuming the vector has been filled continuously
size_t Sketcher::fill_shashVec(   const std::string & seq, size_t start,
                                size_t count, std::vector<uint64_t> & vec) const
{
    size_t i = start;
    while(i < start + count && i < seq.size()){
        size_t slot = i % vec.size(); 
        size_t hash = hasher_(seq.substr(i,s_));
        vec[slot] = hash;
        i++;
    }
    return (i > vec.size()) ? i - vec.size() : 0;
}

//Generates a sketch for an input sequence, where the sketch is the syncmers for
//  the sequence paired with the positions in the sequence at which those syncmers occur
//  The sketch will be empty if the input is smaller than the value of k_
//  The sketch may be empty if the numer of kmers in the input is smaller than k - s
//Input - a string representing a sequence for which to generate a sketch
//Output - a vector of sketch elements
Sketch Sketcher::generate_sketch_impl(const std::string & seq) const {

    Sketch sketch;
    // A sequence must be at least length k to form a single k-mer
    if (seq.length() < k_) {
        return sketch;
    }
    
    // Bounded positional syncmer logic:
    // A k-mer of length k contains (k - s + 1) s-mers.
    // A window of size (k - s + 1) dictates the span where an s-mer can anchor a bounded syncmer.
    // Let's compute syncmers across valid windows.
    
    // Number of k-mers in the sequence
    size_t num_kmers = seq.length() - k_ + 1;
    //size_t num_windows = seq.length() - w_ + 1;
    
    //Sinlge memory allocation to store hashes of s-mers
    size_t num_smers = k_ - s_ + 1;
    std::vector<uint64_t> shashVec(num_smers,0ULL);
    //Start with the first w_ smers,
    // on update fill the oldest slot and update the index
    for(    size_t i = fill_shashVec(seq,0,num_smers,shashVec);
            i < num_kmers;
            i = fill_shashVec(seq,i+num_smers,1,shashVec) )
    {
        // Find the smallest s-mer value inside this kmer
        // Ties broken towards the leftmost s-mer (min_element finds the first occurrence on ties)
        size_t best_smer_local_idx = 0;
        size_t bestSlot = i % num_smers;
        for (size_t j = 1; j < num_smers; ++j) {
            size_t slot = (i + j) % num_smers;
            auto s = std::string(seq.substr(i+j,s_));
            if (shashVec[slot] < shashVec[bestSlot]) {
                bestSlot = slot;
                best_smer_local_idx = j;
            }
        }

        //Test if the kmer meets the conditions
        //  bounded syncmer condition best smer is at start or end
        if(best_smer_local_idx == 0 || best_smer_local_idx == num_smers -1){
            std::string kmer = seq.substr(i, k_);
            uint64_t code = hasher_(kmer);
            //Downsampling condition
            if(code % c_ == 0) {
                sketch.push_back(SketchElement{.hash = code, .position = i+1 });
            }
        }
    }
    return sketch;
}

//Fowler -Noll-Vo hash function
const Sketcher::HashFunction Sketcher::FNVHash = [](const std::string & sv) {
    uint64_t hash = 14695981039346656037ULL;
    for (char c : sv) {
        hash ^= (uint64_t)(c);
        hash *= 1099511628211ULL;
    }
    return hash;
};

uint64_t Sketcher::LexicographicCoder(  const Alphabet & alpha,
                                        const std::string & sv)
{
    //Each digit can be represented by log(|Σ|) / log(2) bits
    double base = alpha.size();
    double bitsPerResidue = std::log( base ) / std::log( 2.0 );
    size_t bits = std::ceil(bitsPerResidue * sv.length());
    if(bits > 64) {
        throw std::invalid_argument("Lexicographic code longer than 64 bits");
    }
    uint64_t code = 0;
    double power = 0.0;
    for(auto it = sv.rbegin(); it != sv.rend(); it++){
        size_t count = alpha.at(*it);
        code += count * std::pow( base, power++);
    }
    return code;
};

const Sketcher::Alphabet Sketcher::DNA_Alphabet = {{'A',0},{'C',1},{'G',2},{'T',3}};
const Sketcher::Alphabet Sketcher::RNA_Alphabet = {{'A',0},{'C',1},{'G',2},{'U',3}};
