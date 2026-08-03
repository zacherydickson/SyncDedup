#include <vector>
#include <string>
#include <string_view>
#include <cstdint>
#include <stdexcept>
#include <functional>
#include <algorithm>
#include <Sketcher.h>


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
    
    // Default hash function if none provided (using a basic Fowler–Noll–Vo or std::hash fallback)
    if (!hasher_) {
        hasher_ = [](std::string_view sv) {
            uint64_t hash = 14695981039346656037ULL;
            for (char c : sv) {
                hash ^= (uint64_t)(c);
                hash *= 1099511628211ULL;
            }
            return hash;
        };
    }
}

void Sketcher::fill_shashVec(   std::string_view seq, size_t start, size_t end,
                                std::vector<uint64_t> & vec) const
{
    for(size_t i = start; i <= end && i < seq.size(); i++){
        size_t slot = i % vec.size(); 
        vec[slot] = hasher_(seq.substr(start,s_));
    }
}

Sketch Sketcher::generate_sketch_impl(std::string_view seq) const {
    // A sequence must be at least length k to form a single k-mer
    if (seq.length() < k_) {
        std::string msg =   "Input sequence is too short to generate a sketch"
                            "(length < k).";
        throw std::out_of_range(msg);
    }

    Sketch sketch;
    
    // Bounded positional syncmer logic:
    // A k-mer of length k contains (k - s + 1) s-mers.
    // A window of size (k - s + 1) dictates the span where an s-mer can anchor a bounded syncmer.
    // Let's compute syncmers across valid windows.
    
    // Number of k-mers in the sequence
    size_t num_kmers = seq.length() - k_ + 1;
    //size_t num_windows = seq.length() - w_ + 1;
    
    //Sinlge memory allocation to store hashes of s-mers
    std::vector<uint64_t> shashVec(w_,0ULL);
    //Start with the first w_ smers
    this->fill_shashVec(seq,0,w_-1,shashVec);

    size_t i;
    while(i < num_kmers){
        // Find the smallest s-mer value inside this kmer
        // Ties broken towards the leftmost s-mer (min_element finds the first occurrence on ties)
        size_t best_smer_local_idx = 0;
        size_t bestSlot = i % w_;
        size_t num_smers = k_ - s_ + 1;
        for (size_t j = 1; j < num_smers; ++j) {
            size_t slot = (i + j) % w_;
            if (shashVec[slot] < shashVec[bestSlot]) {
                bestSlot = slot;
                best_smer_local_idx = j;
            }
        }

        //Test if the kmer meets the conditions
        //  bounded syncmer condition best smer is at start or end
        if(best_smer_local_idx > 0 || best_smer_local_idx < num_smers -1){
            //Downsampling condition
            std::string_view kmer = seq.substr(i, k_);
            uint64_t code = hasher_(kmer);
            if(code % c_ > 0) {
                sketch.push_back(SketchElement{.hash = code, .position = i });
            }
        }

        //Determine how far to advance by the location of the best smer
        if(best_smer_local_idx == num_smers - 1){
            // The best smer is at the end, which precludes the next k -s - 1 kmers, so we advance k-s
            this->fill_shashVec(seq,i+1,i+k_-s_,shashVec);
            i += k_ - s_;
        } else { //The case for best at start, or fails the to pass condition
            //Advance one kmer
            i++;
            this->fill_shashVec(seq,i,i,shashVec);
        }
    }

    return sketch;
}

Sketcher::HashFunction Sketcher::LexicographicCoding = [](std::string_view) {
    uint64_t code = 0;
    return code;
};
