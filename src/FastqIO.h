#ifndef FASTQ_IO_HEADER_GAURD_
#define FASTQ_IO_HEADER_GAURD_

#include <string>
#include <iostream>
#include <memory>

struct FastqSegment_t {
    std::string seq;
    std::string qual;
};

struct FastqTemplate_t {
    std::string name;
    FastqSegment_t r1;
    FastqSegment_t * r2;
    bool hasR2() { return r2; }
};

class FastqIO {
    public:
        FastqIO() = delete;
        FastqIO(const std::string & filepath, bool bInterleaved = false,
                std::ios_base::openmode = std::ios_base::in);
        FastqIO(const std::string & filepath1, const std::string filepath2,
                std::ios_base::openmode = std::ios_base::in);
        FastqIO(std::iostream && file1);
        FastqIO(std::iostream && file1, std::iostream && file2);

        FastqTemplate_t next();
        void write(const FastqTemplate_t & fqtemplate);
    protected:
        const bool bInterleaved_;
        std::iostream file1;
        std::unique_ptr<std::iostream> file2;
};

#endif// FASTQ_IO_HEADER_GAURD_


