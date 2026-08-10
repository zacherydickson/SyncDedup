#ifndef FASTQ_IO_HEADER_GAURD_
#define FASTQ_IO_HEADER_GAURD_

#include <string>
#include <iostream>
#include <memory>
#include <vector>

struct FastqSegment_t {
    std::string seq;
    std::string desc;
    std::string qual;
};

struct FastqTemplate_t {
    std::string name;
    std::vector<FastqSegment_t> segVec;
};

class FastqIO {
    public:
        using istream_ptr = std::unique_ptr<std::istream>;
        using ostream_ptr = std::unique_ptr<std::ostream>;
        using iostream_ptr = std::unique_ptr<std::iostream>;
        enum IO_MODE {
            IO_IN,
            IO_OUT,
            IO_NONE //Used to indicate the Object is in a state where reading/writing is impossible
        };
    protected:
        FastqIO(istream_ptr && in1, istream_ptr && in2,
                ostream_ptr && out1, ostream_ptr && out2, bool bInterleaved);
    public:
        FastqIO() = delete;
        FastqIO(FastqIO &) = delete;
        FastqIO(FastqIO &&);
        FastqIO(const std::string & filepath, bool bInterleaved = false,
                IO_MODE = IO_IN);
        FastqIO(const std::string & filepath1, const std::string filepath2,
                IO_MODE = IO_IN);
        FastqIO(iostream_ptr && file1, bool bInterleaved = false,
                IO_MODE = IO_IN);
        FastqIO(iostream_ptr && file1, iostream_ptr && file2,
                IO_MODE = IO_IN);

        bool next_template(FastqTemplate_t &);
        bool write(const FastqTemplate_t & fqtemplate);
        bool good() { return openmode_ != IO_NONE; } 
        bool bad() { return openmode_ == IO_NONE; } 
        IO_MODE get_mode() { return openmode_; }
        bool isInterleaved() { return bInterleaved_; }
        bool isPaired() { return bPaired_; }
    protected:
        struct FastqReader {
            bool next_template(FastqTemplate_t & fqt); 
            istream_ptr pfile;
        };
        struct FastqWriter {
            bool write(const FastqTemplate_t & fqt);
            ostream_ptr pfile;
        };
        static istream_ptr gzopenpath_in(const std::string & path);
        static ostream_ptr gzopenpath_out(const std::string & path);
    protected:
        const bool bInterleaved_;
        const bool bPaired_;
        IO_MODE openmode_;
        FastqReader reader1_;
        FastqReader reader2_;
        FastqWriter writer1_;
        FastqWriter writer2_;
};

#endif// FASTQ_IO_HEADER_GAURD_


