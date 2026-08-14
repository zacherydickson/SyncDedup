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
    bool operator==(const FastqSegment_t & other) const ;
#ifndef NDEBUG 
    std::string to_string() const;
#endif
};

struct FastqTemplate_t {
    std::string name;
    std::vector<FastqSegment_t> segVec;
    bool operator==(const FastqTemplate_t & other) const;
#ifndef NDEBUG 
    std::string to_string() const;
#endif
};

struct FastqStreamPair_t {
    std::unique_ptr<std::iostream> first;
    std::unique_ptr<std::iostream> second;
};

class FastqIO {
    public:
        using istream_ptr = std::unique_ptr<std::istream>;
        using ostream_ptr = std::unique_ptr<std::ostream>;
        using iostream_ptr = std::unique_ptr<std::iostream>;
        enum IO_FLAGS {
            IO_IN = 0x1,
            IO_OUT = 0x2,
            IO_BAD = 0x4,
            IO_INJECTED = 0x8
        };
        enum READ_RESULT {
            READ_PASS = 0,
            READ_EOF = 1,
            READ_MISSING_LEADER1 = 2,
            READ_MISSING_LEADER2 = 3,
            READ_FAIL = 4,
            READ_MISPAIRED = 5,
            READ_SEQ_QUAL_LEN = 6
        };
    protected:
        FastqIO(istream_ptr && in1, istream_ptr && in2,
                ostream_ptr && out1, ostream_ptr && out2, bool bInterleaved,
                bool bInjected = false);
    public:
        FastqIO() = delete;
        FastqIO(FastqIO &) = delete;
        FastqIO(FastqIO &&);
        FastqIO(const std::string & filepath, IO_FLAGS mode,
                bool bInterleaved = false);
        FastqIO(const std::string & filepath1, const std::string filepath2,
                IO_FLAGS mode);
        FastqIO(iostream_ptr && file1, IO_FLAGS mode,
                bool bInterleaved = false);
        FastqIO(iostream_ptr && file1, iostream_ptr && file2,
                IO_FLAGS mode);
        FastqIO(FastqStreamPair_t && sp, IO_FLAGS mode) :
            FastqIO(std::move(sp.first),std::move(sp.second),mode) {}

        READ_RESULT next_template(FastqTemplate_t &);
        bool write(const FastqTemplate_t & fqtemplate);

        bool isGood() const { return !(flags_ & IO_BAD); }
        bool isBad() const { return (flags_ & IO_BAD); }
        bool isWriter() const { return (flags_ & IO_OUT); }
        bool isReader() const { return (flags_ & IO_IN); }
        bool canWrite() const { return (flags_ & (IO_OUT | IO_BAD)) == IO_OUT; }
        bool canRead() const { return (flags_ & (IO_IN | IO_BAD)) == IO_IN; }
        bool fromInjection() const { return flags_ & IO_INJECTED; }
        bool isInterleaved() const { return bInterleaved_; }
        bool isPaired() const { return bPaired_; }

        READ_RESULT skip_templates(size_t n = 1);
        std::pair<size_t,size_t> tell() const;
        //TODO: Add tests
        bool seek(std::pair<size_t,size_t> posPair);
        //NOTE: Untested
        bool seek(  std::pair<int,int> offPair,
                    std::pair<  std::ios_base::seekdir,
                            std::ios_base::seekdir> dirPair);
        //Convenience versions of known single stream handlers
        bool seek(size_t pos);
        bool seek(int off, std::ios_base::seekdir dir);

        FastqStreamPair_t releaseStreams() &&;
        void close() &&;
    protected:
        struct FastqReader {
            READ_RESULT next_template(FastqTemplate_t & fqt); 
            istream_ptr pfile;
        };
        struct FastqWriter {
            bool write(const FastqTemplate_t & fqt);
            ostream_ptr pfile;
        };
        static istream_ptr gzopenpath_in(const std::string & path);
        static ostream_ptr gzopenpath_out(const std::string & path, int compression);
    protected:
        const bool bInterleaved_;
        const bool bPaired_;
        uint8_t flags_;
        FastqReader reader1_;
        FastqReader reader2_;
        FastqWriter writer1_;
        FastqWriter writer2_;
};

#endif// FASTQ_IO_HEADER_GAURD_


