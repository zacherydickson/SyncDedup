#include "FastqIO.h"

#include <cassert>
#include <limits>
#include <zlib.h>
#include <zstr.hpp>

FastqIO::FastqIO(   FastqIO::istream_ptr && in1, FastqIO::istream_ptr && in2,
                    FastqIO::ostream_ptr && out1, FastqIO::ostream_ptr && out2,
                    bool bInterleaved) 
    :   bInterleaved_(bInterleaved),
        bPaired_(bInterleaved || (in1 && in2) || (out1 && out2)),
        openmode_(in1 ? IO_IN : IO_OUT),
        reader1_({std::move(in1)}), reader2_({std::move(in1)}),
        writer1_({std::move(out1)}), writer2_({std::move(out1)})
{
#ifndef NDEBUG
    //Things that should be true due to the way the code is written
    size_t inCount = bool(reader1_.pfile) + bool(reader2_.pfile);
    size_t outCount = bool(writer1_.pfile) + bool(writer2_.pfile);
    //Require that only readers or only writers are defined
    assert(inCount * outCount == 0);
    //Require at least one reader or writer is defined
    assert(inCount + outCount > 0);
    //Require that the paired file is only present if the main file is
    if(reader2_.pfile) { assert(inCount == 2); }
    if(writer2_.pfile) { assert(outCount == 2); }
#endif
    //Initial checks to make sure the file can be read from / written to
    if(reader1_.pfile && !(*(reader1_.pfile))) { openmode_ = IO_NONE; }
    if(reader2_.pfile && !(*(reader2_.pfile))) { openmode_ = IO_NONE; }
    if(writer1_.pfile && !(*(writer1_.pfile))) { openmode_ = IO_NONE; }
    if(writer2_.pfile && !(*(writer2_.pfile))) { openmode_ = IO_NONE; }
}


FastqIO::FastqIO(FastqIO && other)
    :   bInterleaved_(other.bInterleaved_),
        bPaired_(other.bPaired_),
        openmode_(other.openmode_),
        reader1_(std::move(other.reader1_)), reader2_(std::move(other.reader2_)),
        writer1_(std::move(other.writer1_)), writer2_(std::move(other.writer2_))
{
    other.openmode_ = (IO_NONE);
}


FastqIO::FastqIO(   const std::string & filepath, bool bInterleaved,
                    IO_MODE om) 
    : FastqIO(  (om == IO_IN) ? FastqIO::gzopenpath_in(filepath) : NULL,
                NULL,
                (om == IO_OUT) ? FastqIO::gzopenpath_out(filepath) : NULL,
                NULL,
                bInterleaved ) {}

FastqIO::FastqIO(const std::string & filepath1, const std::string filepath2,
        IO_MODE om)
    : FastqIO(  (om == IO_IN) ? FastqIO::gzopenpath_in(filepath1) : NULL,
                (om == IO_IN) ? FastqIO::gzopenpath_in(filepath2) : NULL,
                (om == IO_OUT) ? FastqIO::gzopenpath_out(filepath1) : NULL,
                (om == IO_OUT) ? FastqIO::gzopenpath_out(filepath2) : NULL,
                false ) {}

FastqIO::FastqIO(   FastqIO::iostream_ptr && file1, bool bInterleaved,
                    IO_MODE om) 
    : FastqIO(  (om == IO_IN) ? std::move(file1) : NULL,
                NULL,
                (om == IO_OUT) ? std::move(file1) : NULL,
                NULL,
                bInterleaved ) {}

FastqIO::FastqIO(   FastqIO::iostream_ptr && file1, FastqIO::iostream_ptr && file2,
                    IO_MODE om)
    : FastqIO(  (om == IO_IN) ? std::move(file1) : NULL,
                (om == IO_IN) ? std::move(file2) : NULL,
                (om == IO_OUT) ? std::move(file1) : NULL,
                (om == IO_OUT) ? std::move(file2) : NULL,
                false ) {}
    
bool FastqIO::FastqReader::next_template(FastqTemplate_t & fqt) {
    FastqSegment_t & seg = fqt.segVec.emplace_back();
    *pfile >> fqt.name;
    if(fqt.name.length() > 1) {
        fqt.name = fqt.name.substr(1);
        size_t descStart = fqt.name.find(' ');
        if(descStart != std::string::npos){
            seg.desc = fqt.name.substr(descStart+1);
            fqt.name = fqt.name.substr(0,descStart);
        } else {
            seg.desc = "";
        }
    }
    *pfile >> seg.seq;
    pfile->ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    *pfile >> seg.qual;
    if(!(*pfile)) { return false; }
    return true;
}

bool FastqIO::next_template(FastqTemplate_t & fqt) {
    if(openmode_ != IO_IN){
        throw std::logic_error("Attempt to read from a FastqIO object not opened for input");
    }
    openmode_ = IO_NONE;
    if(!reader1_.next_template(fqt)){ return false; }
    if(bPaired_){
        FastqTemplate_t read;
        if(bInterleaved_){
            if(!reader1_.next_template(read)) { return false; }
        } else {
            if(!reader2_.next_template(read)) { return false; }
        }
        if(read.name != fqt.name) {
            throw std::invalid_argument("Paired Reads are not properly paired");
        }
        fqt.segVec.push_back(std::move(read.segVec.front()));
        read.segVec.clear();
    }
    openmode_ = IO_IN;
    return true;
}

bool FastqIO::FastqWriter::write(const FastqTemplate_t & fqt) {
    for(const FastqSegment_t & seg : fqt.segVec) {
        *pfile << '@' << fqt.name;
        if(!pfile) { return false; }
        if(seg.desc.length()){
            *pfile << ' ' << seg.desc;
        }
        *pfile << '\n' << seg.seq << "\n+\n" << seg.qual << '\n';
        if(!pfile) { return false; }
    }
    return true;
}

bool FastqIO::write(const FastqTemplate_t & fqtemplate) {
    if(openmode_ != IO_OUT){
        throw std::logic_error("Attempt to write from a FastqIO object not opened for output");
    }
    openmode_ = IO_NONE;
    if(bInterleaved_ || !bPaired_){
        return writer1_.write(fqtemplate); 
    }
    FastqTemplate_t fqt = {fqtemplate.name,{{}}};
    bool bPass = true;
    for(int i = 0; i < 2; i++){
        fqt.segVec[0] = fqtemplate.segVec[i];
        FastqWriter & writer = (i == 0) ? writer1_ : writer2_;
        bPass &= writer.write(fqt);
    }
    if(bPass) { openmode_ = IO_OUT; }
    return bPass;
}

FastqIO::istream_ptr FastqIO::gzopenpath_in(const std::string & path)
{
    FastqIO::istream_ptr result = NULL;
    if(path == "-") {
        result.reset(new zstr::istream(std::cin));
    } else {
        result.reset(new zstr::ifstream(path));
    }
    return result;
}

FastqIO::ostream_ptr FastqIO::gzopenpath_out(const std::string & path)
{
    FastqIO::ostream_ptr result = NULL;
    if(path == "-") {
        result.reset(new zstr::ostream(std::cout));
    } else {
        result.reset(new zstr::ofstream(path));
    }
    return result;
}
