#include "FastqIO.h"

#include <algorithm>
#include <cassert>
#include <limits>
#include <zlib.h>
#include <zstr.hpp>

bool FastqSegment_t::operator==(const FastqSegment_t & other) const {
    return  ( seq == other.seq ) &&
            ( desc == other.desc ) &&
            ( qual == other.qual );
}

#ifndef NDEBUG
std::string FastqSegment_t::to_string() const {
    return desc + ", " + seq + ", " + qual;
}
#endif

bool FastqTemplate_t::operator==(const FastqTemplate_t & other) const {
    return  ( name == other.name ) &&
            std::equal(segVec.begin(),segVec.end(),other.segVec.begin());
}

#ifndef NDEBUG
std::string FastqTemplate_t::to_string() const {
    std::string str = name + ":{";
    for(const auto & seg : segVec){
        str += "[" + seg.to_string() + "] ";
    }
    return str + "}";
}
#endif

FastqIO::FastqIO(   FastqIO::istream_ptr && in1, FastqIO::istream_ptr && in2,
                    FastqIO::ostream_ptr && out1, FastqIO::ostream_ptr && out2,
                    bool bInterleaved, bool bInjected) 
    :   bInterleaved_(bInterleaved),
        bPaired_(bInterleaved || (in1 && in2) || (out1 && out2)),
        flags_(in1 ? IO_IN : IO_OUT),
        reader1_({std::move(in1)}), reader2_({std::move(in2)}),
        writer1_({std::move(out1)}), writer2_({std::move(out2)})
{
    if(bInjected) {
        flags_ |= IO_INJECTED;
    }
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
    if(reader1_.pfile && !(*(reader1_.pfile))) { flags_ |= IO_BAD; }
    if(reader2_.pfile && !(*(reader2_.pfile))) { flags_ |= IO_BAD; }
    if(writer1_.pfile && !(*(writer1_.pfile))) { flags_ |= IO_BAD; }
    if(writer2_.pfile && !(*(writer2_.pfile))) { flags_ |= IO_BAD; }
}


FastqIO::FastqIO(FastqIO && other)
    :   bInterleaved_(other.bInterleaved_),
        bPaired_(other.bPaired_),
        flags_(other.flags_),
        reader1_(std::move(other.reader1_)), reader2_(std::move(other.reader2_)),
        writer1_(std::move(other.writer1_)), writer2_(std::move(other.writer2_))
{
    other.flags_ |= IO_BAD;
}


FastqIO::FastqIO(   const std::string & filepath, IO_FLAGS om, 
                    bool bInterleaved)
    : FastqIO(  (om == IO_IN) ? FastqIO::gzopenpath_in(filepath) : NULL,
                NULL,
                (om == IO_OUT) ? FastqIO::gzopenpath_out(filepath) : NULL,
                NULL,
                bInterleaved ) {}

FastqIO::FastqIO(const std::string & filepath1, const std::string filepath2,
        IO_FLAGS om)
    : FastqIO(  (om == IO_IN) ? FastqIO::gzopenpath_in(filepath1) : NULL,
                (om == IO_IN) ? FastqIO::gzopenpath_in(filepath2) : NULL,
                (om == IO_OUT) ? FastqIO::gzopenpath_out(filepath1) : NULL,
                (om == IO_OUT) ? FastqIO::gzopenpath_out(filepath2) : NULL,
                false )
{
    if(filepath1 == filepath2) {
        throw std::invalid_argument("Paired fastq files must not be the same file");
    }
}

FastqIO::FastqIO(   FastqIO::iostream_ptr && file1, IO_FLAGS om,
                    bool bInterleaved)
    : FastqIO(  (om == IO_IN) ? std::move(file1) : NULL,
                NULL,
                (om == IO_OUT) ? std::move(file1) : NULL,
                NULL,
                bInterleaved, true ) {}

FastqIO::FastqIO(   FastqIO::iostream_ptr && file1, FastqIO::iostream_ptr && file2,
                    IO_FLAGS om)
    : FastqIO(  (om == IO_IN) ? std::move(file1) : NULL,
                (om == IO_IN) ? std::move(file2) : NULL,
                (om == IO_OUT) ? std::move(file1) : NULL,
                (om == IO_OUT) ? std::move(file2) : NULL,
                false, true) {}
    
FastqIO::READ_RESULT FastqIO::FastqReader::next_template(FastqTemplate_t & fqt) {
    FastqSegment_t & seg = fqt.segVec.emplace_back();
    *pfile >> fqt.name;
#ifndef NDEBUG
    if(!(*pfile)) { return pfile->eof() ? READ_EOF : READ_FAIL; }
#endif
    if(fqt.name.length() > 1) {
        if(fqt.name[0] != '@') { //Check for valid formating
            return READ_MISSING_LEADER1;
        }
        fqt.name = fqt.name.substr(1);
    }
    if(pfile->peek() == ' '){
        *pfile >> seg.desc;
#ifndef NDEBUG
        if(!(*pfile)) { return pfile->eof() ? READ_EOF : READ_FAIL; }
#endif
    } else {
        seg.desc = "";
    }
    *pfile >> seg.seq;
#ifndef NDEBUG
    if(!(*pfile)) { return pfile->eof() ? READ_EOF : READ_FAIL; }
#endif
    //Need to get past the whitespace at the end of the sequence 
    pfile->get();
#ifndef NDEBUG
    if(!(*pfile)) { return pfile->eof() ? READ_EOF : READ_FAIL; }
#endif
    if(pfile->peek() != '+') { // Check for valid formating
        return READ_MISSING_LEADER2;
    }
    pfile->ignore(std::numeric_limits<std::streamsize>::max(),'\n');
#ifndef NDEBUG
    if(!(*pfile)) { return pfile->eof() ? READ_EOF : READ_FAIL; }
#endif
    *pfile >> seg.qual;
    if(!(*pfile)) { return pfile->eof() ? READ_EOF : READ_FAIL; }
    if(seg.qual.length() != seg.seq.length()){
        return READ_SEQ_QUAL_LEN;
    }
    return READ_PASS;
}

FastqIO::READ_RESULT FastqIO::next_template(FastqTemplate_t & fqt) {
    if(!(flags_ & IO_IN)){
        throw std::logic_error("Attempt to read from a FastqIO object not opened for input");
    }
    if(flags_ & IO_BAD) {
        throw std::logic_error("Attempt to read from a FastqIO object in a bad state");
    }
    flags_ |= IO_BAD;
    READ_RESULT res;
    res = reader1_.next_template(fqt);
    if(res != READ_PASS){ return res; }
    if(bPaired_){
        FastqTemplate_t read;
        if(bInterleaved_){
            res = reader1_.next_template(read);
        } else {
            res = reader2_.next_template(read);
        }
        if(res != READ_PASS) { return res; }
        if(read.name != fqt.name) { return READ_MISPAIRED; }
        fqt.segVec.push_back(std::move(read.segVec.front()));
        read.segVec.clear();
    }
    uint8_t mask = ~IO_BAD;
    flags_ &= mask;
    return READ_PASS;
}

bool FastqIO::FastqWriter::write(const FastqTemplate_t & fqt) {
    for(const FastqSegment_t & seg : fqt.segVec) {
        *pfile << '@' << fqt.name;
        if(!pfile) { return false; }
        if(seg.desc.length()){
            *pfile << ' ' << seg.desc;
            if(!pfile) { return false; }
        }
        *pfile << '\n' << seg.seq << "\n+\n" << seg.qual << '\n';
        if(!pfile) { return false; }
    }
    return true;
}

bool FastqIO::write(const FastqTemplate_t & fqtemplate) {
    if(!(flags_ & IO_OUT)){
        throw std::logic_error("Attempt to write to a FastqIO object not opened for output");
    }
    if(flags_ & IO_BAD) {
        throw std::logic_error("Attempt to write to a FastqIO object in a bad state");
    }
    if(bPaired_ && fqtemplate.segVec.size() % 2 != 0 ){
        throw std::invalid_argument("Attempt to paired write an odd number of segments");
    }
    flags_ |= IO_BAD;
    bool bPass = true;
    if(bInterleaved_ || !bPaired_){
        bPass = writer1_.write(fqtemplate); 
    } else {
        FastqTemplate_t fqt = {fqtemplate.name,{{}}};
        for(int i = 0; i < 2; i++){
            fqt.segVec[0] = fqtemplate.segVec[i];
            FastqWriter & writer = (i == 0) ? writer1_ : writer2_;
            bPass &= writer.write(fqt);
        }
    }
    if(bPass) { uint8_t mask = ~IO_BAD; flags_ &= mask; }
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

//If the object was instantiated through stream injection
//the underlying streams can be released as they must have originally been 
//iostream objects
//The strategy is to check for injectable initialization, then attempt to
// dynamic cast back to iostream, if successful then ownership is transferred 
// to the return object
//Throws logic_error if called on a non-injected object
//Throws runtime error if dynamic_cast fails
FastqStreamPair_t FastqIO::releaseStreams() && {
    if(!(flags_ & IO_INJECTED)) {
        throw std::logic_error("Call to releaseStreams for a fastq IO handler which was not constructed through stream injection");
    }
    FastqStreamPair_t streamPair;
    if(flags_ & IO_IN) { //Release the reader streams
        std::iostream * p = dynamic_cast<std::iostream*>(reader1_.pfile.get());
        if(!p){
            throw std::runtime_error("Unable to reverse cast reader1");
        }
        reader1_.pfile.release();
        streamPair.first.reset(p);
        p = NULL;
        if(reader2_.pfile) {
            p = dynamic_cast<std::iostream*>(reader2_.pfile.get());
            if(!p){
                throw std::runtime_error("Unable to reverse cast reader2");
            }
            reader2_.pfile.release();
            streamPair.second.reset(p);
            p = NULL;
        }
    } else { //Release the writer streams
        std::iostream * p = dynamic_cast<std::iostream*>(writer1_.pfile.get());
        if(!p){
            throw std::runtime_error("Unable to reverse cast writer1");
        }
        writer1_.pfile.release();
        streamPair.first.reset(p);
        p = NULL;
        if(writer2_.pfile) {
            p = dynamic_cast<std::iostream*>(writer2_.pfile.get());
            if(!p){
                throw std::runtime_error("Unable to reverse cast writer2");
            }
            writer2_.pfile.release();
            streamPair.second.reset(p);
            p = NULL;
        }
    }
    flags_ |= IO_BAD;
    return streamPair;
}
