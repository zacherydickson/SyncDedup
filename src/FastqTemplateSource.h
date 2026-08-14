#ifndef FASTQ_TEMPLATE_SOURCE_HEADER_GAURD_
#define FASTQ_TEMPLATE_SOURCE_HEADER_GAURD_

#include "FastqIO.h"
#include <memory>
#include <vector>

//An interface defining the concept of a Fastq Template Source
//Any type implementing the concept will be able to act a functor, taking a reference to an fqt
//  it will return true if the result was valid, and false otherwise
//  Can throw std::invalid argument if the source encounters a malformed template 
//get_size will return a true value if the fastq source has a defined type, and false otherwise
//  if the source has s defined type the size argument will be set to that size,
//  otherwise the value of size is undefined
//get_block will allow acquisition of multiple templates at once
//  the size of the returned vector will be up to max_n in size
//  Can throw std::invalid argument if the source encounters a malformed template 
struct FastqTemplateSource {
    virtual bool operator()(FastqTemplate_t & fqt) = 0;
    virtual std::vector<FastqTemplate_t> get_block(size_t max_n) = 0;
    virtual bool get_size(size_t & size) const = 0;
};

//TODO: Add Tests
struct FastqIOAsFQTSource : public FastqTemplateSource {
    FastqIOAsFQTSource(std::shared_ptr<FastqIO> handler);
    bool operator()(FastqTemplate_t & fqt) override;
    std::vector<FastqTemplate_t> get_block(size_t max_n) override;
    bool get_size(size_t & size) const override { return false; }
    protected:
    std::shared_ptr<FastqIO> in;
    FastqIO::READ_RESULT getTemplate(FastqTemplate_t & fqt);
};

//TODO: Add Tests
struct VectorAsFQTSource : public FastqTemplateSource {
    VectorAsFQTSource(std::shared_ptr<std::vector<FastqTemplate_t>> vector);
    bool operator()(FastqTemplate_t & fqt) override;
    std::vector<FastqTemplate_t> get_block(size_t max_n) override;
    bool get_size(size_t & size) const override { size = vec->size(); return true; }
    protected:
    std::shared_ptr<std::vector<FastqTemplate_t>> vec;
    size_t pos;
};


#endif // FASTQ_TEMPLATE_SOURCE_HEADER_GAURD_



