#include <FastqTemplateSource.h>

FastqIOAsFQTSource::FastqIOAsFQTSource(std::shared_ptr<FastqIO> handler)
    : in(std::move(handler))
{
    if(!in->isReader()) {
        throw std::invalid_argument("Attempt to use a non reader FastqIO object as a FastQTemplateSource");
    }
}


FastqIO::READ_RESULT FastqIOAsFQTSource::getTemplate(FastqTemplate_t & fqt) {
    FastqIO::READ_RESULT res = in->next_template(fqt);
    if(res != FastqIO::READ_PASS && res != FastqIO::READ_EOF ) {
        throw std::invalid_argument("Malformed fastq entry");
    }
    return res;
}

bool FastqIOAsFQTSource::operator()(FastqTemplate_t & fqt) {
    auto res = getTemplate(fqt); 
    return (res == FastqIO::READ_PASS);
}

std::vector<FastqTemplate_t> FastqIOAsFQTSource::get_block(size_t max_n) {
    FastqIO::READ_RESULT res = (in->isGood()) ? FastqIO::READ_PASS : 
                                                    FastqIO::READ_FAIL;
    std::vector<FastqTemplate_t> block;
    while(block.size() < max_n && res == FastqIO::READ_PASS) {
        FastqTemplate_t fqt;
        res = getTemplate(fqt);
        if(res == FastqIO::READ_PASS) {
            block.push_back(fqt);
        }
    }
    return block;
}



VectorAsFQTSource::VectorAsFQTSource(std::shared_ptr<std::vector<FastqTemplate_t>> vector) 
    : vec(vector), pos(0)
{
}

bool VectorAsFQTSource::operator()(FastqTemplate_t & fqt) {
    if(pos < vec->size()) {
        fqt = vec->at(pos++);
        return true;
    }
    return false;
}

std::vector<FastqTemplate_t> VectorAsFQTSource::get_block(size_t max_n){
    std::vector<FastqTemplate_t> res;
    res.reserve(max_n);
    size_t i = 0;
    while(i < max_n && pos < vec->size()){
        res.push_back(vec->at(pos++));
    }
    return res;
}
