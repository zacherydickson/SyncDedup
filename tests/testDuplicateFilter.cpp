#include <catch2/catch_all.hpp>
#include <Catch2Extensions.hpp>
#include <DuplicateFilter.h>

using namespace DuplicateFilter;

//FastqTemplate_t initFqt1pair = { "Read1",
//    {   { "GACTTTTCCAG", "1:N:0:ACTG", "ABCDEFGHIJK" },
//        { "CTGGAAAAGTC", "2:N:0:ACTG", "KJIHGFEDCBA" } } };
//FastqTemplate_t initFqt1fwd = { initFqt1pair.name,
//    { initFqt1pair.segVec[0] } };
//FastqTemplate_t initFqt1rev = { initFqt1pair.name,
//    { initFqt1pair.segVec[1] } };
//double initMeanQual1fwd = 37;
//double initMeanQual1rev = 37;
//double initMeanQual1pair = 37;
//
//FastqTemplate_t initFqt2pair = { "Read2",
//    {   { "CCCCCCCCC", "1:N:0:ACCC", "IIIIIIIII" },
//        { "GGGGGGGGGG", "2:N:0:ACCC", "[[[[[[[[[[" } } };
//FastqTemplate_t initFqt2fwd = { initFqt2pair.name,
//    { initFqt2pair.segVec[0] } };
//FastqTemplate_t initFqt2rev = { initFqt2pair.name,
//    { initFqt2pair.segVec[1] } };
//double initMeanQual2fwd = 40;
//double initMeanQual2rev = 58;
//double initMeanQual2pair =  (   initMeanQual2fwd *
//                                    initFqt2fwd.segVec[0].seq.length() +
//                                initMeanQual2rev *
//                                    initFqt2rev.segVec[0].seq.length() ) / 
//                            (
//                                initFqt2fwd.segVec[0].seq.length() +
//                                initFqt2rev.segVec[0].seq.length()
//                            );
//
//std::unique_ptr<std::iostream > constructSS(
//        std::vector<FastqTemplate_t> initVec = {}, bool postappend=false,
//        FastqIO::READ_RESULT res = FastqIO::READ_PASS)
//{
//    std::unique_ptr<std::iostream> ss_ptr(new std::stringstream());
//    for(const FastqTemplate_t & fqt : initVec ){
//        for(const FastqSegment_t & seg : fqt.segVec) {
//            if(res != FastqIO::READ_MISSING_LEADER1) {
//                *ss_ptr <<  '@';
//            }
//            *ss_ptr << fqt.name << ' ' << seg.desc << "\n";
//            if(res == FastqIO::READ_EOF) { goto endloop; }
//            *ss_ptr <<seg.seq << "\n";
//            if(res != FastqIO::READ_MISSING_LEADER2) {
//                *ss_ptr << '+';
//            }
//            *ss_ptr << "\n" << seg.qual << "\n";
//        }
//    }
//    //std::cerr << dynamic_cast<std::stringstream*>(ss_ptr.get())->str() << "\n";
//endloop:
//    if(!postappend){
//        ss_ptr->seekg(0);
//        ss_ptr->seekp(0);
//    }
//    return ss_ptr;
//}

SCENARIO ("Template Summaries correctly compare less than","[TemplateSummary_t]") {
    TemplateSummary_t base{30.0,200,123};
    double qualAdj = GENERATE(-1.0, 0.0, 1.0);
    INFO("and Given: qualAdj is " << qualAdj); 
    int lenAdj = GENERATE(-1, 0, 1);
    INFO("and Given: lenAdj is " << lenAdj); 
    int idxAdj = GENERATE(-1, 0, 1);
    INFO("and Given: idxAdj is " << idxAdj); 
    bool exp = (idxAdj == -1);
    if(lenAdj == 1) {exp = true;} else if(lenAdj == -1) {exp = false;}
    if(qualAdj < -TemplateSummary_t::QualTolerance) {
        exp = false;
    } else if (qualAdj > TemplateSummary_t::QualTolerance) {
        exp = true;
    }
    GIVEN("An adjusted template") {
        TemplateSummary_t adj{base.meanQual+qualAdj,base.length+lenAdj,base.idx+idxAdj};
        THEN("Then base template is less than the adjusted when appropriate") {
            REQUIRE( (base < adj) == exp );
        }
    }
}
