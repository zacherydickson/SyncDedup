#include <DuplicateFilter.h>

bool DuplicateFilter::operator()(   const ExtendedFastqTemplate_t & fqt,
                                    const LocalSyncmerMap & lsMap) const
{
    //TODO: Implement Me
    //Find candidate matches for the template
    //Construct synteny blocks from candidates
    //Filter on excessive indels
    //Filter on excessive subs
    return true;
}

