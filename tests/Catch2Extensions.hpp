#ifndef CATCH2_EXT_HEADER_GAURD_
#define CATCH2_EXT_HEADER_GAURD_

#include <iostream>

#include <catch2/catch_test_macros.hpp>

#define CAT2(a,b) a##b
#define CAT(a,b) CAT2(a,b)

#define TODO_TEST_CASE(name, tags, BODY)                  \
    TODO_TEST_CASE_IMPL(name,tags,BODY, __COUNTER__)


#define TODO_TEST_CASE_IMPL(name, tags, BODY, n)          \
    namespace {                                           \
        auto CAT(UNIQUE_BODY_, n) = [] BODY;              \
    }                                                     \
    TEST_CASE(name, tags " [!mayfail]")                   \
    {                                                     \
        CAT(UNIQUE_BODY_, n)();                           \
    }                                                     \
    TEST_CASE(name " (expected failure)",                 \
              tags " [!shouldfail]")                      \
    {                                                     \
        CAT(UNIQUE_BODY_, n)();                           \
    }

#endif //CATCH2_EXT_HEADER_GAURD_
