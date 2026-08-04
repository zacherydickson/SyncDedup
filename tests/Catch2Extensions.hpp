#ifndef CATCH2_EXT_HEADER_GAURD_
#define CATCH2_EXT_HEADER_GAURD_

#include <iostream>

#include <catch2/catch_test_macros.hpp>

#define TODO_TEST_CASE(name, tags, BODY)                  \
    namespace {                                           \
        auto UNIQUE_BODY = [] BODY;                       \
    }                                                     \
    TEST_CASE(name, tags " [!mayfail]")                   \
    {                                                     \
        UNIQUE_BODY();                                    \
    }                                                     \
    TEST_CASE(name " (expected failure)",                 \
              tags " [!shouldfail]")                      \
    {                                                     \
        UNIQUE_BODY();                                    \
    }

#endif //CATCH2_EXT_HEADER_GAURD_
