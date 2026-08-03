#ifndef TODO_TEST_HEADER_GAURD_
#define TODO_TEST_HEADER_GAURD_

#include <iostream>

#include <catch2/catch_test_macros.hpp>
#include <catch2/internal/catch_test_failure_exception.hpp>

template<typename F>
void todo_test(F&& f)
{
    try
    {
        f();
        std::cerr << "CATS\n";
        // The test unexpectedly passed.
        FAIL("TODO test unexpectedly passed. Remove todo_test() wrapper.");
    }
    catch (Catch::TestFailureException const&)
    {
        // Expected failure.
        std::cerr << "DOGS\n";
        SUCCEED("Expected TODO failure.");
    }
}


#endif //TODO_TEST_HEADER_GAURD_
