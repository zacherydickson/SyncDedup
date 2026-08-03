find_package(Catch2 3 REQUIRED)

add_library(todoTest)
target_sources(
    todoTest
    PUBLIC
        FILE_SET HEADERS
        FILES
            todotest.h
)
set_target_properties(todoTest PROPERTIES LINKER_LANGUAGE CXX)

add_executable(TestSketcher)
target_sources(
    TestSketcher
    PRIVATE
    testSketcher.cpp
)
target_compile_features(TestSketcher PRIVATE cxx_std_17)
target_link_libraries(TestSketcher
    PRIVATE 
        Sketcher
        Catch2::Catch2WithMain
        todoTest
)

include(CTest)
include(Catch)
catch_discover_tests(TestSketcher)
