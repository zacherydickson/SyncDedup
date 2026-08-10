find_package(Catch2 3 REQUIRED)

add_library(Catch2Extensions)
target_sources(
    Catch2Extensions 
    PUBLIC
        FILE_SET HEADERS
        FILES
            Catch2Extensions.hpp
)
set_target_properties(Catch2Extensions PROPERTIES LINKER_LANGUAGE CXX)

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
        Catch2Extensions
)

add_executable(TestLocalSyncmerMap)
target_sources(
    TestLocalSyncmerMap
    PRIVATE
    testLocalSyncmerMap.cpp
)
target_compile_features(TestLocalSyncmerMap PRIVATE cxx_std_17)
target_link_libraries(TestLocalSyncmerMap
    PRIVATE 
        LocalSyncmerMap
        Catch2::Catch2WithMain
        Catch2Extensions
)

add_executable(TestFastqIO)
target_sources(TestFastqIO
    PRIVATE
        testFastqIO.cpp
)
target_compile_features(TestFastqIO PRIVATE cxx_std_17)
target_link_libraries(TestFastqIO
    PRIVATE 
        FastqIO
        Catch2::Catch2WithMain
        Catch2Extensions
)

include(CTest)
include(Catch)
catch_discover_tests(TestSketcher)
catch_discover_tests(TestLocalSyncmerMap)
catch_discover_tests(TestFastqIO)
