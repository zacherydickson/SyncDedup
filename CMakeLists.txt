cmake_minimum_required(VERSION 3.25)

project(SyncDedup)
enable_testing()
option(BUILD_TESTING "Enable testing and build tests" OFF)

add_library(Sketcher)
target_sources(Sketcher
    PRIVATE
        src/Sketcher.cpp 
        
    PUBLIC
        FILE_SET sketcherHeaders
        TYPE HEADERS
        BASE_DIRS src
        FILES
            src/Sketcher.h
)

add_library(LocalSyncmerMap)
target_sources(LocalSyncmerMap
    PRIVATE
        src/LocalSyncmerMap.cpp
    PUBLIC
    FILE_SET HEADERS
    BASE_DIRS src
    FILES
        src/LocalSyncmerMap.h
)
set_target_properties(LocalSyncmerMap PROPERTIES LINKER_LANGUAGE CXX)

find_package(ZLIB)
include(FetchContent)
FetchContent_Declare(ZStrGitRepo
  GIT_REPOSITORY    "https://github.com/mateidavid/zstr" # can also be a local filesystem path!
  GIT_TAG           "master"
)
FetchContent_MakeAvailable(ZStrGitRepo) # defines INTERFACE target 'zstr::zstr'

add_library(FastqIO)
target_sources(FastqIO
    PRIVATE
        src/FastqIO.cpp

    PUBLIC
        FILE_SET HEADERS
        BASE_DIRS src
        FILES
            src/FastqIO.h
)
target_link_libraries(FastqIO
    PRIVATE
        ZLIB::ZLIB
        zstr::zstr
)


add_executable(syncdedup)
target_sources(syncdedup
    PRIVATE
        src/syncdedup.cpp
)
#In case we need to require a language feature
#target_compile_features(syncdedup PRIVATE cxx_std_17)
target_link_libraries(syncdedup
    PRIVATE
        Sketcher
        LocalSyncmerMap
        FastqIO
)

if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set_property(CACHE CMAKE_INSTALL_PREFIX PROPERTY VALUE ${CMAKE_CURRENT_SOURCE_DIR})
endif()

install(
    TARGETS syncdedup Sketcher
    FILE_SET sketcherHeaders
)
