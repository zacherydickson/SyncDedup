cmake_minimum_required(VERSION 3.25)

project(SyncDedup)
enable_testing()
option(BUILD_TESTING "Enable testing and build tests" OFF)

###SKETCHER LIBRARY
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


###LOCALSYNCMERMAP LIBRARY
add_library(LocalSyncmerMap)
target_sources(LocalSyncmerMap
    PRIVATE
        src/LocalSyncmerMap.cpp
    PUBLIC
        FILE_SET localSyncmerMapHeaders
        TYPE HEADERS
        BASE_DIRS src
        FILES
            src/LocalSyncmerMap.h
)
set_target_properties(LocalSyncmerMap PROPERTIES LINKER_LANGUAGE CXX)

###ZLIB AND ZSTR LIBS
find_package(ZLIB)
include(FetchContent)
FetchContent_Declare(ZStrGitRepo
  GIT_REPOSITORY    "https://github.com/mateidavid/zstr" # can also be a local filesystem path!
  GIT_TAG           "master"
)
FetchContent_MakeAvailable(ZStrGitRepo) # defines INTERFACE target 'zstr::zstr'

###FASTQIO LIBRARY
add_library(FastqIO)
target_sources(FastqIO
    PRIVATE
        src/FastqIO.cpp
    PUBLIC
        FILE_SET fastqIOHeaders
        TYPE HEADERS
        BASE_DIRS src
        FILES
            src/FastqIO.h
)
target_link_libraries(FastqIO
    PRIVATE
        ZLIB::ZLIB
        zstr::zstr
)

###CTPL LIBRARY
add_library(CTPL)
target_sources(CTPL
    PUBLIC
        FILE_SET ctplHeaders
        TYPE HEADERS
        BASE_DIRS CTPL
        FILES
            CTPL/ctpl.h
            CTPL/ctpl_stl.h
)
set_target_properties(CTPL PROPERTIES LINKER_LANGUAGE CXX)

##SKETCHHASHFACTORY LIBRARY
add_library(SketchHashFactory)
target_sources(SketchHashFactory
    PRIVATE
        src/SketchHashFactory.cpp
    PUBLIC
        FILE_SET fastqIOHeaders
        TYPE HEADERS
        BASE_DIRS src
        FILES
            src/SketchHashFactory.h
)
target_link_libraries(SketchHashFactory
    PRIVATE
        CTPL
        FastqIO
        LocalSyncmerMap
        Sketcher
)


###PROJECT EXECTUTABLE
add_executable(syncdedup)
target_sources(syncdedup
    PRIVATE
        src/syncdedup.cpp
)
#In case we need to require a language feature
#target_compile_features(syncdedup PRIVATE cxx_std_17)
target_link_libraries(syncdedup
    PRIVATE
        SketchHashFactory
)

###TESTING EXECUTABLES
if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

###INSTALL TARGETS
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set_property(CACHE CMAKE_INSTALL_PREFIX PROPERTY VALUE ${CMAKE_CURRENT_SOURCE_DIR})
endif()

install(
    TARGETS syncdedup Sketcher
    FILE_SET sketcherHeaders
)
