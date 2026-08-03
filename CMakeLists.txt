cmake_minimum_required(VERSION 3.23)

project(SyncDedup)


add_library(Sketcher)
target_sources(Sketcher
    PRIVATE
        src/Sketcher.cpp 
        
    PUBLIC
        FILE_SET sketcherHeaders
        TYPE HEADERS
        BASE_DIRS include
        FILES
            include/Sketcher.h
)

add_executable(syncdedup)
target_sources(syncdedup
    PRIVATE
        src/syncdedup.cpp
)
target_link_libraries(syncdedup
    PRIVATE
        Sketcher
)

if(CMAKE_BUILD_TYPE EQUAL "Test")
    add_executable(test)
    target_sources(test
        PRIVATE
            test/test.cpp
    )
endif()
