add_executable(TestSketcher)
target_sources(
    TestSketcher
    PRIVATE
        test.cpp
)

target_link_libraries(TestSketcher PRIVATE Sketcher)

function(SketcherTest op)
    add_test(
        NAME ${op}
        COMMAND TestSketcher ${op}
    )
endfunction()

SketcherTest("Syncmer Length")
