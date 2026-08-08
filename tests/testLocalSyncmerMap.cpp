#include <iostream>
#include <LocalSyncmerMap.h>
#include <string>
#include <catch2/catch_all.hpp>
#include <Catch2Extensions.hpp>

TEST_CASE("Local SyncmerMap object can be constructed", "[LocalSyncmerMap]") {
    REQUIRE_NOTHROW(LocalSyncmerMap());
}

SCENARIO("Map pass-through functions work", "[LocalSyncmerMap]") {
    GIVEN("An empty map, and some elements") {
        LocalSyncmerMap lsMap;
        uint64_t key1 = 0xCE41;
        uint64_t key2 = key1 + 1;
        Location singleHitLocation = {{0,0x01}};
        Location dualHitLocation = {{1,0x02},{2,0x03}};
        Location triHitLocation = {{3,0x04},{4,0x05},{5,0x06}};

        THEN("The map is empty") {
            REQUIRE( lsMap.size() == 0 );
            REQUIRE( lsMap.hits() == 0 );
            AND_WHEN("count() is called for a key not in the map") {
                size_t count = lsMap.count(key1);
                THEN("The count is zero") {
                    REQUIRE( count == 0 );
                }
            }
            AND_WHEN("In a constant context") {
                const LocalSyncmerMap & mapRef = lsMap;
                THEN("at() throws for an invalid key") {
                    REQUIRE_THROWS_AS ( mapRef.at(key1), std::out_of_range);
                }
                AND_THEN("The begin and end iterators are the same") {
                    REQUIRE( mapRef.begin() == mapRef.end() );
                    REQUIRE( mapRef.cbegin() == mapRef.cend() );
                }
                AND_WHEN("find() is called for an invalid key") {
                    auto it = mapRef.find(key1);
                    THEN("The result is the end iterator") {
                        REQUIRE( it == mapRef.end() );
                    }
                }
            }
        }

        WHEN("insert is called with a key and single hit Element pair"){
            lsMap.insert(std::make_pair(key1,singleHitLocation));
            THEN("The map's size increases") {
                REQUIRE( lsMap.size() == 1 );
                REQUIRE( lsMap.hits() == 1 );
                AND_WHEN("count() is called for the key") {
                    size_t count = lsMap.count(key1);
                    THEN("The count is 1") {
                        REQUIRE (count == 1);
                    }
                }
                AND_WHEN("In a constant context") {
                    const LocalSyncmerMap & mapRef = lsMap;
                    THEN("at() does not throw for the key") {
                        REQUIRE_NOTHROW( mapRef.at(key1) );
                    }
                    AND_THEN("The begin and end iterators are not the same") {
                        REQUIRE( mapRef.begin() != mapRef.end() );
                        REQUIRE( mapRef.cbegin() != mapRef.cend() );
                    }
                    AND_WHEN("find() is called for the key") {
                        auto it = mapRef.find(key1);
                        THEN("The result is the begin iterator") {
                            REQUIRE( it == mapRef.begin() );
                        }
                    }
                }
            }
            AND_WHEN("erase is called with the same key"){
                lsMap.erase(key1);
                THEN("The map is now empty again") {
                    REQUIRE( lsMap.size() == 0 );
                    REQUIRE( lsMap.hits() == 0 );
                }
            }
        }

        WHEN("insert() is called with move semantics for a two hit Location") {
            lsMap.insert(std::move(std::make_pair(key1,dualHitLocation)));
            THEN("The map's size and hits increase") {
                REQUIRE( lsMap.size() == 1 );
                REQUIRE( lsMap.hits() == 2 );
            }
            AND_WHEN("insert() is called with the same key for a single hit Location")
            {
                lsMap.insert(std::make_pair(key1,singleHitLocation));
                THEN("The map's hits increase, but the size does not") {
                    REQUIRE( lsMap.size() == 1 );
                    REQUIRE( lsMap.hits() == 3 );
                }
                AND_WHEN("insert() is called with a different key for a tri hit Location")
                {
                    lsMap.insert(std::make_pair(key2,triHitLocation));
                    THEN("The map's size and hits increase"){
                        REQUIRE( lsMap.size() == 2 );
                        REQUIRE( lsMap.hits() == 6 );
                    }
                    AND_WHEN("clear() is called") {
                        lsMap.clear();
                        THEN("The map is now empty") {
                            REQUIRE( lsMap.size() == 0 );
                            REQUIRE( lsMap.hits() == 0 );
                        }
                    }
                }
            }
        }
    }
}

SCENARIO ("Adding and removing Sketch Elements", "[LocalSyncmerMap]") {
    GIVEN("An empty map and two labeled sketch elements") {
        LocalSyncmerMap lsMap;
        uint64_t id1 = 0x01;
        uint64_t id2 = 0x02;
        SketchElement se1 = {0xCE41, 1};
        SketchElement se2 = {0xCE41, 2};

        THEN("No sketch elements can be found") {
            REQUIRE_THROWS_AS( lsMap.at(se1) , std::out_of_range );
            REQUIRE_THROWS_AS( lsMap.at(se2) , std::out_of_range );
            AND_WHEN("count() is called") {
                size_t count = lsMap.count(se1);
                THEN("The count is zero") {
                    REQUIRE( count == 0 );
                }
            }
        }

        WHEN("insert() is called with a the sketch element") {
            lsMap.insert(id1, se1);
            THEN("The map's size and hits increases") {
                REQUIRE ( lsMap.size() == 1 );
                REQUIRE ( lsMap.hits() == 1 );
                    
            }
            AND_WHEN("find() is called with the element's hash") {
                auto it = lsMap.find(se1.hash);
                THEN("The elements's hash is found") {
                    REQUIRE ( it != lsMap.end() );
                    REQUIRE ( it->first == se1.hash );
                }
            }
            AND_WHEN("count() is called") {
                size_t count = lsMap.count(se1);
                THEN("The count is 1") {
                    REQUIRE( count == 1 );
                }
            }
            AND_WHEN("another element with the same hash is added") {
                lsMap.insert(id2, se2);
                THEN("The hits increase, but the size does not") {
                    REQUIRE (lsMap.size() == 1 );
                    REQUIRE (lsMap.hits() == 2 );
                }
                AND_WHEN("The first element is erased") {
                    lsMap.erase(id1, se1);
                    THEN("Only the first element is removed") {
                        REQUIRE( lsMap.size() == 1 ); 
                        REQUIRE( lsMap.hits() == 1 ); 
                        REQUIRE_NOTHROW( lsMap.at(se2) );
                        AND_WHEN("at is called with the second elements hash") {
                            REQUIRE( lsMap.begin()->second.front().id == id2 );
                        }
                    }
                }
            }
        }
    }
}

SCENARIO ("Adding an entire Sketch to a LocalSyncmerMap","[LocalSyncmerMap]") {
    GIVEN("An empty map and a Sketch") {
        LocalSyncmerMap lsMap;
        size_t id = 0x1;
        Sketch sketch = { {0xCE41,1}, {0xCE42,2}, {0xCE43,3}};

        WHEN ("insert() is called for the sketch") {
            lsMap.insert(id, sketch);
            THEN("All elements are added") {
                REQUIRE( lsMap.size() == 3 );
                REQUIRE( lsMap.hits() == 3 );
            }
        }
    }
}
