// Unit tests for the structural FNV-1a 64 EOS hash. The algorithm must stay
// byte-for-byte compatible with the Python and CoolProp C++ implementations
// (dev/scripts/inject_superanc_check_points.py::eos_fnv1a_hex and
// src/Tests/CoolProp-Tests.cpp::TreeHasher) so that CoolProp's
// superancillary-freshness check can compare hashes verbatim.

#include <catch2/catch_test_macros.hpp>

#include "nlohmann/json.hpp"

#include "eos_hash.hpp"

TEST_CASE("FNV-1a 64 structural hash, inline literal", "[eos_hash]") {
    // Expected hash produced by the Python reference implementation
    // (eos_fnv1a_hex in dev/scripts/inject_superanc_check_points.py) over
    // this exact JSON value. This literal covers every branch of the walk:
    // integer, float, bool, null, string, nested array, nested object, and
    // the SUPERANCILLARY key that must be stripped before hashing.
    auto j = nlohmann::json::parse(R"({
        "a": 1,
        "b": 0.5,
        "c": true,
        "d": null,
        "e": "hello",
        "f": [1, 2.5, false, "x"],
        "g": {"nested": [0]},
        "SUPERANCILLARY": {"should_be_stripped": "yes"}
    })");
    REQUIRE(eos_fnv1a_hex(j) == "c3bf8bcc483d1bee");
}

TEST_CASE("SUPERANCILLARY key is stripped before hashing", "[eos_hash]") {
    auto with_super = nlohmann::json::parse(R"({
        "x": 1,
        "SUPERANCILLARY": {"junk": [1,2,3]}
    })");
    auto without_super = nlohmann::json::parse(R"({"x": 1})");
    REQUIRE(eos_fnv1a_hex(with_super) == eos_fnv1a_hex(without_super));
}

TEST_CASE("Integer, float, bool, string tags are distinct", "[eos_hash]") {
    // These four values (0 int, 0.0 float, false, "") would all collapse to
    // the same bytes under a naive encoding; the type-tag prefix keeps them
    // distinct.
    auto h_int    = eos_fnv1a_hex(nlohmann::json::parse(R"({"v": 0})"));
    auto h_float  = eos_fnv1a_hex(nlohmann::json::parse(R"({"v": 0.0})"));
    auto h_false  = eos_fnv1a_hex(nlohmann::json::parse(R"({"v": false})"));
    auto h_empty  = eos_fnv1a_hex(nlohmann::json::parse(R"({"v": ""})"));
    REQUIRE(h_int != h_float);
    REQUIRE(h_int != h_false);
    REQUIRE(h_int != h_empty);
    REQUIRE(h_float != h_false);
    REQUIRE(h_float != h_empty);
    REQUIRE(h_false != h_empty);
}
