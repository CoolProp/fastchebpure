#pragma once

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>

#include "nlohmann/json.hpp"

// ---------------------------------------------------------------------------
// TreeHasher: FNV-1a 64 over a deterministic byte serialization of a parsed
// JSON tree.
//
// This is the fastchebpure side of the structural hash that CoolProp's
// superancillary-freshness test (src/Tests/CoolProp-Tests.cpp,
// "Superancillary eos_hash matches current EOS at bit level") and its Python
// mirror (dev/scripts/inject_superanc_check_points.py::eos_fnv1a_hex) use to
// verify that the fluid JSON's SUPERANCILLARY block was fit against an EOS
// that still matches the EOS shipped in CoolProp.
//
// All three implementations must agree byte-for-byte. See the block comment
// in CoolProp-Tests.cpp for the full byte-stream contract; the short version:
//
//   null      -> 'n'
//   false     -> 'f'
//   true      -> 't'
//   integer   -> 'i' then int64 two's-complement bits as LE u64
//   float     -> 'd' then IEEE-754 bits as LE u64
//   string    -> 's' then LE u64 UTF-8 byte count then UTF-8 bytes
//   array     -> 'a' then LE u64 length then each element walked
//   object    -> 'o' then LE u64 size then, for each key in sorted order:
//                LE u64 UTF-8 byte count, UTF-8 bytes, walked value
//
// Invariants:
//   * nlohmann::json is backed by std::map, so iterating an object yields
//     keys in sorted (lexicographic-by-UTF-8-bytes) order. Do not swap in
//     ordered_json / unordered_json — that would silently diverge from
//     Python's sorted(dict) walk.
//   * is_number_unsigned is lumped with is_number_integer and cast through
//     int64_t; every integer in CoolProp's fluid JSONs fits.
//   * All supported platforms are little-endian, so u64s are encoded LE
//     without an explicit swap. If that ever changes, swap all three impls
//     simultaneously.
//   * FNV-1a is not cryptographic — only used for change detection.
// ---------------------------------------------------------------------------
struct TreeHasher {
    uint64_t h = 0xcbf29ce484222325ULL;

    void mix_u8(uint8_t b) {
        h ^= b;
        h *= 0x100000001b3ULL;
    }
    void mix_bytes(const void* data, std::size_t n) {
        const auto* p = static_cast<const uint8_t*>(data);
        for (std::size_t i = 0; i < n; ++i) mix_u8(p[i]);
    }
    void mix_u64(uint64_t v) {
        for (int i = 0; i < 8; ++i) mix_u8(static_cast<uint8_t>((v >> (i * 8)) & 0xff));
    }

    void walk(const nlohmann::json& j) {
        if (j.is_null()) {
            mix_u8('n');
        } else if (j.is_boolean()) {
            mix_u8(j.get<bool>() ? 't' : 'f');
        } else if (j.is_number_integer() || j.is_number_unsigned()) {
            mix_u8('i');
            mix_u64(static_cast<uint64_t>(j.get<int64_t>()));
        } else if (j.is_number_float()) {
            mix_u8('d');
            uint64_t bits;
            double v = j.get<double>();
            std::memcpy(&bits, &v, 8);
            mix_u64(bits);
        } else if (j.is_string()) {
            const auto& s = j.get_ref<const std::string&>();
            mix_u8('s');
            mix_u64(s.size());
            mix_bytes(s.data(), s.size());
        } else if (j.is_array()) {
            mix_u8('a');
            mix_u64(j.size());
            for (const auto& el : j) walk(el);
        } else if (j.is_object()) {
            mix_u8('o');
            mix_u64(j.size());
            for (auto it = j.begin(); it != j.end(); ++it) {
                const auto& k = it.key();
                mix_u64(k.size());
                mix_bytes(k.data(), k.size());
                walk(it.value());
            }
        }
    }

    std::string hex() const {
        char buf[17];
        std::snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(h));
        return std::string(buf);
    }
};

/// Compute the 16-char lowercase hex FNV-1a 64 structural hash of an EOS
/// JSON subtree, with the SUPERANCILLARY key stripped if present.
inline std::string eos_fnv1a_hex(const nlohmann::json& eos) {
    nlohmann::json stripped = eos;
    if (stripped.is_object() && stripped.contains("SUPERANCILLARY")) {
        stripped.erase("SUPERANCILLARY");
    }
    TreeHasher th;
    th.walk(stripped);
    return th.hex();
}
