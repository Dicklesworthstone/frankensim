//! `philox_normals` admission, draw-index semantics, refusals, and cross-check vectors.
//!
//! Integer fields are bitwise. Box–Muller normals are native `fs_math::det` bits;
//! they are not claimed bitwise-equal to TypeScript host `Math`.

use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use fs_rand::philox::philox4x32_10;
use fs_rand::philox_normals::{
    admit_philox_normals, philox_normals, philox_normals_admitted, TargetKind,
    PHILOX_NORMALS_MAX_COUNT,
};
use fs_rand::{
    STREAM_CHECKPOINT_VERSION, STREAM_SEMANTICS_VERSION, Stream, StreamCheckpoint, StreamKey,
};

const SUITE: &str = "fs-rand/philox-normals";
const BEAD: &str = "am-fs-export-philox-normals-xnv";
const LAST_ACCEPTED_START: u64 = 18_446_744_073_709_551_613;
const FIRST_REFUSED_START: u64 = 18_446_744_073_709_551_614;
const VECTOR_REL: &str = "vectors/philox_cross_check.json";
const SHA256_PLACEHOLDER: &str = "SHA256_PLACEHOLDER_64_HEX_CHARS_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx";

const RANDOM123: [([u32; 4], [u32; 2], [u32; 4]); 3] = [
    (
        [0, 0, 0, 0],
        [0, 0],
        [0x6627_e8d5, 0xe169_c58d, 0xbc57_ac4c, 0x9b00_dbd8],
    ),
    (
        [u32::MAX; 4],
        [u32::MAX; 2],
        [0x408f_276d, 0x41c8_3b0e, 0xa20b_c7c6, 0x6d54_51fd],
    ),
    (
        [0x243f_6a88, 0x85a3_08d3, 0x1319_8a2e, 0x0370_7344],
        [0xa409_3822, 0x299f_31d0],
        [0xd16c_fe09, 0x94fd_cceb, 0x5001_e420, 0x2412_6ea1],
    ),
];

const SEEDS: [u64; 7] = [
    0,
    1,
    4_294_967_295,
    4_294_967_296,
    9_007_199_254_740_992,
    9_007_199_254_740_993,
    u64::MAX,
];
const KERNELS: [u32; 3] = [0, 1, u32::MAX];
const TILES: [u32; 3] = [0, 1, u32::MAX];
const INDICES: [u64; 5] = [0, 1, 4_294_967_295, 4_294_967_296, u64::MAX - 1];
const BELOW_N: [u64; 4] = [1, 6, 4_294_967_297, 18_446_744_073_709_551_557];

fn verdict(case: &str, outcome: &str, extra: &str) {
    println!(
        "{{\"suite\":\"{SUITE}\",\"beadId\":\"{BEAD}\",\"case\":\"{case}\",\"verdict\":\"{outcome}\"{extra}}}"
    );
}

fn key(seed: u64, kernel: u32, tile: u32) -> StreamKey {
    StreamKey { seed, kernel, tile }
}

fn stream_at(key: StreamKey, index: u64) -> Stream {
    Stream::resume(StreamCheckpoint::current(key, index)).expect("current checkpoint")
}

#[test]
fn known_answer_vectors() {
    let t0 = Instant::now();
    for (i, (ctr, k, want)) in RANDOM123.iter().enumerate() {
        let got = philox4x32_10(*ctr, *k);
        assert_eq!(got, *want, "Random123 KAT {i}");
    }
    verdict(
        "known_answer_vectors",
        "pass",
        &format!(
            ",\"detail\":\"3 Random123 vectors\",\"durationMs\":{}",
            t0.elapsed().as_millis()
        ),
    );
}

#[test]
fn philox_normals_count_and_finiteness() {
    let t0 = Instant::now();
    let values = philox_normals(1, 0, 0, 0, 8).expect("admitted");
    assert_eq!(values.len(), 8);
    assert!(values.iter().all(|z| z.is_finite()));
    verdict(
        "philox_normals_count_and_finiteness",
        "pass",
        &format!(",\"durationMs\":{}", t0.elapsed().as_millis()),
    );
}

#[test]
fn resume_positioning() {
    let t0 = Instant::now();
    let seed = 99u64;
    let kernel = 17u32;
    let tile = 29u32;
    let start = 4u64;
    let count = 5usize;
    let from_export = philox_normals(seed, kernel, tile, start, count).unwrap();
    let mut s = stream_at(key(seed, kernel, tile), start);
    let mut manual = Vec::new();
    for _ in 0..count {
        manual.push(s.next_normal());
    }
    assert_eq!(from_export, manual);
    verdict(
        "resume_positioning",
        "pass",
        &format!(",\"durationMs\":{}", t0.elapsed().as_millis()),
    );
}

#[test]
fn draw_indexing_suffix() {
    let t0 = Instant::now();
    let from_zero = philox_normals(1, 0, 0, 0, 8).unwrap();
    let from_six = philox_normals(1, 0, 0, 6, 8).unwrap();
    assert_eq!(&from_six[..5], &from_zero[3..]);
    verdict(
        "draw_indexing_suffix",
        "pass",
        &format!(",\"durationMs\":{}", t0.elapsed().as_millis()),
    );
}

#[test]
fn moments_fixed_seed() {
    let t0 = Instant::now();
    let n = 50_000usize;
    let values = philox_normals(12_345, 1, 2, 0, n).unwrap();
    let mean = values.iter().sum::<f64>() / n as f64;
    let var = values.iter().map(|z| (z - mean) * (z - mean)).sum::<f64>() / n as f64;
    // Precomputed bounds: se(mean) ≈ 1/sqrt(50000) ≈ 0.0045; |mean|<0.03 and
    // var in [0.94, 1.06] are >6 se for a unit normal. Not retried.
    assert!(mean.abs() < 0.03, "mean {mean}");
    assert!((0.94..=1.06).contains(&var), "var {var}");
    verdict(
        "moments_fixed_seed",
        "pass",
        &format!(
            ",\"mean\":{mean},\"variance\":{var},\"n\":{n},\"durationMs\":{}",
            t0.elapsed().as_millis()
        ),
    );
}

#[test]
fn counter_boundary_refusal() {
    let t0 = Instant::now();
    let ok = philox_normals(0, 0, 0, LAST_ACCEPTED_START, 1).expect("last pair admitted");
    assert_eq!(ok.len(), 1);
    assert!(ok[0].is_finite());
    let err = philox_normals(0, 0, 0, FIRST_REFUSED_START, 1).expect_err("would wrap");
    assert_eq!(err.code, "stream-index-overflow");
    assert_eq!(err.target_kind, TargetKind::RefusalCode);
    assert!(err.details.contains("18446744073709551614"));
    verdict(
        "counter_boundary_refusal",
        "pass",
        &format!(",\"durationMs\":{}", t0.elapsed().as_millis()),
    );
}

#[test]
fn refusals() {
    let t0 = Instant::now();
    let zero = admit_philox_normals(0, 0, 0, 0, 0).expect_err("count 0");
    assert_eq!(zero.code, "invalid-parameter");
    assert_eq!(zero.target_kind, TargetKind::RefusalCode);
    assert!(zero.details.contains("\"count\""));

    let budget = admit_philox_normals(0, 0, 0, 0, PHILOX_NORMALS_MAX_COUNT + 1).expect_err("budget");
    assert_eq!(budget.code, "budget-exhausted");
    assert_eq!(budget.target_kind, TargetKind::ExecutionOutcome);

    let empty_is_not_ok = philox_normals(0, 0, 0, 0, 0);
    assert!(empty_is_not_ok.is_err(), "empty Vec is not a refusal stand-in");
    verdict(
        "refusals",
        "pass",
        &format!(",\"durationMs\":{}", t0.elapsed().as_millis()),
    );
}

#[test]
fn admitted_rechecks_forged_counts_via_public_admit() {
    let spec = admit_philox_normals(7, 3, 5, 2, 4).unwrap();
    let a = philox_normals_admitted(&spec).unwrap();
    let b = philox_normals(7, 3, 5, 2, 4).unwrap();
    assert_eq!(a, b);
}

fn hex_u32s(words: [u32; 4]) -> String {
    format!(
        "[\"{:08x}\",\"{:08x}\",\"{:08x}\",\"{:08x}\"]",
        words[0], words[1], words[2], words[3]
    )
}

fn hex_u32_pair(words: [u32; 2]) -> String {
    format!("[\"{:08x}\",\"{:08x}\"]", words[0], words[1])
}

fn hex_f64(z: f64) -> String {
    format!("{:016x}", z.to_bits())
}

fn hex_bytes(bytes: &[u8]) -> String {
    let mut s = String::with_capacity(bytes.len() * 2);
    for b in bytes {
        s.push_str(&format!("{b:02x}"));
    }
    s
}

fn below_record(k: StreamKey, index: u64, n: u64) -> (u64, u64) {
    let mut s = stream_at(k, index);
    let v = s.next_below(n);
    (v, s.index().wrapping_sub(index))
}

fn render_vectors(revision: &str) -> String {
    let mut out = String::with_capacity(1 << 18);
    out.push_str("{\n");
    out.push_str("  \"provenance\": {\n");
    out.push_str(&format!("    \"frankensimRevision\": \"{revision}\",\n"));
    out.push_str(&format!(
        "    \"streamSemanticsVersion\": {STREAM_SEMANTICS_VERSION},\n"
    ));
    out.push_str(&format!(
        "    \"streamCheckpointVersion\": {STREAM_CHECKPOINT_VERSION},\n"
    ));
    out.push_str("    \"generationCommand\": \"cargo test -p fs-rand --test philox_normals emit_cross_check_vectors\",\n");
    out.push_str(&format!("    \"sha256\": \"{SHA256_PLACEHOLDER}\",\n"));
    out.push_str("    \"parity\": {\n");
    out.push_str("      \"block\": \"bitwise\",\n");
    out.push_str("      \"u64\": \"bitwise\",\n");
    out.push_str("      \"f64\": \"bitwise\",\n");
    out.push_str("      \"normalNativeWasm\": \"bitwise when fs_math::det holds on wasm32\",\n");
    out.push_str("      \"normalVsTypeScript\": \"tolerance, not bitwise; TypeScript uses host Math (philox-box-muller-host-v1); FrankenSim uses det ln/cos/sqrt\",\n");
    out.push_str("      \"below\": \"bitwise including rejection draws\",\n");
    out.push_str("      \"checkpoint\": \"bitwise 83-byte StreamCheckpoint frame\"\n");
    out.push_str("    }\n");
    out.push_str("  },\n");
    out.push_str("  \"knownAnswers\": [\n");
    for (i, (ctr, k, want)) in RANDOM123.iter().enumerate() {
        if i > 0 {
            out.push_str(",\n");
        }
        out.push_str("    {\n");
        out.push_str(&format!("      \"counter\": {},\n", hex_u32s(*ctr)));
        out.push_str(&format!("      \"key\": {},\n", hex_u32_pair(*k)));
        out.push_str(&format!("      \"block\": {},\n", hex_u32s(*want)));
        out.push_str("      \"source\": \"Random123 kat_vectors philox4x32 10 rounds\"\n");
        out.push_str("    }");
    }
    out.push_str("\n  ],\n");
    out.push_str("  \"positions\": [\n");
    let mut first_pos = true;
    for seed in SEEDS {
        for kernel in KERNELS {
            for tile in TILES {
                for index in INDICES {
                    if !first_pos {
                        out.push_str(",\n");
                    }
                    first_pos = false;
                    let k = key(seed, kernel, tile);
                    let block = Stream::at(k, index);
                    let mut s = stream_at(k, index);
                    let u = s.next_u64();
                    let mut s = stream_at(k, index);
                    let f = s.next_f64();
                    let mut s = stream_at(k, index);
                    let mut d0 = stream_at(k, index);
                    let n0 = d0.next_u64();
                    let n1 = d0.next_u64();
                    let zn = s.next_normal();
                    let ckpt = StreamCheckpoint::current(k, index).to_canonical_le_bytes();
                    out.push_str("    {\n");
                    out.push_str(&format!("      \"seed\": \"{seed}\",\n"));
                    out.push_str(&format!("      \"streamKernel\": {kernel},\n"));
                    out.push_str(&format!("      \"tile\": {tile},\n"));
                    out.push_str(&format!("      \"index\": \"{index}\",\n"));
                    out.push_str(&format!("      \"block\": {},\n", hex_u32s(block)));
                    out.push_str(&format!("      \"u64\": \"{u}\",\n"));
                    out.push_str(&format!("      \"f64Bits\": \"{}\",\n", hex_f64(f)));
                    out.push_str(&format!("      \"normalBits\": \"{}\",\n", hex_f64(zn)));
                    out.push_str(&format!("      \"normalDraws\": [\"{n0}\", \"{n1}\"],\n"));
                    out.push_str("      \"below\": [\n");
                    for (bi, n) in BELOW_N.iter().enumerate() {
                        let (val, draws) = below_record(k, index, *n);
                        if bi > 0 {
                            out.push_str(",\n");
                        }
                        out.push_str(&format!(
                            "        {{\"n\":\"{n}\",\"value\":\"{val}\",\"drawsConsumed\":{draws}}}"
                        ));
                    }
                    out.push_str("\n      ],\n");
                    out.push_str(&format!("      \"checkpoint\": \"{}\"\n", hex_bytes(&ckpt)));
                    out.push_str("    }");
                }
            }
        }
    }
    out.push_str("\n  ],\n");
    out.push_str("  \"normalSequences\": [\n");
    for (i, start) in [0u64, 6].iter().enumerate() {
        let bits: Vec<String> = philox_normals(1, 0, 0, *start, 8)
            .unwrap()
            .into_iter()
            .map(hex_f64)
            .map(|h| format!("\"{h}\""))
            .collect();
        if i > 0 {
            out.push_str(",\n");
        }
        out.push_str("    {\n");
        out.push_str("      \"seed\": \"1\",\n");
        out.push_str("      \"streamKernel\": 0,\n");
        out.push_str("      \"tile\": 0,\n");
        out.push_str(&format!("      \"startIndex\": \"{start}\",\n"));
        out.push_str("      \"count\": 8,\n");
        out.push_str(&format!("      \"normalBits\": [{}]\n", bits.join(", ")));
        out.push_str("    }");
    }
    out.push_str("\n  ]\n}\n");
    let digest = sha256_hex(out.as_bytes());
    out.replace(SHA256_PLACEHOLDER, &digest)
}

fn vector_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests").join(VECTOR_REL)
}

fn frankensim_revision() -> String {
    if let Ok(text) = fs::read_to_string(vector_path()) {
        if let Some(start) = text.find("\"frankensimRevision\": \"") {
            let rest = &text[start + 23..];
            if let Some(end) = rest.find('"') {
                return rest[..end].to_string();
            }
        }
    }
    std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .current_dir(concat!(env!("CARGO_MANIFEST_DIR"), "/../.."))
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .unwrap_or_else(|| "unknown".to_string())
}

#[test]
fn emit_cross_check_vectors() {
    let t0 = Instant::now();
    let revision = frankensim_revision();
    let rendered = render_vectors(&revision);
    assert_eq!(
        rendered.matches("\"index\":").count(),
        315,
        "expected 315 position records"
    );
    let path = vector_path();
    if std::env::var_os("WRITE_PHILOX_VECTORS").is_some() {
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(&path, &rendered).unwrap();
    }
    let committed = fs::read_to_string(&path).unwrap_or_default();
    assert_eq!(
        committed, rendered,
        "vector file drifted; WRITE_PHILOX_VECTORS=1 cargo test -p fs-rand --test philox_normals emit_cross_check_vectors"
    );
    verdict(
        "emit_cross_check_vectors",
        "pass",
        &format!(
            ",\"positions\":315,\"durationMs\":{}",
            t0.elapsed().as_millis()
        ),
    );
}

#[test]
fn golden_hash() {
    let t0 = Instant::now();
    let text = fs::read_to_string(vector_path()).expect("committed vector file");
    let marker = "\"sha256\": \"";
    let start = text.find(marker).expect("provenance sha256") + marker.len();
    let recorded = &text[start..start + 64];
    let for_hash = text.replacen(recorded, SHA256_PLACEHOLDER, 1);
    let digest = sha256_hex(for_hash.as_bytes());
    assert_eq!(recorded, digest, "provenance sha256 is the hash of the file with the placeholder");
    verdict(
        "golden_hash",
        "pass",
        &format!(
            ",\"sha256\":\"{digest}\",\"durationMs\":{}",
            t0.elapsed().as_millis()
        ),
    );
}

fn sha256_hex(input: &[u8]) -> String {
    let hash = sha256(input);
    hex_bytes(&hash)
}

fn sha256(mut msg: &[u8]) -> [u8; 32] {
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];
    let bit_len = (msg.len() as u64) * 8;
    let mut padded = msg.to_vec();
    padded.push(0x80);
    while (padded.len() % 64) != 56 {
        padded.push(0);
    }
    padded.extend_from_slice(&bit_len.to_be_bytes());
    let mut w = [0u32; 64];
    for chunk in padded.chunks(64) {
        for i in 0..16 {
            w[i] = u32::from_be_bytes(chunk[i * 4..i * 4 + 4].try_into().unwrap());
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let mut a = h[0];
        let mut b = h[1];
        let mut c = h[2];
        let mut d = h[3];
        let mut e = h[4];
        let mut f = h[5];
        let mut g = h[6];
        let mut hh = h[7];
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let t1 = hh
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let t2 = s0.wrapping_add(maj);
            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(t1);
            d = c;
            c = b;
            b = a;
            a = t1.wrapping_add(t2);
        }
        h[0] = h[0].wrapping_add(a);
        h[1] = h[1].wrapping_add(b);
        h[2] = h[2].wrapping_add(c);
        h[3] = h[3].wrapping_add(d);
        h[4] = h[4].wrapping_add(e);
        h[5] = h[5].wrapping_add(f);
        h[6] = h[6].wrapping_add(g);
        h[7] = h[7].wrapping_add(hh);
    }
    let mut out = [0u8; 32];
    for (i, word) in h.iter().enumerate() {
        out[i * 4..i * 4 + 4].copy_from_slice(&word.to_be_bytes());
    }
    out
}

#[test]
fn sha256_self_check() {
    assert_eq!(
        sha256_hex(b""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
    assert_eq!(
        sha256_hex(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
}
