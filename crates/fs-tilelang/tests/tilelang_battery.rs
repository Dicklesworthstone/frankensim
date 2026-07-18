//! fs-tilelang battery (wf9.11): the three reference kernels of the
//! acceptance criteria (batched axpy, stencil apply, SDF-style
//! trilinear grid sample) plus a deterministic-sum reduction kernel —
//! each written ONCE, lowered to scalar + lane variants, checked
//! against hand-written oracles, tier-equivalent bitwise (G0),
//! deterministic on repeat (G5), with intensity metadata logged for
//! the roofline harness. The macro's auto-generated twin tests run
//! alongside (visible as `__twin_tests` in the test list).

use fs_rand::StreamKey;
use fs_tilelang::{
    DeterminismClass, KernelMeta, MAX_KERNEL_NAME_BYTES, MAX_LOG_LABEL_BYTES, MAX_LOG_RECORD_BYTES,
    MAX_METADATA_JSON_BYTES, MetadataRenderError, ReductionKind, kernel,
};
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Eq)]
enum JsonValue {
    Object(BTreeMap<String, JsonValue>),
    String(String),
    Number(String),
}

struct StrictJsonParser<'a> {
    source: &'a str,
    position: usize,
}

impl<'a> StrictJsonParser<'a> {
    fn new(source: &'a str) -> Self {
        Self {
            source,
            position: 0,
        }
    }

    fn parse(mut self) -> Result<JsonValue, String> {
        self.skip_whitespace();
        let value = self.parse_value()?;
        self.skip_whitespace();
        if self.position != self.source.len() {
            let position = self.position;
            return Err(format!("trailing JSON bytes at {position}"));
        }
        Ok(value)
    }

    fn bytes(&self) -> &[u8] {
        self.source.as_bytes()
    }

    fn peek(&self) -> Option<u8> {
        self.bytes().get(self.position).copied()
    }

    fn take(&mut self, expected: u8) -> bool {
        if self.peek() == Some(expected) {
            self.position += 1;
            true
        } else {
            false
        }
    }

    fn expect(&mut self, expected: u8) -> Result<(), String> {
        if self.take(expected) {
            Ok(())
        } else {
            let expected = char::from(expected);
            let position = self.position;
            let found = self.peek().map(char::from);
            Err(format!(
                "expected byte {expected:?} at {position}, found {found:?}"
            ))
        }
    }

    fn skip_whitespace(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\n' | b'\r' | b'\t')) {
            self.position += 1;
        }
    }

    fn parse_value(&mut self) -> Result<JsonValue, String> {
        self.skip_whitespace();
        match self.peek() {
            Some(b'{') => self.parse_object(),
            Some(b'"') => self.parse_string().map(JsonValue::String),
            Some(b'-' | b'0'..=b'9') => self.parse_number().map(JsonValue::Number),
            other => {
                let position = self.position;
                Err(format!("unsupported JSON value at {position}: {other:?}"))
            }
        }
    }

    fn parse_object(&mut self) -> Result<JsonValue, String> {
        self.expect(b'{')?;
        self.skip_whitespace();
        let mut fields = BTreeMap::new();
        if self.take(b'}') {
            return Ok(JsonValue::Object(fields));
        }
        loop {
            self.skip_whitespace();
            let key = self.parse_string()?;
            self.skip_whitespace();
            self.expect(b':')?;
            let value = self.parse_value()?;
            if fields.insert(key.clone(), value).is_some() {
                return Err(format!("duplicate JSON key {key:?}"));
            }
            self.skip_whitespace();
            if self.take(b'}') {
                break;
            }
            self.expect(b',')?;
        }
        Ok(JsonValue::Object(fields))
    }

    fn parse_string(&mut self) -> Result<String, String> {
        self.expect(b'"')?;
        let mut value = String::new();
        loop {
            let byte = self
                .peek()
                .ok_or_else(|| "unterminated JSON string".to_owned())?;
            match byte {
                b'"' => {
                    self.position += 1;
                    return Ok(value);
                }
                b'\\' => {
                    self.position += 1;
                    let escape = self
                        .peek()
                        .ok_or_else(|| "unterminated JSON escape".to_owned())?;
                    self.position += 1;
                    match escape {
                        b'"' => value.push('"'),
                        b'\\' => value.push('\\'),
                        b'/' => value.push('/'),
                        b'b' => value.push('\u{0008}'),
                        b'f' => value.push('\u{000c}'),
                        b'n' => value.push('\n'),
                        b'r' => value.push('\r'),
                        b't' => value.push('\t'),
                        b'u' => value.push(self.parse_unicode_escape()?),
                        _ => return Err(format!("invalid JSON escape byte {escape}")),
                    }
                }
                0x00..=0x1f => {
                    let position = self.position;
                    return Err(format!(
                        "literal control byte {byte} in JSON string at {position}"
                    ));
                }
                _ => {
                    let ch = self.source[self.position..]
                        .chars()
                        .next()
                        .expect("position is within valid UTF-8");
                    value.push(ch);
                    self.position += ch.len_utf8();
                }
            }
        }
    }

    fn parse_hex_quad(&mut self) -> Result<u16, String> {
        let mut value = 0_u16;
        for _ in 0..4 {
            let byte = self
                .peek()
                .ok_or_else(|| "truncated JSON Unicode escape".to_owned())?;
            self.position += 1;
            let digit = match byte {
                b'0'..=b'9' => u16::from(byte - b'0'),
                b'a'..=b'f' => u16::from(byte - b'a' + 10),
                b'A'..=b'F' => u16::from(byte - b'A' + 10),
                _ => return Err(format!("non-hex byte {byte} in JSON Unicode escape")),
            };
            value = (value << 4) | digit;
        }
        Ok(value)
    }

    fn parse_unicode_escape(&mut self) -> Result<char, String> {
        let high = self.parse_hex_quad()?;
        let scalar = if (0xd800..=0xdbff).contains(&high) {
            self.expect(b'\\')?;
            self.expect(b'u')?;
            let low = self.parse_hex_quad()?;
            if !(0xdc00..=0xdfff).contains(&low) {
                return Err(format!("invalid low surrogate {low:#06x}"));
            }
            0x1_0000 + ((u32::from(high) - 0xd800) << 10) + (u32::from(low) - 0xdc00)
        } else if (0xdc00..=0xdfff).contains(&high) {
            return Err(format!("unpaired low surrogate {high:#06x}"));
        } else {
            u32::from(high)
        };
        char::from_u32(scalar).ok_or_else(|| format!("invalid Unicode scalar {scalar:#x}"))
    }

    fn parse_number(&mut self) -> Result<String, String> {
        let start = self.position;
        self.take(b'-');
        match self.peek() {
            Some(b'0') => self.position += 1,
            Some(b'1'..=b'9') => {
                self.position += 1;
                while matches!(self.peek(), Some(b'0'..=b'9')) {
                    self.position += 1;
                }
            }
            _ => return Err(format!("invalid JSON number at {start}")),
        }
        if self.take(b'.') {
            if !matches!(self.peek(), Some(b'0'..=b'9')) {
                let position = self.position;
                return Err(format!("fraction has no digit at {position}"));
            }
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.position += 1;
            }
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            self.position += 1;
            if matches!(self.peek(), Some(b'+' | b'-')) {
                self.position += 1;
            }
            if !matches!(self.peek(), Some(b'0'..=b'9')) {
                let position = self.position;
                return Err(format!("exponent has no digit at {position}"));
            }
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.position += 1;
            }
        }
        Ok(self.source[start..self.position].to_owned())
    }
}

fn object(value: &JsonValue) -> &BTreeMap<String, JsonValue> {
    match value {
        JsonValue::Object(fields) => fields,
        other => panic!("expected JSON object, found {other:?}"),
    }
}

fn string_field<'a>(fields: &'a BTreeMap<String, JsonValue>, key: &str) -> &'a str {
    match fields.get(key) {
        Some(JsonValue::String(value)) => value,
        other => panic!("expected string field {key:?}, found {other:?}"),
    }
}

fn number_field<'a>(fields: &'a BTreeMap<String, JsonValue>, key: &str) -> &'a str {
    match fields.get(key) {
        Some(JsonValue::Number(value)) => value,
        other => panic!("expected number field {key:?}, found {other:?}"),
    }
}

fn reduction_name(kind: ReductionKind) -> &'static str {
    match kind {
        ReductionKind::None => "None",
        ReductionKind::DeterministicSum => "DeterministicSum",
        ReductionKind::FastSum => "FastSum",
    }
}

fn determinism_name(class: DeterminismClass) -> &'static str {
    match class {
        DeterminismClass::BitwiseAllTiers => "BitwiseAllTiers",
        DeterminismClass::PerTier => "PerTier",
    }
}

fn assert_log_record_round_trips(record: &str, case: &str, verdict: &str, meta: &KernelMeta) {
    assert_eq!(record.lines().count(), 1, "one call must emit one record");
    assert!(
        !record
            .chars()
            .any(|ch| ch.is_control() || matches!(ch, '\u{2028}' | '\u{2029}')),
        "record contains a physical-line control: {record:?}"
    );
    let parsed = StrictJsonParser::new(record)
        .parse()
        .unwrap_or_else(|error| panic!("strict JSON refusal: {error}; record={record:?}"));
    let root = object(&parsed);
    assert_eq!(root.len(), 4, "outer schema drifted: {root:?}");
    assert_eq!(string_field(root, "suite"), "fs-tilelang");
    assert_eq!(string_field(root, "case"), case);
    assert_eq!(string_field(root, "verdict"), verdict);
    let detail = object(root.get("detail").expect("nested detail object"));
    assert_eq!(detail.len(), 7, "metadata schema drifted: {detail:?}");
    assert_eq!(string_field(detail, "kernel"), meta.name);
    let expected_flops = meta.flops_per_elem.to_string();
    assert_eq!(
        number_field(detail, "flops_per_elem"),
        expected_flops.as_str()
    );
    let expected_bytes = meta.bytes_per_elem.to_string();
    assert_eq!(
        number_field(detail, "bytes_per_elem"),
        expected_bytes.as_str()
    );
    let expected_intensity = format!("{:.4}", meta.intensity());
    assert_eq!(
        number_field(detail, "intensity"),
        expected_intensity.as_str()
    );
    let expected_halo = meta.halo.to_string();
    assert_eq!(number_field(detail, "halo"), expected_halo.as_str());
    assert_eq!(
        string_field(detail, "reduction"),
        reduction_name(meta.reduction)
    );
    assert_eq!(
        string_field(detail, "determinism"),
        determinism_name(meta.determinism)
    );
}

fn log(case: &str, verdict: &str, meta: &KernelMeta) {
    let record = meta
        .render_log_record(case, verdict)
        .expect("static battery metadata and labels must admit");
    assert_log_record_round_trips(&record, case, verdict, meta);
    println!("{record}");
}

fn rand_vec(n: usize, tile: u32) -> Vec<f64> {
    let mut s = StreamKey {
        seed: 77,
        kernel: 0x71E5,
        tile,
    }
    .stream();
    (0..n).map(|_| 2.0f64.mul_add(s.next_f64(), -1.0)).collect()
}

// Reference kernel 1: batched axpy-like (pure elementwise map).
kernel! {
    name: axpy_k,
    reads: [x, y],
    params: [alpha],
    writes: [out],
    reduction: none,
    body: {
        out = alpha.mul_add(x, y);
    },
}

// Reference kernel 2a: 1D 3-point stencil (literal shifts, halo 1) —
// gets the auto twin tests.
kernel! {
    name: stencil3_k,
    reads: [u],
    params: [c0, c1],
    writes: [out],
    halo: 1,
    reduction: none,
    body: {
        out = c0.mul_add(u, c1 * (shift_sub(u, 1) + shift_add(u, 1)));
    },
}

// Reference kernel 2b: 3D 7-point stencil via stride uparams (halo =
// one xy-plane); uparam kernels drive their own twin checks here.
kernel! {
    name: stencil7_k,
    reads: [u],
    uparams: [nx, nxy],
    params: [c0, c1],
    writes: [out],
    halo: nxy,
    reduction: none,
    body: {
        let ring = shift_sub(u, 1) + shift_add(u, 1)
            + shift_sub(u, nx) + shift_add(u, nx)
            + shift_sub(u, nxy) + shift_add(u, nxy);
        out = c0.mul_add(u, c1 * ring);
    },
}

// Reference kernel 3: SDF-style trilinear grid sample (gather form).
kernel! {
    name: trilinear_k,
    reads: [g, fx, fy, fz],
    index_reads: [ix, iy, iz],
    uparams: [nx, nxy],
    writes: [out],
    reduction: none,
    body: {
        let base = ix + nx * iy + nxy * iz;
        let c00 = gather(g, base) * (1.0 - fx) + gather(g, base + 1) * fx;
        let c10 = gather(g, base + nx) * (1.0 - fx) + gather(g, base + nx + 1) * fx;
        let c01 = gather(g, base + nxy) * (1.0 - fx) + gather(g, base + nxy + 1) * fx;
        let c11 = gather(g, base + nxy + nx) * (1.0 - fx) + gather(g, base + nxy + nx + 1) * fx;
        let c0 = c00 * (1.0 - fy) + c10 * fy;
        let c1 = c01 * (1.0 - fy) + c11 * fy;
        out = c0 * (1.0 - fz) + c1 * fz;
    },
}

// Reduction kernel: deterministic dot product.
kernel! {
    name: dot_k,
    reads: [x, y],
    writes: [],
    reduction: deterministic_sum,
    body: {
        acc = x * y;
    },
}

#[test]
fn axpy_matches_oracle_and_meta() {
    let n = 1537;
    let (x, y) = (rand_vec(n, 1), rand_vec(n, 2));
    let mut out = vec![0.0f64; n];
    axpy_k::run(&x, &y, 1.75, &mut out);
    for i in [0usize, n / 2, n - 1] {
        assert_eq!(
            out[i].to_bits(),
            1.75f64.mul_add(x[i], y[i]).to_bits(),
            "axpy oracle mismatch at {i}"
        );
    }
    assert_eq!(axpy_k::META.flops_per_elem, 2);
    assert_eq!(axpy_k::META.bytes_per_elem, 24);
    assert_eq!(axpy_k::META.reduction, ReductionKind::None);
    assert_eq!(axpy_k::META.determinism, DeterminismClass::BitwiseAllTiers);
    log("axpy", "pass", &axpy_k::META);
}

#[test]
fn stencil3_matches_oracle_and_halo() {
    let n = 513;
    let u = rand_vec(n, 3);
    let mut out = vec![f64::NAN; n];
    stencil3_k::run(&u, 0.5, 0.25, &mut out);
    // Halo untouched (still NaN), interior matches the oracle.
    assert!(
        out[0].is_nan() && out[n - 1].is_nan(),
        "halo must be untouched"
    );
    for i in 1..n - 1 {
        let expect = 0.5f64.mul_add(u[i], 0.25 * (u[i - 1] + u[i + 1]));
        assert_eq!(
            out[i].to_bits(),
            expect.to_bits(),
            "stencil3 mismatch at {i}"
        );
    }
    assert_eq!(stencil3_k::META.halo, 1);
    log("stencil3", "pass", &stencil3_k::META);
}

#[test]
fn stencil7_matches_oracle_and_tier_twins() {
    // 3D grid flattened: nx=12, ny=11, nz=9. Halo of one plane means
    // some non-interior cells get written too (wrap-reads within
    // bounds) — the ORACLE uses identical index arithmetic, so
    // equality is exact everywhere the kernel writes.
    let (nx, ny, nz) = (12usize, 11, 9);
    let nxy = nx * ny;
    let n = nxy * nz;
    let u = rand_vec(n, 4);
    let mut out = vec![0.0f64; n];
    stencil7_k::run(&u, nx, nxy, 0.4, 0.1, &mut out);
    let mut worst = 0u64;
    for i in nxy..n - nxy {
        let ring = u[i - 1] + u[i + 1] + u[i - nx] + u[i + nx] + u[i - nxy] + u[i + nxy];
        let expect = 0.4f64.mul_add(u[i], 0.1 * ring);
        worst = worst.max(out[i].to_bits() ^ expect.to_bits());
    }
    assert_eq!(worst, 0, "stencil7 oracle mismatch");
    // uparam kernels drive their own tier twins (macro can't guess
    // strides): all lane widths bitwise-equal to scalar.
    let mut out_s = vec![0.0f64; n];
    stencil7_k::run_scalar(&u, nx, nxy, 0.4, 0.1, &mut out_s);
    for w in [2usize, 4, 8] {
        let mut out_w = vec![0.0f64; n];
        match w {
            2 => stencil7_k::run_lanes::<2>(&u, nx, nxy, 0.4, 0.1, &mut out_w),
            4 => stencil7_k::run_lanes::<4>(&u, nx, nxy, 0.4, 0.1, &mut out_w),
            _ => stencil7_k::run_lanes::<8>(&u, nx, nxy, 0.4, 0.1, &mut out_w),
        }
        assert!(
            out_s
                .iter()
                .zip(&out_w)
                .all(|(a, b)| a.to_bits() == b.to_bits()),
            "stencil7 lane width {w} diverges from scalar"
        );
    }
    log("stencil7", "pass", &stencil7_k::META);
}

#[test]
fn trilinear_matches_oracle_and_tier_twins() {
    // Grid 8×7×6, 500 query points with in-range bases.
    let (nx, ny, nz) = (8usize, 7, 6);
    let nxy = nx * ny;
    let g = rand_vec(nxy * nz, 5);
    let m = 500usize;
    let mut s = StreamKey {
        seed: 78,
        kernel: 0x71E5,
        tile: 6,
    }
    .stream();
    let ix: Vec<u32> = (0..m)
        .map(|_| u32::try_from(s.next_below(nx as u64 - 1)).expect("small"))
        .collect();
    let iy: Vec<u32> = (0..m)
        .map(|_| u32::try_from(s.next_below(ny as u64 - 1)).expect("small"))
        .collect();
    let iz: Vec<u32> = (0..m)
        .map(|_| u32::try_from(s.next_below(nz as u64 - 1)).expect("small"))
        .collect();
    let fx: Vec<f64> = (0..m).map(|_| s.next_f64()).collect();
    let fy: Vec<f64> = (0..m).map(|_| s.next_f64()).collect();
    let fz: Vec<f64> = (0..m).map(|_| s.next_f64()).collect();
    let mut out = vec![0.0f64; m];
    trilinear_k::run(&g, &fx, &fy, &fz, &ix, &iy, &iz, nx, nxy, &mut out);
    // Hand-written oracle with identical arithmetic order.
    for q in [0usize, 250, 499] {
        let base = ix[q] as usize + nx * (iy[q] as usize) + nxy * (iz[q] as usize);
        let c00 = g[base] * (1.0 - fx[q]) + g[base + 1] * fx[q];
        let c10 = g[base + nx] * (1.0 - fx[q]) + g[base + nx + 1] * fx[q];
        let c01 = g[base + nxy] * (1.0 - fx[q]) + g[base + nxy + 1] * fx[q];
        let c11 = g[base + nxy + nx] * (1.0 - fx[q]) + g[base + nxy + nx + 1] * fx[q];
        let c0 = c00 * (1.0 - fy[q]) + c10 * fy[q];
        let c1 = c01 * (1.0 - fy[q]) + c11 * fy[q];
        let expect = c0 * (1.0 - fz[q]) + c1 * fz[q];
        assert_eq!(
            out[q].to_bits(),
            expect.to_bits(),
            "trilinear mismatch at {q}"
        );
    }
    // Gather kernels drive their own tier twins.
    let mut out_s = vec![0.0f64; m];
    trilinear_k::run_scalar(&g, &fx, &fy, &fz, &ix, &iy, &iz, nx, nxy, &mut out_s);
    let mut out_w = vec![0.0f64; m];
    trilinear_k::run_lanes::<4>(&g, &fx, &fy, &fz, &ix, &iy, &iz, nx, nxy, &mut out_w);
    assert!(
        out_s
            .iter()
            .zip(&out_w)
            .all(|(a, b)| a.to_bits() == b.to_bits()),
        "trilinear lanes diverge from scalar"
    );
    log("trilinear", "pass", &trilinear_k::META);
}

#[test]
fn dot_reduction_deterministic_and_tier_equal() {
    let n = 10_007;
    let (x, y) = (rand_vec(n, 7), rand_vec(n, 8));
    let d1 = dot_k::run(&x, &y);
    let d2 = dot_k::run_scalar(&x, &y);
    let d4 = dot_k::run_lanes::<4>(&x, &y);
    let d8 = dot_k::run_lanes::<8>(&x, &y);
    assert_eq!(d1.to_bits(), d2.to_bits(), "resolved-tier vs scalar");
    assert_eq!(d2.to_bits(), d4.to_bits(), "lanes 4 vs scalar");
    assert_eq!(d2.to_bits(), d8.to_bits(), "lanes 8 vs scalar");
    // Against the runtime's fixed-shape reference combiner.
    let prods: Vec<f64> = x.iter().zip(&y).map(|(a, b)| a * b).collect();
    assert_eq!(
        d2.to_bits(),
        fs_tilelang::deterministic_sum(&prods).to_bits(),
        "kernel reduction must equal the fixed-shape reference"
    );
    // Sanity vs a naive fold (envelope, not bitwise).
    let naive: f64 = prods.iter().sum();
    assert!((d2 - naive).abs() < 1e-9 * naive.abs().max(1.0));
    assert_eq!(dot_k::META.reduction, ReductionKind::DeterministicSum);
    log("dot", "pass", &dot_k::META);
}

#[test]
fn metadata_feeds_the_roofline_table() {
    // The per-kernel variant/intensity table (P6 evidence, logged).
    for meta in [
        axpy_k::META,
        stencil3_k::META,
        stencil7_k::META,
        trilinear_k::META,
        dot_k::META,
    ] {
        assert!(meta.flops_per_elem > 0, "{}: zero flops counted", meta.name);
        assert!(meta.bytes_per_elem > 0);
        assert!(meta.intensity() > 0.0);
        log("roofline-meta", "info", &meta);
    }
}

#[test]
fn metadata_json_escapes_kernel_names() {
    let plain = KernelMeta {
        name: "plain",
        flops_per_elem: 2,
        bytes_per_elem: 24,
        halo: 0,
        reduction: ReductionKind::None,
        determinism: DeterminismClass::BitwiseAllTiers,
    };
    assert_eq!(
        plain.descr(),
        "{\"kernel\":\"plain\",\"flops_per_elem\":2,\"bytes_per_elem\":24,\"intensity\":0.0833,\"halo\":0,\"reduction\":\"None\",\"determinism\":\"BitwiseAllTiers\"}"
    );

    const HOSTILE: &str = "kernel\"\\\u{0008}\u{000c}\n\r\t\u{0000}\u{001f}";
    const ESCAPED: &str = "kernel\\\"\\\\\\b\\f\\n\\r\\t\\u0000\\u001f";
    let hostile = KernelMeta {
        name: HOSTILE,
        ..plain
    };
    let descr = hostile.descr();
    assert_eq!(
        descr,
        format!(
            "{{\"kernel\":\"{ESCAPED}\",\"flops_per_elem\":2,\"bytes_per_elem\":24,\"intensity\":0.0833,\"halo\":0,\"reduction\":\"None\",\"determinism\":\"BitwiseAllTiers\"}}"
        )
    );
    assert_eq!(descr.lines().count(), 1, "{descr}");
    assert!(!descr.chars().any(|ch| ch < ' '), "{descr}");
}

#[test]
fn outer_logger_uses_a_nested_metadata_object_with_exact_stable_bytes() {
    let meta = KernelMeta {
        name: "plain",
        flops_per_elem: 2,
        bytes_per_elem: 24,
        halo: 0,
        reduction: ReductionKind::None,
        determinism: DeterminismClass::BitwiseAllTiers,
    };
    let record = meta
        .render_log_record("roofline-meta", "info")
        .expect("ordinary record admits");
    assert_eq!(
        record,
        "{\"suite\":\"fs-tilelang\",\"case\":\"roofline-meta\",\"verdict\":\"info\",\"detail\":{\"kernel\":\"plain\",\"flops_per_elem\":2,\"bytes_per_elem\":24,\"intensity\":0.0833,\"halo\":0,\"reduction\":\"None\",\"determinism\":\"BitwiseAllTiers\"}}"
    );
    assert!(
        !record.contains("\"detail\":\""),
        "nested JSON must never be quoted as a detail string"
    );
    assert_log_record_round_trips(&record, "roofline-meta", "info", &meta);
}

#[test]
fn exhaustive_hostile_metadata_and_labels_remain_valid_and_round_trip() {
    let mut hostile = "prefix\"\\".to_owned();
    hostile.extend((0_u32..=31).map(|value| char::from_u32(value).expect("C0 scalar")));
    hostile.extend((0x7f_u32..=0x9f).map(|value| char::from_u32(value).expect("C1 scalar")));
    hostile.push('\u{2028}');
    hostile.push('\u{2029}');
    hostile.push('é');
    hostile.push('🦀');
    let static_name: &'static str = Box::leak(hostile.clone().into_boxed_str());
    let meta = KernelMeta {
        name: static_name,
        flops_per_elem: u32::MAX,
        bytes_per_elem: 1,
        halo: u32::MAX,
        reduction: ReductionKind::FastSum,
        determinism: DeterminismClass::PerTier,
    };

    let inner = meta.try_descr().expect("bounded hostile metadata admits");
    let parsed_inner = StrictJsonParser::new(&inner)
        .parse()
        .expect("inner metadata is strict JSON");
    assert_eq!(
        string_field(object(&parsed_inner), "kernel"),
        hostile.as_str()
    );
    assert_eq!(inner.lines().count(), 1);
    for value in 0x7f_u32..=0x9f {
        let escaped = format!("\\u{value:04x}");
        assert_eq!(
            inner.matches(&escaped).count(),
            1,
            "inner metadata must escape C1 scalar U+{value:04X} exactly once"
        );
    }

    let record = meta
        .render_log_record(&hostile, &hostile)
        .expect("bounded hostile record admits");
    assert_log_record_round_trips(&record, &hostile, &hostile, &meta);
    assert_eq!(record.lines().count(), 1);
    assert!(record.len() <= MAX_LOG_RECORD_BYTES);
    for value in 0x7f_u32..=0x9f {
        let escaped = format!("\\u{value:04x}");
        assert_eq!(
            record.matches(&escaped).count(),
            3,
            "kernel, case, and verdict must each escape C1 scalar U+{value:04X}"
        );
    }
}

#[test]
fn metadata_and_log_admission_enforce_every_byte_boundary_before_rendering() {
    let max_name: &'static str = Box::leak("\0".repeat(MAX_KERNEL_NAME_BYTES).into_boxed_str());
    let meta = KernelMeta {
        name: max_name,
        flops_per_elem: u32::MAX,
        bytes_per_elem: 1,
        halo: u32::MAX,
        reduction: ReductionKind::DeterministicSum,
        determinism: DeterminismClass::BitwiseAllTiers,
    };
    let inner = meta.try_descr().expect("inclusive name limit admits");
    assert!(inner.len() <= MAX_METADATA_JSON_BYTES);

    let max_label = "\0".repeat(MAX_LOG_LABEL_BYTES);
    let record = meta
        .render_log_record(&max_label, &max_label)
        .expect("inclusive label limits admit after worst-case escaping");
    assert!(record.len() <= MAX_LOG_RECORD_BYTES);
    assert_log_record_round_trips(&record, &max_label, &max_label, &meta);

    let over_name: &'static str = Box::leak("x".repeat(MAX_KERNEL_NAME_BYTES + 1).into_boxed_str());
    let over_meta = KernelMeta {
        name: over_name,
        ..meta
    };
    assert_eq!(
        over_meta.try_descr(),
        Err(MetadataRenderError::TextTooLong {
            field: "kernel",
            actual: MAX_KERNEL_NAME_BYTES + 1,
            limit: MAX_KERNEL_NAME_BYTES,
        })
    );
    assert_eq!(
        over_meta.render_log_record("case", "pass"),
        Err(MetadataRenderError::TextTooLong {
            field: "kernel",
            actual: MAX_KERNEL_NAME_BYTES + 1,
            limit: MAX_KERNEL_NAME_BYTES,
        })
    );
    let refusal = over_meta.descr();
    assert!(refusal.len() <= MAX_METADATA_JSON_BYTES);
    assert!(
        !refusal.contains(over_name),
        "bounded refusal must not echo rejected attacker text"
    );
    let refusal = StrictJsonParser::new(&refusal)
        .parse()
        .expect("compatibility refusal is strict bounded JSON");
    let refusal = object(&refusal);
    assert_eq!(string_field(refusal, "kernel_metadata"), "refused");
    assert_eq!(string_field(refusal, "field"), "kernel");

    let over_label = "x".repeat(MAX_LOG_LABEL_BYTES + 1);
    assert_eq!(
        meta.render_log_record(&over_label, "pass"),
        Err(MetadataRenderError::TextTooLong {
            field: "case",
            actual: MAX_LOG_LABEL_BYTES + 1,
            limit: MAX_LOG_LABEL_BYTES,
        })
    );
    assert_eq!(
        meta.render_log_record("case", &over_label),
        Err(MetadataRenderError::TextTooLong {
            field: "verdict",
            actual: MAX_LOG_LABEL_BYTES + 1,
            limit: MAX_LOG_LABEL_BYTES,
        })
    );
    assert_eq!(
        KernelMeta { name: "", ..meta }.try_descr(),
        Err(MetadataRenderError::EmptyText { field: "kernel" })
    );
    assert_eq!(
        meta.render_log_record("", "pass"),
        Err(MetadataRenderError::EmptyText { field: "case" })
    );
    assert_eq!(
        meta.render_log_record("case", ""),
        Err(MetadataRenderError::EmptyText { field: "verdict" })
    );
}

#[test]
fn independent_strict_parser_rejects_duplicate_trailing_and_literal_control_attacks() {
    for malformed in [
        "{\"field\":1,\"field\":2}",
        "{} trailing",
        "{\"field\":\"literal\nnewline\"}",
        "{\"field\":01}",
        "{\"field\":\"\\uD800\"}",
    ] {
        assert!(
            StrictJsonParser::new(malformed).parse().is_err(),
            "strict parser admitted {malformed:?}"
        );
    }
}
