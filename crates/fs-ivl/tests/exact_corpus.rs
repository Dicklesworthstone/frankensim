//! Exact rational/algebraic corpus for certified arithmetic
//! (bead frankensim-extreal-program-f85xj.3.2).
//!
//! Arbitrary-precision oracles are still software; exact cases are not.
//! Every row's ground truth is constructed by ALL-INTEGER (i128) dyadic
//! arithmetic or by hand-derivable lattice/algebraic identities, so a
//! failure indicts the implementation (or the corpus) with no third
//! implementation in the loop — the property the MPFR lane (.3.1) cannot
//! give by construction.
//!
//! Shape (mirrors the tracked reference-project pattern): the corpus text
//! is a COMMITTED fixture (`tests/fixtures/exact_corpus.txt`) whose last
//! line carries the row count and an FNV-1a-64 content address; the
//! generator is the ignored `dump_exact_corpus` test in this file; a drift
//! guard pins tracked bytes == generator output; the battery decodes the
//! tracked bytes and executes every row. Decoder refusals (malformed,
//! unknown primitive, duplicate id, bad integer, truncation, count/hash
//! mismatch, oversize) are typed and tested; mutation drills corrupt a
//! numerator, an endpoint order, and an expected classification and must
//! fail. The `end` line's hash is the content-addressed corpus root that
//! .3.7's terminal campaign consumes.
//!
//! Exactness discipline: a dyadic `num·2^k` is admitted into a row only if
//! its f64 conversion reconstructs bit-exactly (checked by integer
//! cross-shifting, no float comparison). Width gates are in ULPs on the
//! ordered-bit line. Cancellation/resource fault classes are N/A by type
//! for this in-process test corpus (no Cx, no I/O beyond one bounded
//! tracked read, no receipts); the elaborate injection points named in the
//! acceptance template belong to the .3.7 terminal harness this corpus
//! feeds. This corpus proves its stated rows only.

use fs_casebook::fnv1a64;
use fs_ivl::predicates::{Sign, incircle, orient2d, orient3d};
use fs_ivl::{Interval, TaylorModel1, expansion, newton};

const CORPUS_REL: &str = "tests/fixtures/exact_corpus.txt";
const MAX_LINES: usize = 4096;
const MAX_LINE_BYTES: usize = 256;

// ---------------------------------------------------------------------------
// Exact dyadic machinery (all-integer)
// ---------------------------------------------------------------------------

/// 2^k as f64 for k in [-1022, 1023] (normal range, exact by construction).
fn pow2(k: i32) -> f64 {
    assert!((-1022..=1023).contains(&k), "pow2 exponent {k}");
    f64::from_bits((u64::try_from(k + 1023).expect("normal")) << 52)
}

/// `num · 2^k` as f64, PROVEN exact by integer reconstruction: the
/// resulting float's (mantissa, exponent) must cross-shift back onto
/// (num, k) in i128 — no float comparisons anywhere.
fn exact_dyadic(num: i128, k: i32) -> f64 {
    assert!(num.unsigned_abs() < (1u128 << 53), "mantissa width: {num}");
    #[allow(clippy::cast_precision_loss)] // |num| < 2^53: exact
    let mut value = num as f64;
    // Apply 2^k in normal-range steps so subnormal targets are reachable;
    // any precision loss on the final subnormal step is caught by the
    // reconstruction check below.
    let mut remaining = k;
    while remaining != 0 {
        let step = remaining.clamp(-1000, 1000);
        value *= pow2(step);
        remaining -= step;
    }
    if num == 0 {
        assert!(value == 0.0);
        return value;
    }
    // Reconstruct: value = m·2^e exactly from the bit pattern.
    let bits = value.to_bits();
    let exp_field = (bits >> 52) & 0x7ff;
    assert!(exp_field != 0x7ff, "dyadic overflowed: {num}*2^{k}");
    let (m, e) = if exp_field == 0 {
        (i128::from(bits & ((1u64 << 52) - 1)), -1074i32)
    } else {
        (
            i128::from((bits & ((1u64 << 52) - 1)) | (1u64 << 52)),
            i32::try_from(exp_field).expect("small") - 1075,
        )
    };
    let m = if bits >> 63 == 1 { -m } else { m };
    // m·2^e == num·2^k  ⇔  m << (e−k) == num (or num << (k−e) == m).
    let equal = if e >= k {
        let shift = u32::try_from(e - k).expect("bounded");
        shift < 74 && m.checked_shl(shift) == Some(num)
    } else {
        let shift = u32::try_from(k - e).expect("bounded");
        shift < 74 && num.checked_shl(shift) == Some(m)
    };
    assert!(equal, "inexact dyadic admission: {num}*2^{k} -> {value:e}");
    value
}

/// Width of an interval in ULP steps on the ordered-bit line.
fn width_ulps(iv: Interval) -> u64 {
    fn ordered(x: f64) -> i64 {
        #[allow(clippy::cast_possible_wrap)] // the wrap IS the ordering trick
        let b = x.to_bits() as i64;
        if b < 0 { i64::MIN ^ b } else { b }
    }
    ordered(iv.hi()).abs_diff(ordered(iv.lo()))
}

// ---------------------------------------------------------------------------
// Corpus format
// ---------------------------------------------------------------------------

/// Typed corpus-load refusals; Display carries stable codes.
#[derive(Debug)]
enum CorpusError {
    Malformed(usize, &'static str),
    UnknownPrimitive(usize, String),
    DuplicateId(usize, String),
    BadInteger(usize, String),
    Truncated,
    CountMismatch { declared: usize, decoded: usize },
    RootMismatch { declared: String, computed: String },
    Oversize(&'static str),
}

impl core::fmt::Display for CorpusError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            CorpusError::Malformed(line, what) => {
                write!(f, "FS-IVL-CORPUS-MALFORMED: line {line}: {what}")
            }
            CorpusError::UnknownPrimitive(line, p) => {
                write!(f, "FS-IVL-CORPUS-UNKNOWN-PRIMITIVE: line {line}: `{p}`")
            }
            CorpusError::DuplicateId(line, id) => {
                write!(f, "FS-IVL-CORPUS-DUPLICATE-ID: line {line}: `{id}`")
            }
            CorpusError::BadInteger(line, field) => {
                write!(f, "FS-IVL-CORPUS-BAD-INTEGER: line {line}: `{field}`")
            }
            CorpusError::Truncated => write!(f, "FS-IVL-CORPUS-TRUNCATED: no `end` line"),
            CorpusError::CountMismatch { declared, decoded } => write!(
                f,
                "FS-IVL-CORPUS-COUNT-MISMATCH: end declares {declared}, decoded {decoded}"
            ),
            CorpusError::RootMismatch { declared, computed } => write!(
                f,
                "FS-IVL-CORPUS-ROOT-MISMATCH: end declares {declared}, computed {computed}"
            ),
            CorpusError::Oversize(what) => write!(f, "FS-IVL-CORPUS-OVERSIZE: {what}"),
        }
    }
}

#[derive(Debug, Clone)]
struct Row {
    id: String,
    primitive: String,
    fields: Vec<i128>,
}

/// Decode the corpus text: comment/blank lines skipped, `row` lines
/// collected, exactly one final `end|count|root` line verified against the
/// preceding bytes. Returns (rows, corpus root).
fn decode(text: &str) -> Result<(Vec<Row>, String), CorpusError> {
    if text.lines().count() > MAX_LINES {
        return Err(CorpusError::Oversize("line count"));
    }
    let mut rows: Vec<Row> = Vec::new();
    let mut consumed = 0usize;
    for (index, line) in text.lines().enumerate() {
        let line_no = index + 1;
        if line.len() > MAX_LINE_BYTES {
            return Err(CorpusError::Oversize("line length"));
        }
        if line.is_empty() || line.starts_with('#') {
            consumed += line.len() + 1;
            continue;
        }
        let parts: Vec<&str> = line.split('|').collect();
        match parts[0] {
            "row" => {
                if parts.len() < 4 {
                    return Err(CorpusError::Malformed(line_no, "row needs id|prim|fields"));
                }
                let id = parts[1].to_string();
                if rows.iter().any(|r| r.id == id) {
                    return Err(CorpusError::DuplicateId(line_no, id));
                }
                let mut fields = Vec::with_capacity(parts.len() - 3);
                for raw in &parts[3..] {
                    fields.push(
                        raw.parse::<i128>()
                            .map_err(|_| CorpusError::BadInteger(line_no, (*raw).to_string()))?,
                    );
                }
                rows.push(Row {
                    id,
                    primitive: parts[2].to_string(),
                    fields,
                });
                consumed += line.len() + 1;
            }
            "end" => {
                if parts.len() != 3 {
                    return Err(CorpusError::Malformed(line_no, "end needs count|root"));
                }
                let declared: usize = parts[1]
                    .parse()
                    .map_err(|_| CorpusError::BadInteger(line_no, parts[1].to_string()))?;
                if declared != rows.len() {
                    return Err(CorpusError::CountMismatch {
                        declared,
                        decoded: rows.len(),
                    });
                }
                let computed = format!("{:016x}", fnv1a64(&text.as_bytes()[..consumed]));
                if computed != parts[2] {
                    return Err(CorpusError::RootMismatch {
                        declared: parts[2].to_string(),
                        computed,
                    });
                }
                return Ok((rows, computed));
            }
            other => {
                return Err(CorpusError::UnknownPrimitive(line_no, other.to_string()));
            }
        }
    }
    Err(CorpusError::Truncated)
}

// ---------------------------------------------------------------------------
// Row execution
// ---------------------------------------------------------------------------

fn sign_of(v: i128) -> Sign {
    match v.cmp(&0) {
        core::cmp::Ordering::Less => Sign::Negative,
        core::cmp::Ordering::Equal => Sign::Zero,
        core::cmp::Ordering::Greater => Sign::Positive,
    }
}

/// Execute one row; Err carries id, primitive, expected and actual
/// relations (bounded rendering, exact rationals spelled `num*2^k`).
#[allow(clippy::too_many_lines)] // one dispatch arm per corpus primitive
fn run_row(row: &Row) -> Result<(), String> {
    let f = &row.fields;
    let need = |n: usize| -> Result<(), String> {
        if f.len() == n {
            Ok(())
        } else {
            Err(format!(
                "{}: {} needs {n} fields, has {}",
                row.id,
                row.primitive,
                f.len()
            ))
        }
    };
    let int_arg = |v: i128| -> Result<i32, String> {
        i32::try_from(v).map_err(|_| format!("{}: exponent field {v} out of range", row.id))
    };
    match row.primitive.as_str() {
        // point-op|a_num|a_k|b_num|b_k|r_num|r_k: [a,a] op [b,b] must
        // contain the exactly representable result with width <= 2 ULPs.
        "add" | "sub" | "mul" | "div" => {
            need(6)?;
            let a = exact_dyadic(f[0], int_arg(f[1])?);
            let b = exact_dyadic(f[2], int_arg(f[3])?);
            let expected = exact_dyadic(f[4], int_arg(f[5])?);
            let (ia, ib) = (Interval::point(a), Interval::point(b));
            let result = match row.primitive.as_str() {
                "add" => ia + ib,
                "sub" => ia - ib,
                "mul" => ia * ib,
                _ => ia / ib,
            };
            if !result.contains(expected) {
                return Err(format!(
                    "{}: {} [{:e},{:e}] excludes exact {}*2^{}",
                    row.id,
                    row.primitive,
                    result.lo(),
                    result.hi(),
                    f[4],
                    f[5]
                ));
            }
            let width = width_ulps(result);
            if width > 2 {
                return Err(format!(
                    "{}: {} width {width} ULPs exceeds the 2-ULP outward gate",
                    row.id, row.primitive
                ));
            }
            Ok(())
        }
        // sqrt-exact|m|k: sqrt of (m^2 · 2^{2k}) must enclose m·2^k tightly.
        "sqrt-exact" => {
            need(2)?;
            let m = f[0];
            let k = int_arg(f[1])?;
            let square = exact_dyadic(m.checked_mul(m).expect("fits"), 2 * k);
            let expected = exact_dyadic(m, k);
            let result = Interval::point(square).sqrt();
            if !result.contains(expected) || width_ulps(result) > 2 {
                return Err(format!(
                    "{}: sqrt-exact [{:e},{:e}] vs {m}*2^{k} (width {})",
                    row.id,
                    result.lo(),
                    result.hi(),
                    width_ulps(result)
                ));
            }
            Ok(())
        }
        // sqrt-witness|n: enclosure of sqrt(n) must satisfy lo^2 <= n <= hi^2
        // checked EXACTLY in integers (the independent algebraic witness).
        "sqrt-witness" => {
            need(1)?;
            let n = f[0];
            #[allow(clippy::cast_precision_loss)] // small integer n
            let result = Interval::point(n as f64).sqrt();
            let square_vs_n = |x: f64| -> core::cmp::Ordering {
                // x = m·2^e exactly; compare m^2·2^{2e} with n·2^0.
                let bits = x.to_bits();
                let exp_field = (bits >> 52) & 0x7ff;
                assert!(exp_field != 0 && exp_field != 0x7ff && bits >> 63 == 0);
                let m = i128::from((bits & ((1u64 << 52) - 1)) | (1u64 << 52));
                let e = i32::try_from(exp_field).expect("small") - 1075;
                let m2 = m.checked_mul(m).expect("106 bits");
                // 2e is negative here (sqrt of a small integer): n << (-2e).
                let shift = u32::try_from(-2 * e).expect("negative exponent");
                m2.cmp(&(n.checked_shl(shift).expect("small n")))
            };
            if square_vs_n(result.lo()) == core::cmp::Ordering::Greater {
                return Err(format!("{}: sqrt-witness lo^2 > {n}", row.id));
            }
            if square_vs_n(result.hi()) == core::cmp::Ordering::Less {
                return Err(format!("{}: sqrt-witness hi^2 < {n}", row.id));
            }
            Ok(())
        }
        // orient2d|ax|ay|bx|by|cx|cy|sign on the integer lattice.
        "orient2d" => {
            need(7)?;
            #[allow(clippy::cast_precision_loss)] // lattice coords < 2^53
            let p = |i: usize| [f[i] as f64, f[i + 1] as f64];
            let got = orient2d(p(0), p(2), p(4));
            let want = sign_of(f[6]);
            if got != want {
                return Err(format!("{}: orient2d {got:?}, expected {want:?}", row.id));
            }
            Ok(())
        }
        // orient3d|12 coords|sign.
        "orient3d" => {
            need(13)?;
            #[allow(clippy::cast_precision_loss)]
            let p = |i: usize| [f[i] as f64, f[i + 1] as f64, f[i + 2] as f64];
            let got = orient3d(p(0), p(3), p(6), p(9));
            let want = sign_of(f[12]);
            if got != want {
                return Err(format!("{}: orient3d {got:?}, expected {want:?}", row.id));
            }
            Ok(())
        }
        // incircle|8 coords|sign.
        "incircle" => {
            need(9)?;
            #[allow(clippy::cast_precision_loss)]
            let p = |i: usize| [f[i] as f64, f[i + 1] as f64];
            let got = incircle(p(0), p(2), p(4), p(6));
            let want = sign_of(f[8]);
            if got != want {
                return Err(format!("{}: incircle {got:?}, expected {want:?}", row.id));
            }
            Ok(())
        }
        // proddiff|a|b|c|d|sign: sign(a*b - c*d) through the expansion path.
        "proddiff" => {
            need(5)?;
            #[allow(clippy::cast_precision_loss)]
            let e = expansion::prod_diff(f[0] as f64, f[1] as f64, f[2] as f64, f[3] as f64);
            let got = expansion::expansion_sign(&e);
            let want = i32::try_from(f[4].signum()).expect("sign");
            if got != want {
                return Err(format!("{}: proddiff sign {got}, expected {want}", row.id));
            }
            Ok(())
        }
        // newton-root|r_num|r_k: f(x) = (x - r)(x + 2) on [-4, 4]; a
        // certified box must contain the exact rational root r.
        "newton-root" => {
            need(2)?;
            let r = exact_dyadic(f[0], int_arg(f[1])?);
            let fr = move |x: Interval| (x - Interval::point(r)) * (x + Interval::point(2.0));
            let fp = move |x: Interval| x * Interval::point(2.0) + Interval::point(2.0 - r);
            let boxes = newton::newton_roots(&fr, &fp, Interval::new(-1.0, 1.5), 1e-12);
            let hit = boxes.iter().any(|b| b.interval().contains(r));
            if !hit {
                return Err(format!(
                    "{}: newton-root no certified box contains {}*2^{} ({} boxes)",
                    row.id,
                    f[0],
                    f[1],
                    boxes.len()
                ));
            }
            Ok(())
        }
        // newton-noroot|c: f(x) = x^2 + c (c > 0) on [0, 2] has no root and
        // must certify none.
        "newton-noroot" => {
            need(1)?;
            let c = exact_dyadic(f[0], 0);
            let fr = move |x: Interval| x * x + Interval::point(c);
            let fp = move |x: Interval| x * Interval::new(2.0, 2.0);
            let boxes = newton::newton_roots(&fr, &fp, Interval::new(0.0, 2.0), 1e-12);
            if !boxes.is_empty() {
                return Err(format!(
                    "{}: newton-noroot returned {} boxes for x^2+{c}",
                    row.id,
                    boxes.len()
                ));
            }
            Ok(())
        }
        // taylor-sq|r_num|r_k: T(x) = x·x − x on [0, 2] evaluated at the
        // exact rational r must enclose r^2 − r (exact dyadic, all-integer).
        "taylor-sq" => {
            need(2)?;
            let r_num = f[0];
            let r_k = int_arg(f[1])?;
            let r = exact_dyadic(r_num, r_k);
            // r^2 − r = 2^{2k}·(r_num² − r_num·2^{−k}); corpus rows keep
            // k ≤ 0 so the inner shift is a plain left shift in i128.
            assert!(r_k <= 0, "taylor-sq rows use non-positive exponents");
            let exact_num = r_num
                .checked_mul(r_num)
                .and_then(|sq| {
                    sq.checked_sub(r_num.checked_shl(u32::try_from(-r_k).expect("k <= 0"))?)
                })
                .expect("fits");
            let expected = exact_dyadic(exact_num, 2 * r_k);
            let x = TaylorModel1::variable(Interval::new(0.0, 2.0), 2).expect("variable model");
            let x_squared = (&x * &x).expect("square model");
            let model = (&x_squared - &x).expect("difference model");
            let result = model.eval_interval(Interval::point(r));
            if !result.contains(expected) {
                return Err(format!(
                    "{}: taylor-sq [{:e},{:e}] excludes exact {exact_num}*2^{}",
                    row.id,
                    result.lo(),
                    result.hi(),
                    2 * r_k
                ));
            }
            Ok(())
        }
        // special-inf|which: extended-real singletons keep their infinite
        // endpoint through addition with a finite point.
        "special-inf" => {
            need(1)?;
            let inf = if f[0] >= 0 {
                f64::INFINITY
            } else {
                f64::NEG_INFINITY
            };
            let result = Interval::point(inf) + Interval::point(1.0);
            let ok = if f[0] >= 0 {
                result.hi() == f64::INFINITY
            } else {
                result.lo() == f64::NEG_INFINITY
            };
            if !ok {
                return Err(format!(
                    "{}: special-inf lost its infinite endpoint",
                    row.id
                ));
            }
            Ok(())
        }
        other => Err(format!(
            "{}: unknown primitive `{other}` at run time",
            row.id
        )),
    }
}

// ---------------------------------------------------------------------------
// The generator (single source of truth for the tracked fixture)
// ---------------------------------------------------------------------------

/// Build the canonical corpus text, `end` line included. Deterministic,
/// all-literal: this IS the provenance of every tracked byte.
#[allow(clippy::too_many_lines)] // a corpus is a list; splitting it obscures it
fn built_corpus() -> String {
    let mut body = String::new();
    body.push_str(
        "# frankensim exact-arithmetic corpus v1 (bead f85xj.3.2).\n\
         # GENERATED by the ignored test `dump_exact_corpus` in\n\
         # crates/fs-ivl/tests/exact_corpus.rs — edit the generator, never\n\
         # this file. Ground truth per row is all-integer dyadic arithmetic\n\
         # or hand-derivable lattice/algebraic identities; the final line is\n\
         # `end|<rows>|<fnv1a64 of preceding bytes>` (the corpus root).\n",
    );
    let mut n = 0usize;
    let mut row = |primitive: &str, fields: &[i128]| {
        n += 1;
        let joined = fields
            .iter()
            .map(i128::to_string)
            .collect::<Vec<_>>()
            .join("|");
        body.push_str(&format!("row|r{n:03}|{primitive}|{joined}\n"));
    };

    // (a) Exact dyadic point chains — reducible spellings included
    // (6·2^-3 = 3·2^-2), signed zeros, subnormal and near-overflow scales.
    row("add", &[3, -2, 1, -2, 1, 0]); // 3/4 + 1/4 = 1
    row("add", &[6, -3, 1, -2, 1, 0]); // reducible 6/8 + 1/4 = 1
    row("add", &[0, 0, 0, 0, 0, 0]); // +0 + +0 = 0
    row("add", &[5, -1074, 3, -1074, 8, -1074]); // subnormal + subnormal
    row("add", &[-7, 40, 7, 40, 0, 0]); // exact cancellation to zero
    row("sub", &[1, 0, 1, -53, 9007199254740991, -53]); // 1 − 2^-53 (Sterbenz-adjacent)
    row("sub", &[3, -2, 3, -2, 0, 0]); // x − x = 0
    row("sub", &[-5, -1074, -5, -1074, 0, 0]); // subnormal self-cancel
    row("mul", &[3, -2, 5, -1, 15, -3]); // 3/4 · 5/2 = 15/8
    row("mul", &[-9, 10, 7, -12, -63, -2]); // sign mix, exponent mix
    row("mul", &[4503599627370497, -52, 2, 52, 4503599627370497, 1]); // full-width mantissa · 2^52
    row("mul", &[1, 511, 1, 511, 1, 1022]); // near-overflow scale, exact
    row("div", &[15, -3, 5, -1, 3, -2]); // (15/8)/(5/2) = 3/4
    row("div", &[-63, -2, 7, -12, -9, 10]); // inverse of the mul row
    row("div", &[1, 0, 1, 10, 1, -10]); // 1/1024
    row("div", &[9, -1074, 3, 0, 3, -1074]); // subnormal / small integer

    // (b) sqrt of perfect squares and integer-square dyadics.
    row("sqrt-exact", &[3, -2]); // sqrt(9/16) = 3/4
    row("sqrt-exact", &[1, -537]); // sqrt(2^-1074) = 2^-537
    row("sqrt-exact", &[67108859, -26]); // near-2^26 mantissa square

    // sqrt witnesses: lo² ≤ n ≤ hi² checked exactly in integers.
    for n_val in [2i128, 3, 5, 7, 10] {
        row("sqrt-witness", &[n_val]);
    }

    // (c) Newton/Krawczyk certified rational roots and refusals.
    row("newton-root", &[3, -2]); // root 3/4
    row("newton-root", &[-1, -3]); // root −1/8
    row("newton-root", &[1023, -10]); // root 1023/1024
    row("newton-noroot", &[1]); // x² + 1 on [0,2]
    row("newton-noroot", &[3]); // x² + 3 on [0,2]

    // Taylor model P(x) = x² − x at exact rationals: remainder must not
    // exclude the exact dyadic value.
    row("taylor-sq", &[3, -2]); // P(3/4) = −3/16
    row("taylor-sq", &[1, -1]); // P(1/2) = −1/4 (the minimum)
    row("taylor-sq", &[7, -2]); // P(7/4) = 21/16

    // (d) Exact predicates on the integer lattice, degeneracies included.
    row("orient2d", &[0, 0, 1, 0, 0, 1, 1]); // CCW
    row("orient2d", &[0, 0, 0, 1, 1, 0, -1]); // CW
    row("orient2d", &[0, 0, 2, 1, 4, 2, 0]); // exactly collinear
    row(
        "orient2d",
        &[1048576, 1048576, 2097152, 2097153, 3145728, 3145730, 0],
    ); // large-lattice collinear
    row(
        "orient2d",
        &[1048576, 1048576, 2097152, 2097153, 3145728, 3145731, 1],
    ); // one off the line
    row("orient3d", &[0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, -1]); // positively oriented tetra (Shewchuk sign convention pinned)
    row("orient3d", &[0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0]); // coplanar
    row("incircle", &[0, 0, 2, 0, 0, 2, 1, 1, 1]); // strictly inside
    row("incircle", &[0, 0, 2, 0, 0, 2, 2, 2, 0]); // exactly cocircular
    row("incircle", &[0, 0, 2, 0, 0, 2, 3, 3, -1]); // strictly outside

    // (e) Expansion identities: sign(a·b − c·d) with exact integer truth.
    row("proddiff", &[3, 7, 21, 1, 0]); // 21 − 21 = 0
    row("proddiff", &[12345, 6789, 83810205, 1, 0]); // exact large product tie
    row("proddiff", &[12345, 6789, 83810204, 1, 1]); // one below the tie
    row("proddiff", &[12345, 6789, 83810206, 1, -1]); // one above the tie
    row("proddiff", &[4503599627370497, 2, 4503599627370496, 2, 1]); // adjacent 53-bit mantissas

    // Extended-real singletons.
    row("special-inf", &[1]);
    row("special-inf", &[-1]);

    let root = fnv1a64(body.as_bytes());
    format!("{body}end|{n}|{root:016x}\n")
}

fn corpus_path() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(CORPUS_REL)
}

/// Generator: writes the tracked fixture. Run explicitly after editing
/// `built_corpus`:
/// `cargo test -p fs-ivl --test exact_corpus -- --ignored dump_exact_corpus`
#[test]
#[ignore = "generator: writes the tracked corpus fixture"]
fn dump_exact_corpus() {
    let path = corpus_path();
    std::fs::create_dir_all(path.parent().expect("fixtures dir")).expect("mkdir");
    std::fs::write(&path, built_corpus()).expect("corpus writes");
    println!("wrote {}", path.display());
}

// ---------------------------------------------------------------------------
// Batteries
// ---------------------------------------------------------------------------

#[test]
fn g0_tracked_corpus_matches_its_generator() {
    let tracked = std::fs::read_to_string(corpus_path()).expect("tracked corpus exists");
    assert_eq!(
        tracked,
        built_corpus(),
        "tracked corpus drifted from its generator; re-run the ignored dump_exact_corpus test"
    );
}

#[test]
fn g0_every_corpus_row_holds_exactly() {
    let tracked = std::fs::read_to_string(corpus_path()).expect("tracked corpus exists");
    let (rows, root) = decode(&tracked).expect("tracked corpus decodes");
    assert!(rows.len() >= 40, "corpus is not vacuous: {}", rows.len());
    let mut failures = Vec::new();
    for r in &rows {
        if let Err(what) = run_row(r) {
            failures.push(what);
        }
    }
    assert!(failures.is_empty(), "corpus rows failed: {failures:#?}");
    println!(
        "{{\"suite\":\"fs-ivl-exact-corpus\",\"case\":\"battery\",\"corpus_root\":\"{root}\",\
         \"rows\":{},\"failures\":0,\"verdict\":\"pass\"}}",
        rows.len()
    );
}

#[test]
fn g3_decoder_refuses_malformed_duplicate_truncated_and_oversized_loads() {
    let good = built_corpus();
    let cases: Vec<(&str, String, &str)> = vec![
        (
            "truncated",
            good.lines()
                .take(good.lines().count() - 1)
                .collect::<Vec<_>>()
                .join("\n"),
            "FS-IVL-CORPUS-TRUNCATED",
        ),
        (
            "unknown-line-tag",
            good.replace("row|r001|add|", "frobnicate|r001|add|"),
            "FS-IVL-CORPUS-UNKNOWN-PRIMITIVE",
        ),
        (
            "duplicate-id",
            good.replace("row|r002|", "row|r001|"),
            "FS-IVL-CORPUS-DUPLICATE-ID",
        ),
        (
            "bad-integer",
            good.replace("row|r001|add|3|", "row|r001|add|three|"),
            "FS-IVL-CORPUS-BAD-INTEGER",
        ),
        (
            "count-mismatch",
            {
                let end_start = good.rfind("end|").expect("end line");
                format!("{}end|1|deadbeefdeadbeef\n", &good[..end_start])
            },
            "FS-IVL-CORPUS-COUNT-MISMATCH",
        ),
        (
            "root-mismatch",
            {
                let end_start = good.rfind("end|").expect("end line");
                let end_line = good[end_start..].trim_end();
                let mut parts: Vec<&str> = end_line.split('|').collect();
                parts[2] = "0000000000000000";
                format!("{}{}\n", &good[..end_start], parts.join("|"))
            },
            "FS-IVL-CORPUS-ROOT-MISMATCH",
        ),
        (
            "oversize-line",
            good.replace(
                "row|r001|add|",
                &format!("row|r001|add|{}", "9".repeat(300)),
            ),
            "FS-IVL-CORPUS-OVERSIZE",
        ),
    ];
    // Every structural refusal fires at its own line, BEFORE the end-line
    // root check — a corrupt structure must never be laundered into a mere
    // hash mismatch.
    for (name, text, code) in &cases {
        let error = decode(text).expect_err("corrupted corpus must refuse");
        let rendered = error.to_string();
        assert!(
            rendered.contains(code),
            "case {name}: got `{rendered}`, wanted `{code}`"
        );
        println!(
            "{{\"suite\":\"fs-ivl-exact-corpus\",\"case\":\"refusal\",\"name\":\"{name}\",\
             \"verdict\":\"refused\"}}"
        );
    }
}

#[test]
fn g3_mutated_rows_fail_the_battery() {
    // Corrupt semantics, not framing: each mutant decodes fine (its text is
    // never re-serialized), and EXECUTION must catch it.
    let (rows, _) = decode(&built_corpus()).expect("canonical corpus decodes");

    // Numerator corruption: r001 is 3/4 + 1/4 = 1; +1 on the expected
    // numerator claims the sum is 2, which the 2-ULP enclosure excludes.
    let mut numerator = rows[0].clone();
    assert_eq!(numerator.primitive, "add");
    numerator.fields[4] += 1;
    let err = run_row(&numerator).expect_err("numerator corruption caught");
    assert!(err.contains("excludes exact"), "{err}");

    // Endpoint/operand direction corruption: swapping operand a with the
    // expected value claims 1 + 1/4 = 3/4.
    let mut swapped = rows[0].clone();
    swapped.fields.swap(0, 4);
    swapped.fields.swap(1, 5);
    let err = run_row(&swapped).expect_err("direction corruption caught");
    assert!(err.contains("excludes exact"), "{err}");

    // Expected-classification flip on an exact predicate row.
    let mut orient = rows
        .iter()
        .find(|r| r.primitive == "orient2d")
        .expect("orient row present")
        .clone();
    orient.fields[6] = if orient.fields[6] > 0 { -1 } else { 1 };
    let err = run_row(&orient).expect_err("classification flip caught");
    assert!(err.contains("expected"), "{err}");
}
