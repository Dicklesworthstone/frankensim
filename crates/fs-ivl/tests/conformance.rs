//! fs-ivl conformance (plan §13.3): random expression-tree containment
//! against the dd oracle, plus the cross-ISA golden hash over interval
//! endpoints — the same evidence discipline as fs-math/fs-fft/fs-sparse.
//! Any reimplementation must pass this suite; the golden-hash case must
//! match bit-for-bit.

use fs_ivl::{AffineCtx, Interval};
use fs_math::dd::Dd;

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
}

/// Expression node over three leaves; op selection is seed-driven so the
/// battery is deterministic and replayable.
#[derive(Debug, Clone, Copy)]
enum Op {
    Add(usize, usize),
    Sub(usize, usize),
    Mul(usize, usize),
    Exp(usize),
    Tanh(usize),
    Sin(usize),
}

fn random_dag(seed: &mut u64, depth: usize) -> Vec<Op> {
    let mut ops = Vec::with_capacity(depth);
    for level in 0..depth {
        let avail = 3 + level; // leaves 0..3, then one node per level
        let pick = |s: &mut u64, n: usize| -> usize {
            *s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (*s >> 33) as usize % n
        };
        let a = pick(seed, avail);
        let b = pick(seed, avail);
        let op = match pick(seed, 6) {
            0 => Op::Add(a, b),
            1 => Op::Sub(a, b),
            2 => Op::Mul(a, b),
            3 => Op::Exp(a),
            4 => Op::Tanh(a),
            _ => Op::Sin(a),
        };
        ops.push(op);
    }
    ops
}

fn eval_interval(leaves: [Interval; 3], ops: &[Op]) -> Interval {
    let mut vals: Vec<Interval> = leaves.to_vec();
    for op in ops {
        let v = match *op {
            Op::Add(a, b) => vals[a] + vals[b],
            Op::Sub(a, b) => vals[a] - vals[b],
            Op::Mul(a, b) => vals[a] * vals[b],
            Op::Exp(a) => clamp_domain(vals[a]).exp(),
            Op::Tanh(a) => vals[a].tanh(),
            Op::Sin(a) => vals[a].sin(),
        };
        vals.push(v);
    }
    *vals.last().unwrap()
}

/// Keep exp's argument in a sane band so the battery doesn't saturate to
/// ±inf (which is still CORRECT, just uninformative for hashing).
fn clamp_domain(x: Interval) -> Interval {
    Interval::new(x.lo().clamp(-30.0, 30.0), x.hi().clamp(-30.0, 30.0))
}

/// Point oracle: dd for arithmetic, strict f64 (fs-math det) for
/// elementary functions — the interval's declared budgets cover the
/// difference, and containment is checked against BOTH the point value and
/// small dd-scale perturbations of it.
fn eval_point(leaves: [f64; 3], ops: &[Op]) -> f64 {
    let mut vals: Vec<Dd> = leaves.iter().map(|&x| Dd::from_f64(x)).collect();
    for op in ops {
        let v = match *op {
            Op::Add(a, b) => vals[a] + vals[b],
            Op::Sub(a, b) => vals[a] - vals[b],
            Op::Mul(a, b) => vals[a] * vals[b],
            Op::Exp(a) => Dd::from_f64(fs_math::det::exp(vals[a].to_f64().clamp(-30.0, 30.0))),
            Op::Tanh(a) => Dd::from_f64(fs_math::det::tanh(vals[a].to_f64())),
            Op::Sin(a) => Dd::from_f64(fs_math::det::sin(vals[a].to_f64())),
        };
        vals.push(v);
    }
    vals.last().unwrap().to_f64()
}

#[test]
fn expression_tree_containment() {
    // The G0 law over random DAGs: interval evaluation must contain the
    // point evaluation at interior samples. Elementary point values are
    // f64-accurate; interval enclosures are ≥ budget wide, so direct
    // containment with a small ULP grace on the comparison.
    let mut seed = 0x0600_DDA6_u64;
    let mut checked = 0u64;
    for _ in 0..4_000 {
        let ops = random_dag(&mut seed, 6);
        let mut leaves_i = [Interval::point(0.0); 3];
        let mut lo = [0.0f64; 3];
        let mut wid = [0.0f64; 3];
        for k in 0..3 {
            let c = lcg(&mut seed) * 4.0;
            let w = (lcg(&mut seed) + 0.5).abs() + 1e-6;
            leaves_i[k] = Interval::new(c - w, c + w);
            lo[k] = c - w;
            wid[k] = 2.0 * w;
        }
        let enc = eval_interval(leaves_i, &ops);
        for t in [0.0, 0.31, 0.5, 0.77, 1.0] {
            // Clamped into each box: the affine sample can round outside at
            // the ends, and containment only speaks for interior points.
            let pts = [
                (lo[0] + t * wid[0]).clamp(leaves_i[0].lo(), leaves_i[0].hi()),
                (lo[1] + (1.0 - t) * wid[1]).clamp(leaves_i[1].lo(), leaves_i[1].hi()),
                (lo[2] + t * 0.5 * wid[2]).clamp(leaves_i[2].lo(), leaves_i[2].hi()),
            ];
            let p = eval_point(pts, &ops);
            if p.is_finite() {
                assert!(
                    enc.contains(p)
                        || enc.contains(fs_math::next_down(p))
                        || enc.contains(fs_math::next_up(p)),
                    "containment violated: {enc:?} vs point {p} (ops {ops:?})"
                );
                checked += 1;
            }
        }
    }
    println!(
        "{{\"suite\":\"fs-ivl\",\"case\":\"dag-containment\",\"verdict\":\"pass\",\"detail\":\"{checked} point checks over 4000 random DAGs, 0 violations\"}}"
    );
}

/// Recorded on aarch64-apple (M4 Pro); must match on x86-64 (trj) — the
/// cross-ISA determinism evidence for interval endpoints and affine
/// collapse.
const GOLDEN_HASH: u64 = 0x3712_a4c1_2d5e_5864;

#[test]
fn golden_hash_of_interval_endpoints() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for b in v.to_bits().to_le_bytes() {
            acc ^= u64::from(b);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let mut seed = 0x001D_BEEF_u64;
    for _ in 0..500 {
        let ops = random_dag(&mut seed, 8);
        let mut leaves = [Interval::point(0.0); 3];
        for l in &mut leaves {
            let c = lcg(&mut seed) * 3.0;
            let w = (lcg(&mut seed) + 0.5).abs() + 1e-9;
            *l = Interval::new(c - w, c + w);
        }
        let enc = eval_interval(leaves, &ops);
        feed(enc.lo());
        feed(enc.hi());
        // Affine collapse folded into the same hash (deterministic symbol
        // ids make this replayable).
        let mut ctx = AffineCtx::new();
        let x = ctx.from_interval(leaves[0]);
        let y = ctx.from_interval(leaves[1]);
        let aff = (&(&x + &y) * &(&x - &y)).to_interval();
        feed(aff.lo());
        feed(aff.hi());
    }
    println!(
        "{{\"suite\":\"fs-ivl\",\"case\":\"golden-hash\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "interval endpoint bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with \
         semantic justification (golden-evidence policy)"
    );
}

// ---------------------------------------------------------------------------
// G0 property adoption (bead frankensim-4nh8): THE interval law —
// containment under arithmetic — generated + shrunk via fs-propcheck.
// The plan (sec 13.1 G0) names "interval containment under random
// rewrites" explicitly; this covers add/mul composition over generated
// endpoints including the signed-zero/negative corners the specials
// generator draws by design. Fixed cases above remain as pins.
// ---------------------------------------------------------------------------

#[test]
fn g0_interval_arithmetic_contains_the_pointwise_result() {
    fs_propcheck::check(
        "interval-add-mul-containment",
        0x1F_1234,
        600,
        |s| {
            // Two intervals as (center, halfwidth>=0) plus interior picks.
            (
                (s.f64_in(-1e6, 1e6), s.f64_in(0.0, 1e3)),
                (s.f64_in(-1e6, 1e6), s.f64_in(0.0, 1e3)),
            )
        },
        |&((c1, w1), (c2, w2))| {
            let a = fs_ivl::Interval::new(c1 - w1, c1 + w1);
            let b = fs_ivl::Interval::new(c2 - w2, c2 + w2);
            // Pointwise members: endpoints and centers are the adversarial
            // picks (extremes of the exact image under monotone ops).
            let xs = [c1 - w1, c1, c1 + w1];
            let ys = [c2 - w2, c2, c2 + w2];
            xs.iter().all(|&x| {
                ys.iter().all(|&y| {
                    (a + b).contains(x + y) && (a * b).contains(x * y)
                })
            })
        },
    );
    println!("{{\"suite\":\"fs-ivl\",\"case\":\"g0-containment\",\"verdict\":\"pass\",\"detail\":\"600 generated interval pairs, add+mul, shrink-armed\"}}");
}
