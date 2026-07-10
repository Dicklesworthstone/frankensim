//! Exact geometric predicates (Shewchuk 1997; plan §6.4): adaptive-precision
//! `orient2d` / `orient3d` / `incircle` / `insphere` whose SIGNS are exact —
//! "the difference between a mesher that works and a mesher that works
//! until Tuesday."
//!
//! Ladder per predicate: a fast floating-point evaluation guarded by a
//! PROVEN static error bound (Shewchuk's stage-A constants, arrangement-
//! specific); when the sign is uncertain, escalation — `orient2d` runs the
//! full faithful adaptive ladder (stages B/C/D); the others go to the exact
//! stage directly (correctness identical; only near-degenerate throughput
//! differs, and the filter-rate tests show stage A absorbs general-position
//! traffic). The exact stage evaluates the standard difference-based
//! determinant over ≤2-component exact difference expansions with
//! [`crate::expansion`] arithmetic — exact by translation invariance of the
//! determinants.
//!
//! Certified domain (Shewchuk's standard caveat, inherited honestly): signs
//! are exact provided no intermediate overflow/underflow occurs — i.e.
//! coordinate differences and their degree-≤5 monomials stay inside the
//! normal f64 range. Inputs violating that are outside the certificate
//! (CONTRACT.md no-claim).
//!
//! Symbolic perturbation: [`orient2d_sos`] and [`orient3d_sos`] implement the
//! Edelsbrunner–Mücke Simulation-of-Simplicity ladder — 2D derived term-by-term
//! in the code; 3D as the sign of the leading κ-coefficient of a moment-curve
//! perturbation `pᵢ + κ·(sᵢ, sᵢ², sᵢ³)` (consistent because moment-curve points
//! are always in general position), with antisymmetry supplied structurally by
//! sorting the points by index and applying the permutation parity. Both are
//! total (never Zero for index-distinct points), antisymmetric, and
//! deterministic — what BRIO Delaunay needs for reproducible tie-breaking (P2).
//!
//! Determinism: straight-line IEEE arithmetic — cross-ISA bit-deterministic
//! by construction.

use crate::expansion::{
    diff_expansion, estimate, expansion_diff, expansion_product, expansion_sign,
    fast_expansion_sum_zeroelim, prod_diff, two_diff,
};

/// The exact sign of a predicate determinant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sign {
    /// Determinant < 0.
    Negative,
    /// Determinant = 0 exactly (degenerate configuration).
    Zero,
    /// Determinant > 0.
    Positive,
}

impl Sign {
    fn of(x: f64) -> Sign {
        if x > 0.0 {
            Sign::Positive
        } else if x < 0.0 {
            Sign::Negative
        } else {
            Sign::Zero
        }
    }

    fn of_i32(x: i32) -> Sign {
        match x.cmp(&0) {
            std::cmp::Ordering::Greater => Sign::Positive,
            std::cmp::Ordering::Less => Sign::Negative,
            std::cmp::Ordering::Equal => Sign::Zero,
        }
    }

    /// The opposite sign.
    #[must_use]
    pub fn flip(self) -> Sign {
        match self {
            Sign::Negative => Sign::Positive,
            Sign::Zero => Sign::Zero,
            Sign::Positive => Sign::Negative,
        }
    }
}

/// Which rung of the ladder resolved the sign (filter-rate telemetry).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stage {
    /// Plain float + static error bound sufficed.
    Filtered,
    /// An adaptive intermediate stage sufficed (orient2d B/C).
    Adaptive,
    /// The exact expansion stage was required.
    Exact,
}

// Shewchuk's epsilon is the HALF-ulp 2⁻⁵³ (f64::EPSILON is 2⁻⁵²). The
// stage-A/B/C constants are the paper's proven forward-error bounds for
// exactly these evaluation arrangements — do not "simplify" the arithmetic
// below without re-deriving them.
const EPS: f64 = f64::EPSILON / 2.0;
const RESULTERRBOUND: f64 = (3.0 + 8.0 * EPS) * EPS;
const CCWERRBOUND_A: f64 = (3.0 + 16.0 * EPS) * EPS;
const CCWERRBOUND_B: f64 = (2.0 + 12.0 * EPS) * EPS;
const CCWERRBOUND_C: f64 = (9.0 + 64.0 * EPS) * EPS * EPS;
const O3DERRBOUND_A: f64 = (7.0 + 56.0 * EPS) * EPS;
const ICCERRBOUND_A: f64 = (10.0 + 96.0 * EPS) * EPS;
const ISPERRBOUND_A: f64 = (16.0 + 224.0 * EPS) * EPS;

// ---------------------------------------------------------------------------
// orient2d — full faithful adaptive ladder
// ---------------------------------------------------------------------------

/// Exact sign of the 2D orientation determinant
/// `|ax−cx ay−cy; bx−cx by−cy|`: Positive iff `a, b, c` wind
/// counterclockwise.
#[must_use]
pub fn orient2d(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2]) -> Sign {
    orient2d_with_stage(pa, pb, pc).0
}

/// [`orient2d`] plus the ladder stage that resolved it.
#[must_use]
pub fn orient2d_with_stage(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2]) -> (Sign, Stage) {
    let detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
    let detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
    let det = detleft - detright;
    let detsum = if detleft > 0.0 {
        if detright <= 0.0 {
            return (Sign::of(det), Stage::Filtered);
        }
        detleft + detright
    } else if detleft < 0.0 {
        if detright >= 0.0 {
            return (Sign::of(det), Stage::Filtered);
        }
        -detleft - detright
    } else {
        // detleft is exactly 0 ⇒ a factor is exactly 0 (in-range inputs),
        // so det = −detright with detright's float sign exact.
        return (Sign::of(det), Stage::Filtered);
    };
    let errbound = CCWERRBOUND_A * detsum;
    if det >= errbound || -det >= errbound {
        return (Sign::of(det), Stage::Filtered);
    }
    orient2d_adapt(pa, pb, pc, detsum)
}

/// Stages B/C/D of Shewchuk's orient2d, ported faithfully.
fn orient2d_adapt(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2], detsum: f64) -> (Sign, Stage) {
    let acx = pa[0] - pc[0];
    let bcx = pb[0] - pc[0];
    let acy = pa[1] - pc[1];
    let bcy = pb[1] - pc[1];

    // Stage B: exact determinant of the ROUNDED differences.
    let b_exp = prod_diff(acx, bcy, acy, bcx);
    let mut det = estimate(&b_exp);
    let errbound_b = CCWERRBOUND_B * detsum;
    if det >= errbound_b || -det >= errbound_b {
        return (Sign::of(det), Stage::Adaptive);
    }

    // Tails of the four differences (exact).
    let acxtail = two_diff(pa[0], pc[0]).1;
    let bcxtail = two_diff(pb[0], pc[0]).1;
    let acytail = two_diff(pa[1], pc[1]).1;
    let bcytail = two_diff(pb[1], pc[1]).1;
    if acxtail == 0.0 && acytail == 0.0 && bcxtail == 0.0 && bcytail == 0.0 {
        // The differences were exact: stage B already was the exact answer.
        return (Sign::of(det), Stage::Adaptive);
    }

    // Stage C: first-order tail correction with its proven bound.
    let errbound_c = CCWERRBOUND_C * detsum + RESULTERRBOUND * det.abs();
    det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
    if det >= errbound_c || -det >= errbound_c {
        return (Sign::of(det), Stage::Adaptive);
    }

    // Stage D: fully exact — fold in every tail cross term.
    let u1 = prod_diff(acxtail, bcy, acytail, bcx);
    let c1 = fast_expansion_sum_zeroelim(&b_exp, &u1);
    let u2 = prod_diff(acx, bcytail, acy, bcxtail);
    let c2 = fast_expansion_sum_zeroelim(&c1, &u2);
    let u3 = prod_diff(acxtail, bcytail, acytail, bcxtail);
    let d = fast_expansion_sum_zeroelim(&c2, &u3);
    (Sign::of_i32(expansion_sign(&d)), Stage::Exact)
}

// ---------------------------------------------------------------------------
// orient3d
// ---------------------------------------------------------------------------

/// Exact sign of the 3D orientation determinant: Positive iff `pd` lies
/// below the plane through `pa, pb, pc` oriented so they appear
/// counterclockwise from above (Shewchuk's convention).
#[must_use]
pub fn orient3d(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3]) -> Sign {
    orient3d_with_stage(pa, pb, pc, pd).0
}

/// [`orient3d`] plus the resolving stage.
#[must_use]
pub fn orient3d_with_stage(
    pa: [f64; 3],
    pb: [f64; 3],
    pc: [f64; 3],
    pd: [f64; 3],
) -> (Sign, Stage) {
    let adx = pa[0] - pd[0];
    let bdx = pb[0] - pd[0];
    let cdx = pc[0] - pd[0];
    let ady = pa[1] - pd[1];
    let bdy = pb[1] - pd[1];
    let cdy = pc[1] - pd[1];
    let adz = pa[2] - pd[2];
    let bdz = pb[2] - pd[2];
    let cdz = pc[2] - pd[2];

    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;
    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;
    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;

    let det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);
    let permanent = (bdxcdy.abs() + cdxbdy.abs()) * adz.abs()
        + (cdxady.abs() + adxcdy.abs()) * bdz.abs()
        + (adxbdy.abs() + bdxady.abs()) * cdz.abs();
    let errbound = O3DERRBOUND_A * permanent;
    if det > errbound || -det > errbound {
        return (Sign::of(det), Stage::Filtered);
    }

    // Exact: the same determinant over exact ≤2-component differences.
    let ax = diff_expansion(pa[0], pd[0]);
    let bx = diff_expansion(pb[0], pd[0]);
    let cx = diff_expansion(pc[0], pd[0]);
    let ay = diff_expansion(pa[1], pd[1]);
    let by = diff_expansion(pb[1], pd[1]);
    let cy = diff_expansion(pc[1], pd[1]);
    let az = diff_expansion(pa[2], pd[2]);
    let bz = diff_expansion(pb[2], pd[2]);
    let cz = diff_expansion(pc[2], pd[2]);

    let m_a = expansion_diff(&expansion_product(&bx, &cy), &expansion_product(&cx, &by));
    let m_b = expansion_diff(&expansion_product(&cx, &ay), &expansion_product(&ax, &cy));
    let m_c = expansion_diff(&expansion_product(&ax, &by), &expansion_product(&bx, &ay));
    let t =
        fast_expansion_sum_zeroelim(&expansion_product(&az, &m_a), &expansion_product(&bz, &m_b));
    let det_e = fast_expansion_sum_zeroelim(&t, &expansion_product(&cz, &m_c));
    (Sign::of_i32(expansion_sign(&det_e)), Stage::Exact)
}

// ---------------------------------------------------------------------------
// incircle
// ---------------------------------------------------------------------------

/// Exact sign of the incircle determinant: Positive iff `pd` lies inside
/// the circle through `pa, pb, pc` (which must wind counterclockwise;
/// clockwise inputs flip the sign — Shewchuk's convention).
#[must_use]
pub fn incircle(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2], pd: [f64; 2]) -> Sign {
    incircle_with_stage(pa, pb, pc, pd).0
}

/// [`incircle`] plus the resolving stage.
#[must_use]
pub fn incircle_with_stage(
    pa: [f64; 2],
    pb: [f64; 2],
    pc: [f64; 2],
    pd: [f64; 2],
) -> (Sign, Stage) {
    let adx = pa[0] - pd[0];
    let bdx = pb[0] - pd[0];
    let cdx = pc[0] - pd[0];
    let ady = pa[1] - pd[1];
    let bdy = pb[1] - pd[1];
    let cdy = pc[1] - pd[1];

    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;
    let alift = adx * adx + ady * ady;
    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;
    let blift = bdx * bdx + bdy * bdy;
    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;
    let clift = cdx * cdx + cdy * cdy;

    let det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);
    let permanent = (bdxcdy.abs() + cdxbdy.abs()) * alift
        + (cdxady.abs() + adxcdy.abs()) * blift
        + (adxbdy.abs() + bdxady.abs()) * clift;
    let errbound = ICCERRBOUND_A * permanent;
    if det > errbound || -det > errbound {
        return (Sign::of(det), Stage::Filtered);
    }

    // Exact over exact differences.
    let ax = diff_expansion(pa[0], pd[0]);
    let bx = diff_expansion(pb[0], pd[0]);
    let cx = diff_expansion(pc[0], pd[0]);
    let ay = diff_expansion(pa[1], pd[1]);
    let by = diff_expansion(pb[1], pd[1]);
    let cy = diff_expansion(pc[1], pd[1]);

    let lift = |x: &[f64], y: &[f64]| -> Vec<f64> {
        fast_expansion_sum_zeroelim(&expansion_product(x, x), &expansion_product(y, y))
    };
    let m_bc = expansion_diff(&expansion_product(&bx, &cy), &expansion_product(&cx, &by));
    let m_ca = expansion_diff(&expansion_product(&cx, &ay), &expansion_product(&ax, &cy));
    let m_ab = expansion_diff(&expansion_product(&ax, &by), &expansion_product(&bx, &ay));
    let t = fast_expansion_sum_zeroelim(
        &expansion_product(&lift(&ax, &ay), &m_bc),
        &expansion_product(&lift(&bx, &by), &m_ca),
    );
    let det_e = fast_expansion_sum_zeroelim(&t, &expansion_product(&lift(&cx, &cy), &m_ab));
    (Sign::of_i32(expansion_sign(&det_e)), Stage::Exact)
}

// ---------------------------------------------------------------------------
// insphere
// ---------------------------------------------------------------------------

/// Exact sign of the insphere determinant: Positive iff `pe` lies inside
/// the sphere through `pa..pd` (which must be positively oriented per
/// [`orient3d`]; negatively oriented inputs flip the sign).
#[must_use]
pub fn insphere(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3], pe: [f64; 3]) -> Sign {
    insphere_with_stage(pa, pb, pc, pd, pe).0
}

/// [`insphere`] plus the resolving stage.
#[must_use]
pub fn insphere_with_stage(
    pa: [f64; 3],
    pb: [f64; 3],
    pc: [f64; 3],
    pd: [f64; 3],
    pe: [f64; 3],
) -> (Sign, Stage) {
    let aex = pa[0] - pe[0];
    let bex = pb[0] - pe[0];
    let cex = pc[0] - pe[0];
    let dex = pd[0] - pe[0];
    let aey = pa[1] - pe[1];
    let bey = pb[1] - pe[1];
    let cey = pc[1] - pe[1];
    let dey = pd[1] - pe[1];
    let aez = pa[2] - pe[2];
    let bez = pb[2] - pe[2];
    let cez = pc[2] - pe[2];
    let dez = pd[2] - pe[2];

    let ab = aex * bey - bex * aey;
    let bc = bex * cey - cex * bey;
    let cd = cex * dey - dex * cey;
    let da = dex * aey - aex * dey;
    let ac = aex * cey - cex * aey;
    let bd = bex * dey - dex * bey;

    let abc = aez * bc - bez * ac + cez * ab;
    let bcd = bez * cd - cez * bd + dez * bc;
    let cda = cez * da + dez * ac + aez * cd;
    let dab = dez * ab + aez * bd + bez * da;

    let alift = aex * aex + aey * aey + aez * aez;
    let blift = bex * bex + bey * bey + bez * bez;
    let clift = cex * cex + cey * cey + cez * cez;
    let dlift = dex * dex + dey * dey + dez * dez;

    let det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

    let ab_p = (aex * bey).abs() + (bex * aey).abs();
    let bc_p = (bex * cey).abs() + (cex * bey).abs();
    let cd_p = (cex * dey).abs() + (dex * cey).abs();
    let da_p = (dex * aey).abs() + (aex * dey).abs();
    let ac_p = (aex * cey).abs() + (cex * aey).abs();
    let bd_p = (bex * dey).abs() + (dex * bey).abs();
    let abc_p = aez.abs() * bc_p + bez.abs() * ac_p + cez.abs() * ab_p;
    let bcd_p = bez.abs() * cd_p + cez.abs() * bd_p + dez.abs() * bc_p;
    let cda_p = cez.abs() * da_p + dez.abs() * ac_p + aez.abs() * cd_p;
    let dab_p = dez.abs() * ab_p + aez.abs() * bd_p + bez.abs() * da_p;
    let permanent = dlift * abc_p + clift * dab_p + blift * cda_p + alift * bcd_p;
    let errbound = ISPERRBOUND_A * permanent;
    if det > errbound || -det > errbound {
        return (Sign::of(det), Stage::Filtered);
    }

    insphere_exact(pa, pb, pc, pd, pe)
}

fn insphere_exact(
    pa: [f64; 3],
    pb: [f64; 3],
    pc: [f64; 3],
    pd: [f64; 3],
    pe: [f64; 3],
) -> (Sign, Stage) {
    let d = |p: [f64; 3], i: usize| diff_expansion(p[i], pe[i]);
    let (ax, ay, az) = (d(pa, 0), d(pa, 1), d(pa, 2));
    let (bx, by, bz) = (d(pb, 0), d(pb, 1), d(pb, 2));
    let (cx, cy, cz) = (d(pc, 0), d(pc, 1), d(pc, 2));
    let (dx, dy, dz) = (d(pd, 0), d(pd, 1), d(pd, 2));

    let minor2 = |px: &[f64], py: &[f64], qx: &[f64], qy: &[f64]| -> Vec<f64> {
        expansion_diff(&expansion_product(px, qy), &expansion_product(qx, py))
    };
    let ab = minor2(&ax, &ay, &bx, &by);
    let bc = minor2(&bx, &by, &cx, &cy);
    let cd = minor2(&cx, &cy, &dx, &dy);
    let da = minor2(&dx, &dy, &ax, &ay);
    let ac = minor2(&ax, &ay, &cx, &cy);
    let bd = minor2(&bx, &by, &dx, &dy);

    let sum3 = |a: &[f64], b: &[f64], c: &[f64]| -> Vec<f64> {
        fast_expansion_sum_zeroelim(&fast_expansion_sum_zeroelim(a, b), c)
    };
    let neg = |e: &[f64]| -> Vec<f64> { e.iter().map(|&x| -x).collect() };

    let abc = sum3(
        &expansion_product(&az, &bc),
        &neg(&expansion_product(&bz, &ac)),
        &expansion_product(&cz, &ab),
    );
    let bcd = sum3(
        &expansion_product(&bz, &cd),
        &neg(&expansion_product(&cz, &bd)),
        &expansion_product(&dz, &bc),
    );
    let cda = sum3(
        &expansion_product(&cz, &da),
        &expansion_product(&dz, &ac),
        &expansion_product(&az, &cd),
    );
    let dab = sum3(
        &expansion_product(&dz, &ab),
        &expansion_product(&az, &bd),
        &expansion_product(&bz, &da),
    );

    let lift = |x: &[f64], y: &[f64], z: &[f64]| -> Vec<f64> {
        sum3(
            &expansion_product(x, x),
            &expansion_product(y, y),
            &expansion_product(z, z),
        )
    };
    let alift = lift(&ax, &ay, &az);
    let blift = lift(&bx, &by, &bz);
    let clift = lift(&cx, &cy, &cz);
    let dlift = lift(&dx, &dy, &dz);

    let det_e = fast_expansion_sum_zeroelim(
        &expansion_diff(
            &expansion_product(&dlift, &abc),
            &expansion_product(&clift, &dab),
        ),
        &expansion_diff(
            &expansion_product(&blift, &cda),
            &expansion_product(&alift, &bcd),
        ),
    );
    (Sign::of_i32(expansion_sign(&det_e)), Stage::Exact)
}

// ---------------------------------------------------------------------------
// Symbolic perturbation (SoS) hooks
// ---------------------------------------------------------------------------

/// Sort three rows by index ascending; returns permutation parity
/// (+1 even, −1 odd).
fn sort3_by_index<T: Copy>(mut rows: [(u64, T); 3]) -> ([(u64, T); 3], i32) {
    let mut parity = 1i32;
    if rows[0].0 > rows[1].0 {
        rows.swap(0, 1);
        parity = -parity;
    }
    if rows[1].0 > rows[2].0 {
        rows.swap(1, 2);
        parity = -parity;
    }
    if rows[0].0 > rows[1].0 {
        rows.swap(0, 1);
        parity = -parity;
    }
    (rows, parity)
}

/// `orient2d` with Edelsbrunner–Mücke Simulation-of-Simplicity
/// tie-breaking: NEVER returns [`Sign::Zero`] for points with distinct
/// indices, is antisymmetric under argument swaps, and resolves identical
/// ties identically (deterministic Delaunay, P2).
///
/// Derivation (indices sorted ascending as rows 1..3, perturbation
/// `ε(i,j) = δ^(2^(2i−j))` so lower index ⇒ larger perturbation): the
/// perturbed determinant's ε-terms in increasing exponent are
/// `D; (x₃−x₂)·ε(1,2); (y₂−y₃)·ε(1,1); (x₁−x₃)·ε(2,2); +1·ε(1,1)ε(2,2)`.
/// The constant final term makes the ladder total. Each coefficient's sign
/// is evaluated EXACTLY (difference expansions), so the hook inherits the
/// predicate's certainty.
#[must_use]
pub fn orient2d_sos(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2], ia: u64, ib: u64, ic: u64) -> Sign {
    debug_assert!(
        ia != ib && ib != ic && ia != ic,
        "SoS requires distinct indices"
    );
    let base = orient2d(pa, pb, pc);
    if base != Sign::Zero {
        return base;
    }
    let ([(_, p1), (_, p2), (_, p3)], parity) = sort3_by_index([(ia, pa), (ib, pb), (ic, pc)]);
    let apply = |s: Sign| if parity > 0 { s } else { s.flip() };
    let terms = [
        diff_expansion(p3[0], p2[0]), // x3 − x2
        diff_expansion(p2[1], p3[1]), // y2 − y3
        diff_expansion(p1[0], p3[0]), // x1 − x3
    ];
    for t in &terms {
        let s = expansion_sign(t);
        if s != 0 {
            return apply(Sign::of_i32(s));
        }
    }
    apply(Sign::Positive) // the +1·ε(1,1)ε(2,2) constant term
}

/// `orient3d` with Edelsbrunner–Mücke Simulation-of-Simplicity tie-breaking:
/// NEVER returns [`Sign::Zero`] for points with distinct indices, and is
/// ANTISYMMETRIC under argument swaps — reproducible Delaunay tie-breaking (P2).
///
/// On an exactly coplanar input the sign is the limit, as `κ → 0⁺`, of
/// `orient3d` of the points moved along a MOMENT CURVE:
/// `pᵢ ↦ pᵢ + κ·(sᵢ, sᵢ², sᵢ³)` with distinct `sᵢ`. Moment-curve points are
/// always in general position, so this perturbation is CONSISTENT (realizable)
/// — the SoS answer is a genuine orientation, not an ad-hoc tie. Antisymmetry
/// is STRUCTURAL: the four points are sorted by index (parity tracked), the
/// perturbation is keyed to sorted RANK so the coefficient ladder is a pure
/// function of the point SET, and the permutation parity supplies the swap
/// sign. (The previous projection heuristic ignored the 4th point's position
/// in its exact branch and so was NOT antisymmetric — bead wa8i V1.)
///
/// `D(κ) = det(pᵢ + κ·vᵢ, 1) = Σₘ κᵐ Cₘ` with `C₀ = 0` (coplanar) and `C₄ = 0`
/// (the perturbation rows' constant column is 0); the sign of the first nonzero
/// `Cₘ` (`m = 1, 2, 3`) decides. Expanding each `det` along the constant column
/// reduces every `Cₘ` to a signed sum of EXACT 3×3 determinants.
#[must_use]
pub fn orient3d_sos(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3], idx: [u64; 4]) -> Sign {
    let [ia, ib, ic, id] = idx;
    debug_assert!(
        ia != ib && ia != ic && ia != id && ib != ic && ib != id && ic != id,
        "SoS requires distinct indices"
    );
    let base = orient3d(pa, pb, pc, pd);
    if base != Sign::Zero {
        return base;
    }
    // Coplanar: sort by index (parity) so the tie-break is order-independent;
    // parity then supplies antisymmetry under any argument swap.
    let (rows, parity) = sort4_by_index([(ia, pa), (ib, pb), (ic, pc), (id, pd)]);
    let q = [rows[0].1, rows[1].1, rows[2].1, rows[3].1];
    // Moment-curve directions keyed to SORTED rank (distinct sᵢ ⇒ general
    // position ⇒ consistent). sᵢ = i+1 keeps every direction nonzero.
    let v: [[f64; 3]; 4] = core::array::from_fn(|j| {
        let s = (j + 1) as f64;
        [s, s * s, s * s * s]
    });
    let apply = |s: i32| {
        let sign = Sign::of_i32(s);
        if parity > 0 { sign } else { sign.flip() }
    };
    for m in 1..=3 {
        let s = moment_coeff_sign(&q, &v, m);
        if s != 0 {
            return apply(s);
        }
    }
    // Distinct indices ⇒ some Cₘ is nonzero (moment-curve genericity); this
    // only guards an impossible all-degenerate residue.
    apply(1)
}

/// Sort four `(index, row)` pairs by index ascending; returns the sorted rows
/// and the permutation parity (+1 even, −1 odd).
fn sort4_by_index<T: Copy>(mut rows: [(u64, T); 4]) -> ([(u64, T); 4], i32) {
    let mut parity = 1i32;
    for i in 0..4 {
        for j in 0..3 - i {
            if rows[j].0 > rows[j + 1].0 {
                rows.swap(j, j + 1);
                parity = -parity;
            }
        }
    }
    (rows, parity)
}

/// The sign of the `κᵐ` coefficient of `det(qᵢ + κ·vᵢ, 1)`, computed EXACTLY.
/// Expanding `det(M_S)` — the rows in a size-`m` subset `S` replaced by their
/// perturbation rows `vᵢ` (constant column 0), the rest kept as `(qᵢ, 1)` —
/// along the constant column leaves only the `(qᵢ, 1)` rows, so
/// `det(M_S) = Σ_{r∉S} (−1)^{r+1} · det3(rows ≠ r)` (row `j` is `v[j]` if
/// `j ∈ S` else `q[j]`). `Cₘ` sums that over every size-`m` subset.
fn moment_coeff_sign(q: &[[f64; 3]; 4], v: &[[f64; 3]; 4], m: usize) -> i32 {
    let mut acc: Vec<f64> = Vec::new();
    for mask in 0u32..16 {
        if mask.count_ones() as usize != m {
            continue;
        }
        for r in 0..4 {
            if mask & (1 << r) != 0 {
                continue; // r must be OUTSIDE S (its (q, 1) row survives the expansion)
            }
            let mut three = [[0.0f64; 3]; 3];
            let mut t = 0;
            for j in 0..4 {
                if j == r {
                    continue;
                }
                three[t] = if mask & (1 << j) != 0 { v[j] } else { q[j] };
                t += 1;
            }
            let mut d3 = exact_det3(three[0], three[1], three[2]);
            if r % 2 == 0 {
                // (−1)^{r+1} = −1 for even r.
                for x in &mut d3 {
                    *x = -*x;
                }
            }
            acc = fast_expansion_sum_zeroelim(&acc, &d3);
        }
    }
    expansion_sign(&acc)
}

/// The exact 3×3 determinant of three row vectors, as an expansion:
/// `a0(b1c2 − b2c1) − a1(b0c2 − b2c0) + a2(b0c1 − b1c0)`.
fn exact_det3(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Vec<f64> {
    let m0 = prod_diff(b[1], c[2], b[2], c[1]);
    let m1 = prod_diff(b[0], c[2], b[2], c[0]);
    let m2 = prod_diff(b[0], c[1], b[1], c[0]);
    let t0 = expansion_product(&[a[0]], &m0);
    let t1 = expansion_product(&[a[1]], &m1);
    let t2 = expansion_product(&[a[2]], &m2);
    fast_expansion_sum_zeroelim(&expansion_diff(&t0, &t1), &t2)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- i128 lattice oracle -------------------------------------------

    fn o2d_i128(a: [i64; 2], b: [i64; 2], c: [i64; 2]) -> i128 {
        let (ax, ay) = (i128::from(a[0] - c[0]), i128::from(a[1] - c[1]));
        let (bx, by) = (i128::from(b[0] - c[0]), i128::from(b[1] - c[1]));
        ax * by - ay * bx
    }

    fn o3d_i128(a: [i64; 3], b: [i64; 3], c: [i64; 3], d: [i64; 3]) -> i128 {
        let v = |p: [i64; 3]| {
            [
                i128::from(p[0] - d[0]),
                i128::from(p[1] - d[1]),
                i128::from(p[2] - d[2]),
            ]
        };
        let (a, b, c) = (v(a), v(b), v(c));
        a[2] * (b[0] * c[1] - c[0] * b[1])
            + b[2] * (c[0] * a[1] - a[0] * c[1])
            + c[2] * (a[0] * b[1] - b[0] * a[1])
    }

    fn icc_i128(a: [i64; 2], b: [i64; 2], c: [i64; 2], d: [i64; 2]) -> i128 {
        let v = |p: [i64; 2]| {
            let x = i128::from(p[0] - d[0]);
            let y = i128::from(p[1] - d[1]);
            (x, y, x * x + y * y)
        };
        let (a, b, c) = (v(a), v(b), v(c));
        a.2 * (b.0 * c.1 - c.0 * b.1)
            + b.2 * (c.0 * a.1 - a.0 * c.1)
            + c.2 * (a.0 * b.1 - b.0 * a.1)
    }

    fn isp_i128(a: [i64; 3], b: [i64; 3], c: [i64; 3], d: [i64; 3], e: [i64; 3]) -> i128 {
        let v = |p: [i64; 3]| {
            let x = i128::from(p[0] - e[0]);
            let y = i128::from(p[1] - e[1]);
            let z = i128::from(p[2] - e[2]);
            (x, y, z, x * x + y * y + z * z)
        };
        let (a, b, c, d) = (v(a), v(b), v(c), v(d));
        let m2 = |p: (i128, i128, i128, i128), q: (i128, i128, i128, i128)| p.0 * q.1 - q.0 * p.1;
        let (ab, bc, cd, da, ac, bd) = (m2(a, b), m2(b, c), m2(c, d), m2(d, a), m2(a, c), m2(b, d));
        let abc = a.2 * bc - b.2 * ac + c.2 * ab;
        let bcd = b.2 * cd - c.2 * bd + d.2 * bc;
        let cda = c.2 * da + d.2 * ac + a.2 * cd;
        let dab = d.2 * ab + a.2 * bd + b.2 * da;
        (d.3 * abc - c.3 * dab) + (b.3 * cda - a.3 * bcd)
    }

    fn sgn(x: i128) -> Sign {
        match x.cmp(&0) {
            std::cmp::Ordering::Greater => Sign::Positive,
            std::cmp::Ordering::Less => Sign::Negative,
            std::cmp::Ordering::Equal => Sign::Zero,
        }
    }

    fn f2(p: [i64; 2]) -> [f64; 2] {
        [p[0] as f64, p[1] as f64]
    }

    fn f3(p: [i64; 3]) -> [f64; 3] {
        [p[0] as f64, p[1] as f64, p[2] as f64]
    }

    fn lcg(seed: &mut u64) -> i64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 24).cast_signed() % 41) - 20 // small lattice band forces degeneracy
    }

    #[test]
    fn lattice_battery_matches_i128_oracle() {
        let mut seed = 0x5EED_09E0_0000_0014u64;
        let mut zeros = 0usize;
        for round in 0..4000 {
            let p = |s: &mut u64| [lcg(s), lcg(s)];
            let (a, b, c, d) = (p(&mut seed), p(&mut seed), p(&mut seed), p(&mut seed));
            assert_eq!(
                orient2d(f2(a), f2(b), f2(c)),
                sgn(o2d_i128(a, b, c)),
                "orient2d vs oracle, round {round}: {a:?} {b:?} {c:?}"
            );
            let s = incircle(f2(a), f2(b), f2(c), f2(d));
            assert_eq!(s, sgn(icc_i128(a, b, c, d)), "incircle round {round}");
            if s == Sign::Zero {
                zeros += 1;
            }
        }
        assert!(zeros > 0, "battery too easy: no exact degeneracies hit");
    }

    #[test]
    fn lattice_battery_3d_matches_i128_oracle() {
        let mut seed = 0x5EED_09E0_0000_003Du64;
        for round in 0..2000 {
            let p = |s: &mut u64| [lcg(s), lcg(s), lcg(s)];
            let (a, b, c, d, e) = (
                p(&mut seed),
                p(&mut seed),
                p(&mut seed),
                p(&mut seed),
                p(&mut seed),
            );
            assert_eq!(
                orient3d(f3(a), f3(b), f3(c), f3(d)),
                sgn(o3d_i128(a, b, c, d)),
                "orient3d vs oracle, round {round}"
            );
            assert_eq!(
                insphere(f3(a), f3(b), f3(c), f3(d), f3(e)),
                sgn(isp_i128(a, b, c, d, e)),
                "insphere vs oracle, round {round}"
            );
        }
    }

    #[test]
    fn exact_degeneracies_are_zero() {
        // Collinear grid points.
        assert_eq!(orient2d([0.0, 0.0], [1.0, 1.0], [2.0, 2.0]), Sign::Zero);
        // Cocircular lattice points on x² + y² = 25.
        let (a, b, c, d) = ([5.0, 0.0], [3.0, 4.0], [-3.0, 4.0], [0.0, -5.0]);
        assert_eq!(incircle(a, b, c, d), Sign::Zero);
        // Coplanar 3D.
        assert_eq!(
            orient3d([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]),
            Sign::Zero
        );
        // Cospherical lattice points on x² + y² + z² = 9.
        let s = [
            [3.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 3.0],
            [-3.0, 0.0, 0.0],
            [0.0, -3.0, 0.0],
        ];
        assert_eq!(insphere(s[0], s[1], s[2], s[3], s[4]), Sign::Zero);
    }

    #[test]
    fn one_ulp_perturbations_resolve_with_known_truth() {
        // Exactly representable dyadic scaling: coords ≤ 2¹⁰, factor
        // 1 + 2⁻⁴⁰ ⇒ products have ≤ 51 mantissa bits (EXACT), so the
        // scaled point is exactly OUTSIDE its sphere/circle.
        let grow = 1.0 + (2.0f64).powi(-40); // det-ok: base 2, exact (4xnt)
        let (a, b, c) = ([5.0, 0.0], [3.0, 4.0], [-3.0, 4.0]);
        let d_on = [0.0, -5.0];
        let d_out = [0.0, -5.0 * grow];
        let d_in = [0.0, -5.0 / 1.000_000_1];
        assert_eq!(incircle(a, b, c, d_on), Sign::Zero);
        assert_eq!(incircle(a, b, c, d_out), Sign::Negative, "outside");
        assert_eq!(incircle(a, b, c, d_in), Sign::Positive, "inside");
        // 3D: positively oriented base, then on/out/in.
        let (pa, pb, pc, pd) = (
            [3.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 3.0],
            [-3.0, 0.0, 0.0],
        );
        assert_eq!(
            orient3d(pa, pb, pc, pd),
            Sign::Positive,
            "base orientation premise"
        );
        assert_eq!(insphere(pa, pb, pc, pd, [0.0, -3.0, 0.0]), Sign::Zero);
        assert_eq!(
            insphere(pa, pb, pc, pd, [0.0, -3.0 * grow, 0.0]),
            Sign::Negative
        );
        assert_eq!(
            insphere(pa, pb, pc, pd, [0.0, -2.999_999, 0.0]),
            Sign::Positive
        );
    }

    #[test]
    fn kettner_class_breakdown_naive_float_disagrees_somewhere() {
        // Sweep a near-collinear family; the NAIVE float determinant must
        // misjudge at least one case that the exact predicate gets right
        // (the Kettner-class demonstration), while the exact predicate
        // matches the analytically known sign everywhere.
        // Kettner et al.'s classic grid: p walks a half-ulp lattice around
        // (0.5, 0.5) with INDEPENDENT x/y offsets; q = (12,12), r = (24,24)
        // sit far along the diagonal. The naive float determinant misjudges
        // many cells; the exact predicate must stay internally consistent
        // (antisymmetric) everywhere and match the i128 truth of the
        // scaled-integer reformulation.
        let ulp = f64::EPSILON * 0.5; // exact steps: ulp(0.5) = 2⁻⁵³
        let q = [12.0, 12.0];
        let r = [24.0, 24.0];
        let mut naive_wrong = 0usize;
        for i in 0..17u32 {
            for k in 0..17u32 {
                let p = [0.5 + f64::from(i) * ulp, 0.5 + f64::from(k) * ulp];
                let exact = orient2d(p, q, r);
                // Ground truth via exact integer reformulation: scale by
                // 2⁵⁴ (all values integral: 0.5·2⁵⁴ + i·2, 12·2⁵⁴, 24·2⁵⁴).
                // det-ok: base 2, exact (4xnt)
                let s = |x: f64| (x * (2.0f64).powi(54)) as i128;
                let truth = ((s(p[0]) - s(r[0])) * (s(q[1]) - s(r[1]))
                    - (s(p[1]) - s(r[1])) * (s(q[0]) - s(r[0])))
                .signum();
                assert_eq!(
                    exact,
                    Sign::of_i32(truth as i32),
                    "exact predicate wrong at i={i}, k={k}"
                );
                // Internal consistency: antisymmetry under swaps.
                assert_eq!(orient2d(q, p, r), exact.flip(), "antisymmetry i={i} k={k}");
                let naive = Sign::of((p[0] - r[0]) * (q[1] - r[1]) - (p[1] - r[1]) * (q[0] - r[0]));
                if naive != exact {
                    naive_wrong += 1;
                }
            }
        }
        assert!(
            naive_wrong > 0,
            "sweep failed to produce a naive-float breakdown — weak battery"
        );
    }

    #[test]
    fn sos_is_total_antisymmetric_and_consistent() {
        // Degenerate family: collinear, coincident, and mixed points.
        let pts: [[f64; 2]; 4] = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [1.0, 1.0]];
        for (i, &a) in pts.iter().enumerate() {
            for (j, &b) in pts.iter().enumerate() {
                for (k, &c) in pts.iter().enumerate() {
                    if i == j || j == k || i == k {
                        continue;
                    }
                    let (ia, ib, ic) = (i as u64, j as u64, k as u64);
                    let s = orient2d_sos(a, b, c, ia, ib, ic);
                    assert_ne!(s, Sign::Zero, "SoS must be total ({i},{j},{k})");
                    // Antisymmetry under one swap.
                    assert_eq!(
                        orient2d_sos(b, a, c, ib, ia, ic),
                        s.flip(),
                        "SoS antisymmetry ({i},{j},{k})"
                    );
                    // Determinism.
                    assert_eq!(orient2d_sos(a, b, c, ia, ib, ic), s);
                }
            }
        }
        // Agreement with the plain predicate when nonzero.
        let s = orient2d_sos([0.0, 0.0], [1.0, 0.0], [0.0, 1.0], 0, 1, 2);
        assert_eq!(s, Sign::Positive);
        // 3D cascade: total + antisymmetric on a coplanar quadruple.
        let q: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        let s3 = orient3d_sos(q[0], q[1], q[2], q[3], [0, 1, 2, 3]);
        assert_ne!(s3, Sign::Zero);
        assert_eq!(
            orient3d_sos(q[1], q[0], q[2], q[3], [1, 0, 2, 3]),
            s3.flip()
        );
    }

    #[test]
    fn orient3d_sos_full_em_ladder_antisymmetric_total_consistent() {
        // Coplanar quadruples spanning generic, tilted, collinear, and
        // coincident degeneracies (indices always distinct).
        let configs: [[[f64; 3]; 4]; 5] = [
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
            ], // z = 0
            [
                [0.0, 0.0, 0.0],
                [4.0, 4.0, 0.0],
                [4.0, 0.0, 4.0],
                [2.0, 1.0, 1.0],
            ], // x = y+z (V1)
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 1.0],
                [0.0, 1.0, 2.0],
                [1.0, 1.0, 3.0],
            ], // z = x+2y
            [
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
                [2.0, 2.0, 2.0],
                [3.0, 3.0, 3.0],
            ], // collinear
            [
                [1.0, 1.0, 1.0],
                [1.0, 1.0, 1.0],
                [2.0, 3.0, 4.0],
                [5.0, 6.0, 7.0],
            ], // coincident
        ];
        let parity = |p: [usize; 4]| -> i32 {
            let mut inv = 0;
            for i in 0..4 {
                for j in (i + 1)..4 {
                    if p[i] > p[j] {
                        inv += 1;
                    }
                }
            }
            if inv % 2 == 0 { 1 } else { -1 }
        };
        for cfg in &configs {
            // Must actually exercise SoS: the base predicate is exactly Zero.
            assert_eq!(
                orient3d(cfg[0], cfg[1], cfg[2], cfg[3]),
                Sign::Zero,
                "config is not coplanar"
            );
            let s0 = orient3d_sos(cfg[0], cfg[1], cfg[2], cfg[3], [0, 1, 2, 3]);
            assert_ne!(s0, Sign::Zero, "SoS must be total");
            // Antisymmetry + totality + determinism over ALL 24 permutations.
            for a in 0..4 {
                for b in 0..4 {
                    for c in 0..4 {
                        for d in 0..4 {
                            let p = [a, b, c, d];
                            if a == b || a == c || a == d || b == c || b == d || c == d {
                                continue;
                            }
                            let idx = [a as u64, b as u64, c as u64, d as u64];
                            let s = orient3d_sos(cfg[a], cfg[b], cfg[c], cfg[d], idx);
                            assert_ne!(s, Sign::Zero, "SoS total under perm {p:?}");
                            let want = if parity(p) > 0 { s0 } else { s0.flip() };
                            assert_eq!(s, want, "antisymmetry under perm {p:?}");
                            // Determinism.
                            assert_eq!(orient3d_sos(cfg[a], cfg[b], cfg[c], cfg[d], idx), s);
                        }
                    }
                }
            }
            // CONSISTENCY: the SoS sign must equal orient3d of the points moved
            // along the SAME moment curve for a small concrete κ — proof the
            // exact ladder computes the true perturbation limit, not an ad-hoc
            // tie. (Canonical order ⇒ sorted rank = position, parity = +1.)
            for &kappa in &[1e-3_f64, 1e-4, 1e-5] {
                let pert: [[f64; 3]; 4] = core::array::from_fn(|i| {
                    let s = (i + 1) as f64;
                    [
                        kappa.mul_add(s, cfg[i][0]),
                        kappa.mul_add(s * s, cfg[i][1]),
                        kappa.mul_add(s * s * s, cfg[i][2]),
                    ]
                });
                assert_eq!(
                    orient3d(pert[0], pert[1], pert[2], pert[3]),
                    s0,
                    "SoS inconsistent with concrete moment perturbation (κ={kappa})"
                );
            }
        }
    }

    #[test]
    fn filter_rates_general_position_resolves_fast() {
        let mut seed = 0x5EED_F117_0000_0AAAu64;
        let mut rnd = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((seed >> 11) as f64) / (1u64 << 53) as f64 * 100.0 - 50.0
        };
        let (mut filtered, mut total) = (0usize, 0usize);
        for _ in 0..2000 {
            let p = |r: &mut dyn FnMut() -> f64| [r(), r()];
            let (a, b, c, d) = (p(&mut rnd), p(&mut rnd), p(&mut rnd), p(&mut rnd));
            for stage in [
                orient2d_with_stage(a, b, c).1,
                incircle_with_stage(a, b, c, d).1,
            ] {
                total += 1;
                if stage == Stage::Filtered {
                    filtered += 1;
                }
            }
        }
        let rate = filtered as f64 / total as f64;
        println!(
            "{{\"suite\":\"fs-ivl/predicates\",\"metric\":\"filter_rate\",\
             \"value\":{rate:.4},\"n\":{total}}}"
        );
        assert!(
            rate > 0.99,
            "general-position inputs must resolve at stage A: {rate}"
        );
    }
}
