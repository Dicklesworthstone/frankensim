//! fs-feec battery (tfz.5): dd = 0 exactly (integer path AND the f64
//! CSR materialization) on the fixture zoo, Betti rank–nullity
//! bookkeeping, de-Rham commutation R∘∇/curl/div = d∘R (the test that
//! pins every orientation sign), Whitney patch tests (constants and
//! affine fields reproduced exactly), mass-matrix partition-of-unity
//! and SPD checks, G1 MMS convergence for primal Poisson through the
//! FEEC stiffness composition, Hodge star positivity, and the
//! cross-ISA golden hash over assembled operators.

use fs_feec::{
    Cochain, betti_numbers, deram0, deram1, deram2, deram3, element_geometry, galerkin_star,
    hodge_diagonal_barycentric, incidence_to_csr, kuhn_cube, mass_matrix, on_unit_cube_boundary,
    single_tet, stiffness, two_tets,
};
use fs_qty::Dims;
use fs_rand::StreamKey;
use fs_rep_mesh::TetComplex;

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-feec\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

fn zoo() -> Vec<(&'static str, TetComplex, Vec<[f64; 3]>)> {
    let (c1, p1) = single_tet();
    let (c2, p2) = two_tets();
    let (c3, p3) = kuhn_cube(1);
    let (c4, p4) = kuhn_cube(2);
    let (c5, p5) = kuhn_cube(3);
    vec![
        ("single-tet", c1, p1),
        ("two-tets", c2, p2),
        ("kuhn-1", c3, p3),
        ("kuhn-2", c4, p4),
        ("kuhn-3", c5, p5),
    ]
}

fn cycle_integer_cochain(seed_values: &[i64], len: usize) -> Vec<i64> {
    if seed_values.is_empty() {
        vec![0; len]
    } else {
        (0..len)
            .map(|index| seed_values[index % seed_values.len()])
            .collect()
    }
}

#[test]
fn dd_is_zero_exactly_on_the_zoo() {
    let mut stream = StreamKey {
        seed: 5,
        kernel: 0xFEEC,
        tile: 1,
    }
    .stream();
    for (name, complex, _) in zoo() {
        let (d0, d1, d2) = (complex.d0(), complex.d1(), complex.d2());
        // Integer path: d(d(x)) == 0 for random integer cochains.
        for _ in 0..3 {
            let x0: Vec<i64> = (0..complex.vertex_count)
                .map(|_| i64::try_from(stream.next_below(2001)).expect("small") - 1000)
                .collect();
            assert!(
                d1.apply(&d0.apply(&x0)).iter().all(|&v| v == 0),
                "{name}: d1 d0 != 0 (integer)"
            );
            let x1: Vec<i64> = (0..complex.edges.len())
                .map(|_| i64::try_from(stream.next_below(2001)).expect("small") - 1000)
                .collect();
            assert!(
                d2.apply(&d1.apply(&x1)).iter().all(|&v| v == 0),
                "{name}: d2 d1 != 0 (integer)"
            );
        }
        // Materialized f64 path: the spgemm product must be EXACTLY
        // zero-valued (sums of ±1, exact in f64).
        let dd10 = fs_sparse::ops::spgemm(&incidence_to_csr(&d1), &incidence_to_csr(&d0));
        let dd21 = fs_sparse::ops::spgemm(&incidence_to_csr(&d2), &incidence_to_csr(&d1));
        let dense10 = dd10.to_dense();
        let dense21 = dd21.to_dense();
        assert!(dense10.iter().all(|&v| v == 0.0), "{name}: CSR d1 d0 != 0");
        assert!(dense21.iter().all(|&v| v == 0.0), "{name}: CSR d2 d1 != 0");
        log("dd-zero", "pass", name);
    }
}

/// G0 generated exact-sequence battery (bead frankensim-4nh8). The existing
/// three fixed Philox trials per fixture remain unchanged; this harness adds
/// replay seeds and shrinking over both fixture selection and cochain values.
#[test]
fn generated_integer_cochains_satisfy_dd_zero_on_the_zoo() {
    let fixtures = zoo();
    fs_propcheck::check(
        "feec-integer-dd-zero-on-fixture-zoo",
        0xFEEC_4A48_0001,
        600,
        |s| {
            (
                s.next_u64(),
                s.vec_of(32, |s| s.int_in(-1_000, 1_000)),
                s.vec_of(32, |s| s.int_in(-1_000, 1_000)),
            )
        },
        |(fixture_index, vertex_seed, edge_seed)| {
            let fixture_index = usize::try_from(
                *fixture_index % u64::try_from(fixtures.len()).expect("fixture count fits u64"),
            )
            .expect("reduced fixture index fits usize");
            let (_, complex, _) = &fixtures[fixture_index];
            let x0 = cycle_integer_cochain(vertex_seed, complex.vertex_count);
            let x1 = cycle_integer_cochain(edge_seed, complex.edges.len());
            let (d0, d1, d2) = (complex.d0(), complex.d1(), complex.d2());
            d1.apply(&d0.apply(&x0)).iter().all(|&value| value == 0)
                && d2.apply(&d1.apply(&x1)).iter().all(|&value| value == 0)
        },
    );
    log(
        "dd-zero-propcheck",
        "pass",
        "600 generated fixture/cochain cases, shrink-armed",
    );
}

#[test]
fn kuhn_fixture_is_positively_oriented_and_conforming() {
    for n in [1usize, 2, 3] {
        let (complex, positions) = kuhn_cube(n);
        let geo = element_geometry(&complex, &positions);
        assert!(
            geo.vol_signed.iter().all(|&v| v > 0.0),
            "kuhn({n}): non-positive stored-order volume"
        );
        let total: f64 = geo.vol_signed.iter().sum();
        assert!(
            (total - 1.0).abs() < 1e-12,
            "kuhn({n}): volumes sum to {total}, expected 1"
        );
        // Conforming: interior faces shared by exactly 2 tets, boundary
        // by 1 — count via d2 column occupancy.
        let d2 = complex.d2();
        let mut face_use = vec![0usize; complex.faces.len()];
        for row in &d2.rows {
            for &(f, _) in row {
                face_use[f] += 1;
            }
        }
        assert!(
            face_use.iter().all(|&c| c == 1 || c == 2),
            "kuhn({n}): non-conforming face incidence"
        );
        log(
            "kuhn-fixture",
            "pass",
            &format!("n={n} tets={} vol=1 exact-ish", complex.tets.len()),
        );
    }
}

#[test]
fn betti_numbers_of_ball_fixtures() {
    for (name, complex, _) in zoo() {
        let b = betti_numbers(&complex);
        assert_eq!(b, [1, 0, 0, 0], "{name}: Betti {b:?}");
        log("betti", "pass", &format!("{name} -> {b:?}"));
    }
}

/// Quadratic scalar field and its exact gradient.
fn f_quad(p: [f64; 3]) -> f64 {
    let q = p[0] * p[1] + 0.5 * p[2] * p[2];
    2.0f64.mul_add(p[0], q) - 0.7 * p[1]
}
fn grad_f_quad(p: [f64; 3]) -> [f64; 3] {
    [p[1] + 2.0, p[0] - 0.7, p[2]]
}

/// Affine vector field, its curl (constant) and divergence (constant).
fn a_affine(p: [f64; 3]) -> [f64; 3] {
    [
        0.5f64.mul_add(p[1], 1.0) - 0.25 * p[2],
        2.0f64.mul_add(p[2], -p[0]) + 0.5,
        0.3f64.mul_add(p[0], p[1]) - 2.0,
    ]
}
const CURL_A: [f64; 3] = [1.0 - 2.0, -0.25 - 0.3, -1.0 - 0.5];
const DIV_A: f64 = 0.0 + 0.0 + 0.0; // components have no self-derivative
fn b_affine(p: [f64; 3]) -> [f64; 3] {
    [
        1.5f64.mul_add(p[0], 0.2 * p[1]),
        0.7f64.mul_add(p[1], -0.4 * p[2]) + 1.0,
        2.0f64.mul_add(p[2], 0.1 * p[0]) - 0.5,
    ]
}
const DIV_B: f64 = 1.5 + 0.7 + 2.0;

#[test]
fn deram_maps_commute_with_d() {
    // THE orientation test: R∘(grad/curl/div) = d∘R pins every sign
    // convention (edge direction, face circulation, cell parity).
    let (complex, positions) = kuhn_cube(2);
    let geo = element_geometry(&complex, &positions);
    // Gradient: d0 R0(f) == R1(∇f) for quadratic f (Simpson exact).
    let r0 = deram0(&positions, &f_quad);
    let d0 = complex.d0();
    let mut worst = 0.0f64;
    let r1_grad = deram1(&complex, &positions, &grad_f_quad);
    for (e, row) in d0.rows.iter().enumerate() {
        let mut acc = 0.0f64;
        for &(v, s) in row {
            acc += f64::from(s) * r0[v];
        }
        worst = worst.max((acc - r1_grad[e]).abs());
    }
    assert!(worst < 1e-13, "grad commutation residual {worst:.3e}");
    // Curl: d1 R1(A) == R2(curl A) for affine A.
    let r1 = deram1(&complex, &positions, &a_affine);
    let r2_curl = deram2(&complex, &positions, &|_| CURL_A);
    let d1 = complex.d1();
    let mut worst_c = 0.0f64;
    for (f, row) in d1.rows.iter().enumerate() {
        let mut acc = 0.0f64;
        for &(e, s) in row {
            acc += f64::from(s) * r1[e];
        }
        worst_c = worst_c.max((acc - r2_curl[f]).abs());
    }
    assert!(worst_c < 1e-13, "curl commutation residual {worst_c:.3e}");
    // Divergence: d2 R2(B) == R3(div B) for affine B; and the affine A
    // above happens to be divergence-free — both are exercised.
    let d2 = complex.d2();
    for (name, field, divv) in [
        ("B", &b_affine as &dyn Fn([f64; 3]) -> [f64; 3], DIV_B),
        ("A", &a_affine as &dyn Fn([f64; 3]) -> [f64; 3], DIV_A),
    ] {
        let r2 = deram2(&complex, &positions, &|p| field(p));
        let r3_div = deram3(&complex, &positions, &geo, &|_| divv);
        let mut worst_d = 0.0f64;
        for (t, row) in d2.rows.iter().enumerate() {
            let mut acc = 0.0f64;
            for &(f, s) in row {
                acc += f64::from(s) * r2[f];
            }
            worst_d = worst_d.max((acc - r3_div[t]).abs());
        }
        assert!(
            worst_d < 1e-13,
            "div({name}) commutation residual {worst_d:.3e}"
        );
    }
    log(
        "commutation",
        "pass",
        "grad/curl/div all commute with d (signs pinned)",
    );
}

#[test]
fn whitney_masses_reproduce_constants_exactly() {
    // Λᵏ patch tests: Whitney interpolation of a CONSTANT field has
    // exactly the right energy — ⟨R c, M_k R c⟩ = |c|²·V — because
    // lowest-order Whitney spaces contain constants.
    let (complex, positions) = kuhn_cube(2);
    let geo = element_geometry(&complex, &positions);
    let energy = |m: &fs_sparse::Csr, x: &[f64]| -> f64 {
        let mut y = vec![0.0f64; x.len()];
        m.spmv(x, &mut y);
        x.iter().zip(&y).map(|(a, b)| a * b).sum()
    };
    // k = 0: constant scalar 3.
    let m0 = mass_matrix(&complex, &geo, 0);
    let r0 = vec![3.0f64; complex.vertex_count];
    assert!((energy(&m0, &r0) - 9.0).abs() < 1e-12, "M0 constant energy");
    // k = 1: constant vector c.
    let c = [0.8f64, -1.1, 0.6];
    let m1 = mass_matrix(&complex, &geo, 1);
    let r1 = deram1(&complex, &positions, &|_| c);
    let c2 = c[0].mul_add(c[0], c[1].mul_add(c[1], c[2] * c[2]));
    let e1 = energy(&m1, &r1);
    assert!((e1 - c2).abs() < 1e-12, "M1 constant energy: {e1} vs {c2}");
    // k = 2: constant vector b.
    let b = [0.4f64, 1.3, -0.9];
    let m2 = mass_matrix(&complex, &geo, 2);
    let r2 = deram2(&complex, &positions, &|_| b);
    let b2 = b[0].mul_add(b[0], b[1].mul_add(b[1], b[2] * b[2]));
    let e2 = energy(&m2, &r2);
    assert!((e2 - b2).abs() < 1e-12, "M2 constant energy: {e2} vs {b2}");
    // k = 3: constant density.
    let m3 = mass_matrix(&complex, &geo, 3);
    let r3 = deram3(&complex, &positions, &geo, &|_| 2.5);
    let e3 = energy(&m3, &r3);
    assert!((e3 - 6.25).abs() < 1e-12, "M3 constant energy: {e3}");
    // SPD: positive energy on deterministic pseudo-random vectors.
    let mut stream = StreamKey {
        seed: 9,
        kernel: 0xFEEC,
        tile: 2,
    }
    .stream();
    for (k, m) in [(0u8, &m0), (1, &m1), (2, &m2), (3, &m3)] {
        let x: Vec<f64> = (0..m.nrows())
            .map(|_| 2.0f64.mul_add(stream.next_f64(), -1.0))
            .collect();
        assert!(
            energy(m, &x) > 0.0,
            "M{k} not positive definite on random vector"
        );
    }
    log(
        "whitney-patch",
        "pass",
        "constant fields exact for k=0..3, masses PD",
    );
}

#[test]
fn stiffness_composition_kills_affine_fields() {
    // K0 = d0ᵀ·M1·d0: interior rows annihilate affine vertex data (the
    // classical patch test, derived from the complex).
    let (complex, positions) = kuhn_cube(2);
    let geo = element_geometry(&complex, &positions);
    let k0 = stiffness(
        &incidence_to_csr(&complex.d0()),
        &mass_matrix(&complex, &geo, 1),
    );
    let u: Vec<f64> = positions
        .iter()
        .map(|&p| 0.3f64.mul_add(p[0], 1.2f64.mul_add(p[1], -0.5 * p[2])) + 2.0)
        .collect();
    let mut ku = vec![0.0f64; u.len()];
    k0.spmv(&u, &mut ku);
    let mut worst = 0.0f64;
    for (v, &r) in ku.iter().enumerate() {
        if !on_unit_cube_boundary(positions[v]) {
            worst = worst.max(r.abs());
        }
    }
    assert!(worst < 1e-12, "affine patch residual {worst:.3e}");
    // Symmetry of the composition.
    let dense = k0.to_dense();
    let n = k0.nrows();
    let mut asym = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            asym = asym.max((dense[i * n + j] - dense[j * n + i]).abs());
        }
    }
    assert!(asym < 1e-14, "K0 asymmetry {asym:.3e}");
    log(
        "stiffness-patch",
        "pass",
        &format!("affine residual {worst:.2e}, asym {asym:.2e}"),
    );
}

#[test]
fn mms_poisson_primal_converges_at_second_order() {
    // G1: −Δu = f, u = sin(πx)sin(πy)sin(πz) (test-side std oracle),
    // homogeneous Dirichlet; FEEC stiffness d0ᵀM1d0; M0-weighted L2
    // error must fall at O(h²).
    let pi = std::f64::consts::PI;
    let u_exact = |p: [f64; 3]| (pi * p[0]).sin() * (pi * p[1]).sin() * (pi * p[2]).sin();
    let f_rhs = move |p: [f64; 3]| 3.0 * pi * pi * u_exact(p);
    let mut errs = Vec::new();
    for n in [4usize, 8, 16] {
        let (complex, positions) = kuhn_cube(n);
        let geo = element_geometry(&complex, &positions);
        let m0 = mass_matrix(&complex, &geo, 0);
        let k0 = stiffness(
            &incidence_to_csr(&complex.d0()),
            &mass_matrix(&complex, &geo, 1),
        );
        // RHS b = M0 · R0(f).
        let r0f = deram0(&positions, &f_rhs);
        let mut b = vec![0.0f64; r0f.len()];
        m0.spmv(&r0f, &mut b);
        // Reduce to interior vertices and PCG-solve (SPD stiffness;
        // solver tolerance far below discretization error).
        let interior: Vec<usize> = (0..positions.len())
            .filter(|&v| !on_unit_cube_boundary(positions[v]))
            .collect();
        let ni = interior.len();
        let mut slot = vec![usize::MAX; positions.len()];
        for (i, &v) in interior.iter().enumerate() {
            slot[v] = i;
        }
        let mut red = fs_sparse::Coo::new(ni, ni);
        for (i, &v) in interior.iter().enumerate() {
            let (cols, vals) = k0.row(v);
            for (&c, &val) in cols.iter().zip(vals) {
                if slot[c] != usize::MAX {
                    red.push(i, slot[c], val);
                }
            }
        }
        let a = red.assemble();
        let rhs: Vec<f64> = interior.iter().map(|&v| b[v]).collect();
        let mut x = vec![0.0f64; ni];
        let report = fs_sparse::precond::pcg(
            &a,
            &rhs,
            &mut x,
            &fs_sparse::precond::IdentityPrecond,
            1e-12,
            10_000,
        );
        assert!(report.converged, "PCG failed at n={n}: {report:?}");
        // M0-weighted L2 error over ALL vertices (boundary exact).
        let mut e = vec![0.0f64; positions.len()];
        for (i, &v) in interior.iter().enumerate() {
            e[v] = x[i] - u_exact(positions[v]);
        }
        let mut me = vec![0.0f64; e.len()];
        m0.spmv(&e, &mut me);
        let l2: f64 = e.iter().zip(&me).map(|(a, b)| a * b).sum::<f64>().sqrt();
        errs.push(l2);
        log("mms-poisson", "info", &format!("n={n} L2={l2:.4e}"));
    }
    let (o1, o2) = ((errs[0] / errs[1]).log2(), (errs[1] / errs[2]).log2());
    assert!(
        (o1 - 2.0).abs() < 0.5 && (o2 - 2.0).abs() < 0.35,
        "MMS orders {o1:.2}, {o2:.2} (errors {errs:?})"
    );
    log("mms-order", "pass", &format!("orders {o1:.2}, {o2:.2}"));
}

#[test]
fn hodge_stars_positive_and_consistent() {
    let (complex, positions) = kuhn_cube(2);
    let geo = element_geometry(&complex, &positions);
    for k in 0..=3u8 {
        let diag = hodge_diagonal_barycentric(&complex, &positions, &geo, k);
        assert!(
            diag.iter().all(|&v| v > 0.0),
            "diagonal star {k} not positive"
        );
        let gal = galerkin_star(&complex, &geo, k);
        assert_eq!(gal.nrows(), diag.len());
        log(
            "hodge",
            "pass",
            &format!(
                "k={k} diag range [{:.3e}, {:.3e}]",
                diag.iter().copied().fold(f64::INFINITY, f64::min),
                diag.iter().copied().fold(0.0f64, f64::max)
            ),
        );
    }
    // Total dual volume at k=0 equals the domain volume.
    let d0 = hodge_diagonal_barycentric(&complex, &positions, &geo, 0);
    let total: f64 = d0.iter().sum();
    assert!(
        (total - 1.0).abs() < 1e-12,
        "dual volumes must sum to |Ω|: {total}"
    );
}

#[test]
fn cochain_container_semantics() {
    let (complex, _positions) = kuhn_cube(2);
    let vals: Vec<f64> = (0..complex.vertex_count).map(|v| v as f64 * 0.5).collect();
    let dims = Dims([0, 0, 1, 0, 0]); // seconds, say
    let c = Cochain::from_values(&complex, 0, &vals, dims);
    assert_eq!(c.degree(), 0);
    assert_eq!(c.dims(), dims);
    let dc = c.d(&complex);
    assert_eq!(dc.degree(), 1);
    assert_eq!(dc.dims(), dims, "d preserves the dimension tag");
    assert_eq!(dc.values().len(), complex.edges.len());
    // dd = 0 through the container too.
    let ddc = dc.d(&complex);
    assert!(ddc.values().iter().all(|&v| v == 0.0), "container dd != 0");
    let view = c.view();
    assert_eq!(view.name, "cochain0");
    assert_eq!(view.addr % 128, 0, "cochain buffer misaligned");
    log("cochain", "pass", "dims tag + view + container dd=0");
}

const GOLDEN_HASH: u64 = 0xa973_ca6b_07c3_9639; // recorded at tfz.5 landing, frozen

#[test]
fn feec_golden_hash() {
    // Hash over PRODUCTION-PATH values only: assembled operator
    // entries and Betti counts on kuhn(2). Positions are rational,
    // assembly is +/−/×/÷ — no libm anywhere in the hashed pipeline.
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let (complex, positions) = kuhn_cube(2);
    let geo = element_geometry(&complex, &positions);
    for k in 0..=3u8 {
        let m = mass_matrix(&complex, &geo, k);
        let dense = m.to_dense();
        for v in dense.iter().step_by(7) {
            feed(*v);
        }
    }
    let k0 = stiffness(
        &incidence_to_csr(&complex.d0()),
        &mass_matrix(&complex, &geo, 1),
    );
    for v in k0.to_dense().iter().step_by(11) {
        feed(*v);
    }
    for b in betti_numbers(&complex) {
        feed(f64::from(u32::try_from(b).expect("small")));
    }
    for v in geo.vol_signed.iter().take(16) {
        feed(*v);
    }
    log("feec-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "feec bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}

#[test]
fn on_unit_cube_boundary_detects_reconstructed_far_faces() {
    // Regression: `kuhn_cube(n)` reconstructs the far face as `n * (1.0/n)`,
    // which is 1 ULP below 1.0 for many n (first at n = 49). A bit-exact
    // boundary check treated those face vertices as INTERIOR — under Dirichlet
    // pinning that yields a near-singular reduced system. `far` below is EXACTLY
    // the coordinate `kuhn_cube` stores for the i = n plane.
    for &n in &[49usize, 98, 103, 107] {
        let far = n as f64 * (1.0 / n as f64);
        assert!(
            far < 1.0,
            "n={n}: far-face reconstruction must be < 1.0 (the bug trigger)"
        );
        assert!(
            on_unit_cube_boundary([far, 0.5, 0.5]),
            "n={n}: x = {far} (the i=n plane) must be on the cube boundary"
        );
        assert!(on_unit_cube_boundary([0.5, far, 0.5]), "n={n}: y = {far}");
        assert!(on_unit_cube_boundary([0.5, 0.5, far]), "n={n}: z = {far}");
    }
    // Exact endpoints/origin stay boundary; strictly-interior vertices (≥ 1/n
    // from every face) stay interior — no false positives.
    assert!(on_unit_cube_boundary([0.0, 0.5, 0.5]));
    assert!(on_unit_cube_boundary([1.0, 0.5, 0.5]));
    assert!(!on_unit_cube_boundary([0.5, 0.5, 0.5]));
    assert!(!on_unit_cube_boundary([1.0 / 49.0, 0.5, 0.5])); // one lattice step in
    log(
        "unit-cube-boundary",
        "pass",
        "reconstructed far faces detected; interior unaffected",
    );
}
