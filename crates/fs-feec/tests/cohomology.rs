//! Cohomology conformance (the tfz.7 bead). Acceptance: harmonic
//! dimensions match Betti numbers across the fixture zoo (the
//! geometry/physics internal consistency check); decomposition
//! components are M-orthogonal and re-sum exactly; the circulation
//! payoff recovers Γ (and Kutta–Joukowski lift) from the harmonic
//! component of a sampled vortex flow within discretization error,
//! stable under refinement; harmonic deflation kills the circulation a
//! saddle solver must not see.

use fs_feec::whitney::element_geometry;
use fs_feec::{
    betti_numbers, circulation, deflate_harmonics, harmonic_basis, hodge_decompose, kuhn_cube,
    masked_cube_grid,
};
use fs_rep_mesh::TetComplex;

const SUITE: &str = "fs-feec/cohomology";
const FIXED_INPUT_SEED: u64 = 0;

fn verdict(case: &str, detail: &str) {
    let mut emitter = fs_obs::Emitter::new(SUITE, case);
    let event = emitter.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass: true,
            detail: detail.to_string(),
            seed: FIXED_INPUT_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("cohomology verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("cohomology verdict must use the fs-obs wire schema");
    println!("{line}");
}

fn measurement(case: &str, name: &str, json: String) {
    let identity = format!("{case}/measurement");
    let mut emitter = fs_obs::Emitter::new(SUITE, &identity);
    let event = emitter.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::Custom {
            name: name.to_string(),
            json,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("cohomology measurement must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("cohomology measurement must use the fs-obs wire schema");
    println!("{line}");
}

fn finite_json(value: f64) -> String {
    if value.is_finite() {
        value.to_string()
    } else {
        "null".to_string()
    }
}

/// Ring slab: nx × nx × 1 cells with a centered hole_w × hole_w hole.
fn ring(nx: usize, hole_w: usize) -> (TetComplex, Vec<[f64; 3]>) {
    let lo = (nx - hole_w) / 2;
    let hi = lo + hole_w;
    masked_cube_grid(nx, nx, 1, &|i, j, _| {
        !(i >= lo && i < hi && j >= lo && j < hi)
    })
}

#[test]
fn ch_001_harmonic_dimensions_match_betti() {
    type Fixture = (TetComplex, Vec<[f64; 3]>);
    // The zoo: ball, ring (1 handle), two-hole slab (2 handles),
    // hollow cube (1 void).
    let ball = kuhn_cube(2);
    let ring1 = ring(4, 2);
    let two_hole = masked_cube_grid(5, 3, 1, &|i, j, _| !((i == 1 || i == 3) && j == 1));
    let hollow = masked_cube_grid(3, 3, 3, &|i, j, k| !(i == 1 && j == 1 && k == 1));
    let zoo: [(&str, &Fixture, [usize; 3]); 4] = [
        ("ball", &ball, [1, 0, 0]),
        ("ring", &ring1, [1, 1, 0]),
        ("two-hole", &two_hole, [1, 2, 0]),
        ("hollow", &hollow, [1, 0, 1]),
    ];
    for (name, (complex, positions), expect) in zoo {
        let betti = betti_numbers(complex);
        assert_eq!(
            [betti[0], betti[1], betti[2]],
            expect,
            "{name}: integer Betti"
        );
        let geo = element_geometry(complex, positions);
        for degree in 0..3u8 {
            let basis = harmonic_basis(
                complex,
                positions,
                &geo,
                degree,
                expect[degree as usize] + 2,
            );
            assert_eq!(
                basis.len(),
                expect[degree as usize],
                "{name}: harmonic dim at degree {degree} must equal b_{degree}"
            );
        }
    }
    verdict(
        "ch-001",
        "harmonic dimensions equal Betti numbers on ball/ring/two-hole/hollow (degrees \
         0..2) — geometry and physics agree",
    );
}

#[test]
fn ch_002_decomposition_orthogonal_and_resums() {
    let (complex, positions) = ring(4, 2);
    let geo = element_geometry(&complex, &positions);
    let n_edges = complex.edges.len();
    // A deterministic messy 1-cochain.
    let x: Vec<f64> = (0..n_edges)
        .map(|k| ((k * 2654435761) % 1000) as f64 / 500.0 - 1.0)
        .collect();
    let parts = hodge_decompose(&complex, &positions, &geo, 1, &x);
    // Re-sum EXACTLY (harmonic is defined as the remainder).
    for (i, ((e, c), h)) in parts
        .exact
        .iter()
        .zip(&parts.coexact)
        .zip(&parts.harmonic)
        .enumerate()
    {
        assert!((e + c + h - x[i]).abs() < 1e-12, "resum at {i}");
    }
    // M-orthogonality to solver tolerance.
    for (k, r) in parts.ortho_residuals.iter().enumerate() {
        assert!(*r < 1e-8, "orthogonality residual {k}: {r}");
    }
    // Projection idempotence: decomposing the exact part returns it.
    let again = hodge_decompose(&complex, &positions, &geo, 1, &parts.exact);
    let mut drift = 0.0f64;
    for i in 0..n_edges {
        drift = drift.max((again.exact[i] - parts.exact[i]).abs());
        drift = drift.max(again.harmonic[i].abs());
    }
    assert!(drift < 1e-8, "idempotence drift: {drift}");
    verdict(
        "ch-002",
        "exact+coexact+harmonic re-sums to 1e-12; orthogonality residuals < 1e-8; \
         projection idempotent",
    );
}

/// Sample a velocity field as a 1-cochain by midpoint line integrals.
fn sample_flow(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    u: &dyn Fn([f64; 3]) -> [f64; 3],
) -> Vec<f64> {
    complex
        .edges
        .iter()
        .map(|&[a, b]| {
            let (pa, pb) = (positions[a as usize], positions[b as usize]);
            let mid = [
                f64::midpoint(pa[0], pb[0]),
                f64::midpoint(pa[1], pb[1]),
                f64::midpoint(pa[2], pb[2]),
            ];
            let vel = u(mid);
            vel[0] * (pb[0] - pa[0]) + vel[1] * (pb[1] - pa[1]) + vel[2] * (pb[2] - pa[2])
        })
        .collect()
}

/// An axis-aligned square edge loop at height z = 0, corners
/// (lo, lo) → (hi, lo) → (hi, hi) → (lo, hi), counter-clockwise.
fn square_loop(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    lo: f64,
    hi: f64,
) -> Vec<(usize, f64)> {
    let vert_at = |x: f64, y: f64| -> u32 {
        positions
            .iter()
            .position(|p| (p[0] - x).abs() < 1e-9 && (p[1] - y).abs() < 1e-9 && p[2].abs() < 1e-9)
            .map(|i| u32::try_from(i).expect("fits"))
            .expect("grid vertex present")
    };
    let edge_of = |a: u32, b: u32| -> (usize, f64) {
        let key = [a.min(b), a.max(b)];
        let idx = complex
            .edges
            .binary_search(&key)
            .expect("grid edge present");
        (idx, if a < b { 1.0 } else { -1.0 })
    };
    let mut cycle = Vec::new();
    let (loi, hii) = (lo as i64, hi as i64);
    let mut walk = |x0: i64, y0: i64, x1: i64, y1: i64| {
        #[allow(clippy::cast_precision_loss)]
        let (a, b) = (vert_at(x0 as f64, y0 as f64), vert_at(x1 as f64, y1 as f64));
        cycle.push(edge_of(a, b));
    };
    for x in loi..hii {
        walk(x, loi, x + 1, loi);
    }
    for y in loi..hii {
        walk(hii, y, hii, y + 1);
    }
    for x in (loi..hii).rev() {
        walk(x + 1, hii, x, hii);
    }
    for y in (loi..hii).rev() {
        walk(loi, y + 1, loi, y);
    }
    cycle
}

fn gamma_recovery(nx: usize, hole_w: usize, loop_lo: f64, loop_hi: f64) -> (f64, f64) {
    let (complex, positions) = ring(nx, hole_w);
    let geo = element_geometry(&complex, &positions);
    #[allow(clippy::cast_precision_loss)]
    let c = nx as f64 / 2.0;
    let gamma = 2.5f64;
    let v_inf = 1.0f64;
    // Uniform flow + a Γ-vortex around the hole axis.
    let u = move |p: [f64; 3]| -> [f64; 3] {
        let (dx, dy) = (p[0] - c, p[1] - c);
        let r2 = (dx * dx + dy * dy).max(1e-12);
        let k = gamma / (2.0 * std::f64::consts::PI * r2);
        [v_inf - k * dy, k * dx, 0.0]
    };
    let x = sample_flow(&complex, &positions, &u);
    let parts = hodge_decompose(&complex, &positions, &geo, 1, &x);
    let cycle = square_loop(&complex, &positions, loop_lo, loop_hi);
    let rec = circulation(&parts.harmonic, &cycle);
    (rec, gamma)
}

#[test]
fn ch_003_circulation_to_lift_g1() {
    // Coarse ring and a refined ring: Γ extracted from the HARMONIC
    // component within discretization error, improving (or holding)
    // under refinement — then Kutta–Joukowski turns it into lift.
    let (rec_coarse, gamma) = gamma_recovery(6, 2, 1.0, 5.0);
    let err_coarse = (rec_coarse - gamma).abs() / gamma;
    let (rec_fine, _) = gamma_recovery(12, 4, 2.0, 10.0);
    let err_fine = (rec_fine - gamma).abs() / gamma;
    assert!(
        err_coarse < 0.08,
        "coarse Γ recovery within 8%: {rec_coarse} vs {gamma} ({err_coarse})"
    );
    assert!(
        err_fine < 0.04,
        "refined Γ recovery within 4%: {rec_fine} vs {gamma} ({err_fine})"
    );
    assert!(
        err_fine <= err_coarse + 1e-9,
        "refinement does not regress: {err_fine} vs {err_coarse}"
    );
    // The payoff line: L' = ρ V Γ.
    let (rho, v_inf) = (1.225, 1.0);
    let lift = rho * v_inf * rec_fine;
    let lift_ref = rho * v_inf * gamma;
    measurement(
        "ch-003",
        "circulation-lift",
        format!(
            "{{\"gamma_ref\":{},\"gamma_coarse\":{},\"gamma_fine\":{},\"lift\":{},\"lift_ref\":{}}}",
            finite_json(gamma),
            finite_json(rec_coarse),
            finite_json(rec_fine),
            finite_json(lift),
            finite_json(lift_ref)
        ),
    );
    assert!((lift - lift_ref).abs() / lift_ref < 0.04);
    verdict(
        "ch-003",
        "harmonic circulation recovers Γ to 4% on the refined ring (8% coarse, \
         monotone) and Kutta–Joukowski lift follows",
    );
}

#[test]
fn ch_004_deflation_makes_saddles_well_posed() {
    let (complex, positions) = ring(4, 2);
    let geo = element_geometry(&complex, &positions);
    let basis = harmonic_basis(&complex, &positions, &geo, 1, 3);
    assert_eq!(basis.len(), 1, "one handle, one harmonic 1-cochain");
    // A cochain with a deliberate harmonic component.
    let n = complex.edges.len();
    let mut x: Vec<f64> = (0..n).map(|k| ((k * 48271) % 100) as f64 / 100.0).collect();
    for (xi, hi) in x.iter_mut().zip(&basis[0]) {
        *xi += 3.0 * hi;
    }
    let cycle = square_loop(&complex, &positions, 0.0, 4.0);
    let before = hodge_decompose(&complex, &positions, &geo, 1, &x);
    let circ_before = circulation(&before.harmonic, &cycle).abs();
    assert!(circ_before > 0.1, "the seeded handle mode is visible");
    // Deflate: the component a saddle solver must not see is GONE.
    let deflated = deflate_harmonics(&complex, &positions, &geo, 1, &basis, &x);
    let after = hodge_decompose(&complex, &positions, &geo, 1, &deflated);
    let circ_after = circulation(&after.harmonic, &cycle).abs();
    assert!(
        circ_after < 1e-8 * circ_before.max(1.0),
        "deflation kills the harmonic circulation: {circ_after}"
    );
    verdict(
        "ch-004",
        "harmonic deflation removes the handle mode (circulation 0.1+ -> ~0): the \
         well-posedness constraint rows work",
    );
}
