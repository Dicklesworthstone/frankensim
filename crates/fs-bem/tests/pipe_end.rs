//! The open end of a thin unflanged pipe on the open P1 screen (RULING .100
//! (2)): a measured ladder against Levine-Schwinger, cup-length independence
//! and two planted negatives. No accuracy gate: GreenPuma owns the gate and
//! its declaration. Run release, --test-threads=1.
//!
//! The ladder refines by about 1.27 per rung (36, 46, 58 vertices per ring),
//! so three rungs fit under the dense cap of 8192 facets (3114, 5060, 8062).
//! A ratio that close to one makes the observed-order estimate weak: its
//! changes are small and their ratio is noisy.
use std::time::Instant;

use fs_bem::helmholtz::Medium;
use fs_bem::reciprocal_radiation::open_screen::OpenOptions;
use fs_bem::reciprocal_radiation::pipe_end::{
    NODE_RING, PipeEndRow, PipeEndSpec, pipe_end_reflection, standing_wave_fit,
};
use fs_duct::levine_schwinger::{LsQuadrature, LsSolution, levine_schwinger};
use fs_math::c64::C64;

const KA: [f64; 7] = [0.05, 0.5, 1.0, 2.0, 3.0, 3.5, 3.8];
/// (circumferential, axial) per rung.
const RUNGS: [(usize, usize); 3] = [(36, 36), (46, 46), (58, 58)];
/// The mean refinement ratio of the rungs, sqrt(58 / 36).
const RATIO: f64 = 1.269_295_517_643_492;

fn options(level: u32) -> OpenOptions {
    OpenOptions {
        edge_level: level,
        maximum_edge_level: 8,
        triangle_order: 4,
        angular_refinement: 2,
        maximum_unknowns: 8192,
        maximum_cache_bytes: 8_000_000_000,
        maximum_work: u64::MAX,
        workers: 16,
    }
}

/// The declared cup: the sources and the window sit 1.5 radii from the
/// sources, the cap and the mouth, for every ka. The source rings sit on the
/// node ring, eight points each.
fn spec(length: f64, rung: (usize, usize), ring: f64) -> PipeEndSpec {
    PipeEndSpec {
        radius: 1.0,
        length,
        circumferential: rung.0,
        axial: rung.1,
        source: 0.3,
        source_spacing: 0.4,
        source_radius: NODE_RING,
        source_points: 8,
        window: (2.2, length - 1.5),
        stations: 16,
        ring_radius: ring,
    }
}

fn ls(ka: f64) -> LsSolution {
    levine_schwinger(ka, LsQuadrature { nodes: 1024 }).unwrap()
}

fn solve(s: &PipeEndSpec, ka: &[f64], level: u32) -> (Vec<PipeEndRow>, f64) {
    let start = Instant::now();
    let rows = pipe_end_reflection(s, ka, options(level), Medium::air(), || false).unwrap();
    (rows, start.elapsed().as_secs_f64())
}

fn print(case: &str, rung: usize, row: &PipeEndRow, seconds: f64) {
    let oracle = ls(row.ka);
    let r_ls = oracle.reflection();
    println!(
        "{{\"suite\":\"fs-bem\",\"case\":\"{case}\",\"rung\":{rung},\"edge_level\":{},\"unknowns\":{},\"ka\":{},\"R\":[{:.10e},{:.10e}],\"R_abs\":{:.8},\"R_LS\":[{:.10e},{:.10e}],\"R_err\":{:.4e},\"l_over_a\":{:.8},\"l_over_a_zero_crossing\":{:.8},\"l_over_a_LS\":{:.8},\"l_err\":{:.4e},\"l_err_zero_crossing\":{:.4e},\"fit_residual\":{:.3e},\"fit_condition\":{:.3},\"k_span\":{:.4},\"backward_residual\":{:.2e},\"seconds\":{seconds:.1}}}",
        row.edge_level,
        row.unknowns,
        row.ka,
        row.reflection.re,
        row.reflection.im,
        row.reflection.abs(),
        r_ls.re,
        r_ls.im,
        (row.reflection - r_ls).abs(),
        row.end_correction,
        row.zero_crossing_end_correction,
        oracle.l_over_a,
        (row.end_correction - oracle.l_over_a).abs(),
        (row.zero_crossing_end_correction - oracle.l_over_a).abs(),
        row.fit_residual,
        row.fit_condition,
        row.k_span,
        row.backward_residual,
    );
}

#[test]
fn open_pipe_end_ladder_against_levine_schwinger() {
    let mut ladder: Vec<Vec<PipeEndRow>> = Vec::new();
    for (rung, &mesh) in RUNGS.iter().enumerate() {
        let (rows, seconds) = solve(&spec(4.6, mesh, NODE_RING), &KA, rung as u32);
        for row in &rows {
            print("pipe-end-ladder", rung, row, seconds);
            assert!(row.reflection.re.is_finite() && row.reflection.im.is_finite());
        }
        ladder.push(rows);
    }
    // Signed Aitken limit of the three rungs (exact for v + c q^n with complex
    // q), its order for the mean ratio, and the limit against LS.
    for (j, &ka) in KA.iter().enumerate() {
        let [a, b, c] = [0, 1, 2].map(|r| ladder[r][j].reflection);
        let q = (c - b) / (b - a);
        let limit = c + (c - b) * q / (C64::ONE - q);
        let order = -q.abs().ln() / RATIO.ln();
        let [la, lb, lc] = [0, 1, 2].map(|r| ladder[r][j].end_correction);
        let ql = (lc - lb) / (lb - la);
        let l_limit = lc + (lc - lb) * ql / (1.0 - ql);
        let oracle = ls(ka);
        println!(
            "{{\"suite\":\"fs-bem\",\"case\":\"pipe-end-ladder-limit\",\"ka\":{ka},\"q\":[{:.4e},{:.4e}],\"observed_order\":{order:.3},\"ratio\":{RATIO:.4},\"order_note\":\"weak at a ratio this close to 1\",\"R_limit\":[{:.10e},{:.10e}],\"R_limit_err\":{:.4e},\"l_limit\":{l_limit:.8},\"l_limit_err\":{:.4e},\"finest_R_err\":{:.4e}}}",
            q.re,
            q.im,
            limit.re,
            limit.im,
            (limit - oracle.reflection()).abs(),
            (l_limit - oracle.l_over_a).abs(),
            (c - oracle.reflection()).abs(),
        );
    }
}

#[test]
fn open_pipe_end_does_not_depend_on_cup_length() {
    // The cap and sources change only the incident amplitude, so R must not
    // move with the cup's length. The axial count scales with the length to
    // keep the spacing.
    let ka = [1.0, 3.0];
    let (short, _) = solve(&spec(4.6, RUNGS[1], NODE_RING), &ka, 1);
    let (long, seconds) = solve(&spec(5.6, (RUNGS[1].0, 56), NODE_RING), &ka, 1);
    for (s, l) in short.iter().zip(&long) {
        print("pipe-end-length-4.6", 1, s, 0.0);
        print("pipe-end-length-5.6", 1, l, seconds);
        println!(
            "{{\"suite\":\"fs-bem\",\"case\":\"pipe-end-length-independence\",\"ka\":{},\"lengths\":[4.6,5.6],\"R_change\":{:.4e},\"l_change\":{:.4e}}}",
            s.ka,
            (s.reflection - l.reflection).abs(),
            (s.end_correction - l.end_correction).abs()
        );
    }
}

#[test]
fn open_pipe_end_planted_negatives() {
    // (a) Fitting under exp(+i omega t) swaps the waves: R -> 1/R, so l/a
    // changes sign and the LS comparison is off by orders of magnitude.
    let (rows, _) = solve(&spec(4.6, RUNGS[0], NODE_RING), &[0.5, 1.0], 0);
    for row in &rows {
        let k = row.ka;
        let honest = (row.reflection - ls(row.ka).reflection()).abs();
        let flipped = standing_wave_fit(-k, 1.0, 4.6, &row.stations, &row.pressures).unwrap();
        let wrong = (flipped.reflection - ls(row.ka).reflection()).abs();
        println!(
            "{{\"suite\":\"fs-bem\",\"case\":\"pipe-end-conjugate-convention\",\"ka\":{},\"l\":{:.6},\"l_conjugate\":{:.6},\"R_err\":{honest:.4e},\"R_err_conjugate\":{wrong:.4e}}}",
            row.ka, row.end_correction, flipped.end_correction
        );
        assert!(
            (flipped.end_correction + row.end_correction).abs() < 1e-9,
            "l/a did not change sign"
        );
        assert!(
            wrong > 0.03 && wrong > 30.0 * honest.max(1e-3),
            "the conjugate fit was not refused: {wrong:e} against {honest:e}"
        );
    }
    // (b) Receivers on the axis instead of the node ring: the first radial
    // mode is barely evanescent at ka >= 3.5 and must show in the residual.
    // Axial sources, which excite that mode, make the check physical.
    let ka = [3.6];
    let axial = |ring: f64| PipeEndSpec {
        source_radius: 0.0,
        source_points: 1,
        ..spec(4.6, RUNGS[0], ring)
    };
    let (node, _) = solve(&axial(NODE_RING), &ka, 0);
    let (axis, _) = solve(&axial(0.0), &ka, 0);
    print("pipe-end-node-ring", 0, &node[0], 0.0);
    print("pipe-end-axis", 0, &axis[0], 0.0);
    assert!(
        axis[0].fit_residual > 10.0 * node[0].fit_residual,
        "the axis residual {:e} does not rise over the node ring's {:e}",
        axis[0].fit_residual,
        node[0].fit_residual
    );
}

#[test]
fn open_pipe_end_node_ring_sources_at_the_radial_cut_on() {
    // Prediction, declared before the run: at ka 3.8 the axial sources feed
    // the barely evanescent first radial mode to the mouth, where it converts
    // into the reflected plane wave (R missed LS by 3.6-3.7e-2 on every rung).
    // Node-ring sources do not excite it, so they must cut the miss at least
    // threefold on the same rung.
    let ka = [3.5, 3.8];
    let ring_spec = spec(4.6, RUNGS[1], NODE_RING);
    let axial_spec = PipeEndSpec {
        source_radius: 0.0,
        source_points: 1,
        ..ring_spec
    };
    let (ring, seconds) = solve(&ring_spec, &ka, 1);
    let (axial, _) = solve(&axial_spec, &ka, 1);
    for (r, a) in ring.iter().zip(&axial) {
        print("pipe-end-ring-sources", 1, r, seconds);
        print("pipe-end-axial-sources", 1, a, 0.0);
    }
    let miss = |row: &PipeEndRow| (row.reflection - ls(row.ka).reflection()).abs();
    assert!(
        miss(&ring[1]) < miss(&axial[1]) / 3.0,
        "node-ring sources miss LS by {:e} at ka 3.8, axial ones by {:e}",
        miss(&ring[1]),
        miss(&axial[1])
    );
}
