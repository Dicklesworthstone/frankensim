//! M1.2-W3 (bead pexfp.52): shell-beam ties. This file holds the Timoshenko beam element gate,
//! `g1_timoshenko_beam_element`:
//! (a) as G -> infinity the free-free bending frequencies converge to those of the Euler-Bernoulli beam,
//!     with a difference that falls as 1/G (slope -1 +- 0.2 per decade of G);
//! (b) the free-free bending frequencies of a refined chain lie within 0.1 % of the roots of the
//!     Timoshenko frequency equation (shear deformation and rotary inertia), set up independently here
//!     from the two coupled Timoshenko equations (the theory of Han, Benaroya and Wei 1999, J. Sound
//!     Vib. 225:935): a 4x4 boundary determinant on the cosh/sinh/cos/sin basis.
//!
//! Fixture: a straight wooden bar 300 mm long, 20 mm (y) x 30 mm (z), E_L 12 GPa, G 0.75 GPa,
//! 450 kg/m3. The shear areas are 5/6 A, the classical value of a rectangle, as section data.
//! Bending is in the x-y plane (v, theta_z), which a straight bar decouples exactly from the rest.
use fs_modal::eigh_gen_dense;
use fs_plate::shell::stiffened::BeamSection;
use fs_plate::shell::stiffened::timoshenko::{BeamKind, ShearAreas, assemble_beam_chain};

const L: f64 = 0.3;
const HY: f64 = 0.02;
const HZ: f64 = 0.03;
const E: f64 = 12e9;
const G: f64 = 0.75e9;
const RHO: f64 = 450.0;

fn section(shear_pa: f64) -> BeamSection {
    let a = HY * HZ;
    BeamSection {
        young_pa: E,
        shear_pa,
        density_kg_m3: RHO,
        area_m2: a,
        iy_m4: HY * HZ.powi(3) / 12.0,
        iz_m4: HZ * HY.powi(3) / 12.0,
        // Saint-Venant constant of a 30 x 20 mm rectangle, beta a b^3 with beta = 0.196 (a/b = 1.5);
        // torsion does not enter the x-y bending block.
        torsion_m4: 0.196 * HZ * HY.powi(3),
    }
}

fn shear_areas() -> ShearAreas {
    let a = HY * HZ;
    ShearAreas {
        y_m2: 5.0 / 6.0 * a,
        z_m2: 5.0 / 6.0 * a,
    }
}

/// Free-free x-y bending frequencies [Hz] of a chain of `n` equal segments, rigid pair removed.
fn chain_bending_hz(n: usize, shear_pa: f64, kind: BeamKind, count: usize) -> Vec<f64> {
    let nodes: Vec<[f64; 3]> = (0..=n)
        .map(|i| [L * i as f64 / n as f64, 0.0, 0.0])
        .collect();
    let (k, m) = assemble_beam_chain(&nodes, [0.0, 0.0, 1.0], section(shear_pa), kind).unwrap();
    let dofs: Vec<usize> = (0..=n).flat_map(|i| [6 * i + 1, 6 * i + 5]).collect();
    let size = dofs.len();
    let dense = |a: &fs_sparse::Csr| -> Vec<f64> {
        let mut out = vec![0.0; size * size];
        for (r, &row) in dofs.iter().enumerate() {
            let (cols, vals) = a.row(row);
            for (&col, &val) in cols.iter().zip(vals) {
                if let Some(c) = dofs.iter().position(|&d| d == col) {
                    out[r * size + c] = val;
                }
            }
        }
        out
    };
    let modes = eigh_gen_dense(&dense(&k), &dense(&m), size).unwrap();
    // Rigid modes sit at roundoff (about eps x the largest eigenvalue, tens of rad^2/s^2 here); the first
    // elastic mode is near 1 kHz. A 10 Hz threshold separates them by orders of magnitude both ways.
    let rigid_limit = (2.0 * std::f64::consts::PI * 10.0).powi(2);
    let elastic: Vec<f64> = modes
        .iter()
        .filter(|p| p.lambda > rigid_limit)
        .map(|p| p.lambda.sqrt() / (2.0 * std::f64::consts::PI))
        .collect();
    assert_eq!(
        modes.len() - elastic.len(),
        2,
        "two rigid modes (translation and rotation) expected"
    );
    elastic.into_iter().take(count).collect()
}

/// 4x4 free-free boundary determinant of the Timoshenko beam at nondimensional frequency b,
/// b^2 = rho A L^4 omega^2 / (E I). Rows: moment (psi') and shear (w' - psi) at xi = 0 and xi = 1.
fn timoshenko_det(b: f64, r2: f64, s2: f64) -> f64 {
    let root = b * (b * b * (r2 - s2).powi(2) + 4.0).sqrt();
    let p = ((-b * b * (r2 + s2) + root) / 2.0).sqrt();
    let q = ((b * b * (r2 + s2) + root) / 2.0).sqrt();
    let bs = b * b * s2;
    let ap = (p * p + bs) / p;
    let bq = (bs - q * q) / q;
    // Each basis function: (w, w', psi, psi') at xi.
    let basis = |xi: f64| -> [[f64; 4]; 4] {
        let (ch, sh) = ((p * xi).cosh(), (p * xi).sinh());
        let (c, s) = ((q * xi).cos(), (q * xi).sin());
        [
            [ch, p * sh, ap * sh, ap * p * ch],
            [sh, p * ch, ap * ch, ap * p * sh],
            [c, -q * s, bq * s, bq * q * c],
            [s, q * c, -bq * c, bq * q * s],
        ]
    };
    let (f0, f1) = (basis(0.0), basis(1.0));
    let mut a = [[0.0; 4]; 4];
    for j in 0..4 {
        a[0][j] = f0[j][3];
        a[1][j] = f0[j][1] - f0[j][2];
        a[2][j] = f1[j][3];
        a[3][j] = f1[j][1] - f1[j][2];
    }
    det4(a)
}

fn det4(mut a: [[f64; 4]; 4]) -> f64 {
    let mut det = 1.0;
    for col in 0..4 {
        let pivot = (col..4)
            .max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))
            .unwrap();
        if a[pivot][col] == 0.0 {
            return 0.0;
        }
        if pivot != col {
            a.swap(pivot, col);
            det = -det;
        }
        det *= a[col][col];
        for row in col + 1..4 {
            let f = a[row][col] / a[col][col];
            for k in col..4 {
                a[row][k] -= f * a[col][k];
            }
        }
    }
    det
}

/// The first `count` roots b of the frequency equation, by a fine scan and bisection, as Hz.
fn timoshenko_bending_hz(count: usize) -> Vec<f64> {
    let s = section(G);
    let r2 = s.iz_m4 / (s.area_m2 * L * L);
    let s2 = E * s.iz_m4 / (shear_areas().y_m2 * G * L * L);
    let mut roots = Vec::new();
    let (mut b0, step) = (0.5, 1e-3);
    let mut d0 = timoshenko_det(b0, r2, s2);
    while roots.len() < count {
        let b1 = b0 + step;
        // Stay below the second spectrum (b^2 r^2 s^2 < 1), where the basis above is the right one.
        assert!(b1 * b1 * r2 * s2 < 1.0, "scan reached the second spectrum");
        let d1 = timoshenko_det(b1, r2, s2);
        if d0 == 0.0 || d0.signum() != d1.signum() {
            let (mut lo, mut hi, mut dlo) = (b0, b1, d0);
            for _ in 0..200 {
                let mid = 0.5 * (lo + hi);
                let dm = timoshenko_det(mid, r2, s2);
                if dm.signum() == dlo.signum() {
                    lo = mid;
                    dlo = dm;
                } else {
                    hi = mid;
                }
            }
            roots.push(0.5 * (lo + hi));
        }
        b0 = b1;
        d0 = d1;
    }
    // b^2 = rho A L^4 omega^2 / (E I), so omega = b sqrt(E I / (rho A L^4)); in the Euler-Bernoulli limit the
    // roots are b = (beta L)^2, 22.37, 61.67, ... for a free-free bar.
    let scale = (E * s.iz_m4 / (RHO * s.area_m2 * L.powi(4))).sqrt();
    roots
        .into_iter()
        .map(|b| b * scale / (2.0 * std::f64::consts::PI))
        .collect()
}

#[test]
fn g1_timoshenko_beam_element() {
    const MODES: usize = 4;
    // (b) refined chains against the frequency equation.
    let exact = timoshenko_bending_hz(MODES);
    let euler = chain_bending_hz(320, G, BeamKind::EulerBernoulli, MODES);
    println!(
        "Timoshenko frequency equation: {exact:.3?} Hz (Euler-Bernoulli chain, 320 segments: {euler:.3?} Hz)"
    );
    let mut previous: Option<Vec<f64>> = None;
    for n in [80, 160, 320] {
        let fe = chain_bending_hz(n, G, BeamKind::Timoshenko(shear_areas()), MODES);
        let dev: Vec<f64> = fe.iter().zip(&exact).map(|(f, x)| f / x - 1.0).collect();
        let dev_text: Vec<String> = dev.iter().map(|d| format!("{d:.2e}")).collect();
        println!("{n:>4} segments: {fe:.3?} Hz, deviation {dev_text:?}");
        if let Some(prev) = &previous {
            for (i, (now, before)) in dev.iter().zip(prev).enumerate() {
                assert!(
                    now.abs() < before.abs(),
                    "mode {i}: the deviation must fall with refinement"
                );
            }
        }
        if n == 320 {
            for (i, d) in dev.iter().enumerate() {
                assert!(
                    d.abs() < 1e-3,
                    "mode {i}: {d:.2e} against the Timoshenko frequency equation"
                );
            }
            // The fixture is meant to need the shear theory: Euler-Bernoulli is far off on every mode.
            for (i, (e, x)) in euler.iter().zip(&exact).enumerate() {
                assert!(
                    e / x - 1.0 > 1e-2,
                    "mode {i}: Euler-Bernoulli only {:.2e} above",
                    e / x - 1.0
                );
            }
        }
        previous = Some(dev);
    }
    // (a) G -> infinity: the Timoshenko chain approaches the Euler-Bernoulli chain as 1/G.
    let n = 40;
    let mut gaps = Vec::new();
    for g in [G * 1e1, G * 1e2, G * 1e3, G * 1e4] {
        let t = chain_bending_hz(n, g, BeamKind::Timoshenko(shear_areas()), MODES);
        let e = chain_bending_hz(n, g, BeamKind::EulerBernoulli, MODES);
        gaps.push(
            t.iter()
                .zip(&e)
                .map(|(a, b)| (a / b - 1.0).abs())
                .fold(0.0_f64, f64::max),
        );
    }
    let gap_text: Vec<String> = gaps.iter().map(|g| format!("{g:.3e}")).collect();
    println!("max relative gap to Euler-Bernoulli at G x 10, 100, 1000, 10000: {gap_text:?}");
    for pair in gaps.windows(2) {
        let slope = (pair[1] / pair[0]).log10();
        assert!(
            (slope + 1.0).abs() <= 0.2,
            "slope {slope:.3} per decade of G"
        );
    }
}
