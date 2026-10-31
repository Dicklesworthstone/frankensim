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
//!
//! It also holds the off-node beam-to-shell tie (`tie_beams_to_shells`, the linear 1-D dual mortar from
//! the beam's own nodes to a shell skin line, with the arm u + theta x d):
//! - `g1_matching_beam_tie_reproduces_on_node_beam`: beam nodes on the shell nodes (C = I) reproduce the
//!   bonded eccentric `ShellBeam` model (CST/DKT shell, Euler-Bernoulli beam) to 1e-9;
//! - `g3_on_node_beam_tie_matches_legacy_to_rounding`: the same on-node tie agrees with the consistent
//!   bonded `assemble_stiffened_shell_with_mass` pencil entry by entry within a derived rounding bound.
//!   The plan's "bitwise" form of this claim is falsified (see the test);
//! - `g3_legacy_bonded_beam_pencil_golden`: the bits of the legacy bonded pencil the piano boards consume,
//!   pinned under docs/GOLDEN_POLICY.md;
//! - `g1_off_node_brace_converges`: a straight and a curved Timoshenko brace on an MITC3+ plate, each
//!   with its own beam nodes off the shell nodes, converge to the on-node tie (beam nodes on the shell
//!   nodes). One-sided gate of the orchestrator's ruling (mail 6174): the relative frequency gap between
//!   the two falls strictly at every rung n = 16, 32, 64, 128, and its fitted order over n = 32..128 is at
//!   least 1.8 (faster than 2 to the same limit is not a defect). A planted negative, the same beam tied
//!   12.5 mm off its station (one element of the coarsest rung, fixed in space), must fail that gate.
//!   MITC3+ (Reissner-Mindlin) with a Timoshenko beam is the consistent pairing for the linear transfer;
//!   on Kirchhoff CST/DKT shells an off-node Euler-Bernoulli beam converges to a stiffer limit (see
//!   `tie_beams_to_shells`);
//! - the refusals: a beam path off the surface, a one-node beam, and the tie budgets at cap and cap+1.
//!
//! Brace fixture: a cantilevered plate 400 x 200 x 3 mm (clamped along x = 0), isotropic E 10 GPa,
//! nu 0.3, 450 kg/m3, structured mesh of 2n x n cells; a glued 6 x 10 mm spruce brace (E_L 13.8 GPa,
//! G 0.8 GPa, 460 kg/m3, shear areas 5/6 A) along the plate's mid line y = W/2 + c sin(pi x / L), with
//! c = 0 (straight) or 20 mm (curved: the whole mesh follows the same bump, so the brace path is a mesh
//! line). The off-node beam has 3 nodes per 2 shell edges, so only its ends are shell nodes.
use fs_mms::{RefinementLadder, fit_order};
use fs_modal::{SliceOptions, eigh_gen_dense, slice_window};
use fs_plate::PlateSection;
use fs_plate::shell::junction::{
    BeamOnShell, ShellComponent, TieBudget, TieSide, tie_beams_to_shells, tie_shells,
};
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, NodalDirectors, ShearShellSection, ShellConstraint, ShellMass,
    assemble_mitc3p,
};
use fs_plate::shell::stiffened::BeamSection;
use fs_plate::shell::stiffened::timoshenko::{
    BeamKind, ShearAreas, assemble_beam_chain, segment_matrices,
};
use fs_plate::shell::stiffened::{
    BeamMass, ShellBeam, assemble_stiffened_shell, assemble_stiffened_shell_with_mass,
};
use fs_plate::shell::{ShellMesh, ShellModel, ShellSupport, assemble_shell_sections};

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

// ---------------------------------------------------------------- brace fixture
const PL_L: f64 = 0.4;
const PL_W: f64 = 0.2;
const PL_H: f64 = 0.003;
const PL_E: f64 = 10e9;
const PL_NU: f64 = 0.3;
const PL_RHO: f64 = 450.0;
const BR_B: f64 = 0.006;
const BR_H: f64 = 0.010;

fn brace_section() -> BeamSection {
    BeamSection {
        young_pa: 13.8e9,
        shear_pa: 0.8e9,
        density_kg_m3: 460.0,
        area_m2: BR_B * BR_H,
        iy_m4: BR_B * BR_H.powi(3) / 12.0,
        iz_m4: BR_H * BR_B.powi(3) / 12.0,
        // Saint-Venant constant of a 10 x 6 mm rectangle, beta a b^3 with beta = 0.21 (a/b = 1.67).
        torsion_m4: 0.21 * BR_H * BR_B.powi(3),
    }
}

fn brace_timoshenko() -> BeamKind {
    BeamKind::Timoshenko(ShearAreas {
        y_m2: 5.0 / 6.0 * BR_B * BR_H,
        z_m2: 5.0 / 6.0 * BR_B * BR_H,
    })
}

fn budget() -> TieBudget {
    TieBudget {
        max_segments: 1_000_000,
        max_slave_dofs: 1_000_000,
        gap_tolerance_m: 1e-3,
        angle_tolerance_rad: 0.3,
    }
}

fn mid_line_y(x: f64, bump: f64) -> f64 {
    0.5 * PL_W + bump * (std::f64::consts::PI * x / PL_L).sin()
}

struct Plate {
    mesh: ShellMesh,
    clamped: Vec<usize>,
    line: Vec<usize>,
    /// The grid line 12.5 mm off the mid line, one element of the coarsest rung: the planted
    /// negative's wrong station.
    wrong_line: Vec<usize>,
}

/// 2n x n structured cells; every grid line j = const follows the bump, so the mid line is a mesh line.
fn plate(n: usize, bump: f64) -> Plate {
    let (nx, ny) = (2 * n, n);
    let id = |i: usize, j: usize| j * (nx + 1) + i;
    let mut nodes = Vec::new();
    for j in 0..=ny {
        for i in 0..=nx {
            let x = PL_L * i as f64 / nx as f64;
            let y = PL_W * j as f64 / ny as f64 + bump * (std::f64::consts::PI * x / PL_L).sin();
            nodes.push([x, y, 0.0]);
        }
    }
    let mut tris = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
            if (i + j) % 2 == 0 {
                tris.push([a, b, c]);
                tris.push([a, c, d]);
            } else {
                tris.push([a, b, d]);
                tris.push([b, c, d]);
            }
        }
    }
    Plate {
        mesh: ShellMesh::new(nodes, tris).unwrap(),
        clamped: (0..=ny).map(|j| id(0, j)).collect(),
        line: (0..=nx).map(|i| id(i, ny / 2)).collect(),
        wrong_line: (0..=nx).map(|i| id(i, ny / 2 + ny / 16)).collect(),
    }
}

fn dkt_sections(p: &Plate) -> Vec<PlateSection> {
    vec![PlateSection::isotropic(PL_E, PL_NU, PL_H, PL_RHO).unwrap(); p.mesh.tris.len()]
}

/// The plate as CST/DKT (Kirchhoff), without beams.
fn dkt_plate(p: &Plate) -> ShellModel {
    assemble_shell_sections(&p.mesh, &dkt_sections(p), &p.clamped, ShellSupport::Clamped).unwrap()
}

/// The plate as MITC3+ (Reissner-Mindlin), clamped along x = 0 in global frames.
fn mitc_plate(p: &Plate) -> ShellModel {
    let g = PL_E / (2.0 * (1.0 + PL_NU));
    let section = ShearShellSection::from_orthotropic_solid(
        PL_E,
        PL_E,
        PL_NU,
        g,
        g,
        g,
        PL_H,
        PL_RHO,
        0.0,
        5.0 / 6.0,
    )
    .unwrap();
    let directors =
        NodalDirectors::new(&p.mesh, vec![[0.0, 0.0, 1.0]; p.mesh.nodes.len()], 0.9).unwrap();
    let constraints: Vec<ShellConstraint> = p
        .clamped
        .iter()
        .map(|&node| ShellConstraint {
            node,
            mask: [true; 6],
            frame: ConstraintFrame::Global,
        })
        .collect();
    assemble_mitc3p(
        &p.mesh,
        &directors,
        &vec![section; p.mesh.tris.len()],
        &constraints,
        ShellMass::Lumped,
        MitcBudget {
            max_nodes: p.mesh.nodes.len(),
            max_triangles: p.mesh.tris.len(),
        },
        1e-3,
        || false,
    )
    .unwrap()
    .pencil()
    .clone()
}

fn centroid_arm() -> f64 {
    0.5 * PL_H + 0.5 * BR_H
}

/// The bonded eccentric reference on the CST/DKT plate: one `ShellBeam` per shell edge on the mid line.
fn on_node_shell_beam(p: &Plate) -> (fs_sparse::Csr, fs_sparse::Csr) {
    let model = on_node_model(p);
    (model.k, model.m)
}

/// The same bonded eccentric model, whole.
fn on_node_model(p: &Plate) -> ShellModel {
    on_node_model_with_mass(p, BeamMass::Consistent)
}

fn on_node_model_with_mass(p: &Plate, mass: BeamMass) -> ShellModel {
    let beams: Vec<ShellBeam> = p
        .line
        .windows(2)
        .map(|w| ShellBeam {
            nodes: [w[0], w[1]],
            offsets_m: [[0.0, 0.0, centroid_arm()]; 2],
            section_z: [0.0, 0.0, 1.0],
            section: brace_section(),
        })
        .collect();
    match mass {
        BeamMass::Lumped => assemble_stiffened_shell(
            &p.mesh, &dkt_sections(p), &p.clamped, ShellSupport::Clamped, &beams,
        ),
        BeamMass::Consistent => assemble_stiffened_shell_with_mass(
            &p.mesh, &dkt_sections(p), &p.clamped, ShellSupport::Clamped, &beams, mass,
        ),
    }
    .unwrap()
}

/// `base` (the plate without beams) tied to a beam chain through `beam_nodes` (centroids) along the mid line.
fn tie_beam(
    p: &Plate,
    base: &ShellModel,
    beam_nodes: &[[f64; 3]],
    kind: BeamKind,
    budget: TieBudget,
) -> Result<fs_plate::shell::junction::BeamTiedShellModel, fs_plate::PlateError> {
    tie_beam_on(p, base, &p.line, beam_nodes, kind, budget)
}

/// The same, with the master polyline given.
fn tie_beam_on(
    p: &Plate,
    base: &ShellModel,
    master_line: &[usize],
    beam_nodes: &[[f64; 3]],
    kind: BeamKind,
    budget: TieBudget,
) -> Result<fs_plate::shell::junction::BeamTiedShellModel, fs_plate::PlateError> {
    let thickness = vec![PL_H; p.mesh.nodes.len()];
    let components = [ShellComponent {
        id: 0,
        mesh: &p.mesh,
        model: base,
        nodal_thickness_m: &thickness,
    }];
    let tied = tie_shells(&components, &[], budget, || false)?;
    let beam = BeamOnShell {
        beam_nodes_m: beam_nodes,
        section_z: [0.0, 0.0, 1.0],
        section: brace_section(),
        kind,
        master: TieSide {
            component: 0,
            nodes: master_line.to_vec(),
        },
        side: 1,
        centroid_height_m: 0.5 * BR_H,
    };
    tie_beams_to_shells(&tied, &components, &[beam], budget, || false)
}

/// Beam centroids on the shell line nodes (the on-node beam).
fn on_line_nodes(p: &Plate) -> Vec<[f64; 3]> {
    p.line
        .iter()
        .map(|&n| {
            let x = p.mesh.nodes[n];
            [x[0], x[1], centroid_arm()]
        })
        .collect()
}

/// Beam centroids on the mid line, `segments` equal steps in x (off the shell nodes).
fn beam_nodes(segments: usize, bump: f64) -> Vec<[f64; 3]> {
    (0..=segments)
        .map(|k| {
            let x = PL_L * k as f64 / segments as f64;
            [x, mid_line_y(x, bump), centroid_arm()]
        })
        .collect()
}

fn lowest_hz(k: &fs_sparse::Csr, m: &fs_sparse::Csr, hi: f64, count: usize) -> Vec<f64> {
    let report = slice_window(k, m, (0.0, hi), &SliceOptions::default()).unwrap();
    assert!(
        report.modes.len() >= count,
        "window holds {} modes",
        report.modes.len()
    );
    report
        .modes
        .iter()
        .take(count)
        .map(|p| p.lambda.sqrt() / (2.0 * std::f64::consts::PI))
        .collect()
}

/// Twice the sixth eigenvalue of a coarse (n = 8) model, from a dense solve.
fn window_from(k: &fs_sparse::Csr, m: &fs_sparse::Csr) -> f64 {
    let n = k.nrows();
    let dense = |a: &fs_sparse::Csr| -> Vec<f64> {
        let mut out = vec![0.0; n * n];
        for r in 0..n {
            let (cols, vals) = a.row(r);
            for (&c, &v) in cols.iter().zip(vals) {
                out[r * n + c] = v;
            }
        }
        out
    };
    let modes = eigh_gen_dense(&dense(k), &dense(m), n).unwrap();
    2.0 * modes[5].lambda
}

#[test]
fn g1_matching_beam_tie_reproduces_on_node_beam() {
    for bump in [0.0, 0.02] {
        let p = plate(8, bump);
        let (k, m) = on_node_shell_beam(&p);
        let hi = window_from(&k, &m);
        let tied = tie_beam(
            &p,
            &dkt_plate(&p),
            &on_line_nodes(&p),
            BeamKind::EulerBernoulli,
            budget(),
        )
        .unwrap();
        let a = lowest_hz(&k, &m, hi, 6);
        let b = lowest_hz(&tied.model.k, &tied.model.m, hi, 6);
        println!(
            "matching nodes, bump {bump}: on-node {a:.4?} Hz, tied {b:.4?} Hz, gap {:.1e} m",
            tied.receipt.max_gap_m
        );
        for (x, y) in a.iter().zip(&b) {
            assert!(
                (y / x - 1.0).abs() < 1e-9,
                "matching tie {y} against on-node {x}"
            );
        }
    }
}

fn dense_mul(a: &[f64; 144], b: &[f64; 144]) -> [f64; 144] {
    let mut out = [0.0; 144];
    for i in 0..12 {
        for k in 0..12 {
            for j in 0..12 {
                out[12 * i + j] += a[12 * i + k] * b[12 * k + j];
            }
        }
    }
    out
}

fn dense_transpose(a: &[f64; 144]) -> [f64; 144] {
    std::array::from_fn(|n| a[12 * (n % 12) + n / 12])
}

fn unit3(v: [f64; 3]) -> [f64; 3] {
    let length = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    v.map(|x| x / length)
}

/// Per-entry magnitude of the summands that both the legacy `ShellBeam` assembly and the on-node tie add,
/// dense on the legacy free coordinates: S = |P| + sum_e |T_e|^T |A_e| |T_e|. P is the plate pencil; for
/// each segment e, A_e = R_e A_e^glob R_e^T is its matrix in its chord axes R_e (x along the chord, z the
/// normal part of the section hint, y = z x x, as the legacy element builds them), and T_e = R_e T_arm maps
/// the shell DOFs of its two nodes to them, T_arm being the rigid arm (u + theta x d, theta) with d the
/// centroid offset. Both paths evaluate exactly this sum of terms, each in its own order.
fn summand_magnitude(
    p: &Plate,
    plate: &fs_sparse::Csr,
    dof_map: &[Option<usize>],
    mass: bool,
) -> Vec<f64> {
    let n = plate.nrows();
    let mut s: Vec<f64> = plate.to_dense().iter().map(|v| v.abs()).collect();
    let arm = centroid_arm();
    for w in p.line.windows(2) {
        let c: [[f64; 3]; 2] = std::array::from_fn(|i| {
            let x = p.mesh.nodes[w[i]];
            [x[0], x[1], x[2] + arm]
        });
        let (k, m) = segment_matrices(
            c[0],
            c[1],
            [0.0, 0.0, 1.0],
            brace_section(),
            BeamKind::EulerBernoulli,
        )
        .unwrap();
        let global = if mass { m } else { k };
        let x = unit3(std::array::from_fn(|i| c[1][i] - c[0][i]));
        let z = unit3(std::array::from_fn(|i| [0.0, 0.0, 1.0][i] - x[2] * x[i]));
        let y = [
            z[1] * x[2] - z[2] * x[1],
            z[2] * x[0] - z[0] * x[2],
            z[0] * x[1] - z[1] * x[0],
        ];
        let axes = [x, y, z];
        let mut r = [0.0; 144];
        for block in 0..4 {
            for a in 0..3 {
                for j in 0..3 {
                    r[12 * (3 * block + a) + 3 * block + j] = axes[a][j];
                }
            }
        }
        let local = dense_mul(&dense_mul(&r, &global), &dense_transpose(&r));
        // T_arm per node: u + theta x d with d = (0, 0, arm), i.e. u_x += arm theta_y, u_y -= arm theta_x.
        let mut t = [0.0; 144];
        for node in 0..2 {
            let o = 6 * node;
            for i in 0..6 {
                t[12 * (o + i) + o + i] = 1.0;
            }
            t[12 * o + o + 4] = arm;
            t[12 * (o + 1) + o + 3] = -arm;
        }
        let te = dense_mul(&r, &t).map(f64::abs);
        let contribution = dense_mul(&dense_mul(&dense_transpose(&te), &local.map(f64::abs)), &te);
        for i in 0..12 {
            let Some(row) = dof_map[6 * w[i / 6] + i % 6] else {
                continue;
            };
            for j in 0..12 {
                let Some(col) = dof_map[6 * w[j / 6] + j % 6] else {
                    continue;
                };
                s[n * row + col] += contribution[12 * i + j];
            }
        }
    }
    s
}

/// The tie's transfer from the reduced shell coordinates to the beam DOFs (rows, global axes), read
/// column by column through `expand_beam`.
fn tie_transfer(tied: &fs_plate::shell::junction::BeamTiedShellModel) -> Vec<Vec<f64>> {
    let n = tied.model.free;
    let columns: Vec<Vec<f64>> = (0..n)
        .map(|q| {
            let mut e = vec![0.0; n];
            e[q] = 1.0;
            tied.expand_beam(0, &e).unwrap()
        })
        .collect();
    (0..columns[0].len())
        .map(|a| columns.iter().map(|c| c[a]).collect())
        .collect()
}

/// The legacy element's exact transfer in global axes: beam node k follows shell node `line[k]` rigidly,
/// u_b = u + theta x d and theta_b = theta, with d = (0, 0, arm).
fn arm_transfer(p: &Plate, dof_map: &[Option<usize>], n: usize) -> Vec<Vec<f64>> {
    let arm = centroid_arm();
    let mut t = vec![vec![0.0; n]; 6 * p.line.len()];
    for (k, &node) in p.line.iter().enumerate() {
        for c in 0..6 {
            if let Some(q) = dof_map[6 * node + c] {
                t[6 * k + c][q] = 1.0;
            }
        }
        if let Some(q) = dof_map[6 * node + 4] {
            t[6 * k][q] = arm;
        }
        if let Some(q) = dof_map[6 * node + 3] {
            t[6 * k + 1][q] = -arm;
        }
    }
    t
}

/// sum_e A_e^T G_e B_e over the beam segments (G_e the segment's global matrix, A_e and B_e the rows of
/// transfers A and B for its two nodes), dense n x n; with `magnitude`, of the entries' absolute values.
fn segment_products(
    p: &Plate,
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    mass: bool,
    n: usize,
    magnitude: bool,
) -> Vec<f64> {
    let v = |x: f64| if magnitude { x.abs() } else { x };
    let arm = centroid_arm();
    let mut out = vec![0.0; n * n];
    for (e, w) in p.line.windows(2).enumerate() {
        let c: [[f64; 3]; 2] = std::array::from_fn(|i| {
            let x = p.mesh.nodes[w[i]];
            [x[0], x[1], x[2] + arm]
        });
        let (k, m) = segment_matrices(
            c[0],
            c[1],
            [0.0, 0.0, 1.0],
            brace_section(),
            BeamKind::EulerBernoulli,
        )
        .unwrap();
        let g = if mass { m } else { k };
        let support = |t: &[Vec<f64>]| -> Vec<usize> {
            (0..n)
                .filter(|&q| (0..12).any(|r| t[6 * e + r][q] != 0.0))
                .collect()
        };
        let (rows_a, columns_b) = (support(a), support(b));
        for &j in &columns_b {
            let x: [f64; 12] = std::array::from_fn(|r| {
                (0..12).map(|s| v(g[12 * r + s]) * v(b[6 * e + s][j])).sum()
            });
            for &i in &rows_a {
                out[n * i + j] += (0..12).map(|r| v(a[6 * e + r][i]) * x[r]).sum::<f64>();
            }
        }
    }
    out
}

/// The plan's original G3 claim, "an on-node beam tie reproduces `assemble_stiffened_shell` K and M bitwise", is
/// FALSIFIED (measured 2026-10-23; ruling 6285): on this fixture 3212 of the 864^2 K entries differ in their
/// bits, by at most 9.6e-17 of the largest entry. The legacy path itself is untouched, so the piano boards
/// are not affected (their pencil is pinned by `g3_legacy_bonded_beam_pencil_golden`). What holds is
/// agreement to rounding, gated by a DERIVED bound rather than a chosen one. Both sides now select
/// consistent beam mass; the separate legacy golden below still calls the frozen piano entry point.
/// The shared rounded local mass integral is identical on both sides, so its quadrature evaluation
/// precedes the comparison's rounding chain. The two paths differ in three
/// ways, and each is accounted for:
/// - Arithmetic order. Both evaluate a sum of the plate entry P and each segment's congruence terms. Along a
///   chain of at most m roundings a computed sum of products is within gamma_m of the sum of its terms'
///   magnitudes (gamma_m = m u / (1 - m u), u = 2^-53; Higham, Accuracy and Stability of Numerical
///   Algorithms, 2nd ed., section 3.1). The legacy terms, in the segment's chord axes R_e, have magnitude
///   S_leg = |P| + sum_e |R_e T_arm|^T |A_e| |R_e T_arm| (`summand_magnitude`); the tie's, in global axes,
///   S_tie = |P| + sum_e |T|^T |G_e| |T|. m = 64 bounds either chain: the legacy element is at most 41 deep
///   (7 for its local entry, 6 for the transform, two 12-term inner products, 2 assembly additions) and the
///   tie at most 29 (8 local, 8 rotation, 11 for the congruence, 2 assembly).
/// - The transfer. The legacy element uses the exact rigid arm T_arm; the tie uses its dual mortar of matching
///   nodes, which equals T_arm only up to the rounding of its own quadrature, T = T_arm + dT. That adds
///   exactly E = T^T G T - T_arm^T G T_arm = dT^T G T + T_arm^T G dT, which reaches entries that the legacy
///   element leaves exactly zero. The test computes E (signed) and subtracts it; computing it is itself
///   within gamma_64 (|dT|^T |G| |T| + |T_arm|^T |G| |dT|).
/// - Symmetry. The tie builds its pencil from the UPPER triangle of each source matrix and mirrors it
///   (`accumulate_congruence`), so it is symmetric bit for bit; the legacy assembly copies the assembled plate
///   pencil P, which is symmetric only to rounding (its asymmetry max |P_ij - P_ji| is printed). The entries are
///   therefore compared on and above the diagonal, where both read the same plate entry, and the tie's own
///   bit symmetry is asserted.
///
/// So for i <= j, |K_tie - K_legacy - E|_ij <= gamma_64 (S_leg + S_tie + |dT|^T |G| |T| + |T_arm|^T |G| |dT|)_ij,
/// and an entry with no term at all must agree exactly. The same holds for M.
#[test]
fn g3_on_node_beam_tie_matches_legacy_to_rounding() {
    let depth = 64.0;
    let u = f64::EPSILON / 2.0;
    let gamma = depth * u / (1.0 - depth * u);
    for bump in [0.0, 0.02] {
        let p = plate(8, bump);
        let legacy = on_node_model(&p);
        let base = dkt_plate(&p);
        let tied = tie_beam(
            &p,
            &base,
            &on_line_nodes(&p),
            BeamKind::EulerBernoulli,
            budget(),
        )
        .unwrap();
        let n = legacy.free;
        let t_tie = tie_transfer(&tied);
        let t_arm = arm_transfer(&p, &legacy.dof_map, n);
        let delta: Vec<Vec<f64>> = t_tie
            .iter()
            .zip(&t_arm)
            .map(|(a, b)| a.iter().zip(b).map(|(x, y)| x - y).collect())
            .collect();
        let delta_max = delta.iter().flatten().fold(0.0_f64, |s, v| s.max(v.abs()));
        let delta_count = delta.iter().flatten().filter(|v| **v != 0.0).count();
        for (name, a, b, plate_part, mass) in [
            ("K", &legacy.k, &tied.model.k, &base.k, false),
            ("M", &legacy.m, &tied.model.m, &base.m, true),
        ] {
            assert_eq!((a.nrows(), a.ncols()), (b.nrows(), b.ncols()));
            let s_leg = summand_magnitude(&p, plate_part, &legacy.dof_map, mass);
            let plate = plate_part.to_dense();
            let plate_abs: Vec<f64> = plate.iter().map(|v| v.abs()).collect();
            let s_tie = segment_products(&p, &t_tie, &t_tie, mass, n, true);
            let d1 = segment_products(&p, &delta, &t_tie, mass, n, true);
            let d2 = segment_products(&p, &t_arm, &delta, mass, n, true);
            let e1 = segment_products(&p, &delta, &t_tie, mass, n, false);
            let e2 = segment_products(&p, &t_arm, &delta, mass, n, false);
            let (x, y) = (a.to_dense(), b.to_dense());
            let (mut bits, mut worst, mut transfer_only) = (0usize, 0.0_f64, 0usize);
            let mut transfer_max = 0.0_f64;
            let mut asymmetry = 0.0_f64;
            for index in 0..n * n {
                let (row, column) = (index / n, index % n);
                let mirror = n * column + row;
                assert!(
                    y[index].to_bits() == y[mirror].to_bits(),
                    "bump {bump}, {name}: the tied pencil is not symmetric at ({row}, {column})"
                );
                asymmetry = asymmetry.max((plate[index] - plate[mirror]).abs());
                if row > column {
                    continue;
                }
                let (l, t) = (x[index], y[index]);
                bits += usize::from(l.to_bits() != t.to_bits());
                let transfer = e1[index] + e2[index];
                transfer_max = transfer_max.max(transfer.abs());
                if s_leg[index] == 0.0 && d1[index] + d2[index] > 0.0 {
                    transfer_only += 1;
                }
                let bound = gamma
                    * (s_leg[index] + plate_abs[index] + s_tie[index] + d1[index] + d2[index]);
                let residual = t - l - transfer;
                if bound == 0.0 {
                    assert!(
                        residual == 0.0,
                        "bump {bump}, {name}: entry ({row}, {column}) has no term but differs: legacy {l:e}, tie {t:e}"
                    );
                } else {
                    worst = worst.max(residual.abs() / bound);
                }
            }
            println!(
                "bump {bump}, {name}: {n} x {n} entries, {bits} on or above the diagonal differ in bits (legacy plate asymmetry {asymmetry:.2e}); the tie's transfer departs from the exact arm in {delta_count} entries, by at most {delta_max:.2e}, and that alone reaches {transfer_only} entries the legacy element leaves zero (|E| up to {transfer_max:.2e}); worst |tie - legacy - E| is {worst:.3e} of the derived bound"
            );
            assert!(
                worst <= 1.0,
                "bump {bump}, {name}: the on-node tie departs from the legacy pencil beyond its derived bound ({worst:.3e})"
            );
        }
    }
}

const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

fn fold(acc: u64, v: u64) -> u64 {
    v.to_le_bytes()
        .iter()
        .fold(acc, |a, &b| (a ^ u64::from(b)).wrapping_mul(FNV_PRIME))
}

/// A pencil matrix's row structure and value bits.
fn fold_matrix(mut acc: u64, a: &fs_sparse::Csr) -> u64 {
    acc = fold(acc, a.nrows() as u64);
    for r in 0..a.nrows() {
        let (columns, values) = a.row(r);
        acc = fold(acc, columns.len() as u64);
        for (&c, &v) in columns.iter().zip(values) {
            acc = fold(fold(acc, c as u64), v.to_bits());
        }
    }
    acc
}

/// Golden (docs/GOLDEN_POLICY.md; `golden-couplings.json` row `fs-plate:legacy-bonded-beam-pencil`, surface
/// `fs-plate:stiffened-shell-bits`): the bits of the legacy bonded `ShellBeam` pencil of
/// `assemble_stiffened_shell`, which the piano boards consume (`crowned_board.rs`, `downbearing.rs`). Frozen
/// on x86_64 Linux in debug and release alike. It is not a cross-ISA claim: the curved fixture's
/// coordinates use the platform `sin`.
const LEGACY_BONDED_PENCIL_GOLDEN: u64 = 0x24de_da52_df52_7693;

#[test]
fn g3_legacy_bonded_beam_pencil_golden() {
    let mut acc = FNV_OFFSET;
    for bump in [0.0, 0.02] {
        let model = on_node_model_with_mass(&plate(8, bump), BeamMass::Lumped);
        acc = fold_matrix(fold_matrix(acc, &model.k), &model.m);
    }
    println!(
        "legacy bonded-beam pencil digest {acc:#018x} ({} build)",
        if cfg!(debug_assertions) {
            "debug"
        } else {
            "release"
        }
    );
    assert_eq!(
        acc, LEGACY_BONDED_PENCIL_GOLDEN,
        "the legacy bonded-beam pencil's bits moved: {acc:#018x} vs {LEGACY_BONDED_PENCIL_GOLDEN:#018x}. This is a golden event \
         (docs/GOLDEN_POLICY.md): name the semantic change, bump STIFFENED_SHELL_BIT_SEMANTICS_VERSION and re-freeze \
         in debug and release; never re-pin silently. The piano boards consume this pencil."
    );
}

/// Local orders between consecutive rungs (each rung halves h).
fn local_orders(gaps: &[f64]) -> Vec<f64> {
    gaps.windows(2).map(|w| (w[0] / w[1]).log2()).collect()
}

/// The one-sided gate of the orchestrator's ruling (mail 6174): the gap falls strictly at every rung, and
/// its fitted order over the three finest rungs is at least 1.8. Returns the fit, or why it fails.
fn falls_at_order_two(rungs: &[usize], gaps: &[f64]) -> Result<f64, String> {
    if gaps.windows(2).any(|w| w[1] >= w[0]) {
        return Err(format!("the gap does not fall at every rung: {gaps:?}"));
    }
    let hs = rungs[1..].iter().map(|&n| 1.0 / n as f64).collect();
    let fit = fit_order(&RefinementLadder::new(hs, gaps[1..].to_vec()).map_err(|e| e.to_string())?);
    if fit.observed < 1.8 {
        return Err(format!("fitted order {:.3} < 1.8", fit.observed));
    }
    Ok(fit.observed)
}

#[test]
fn g1_off_node_brace_converges() {
    const MODES: usize = 4;
    let rungs = [16usize, 32, 64, 128];
    for (name, bump) in [("straight", 0.0), ("curved", 0.02)] {
        // The window: twice the sixth eigenvalue of the coarse CST/DKT on-node model (a dense solve);
        // a 3 mm plate has nearly the same low modes in both elements.
        let (k8, m8) = on_node_shell_beam(&plate(8, bump));
        let hi = window_from(&k8, &m8);
        let mut on = Vec::new();
        let mut off = Vec::new();
        let mut wrong = Vec::new();
        for &n in &rungs {
            let p = plate(n, bump);
            let base = mitc_plate(&p);
            let on_tie =
                tie_beam(&p, &base, &on_line_nodes(&p), brace_timoshenko(), budget()).unwrap();
            on.push(lowest_hz(&on_tie.model.k, &on_tie.model.m, hi, MODES));
            // The beam's own nodes: 3 per 2 shell edges, so only the ends coincide with shell nodes.
            let nodes = beam_nodes(3 * n, bump);
            let off_tie = tie_beam(&p, &base, &nodes, brace_timoshenko(), budget()).unwrap();
            off.push(lowest_hz(&off_tie.model.k, &off_tie.model.m, hi, MODES));
            println!(
                "{name} n={n}: on-node {:.4?} Hz, off-node {:.4?} Hz, beam gap {:.2e} m",
                on.last().unwrap(),
                off.last().unwrap(),
                off_tie.receipt.max_gap_m
            );
            if bump == 0.0 {
                // Planted negative (ruling 6174): the same beam tied to the wrong station, the grid line one
                // coarsest-rung element (12.5 mm) off the mid line at every rung. The distance is fixed in
                // space, so the error does not refine away and the gap should stop falling.
                let wide = TieBudget {
                    gap_tolerance_m: 0.02,
                    ..budget()
                };
                let w = tie_beam_on(&p, &base, &p.wrong_line, &nodes, brace_timoshenko(), wide)
                    .unwrap();
                wrong.push(lowest_hz(&w.model.k, &w.model.m, hi, MODES));
            }
        }
        // The quantity: the off-node beam's distance to the on-node beam at the same rung. The plate's own
        // discretization error is common to both and cancels; a constant gap (a wrong limit) would stop
        // falling. Modes within 1 % of each other on the finest rung are gated as their cluster mean (a near
        // pair trades character between rungs); a cluster whose gap at the coarsest rung is below 1e-5 carries
        // almost no brace energy and sits at the eigensolver's accuracy, so it is printed and not gated.
        let finest = on.last().unwrap();
        let mut clusters: Vec<Vec<usize>> = Vec::new();
        for mode in 0..MODES {
            let joins = clusters
                .last()
                .is_some_and(|c| finest[mode] / finest[*c.last().unwrap()] - 1.0 < 0.01);
            if joins {
                clusters.last_mut().unwrap().push(mode);
            } else {
                clusters.push(vec![mode]);
            }
        }
        let mean = |r: &[f64], c: &[usize]| c.iter().map(|&i| r[i]).sum::<f64>() / c.len() as f64;
        let gaps_of = |series: &[Vec<f64>], c: &[usize]| -> Vec<f64> {
            (0..rungs.len())
                .map(|i| (mean(&series[i], c) / mean(&on[i], c) - 1.0).abs())
                .collect()
        };
        let mut gated = 0;
        for cluster in &clusters {
            let gaps = gaps_of(&off, cluster);
            let text: Vec<String> = gaps.iter().map(|g| format!("{g:.2e}")).collect();
            let orders: Vec<String> = local_orders(&gaps)
                .iter()
                .map(|p| format!("{p:.2}"))
                .collect();
            if gaps[0] < 1e-5 {
                println!(
                    "{name} modes {cluster:?}: relative gap {text:?}, local orders {orders:?}; {}",
                    "below 1e-5 at the coarsest rung: no brace energy, not gated"
                );
                continue;
            }
            let fit = falls_at_order_two(&rungs, &gaps)
                .unwrap_or_else(|e| panic!("{name} brace modes {cluster:?}: {e}"));
            println!(
                "{name} modes {cluster:?}: relative gap {text:?}, local orders {orders:?}, {}",
                format!("fitted order over n = 32..128 {fit:.3} (>= 1.8)")
            );
            gated += 1;
            if bump == 0.0 {
                let planted = gaps_of(&wrong, cluster);
                let planted_text: Vec<String> =
                    planted.iter().map(|g| format!("{g:.2e}")).collect();
                let verdict = falls_at_order_two(&rungs, &planted);
                println!(
                    "{name} modes {cluster:?}, planted wrong station: relative gap {planted_text:?}, gate {verdict:?}"
                );
                assert!(
                    verdict.is_err(),
                    "the gate must catch a beam tied 12.5 mm off its station"
                );
            }
        }
        assert!(gated > 0, "{name}: no brace-carried mode to gate");
    }
}

#[test]
fn g0_beam_path_off_surface_refuses() {
    let p = plate(8, 0.0);
    let base = mitc_plate(&p);
    let mut nodes = beam_nodes(24, 0.0);
    assert!(tie_beam(&p, &base, &nodes, brace_timoshenko(), budget()).is_ok());
    // 5 mm above where a beam standing on the skin can be: beyond the 1 mm gap tolerance.
    for node in &mut nodes {
        node[2] += 0.005;
    }
    assert!(tie_beam(&p, &base, &nodes, brace_timoshenko(), budget()).is_err());
    // A one-node beam refuses too.
    assert!(tie_beam(&p, &base, &nodes[..1], brace_timoshenko(), budget()).is_err());
}

#[test]
fn g0_beam_tie_budget_cap_and_cap_plus_one() {
    let p = plate(8, 0.0);
    let base = mitc_plate(&p);
    let nodes = beam_nodes(24, 0.0);
    let full = tie_beam(&p, &base, &nodes, brace_timoshenko(), budget()).unwrap();
    let dofs = 6 * nodes.len();
    assert_eq!(full.receipt.slave_dofs, dofs);
    let segments = full.receipt.segments;
    let at = |max_segments: usize, max_slave_dofs: usize| {
        tie_beam(
            &p,
            &base,
            &nodes,
            brace_timoshenko(),
            TieBudget {
                max_segments,
                max_slave_dofs,
                ..budget()
            },
        )
    };
    assert!(at(segments, dofs).is_ok());
    assert!(at(segments - 1, dofs).is_err());
    assert!(at(segments, dofs - 1).is_err());
}
