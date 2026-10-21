//! Bar-S1 (.266), G0/G1/G2: existing TET10 solid (P1 kept as a reported
//! locking diagnostic), certified modal slicer, Burton--Miller radiation and
//! material Hertz owners composed as fixture data.
//!
//! This is module evidence, never a reconstruction of the Iowa specimen.
//! Geometry and longitudinal material data reuse `xylophone_bar.rs`. Unknown
//! rosewood transverse constants are an explicitly bounded, correlated set of
//! three hardwood-ratio scenarios, NOT measurements of rosewood or probability
//! draws. Material loss, supports, resonator, acoustic back-loading and the
//! mallet shaft are absent. The lossless render is a diagnostic scenario;
//! neither its decay nor its modal/narrow-band truncation is qualified.
//! Native execution and grading belong to the central Linux lane.

use fs_bem::helmholtz::{Formulation, Medium, exterior_pressure_at_points, solve_radiation_batch};
use fs_bem::panel3d::SpherePanels;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel, ModalAcousticWorkspace,
};
use fs_exec::Cx;
use fs_material::OrthotropicElastic;
use fs_math::{c64::C64, det};
use fs_modal::{ModePair, SliceOptions, SliceReport, eigh_gen_dense, slice_window};
use fs_solid::linear3::tet10::{
    BoundaryProjection, QuadraticTetProblem, TetQuadrature, promote_p1_to_p2,
};
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticError, TetElasticMaterial,
    TetLinearElasticProblem, TetMaterialField,
};
use std::path::{Path, PathBuf};
use std::sync::{Mutex, OnceLock};
use std::time::Instant;

const L: f64 = 0.38;
const W: f64 = 0.038;
const T: f64 = 0.018;
const RHO: f64 = 750.0; // Authored Estimate, NOT SG-derived 12% MC density.
const TAU: f64 = core::f64::consts::TAU;
const BEAM_TOL: f64 = 0.01;
const ET_DB5_HZ: f64 = 554.365_261_953_744_2;
const RECEIVER: [f64; 3] = [0.114, 0.20, 1.0];
const RATE: u32 = 48_000;
const FULL_SCALE_PA: f64 = 0.1;
// Four rungs allow two independently observed orders. No historical rung is
// substituted. Long native qualification: the last has 147456 P1 tetrahedra.
const RUNGS: [usize; 4] = [1, 2, 4, 8];
// DEF-08: P1 tetrahedra lock in bending, so the beam cross-check is gated on
// TET10 (M1.3-W1a). The finest rung has 18432 quadratic tetrahedra; rung 3 is
// not nested in rung 2, so orders use the actual h ratios.
const TET10_RUNGS: [usize; 4] = [1, 2, 3, 4];
// A mode whose whole bound |p| + u stays below this fraction of the strongest
// mode's receiver pressure (per unit modal velocity, i.e. per unit kinetic
// energy) is a near-null radiator: relative error is undefined as |p| -> 0,
// so it is BOUNDED against the dominant scale, never claimed resolved. The
// prism's torsion-1 is a near-cancelled quadrupole: |p| falls 3x per rung.
const NEAR_NULL_FRACTION: f64 = 0.1;
// Formal order 1: separation is detected on the substep grid and the mallet
// uses explicit symplectic Euler; the modal body step itself is exact.
const SUBSTEPS: [u32; 4] = [2, 4, 8, 16];
const FMAX: f64 = 4_500.0;
const KINDS: [&str; 5] = [
    "vertical-1",
    "vertical-2",
    "vertical-3",
    "torsion-1",
    "lateral-1",
];
static SERIAL: Mutex<()> = Mutex::new(());
static UNIFORM: OnceLock<Ladder> = OnceLock::new();
static UNIFORM_TET10: OnceLock<Ladder> = OnceLock::new();
static UNDERCUT: OnceLock<Ladder> = OnceLock::new();
static RADIATION: OnceLock<RadiationLadder> = OnceLock::new();

const ROSEWOOD: &str =
    include_str!("../../../data/matdb/seed-v1/rosewood-indian-fpl-gtr282/properties.tsv");
const HONDURAS: &str =
    include_str!("../../../data/matdb/seed-v1/mahogany-honduras-fpl-gtr282/properties.tsv");
const AFRICAN: &str =
    include_str!("../../../data/matdb/seed-v1/mahogany-african-fpl-gtr282/properties.tsv");
const WALNUT: &str =
    include_str!("../../../data/matdb/seed-v1/walnut-black-fpl-gtr282/properties.tsv");

/// Measured rosewood-bar loss constants, one `scalar` row each (Chaigne and
/// Doutaut 1997, JASA 101(1) 539, Table I and Eqs. 18-19).
const ROSEWOOD_LOSS: &str = include_str!(
    "../../../data/instruments/descriptions/bars/rosewood-loss-chaigne-doutaut1997.tsv"
);

fn scalar(pack: &str, property: &str) -> f64 {
    let row = pack.lines().find_map(|line| {
        let fields: Vec<_> = line.split('\t').collect();
        (fields.first() == Some(&"scalar") && fields.get(3) == Some(&property))
            .then(|| fields[4].parse::<f64>().expect("source scalar"))
    });
    row.unwrap_or_else(|| panic!("required source property missing: {property}"))
}

/// Chaigne and Doutaut measured the free decay of rosewood-bar partials with a
/// laser vibrometer and fitted the amplitude decay rate 1/tau_d = a0 + a2 f^2
/// (Eq. 18), with a0 = gamma_B / 2 and a2 = 2 pi^2 eta (Eq. 19). For a partial
/// of angular frequency w the equivalent modal damping ratio is
/// zeta = (a0 + a2 f^2) / w. The constants belong to THEIR bar (E 21.3 GPa,
/// 1015 kg/m3), so applying them to this one is a declared scenario.
fn modal_damping_ratio(f_hz: f64) -> f64 {
    let pi = core::f64::consts::PI;
    let a0 = scalar(ROSEWOOD_LOSS, "fluid_damping_gamma_b") / 2.0;
    let a2 = 2.0 * pi * pi * scalar(ROSEWOOD_LOSS, "viscoelastic_constant_eta");
    (a0 + a2 * f_hz * f_hz) / (TAU * f_hz)
}

fn young() -> f64 {
    scalar(ROSEWOOD, "young_modulus_longitudinal") * 1e6
}

fn law(pack: &str) -> OrthotropicElastic {
    let e = young();
    OrthotropicElastic::new(
        [
            e,
            e * scalar(pack, "er_over_el"),
            e * scalar(pack, "et_over_el"),
        ],
        [
            scalar(pack, "nu_lr"),
            scalar(pack, "nu_lt"),
            scalar(pack, "nu_rt"),
        ],
        [
            e * scalar(pack, "glr_over_el"),
            e * scalar(pack, "grt_over_el"),
            e * scalar(pack, "glt_over_el"),
        ],
        1e-3, // Declared small-strain domain, not a measured failure strain.
    )
    .expect("complete source-ratio scenario must be SPD")
}

fn material(pack: &str) -> TetElasticMaterial {
    let source = format!("{ROSEWOOD}\n{pack}\ndensity-authored-estimate={RHO}");
    let identity = fs_blake3::hash_domain(
        "org.frankensim.bar-s1.material-scenario.v1",
        source.as_bytes(),
    );
    TetElasticMaterial::try_new_oriented_orthotropic(
        RHO,
        &law(pack),
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        identity.0,
        fs_blake3::hash_domain("org.frankensim.bar-s1.axes.LRT.v1", b"x=L,y=R,z=T").0,
    )
    .expect("oriented source scenario")
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 266,
                kernel_id: 266,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// Same record shape for each numerical check; detail is a JSON object made
/// only from finite measured values, fixture data and explicit refusal codes.
fn row(case: &str, axis: &str, verdict: &str, detail: &str, started: Instant) {
    println!(
        "{{\"suite\":\"bar-s1\",\"case\":\"{case}\",\"axis\":\"{axis}\",\"fixture\":\"rosewood-prism-plus-declared-undercut\",\"source\":\"FPL-GTR-282 ch5; xylophone_bar.rs\",\"verdict\":\"{verdict}\",\"elapsed_s\":{},\"detail\":{detail}}}",
        started.elapsed().as_secs_f64()
    );
}

#[derive(Clone, Copy, Debug)]
struct Design {
    depth: f64,
    scale: f64,
}

impl Design {
    const PRISM: Self = Self {
        depth: 0.0,
        scale: 1.0,
    };

    /// A smooth, declared maker's cut on the lower face, over the middle 60%.
    /// Depth is solved on the 3-D pencil; no recording enters construction.
    fn point(self, u: [f64; 3]) -> [f64; 3] {
        let a = ((u[0] - 0.5) / 0.3).abs();
        let cut = if a < 1.0 {
            self.depth * (1.0 - a * a).powi(2)
        } else {
            0.0
        };
        [
            L * u[0] * self.scale,
            W * (u[1] - 0.5) * self.scale,
            T * (-0.5 + cut + u[2] * (1.0 - cut)) * self.scale,
        ]
    }
}

struct Mesh {
    nodes: Vec<[f64; 3]>,
    tets: Vec<[usize; 4]>,
    counts: [usize; 3],
}

fn mesh(counts: [usize; 3], design: Design) -> Mesh {
    let (complex, nodes) = fs_conduction::fixtures::structured_tets(counts, &|i, j, k| {
        design.point([
            i as f64 / counts[0] as f64,
            j as f64 / counts[1] as f64,
            k as f64 / counts[2] as f64,
        ])
    });
    Mesh {
        nodes,
        tets: complex.tets.iter().map(|t| t.map(|i| i as usize)).collect(),
        counts,
    }
}

fn assemble(mesh: &Mesh, mat: &TetElasticMaterial, cx: &Cx<'_>) -> TetElasticAssembly {
    TetLinearElasticProblem {
        nodes_m: &mesh.nodes,
        tetrahedra: &mesh.tets,
        materials: TetMaterialField::Uniform(mat),
        fixed_dofs: &[],
        budget: TetAssemblyBudget::standard(),
    }
    .assemble(cx)
    .expect("free-free P1 assembly")
}

fn options() -> SliceOptions {
    SliceOptions {
        max_lanczos: 2048,
        max_restarts: 48,
        ritz_tol: 1e-13,
        mass_diagonal_equilibration: true,
        ..SliceOptions::default()
    }
}

fn slice(a: &TetElasticAssembly, low_hz: f64, high_hz: f64) -> SliceReport {
    slice_window(
        &a.stiffness,
        &a.mass,
        ((TAU * low_hz).powi(2), (TAU * high_hz).powi(2)),
        &options(),
    )
    .expect("certified complete frequency window")
}

fn hz(m: &ModePair) -> f64 {
    m.lambda.max(0.0).sqrt() / TAU
}

fn hz_error(m: &ModePair) -> f64 {
    let f = hz(m);
    ((m.interval.1.max(0.0).sqrt() / TAU - f).abs())
        .max((f - m.interval.0.max(0.0).sqrt() / TAU).abs())
}

/// Cross-section translation distinguishes the two bending families; torsion
/// is checked independently by its normalized overlap with a rigid twist of
/// the end section. Frequencies alone never assign the mode labels.
fn classify(mesh: &Mesh, modes: &[ModePair]) -> [usize; 5] {
    let [_, ny, nz] = mesh.counts;
    let end = (ny + 1) * (nz + 1);
    let mut vertical = Vec::new();
    let mut lateral = Vec::new();
    let mut torsion = Vec::new();
    for (i, m) in modes.iter().enumerate() {
        let mut mean = [0.0; 3];
        let mut rms = 0.0;
        let mut twist = 0.0;
        let mut twist_norm = 0.0;
        let mut weight = 0.0;
        for node in 0..end {
            let j = node / (nz + 1);
            let k = node % (nz + 1);
            let w = (if j == 0 || j == ny { 0.5 } else { 1.0 })
                * (if k == 0 || k == nz { 0.5 } else { 1.0 });
            let u: [f64; 3] = core::array::from_fn(|c| m.phi[3 * node + c]);
            for c in 0..3 {
                mean[c] += w * u[c];
            }
            let p = mesh.nodes[node];
            twist += w * (-p[2] * u[1] + p[1] * u[2]);
            twist_norm += w * (p[1] * p[1] + p[2] * p[2]);
            rms += w * u.iter().map(|v| v * v).sum::<f64>();
            weight += w;
        }
        mean = mean.map(|v| v / weight);
        let amplitude = (rms / weight).sqrt();
        if mean[2].abs() > 0.3 * amplitude && mean[2].abs() > 2.0 * mean[1].abs() {
            vertical.push(i);
        } else if mean[1].abs() > 0.3 * amplitude && mean[1].abs() > 2.0 * mean[2].abs() {
            lateral.push(i);
        } else if twist.abs() > 0.6 * (twist_norm * rms).sqrt() {
            torsion.push(i);
        }
    }
    assert!(
        vertical.len() >= 3 && !lateral.is_empty() && !torsion.is_empty(),
        "window lacks shape-identified families: vertical={vertical:?}, lateral={lateral:?}, torsion={torsion:?}"
    );
    [
        vertical[0],
        vertical[1],
        vertical[2],
        torsion[0],
        lateral[0],
    ]
}

struct Rung {
    h: f64,
    n_tets: usize,
    frequencies: [f64; 5],
    eigen_errors: [f64; 5],
}

struct Ladder {
    design: Design,
    rungs: Vec<Rung>,
    mesh: Mesh,
    modes: Vec<ModePair>,
    selected: [usize; 5],
    rigid_projectors: Vec<f64>,
}

fn strike_node(mesh: &Mesh) -> usize {
    let [nx, ny, nz] = mesh.counts;
    assert_eq!(nx % 4, 0);
    assert_eq!(ny % 2, 0);
    ((nx / 4) * (ny + 1) + ny / 2) * (nz + 1) + nz
}

/// Canonical phase for reproducible pressure comparisons and the render.
fn canonical_phase(modes: &mut [ModePair]) {
    for m in modes {
        let at = m
            .phi
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.abs().total_cmp(&b.1.abs()))
            .unwrap()
            .0;
        if m.phi[at] < 0.0 {
            for v in &mut m.phi {
                *v = -*v;
            }
        }
    }
}

/// Preserve all six free-body coordinates in contact. Dropping them would
/// spuriously anchor a nominally free-free bar. The analytic affine fields
/// are M-orthonormalized on the actual assembled body, including the cut.
/// `nodes` are the assembly's nodes (the TET10 mid-edge nodes included).
fn rigid_projectors(mesh: &Mesh, nodes: &[[f64; 3]], a: &TetElasticAssembly) -> Vec<f64> {
    let mut basis: Vec<Vec<f64>> = Vec::new();
    let mut mass_basis: Vec<Vec<f64>> = Vec::new();
    let mut scratch = vec![0.0; a.free_dofs.len()];
    let mut projections = Vec::new();
    for axis in 0..6 {
        let mut field: Vec<f64> = nodes
            .iter()
            .flat_map(|&[x, y, z]| match axis {
                0 => [1.0, 0.0, 0.0],
                1 => [0.0, 1.0, 0.0],
                2 => [0.0, 0.0, 1.0],
                3 => [0.0, -z, y],
                4 => [z, 0.0, -x],
                _ => [-y, x, 0.0],
            })
            .collect();
        for _ in 0..2 {
            for (b, mb) in basis.iter().zip(&mass_basis) {
                let dot: f64 = field.iter().zip(mb).map(|(x, y)| x * y).sum();
                for (x, y) in field.iter_mut().zip(b) {
                    *x -= dot * y;
                }
            }
        }
        a.mass.spmv(&field, &mut scratch);
        let norm = field
            .iter()
            .zip(&scratch)
            .map(|(x, y)| x * y)
            .sum::<f64>()
            .sqrt();
        assert!(norm.is_finite() && norm > 0.0);
        for x in &mut field {
            *x /= norm;
        }
        a.mass.spmv(&field, &mut scratch);
        let check: f64 = field.iter().zip(&scratch).map(|(x, y)| x * y).sum();
        assert!((check - 1.0).abs() < 1e-10);
        projections.push(-field[3 * strike_node(mesh) + 2]);
        basis.push(field);
        mass_basis.push(scratch.clone());
    }
    projections
}

fn calculate_ladder(design: Design, bar: &str) -> Ladder {
    with_cx(|cx| {
        let mat = material(HONDURAS);
        let mut rows = Vec::new();
        let mut last = None;
        for r in RUNGS {
            let started = Instant::now();
            cx.checkpoint().expect("cancel before rung");
            let mesh = mesh([24 * r, 2 * r, r], design);
            let a = assemble(&mesh, &mat, cx);
            let mut report = slice(&a, 20.0, FMAX / design.scale);
            canonical_phase(&mut report.modes);
            let selected = classify(&mesh, &report.modes);
            let frequencies = selected.map(|i| hz(&report.modes[i]));
            let eigen_errors = selected.map(|i| hz_error(&report.modes[i]));
            let h = (L / (24 * r) as f64)
                .max(W / (2 * r) as f64)
                .max(T / r as f64)
                * design.scale;
            row(
                "structural-rung",
                "h",
                "pass",
                &format!(
                    "{{\"bar_id\":\"{bar}\",\"rung_h_m\":{h},\"n_tets\":{},\"free_dofs\":{},\"f_modes_hz\":{frequencies:?},\"eigen_uncertainty_hz\":{eigen_errors:?},\"factor_peak_bytes\":{},\"factor_nnz_l\":{},\"window_hz\":[20,{}],\"inertia_count\":{}}}",
                    mesh.tets.len(),
                    a.free_dofs.len(),
                    report.stats.factor_peak_bytes,
                    report.stats.factor_nnz_l,
                    FMAX / design.scale,
                    report.expected
                ),
                started,
            );
            rows.push(Rung {
                h,
                n_tets: mesh.tets.len(),
                frequencies,
                eigen_errors,
            });
            let rigid = if r == RUNGS[3] {
                rigid_projectors(&mesh, &mesh.nodes, &a)
            } else {
                Vec::new()
            };
            last = Some((mesh, report.modes, selected, rigid));
        }
        let (mesh, modes, selected, rigid_projectors) = last.unwrap();
        Ladder {
            design,
            rungs: rows,
            mesh,
            modes,
            selected,
            rigid_projectors,
        }
    })
}

fn uniform() -> &'static Ladder {
    UNIFORM.get_or_init(|| calculate_ladder(Design::PRISM, "A-uniform"))
}

/// Boundary mid-edge nodes land on the exact declared surface. The undercut
/// face is curved, so straight midpoints would leave a P1 geometry error
/// under the quadratic field; on the prism this is the identity.
struct DesignSurface(Design);

impl Design {
    fn reference(self, p: [f64; 3]) -> [f64; 3] {
        let u0 = p[0] / (L * self.scale);
        let a = ((u0 - 0.5) / 0.3).abs();
        let cut = if a < 1.0 {
            self.depth * (1.0 - a * a).powi(2)
        } else {
            0.0
        };
        [
            u0,
            p[1] / (W * self.scale) + 0.5,
            (p[2] / (T * self.scale) + 0.5 - cut) / (1.0 - cut),
        ]
    }
}

impl BoundaryProjection for DesignSurface {
    /// x and y are linear in the reference coordinates, so a boundary-edge
    /// midpoint keeps the exact reference average along the face; only its
    /// face coordinate moves off the face (by the cut's sag, ~1e-3 in u),
    /// and is snapped back to the nearest face of the reference box.
    fn project(&self, p: [f64; 3]) -> Result<[f64; 3], TetElasticError> {
        let mut u = self.0.reference(p);
        let axis = (0..3)
            .min_by(|&a, &b| u[a].min(1.0 - u[a]).total_cmp(&u[b].min(1.0 - u[b])))
            .unwrap();
        u[axis] = if u[axis] < 0.5 { 0.0 } else { 1.0 };
        Ok(self.0.point(u))
    }
}

/// The P1 mesh promoted to TET10 on the existing quadratic owner. Promotion
/// keeps the P1 node IDs, so `classify`, the strike node and the BEM vertex
/// map read the same vertex DOFs as on P1.
fn tet10_assemble(
    mesh: &Mesh,
    design: Design,
    mat: &TetElasticMaterial,
    cx: &Cx<'_>,
) -> (Vec<[f64; 3]>, usize, TetElasticAssembly) {
    let surface = DesignSurface(design);
    let (nodes, elements) = promote_p1_to_p2(&mesh.nodes, &mesh.tets, Some(&surface), cx)
        .expect("P1 to TET10 promotion");
    let a = QuadraticTetProblem {
        nodes_m: &nodes,
        elements: &elements,
        materials: TetMaterialField::Uniform(mat),
        fixed_dofs: &[],
        budget: TetAssemblyBudget::standard(),
        // Affine minima: stiffness degree 2, consistent mass degree 4.
        quadrature: TetQuadrature {
            stiffness_points_per_dir: 3,
            mass_points_per_dir: 4,
        },
    }
    .assemble(cx)
    .expect("free-free TET10 assembly");
    (nodes, elements.len(), a)
}

/// TET10 counterpart of [`calculate_ladder`] on rungs [`TET10_RUNGS`].
fn calculate_tet10_ladder(design: Design, bar: &str) -> Ladder {
    with_cx(|cx| {
        let mat = material(HONDURAS);
        let mut rows = Vec::new();
        let mut last = None;
        for r in TET10_RUNGS {
            let started = Instant::now();
            cx.checkpoint().expect("cancel before rung");
            let mesh = mesh([24 * r, 2 * r, r], design);
            let (nodes, n_tets, a) = tet10_assemble(&mesh, design, &mat, cx);
            let mut report = slice(&a, 20.0, FMAX / design.scale);
            canonical_phase(&mut report.modes);
            let selected = classify(&mesh, &report.modes);
            let frequencies = selected.map(|i| hz(&report.modes[i]));
            let eigen_errors = selected.map(|i| hz_error(&report.modes[i]));
            let h = (L / (24 * r) as f64)
                .max(W / (2 * r) as f64)
                .max(T / r as f64)
                * design.scale;
            row(
                "structural-rung",
                "h",
                "pass",
                &format!(
                    "{{\"bar_id\":\"{bar}\",\"element\":\"TET10\",\"rung_h_m\":{h},\"n_tets\":{n_tets},\"free_dofs\":{},\"f_modes_hz\":{frequencies:?},\"eigen_uncertainty_hz\":{eigen_errors:?},\"factor_peak_bytes\":{},\"factor_nnz_l\":{},\"window_hz\":[20,{}],\"inertia_count\":{}}}",
                    a.free_dofs.len(),
                    report.stats.factor_peak_bytes,
                    report.stats.factor_nnz_l,
                    FMAX / design.scale,
                    report.expected
                ),
                started,
            );
            rows.push(Rung {
                h,
                n_tets,
                frequencies,
                eigen_errors,
            });
            let rigid = if r == *TET10_RUNGS.last().unwrap() {
                rigid_projectors(&mesh, &nodes, &a)
            } else {
                Vec::new()
            };
            last = Some((mesh, report.modes, selected, rigid));
        }
        let (mesh, modes, selected, rigid_projectors) = last.unwrap();
        Ladder {
            design,
            rungs: rows,
            mesh,
            modes,
            selected,
            rigid_projectors,
        }
    })
}

fn uniform_tet10() -> &'static Ladder {
    UNIFORM_TET10.get_or_init(|| calculate_tet10_ladder(Design::PRISM, "A-uniform"))
}

/// Exact free-free Timoshenko characteristic equations, split into parity
/// about the midpoint. Dimensionless lambda = rho*A*omega²*L⁴/(E*I),
/// r=I/(A*L²), s=E*I/(kappa*G*A*L²). The coupled ODE gives
/// y²+lambda*(r+s)*y+r*s*lambda²-lambda=0, y=alpha²,-beta².
/// M=0 and V=0 at both ends give the determinants below. kappa=5/6
/// is the declared rectangular-section energy correction, not a fit.
fn beam_frequencies(shear: f64, rotary: f64) -> [f64; 3] {
    let characteristic = |q: f64, even: bool| {
        let lambda = q.powi(4);
        let sum = lambda * (rotary + shear);
        let disc = (sum * sum + 4.0 * (lambda - rotary * shear * lambda * lambda)).sqrt();
        let alpha = (2.0 * (lambda - rotary * shear * lambda * lambda) / (disc + sum)).sqrt();
        let beta = (0.5 * (disc + sum)).sqrt();
        let a = alpha * (beta * beta - rotary * lambda);
        let b = beta * (alpha * alpha + rotary * lambda);
        if even {
            a * det::sin(beta / 2.0) + b * det::cos(beta / 2.0) * det::tanh(alpha / 2.0)
        } else {
            a * det::cos(beta / 2.0) * det::tanh(alpha / 2.0) - b * det::sin(beta / 2.0)
        }
    };
    let mut roots = Vec::new();
    for even in [true, false] {
        let mut lo = 0.2;
        for j in 1..=1800 {
            let hi = 0.2 + j as f64 * 0.01;
            if characteristic(lo, even) * characteristic(hi, even) < 0.0 {
                let (mut a, mut b) = (lo, hi);
                for _ in 0..64 {
                    let mid = (a + b) / 2.0;
                    if characteristic(a, even) * characteristic(mid, even) <= 0.0 {
                        b = mid;
                    } else {
                        a = mid;
                    }
                }
                roots.push((a + b) / 2.0);
            }
            lo = hi;
        }
    }
    roots.sort_by(f64::total_cmp);
    assert!(roots.len() >= 3);
    core::array::from_fn(|i| {
        roots[i].powi(2) * (young() * T * T / (12.0 * RHO)).sqrt() / (TAU * L * L)
    })
}

fn beam_for(pack: &str) -> [f64; 3] {
    let mat = law(pack);
    let rotary = T * T / (12.0 * L * L);
    beam_frequencies(young() * rotary / ((5.0 / 6.0) * mat.g[2]), rotary)
}

struct Estimate {
    orders: [f64; 2],
    limit_hz: f64,
    uncertainty_hz: f64,
    asymptotic: bool,
}

fn estimate(ladder: &Ladder, mode: usize) -> Estimate {
    let f: Vec<_> = ladder.rungs.iter().map(|r| r.frequencies[mode]).collect();
    let differences = [f[0] - f[1], f[1] - f[2], f[2] - f[3]];
    let orders = core::array::from_fn(|i| {
        det::ln((differences[i] / differences[i + 1]).abs()) / det::ln(2.0)
    });
    // Only the final THREE rungs support the asymptotic claim. The coarse
    // diagnostic order is retained so P1 locking stays visible in the log.
    let asymptotic = differences.iter().all(|d| *d > 0.0)
        && orders[1].is_finite()
        && (1.7..=2.3).contains(&orders[1])
        && ladder.rungs[3].eigen_errors[mode] + ladder.rungs[2].eigen_errors[mode]
            < 0.1 * differences[2];
    let correction = differences[2] / 3.0;
    let eig = ladder.rungs[3].eigen_errors[mode] + ladder.rungs[2].eigen_errors[mode] / 3.0;
    Estimate {
        orders,
        limit_hz: f[3] - correction,
        // Conservative diagnostic Richardson allowance, not a certified bound.
        uncertainty_hz: 2.0 * correction.abs() + eig,
        asymptotic,
    }
}

fn h_pow(h: f64, p: f64) -> f64 {
    det::exp(p * det::ln(h))
}

/// Observed order from three rungs with arbitrary refinement ratios: the root
/// p of (f0-f1)/(f1-f2) = (h0^p-h1^p)/(h1^p-h2^p). The right side increases
/// with p for h0 > h1 > h2; outside the [0.25, 12] bracket there is no order.
fn observed_order(h: [f64; 3], f: [f64; 3]) -> f64 {
    let target = (f[0] - f[1]) / (f[1] - f[2]);
    let ratio = |p: f64| (h_pow(h[0], p) - h_pow(h[1], p)) / (h_pow(h[1], p) - h_pow(h[2], p));
    let (mut lo, mut hi) = (0.25, 12.0);
    if !(target.is_finite() && ratio(lo) <= target && target <= ratio(hi)) {
        return f64::NAN;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if ratio(mid) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// TET10 counterpart of [`estimate`]: formal eigenfrequency order 2p = 4,
/// the canonical owner's admitted band [3.5, 4.3] (fs-solid tet10 benchmarks),
/// and Richardson at the formal order with the actual final h ratio.
fn tet10_estimate(rungs: &[Rung], mode: usize) -> Estimate {
    let f: Vec<_> = rungs.iter().map(|r| r.frequencies[mode]).collect();
    let h: Vec<_> = rungs.iter().map(|r| r.h).collect();
    let orders = [
        observed_order([h[0], h[1], h[2]], [f[0], f[1], f[2]]),
        observed_order([h[1], h[2], h[3]], [f[1], f[2], f[3]]),
    ];
    let last = f[2] - f[3];
    let weight = h_pow(h[3], 4.0) / (h_pow(h[2], 4.0) - h_pow(h[3], 4.0));
    let correction = last * weight;
    let [e2, e3] = [rungs[2].eigen_errors[mode], rungs[3].eigen_errors[mode]];
    let asymptotic = f.windows(2).all(|w| w[0] > w[1])
        && (3.5..=4.3).contains(&orders[1])
        && e2 + e3 < 0.1 * last;
    Estimate {
        orders,
        limit_hz: f[3] - correction,
        // Same diagnostic Richardson allowance as P1, not a certified bound.
        uncertainty_hz: 2.0 * correction.abs() + e3 + (e2 + e3) * weight,
        asymptotic,
    }
}

#[test]
fn uniform_prism_3d_ladder_converges_to_free_free_beam_theory() {
    let _guard = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let p1 = uniform();
    let tet10 = uniform_tet10();
    let timoshenko = beam_for(HONDURAS);
    let eb = beam_frequencies(0.0, 0.0);
    for mode in 0..3 {
        let started = Instant::now();
        // The P1 ladder is reported against its formal order 2, never gated:
        // at 147456 tetrahedra it is still locking (DEF-08).
        let locked = estimate(p1, mode);
        row(
            "uniform-prism-P1-locking-diagnostic",
            KINDS[mode],
            if locked.asymptotic { "pass" } else { "refused" },
            &format!(
                "{{\"element\":\"P1\",\"observed_order\":{:?},\"formal_order\":2,\"extrapolated_hz\":{},\"oracle_hz\":{},\"route\":\"TET10:M1.3-W1a\"}}",
                locked.orders, locked.limit_hz, timoshenko[mode]
            ),
            started,
        );
        let est = tet10_estimate(&tet10.rungs, mode);
        let error = (est.limit_hz / timoshenko[mode] - 1.0).abs();
        row(
            "uniform-prism-vs-free-free-Timoshenko",
            KINDS[mode],
            if est.asymptotic && error <= BEAM_TOL {
                "pass"
            } else {
                "fail"
            },
            &format!(
                "{{\"element\":\"TET10\",\"rungs\":{TET10_RUNGS:?},\"f_hz\":{:?},\"observed_order\":{:?},\"formal_order\":4,\"order_band\":[3.5,4.3],\"extrapolated_hz\":{},\"oracle_hz\":{},\"euler_bernoulli_hz\":{},\"relative_error\":{error},\"tolerance\":{BEAM_TOL},\"numerical_uncertainty_hz\":{}}}",
                tet10
                    .rungs
                    .iter()
                    .map(|r| r.frequencies[mode])
                    .collect::<Vec<_>>(),
                est.orders,
                est.limit_hz,
                timoshenko[mode],
                eb[mode],
                est.uncertainty_hz
            ),
            started,
        );
        assert!(
            timoshenko[mode] < eb[mode],
            "shear/rotary inertia LOWER the bending branch; bead's 'above' is a sign typo"
        );
        assert!(
            est.asymptotic,
            "TET10 pre-asymptotic orders {:?}; no convergence claim",
            est.orders
        );
        assert!(
            error <= BEAM_TOL,
            "vertical {}: 3-D {} vs analytic Timoshenko {} (frozen 1%)",
            mode + 1,
            est.limit_hz,
            timoshenko[mode]
        );
    }
}

#[test]
fn timoshenko_oracle_has_euler_bernoulli_limit_and_live_shear_rotary_axes() {
    let started = Instant::now();
    let r = T * T / (12.0 * L * L);
    let s = young() * r / ((5.0 / 6.0) * law(HONDURAS).g[2]);
    let eb = beam_frequencies(0.0, 0.0);
    let shear = beam_frequencies(s, 0.0);
    let rotary = beam_frequencies(0.0, r);
    let both = beam_frequencies(s, r);
    let scale = (young() * T * T / (12.0 * RHO)).sqrt() / (TAU * L * L);
    for mode in 0..3 {
        let q = (eb[mode] / scale).sqrt();
        let cosh = (det::exp(q) + det::exp(-q)) / 2.0;
        assert!((det::cos(q) - 1.0 / cosh).abs() < 1e-12);
        assert!(shear[mode] < eb[mode] && rotary[mode] < eb[mode]);
        assert!(both[mode] < shear[mode].min(rotary[mode]));
    }
    row(
        "analytic-beam-oracle-controls",
        "shear-and-rotary",
        "pass",
        &format!(
            "{{\"euler_bernoulli_hz\":{eb:?},\"shear_only_hz\":{shear:?},\"rotary_only_hz\":{rotary:?},\"Timoshenko_hz\":{both:?},\"characteristic_residual_tolerance\":1e-12}}"
        ),
        started,
    );
}

#[test]
fn free_free_bar_has_six_rigid_modes_below_1_hz() {
    let _guard = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let started = Instant::now();
    with_cx(|cx| {
        let mesh = mesh([6, 1, 1], Design::PRISM);
        let a = assemble(&mesh, &material(HONDURAS), cx);
        let modes = eigh_gen_dense(
            &a.stiffness.to_dense(),
            &a.mass.to_dense(),
            a.free_dofs.len(),
        )
        .unwrap();
        assert_eq!(
            modes
                .iter()
                .filter(|m| m.lambda.abs().sqrt() / TAU < 1.0)
                .count(),
            6
        );
        assert!(hz(&modes[6]) > 100.0);
        row(
            "six-rigid-modes",
            "boundary-condition",
            "pass",
            "{\"observed\":6,\"oracle\":6,\"tolerance_hz\":1}",
            started,
        );
    });
}

#[test]
fn slice_window_matches_dense_oracle_on_coarse_rung() {
    let _guard = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let started = Instant::now();
    with_cx(|cx| {
        let mesh = mesh([6, 1, 1], Design::PRISM);
        let a = assemble(&mesh, &material(HONDURAS), cx);
        let dense = eigh_gen_dense(
            &a.stiffness.to_dense(),
            &a.mass.to_dense(),
            a.free_dofs.len(),
        )
        .unwrap();
        // Compare the first ten elastic modes. The six near-zero eigenvalues
        // are checked separately: relative error at zero is undefined.
        let high = 0.5 * (hz(&dense[15]) + hz(&dense[16]));
        let actual = slice(&a, 1.0, high);
        assert_eq!(actual.modes.len(), 10);
        for (i, (got, oracle)) in actual.modes.iter().zip(&dense[6..16]).enumerate() {
            let error = (got.lambda / oracle.lambda - 1.0).abs();
            assert!(error <= 1e-9, "elastic mode {i}: {error:e}");
        }
        row(
            "sparse-vs-dense-first-ten-elastic",
            "eigensolver",
            "pass",
            "{\"mode_count\":10,\"tolerance_relative_eigenvalue\":1e-9}",
            started,
        );
    });
}

#[test]
fn mode_window_cap_is_live() {
    let _guard = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let started = Instant::now();
    with_cx(|cx| {
        let mesh = mesh([6, 1, 1], Design::PRISM);
        let a = assemble(&mesh, &material(HONDURAS), cx);
        let dense = eigh_gen_dense(
            &a.stiffness.to_dense(),
            &a.mass.to_dense(),
            a.free_dofs.len(),
        )
        .unwrap();
        let lower = 0.5 * (hz(&dense[6]) + hz(&dense[7]));
        let upper = 0.5 * (hz(&dense[7]) + hz(&dense[8]));
        assert_eq!(slice(&a, 1.0, lower).expected, 1);
        assert_eq!(slice(&a, 1.0, upper).expected, 2);
        // SliceOptions has a Lanczos WORK cap, not a number-of-modes cap.
        // Both an exactly exhausted one-step budget and the additional mode
        // refuse, rather than publishing a shortened certified window.
        let limited = SliceOptions {
            max_lanczos: 1,
            ..options()
        };
        for high in [lower, upper] {
            assert!(matches!(
                slice_window(
                    &a.stiffness,
                    &a.mass,
                    ((TAU * 1.0).powi(2), (TAU * high).powi(2)),
                    &limited
                ),
                Err(fs_modal::ModalError::WindowUnresolved { .. })
            ));
        }
        row(
            "window-and-work-cap-live",
            "mode-window",
            "pass",
            "{\"window_mode_counts\":[1,2],\"max_lanczos\":1,\"refused_counts\":[1,2]}",
            started,
        );
    });
}

/// The design is solved on the locking-free TET10 owner at rung 3. On the
/// coarse P1 grid that solved it before, shear locking lifted the
/// fundamental, so the converged bar sounded about 330 cents under Db5.
const DESIGN_GRID: [usize; 3] = [72, 6, 3];

/// A design solved by an earlier run of the bisection below, as
/// `depth=<f64>` and `scale=<f64>` lines in the file named by
/// BAR_S1_DESIGN_CACHE. Never trusted unchecked: see `tune_design`.
fn cached_design() -> Option<(f64, f64)> {
    let path = std::env::var_os("BAR_S1_DESIGN_CACHE")?;
    let text = std::fs::read_to_string(path).expect("readable design cache");
    let value = |key: &str| -> f64 {
        text.lines()
            .find_map(|line| line.strip_prefix(key))
            .unwrap_or_else(|| panic!("design cache lacks {key}"))
            .trim()
            .parse()
            .expect("design cache number")
    };
    Some((value("depth="), value("scale=")))
}

fn tune_design() -> Design {
    with_cx(|cx| {
        let mat = material(HONDURAS);
        let evaluate = |depth: f64| {
            let design = Design { depth, scale: 1.0 };
            let mesh = mesh(DESIGN_GRID, design);
            let (_, _, a) = tet10_assemble(&mesh, design, &mat, cx);
            let report = slice(&a, 20.0, FMAX);
            let ids = classify(&mesh, &report.modes);
            (
                hz(&report.modes[ids[1]]) / hz(&report.modes[ids[0]]),
                hz(&report.modes[ids[0]]),
            )
        };
        if let Some((depth, scale)) = cached_design() {
            // Accepted only after ONE fresh TET10 solve reproduces the 1:3
            // ratio and the Db5 fundamental that fixed `scale`.
            let started = Instant::now();
            let (ratio, fundamental) = evaluate(depth);
            assert!(
                (ratio / 3.0 - 1.0).abs() < 1e-4,
                "cached design no longer gives 1:3: {ratio}"
            );
            assert!(
                (fundamental / ET_DB5_HZ / scale - 1.0).abs() < 1e-9,
                "cached scale no longer reproduces the Db5 fundamental: {fundamental}"
            );
            row(
                "synthetic-undercut-construction",
                "cut-depth",
                "pass",
                &format!(
                    "{{\"bar_id\":\"B-undercut\",\"depth_fraction\":{depth},\"uniform_scale\":{scale},\"design_cache\":\"reused; re-verified by one TET10 solve\",\"achieved_design_grid_ratio\":{ratio},\"design_grid\":{DESIGN_GRID:?},\"same_specimen\":false}}"
                ),
                started,
            );
            return Design { depth, scale };
        }
        let (mut lo, mut hi) = (0.0, 0.85);
        let left = evaluate(lo).0 - 3.0;
        let right = evaluate(hi).0 - 3.0;
        assert!(
            left * right < 0.0,
            "declared cut family must bracket 1:3: {left}, {right}"
        );
        for _ in 0..18 {
            cx.checkpoint().expect("cancel between construction solves");
            let mid = (lo + hi) / 2.0;
            if (evaluate(mid).0 - 3.0).signum() == left.signum() {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let depth = (lo + hi) / 2.0;
        let (ratio, fundamental) = evaluate(depth);
        assert!((ratio / 3.0 - 1.0).abs() < 1e-4, "model-based design solve");
        // Uniform geometric scaling changes every eigenfrequency by 1/scale,
        // preserving ratios and orthotropic material axes exactly.
        let design = Design {
            depth,
            scale: fundamental / ET_DB5_HZ,
        };
        row(
            "synthetic-undercut-construction",
            "cut-depth",
            "pass",
            &format!(
                "{{\"bar_id\":\"B-undercut\",\"profile\":\"(1-((x/L-0.5)/0.3)^2)^2 on central 60 percent\",\"depth_fraction\":{depth},\"uniform_scale\":{},\"targeted_ratio\":3,\"achieved_design_grid_ratio\":{ratio},\"targeted_et_note\":\"Db5\",\"target_hz\":{ET_DB5_HZ},\"design_grid\":{DESIGN_GRID:?},\"design_element\":\"TET10\",\"target_is_not_validation\":true,\"same_specimen\":false}}",
                design.scale
            ),
            Instant::now(),
        );
        design
    })
}

fn undercut() -> &'static Ladder {
    UNDERCUT.get_or_init(|| calculate_tet10_ladder(tune_design(), "B-undercut"))
}

#[test]
fn undercut_bar_ladder_reports_modes_with_uncertainty() {
    let _guard = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let ladder = undercut();
    let mut qualified = true;
    for (mode, kind) in KINDS.iter().enumerate() {
        let started = Instant::now();
        let est = tet10_estimate(&ladder.rungs, mode);
        let f = ladder.rungs.last().unwrap().frequencies[mode];
        let cents = 1200.0 * det::ln(1.0 + est.uncertainty_hz / f) / det::ln(2.0);
        let admitted = est.asymptotic && cents <= 1.0;
        if mode < 2 {
            qualified &= admitted;
        }
        row(
            "undercut-modal-uncertainty",
            kind,
            if admitted { "pass" } else { "refused" },
            &format!(
                "{{\"bar_id\":\"B-undercut\",\"element\":\"TET10\",\"f_hz\":{f},\"extrapolated_hz\":{},\"observed_order\":{:?},\"formal_order\":4,\"order_band\":[3.5,4.3],\"uncertainty_hz\":{},\"uncertainty_cents\":{cents},\"tolerance_cents\":1,\"uncertainty_kind\":\"Richardson estimate plus eigen residual; not certified\",\"asymptotic\":{},\"route\":\"{}\"}}",
                est.limit_hz,
                est.orders,
                est.uncertainty_hz,
                est.asymptotic,
                if admitted {
                    "TET10-qualified-numerical-band-only"
                } else {
                    "finer TET10 rung"
                }
            ),
            started,
        );
        assert!(est.uncertainty_hz.is_finite() && est.uncertainty_hz > 0.0);
        // A rejected one-cent band is explicitly allowed by DONE(b). It must
        // never be represented as a qualified render by the script/metadata.
        if !admitted {
            assert!(cents > 1.0 || !est.asymptotic);
        }
    }
    // The design target is reported against the extrapolated fundamental; it
    // is a construction target, never a validation.
    let fundamental = tet10_estimate(&ladder.rungs, 0).limit_hz;
    row(
        "undercut-band",
        "h",
        if qualified { "pass" } else { "refused" },
        &format!(
            "{{\"targeted_partials\":[1,2],\"ratio_target\":3,\"finest_ratio\":{},\"extrapolated_fundamental_hz\":{fundamental},\"target_hz\":{ET_DB5_HZ},\"fundamental_cents_vs_target\":{},\"route\":\"finer TET10 rung when refused\"}}",
            ladder.rungs[3].frequencies[1] / ladder.rungs[3].frequencies[0],
            1200.0 * det::ln(fundamental / ET_DB5_HZ) / det::ln(2.0)
        ),
        Instant::now(),
    );
}

#[test]
fn undercut_ratio_disclosed_in_fixture_metadata() {
    let _guard = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let started = Instant::now();
    assert_eq!(young(), 13.53e9);
    assert!(
        !ROSEWOOD.contains("\ter_over_el\t"),
        "if rosewood data land, replace this scenario policy explicitly"
    );
    // Bounds are the extrema of this predeclared correlated scenario SET.
    // They are not independent Cartesian parameter draws or species coverage.
    let properties = [
        "er_over_el",
        "et_over_el",
        "glr_over_el",
        "grt_over_el",
        "glt_over_el",
        "nu_lr",
        "nu_lt",
        "nu_rt",
    ];
    for property in properties {
        let values = [HONDURAS, AFRICAN, WALNUT].map(|pack| scalar(pack, property));
        let min = values.into_iter().fold(f64::INFINITY, f64::min);
        let max = values.into_iter().fold(f64::NEG_INFINITY, f64::max);
        assert!(min.is_finite() && min > 0.0 && max >= min);
        row(
            "transverse-material-scenario-envelope",
            property,
            "pass",
            &format!(
                "{{\"scenario_values\":{values:?},\"bounds\":[{min},{max}],\"source\":\"FPL-GTR-282 tables5-1,5-2: Honduras mahogany, African mahogany, black walnut\",\"rosewood_coverage_claim\":false,\"unit\":\"1\"}}"
            ),
            started,
        );
    }
    let beam = [HONDURAS, AFRICAN, WALNUT].map(beam_for);
    row(
        "material-uncertainty-not-random",
        "transverse-constants",
        "pass",
        &format!(
            "{{\"beam_hz_by_correlated_scenario\":{beam:?},\"nominal_scenario\":\"Honduras mahogany ratios applied to Indian rosewood E_L\",\"density_kg_m3\":{RHO},\"density_authority\":\"authored Estimate\",\"wood_loss\":\"unknown:M8-W2c\",\"receiver_m\":{RECEIVER:?},\"receiver_authority\":\"single declared observation prior; not measured Iowa microphone\",\"targeted_ratio\":3,\"same_specimen\":false}}"
        ),
        started,
    );
    for pack in [HONDURAS, AFRICAN, WALNUT] {
        material(pack);
    }
}

/// The six parametric boundary faces of the same Kuhn mesh. Each quad uses
/// its low/low--high/high diagonal. Orientation is checked using the fixed
/// face-normal component (also valid on the curved lower undercut face).
fn boundary(mesh: &Mesh) -> Vec<[usize; 3]> {
    let counts = mesh.counts;
    let index = |q: [usize; 3]| (q[0] * (counts[1] + 1) + q[1]) * (counts[2] + 1) + q[2];
    let mut faces = Vec::new();
    for axis in 0..3 {
        let other: Vec<_> = (0..3).filter(|&a| a != axis).collect();
        for high in [false, true] {
            for i in 0..counts[other[0]] {
                for j in 0..counts[other[1]] {
                    let mut q = [0; 3];
                    q[axis] = if high { counts[axis] } else { 0 };
                    q[other[0]] = i;
                    q[other[1]] = j;
                    let a = index(q);
                    q[other[0]] += 1;
                    let b = index(q);
                    q[other[1]] += 1;
                    let c = index(q);
                    q[other[0]] -= 1;
                    let d = index(q);
                    for mut face in [[a, b, c], [a, c, d]] {
                        let p = face.map(|v| mesh.nodes[v]);
                        let ab: [f64; 3] = core::array::from_fn(|k| p[1][k] - p[0][k]);
                        let ac: [f64; 3] = core::array::from_fn(|k| p[2][k] - p[0][k]);
                        let normal_component = ab[(axis + 1) % 3] * ac[(axis + 2) % 3]
                            - ab[(axis + 2) % 3] * ac[(axis + 1) % 3];
                        if (normal_component > 0.0) != high {
                            face.swap(1, 2);
                        }
                        faces.push(face);
                    }
                }
            }
        }
    }
    faces
}

/// Restrict the finest P1 boundary displacement to nested boundary nodes.
/// This keeps the structural mode fixed while the BEM axis changes. Every
/// selected coarse node has an exact corresponding finest node; no nearest
/// point search or shape-fit can fabricate pressure convergence.
fn surface_and_shapes(ladder: &Ladder, refinement: usize) -> (SpherePanels, Vec<Vec<f64>>) {
    let counts = [12 * refinement, 2 * refinement, refinement];
    let coarse = mesh(counts, ladder.design);
    let faces = boundary(&coarse);
    let surface =
        SpherePanels::from_triangles(faces.iter().map(|f| f.map(|v| coarse.nodes[v])).collect())
            .unwrap();
    let fine = ladder.mesh.counts;
    let map: Vec<_> = (0..coarse.nodes.len())
        .map(|node| {
            let k = node % (counts[2] + 1);
            let j = node / (counts[2] + 1) % (counts[1] + 1);
            let i = node / ((counts[1] + 1) * (counts[2] + 1));
            assert!(fine.iter().zip(counts).all(|(a, b)| a % b == 0));
            ((i * fine[0] / counts[0]) * (fine[1] + 1) + j * fine[1] / counts[1]) * (fine[2] + 1)
                + k * fine[2] / counts[2]
        })
        .collect();
    let shapes = ladder
        .selected
        .iter()
        .map(|&m| {
            faces
                .iter()
                .zip(surface.normals())
                .map(|(face, normal)| {
                    face.iter()
                        .map(|&v| {
                            (0..3)
                                .map(|c| ladder.modes[m].phi[3 * map[v] + c] * normal[c])
                                .sum::<f64>()
                        })
                        .sum::<f64>()
                        / 3.0
                })
                .collect()
        })
        .collect();
    (surface, shapes)
}

struct RadiationLadder {
    pressures: Vec<[C64; 5]>,
    panel_counts: Vec<usize>,
}

fn radiation() -> &'static RadiationLadder {
    RADIATION.get_or_init(|| {
        let ladder = undercut();
        let mut pressures = Vec::new();
        let mut panel_counts = Vec::new();
        for refinement in [1, 2, 4] {
            let started = Instant::now();
            let (surface, shapes) = surface_and_shapes(ladder, refinement);
            let mut p = [C64::ZERO; 5];
            for mode in 0..5 {
                let omega = TAU * ladder.rungs[3].frequencies[mode];
                // Unit generalized modal velocity, NOT displacement. H has
                // units Pa s/(m sqrt(kg)); physical exp(-iwt) velocity of a
                // unit displacement is -i*w*phi, and the corresponding
                // displacement transfer is therefore -i*w*H.
                let velocity: Vec<_> = shapes[mode].iter().map(|&phi| C64::from_re(phi)).collect();
                let solved = solve_radiation_batch(&surface, omega / Medium::air().sound_speed,
                    Medium::air(), &[velocity.as_slice()], Formulation::BurtonMiller).unwrap();
                p[mode] = exterior_pressure_at_points(&surface, &solved[0], Medium::air(), &[RECEIVER]).unwrap()[0];
                assert!(p[mode].re.is_finite() && p[mode].im.is_finite());
                let displacement = p[mode] * C64::new(0.0, -omega);
                row("Burton-Miller-modal-receiver", KINDS[mode], "pass", &format!(
                    "{{\"bar_id\":\"B-undercut\",\"panel_count\":{},\"frequency_hz\":{},\"p_receiver_pa_per_unit_modal_velocity\":[{},{}],\"p_receiver_pa_per_unit_modal_displacement\":[{},{}],\"receiver_m\":{RECEIVER:?},\"formulation\":\"BurtonMiller\",\"panels_per_wavelength\":{}}}",
                    surface.centroids().len(), omega / TAU, p[mode].re, p[mode].im,
                    displacement.re, displacement.im, solved[0].panels_per_wavelength), started);
            }
            panel_counts.push(surface.centroids().len());
            pressures.push(p);
        }
        RadiationLadder { pressures, panel_counts }
    })
}

#[test]
fn bem_panel_ladder_converges_modal_radiated_pressure() {
    let _guard = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let radiation = radiation();
    let finest = radiation.pressures.last().unwrap();
    let dominant = finest.iter().map(|p| p.abs()).fold(0.0, f64::max);
    for (mode, kind) in KINDS.iter().enumerate() {
        let started = Instant::now();
        let p: Vec<_> = radiation.pressures.iter().map(|r| r[mode]).collect();
        let differences: Vec<f64> = p.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
        let [d0, d1] = [differences[0], differences[1]];
        let order = det::ln(d0 / d1) / det::ln(2.0);
        assert!(
            d0.is_finite() && d1.is_finite() && d0 > d1 && d1 > 0.0,
            "{kind} complex receiver pressure must converge: {p:?}; differences {differences:?}"
        );
        let uncertainty = d1 / (d0 / d1 - 1.0);
        let bound = finest[mode].abs() + uncertainty;
        // A mode that is not resolved to 10% of itself may still be claimed
        // as a near-null radiator, only if its whole bound |p| + u lies under
        // the near-null threshold; the same 10% band then applies to that
        // threshold scale. A resolved mode is never downgraded to a bound.
        let resolved = uncertainty <= 0.1 * finest[mode].abs();
        let near_null = !resolved && bound <= NEAR_NULL_FRACTION * dominant;
        let (scale, claim) = if near_null {
            (
                NEAR_NULL_FRACTION * dominant,
                "near-null radiator: bounded, not resolved",
            )
        } else {
            (finest[mode].abs(), "resolved")
        };
        let relative = uncertainty / scale;
        row(
            "panel-pressure-convergence",
            kind,
            if relative <= 0.1 { "pass" } else { "fail" },
            &format!(
                "{{\"panel_counts\":{:?},\"differences_pa_per_unit_modal_velocity\":{differences:?},\"observed_order\":{order},\"uncertainty_pa_per_unit_modal_velocity\":{uncertainty},\"claim\":\"{claim}\",\"gate_scale_pa_per_unit_modal_velocity\":{scale},\"dominant_pa_per_unit_modal_velocity\":{dominant},\"near_null_fraction\":{NEAR_NULL_FRACTION},\"bound_pa_per_unit_modal_velocity\":{bound},\"relative_uncertainty\":{relative},\"tolerance_relative\":0.1,\"uncertainty_kind\":\"three-rung Richardson estimate; phase retained\"}}",
                radiation.panel_counts
            ),
            started,
        );
        assert!(
            relative <= 0.1,
            "{kind} ({claim}) receiver uncertainty {relative} of {scale}; refine panels, do not loosen"
        );
    }
}

/// Ideal elastic silicone head. Johnston et al., JMM24(2014)035017,
/// doi:10.1088/0960-1317/24/3/035017, table1: 25C-cured Sylgard184
/// E=1.32 +/-0.07 MPa (reported spread, not a confidence interval).
/// https://researchprofiles.herts.ac.uk/files/13511316/906746.pdf
/// nu=0.49 in the DECLARED [0.45,0.499] near-incompressible design
/// interval; the paper does not measure nu. Radius, mass, speed are authored
/// fixture inputs, not an Iowa mallet identification. Rate loss is unknown.
const HEAD_E: f64 = 1.32e6;
const HEAD_NU: f64 = 0.49;
const HEAD_R: f64 = 0.008;
const HEAD_MASS: f64 = 0.002;
const STROKE: f64 = 0.10;

fn hertz_interface() -> fs_tribo::InterfaceSystemRef {
    fs_tribo::InterfaceSystemRef::new(
        "bar-s1:ideal-silicone-head/target",
        "first-strike",
        "Johnston2014:table1-E; declared nu interval; elastic-only diagnostic",
        fs_tribo::InputAuthority::Estimated,
        fs_tribo::InterfaceMedium::Dry,
    )
    .unwrap()
}

fn hertz_law(e: f64, nu: f64) -> fs_tribo::HertzSpherePlane {
    fs_tribo::HertzSpherePlane {
        effective_radius: HEAD_R,
        reduced_modulus: e / (1.0 - nu * nu),
    }
}

fn head_modulus_relative_uncertainty() -> f64 {
    let nominal = hertz_law(HEAD_E, HEAD_NU).reduced_modulus;
    [1.25e6, 1.39e6]
        .into_iter()
        .flat_map(|e| {
            [0.45, 0.499].map(move |nu| (hertz_law(e, nu).reduced_modulus / nominal - 1.0).abs())
        })
        .fold(0.0, f64::max)
}

fn force(
    law: fs_tribo::HertzSpherePlane,
    interface: &fs_tribo::InterfaceSystemRef,
    indentation: f64,
) -> f64 {
    law.response(interface, indentation.max(0.0))
        .expect("material-derived Hertz force")
        .normal_force_n
}

fn closed_contact_time(law: fs_tribo::HertzSpherePlane) -> (f64, f64) {
    let delta = det::pow(
        15.0 * HEAD_MASS * STROKE * STROKE / (16.0 * law.reduced_modulus * HEAD_R.sqrt()),
        0.4,
    );
    (2.94 * delta / STROKE, delta)
}

/// Independent time integration of the existing force owner. The oracle
/// only selects the work grid; it never supplies a force or a trajectory.
fn rigid_impact(law: fs_tribo::HertzSpherePlane, steps: usize) -> (f64, f64) {
    let interface = hertz_interface();
    let (oracle, _) = closed_contact_time(law);
    let h = oracle / steps as f64;
    let (mut x, mut v) = (0.0, STROKE);
    let mut f = force(law, &interface, x);
    for step in 0..2 * steps {
        let old = x;
        v -= 0.5 * h * f / HEAD_MASS;
        x += h * v;
        f = force(law, &interface, x);
        v -= 0.5 * h * f / HEAD_MASS;
        if x <= 0.0 && step > 0 {
            let tc = h * (step as f64 + old / (old - x));
            return (tc, (v * v / (STROKE * STROKE) - 1.0).abs());
        }
    }
    panic!("Hertz impact failed to separate within its declared work cap");
}

#[test]
fn material_hertz_mallet_contact_time_matches_closed_form() {
    let _guard = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    use fs_contact::normal_patch::*;
    let started = Instant::now();
    let law = hertz_law(HEAD_E, HEAD_NU);
    let (oracle, delta) = closed_contact_time(law);
    let request = NormalPatchRequest {
        identity: NormalPatchIdentity {
            model_id: "bar-s1-Hertz".into(),
            source_id: "Johnston2014:table1; authored ideal head".into(),
            state_id: "elastic-only".into(),
        },
        interface: hertz_interface(),
        law: NormalPatchLaw::HertzSpherePlane {
            effective_radius_m: HEAD_R,
            reduced_modulus_pa: law.reduced_modulus,
        },
        geometry: NormalPatchGeometry::SpherePlane,
        indentation_m: delta,
        indentation_rate_m_per_s: 0.0,
        step_s: oracle / 2048.0,
        line_load_n_per_m: 0.0,
        applicability: ApplicabilityInput {
            half_space_depth_m: 1.0,
            layer_thickness_m: 1.0,
            yield_strength_pa: 3.0e6,
            characteristic_rate_m_per_s: STROKE,
            temperature_k: 293.15,
            adhesion_energy_j_per_m2: 0.0,
        },
        limits: ApplicabilityLimits {
            max_patch_to_radius: 0.2,
            max_strain: 0.2,
            max_patch_to_depth: 0.1,
            max_patch_to_layer: 0.1,
            max_pressure_to_yield: 0.1,
            max_rate_ratio: 1.0,
            min_temperature_k: 290.0,
            max_temperature_k: 300.0,
        },
        uncertainty: InputUncertainty {
            radius_relative: 0.0,
            modulus_relative: head_modulus_relative_uncertainty(),
            load_relative: 0.0,
        },
    };
    // The applicability numbers describe this ideal rigid-target exercise,
    // not a claim of a measured yield surface or an actual bar contact patch.
    let NormalPatchReceipt::Point(receipt) = request.evaluate().unwrap() else {
        panic!("point law")
    };
    let reference_force = force(law, &hertz_interface(), delta);
    assert!((receipt.normal_force_n / reference_force - 1.0).abs() < 1e-12);
    for steps in [512, 1024, 2048] {
        let (actual, energy_error) = rigid_impact(law, steps);
        let relative = (actual / oracle - 1.0).abs();
        row(
            "Hertz-rigid-target-contact-time",
            "time-step",
            if relative <= 0.01 { "pass" } else { "fail" },
            &format!(
                "{{\"steps_per_contact\":{steps},\"contact_time_s\":{actual},\"oracle_s\":{oracle},\"delta_max_m\":{delta},\"relative_error\":{relative},\"tolerance\":0.01,\"rebound_energy_relative_error\":{energy_error},\"E_pa\":{HEAD_E},\"E_reported_spread_pa\":70000,\"nu\":{HEAD_NU},\"nu_declared_bounds\":[0.45,0.499],\"radius_m\":{HEAD_R},\"mass_kg\":{HEAD_MASS},\"stroke_m_s\":{STROKE},\"stroke_authority\":\"declared\"}}"
            ),
            started,
        );
        assert!(
            relative <= 0.01,
            "frozen 1% Hertz contact-time gate: {relative}"
        );
        assert!(energy_error < 2e-5, "elastic rebound energy {energy_error}");
    }
    // Retain epistemic input envelope, separate from numerical error.
    let times: Vec<_> = [1.25e6, 1.39e6]
        .into_iter()
        .flat_map(|e| [0.45, 0.499].map(move |nu| closed_contact_time(hertz_law(e, nu)).0))
        .collect();
    row(
        "Hertz-input-envelope",
        "material",
        "pass",
        &format!(
            "{{\"contact_time_s_scenarios\":{times:?},\"source\":\"Johnston et al 2014 doi:10.1088/0960-1317/24/3/035017 table1\",\"elastic_idealization\":true,\"viscoelastic_loss\":\"unknown\"}}"
        ),
        started,
    );
}

fn projector(ladder: &Ladder) -> Vec<f64> {
    // Fixed maker/player station, exactly a grid node on every rung.
    let node = strike_node(&ladder.mesh);
    ladder
        .selected
        .iter()
        .map(|&m| -ladder.modes[m].phi[3 * node + 2])
        .collect()
}

/// SPD comparison floor: lambda_min(C) >= 1/||S||_infinity for symmetric
/// compliance S in the work-conjugate Mandel basis. An isotropic comparison
/// half-space with C=floor*I has nu=0, E=floor. Variational monotonicity then
/// bounds local half-space compliance by 1/floor. This is NOT a proof about
/// finite-body truncation, head viscoelasticity, or the actual Iowa contact.
fn local_wood_compliance_ratio() -> f64 {
    let mut ratio: f64 = 0.0;
    let head_max = hertz_law(1.39e6, 0.499).reduced_modulus;
    for pack in [HONDURAS, AFRICAN, WALNUT] {
        let m = law(pack);
        let [a, b, c] = m.e;
        let [ab, ac, bc] = m.nu;
        let norm = [
            (1.0 + ab.abs() + ac.abs()) / a,
            (1.0 + bc.abs()) / b + ab.abs() / a,
            1.0 / c + ac.abs() / a + bc.abs() / b,
            1.0 / (2.0 * m.g[0]),
            1.0 / (2.0 * m.g[1]),
            1.0 / (2.0 * m.g[2]),
        ]
        .into_iter()
        .fold(0.0, f64::max);
        ratio = ratio.max(head_max * norm);
    }
    ratio
}

fn render(
    ladder: &Ladder,
    transfer: [C64; 5],
    substeps: u32,
    local_compliance_ratio: f64,
) -> (Vec<f64>, f64) {
    render_with(ladder, transfer, substeps, local_compliance_ratio, false)
}

/// `rosewood_loss` selects the measured rosewood-bar decay law for the five
/// elastic modes; `false` is the explicit lossless scenario.
fn render_with(
    ladder: &Ladder,
    transfer: [C64; 5],
    substeps: u32,
    local_compliance_ratio: f64,
    rosewood_loss: bool,
) -> (Vec<f64>, f64) {
    let mut modes: Vec<_> = (0..5)
        .map(|mode| ModalAcousticMode {
            angular_frequency_rad_s: TAU * ladder.rungs[3].frequencies[mode],
            // Lossless unless the caller selects the sourced rosewood law; head
            // loss stays UNKNOWN either way. Do not promote the zero to measured
            // damping or qualified decay.
            damping_ratio: if rosewood_loss {
                modal_damping_ratio(ladder.rungs[3].frequencies[mode])
            } else {
                0.0
            },
            pressure_per_modal_velocity: transfer[mode],
        })
        .collect();
    modes.extend(ladder.rigid_projectors.iter().map(|_| ModalAcousticMode {
        angular_frequency_rad_s: 0.0,
        damping_ratio: 0.0,
        pressure_per_modal_velocity: C64::ZERO,
    }));
    let mut body = ModalAcousticTimeModel::try_new_with_free_coordinates(
        RATE * substeps,
        modes,
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let mut phi = projector(ladder);
    phi.extend_from_slice(&ladder.rigid_projectors);
    let mut workspace = ModalAcousticWorkspace::new(&body);
    let mut forces = vec![0.0; phi.len()];
    let mut law = hertz_law(HEAD_E, HEAD_NU);
    law.reduced_modulus /= 1.0 + local_compliance_ratio;
    let interface = hertz_interface();
    let h = 1.0 / f64::from(RATE * substeps);
    let (mut position, mut speed) = (0.0, STROKE);
    let mut detached = false;
    let mut contact_time = 0.0;
    let mut samples = Vec::with_capacity(6 * RATE as usize);
    for sample in 0..6 * RATE as usize {
        let mut pressure = 0.0;
        for substep in 0..substeps {
            let displacement: f64 = body
                .states()
                .iter()
                .zip(&phi)
                .map(|(s, p)| s.displacement_m_sqrt_kg * p)
                .sum();
            let overlap = position - displacement;
            if !detached && overlap < 0.0 && speed < 0.0 {
                detached = true;
                contact_time = (sample as f64 * f64::from(substeps) + f64::from(substep)) * h;
            }
            let load = if detached {
                0.0
            } else {
                force(law, &interface, overlap)
            };
            for (f, p) in forces.iter_mut().zip(&phi) {
                *f = p * load;
            }
            let frame = body
                .step_into(&forces, &mut workspace)
                .expect("material strike on 3-D modes");
            pressure = frame.observer_pressure_pa;
            speed -= h * load / HEAD_MASS;
            position += h * speed;
        }
        samples.push(pressure);
    }
    assert!(detached && contact_time > 0.0);
    (samples, contact_time)
}

fn artifact_dir() -> Option<PathBuf> {
    std::env::var_os("BAR_S1_OUT_DIR").map(PathBuf::from)
}

fn write_artifacts(
    out: &Path,
    ladder: &Ladder,
    transfer: [C64; 5],
    signal: &[f64],
    contact_time: f64,
) {
    std::fs::create_dir_all(out).unwrap();
    // A fixed physical pressure scale, never inferred from the signal peak.
    let (wav, clipped) = fs_couple::pcm_wav::encode_pcm16_wav(signal, RATE, FULL_SCALE_PA).unwrap();
    assert_eq!(
        clipped, 0,
        "increase declared physical range explicitly if needed"
    );
    assert!(
        wav[44..].iter().any(|b| *b != 0),
        "pressure must survive fixed-scale PCM quantization"
    );
    let hash = fs_blake3::hash_domain("org.frankensim.music-render.wav.v1", &wav).to_hex();
    std::fs::write(out.join("bar-s1-after.wav"), wav).unwrap();
    // The bare near-null radiator is quiet: report how many PCM16 steps the
    // signal spans, and keep the same pressure history at float32 precision.
    let lsb = FULL_SCALE_PA / f64::from(i16::MAX);
    let rms = (signal.iter().map(|p| p * p).sum::<f64>() / signal.len() as f64).sqrt();
    let peak = signal.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
    println!(
        "{{\"suite\":\"bar-s1-resolution\",\"pcm16_lsb_pa\":{lsb},\"rms_pa\":{rms},\"peak_pa\":{peak},\"rms_lsb\":{},\"peak_lsb\":{}}}",
        rms / lsb,
        peak / lsb
    );
    let f32_wav = fs_couple::pcm_wav::high_resolution::encode_pressure_wav(
        signal,
        RATE,
        1,
        FULL_SCALE_PA,
        fs_couple::pcm_wav::high_resolution::Encoding::Float32,
    )
    .unwrap();
    std::fs::write(out.join("bar-s1-after-f32.wav"), f32_wav).unwrap();
    let finest = &ladder.rungs[3];
    let est = tet10_estimate(&ladder.rungs, 0);
    let p: Vec<_> = transfer.iter().map(|p| [p.re, p.im]).collect();
    let band = (0..2).all(|mode| {
        let e = tet10_estimate(&ladder.rungs, mode);
        e.asymptotic
            && 1200.0 * det::ln(1.0 + e.uncertainty_hz / finest.frequencies[mode]) / det::ln(2.0)
                <= 1.0
    });
    let metadata = format!(
        "{{\"bar_id\":\"B-undercut\",\"rung_h_m\":{},\"n_tets\":{},\"f_modes_hz\":{:?},\"observed_order\":{:?},\"panel_count\":{},\"p_receiver_pa\":{p:?},\"contact_time_s\":{contact_time},\"wav_blake3\":\"{hash}\",\"sample_rate_hz\":{RATE},\"full_scale_pa\":{FULL_SCALE_PA},\"clipped_samples\":0,\"rigid_contact_coordinates\":6,\"direct_rigid_acceleration_radiation\":\"not modelled\",\"grader_json_path\":\"data/instruments/grades/bars/bar-s1-after.json\",\"before_db\":28.4,\"after_db\":null,\"targeted_ratio\":3,\"targeted_et_note\":\"Db5\",\"cut_depth_fraction\":{},\"uniform_geometry_scale\":{},\"same_specimen\":false,\"element\":\"TET10\",\"numerical_one_cent_band\":{},\"band_refusal_route\":\"finer TET10 rung\",\"wood_loss\":\"unknown:M8-W2c\",\"decay_prediction\":\"unknown\",\"radiation_feedback\":false,\"supports\":false,\"resonator\":false,\"mallet\":\"ideal elastic silicone sphere; no shaft; unknown rate loss\",\"material_scenario\":\"Indian rosewood E_L with bounded correlated hardwood transverse-ratio scenarios\",\"receiver_m\":{RECEIVER:?},\"observer\":\"narrow-band modal transfer at five retained elastic poles\",\"acceptance\":false}}\n",
        finest.h,
        finest.n_tets,
        finest.frequencies,
        est.orders,
        radiation().panel_counts.last().unwrap(),
        ladder.design.depth,
        ladder.design.scale,
        band
    );
    std::fs::write(out.join("bar-s1-render.json"), metadata).unwrap();
}

#[test]
fn bar_s1_material_strike_radiated_render_replays_and_time_axis_is_live() {
    let _guard = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let started = Instant::now();
    let ladder = undercut();
    let transfer = *radiation().pressures.last().unwrap();
    let mut last = None;
    let mut differences = Vec::new();
    for substeps in SUBSTEPS {
        let (signal, contact) = render(ladder, transfer, substeps, 0.0);
        if let Some((_, previous)) = &last {
            let energy: f64 = signal.iter().map(|p| p * p).sum();
            let error: f64 = signal
                .iter()
                .zip(previous)
                .map(|(a, b): (&f64, &f64)| (a - b).powi(2))
                .sum();
            let relative = (error / energy).sqrt();
            assert!(
                relative.is_finite() && relative > 0.0,
                "time axis must be live"
            );
            differences.push(relative);
        }
        row(
            "coupled-mallet-modal-pressure-render",
            "time-substep",
            "pass",
            &format!(
                "{{\"substeps\":{substeps},\"samples\":{},\"contact_time_s\":{contact},\"decay\":\"unknown; lossless diagnostic scenario\",\"acceptance\":false}}",
                signal.len()
            ),
            started,
        );
        last = Some((contact, signal));
    }
    let orders: Vec<f64> = differences
        .windows(2)
        .map(|d| det::ln(d[0] / d[1]) / det::ln(2.0))
        .collect();
    let finest_order = *orders.last().unwrap();
    assert!(
        differences.windows(2).all(|d| d[1] < d[0])
            && *differences.last().unwrap() < 0.01
            && (0.7..=1.3).contains(&finest_order),
        "strike time ladder needs refinement: {differences:?}, orders {orders:?} vs formal 1"
    );
    let (contact, signal) = last.unwrap();
    let finest = *SUBSTEPS.last().unwrap();
    let (replay, replay_contact) = render(ladder, transfer, finest, 0.0);
    assert_eq!(contact.to_bits(), replay_contact.to_bits());
    assert!(
        signal
            .iter()
            .zip(&replay)
            .all(|(a, b)| a.to_bits() == b.to_bits()),
        "same-host exact replay"
    );
    assert!(
        signal.iter().any(|v| v.abs() > 1e-7),
        "radiation is an actual pressure observation"
    );
    let local_ratio = local_wood_compliance_ratio();
    assert!(
        local_ratio < 0.01,
        "local half-space rigidity needs a finer contact model"
    );
    let change = |other: &[f64]| {
        (signal
            .iter()
            .zip(other)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            / signal.iter().map(|p| p * p).sum::<f64>())
        .sqrt()
    };
    let (softest, _) = render(ladder, transfer, finest, local_ratio);
    let contact_material_effect = change(&softest);
    // The rigid-wood idealization is admitted only while its bound stays an
    // order of magnitude under the receiver change that the DECLARED head
    // input envelope (Johnston E spread x nu interval) already carries. An
    // absolute waveform band would demand more of this omission than of the
    // inputs: the lossless six-second receiver is phase-sensitive at 2.6 kHz.
    let nominal = hertz_law(HEAD_E, HEAD_NU).reduced_modulus;
    let head_envelope_effects = [(1.25e6, 0.45), (1.39e6, 0.499)].map(|(e, nu)| {
        let softening = nominal / hertz_law(e, nu).reduced_modulus - 1.0;
        change(&render(ladder, transfer, finest, softening).0)
    });
    let envelope = head_envelope_effects[0].min(head_envelope_effects[1]);
    let dominance = contact_material_effect / envelope;
    // Liveness: ten times the compliance bound must fail the same gate.
    let planted = change(&render(ladder, transfer, finest, 10.0 * local_ratio).0) / envelope;
    assert!(planted > 0.1, "compliance gate is not live: {planted}");
    let verdict = if dominance <= 0.1 { "pass" } else { "fail" };
    row(
        "local-contact-compliance-endpoint",
        "material",
        verdict,
        &format!(
            "{{\"local_half_space_compliance_ratio_bound\":{local_ratio},\"receiver_relative_l2_endpoint_change\":{contact_material_effect},\"head_envelope_receiver_relative_l2_changes\":{head_envelope_effects:?},\"omission_to_input_envelope_ratio\":{dominance},\"tolerance_ratio\":0.1,\"planted_10x_ratio\":{planted},\"finite_body_contact_qualified\":false}}"
        ),
        started,
    );
    assert!(
        dominance <= 0.1,
        "finite-head/wood compliance endpoint needs a fuller contact model: {dominance}"
    );
    if let Some(out) = artifact_dir() {
        write_artifacts(&out, ladder, transfer, &signal, contact);
    }
    row(
        "render-replay-and-time-convergence",
        "time",
        "pass",
        &format!(
            "{{\"substeps\":{SUBSTEPS:?},\"relative_l2_differences\":{differences:?},\"observed_orders\":{orders:?},\"formal_order\":1,\"order_band\":[0.7,1.3],\"tolerance_relative_l2\":0.01,\"replay_bits_equal\":true,\"contact_time_s\":{contact},\"acceptance\":false}}"
        ),
        started,
    );
}

/// Same strike, same radiation, but the five elastic modes decay by the
/// measured rosewood law. Float32 output only: the bar radiates at a few LSB of
/// the fixed PCM16 scale used by the lossless render.
#[test]
fn bar_s1_rosewood_loss_scenario_render_replays_and_decays() {
    let _guard = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let started = Instant::now();
    let ladder = undercut();
    let transfer = *radiation().pressures.last().unwrap();
    let finest = *SUBSTEPS.last().unwrap();
    let (signal, contact) = render_with(ladder, transfer, finest, 0.0, true);
    let (replay, replay_contact) = render_with(ladder, transfer, finest, 0.0, true);
    let (coarser, _) = render_with(ladder, transfer, SUBSTEPS[2], 0.0, true);
    let energy: f64 = signal.iter().map(|p| p * p).sum();
    let change: f64 = signal
        .iter()
        .zip(&coarser)
        .map(|(a, b)| (a - b) * (a - b))
        .sum();
    let time_axis = (change / energy).sqrt();
    let rate = RATE as usize;
    let rms = |from: usize, to: usize| {
        (signal[from..to].iter().map(|p| p * p).sum::<f64>() / (to - from) as f64).sqrt()
    };
    let (early, late) = (rms(0, rate / 4), rms(rate, signal.len()));
    let frequencies = ladder.rungs[3].frequencies;
    let zetas: Vec<f64> = frequencies
        .iter()
        .map(|&f| modal_damping_ratio(f))
        .collect();
    let rates: Vec<f64> = frequencies
        .iter()
        .zip(&zetas)
        .map(|(&f, &zeta)| zeta * TAU * f)
        .collect();
    if let Some(out) = artifact_dir() {
        std::fs::create_dir_all(&out).unwrap();
        let wav = fs_couple::pcm_wav::high_resolution::encode_pressure_wav(
            &signal,
            RATE,
            1,
            FULL_SCALE_PA,
            fs_couple::pcm_wav::high_resolution::Encoding::Float32,
        )
        .unwrap();
        let hash = fs_blake3::hash_domain("org.frankensim.music-render.wav.v1", &wav).to_hex();
        std::fs::write(out.join("bar-s1-loss-f32.wav"), wav).unwrap();
        let cache =
            std::env::var("BAR_S1_DESIGN_CACHE").map_or("null".to_owned(), |p| format!("\"{p}\""));
        let metadata = format!(
            "{{\"bar_id\":\"B-undercut\",\"f_modes_hz\":{frequencies:?},\"damping_ratio\":{zetas:?},\"amplitude_decay_rate_per_s\":{rates:?},\"loss_law\":\"1/tau_d = a0 + a2 f^2, a0 = gamma_B/2, a2 = 2 pi^2 eta; Chaigne and Doutaut 1997 JASA 101(1) 539 Table I rosewood bar, Eqs. 18-19\",\"loss_is_scenario\":\"constants measured on another rosewood bar (E 21.3 GPa, 1015 kg/m3); mallet head loss unknown\",\"design_depth\":{},\"design_scale\":{},\"design_cache\":{cache},\"contact_time_s\":{contact},\"time_axis_relative_l2_substeps_8_to_16\":{time_axis},\"early_rms_pa_first_quarter_s\":{early},\"late_rms_pa_after_1_s\":{late},\"wav_blake3\":\"{hash}\",\"sample_rate_hz\":{RATE},\"full_scale_pa\":{FULL_SCALE_PA},\"encoding\":\"float32\",\"receiver_m\":{RECEIVER:?},\"observer\":\"narrow-band modal transfer at five retained elastic poles, Burton-Miller\",\"supports\":false,\"resonator\":false,\"same_specimen\":false,\"acceptance\":false}}\n",
            ladder.design.depth, ladder.design.scale
        );
        std::fs::write(out.join("bar-s1-loss-render.json"), metadata).unwrap();
    }
    row(
        "rosewood-loss-scenario-render",
        "loss",
        "pass",
        &format!(
            "{{\"damping_ratio\":{zetas:?},\"amplitude_decay_rate_per_s\":{rates:?},\"time_axis_relative_l2\":{time_axis},\"early_rms_pa\":{early},\"late_rms_pa\":{late},\"contact_time_s\":{contact},\"acceptance\":false}}"
        ),
        started,
    );
    assert_eq!(contact.to_bits(), replay_contact.to_bits());
    assert!(
        signal
            .iter()
            .zip(&replay)
            .all(|(a, b)| a.to_bits() == b.to_bits()),
        "same-host exact replay"
    );
    assert!(
        time_axis < 0.01,
        "strike time step needs refinement: {time_axis}"
    );
    assert!(
        early > 0.0 && late < 1e-3 * early,
        "loss law must damp the note"
    );
}
