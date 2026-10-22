//! IRP M1.2-W1c (pexfp.49): thin-wall model-form cross-check.
//!
//! A free-free, quarter-sawn Sitka spruce plate, solved three ways:
//! - **Material.** L along x, R along y, T through the thickness. The full
//!   tensor comes from the FPL-GTR-282 seed pack through the M0.2-W1
//!   reference-ratio adapter; no modulus is authored here.
//! - **3-D reference.** TET10 on an eighth model: x and y mirror planes, one
//!   symmetry class at a time, and the mid-plane with `u_x = u_y = 0`. That
//!   last condition keeps exactly the flexural subspace of this z-symmetric
//!   plate. It is an h-ladder with 2, 4 and 8 element layers through the full
//!   thickness (1, 2 and 4 in the half) and paired in-plane refinement; the
//!   3-D factorizations use geometric nested dissection (see [`RUNGS_3D`]).
//!   A true through-thickness p-FEM does not exist yet, so this is an
//!   h-ladder of quadratic elements.
//! - **Shell.** The MITC3+ shear-deformable shell, with κ = 5/6 as input data.
//! - **Kirchhoff.** The DKT plate.
//!
//! Modes are paired by the MAC of the mid-surface `w` on the quarter's coarse
//! vertex grid, within one symmetry class; frequencies never pair them.
//! In-plane modes (`w` energy share < 1/5) are excluded from the comparison
//! on every model. A mixed mode refuses.

use fs_evidence::ValidityDomain;
use fs_exec::Cx;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::OrthotropicElastic;
use fs_material::state_point::*;
use fs_math::det;
use fs_modal::{ModePair, SliceOptions, WindowPlan, slice_window, slice_windows};
use fs_plate::shell::mitc::{
    MitcBudget, NodalDirectors, ShearShellSection, ShellMass, assemble_mitc3p,
};
use fs_plate::shell::survey::material_axis_angle;
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh, PlateSection, ShellMesh};
use fs_qty::{Density, Dims, Pressure};
use fs_solid::linear3::tet10::{QuadraticTetProblem, TetQuadrature, promote_p1_to_p2};
use fs_solid::linear3::{TetAssemblyBudget, TetElasticMaterial, TetMaterialField};
use fs_sparse::{DirectOrdering, NdOptions};
use std::collections::BTreeMap;
use std::sync::{Mutex, OnceLock};
use std::time::Instant;

const A: f64 = 0.16; // extent along x (L) [m]
const B: f64 = 0.10; // extent along y (R) [m]
const H: f64 = 0.003; // thickness (T) [m]
const TAU: f64 = core::f64::consts::TAU;
/// The compared band [Hz].
const BAND: (f64, f64) = (2000.0, 8000.0);
/// The solved band [Hz]: wider, so modes at the band edges find partners.
/// Starts below the first elastic mode, so the DKT 10 c crossing is found
/// wherever it lies (a 1.5 kHz start reported the first band mode instead).
const SOLVE: (f64, f64) = (20.0, 9500.0);
/// Model-form tolerance for the shear shell against the 3-D reference.
const TOL: f64 = 0.01;
const MAC_PAIR: f64 = 0.9;
/// Below this `w` energy share a mode is in-plane; above `1 − ` it is flexural.
const IN_PLANE_SHARE: f64 = 0.2;
/// Eighth-model cells at the coarsest 3-D rung (x, y, half-thickness); the
/// in-plane counts are also the MAC sample grid.
const BASE: [usize; 3] = [6, 8, 1];
/// Paired refinement factors: 2, 4 and 8 element layers through the full
/// thickness, as the bead requires. With AMD the 8-layer rung (about 78k DOF
/// per class) did not finish within 20 minutes; measured 2026-10-22, the x3
/// rung alone took 87.6 s per class (34k DOF, 3 factorizations, nnz(L) 19.5M),
/// all of it in the factorizations. Geometric nested dissection on the TET10
/// node coordinates makes the x4 rung affordable (see [`solid_modes`]). The
/// order uses the actual h ratios.
const RUNGS_3D: [usize; 3] = [1, 2, 4];
/// In-plane refinement factors for the thin isotropic sanity plate.
/// Measured 2026-10-22: x2/x3/x4 left the 3-D reference pre-asymptotic
/// (orders 2.2-2.4, last increment up to 0.22 %, spread 0.274 % > 0.2 %).
const RUNGS_THIN: [usize; 3] = [4, 6, 8];
/// Full-plate shell cells are `2·BASE·r`, so shell nodes contain the samples.
const RUNGS_SHELL: [usize; 3] = [2, 4, 8];
const KAPPA: f64 = 5.0 / 6.0;
/// The spruce pack's only validity state.
const MC_PCT: f64 = 12.0;

const SPRUCE: &str =
    include_str!("../../../data/matdb/seed-v1/spruce-sitka-fpl-gtr282/properties.tsv");
const ALUMINIUM: &str =
    include_str!("../../../data/matdb/seed-v1/aluminum-2024-t3-sheet-mil-hdbk-5j/properties.tsv");

static SERIAL: Mutex<()> = Mutex::new(());
static SPRUCE_STUDY: OnceLock<Study> = OnceLock::new();

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 49,
                kernel_id: 49,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx)
    })
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

fn window(band: (f64, f64)) -> (f64, f64) {
    ((TAU * band.0).powi(2), (TAU * band.1).powi(2))
}

/// One window per 3-D solve: measured, the cost is factorization-bound, and
/// split windows only add endpoint factorizations (12.9 s vs 5.2 s at 10k DOF).
const SUB_WINDOWS: usize = 1;

fn sub_windows(band: (f64, f64), count: usize) -> Vec<(f64, f64)> {
    let ratio = det::ln(band.1 / band.0) / count as f64;
    let edge = |i: usize| match i {
        0 => band.0,
        i if i == count => band.1,
        i => band.0 * det::exp(ratio * i as f64),
    };
    (0..count).map(|i| window((edge(i), edge(i + 1)))).collect()
}

fn hz(m: &ModePair) -> f64 {
    m.lambda.max(0.0).sqrt() / TAU
}

fn cents(f: f64, reference: f64) -> f64 {
    1200.0 * det::ln(f / reference) / core::f64::consts::LN_2
}

/// A receipt number: `null` when non-finite, since bare `NaN` is not JSON.
fn json_f64(x: f64, digits: usize) -> String {
    if x.is_finite() {
        format!("{x:.digits$}")
    } else {
        "null".to_string()
    }
}

// ---------------------------------------------------------------------------
// Material from the seed pack
// ---------------------------------------------------------------------------

/// Every `scalar` row of a seed pack, with its validity rows, as claims.
fn pack_card(tsv: &str, label: &str) -> MaterialCard {
    let mut validity: BTreeMap<&str, (&str, f64, f64)> = BTreeMap::new();
    for line in tsv.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.first() == Some(&"validity") {
            validity.insert(f[1], (f[2], f[3].parse().unwrap(), f[4].parse().unwrap()));
        }
    }
    let mut claims = ClaimSet::new();
    for line in tsv.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.first() != Some(&"scalar") {
            continue;
        }
        let raw: f64 = f[4].parse().expect("source scalar");
        let (dims, value) = match f[5] {
            "1" => (Dims::NONE, raw),
            "MPa" => (Pressure::DIMS, raw * 1e6),
            "GPa" => (Pressure::DIMS, raw * 1e9),
            "kg/m3" => (Density::DIMS, raw),
            unit => panic!("unexpected source unit {unit} on {}", f[3]),
        };
        let mut domain = ValidityDomain::unconstrained();
        if let Some(&(axis, lo, hi)) = validity.get(f[1]) {
            domain = domain.with(axis, lo, hi);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(f[3], dims),
                value: PropertyValue::Scalar { value, dims },
                validity: domain,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: format!("{label} seed-v1 row {}", f[1]),
                    license: "Work-of-the-US-Government-Public-Use-Permitted".into(),
                    artifact: None,
                },
            })
            .expect("pack claim");
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: label.into(),
            phase: "clear wood".into(),
            process: "seed-v1 pack rows".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .expect("pack card")
}

struct Spruce {
    state: OrthotropicElasticStatePoint,
    law: OrthotropicElastic,
    rho: f64,
}

fn spruce() -> Spruce {
    let card = pack_card(SPRUCE, "spruce-sitka-fpl-gtr282");
    let resolved = resolve_orthotropic_elastic_state_point_from_reference_ratios(
        &card,
        &QueryPoint::new().with("moisture_content", MC_PCT).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        1e-3, // declared small-strain domain, not a measured failure strain
        OrthotropicRequirementSet::Full3d,
        PoissonReciprocityRule::MajorRatios,
    )
    .expect("full 3-D tensor from the seed pack");
    let rho = resolved.density_kg_m3();
    let law = resolved.law().expect("full 3-D law").clone();
    let ReferenceRatioOrthotropicStatePoint::Full3d(state) = resolved else {
        panic!("a Full3d request returned a plane state");
    };
    Spruce { state, law, rho }
}

/// The pack's isotropic scalars, read directly (sanity fixture only).
fn scalar(tsv: &str, property: &str) -> f64 {
    tsv.lines()
        .find_map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (f.first() == Some(&"scalar") && f.get(3) == Some(&property)).then(|| {
                let v: f64 = f[4].parse().unwrap();
                if f[5] == "GPa" { v * 1e9 } else { v }
            })
        })
        .unwrap_or_else(|| panic!("source property missing: {property}"))
}

// ---------------------------------------------------------------------------
// Symmetry classes and mode records
// ---------------------------------------------------------------------------

/// Parity of `w` across a mirror plane through the plate centre.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Parity {
    Even,
    Odd,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Class {
    x: Parity,
    y: Parity,
}

const CLASSES: [Class; 4] = [
    Class {
        x: Parity::Even,
        y: Parity::Even,
    },
    Class {
        x: Parity::Even,
        y: Parity::Odd,
    },
    Class {
        x: Parity::Odd,
        y: Parity::Even,
    },
    Class {
        x: Parity::Odd,
        y: Parity::Odd,
    },
];

impl Class {
    fn label(self) -> &'static str {
        match (self.x, self.y) {
            (Parity::Even, Parity::Even) => "EE",
            (Parity::Even, Parity::Odd) => "EO",
            (Parity::Odd, Parity::Even) => "OE",
            (Parity::Odd, Parity::Odd) => "OO",
        }
    }
}

/// One mode as the comparison sees it.
#[derive(Clone, Debug)]
struct Mode {
    hz: f64,
    /// Share of the translational energy in `w`.
    w_share: f64,
    /// Mid-surface `w` at the quarter's coarse vertex grid.
    samples: Vec<f64>,
}

#[derive(Debug, PartialEq, Eq)]
enum PairRefusal {
    Unpaired,
    Ambiguous(Vec<usize>),
}

fn mac(a: &[f64], b: &[f64]) -> f64 {
    let dot = |u: &[f64], v: &[f64]| u.iter().zip(v).map(|(x, y)| x * y).sum::<f64>();
    let ab = dot(a, b);
    ab * ab / (dot(a, a) * dot(b, b))
}

/// The unique candidate whose MAC with `target` reaches [`MAC_PAIR`].
fn pair(target: &[f64], candidates: &[&[f64]]) -> Result<usize, PairRefusal> {
    let hits: Vec<usize> = candidates
        .iter()
        .enumerate()
        .filter(|(_, c)| mac(target, c) >= MAC_PAIR)
        .map(|(i, _)| i)
        .collect();
    match hits.as_slice() {
        [i] => Ok(*i),
        [] => Err(PairRefusal::Unpaired),
        _ => Err(PairRefusal::Ambiguous(hits)),
    }
}

/// Share of `v` in span{x, y}: `vᵀPv / vᵀv` for the orthogonal projector P.
fn span_fraction(v: &[f64], x: &[f64], y: &[f64]) -> f64 {
    let dot = |u: &[f64], w: &[f64]| u.iter().zip(w).map(|(a, b)| a * b).sum::<f64>();
    let (xx, xy, yy, xv, yv) = (dot(x, x), dot(x, y), dot(y, y), dot(x, v), dot(y, v));
    let det = xx * yy - xy * xy;
    if !(det > 0.0) {
        return 0.0;
    }
    let (a, b) = ((yy * xv - xy * yv) / det, (xx * yv - xy * xv) / det);
    (a * xv + b * yv) / dot(v, v)
}

/// Frequency of `target`'s partner on another rung of one model.
/// - **Unique.** A unique MAC match.
/// - **Two-mode cluster.** A coarser rung rotates a near-degenerate pair into
///   a different basis (measured: MACs 0.50 and 0.55). The target and exactly
///   one other fine mode must both lie at least [`MAC_PAIR`] in the span of the
///   same two coarse modes; within that certified cluster they pair by
///   frequency rank. Global sorted order is never used.
fn track(target: usize, fine: &[Mode], coarse: &[Mode]) -> Result<f64, PairRefusal> {
    let candidates: Vec<&[f64]> = coarse.iter().map(|c| c.samples.as_slice()).collect();
    match pair(&fine[target].samples, &candidates) {
        Ok(k) => return Ok(coarse[k].hz),
        Err(PairRefusal::Ambiguous(hits)) => return Err(PairRefusal::Ambiguous(hits)),
        Err(PairRefusal::Unpaired) => {}
    }
    let mut cluster = None;
    for a in 0..coarse.len() {
        for b in a + 1..coarse.len() {
            let (x, y) = (&coarse[a].samples, &coarse[b].samples);
            if span_fraction(&fine[target].samples, x, y) < MAC_PAIR {
                continue;
            }
            let partners: Vec<usize> = (0..fine.len())
                .filter(|&j| j != target && span_fraction(&fine[j].samples, x, y) >= MAC_PAIR)
                .collect();
            if let [j] = partners[..] {
                if cluster.is_some() {
                    return Err(PairRefusal::Ambiguous(vec![a, b]));
                }
                cluster = Some((a, b, j));
            }
        }
    }
    let (a, b, j) = cluster.ok_or(PairRefusal::Unpaired)?;
    let (lo, hi) = if coarse[a].hz <= coarse[b].hz {
        (coarse[a].hz, coarse[b].hz)
    } else {
        (coarse[b].hz, coarse[a].hz)
    };
    Ok(if fine[target].hz <= fine[j].hz {
        lo
    } else {
        hi
    })
}

fn flexural(modes: &[Mode]) -> Vec<Mode> {
    for m in modes {
        assert!(
            m.w_share < IN_PLANE_SHARE || m.w_share > 1.0 - IN_PLANE_SHARE,
            "mixed in-plane/flexural mode at {} Hz (w share {})",
            m.hz,
            m.w_share
        );
    }
    modes
        .iter()
        .filter(|m| m.w_share > 1.0 - IN_PLANE_SHARE)
        .cloned()
        .collect()
}

// ---------------------------------------------------------------------------
// 3-D TET10 quarter model
// ---------------------------------------------------------------------------

struct Quarter {
    nodes: Vec<[f64; 3]>,
    tets: Vec<[usize; 4]>,
    counts: [usize; 3],
}

fn quarter(counts: [usize; 3], h: f64) -> Quarter {
    let (complex, nodes) = fs_conduction::fixtures::structured_tets(counts, &|i, j, k| {
        [
            i as f64 * 0.5 * A / counts[0] as f64,
            j as f64 * 0.5 * B / counts[1] as f64,
            -0.5 * h + k as f64 * h / counts[2] as f64,
        ]
    });
    Quarter {
        nodes,
        tets: complex.tets.iter().map(|t| t.map(|i| i as usize)).collect(),
        counts,
    }
}

/// The half-thickness eighth, z in [0, h/2], with the mid-plane at k = 0.
fn eighth(counts: [usize; 3], h: f64) -> Quarter {
    let (complex, nodes) = fs_conduction::fixtures::structured_tets(counts, &|i, j, k| {
        [
            i as f64 * 0.5 * A / counts[0] as f64,
            j as f64 * 0.5 * B / counts[1] as f64,
            k as f64 * 0.5 * h / counts[2] as f64,
        ]
    });
    Quarter {
        nodes,
        tets: complex.tets.iter().map(|t| t.map(|i| i as usize)).collect(),
        counts,
    }
}

/// P1 index of structured vertex (i, j, k), as `structured_tets` numbers it.
fn vertex(counts: [usize; 3], i: usize, j: usize, k: usize) -> usize {
    i * (counts[1] + 1) * (counts[2] + 1) + j * (counts[2] + 1) + k
}

/// The four quarters mirrored into the full plate, with exactly mirrored
/// coordinates; reflected tetrahedra are reoriented by the P2 promotion.
fn mirrored_full(q: &Quarter) -> Quarter {
    let [nx, ny, nz] = q.counts;
    let full = [2 * nx, 2 * ny, nz];
    let index = |i: isize, j: isize, k: usize| -> usize {
        vertex(
            full,
            (i + nx as isize) as usize,
            (j + ny as isize) as usize,
            k,
        )
    };
    let mut nodes = vec![[0.0; 3]; (2 * nx + 1) * (2 * ny + 1) * (nz + 1)];
    for i in 0..=nx {
        for j in 0..=ny {
            for k in 0..=nz {
                let p = q.nodes[vertex(q.counts, i, j, k)];
                for (sx, sy) in [(1.0, 1.0), (-1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)] {
                    let (ii, jj) = (sx as isize * i as isize, sy as isize * j as isize);
                    nodes[index(ii, jj, k)] = [sx * p[0], sy * p[1], p[2]];
                }
            }
        }
    }
    // Recover each quarter vertex's grid index from its P1 number.
    let grid = |v: usize| {
        let k = v % (nz + 1);
        let j = (v / (nz + 1)) % (ny + 1);
        let i = v / ((nz + 1) * (ny + 1));
        (i as isize, j as isize, k)
    };
    let mut tets = Vec::with_capacity(4 * q.tets.len());
    for (sx, sy) in [(1, 1), (-1, 1), (1, -1), (-1, -1)] {
        for t in &q.tets {
            tets.push(t.map(|v| {
                let (i, j, k) = grid(v);
                index(sx * i, sy * j, k)
            }));
        }
    }
    Quarter {
        nodes,
        tets,
        counts: full,
    }
}

/// Homogeneous symmetry conditions of `class` on the planes x = 0 and y = 0,
/// plus `u_x = u_y = 0` on z = 0 when only the flexural subspace is kept.
fn symmetry_dofs(nodes: &[[f64; 3]], class: Class, flexural: bool) -> Vec<usize> {
    let mut fixed = Vec::new();
    for (n, p) in nodes.iter().enumerate() {
        if flexural && p[2] == 0.0 {
            fixed.extend([3 * n, 3 * n + 1]);
        }
        if p[0].abs() <= 1e-12 * A {
            match class.x {
                Parity::Even => fixed.push(3 * n),
                Parity::Odd => fixed.extend([3 * n + 1, 3 * n + 2]),
            }
        }
        if p[1].abs() <= 1e-12 * B {
            match class.y {
                Parity::Even => fixed.push(3 * n + 1),
                Parity::Odd => fixed.extend([3 * n, 3 * n + 2]),
            }
        }
    }
    fixed.sort_unstable();
    fixed.dedup();
    fixed
}

/// The certified modes of one 3-D window.
struct Solved {
    modes: Vec<Mode>,
    /// Sylvester inertia below the window: the rigid modes, plus any elastic
    /// mode under its lower edge.
    below: usize,
}

/// TET10 modes of `mesh` under `fixed`; `sample` maps the coarse sample grid
/// to P1 vertex numbers on the mid-plane.
fn solid_modes(
    mesh: &Quarter,
    material: &TetElasticMaterial,
    class: Option<Class>,
    flexural: bool,
    band: (f64, f64),
    windows: usize,
    sample: &[usize],
    cx: &Cx<'_>,
) -> Solved {
    let started = Instant::now();
    let (nodes, tets) = promote_p1_to_p2(&mesh.nodes, &mesh.tets, None, cx).expect("P2 mesh");
    let fixed = class.map_or_else(Vec::new, |c| symmetry_dofs(&nodes, c, flexural));
    let assembly = QuadraticTetProblem {
        nodes_m: &nodes,
        elements: &tets,
        materials: TetMaterialField::Uniform(material),
        fixed_dofs: &fixed,
        budget: TetAssemblyBudget::standard(),
        quadrature: TetQuadrature {
            stiffness_points_per_dir: 3,
            mass_points_per_dir: 4,
        },
    }
    .assemble(cx)
    .expect("TET10 assembly");
    let assembled_s = started.elapsed().as_secs_f64();
    // The cost is factorization-bound, so the ordering matters more than the
    // windows: nested dissection instead of AMD, on the same certified path.
    // Counts stay inertia-certified per window and the half-open windows tile
    // the band.
    let mut opts = options();
    opts.ordering = DirectOrdering::NestedDissection(NdOptions {
        leaf_size: 64,
        refinement_passes: 4,
    });
    let reports = if windows == 1 {
        vec![
            slice_window(&assembly.stiffness, &assembly.mass, window(band), &opts)
                .expect("certified 3-D window"),
        ]
    } else {
        let plan = WindowPlan {
            windows: sub_windows(band, windows),
            threads: 1,
            maximum_windows: windows,
            maximum_bytes: 64 << 30,
            maximum_relative_halfwidth: 1e-8,
            coordinates: None,
        };
        opts.ldlt.max_factor_bytes = Some(32 << 30);
        slice_windows(&assembly.stiffness, &assembly.mass, &plan, &opts, cx)
            .expect("certified 3-D windows")
            .reports
    };
    let modes: Vec<&ModePair> = reports.iter().flat_map(|r| r.modes.iter()).collect();
    let (factorizations, lanczos, nnz_l) = reports.iter().fold((0, 0, 0), |acc, r| {
        (
            acc.0 + r.stats.factorizations,
            acc.1 + r.stats.lanczos_iters,
            acc.2.max(r.stats.factor_nnz_l),
        )
    });
    eprintln!(
        "pexfp49 3-D {:?} class {} dofs {} modes {} factorizations {factorizations} lanczos {lanczos} nnz_l {nnz_l} assemble {assembled_s:.1}s total {:.1}s",
        mesh.counts,
        class.map_or("full", Class::label),
        assembly.free_dofs.len(),
        modes.len(),
        started.elapsed().as_secs_f64()
    );
    let vertices = mesh.nodes.len();
    let modes = modes
        .into_iter()
        .map(|m| {
            let mut u = vec![0.0; 3 * nodes.len()];
            for (r, &d) in assembly.free_dofs.iter().enumerate() {
                u[d] = m.phi[r];
            }
            let (mut w2, mut all2) = (0.0, 0.0);
            for n in 0..vertices {
                w2 += u[3 * n + 2] * u[3 * n + 2];
                all2 +=
                    u[3 * n] * u[3 * n] + u[3 * n + 1] * u[3 * n + 1] + u[3 * n + 2] * u[3 * n + 2];
            }
            Mode {
                hz: hz(m),
                w_share: w2 / all2,
                samples: sample.iter().map(|&n| u[3 * n + 2]).collect(),
            }
        })
        .collect();
    Solved {
        modes,
        below: reports.iter().map(|r| r.below_low).min().unwrap(),
    }
}

/// Coarse-grid mid-plane (k = 0) vertices of an eighth at refinement `r`.
fn mid_plane_samples(base: [usize; 3], counts: [usize; 3], r: usize) -> Vec<usize> {
    let mut out = Vec::new();
    for ic in 0..=base[0] {
        for jc in 0..=base[1] {
            out.push(vertex(counts, ic * r, jc * r, 0));
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Full-plate shells
// ---------------------------------------------------------------------------

struct PlateModes {
    class: Vec<Class>,
    modes: Vec<Mode>,
}

/// Full-plate grid of `2·BASE·r` cells and the sample nodes inside it.
fn plate_grid(base: [usize; 3], r: usize) -> (PlateMesh, Vec<usize>, [usize; 2]) {
    let (nx, ny) = (2 * base[0] * r, 2 * base[1] * r);
    let (qx, qy) = (nx / 2, ny / 2);
    // One quarter triangulation mirrored into all four quadrants, so the
    // discrete operator commutes with both reflections and every mode of a
    // simple eigenvalue lies in exactly one class (a one-diagonal grid mixes
    // near-degenerate modes of different classes).
    let (dx, dy) = (0.5 * A / qx as f64, 0.5 * B / qy as f64);
    let id = |i: usize, j: usize| j * (nx + 1) + i;
    let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1));
    for j in 0..=ny {
        for i in 0..=nx {
            nodes.push((
                (i as f64 - qx as f64) * dx + 0.5 * A,
                (j as f64 - qy as f64) * dy + 0.5 * B,
            ));
        }
    }
    let mut tris = Vec::with_capacity(2 * nx * ny);
    for (sx, sy) in [(1_isize, 1_isize), (-1, 1), (1, -1), (-1, -1)] {
        let at = |a: usize, b: usize| {
            id(
                (qx as isize + sx * a as isize) as usize,
                (qy as isize + sy * b as isize) as usize,
            )
        };
        for b in 0..qy {
            for a in 0..qx {
                let (p00, p10, p01, p11) = (at(a, b), at(a + 1, b), at(a, b + 1), at(a + 1, b + 1));
                for t in [[p00, p10, p11], [p00, p11, p01]] {
                    // A single reflection reverses orientation.
                    tris.push(if sx * sy > 0 { t } else { [t[0], t[2], t[1]] });
                }
            }
        }
    }
    let mesh = PlateMesh::from_unstructured(nodes, tris).expect("mirrored plate mesh");
    let mut samples = Vec::new();
    for ic in 0..=base[0] {
        for jc in 0..=base[1] {
            let (i, j) = (base[0] * r + ic * r, base[1] * r + jc * r);
            samples.push(j * (nx + 1) + i);
        }
    }
    (mesh, samples, [nx, ny])
}

/// Parity class of a full-plate `w` field on the structured grid.
fn classify(w: &[f64], [nx, ny]: [usize; 2]) -> Class {
    let at = |i: usize, j: usize| w[j * (nx + 1) + i];
    let norm: f64 = w.iter().map(|v| v * v).sum();
    let (mut px, mut py) = (0.0, 0.0);
    for j in 0..=ny {
        for i in 0..=nx {
            px += at(i, j) * at(nx - i, j);
            py += at(i, j) * at(i, ny - j);
        }
    }
    let parity = |p: f64| {
        let p = p / norm;
        assert!(p.abs() > 0.9, "mode of no clear parity ({p})");
        if p > 0.0 { Parity::Even } else { Parity::Odd }
    };
    Class {
        x: parity(px),
        y: parity(py),
    }
}

fn shell_modes(
    base: [usize; 3],
    r: usize,
    section: &ShearShellSection,
    band: (f64, f64),
) -> PlateModes {
    let (plate, samples, grid) = plate_grid(base, r);
    let mesh = ShellMesh::new(
        plate.nodes.iter().map(|&(x, y)| [x, y, 0.0]).collect(),
        plate.tris.clone(),
    )
    .expect("shell mesh");
    let n = mesh.nodes.len();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; n], 0.5).expect("directors");
    let sections: Vec<ShearShellSection> = (0..mesh.tris.len())
        .map(|e| {
            let angle = material_axis_angle(&mesh.facet(e).unwrap(), [1.0, 0.0, 0.0]).unwrap();
            section.with_material_angle_rad(angle).unwrap()
        })
        .collect();
    let model = assemble_mitc3p(
        &mesh,
        &directors,
        &sections,
        &[],
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: n,
            max_triangles: mesh.tris.len(),
        },
        1e-3,
        || false,
    )
    .expect("MITC3+ assembly");
    let pencil = model.pencil();
    let report =
        slice_window(&pencil.k, &pencil.m, window(band), &options()).expect("shell window");
    let mut out = PlateModes {
        class: Vec::new(),
        modes: Vec::new(),
    };
    for m in &report.modes {
        let w: Vec<f64> = (0..n).map(|v| m.phi[6 * v + 2]).collect();
        let w2: f64 = w.iter().map(|v| v * v).sum();
        let all2: f64 = (0..n)
            .map(|v| {
                (0..3)
                    .map(|c| m.phi[6 * v + c] * m.phi[6 * v + c])
                    .sum::<f64>()
            })
            .sum();
        let mode = Mode {
            hz: hz(m),
            w_share: w2 / all2,
            samples: samples.iter().map(|&v| w[v]).collect(),
        };
        if mode.w_share > 1.0 - IN_PLANE_SHARE {
            out.class.push(classify(&w, grid));
        } else {
            assert!(
                mode.w_share < IN_PLANE_SHARE,
                "mixed shell mode at {} Hz",
                mode.hz
            );
            continue;
        }
        out.modes.push(mode);
    }
    out
}

fn dkt_modes(base: [usize; 3], r: usize, section: &PlateSection, band: (f64, f64)) -> PlateModes {
    let (plate, samples, grid) = plate_grid(base, r);
    let model = fs_plate::assemble(
        &plate,
        section,
        &[],
        &[],
        &AssemblyOptions {
            pretension: 0.0,
            support: EdgeSupport::SimplySupported, // ignored: no boundary nodes
        },
    )
    .expect("DKT assembly");
    let report = fs_plate::modes(&model, window(band), &options()).expect("DKT window");
    let n = plate.nodes.len();
    let mut out = PlateModes {
        class: Vec::new(),
        modes: Vec::new(),
    };
    for m in &report.modes {
        let w: Vec<f64> = (0..n)
            .map(|v| model.dof_map[3 * v].map_or(0.0, |d| m.phi[d]))
            .collect();
        out.class.push(classify(&w, grid));
        out.modes.push(Mode {
            hz: hz(m),
            w_share: 1.0,
            samples: samples.iter().map(|&v| w[v]).collect(),
        });
    }
    out
}

// ---------------------------------------------------------------------------
// The spruce study (shared by the three spruce tests)
// ---------------------------------------------------------------------------

/// Observed order, reference and bound from three rungs with refinement
/// factors `factors` (h = 1/factor). The order solves
/// `(h1^p − h2^p)/(h2^p − h3^p) = (f1 − f2)/(f2 − f3)` by bisection; it is
/// NaN for a non-monotone ladder. Richardson applies only when the order is
/// asymptotic (3..=5); otherwise the declared bound is the last increment.
fn ladder(f: [f64; 3], factors: [usize; 3]) -> (f64, f64, f64) {
    let h = factors.map(|k| 1.0 / k as f64);
    let (d1, d2) = (f[0] - f[1], f[1] - f[2]);
    let pow = |x: f64, p: f64| fs_math::det::pow(x, p);
    let g = |p: f64| (pow(h[0], p) - pow(h[1], p)) / (pow(h[1], p) - pow(h[2], p));
    let ratio = d1 / d2;
    let order = if d1 > 0.0 && d2 > 0.0 && (g(0.05)..=g(12.0)).contains(&ratio) {
        let (mut lo, mut hi) = (0.05, 12.0);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if g(mid) < ratio {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    } else {
        f64::NAN
    };
    if (3.0..=5.0).contains(&order) {
        let tail = d2 * pow(h[2], order) / (pow(h[1], order) - pow(h[2], order));
        (order, f[2] - tail, tail.abs())
    } else {
        (order, f[2], d2.abs())
    }
}

/// One 3-D reference mode with its ladder.
#[derive(Clone, Debug)]
struct Reference {
    class: Class,
    /// Frequencies on the 3-D rungs, coarse to fine.
    rungs_hz: [f64; 3],
    order: f64,
    /// Reference frequency: Richardson when asymptotic, else the finest rung.
    hz: f64,
    /// Error bound on `hz`.
    bound_hz: f64,
    samples: Vec<f64>,
}

struct Study {
    references: Vec<Reference>,
    shell: Vec<PlateModes>,
    dkt: Vec<PlateModes>,
    elapsed_s: f64,
}

fn spruce_study() -> &'static Study {
    SPRUCE_STUDY.get_or_init(|| {
        let started = Instant::now();
        let s = spruce();
        let material = TetElasticMaterial::try_from_resolved_orthotropic_state(
            &s.state,
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            fs_blake3::hash_domain("org.frankensim.pexfp49.axes.v1", b"x=L,y=R,z=T"),
        )
        .expect("oriented spruce tensor");
        let mut references = Vec::new();
        with_cx(|cx| {
            for class in CLASSES {
                let ladders: Vec<Vec<Mode>> = RUNGS_3D
                    .iter()
                    .map(|&r| {
                        let mesh = eighth(BASE.map(|c| c * r), H);
                        let samples = mid_plane_samples(BASE, mesh.counts, r);
                        flexural(
                            &solid_modes(
                                &mesh,
                                &material,
                                Some(class),
                                true,
                                SOLVE,
                                SUB_WINDOWS,
                                &samples,
                                cx,
                            )
                            .modes,
                        )
                    })
                    .collect();
                let finest = &ladders[2];
                for (index, m) in finest.iter().enumerate() {
                    // Above the band a mode only partners the band edge: its own
                    // coarse partner may lie above SOLVE (measured: EE 9109 Hz).
                    // Modes below the band stay, for the DKT refusal band.
                    if m.hz > BAND.1 {
                        continue;
                    }
                    let mut rungs = [0.0; 3];
                    rungs[2] = m.hz;
                    for rung in 0..2 {
                        rungs[rung] = track(index, finest, &ladders[rung]).unwrap_or_else(|e| {
                            panic!(
                                "{} mode at {} Hz untracked on rung {rung}: {e:?}",
                                class.label(),
                                m.hz
                            )
                        });
                    }
                    let (order, reference, bound) = ladder(rungs, RUNGS_3D);
                    references.push(Reference {
                        class,
                        rungs_hz: rungs,
                        order,
                        hz: reference,
                        bound_hz: bound,
                        samples: m.samples.clone(),
                    });
                }
            }
        });
        let s_section = ShearShellSection::from_orthotropic_solid(
            s.law.e[0],
            s.law.e[1],
            s.law.nu[0],
            s.law.g[0],
            s.law.g[2],
            s.law.g[1],
            H,
            s.rho,
            0.0,
            KAPPA,
        )
        .expect("shear shell section");
        let d_section = PlateSection::orthotropic(&s.law, H, s.rho).expect("DKT section");
        let shell = RUNGS_SHELL
            .iter()
            .map(|&r| shell_modes(BASE, r, &s_section, SOLVE))
            .collect();
        let dkt = RUNGS_SHELL
            .iter()
            .map(|&r| dkt_modes(BASE, r, &d_section, SOLVE))
            .collect();
        Study {
            references,
            shell,
            dkt,
            elapsed_s: started.elapsed().as_secs_f64(),
        }
    })
}

fn in_band(f: f64) -> bool {
    (BAND.0..=BAND.1).contains(&f)
}

/// The partner of `reference` among a plate model's modes of the same class.
fn partner<'a>(reference: &Reference, model: &'a PlateModes) -> &'a Mode {
    try_partner(reference, model).unwrap_or_else(|(e, best)| {
        panic!(
            "{} reference at {} Hz has no plate partner: {e:?} (best MAC {best})",
            reference.class.label(),
            reference.hz
        )
    })
}

/// The partner, or the refusal together with the best MAC that was found.
fn try_partner<'a>(
    reference: &Reference,
    model: &'a PlateModes,
) -> Result<&'a Mode, (PairRefusal, f64)> {
    let (indices, candidates): (Vec<usize>, Vec<&[f64]>) = model
        .modes
        .iter()
        .enumerate()
        .filter(|(i, _)| model.class[*i] == reference.class)
        .map(|(i, m)| (i, m.samples.as_slice()))
        .unzip();
    match pair(&reference.samples, &candidates) {
        Ok(k) => Ok(&model.modes[indices[k]]),
        Err(e) => {
            let best = candidates
                .iter()
                .map(|c| mac(&reference.samples, c))
                .fold(0.0, f64::max);
            Err((e, best))
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn g0_mac_pairing_refuses_ambiguous_pairs() {
    let target = [1.0, 0.5, -0.25, 0.0];
    let near_a = [1.0, 0.52, -0.24, 0.01];
    let near_b = [0.98, 0.49, -0.26, -0.02];
    let other = [0.0, -0.25, 0.5, 1.0];
    assert_eq!(pair(&target, &[&near_a, &other]), Ok(0));
    assert_eq!(
        pair(&target, &[&near_a, &other, &near_b]),
        Err(PairRefusal::Ambiguous(vec![0, 2]))
    );
    assert_eq!(pair(&target, &[&other]), Err(PairRefusal::Unpaired));
    // A near-degenerate pair that a coarser rung rotates by 45°: no unique MAC
    // (each 0.5), so the certified two-mode cluster pairs by rank within it.
    let mode = |hz: f64, samples: Vec<f64>| Mode {
        hz,
        w_share: 1.0,
        samples,
    };
    let fine = [
        mode(100.0, vec![1.0, 0.0, 0.0]),
        mode(101.0, vec![0.0, 1.0, 0.0]),
        mode(300.0, vec![0.0, 0.0, 1.0]),
    ];
    let s = core::f64::consts::FRAC_1_SQRT_2;
    let coarse = [
        mode(110.0, vec![s, s, 0.0]),
        mode(108.0, vec![s, -s, 0.0]),
        mode(310.0, vec![0.0, 0.0, 1.0]),
    ];
    assert_eq!(track(0, &fine, &coarse), Ok(108.0));
    assert_eq!(track(1, &fine, &coarse), Ok(110.0));
    assert_eq!(track(2, &fine, &coarse), Ok(310.0));
    // A third fine mode in the same span makes the cluster unresolvable.
    let crowded = [
        mode(100.0, vec![1.0, 0.0, 0.0]),
        mode(101.0, vec![0.0, 1.0, 0.0]),
        mode(102.0, vec![s, s, 0.0]),
    ];
    assert!(track(0, &crowded, &coarse).is_err());
    println!(
        "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W1c\",\"test\":\"g0_mac_pairing_refuses_ambiguous_pairs\",\"tier\":\"G0\",\"verdict\":\"pass\"}}"
    );
}

#[test]
fn g0_symmetry_class_merge_recovers_full_model_spectrum() {
    let _serial = SERIAL
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let s = spruce();
    let material = TetElasticMaterial::try_from_resolved_orthotropic_state(
        &s.state,
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        fs_blake3::hash_domain("org.frankensim.pexfp49.axes.v1", b"x=L,y=R,z=T"),
    )
    .unwrap();
    let band = (100.0, 12_000.0);
    with_cx(|cx| {
        let q = quarter([3, 4, 2], H);
        let full = mirrored_full(&q);
        let classes: Vec<Solved> = CLASSES
            .iter()
            .map(|&c| solid_modes(&q, &material, Some(c), false, band, 1, &[], cx))
            .collect();
        let mut merged: Vec<f64> = classes
            .iter()
            .flat_map(|s| s.modes.iter().map(|m| m.hz))
            .collect();
        merged.sort_by(f64::total_cmp);
        let full_model = solid_modes(&full, &material, None, false, band, 1, &[], cx);
        // The six free-free rigid modes are excluded by inertia count: exactly
        // six eigenvalues lie below the window, and the classes split them.
        let class_rigid: usize = classes.iter().map(|s| s.below).sum();
        assert_eq!(
            full_model.below, 6,
            "full-model eigenvalues below the window"
        );
        assert_eq!(class_rigid, 6, "class eigenvalues below the window");
        let whole: Vec<f64> = full_model.modes.iter().map(|m| m.hz).collect();
        assert_eq!(merged.len(), whole.len(), "class counts must add up");
        let worst = merged
            .iter()
            .zip(&whole)
            .map(|(a, b)| (a / b - 1.0).abs())
            .fold(0.0, f64::max);
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W1c\",\"test\":\"g0_symmetry_class_merge_recovers_full_model_spectrum\",\"tier\":\"G0\",\"observed\":{{\"modes\":{},\"max_rel_diff\":{worst:e},\"rigid_below_window\":{},\"class_rigid_below_window\":{class_rigid}}},\"tolerance\":1e-10}}",
            whole.len(),
            full_model.below
        );
        assert!(worst <= 1e-10, "quarter classes vs full model: {worst:e}");
        assert!(whole.len() >= 10, "the window must hold a real spectrum");
    });
}

#[test]
fn g1_3d_reference_ladder_converges() {
    let _serial = SERIAL
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let study = spruce_study();
    let band: Vec<&Reference> = study.references.iter().filter(|r| in_band(r.hz)).collect();
    assert!(band.len() >= 15, "only {} modes in 2-8 kHz", band.len());
    let mut worst_bound = 0.0_f64;
    let mut lowest_order = f64::INFINITY;
    for r in band.iter().copied() {
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W1c\",\"test\":\"g1_3d_reference_ladder_converges\",\"tier\":\"G1\",\"axis\":\"h\",\"inputs\":{{\"pack\":\"spruce-sitka-fpl-gtr282\",\"pack_hash\":\"{}\",\"mc_pct\":{MC_PCT},\"h_m\":{H},\"class\":\"{}\"}},\"rungs\":{{\"layers\":{:?}}},\"observed\":{{\"f_hz\":{:?},\"order\":{},\"reference_hz\":{:.4},\"bound_hz\":{:e}}}}}",
            fs_blake3::hash_bytes(SPRUCE.as_bytes()).to_hex(),
            r.class.label(),
            RUNGS_3D.map(|r| 2 * r),
            r.rungs_hz,
            json_f64(r.order, 3),
            r.hz,
            r.bound_hz
        );
        // A non-finite order (non-monotone ladder) fails here: NaN >= 3 is false.
        assert!(
            r.order >= 3.0,
            "{} mode at {} Hz: 3-D frequency order {} < 3",
            r.class.label(),
            r.hz,
            r.order
        );
        worst_bound = worst_bound.max(r.bound_hz / r.hz);
        lowest_order = lowest_order.min(r.order);
    }
    println!(
        "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W1c\",\"test\":\"g1_3d_reference_ladder_converges\",\"observed\":{{\"modes\":{},\"min_order\":{lowest_order:.3},\"max_rel_bound\":{worst_bound:e},\"study_s\":{:.1}}},\"tolerance\":{{\"min_order\":3.0,\"max_rel_bound\":{}}}}}",
        band.len(),
        study.elapsed_s,
        TOL / 3.0
    );
    assert!(
        worst_bound <= TOL / 3.0,
        "3-D finest-rung bound {worst_bound:e} > TOL/3"
    );
}

#[test]
fn g2_shear_shell_matches_3d_reference_2_to_8_khz() {
    let _serial = SERIAL
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let study = spruce_study();
    let mut worst = 0.0_f64;
    let mut pairs = 0usize;
    for r in study.references.iter().filter(|r| in_band(r.hz)) {
        pairs += 1;
        let ladder: Vec<f64> = study.shell.iter().map(|m| partner(r, m).hz).collect();
        let f = *ladder.last().unwrap();
        let dev = (f - r.hz) / r.hz;
        worst = worst.max(dev.abs());
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W1c\",\"test\":\"g2_shear_shell_matches_3d_reference_2_to_8_khz\",\"tier\":\"G2\",\"inputs\":{{\"class\":\"{}\",\"kappa\":{KAPPA}}},\"rungs\":{{\"shell_cells\":{:?}}},\"observed\":{{\"f_3d\":{:.4},\"f_shell\":{ladder:?},\"dev_shell\":{dev:e}}},\"tolerance\":{TOL}}}",
            r.class.label(),
            RUNGS_SHELL.map(|k| [2 * BASE[0] * k, 2 * BASE[1] * k]),
            r.hz
        );
    }
    println!(
        "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W1c\",\"test\":\"g2_shear_shell_matches_3d_reference_2_to_8_khz\",\"observed\":{{\"pairs\":{},\"max_abs_dev\":{worst:e}}},\"tolerance\":{TOL},\"verdict\":\"{}\"}}",
        pairs,
        if worst <= TOL { "pass" } else { "fail" }
    );
    assert!(
        worst <= TOL,
        "shear shell deviates {worst:e} from the 3-D reference"
    );
}

#[test]
fn g2_dkt_deviation_table_and_refusal_band() {
    let _serial = SERIAL
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let study = spruce_study();
    let shell = study.shell.last().unwrap();
    let dkt = study.dkt.last().unwrap();
    // (f_3d, dev_shell, dev_dkt) in cents. A DKT mode that no Kirchhoff mode
    // matches (MAC < 0.9) is a model-form break: its deviation is recorded as
    // infinite, with the best MAC, and never paired by frequency.
    let mut table: Vec<(f64, f64, f64)> = Vec::new();
    let mut dkt_unpaired = 0usize;
    for r in &study.references {
        let fs = partner(r, shell).hz;
        let cs = cents(fs, r.hz);
        let (fd, cd, best_mac) = match try_partner(r, dkt) {
            Ok(m) => (m.hz, cents(m.hz, r.hz), mac(&r.samples, &m.samples)),
            Err((_, best)) => {
                dkt_unpaired += 1;
                (f64::NAN, f64::INFINITY, best)
            }
        };
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W1c\",\"test\":\"g2_dkt_deviation_table_and_refusal_band\",\"tier\":\"G2\",\"inputs\":{{\"class\":\"{}\"}},\"observed\":{{\"f_3d\":{:.4},\"f_shell\":{fs:.4},\"f_dkt\":{},\"dev_shell_cents\":{cs:.3},\"dev_dkt_cents\":{},\"dkt_mac\":{best_mac:.4}}}}}",
            r.class.label(),
            r.hz,
            json_f64(fd, 4),
            if cd.is_finite() {
                format!("{cd:.3}")
            } else {
                "\"unpaired\"".to_string()
            }
        );
        table.push((r.hz, cs, cd));
    }
    table.sort_by(|a, b| a.0.total_cmp(&b.0));
    // The refusal band: the lowest frequency above which every paired mode's
    // DKT deviation exceeds 10 c. Computed, never asserted.
    let mut refusal_hz = None;
    for (i, &(f, _, _)) in table.iter().enumerate() {
        if table[i..].iter().all(|&(_, _, d)| d.abs() > 10.0) {
            refusal_hz = Some(f);
            break;
        }
    }
    let near_4k = table
        .iter()
        .filter(|row| row.2.is_finite())
        .min_by(|a, b| (a.0 - 4000.0).abs().total_cmp(&(b.0 - 4000.0).abs()))
        .unwrap();
    let dev_4k = fs_math::det::pow(2.0, near_4k.2 / 1200.0) - 1.0;
    let estimate = if (0.02..=0.10).contains(&dev_4k) {
        "consistent"
    } else {
        "falsified"
    };
    println!(
        "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W1c\",\"test\":\"g2_dkt_deviation_table_and_refusal_band\",\"observed\":{{\"pairs\":{},\"dkt_unpaired\":{dkt_unpaired},\"dkt_refusal_hz\":{},\"mode_near_4khz\":{:.1},\"dkt_dev_near_4khz\":{dev_4k:e},\"plan_estimate\":[0.04,0.05],\"plan_estimate_status\":\"{estimate}\"}}}}",
        table.len(),
        refusal_hz.map_or("null".to_string(), |f| format!("{f:.1}")),
        near_4k.0
    );
    for &(f, s, d) in table.iter().filter(|row| in_band(row.0)) {
        assert!(
            d.abs() > s.abs(),
            "at {f} Hz DKT ({d:.3} c) does not deviate more than the shear shell ({s:.3} c)"
        );
    }
}

#[test]
fn g1_isotropic_thin_plate_three_way_agreement() {
    let _serial = SERIAL
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let (e, nu, rho) = (
        scalar(ALUMINIUM, "young_modulus"),
        scalar(ALUMINIUM, "poisson_ratio"),
        scalar(ALUMINIUM, "density"),
    );
    let g = e / (2.0 * (1.0 + nu)); // isotropy, not the table's separate G
    let h = 0.0005;
    // Solved well above the first ten modes, so that the stiffer coarse rungs
    // still hold every partner.
    let band = (20.0, 5000.0);
    let material = TetElasticMaterial::try_new(
        rho,
        e,
        nu,
        fs_blake3::hash_domain("org.frankensim.pexfp49.al2024.v1", ALUMINIUM.as_bytes()).0,
    )
    .unwrap();
    // At h/λ ≤ 0.01 the 3-D error is in-plane: one quadratic layer per half
    // already holds the bending field's z-dependence (linear, plus the
    // quadratic Poisson thickness change). Paired x1 rungs would put element
    // aspect near 50 and leave the coarse rung far from asymptotic (measured:
    // best MAC 0.69). So this sanity ladder refines in-plane only.
    let mut references = Vec::new();
    // Inertia below the window per rung, summed over the classes. The
    // flexural eighth keeps three of the six rigid modes: w and two tilts.
    let mut rigid = [0usize; 3];
    with_cx(|cx| {
        for class in CLASSES {
            let ladders: Vec<Vec<Mode>> = RUNGS_THIN
                .iter()
                .zip(rigid.iter_mut())
                .map(|(&r, rigid)| {
                    let mesh = eighth([BASE[0] * r, BASE[1] * r, 1], h);
                    let solved = solid_modes(
                        &mesh,
                        &material,
                        Some(class),
                        true,
                        band,
                        1,
                        &mid_plane_samples(BASE, mesh.counts, r),
                        cx,
                    );
                    *rigid += solved.below;
                    flexural(&solved.modes)
                })
                .collect();
            // Only the low modes the comparison uses are tracked (the first ten
            // lie far below 1.5 kHz); higher modes are pre-asymptotic on the
            // coarse rung (measured: best MAC 0.81 at 2.57 kHz) and a stiffer
            // rung can push a band-edge partner out of the window.
            for (index, m) in ladders[2].iter().enumerate() {
                if m.hz > 0.3 * band.1 {
                    continue;
                }
                let mut rungs = [0.0; 3];
                rungs[2] = m.hz;
                for rung in 0..2 {
                    rungs[rung] = track(index, &ladders[2], &ladders[rung]).unwrap_or_else(|e| {
                        let near: Vec<(f64, f64)> = ladders[rung]
                            .iter()
                            .map(|c| (c.hz, mac(&m.samples, &c.samples)))
                            .filter(|(_, q)| *q > 0.05)
                            .collect();
                        panic!(
                            "{} mode at {} Hz untracked on rung {rung}: {e:?}; candidates (hz, MAC) {near:?}; finest rung {:?}",
                            class.label(),
                            m.hz,
                            ladders[2].iter().map(|c| c.hz).collect::<Vec<_>>()
                        )
                    });
                }
                let (order, reference, bound) = ladder(rungs, RUNGS_THIN);
                references.push(Reference {
                    class,
                    rungs_hz: rungs,
                    order,
                    hz: reference,
                    bound_hz: bound,
                    samples: m.samples.clone(),
                });
            }
        }
    });
    references.sort_by(|a, b| a.hz.total_cmp(&b.hz));
    assert_eq!(rigid, [3; 3], "rigid modes below the window, per rung");
    assert!(
        references.len() >= 10,
        "only {} flexural modes",
        references.len()
    );
    let iso = OrthotropicElastic::new([e; 3], [nu; 3], [g; 3], 1e-3).unwrap();
    let shell_section =
        ShearShellSection::from_orthotropic_solid(e, e, nu, g, g, g, h, rho, 0.0, KAPPA).unwrap();
    // Shells on the same full-plate grid as the spruce finest rung.
    let shell = shell_modes(BASE, 8, &shell_section, band);
    let dkt = dkt_modes(
        BASE,
        8,
        &PlateSection::orthotropic(&iso, h, rho).unwrap(),
        band,
    );
    let mut worst = 0.0_f64;
    for r in references.iter().take(10) {
        let (fs, fd) = (partner(r, &shell).hz, partner(r, &dkt).hz);
        let spread = [r.hz, fs, fd];
        let hi = spread.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let lo = spread.iter().copied().fold(f64::INFINITY, f64::min);
        worst = worst.max(hi / lo - 1.0);
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W1c\",\"test\":\"g1_isotropic_thin_plate_three_way_agreement\",\"tier\":\"G1\",\"inputs\":{{\"pack\":\"aluminum-2024-t3-sheet-mil-hdbk-5j\",\"pack_hash\":\"{}\",\"h_m\":{h},\"class\":\"{}\"}},\"rungs\":{{\"in_plane\":{RUNGS_THIN:?},\"half_layers\":1,\"shell_cells\":{:?}}},\"observed\":{{\"f_3d\":{:.4},\"f_3d_rungs\":{:?},\"f_3d_order\":{},\"f_shell\":{fs:.4},\"f_dkt\":{fd:.4}}}}}",
            fs_blake3::hash_bytes(ALUMINIUM.as_bytes()).to_hex(),
            r.class.label(),
            [2 * BASE[0] * 8, 2 * BASE[1] * 8],
            r.hz,
            r.rungs_hz,
            json_f64(r.order, 2)
        );
    }
    println!(
        "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W1c\",\"test\":\"g1_isotropic_thin_plate_three_way_agreement\",\"observed\":{{\"max_spread\":{worst:e},\"rigid_below_window\":{rigid:?}}},\"tolerance\":0.002}}"
    );
    assert!(worst <= 0.002, "three-way spread {worst:e} > 0.2 %");
}
