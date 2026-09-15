//! The compliance descent loop — THE marquee coupling: physics by
//! fs-cutfem's canonical elasticity operator DIRECTLY on the evolving [`GridSdf`]
//! (zero meshing, ever), shape velocity `v_n = w − ℓ` from the energy
//! density (grow where strain energy exceeds the volume multiplier,
//! shrink where it doesn't), Sobolev-smoothed through fs-adjoint's
//! Riesz step, extended off the interface, advected by WENO5 normal
//! flow, redistanced with drift audits, and — on schedule — hole
//! nucleation by topological derivative. Volume rides an augmented
//! Lagrangian multiplier; every iteration ledgers compliance, volume,
//! multiplier, audits, events, and an FNV snapshot hash of φ.

use crate::fim::{RedistanceAudit, redistance};
use crate::gridsdf::GridSdf;
use crate::topder::{NucleationEvent, nucleate, topological_derivative};
use crate::veloext::extend_velocity;
use crate::weno::{Velocity, advect, build_band};
use fs_cutfem::quad::cut_cell_rules;
use fs_cutfem::{
    BoundaryTraction, CutElasticity, CutElasticitySolution, CutFemError, CutSdf, DesignBoxEdge,
    EdgeBand, MAX_PLANE_STRAIN_STIFFNESS_RATIO, Quadtree,
};
use fs_material::IsotropicElastic;
use std::fmt::Write as _;

const MATERIAL_STRAIN_LIMIT: f64 = 1.0;
const SOLVER_TOL: f64 = 1e-12;
const SOLVER_MAX_ITERS: usize = 60_000;

/// Optimizer controls.
#[derive(Debug, Clone, Copy)]
pub struct OptimizeSettings {
    /// Grid level (cells per side = 2^level; GridSdf n must match).
    pub level: u32,
    /// Target volume fraction of the design box.
    pub volfrac: f64,
    /// Descent iterations.
    pub iterations: usize,
    /// Narrow band half-width in cells.
    pub band_cells: f64,
    /// Interface travel per iteration, in cells.
    pub move_cells: f64,
    /// Initial volume multiplier ℓ.
    pub ell0: f64,
    /// Augmented-Lagrangian multiplier gain.
    pub mu_al: f64,
    /// Sobolev smoothing α (≈ h²·scale).
    pub sobolev_alpha: f64,
    /// Nucleation period (0 disables).
    pub nucleation_period: usize,
    /// Nucleation hole radius (in cells).
    pub hole_radius_cells: f64,
    /// Young's modulus; must be finite and positive.
    pub youngs: f64,
    /// Poisson ratio in the canonical certified plane-strain regime:
    /// `(lambda + 2*mu) / mu <= 4` (equivalently `nu <= 1/3`).
    pub poisson: f64,
}

impl Default for OptimizeSettings {
    fn default() -> Self {
        OptimizeSettings {
            level: 5,
            volfrac: 0.5,
            iterations: 25,
            band_cells: 6.0,
            move_cells: 0.5,
            ell0: 0.0,
            mu_al: 4.0,
            sobolev_alpha: 2.0,
            nucleation_period: 8,
            hole_radius_cells: 2.5,
            youngs: 1.0,
            poisson: 0.3,
        }
    }
}

/// The ledgered trajectory.
#[derive(Debug, Clone, Default)]
pub struct OptimizeReport {
    /// Compliance of the successfully evaluated post-evolution geometry per iteration.
    pub compliance: Vec<f64>,
    /// Material volume per iteration.
    pub volume: Vec<f64>,
    /// Multiplier per iteration.
    pub ell: Vec<f64>,
    /// Redistancing audits.
    pub audits: Vec<RedistanceAudit>,
    /// Nucleation events.
    pub events: Vec<NucleationEvent>,
    /// FNV-64 hashes of the φ bits per iteration (evolution snapshots).
    pub snapshots: Vec<u64>,
    /// Nodal assignments made to retain the non-design load pad each
    /// iteration. This is intervention evidence, not a convergence metric.
    pub load_pad_nodes: Vec<usize>,
    /// Ledger rows.
    pub rows: Vec<String>,
}

fn fnv(phi: &GridSdf) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for v in phi.nodes() {
        for b in v.to_bits().to_le_bytes() {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    h
}

/// Material area of `{φ < 0}` by certified cut quadrature.
#[must_use]
pub fn material_volume(grid: &Quadtree, phi: &GridSdf) -> f64 {
    let mut vol = 0.0;
    for c in grid.leaves() {
        let (lo, hi) = grid.rect(c);
        let iv = phi.enclose(lo, hi);
        if iv.hi() < 0.0 {
            vol += (hi[0] - lo[0]) * (hi[1] - lo[1]);
        } else if iv.lo() <= 0.0 {
            vol += cut_cell_rules(phi, lo, hi, 2)
                .bulk
                .iter()
                .map(|&(_, w)| w)
                .sum::<f64>();
        }
    }
    vol
}

/// The cantilever fixture: clamp on the left box edge, downward
/// traction band on the right edge around mid-height.
#[derive(Debug, Clone, Copy)]
pub struct Cantilever {
    /// Finite strictly positive traction magnitude.
    pub load: f64,
    /// Finite load-band half-width in `[0, 0.5]`.
    pub band: f64,
}

fn invalid_input(what: impl Into<String>) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: what.into() }
}

fn cantilever_support(fixture: Cantilever) -> Result<EdgeBand, CutFemError> {
    if !(fixture.load.is_finite() && fixture.load > 0.0) {
        return Err(invalid_input(format!(
            "cantilever load magnitude {} must be finite and strictly positive",
            fixture.load
        )));
    }
    EdgeBand::new(DesignBoxEdge::Right, 0.5 - fixture.band, 0.5 + fixture.band).map_err(|error| {
        invalid_input(format!(
            "cantilever load-band half-width {} must be finite and lie in [0, 0.5]: {error}",
            fixture.band
        ))
    })
}

/// Retain a two-cell-deep non-design material pad around the checked load.
///
/// The pad extends one cell past each support endpoint and one cell beyond the
/// design box, so the supported right-edge segment is strictly inside material
/// rather than coincident with the pad boundary. Only the final three lattice
/// columns can change. The caller still receives a fail-closed refusal when its
/// initial geometry cuts the support; this policy repairs only geometry changes
/// made by the optimizer after a successful canonical solve.
fn retain_cantilever_load_pad(phi: &mut GridSdf, support: EdgeBand) -> usize {
    let n = phi.n();
    let h = phi.h();
    let x_min = 1.0 - 2.0 * h;
    let x_max = 1.0 + h;
    let y_min = support.start() - h;
    let y_max = support.end() + h;
    let mut changed = 0usize;
    for j in 0..=n {
        for i in n.saturating_sub(2)..=n {
            let [x, y] = phi.pos(i, j);
            let pad = (x_min - x).max(x - x_max).max(y_min - y).max(y - y_max);
            if pad <= 0.0 && pad < phi.node(i, j) {
                *phi.node_mut(i, j) = pad;
                changed += 1;
            }
        }
    }
    changed
}

fn validated_plane_strain_material(
    settings: OptimizeSettings,
) -> Result<(IsotropicElastic, f64, f64), CutFemError> {
    if !(settings.youngs.is_finite() && settings.youngs > 0.0) {
        return Err(invalid_input(format!(
            "optimizer Young's modulus {} must be finite and positive",
            settings.youngs
        )));
    }
    if !(settings.poisson.is_finite() && settings.poisson > -1.0 && settings.poisson < 0.5) {
        return Err(invalid_input(format!(
            "optimizer Poisson ratio {} must lie in (-1, 0.5)",
            settings.poisson
        )));
    }
    let material = IsotropicElastic::new(settings.youngs, settings.poisson, MATERIAL_STRAIN_LIMIT)
        .map_err(|error| invalid_input(format!("optimizer material card was refused: {error}")))?;
    let (lambda, mu) = material.lame();
    let bulk_2d = lambda + mu;
    let stiffness_ratio = (lambda + 2.0 * mu) / mu;
    if !(lambda.is_finite()
        && mu.is_finite()
        && mu > 0.0
        && bulk_2d.is_finite()
        && bulk_2d > 0.0
        && stiffness_ratio.is_finite())
    {
        return Err(invalid_input(
            "optimizer material does not define a finite coercive plane-strain law",
        ));
    }
    if stiffness_ratio > MAX_PLANE_STRAIN_STIFFNESS_RATIO {
        return Err(invalid_input(format!(
            "optimizer plane-strain stiffness ratio (lambda + 2*mu)/mu = {stiffness_ratio} exceeds the certified limit {MAX_PLANE_STRAIN_STIFFNESS_RATIO}"
        )));
    }
    Ok((material, lambda, mu))
}

/// Uniform Q1 mass/stiffness on the full node lattice (Sobolev step).
fn mass_stiffness(n: usize) -> (fs_sparse::Csr, fs_sparse::Csr) {
    #[allow(clippy::cast_precision_loss)]
    let h = 1.0 / n as f64;
    let stride = n + 1;
    let nn = stride * stride;
    let mut mc = fs_sparse::Coo::new(nn, nn);
    let mut kc = fs_sparse::Coo::new(nn, nn);
    // Q1 element matrices on a square of side h (standard closed
    // forms): lumped mass h²/4 per corner; stiffness pattern of the
    // Laplacian.
    let ke = [
        [2.0 / 3.0, -1.0 / 6.0, -1.0 / 3.0, -1.0 / 6.0],
        [-1.0 / 6.0, 2.0 / 3.0, -1.0 / 6.0, -1.0 / 3.0],
        [-1.0 / 3.0, -1.0 / 6.0, 2.0 / 3.0, -1.0 / 6.0],
        [-1.0 / 6.0, -1.0 / 3.0, -1.0 / 6.0, 2.0 / 3.0],
    ];
    for cj in 0..n {
        for ci in 0..n {
            let ids = [
                ci + cj * stride,
                ci + 1 + cj * stride,
                ci + 1 + (cj + 1) * stride,
                ci + (cj + 1) * stride,
            ];
            for (a, &ia) in ids.iter().enumerate() {
                mc.push(ia, ia, 0.25 * h * h);
                for (b, &ib) in ids.iter().enumerate() {
                    kc.push(ia, ib, ke[a][b]);
                }
            }
        }
    }
    (mc.assemble(), kc.assemble())
}

mod engine;
pub use engine::optimize_compliance;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn retained_load_pad_is_strict_material_and_replays_without_global_fill() {
        let support =
            EdgeBand::new(DesignBoxEdge::Right, 0.375, 0.625).expect("valid right-edge support");
        let initial = GridSdf::from_fn(8, &|_, y| y - 0.5);
        let mut first = initial.clone();
        let mut replay = initial.clone();

        let changed = retain_cantilever_load_pad(&mut first, support);
        let replay_changed = retain_cantilever_load_pad(&mut replay, support);
        assert!(changed > 0, "fixture must exercise load-pad retention");
        assert_eq!(changed, replay_changed);
        assert_eq!(first.nodes(), replay.nodes());
        assert!(first.value_at([1.0, support.start()]) < 0.0);
        assert!(first.value_at([1.0, support.end()]) < 0.0);
        assert!(
            first
                .enclose([1.0, support.start()], [1.0, support.end()])
                .hi()
                < 0.0,
            "the complete loaded edge band must be strictly inside material"
        );
        assert_eq!(
            first.node(0, 8).to_bits(),
            initial.node(0, 8).to_bits(),
            "load-pad retention must not fill unrelated lattice columns"
        );
    }

    fn beam(level: u32) -> GridSdf {
        GridSdf::from_fn(1usize << level, &|_, y| (y - 0.5).abs() - 0.35)
    }

    fn oracle_compliance(phi: &GridSdf, level: u32) -> f64 {
        let grid = Quadtree::uniform(level);
        let material = IsotropicElastic::new(1.0, 0.3, 1.0).expect("fixture material");
        let clamp = |x: f64, _: f64| x < 1e-9;
        let traction = |_: f64, _: f64| [0.0, -1.0];
        CutElasticity {
            grid: &grid, sdf: phi, material: &material,
            nitsche_beta: 20.0, ghost_gamma: 0.5,
            stabilization_scaling: fs_cutfem::CutStabilizationScaling::LongitudinalModulus,
            quad_depth: 2, clamp: Some(&clamp), boundary_traction: None,
            traction_free_interface: true, solver_tol: 1e-12, solver_max_iters: 60_000,
        }.solve_with_boundary_traction(
            &|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0],
            BoundaryTraction::EdgeBand {
                support: EdgeBand::new(DesignBoxEdge::Right, 0.375, 0.625).expect("fixture support"),
                value: &traction,
            },
        ).expect("independent final-geometry solve").compliance()
    }

    #[test]
    fn reported_compliance_and_snapshot_belong_to_the_returned_geometry() {
        let settings = OptimizeSettings {
            level: 4, iterations: 2, move_cells: 0.1, nucleation_period: 0,
            ..OptimizeSettings::default()
        };
        let fixture = Cantilever { load: 1.0, band: 0.125 };
        let mut phi = beam(settings.level);
        let original_snapshot = fnv(&phi);
        let report = optimize_compliance(&mut phi, fixture, settings).expect("real two-step run");
        assert_eq!(report.compliance.len(), 2);
        assert_ne!(fnv(&phi), original_snapshot, "the fixture must actually evolve");
        assert_eq!(report.snapshots.last(), Some(&fnv(&phi)));
        assert_eq!(report.compliance[1].to_bits(), oracle_compliance(&phi, settings.level).to_bits());
        assert_eq!(report.volume[1].to_bits(), material_volume(&Quadtree::uniform(settings.level), &phi).to_bits());
    }

    #[test]
    fn nucleation_volume_and_compliance_are_measured_after_the_holes() {
        let settings = OptimizeSettings {
            level: 5, iterations: 2, move_cells: 0.0, ell0: 1e6,
            nucleation_period: 1, hole_radius_cells: 1.25,
            ..OptimizeSettings::default()
        };
        let mut phi = beam(settings.level);
        let report = optimize_compliance(&mut phi, Cantilever { load: 1.0, band: 0.125 }, settings)
            .expect("real nucleation run");
        assert!(!report.events.is_empty(), "the regression must punch actual holes");
        assert!(report.volume[1] < report.volume[0]);
        assert_eq!(report.volume[1].to_bits(), material_volume(&Quadtree::uniform(settings.level), &phi).to_bits());
        assert_eq!(report.compliance[1].to_bits(), oracle_compliance(&phi, settings.level).to_bits());
    }

    #[test]
    fn bad_evolution_controls_refuse_without_mutating_the_design() {
        let original = beam(3);
        let settings = OptimizeSettings { level: 3, iterations: 0, ..OptimizeSettings::default() };
        for bad in [
            OptimizeSettings { level: u32::MAX, ..settings },
            OptimizeSettings { volfrac: 0.0, ..settings },
            OptimizeSettings { volfrac: f64::NAN, ..settings },
            OptimizeSettings { move_cells: f64::INFINITY, ..settings },
            OptimizeSettings { band_cells: 0.0, ..settings },
            OptimizeSettings { mu_al: 0.0, ..settings },
            OptimizeSettings { sobolev_alpha: -1.0, ..settings },
        ] {
            let mut phi = original.clone();
            assert!(optimize_compliance(&mut phi, Cantilever { load: 1.0, band: 0.125 }, bad).is_err());
            assert_eq!(phi.nodes(), original.nodes());
        }
    }

}
