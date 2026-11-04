//! Held affine transport and positive viscous evolution on one adaptive clock.

use super::{
    AdaptiveError, AdaptiveRun, AmbientAdaptiveBudget, CompactNeighborBudget,
    ViscousAdaptiveOptions, viscous_run,
};
use crate::conservative::ambient::{AffineFlow2, AmbientMomentChange, HalfFlow};
use crate::diffusion::PseDiffusion;
use crate::diffusion::moving::MovingDiffusion;
use crate::{AdmittedBudget, Cx, VortexParticle, VpmError, WorkMeter};
use fs_math::dd::Dd;

/// Actual accepted PSE vortex-moment changes and native representation defects.
///
/// These are discrete particle moments about the world origin, not a closed
/// finite-volume fluid/body momentum or viscous-work account. On an arbitrary
/// finite cloud the original positive PSE law conserves circulation but does
/// not enforce continuum impulse identities. Its commanded pair transfers
/// therefore generally change first/quadratic moments. Those changes are
/// retained, separately from floating-point update loss and the held source.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct DiffusionMomentChange {
    /// Actual change in `sum Gamma*x` and `sum Gamma*y` [m^3/s].
    pub first_moment_m3_s: [f64; 2],
    /// Actual change in `sum Gamma*|x|^2` [m^4/s].
    pub quadratic_moment_m4_s: f64,
    /// Actual circulation change; each commanded pair has zero sum [m^2/s].
    pub circulation_residual_m2_s: f64,
    /// Actual minus commanded pair-transfer first moments [m^3/s].
    pub first_moment_residual_m3_s: [f64; 2],
    /// Actual minus commanded pair-transfer quadratic moment [m^4/s].
    pub quadratic_moment_residual_m4_s: f64,
}

/// Accepted native viscous trajectory with disjoint source and PSE accounts.
#[derive(Debug, Clone, PartialEq)]
pub struct AmbientViscousAdaptiveRun {
    /// The sole published particle inventory and original adaptive work receipt.
    pub run: AdaptiveRun,
    /// Actual vortex-moment increments during accepted affine maps only.
    pub ambient_moment_change: AmbientMomentChange,
    /// Actual accepted PSE increments and their pair-update representation loss.
    pub diffusion_moment_change: DiffusionMomentChange,
    /// Computational source visits, including coarse, rejected and witness work.
    pub source_particle_updates: u64,
    /// Completed compact geometry rebuilds, including rejected trials; zero
    /// for the dense path or inactive diffusion.
    pub graph_rebuilds: usize,
    /// Actual compact distance tests across all geometry rebuilds; zero for
    /// the dense path or inactive diffusion.
    pub candidate_visits: u64,
}

/// Adapt a moving viscous cloud in one held incompressible affine background.
///
/// Each native second-order step is `F(h/2) D(h/2) A(h) D(h/2) F(h/2)`.
/// F is the original exact held-source map, A the original exact-pair advection,
/// and D the original positive PSE exchange at its CURRENT substage geometry.
/// Particle areas stay attached because F is incompressible. A coarse step and
/// two half steps must satisfy both original vorticity and position tolerances;
/// only the fine state is accepted, without extrapolation or moment repair.
///
/// Positive viscosity retains the independent reverse-order stiffness witness
/// `F(h/2) A(h/2) D(h) A(h/2) F(h/2)`. This exposes initial induced transport
/// that rapid diffusion can erase in both ordinary candidates. The witness is
/// an additional local estimator, not a proof for every stiff source transient.
/// Rejected/coarse/witness evolution adds computational work but no physical
/// source or diffusion moment to the accepted-path accounts.
///
/// Adjacent fine source maps coalesce. With P=N(N-1), a positive-viscosity trial
/// with N>=2 uses 12P pair visits, 7P diffusion-kernel evaluations and 7N source
/// updates; its total work is 12P+11N+5. Zero viscosity or a singleton uses 3P
/// pairs, 5N source updates and 3P+7N+4 work. Zero fields add no source work.
/// Three particle banks are reserved once. Moment increments are observed
/// inside existing fine-candidate pair visits; only accepted increments are
/// retained, with constant stack storage and no snapshot cloud/pass.
///
/// All whole-window attempt, pair, source, work and memory limits include
/// rejected trials. Both source maps are prepared before candidate copy work.
/// Source-active windows share the original inviscid compensated clock, so
/// rounded accumulated time cannot invent an unresolvable final source map.
/// Zero-source native trajectories and execution counts match
/// [`super::simulate_viscous_with_cx`]; this entry additionally admits the finite
/// moment receipt. Empty/zero-duration cases retain the original identity path;
/// a nonzero held source transports a singleton even without pair interactions.
///
/// This is temporal error control of the stated finite-cloud model. Kernel
/// overlap, padding and quadrature still require spatial refinement. PSE moment
/// drift is neither labeled as an external source nor projected away. The held
/// background adds no finite carrier recoil, boundary work or viscous heating.
/// Any numerical, moment-representation, budget or cancellation refusal returns
/// no partial cloud and leaves all borrowed inputs unchanged.
#[allow(clippy::too_many_arguments)]
pub fn simulate_viscous_with_affine_flow_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    areas_m2: &[f64],
    duration_s: f64,
    core: f64,
    law: PseDiffusion,
    flow: AffineFlow2,
    options: ViscousAdaptiveOptions,
    budget: AmbientAdaptiveBudget,
) -> Result<AmbientViscousAdaptiveRun, AdaptiveError> {
    viscous_run::<true>(
        cx,
        particles,
        areas_m2,
        duration_s,
        core,
        law,
        (!flow.is_zero()).then_some(flow),
        options,
        budget,
        None,
    )
}

/// The same source-coupled viscous trajectory with bounded compact neighbors.
///
/// This selects the original `MovingDiffusion` backend in the same adaptive
/// owner. Only compact laws are admitted; Gaussian tails are never truncated.
/// Every active source map invalidates the geometry, as do native advection
/// and candidate resets. Fractions retain the original current-geometry
/// evaluation and shared forward/reverse exchange cache. The accepted state
/// and physical source/PSE accounts match the dense compact law bit for bit.
///
/// All graph rebuild/search/sort/fraction work, source visits, witness work
/// and rejected trials use the same whole-window caps and execution account.
/// The complete neighbor workspace is reserved once alongside three particle
/// banks. Work with data-dependent edge counts is admitted before its bounded
/// operations. Direct advection remains O(N^2); no measured speedup is claimed.
#[allow(clippy::too_many_arguments)]
pub fn simulate_viscous_compact_with_affine_flow_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    areas_m2: &[f64],
    duration_s: f64,
    core: f64,
    law: PseDiffusion,
    flow: AffineFlow2,
    options: ViscousAdaptiveOptions,
    budget: AmbientAdaptiveBudget,
    neighbors: CompactNeighborBudget,
) -> Result<AmbientViscousAdaptiveRun, AdaptiveError> {
    viscous_run::<true>(
        cx,
        particles,
        areas_m2,
        duration_s,
        core,
        law,
        (!flow.is_zero()).then_some(flow),
        options,
        budget,
        Some(neighbors),
    )
}

pub(super) struct SourceStep {
    pub half: HalfFlow,
    pub quarter: HalfFlow,
}
impl SourceStep {
    pub fn new(flow: AffineFlow2, h: f64) -> Result<Self, VpmError> {
        Ok(Self {
            half: HalfFlow::new(flow, h)?,
            // Keep the physical interval and its fraction separate, including
            // when h/2 itself would round to zero before a finite source rate.
            quarter: HalfFlow::weighted(flow, h, 0.5)?,
        })
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn apply_source(
    backend: &mut MovingDiffusion,
    source: Option<&HalfFlow>,
    state: &mut [VortexParticle],
    moments: &mut [Dd; 3],
    updates: &mut u64,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Result<(), VpmError> {
    if let Some(source) = source {
        backend.admit_work(state.len() as u64, meter, accountant)?;
        source.apply(state, moments, meter, accountant, cx)?;
        // Source motion can create or remove supported pairs even between
        // adjacent fine diffusion stages. Never reuse their old geometry.
        backend.invalidate();
        // The complete trial's checked plan admits this count before any work.
        *updates += state.len() as u64;
    }
    Ok(())
}

pub(super) fn add_source(total: &mut [Dd; 3], increment: [Dd; 3]) -> Result<(), VpmError> {
    for i in 0..3 {
        total[i] = total[i] + increment[i];
        if !finite(total[i]) {
            return Err(range_error(0));
        }
    }
    Ok(())
}

#[derive(Default)]
pub(super) struct DiffusionMoments {
    commanded: [Dd; 3],
    residual: [Dd; 4],
}
impl DiffusionMoments {
    pub fn observe(
        &mut self,
        particle: usize,
        before: [VortexParticle; 2],
        after: [VortexParticle; 2],
        transfer: f64,
    ) -> Result<(), VpmError> {
        if transfer == 0.0 {
            return Ok(());
        }
        let t = dd(transfer);
        let first = before[0].pos.map(dd);
        let second = before[1].pos.map(dd);
        let difference = [first[0] - second[0], first[1] - second[1]];
        let sum = [first[0] + second[0], first[1] + second[1]];
        let quadratic = |strength: Dd| {
            (strength * difference[0]) * sum[0] + (strength * difference[1]) * sum[1]
        };
        let actual = [0, 1].map(|i| dd(after[i].circulation) - dd(before[i].circulation));
        let defect = [actual[0] - t, actual[1] + t];
        let circulation_defect = defect[0] + defect[1];
        // Polarized relative-position products avoid subtracting large whole
        // moments or large squared positions when a pair exchanges little.
        let commanded = [t * difference[0], t * difference[1], quadratic(t)];
        let origin_quadratic = if circulation_defect == Dd::ZERO {
            Dd::ZERO
        } else {
            (circulation_defect * first[0]) * first[0] + (circulation_defect * first[1]) * first[1]
        };
        let residual = [
            circulation_defect,
            circulation_defect * first[0] - defect[1] * difference[0],
            circulation_defect * first[1] - defect[1] * difference[1],
            origin_quadratic - quadratic(defect[1]),
        ];
        for i in 0..3 {
            self.commanded[i] = self.commanded[i] + commanded[i];
        }
        for i in 0..4 {
            self.residual[i] = self.residual[i] + residual[i];
        }
        self.admit(particle)
    }

    pub fn add(&mut self, other: &Self) -> Result<(), VpmError> {
        for i in 0..3 {
            self.commanded[i] = self.commanded[i] + other.commanded[i];
        }
        for i in 0..4 {
            self.residual[i] = self.residual[i] + other.residual[i];
        }
        self.admit(0)
    }

    fn admit(&self, particle: usize) -> Result<(), VpmError> {
        if self
            .commanded
            .iter()
            .chain(&self.residual)
            .any(|&v| !finite(v))
            || (0..3).any(|i| !finite(self.commanded[i] + self.residual[i + 1]))
        {
            return Err(range_error(particle));
        }
        Ok(())
    }

    pub fn finish(self) -> DiffusionMomentChange {
        DiffusionMomentChange {
            first_moment_m3_s: [
                (self.commanded[0] + self.residual[1]).to_f64(),
                (self.commanded[1] + self.residual[2]).to_f64(),
            ],
            quadratic_moment_m4_s: (self.commanded[2] + self.residual[3]).to_f64(),
            circulation_residual_m2_s: self.residual[0].to_f64(),
            first_moment_residual_m3_s: [self.residual[1].to_f64(), self.residual[2].to_f64()],
            quadratic_moment_residual_m4_s: self.residual[3].to_f64(),
        }
    }
}

fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}
fn finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite() && value.to_f64().is_finite()
}
fn range_error(particle: usize) -> VpmError {
    VpmError::NonFiniteResult {
        phase: "fs-vpm.viscous-adaptive-moments",
        particle,
    }
}
