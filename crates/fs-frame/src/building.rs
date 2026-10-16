//! Stage 3, FULL TIER: distributed-plasticity BUILDING frames (plan §15.2
//! step 3) on fs-solid's planar force-based fiber engine — the recorded
//! successor of the one-story hinge model in [`crate::history`].
//!
//! A regular moment frame of `stories × bays` reinforced-concrete members:
//! every column and beam is ONE force-based element with Gauss–Lobatto fiber
//! sections (Mander confined/unconfined concrete, Menegotto–Pinto steel), so
//! plasticity spreads along members instead of living in a calibrated hinge
//! spring. Geometry is corotational (P-Δ and large drift are exact chord
//! kinematics). The building is preloaded by its own floor weight, its modes
//! are extracted from the gravity-loaded tangent, Rayleigh damping is fitted
//! to the first and last retained lateral modes, and the base record drives
//! an implicit Newmark history with an energy ledger.
//!
//! [`building_fragility`] reruns the stage-4 anytime-valid study at this tier:
//! an e-stopped confidence sequence on the exceedance indicator of the PEAK
//! INTER-STORY DRIFT RATIO, plus a MULTI-FIDELITY MLMC on the peak drift
//! ratio whose level 0 is the cheap story model and whose level 1 is the
//! building-minus-story correction on the SAME ground motion (telescoping, so
//! unbiased for the building whatever the story model's bias).
//!
//! No-claim: the frame is planar with rigid joints, fixed bases, no shear
//! failure, bond-slip or joint-panel flexibility, and floor mass lumped at
//! beam–column nodes; the confidence sequence is valid for the ensemble's
//! synthetic Kanai–Tajimi motions only (recorded-motion suites are the next
//! successor), and the fragility is of the model, not of a real building.

use fs_eproc::GaussianMixtureCs;
use fs_scenario::ensemble::StochasticEnsemble;
use fs_solid::fiber::rc_section;
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel, SolidError};
use fs_uq::mlmc::{MlmcReport, mlmc_estimate};

use crate::assert_ground_motion_ensemble;
use crate::history::{StoryFrame, StoryParams, peak_drift};

/// Standard gravity (m/s²).
const G: f64 = 9.806_65;

/// A rectangular RC member section: depth, width, steel area per face (SI).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RcMember {
    /// Section depth (m), in the bending plane.
    pub depth: f64,
    /// Section width (m).
    pub width: f64,
    /// Longitudinal steel area per face (m²).
    pub steel_area: f64,
}

/// A regular RC moment-frame building.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BuildingSpec {
    /// Number of stories (≥ 1).
    pub stories: usize,
    /// Number of bays (≥ 1).
    pub bays: usize,
    /// Story height (m).
    pub story_height: f64,
    /// Bay width (m).
    pub bay_width: f64,
    /// Seismic mass per floor (kg), lumped equally at the floor's nodes.
    pub floor_mass: f64,
    /// Column section.
    pub column: RcMember,
    /// Beam section.
    pub beam: RcMember,
    /// Column fiber-area scale (the CVaR sizing variable of stage 5).
    pub scale: f64,
    /// Modal damping ratio for the Rayleigh fit.
    pub zeta: f64,
    /// Gauss–Lobatto sections per member (3..=7).
    pub sections_per_member: usize,
}

impl Default for BuildingSpec {
    fn default() -> Self {
        BuildingSpec {
            stories: 3,
            bays: 2,
            story_height: 3.2,
            bay_width: 6.0,
            floor_mass: 1.2e5,
            column: RcMember {
                depth: 0.5,
                width: 0.5,
                steel_area: 1.5e-3,
            },
            beam: RcMember {
                depth: 0.6,
                width: 0.35,
                steel_area: 1.2e-3,
            },
            scale: 1.0,
            zeta: 0.03,
            sections_per_member: 5,
        }
    }
}

/// A building model ready for time histories.
#[derive(Debug, Clone)]
pub struct BuildingFrame {
    frame: Frame2d,
    spec: BuildingSpec,
    /// Node ids per floor level (level 0 = base), left to right.
    floors: Vec<Vec<usize>>,
    periods: Vec<f64>,
    damping: Rayleigh,
}

/// One building time-history record.
#[derive(Debug, Clone, PartialEq)]
pub struct BuildingResponse {
    /// Peak inter-story drift ratio per story (bottom first).
    pub peak_interstory_drift: Vec<f64>,
    /// The largest of them (the fragility QoI).
    pub peak_drift_ratio: f64,
    /// Roof displacement relative to the base, per step.
    pub roof_drift: Vec<f64>,
    /// Peak base shear magnitude (N).
    pub peak_base_shear: f64,
    /// Final cumulative internal work (J) — strain energy plus hysteretic
    /// dissipation; at rest this is the energy the frame absorbed.
    pub internal_work: f64,
    /// Relative energy-balance error of the integration.
    pub energy_balance_error: f64,
}

fn scaled(member: RcMember, scale: f64) -> SectionModel {
    let mut s = rc_section(member.depth, member.width, 16, member.steel_area);
    for f in &mut s.fibers {
        f.area *= scale;
    }
    SectionModel::Fiber(s)
}

impl BuildingFrame {
    /// Build the frame, apply gravity, extract modes, fit damping.
    ///
    /// # Errors
    /// [`SolidError::InvalidInput`] for a degenerate spec; solver errors if
    /// the gravity preload or modal extraction fails.
    #[allow(clippy::needless_range_loop)] // floor/bay grid indexing
    pub fn new(spec: BuildingSpec) -> Result<BuildingFrame, SolidError> {
        if spec.stories == 0
            || spec.bays == 0
            || !(spec.story_height > 0.0 && spec.bay_width > 0.0 && spec.floor_mass > 0.0)
            || !(spec.scale > 0.0 && spec.zeta >= 0.0 && spec.zeta < 1.0)
        {
            return Err(SolidError::InvalidInput {
                what: "building spec needs ≥1 story and bay, positive geometry/mass/scale, \
                       and zeta in [0, 1)"
                    .into(),
            });
        }
        let mut frame = Frame2d::new(Geometry::Corotational);
        let mut floors = Vec::with_capacity(spec.stories + 1);
        for s in 0..=spec.stories {
            let row: Vec<usize> = (0..=spec.bays)
                .map(|b| frame.add_node(spec.bay_width * b as f64, spec.story_height * s as f64))
                .collect();
            floors.push(row);
        }
        for &n in &floors[0] {
            frame.fix(n, [true; 3]);
        }
        let np = spec.sections_per_member;
        let node_mass = spec.floor_mass / (spec.bays + 1) as f64;
        let column = || scaled(spec.column, spec.scale);
        let beam = || scaled(spec.beam, 1.0);
        for s in 1..=spec.stories {
            for b in 0..=spec.bays {
                frame.add_element(floors[s - 1][b], floors[s][b], np, 0.0, &column)?;
                frame.add_mass(floors[s][b], node_mass);
            }
            for b in 0..spec.bays {
                frame.add_element(floors[s][b], floors[s][b + 1], np, 0.0, &beam)?;
            }
        }
        // Gravity preload: the floor weight on its nodes.
        let mut gravity = vec![0.0; frame.ndof()];
        for row in &floors[1..] {
            for &n in row {
                gravity[3 * n + 1] = -node_mass * G;
            }
        }
        frame.static_load(&gravity, 4)?;
        let w = frame.natural_frequencies()?;
        let lateral: Vec<f64> = w.iter().copied().filter(|x| *x > 0.0).collect();
        let w1 = *lateral
            .first()
            .ok_or_else(|| SolidError::InternalInvariant {
                what: "building has no positive natural frequency".into(),
            })?;
        let w2 = lateral
            .get(spec.stories.min(lateral.len()) - 1)
            .copied()
            .unwrap_or(w1);
        let damping = if w2 > w1 {
            Rayleigh::from_modes(spec.zeta, w1, w2)
        } else {
            Rayleigh {
                a0: 2.0 * spec.zeta * w1,
                a1: 0.0,
            }
        };
        let periods = lateral
            .iter()
            .map(|x| 2.0 * std::f64::consts::PI / x)
            .collect();
        Ok(BuildingFrame {
            frame,
            spec,
            floors,
            periods,
            damping,
        })
    }

    /// Natural periods of the gravity-loaded frame (s), descending.
    #[must_use]
    pub fn periods(&self) -> &[f64] {
        &self.periods
    }

    /// The fitted Rayleigh damping.
    #[must_use]
    pub fn damping(&self) -> Rayleigh {
        self.damping
    }

    /// The building spec.
    #[must_use]
    pub fn spec(&self) -> &BuildingSpec {
        &self.spec
    }

    /// Run a horizontal base-acceleration record `ag` sampled every `dt`.
    /// Consumes the gravity-loaded state; clone the frame to rerun.
    ///
    /// # Errors
    /// Solver errors if a step fails to converge (reduce `dt`).
    pub fn run(&mut self, ag: &[f64], dt: f64) -> Result<BuildingResponse, SolidError> {
        let u_grav = self.frame.displacements().to_vec();
        let hist = self.frame.newmark(ag, dt, (1.0, 0.0), self.damping)?;
        let floor_ux = |u: &[f64], level: usize| -> f64 {
            let row = &self.floors[level];
            row.iter().map(|&n| u[3 * n] - u_grav[3 * n]).sum::<f64>() / row.len() as f64
        };
        let ns = self.spec.stories;
        let mut peak = vec![0.0f64; ns];
        let mut roof = Vec::with_capacity(hist.u.len());
        for u in &hist.u {
            let mut below = 0.0;
            for s in 1..=ns {
                let here = floor_ux(u, s);
                peak[s - 1] = peak[s - 1].max(((here - below) / self.spec.story_height).abs());
                below = here;
            }
            roof.push(below);
        }
        Ok(BuildingResponse {
            peak_drift_ratio: peak.iter().copied().fold(0.0, f64::max),
            peak_interstory_drift: peak,
            roof_drift: roof,
            peak_base_shear: hist.base_shear.iter().fold(0.0f64, |m, v| m.max(v.abs())),
            internal_work: hist.internal_work.last().copied().unwrap_or(0.0),
            energy_balance_error: hist.energy_balance_error(),
        })
    }
}

/// The building-tier fragility record.
#[derive(Debug, Clone)]
pub struct BuildingFragility {
    /// Members consumed before the e-stop.
    pub members_used: u32,
    /// Exceedance estimate (confidence-sequence center).
    pub p_hat: f64,
    /// Confidence-sequence radius at the stop.
    pub radius: f64,
    /// Stopped before exhausting the ensemble?
    pub stopped_early: bool,
    /// Exceedances observed.
    pub exceedances: u32,
    /// Peak inter-story drift ratio of every consumed member.
    pub peak_drifts: Vec<f64>,
    /// Multi-fidelity MLMC on the peak drift ratio (level 0: story model;
    /// level 1: building − story on the same motion).
    pub mlmc: MlmcReport,
}

/// Anytime-valid fragility of the BUILDING: `P(peak inter-story drift ratio
/// > drift_limit)` over the ensemble, e-stopped once the confidence-sequence
/// radius is `≤ margin` (after at least 8 members), plus the multi-fidelity
/// MLMC report. `story` parameterizes the level-0 story model.
///
/// # Errors
/// Solver errors from any building time history.
///
/// # Panics
/// On ensemble realization errors or non-ground-motion ensembles (spec
/// defects, as in [`crate::fragility::e_stopped_fragility`]).
pub fn building_fragility(
    ensemble: &StochasticEnsemble,
    spec: BuildingSpec,
    story: StoryParams,
    drift_limit: f64,
    alpha: f64,
    margin: f64,
) -> Result<BuildingFragility, SolidError> {
    assert_ground_motion_ensemble(ensemble);
    let dt = ensemble.dt.value;
    let base = BuildingFrame::new(spec)?;
    let mut cs = GaussianMixtureCs::new(0.5, 8.0, alpha);
    let (mut used, mut exceedances, mut stopped_early) = (0u32, 0u32, false);
    let mut peak_drifts = Vec::new();
    for member in 0..ensemble.members {
        let real = ensemble.realize(member).expect("ensemble realizes");
        let pd = base.clone().run(&real.values, dt)?.peak_drift_ratio;
        peak_drifts.push(pd);
        let x = if pd > drift_limit { 1.0 } else { 0.0 };
        if x > 0.5 {
            exceedances += 1;
        }
        cs.observe(x);
        used = member + 1;
        if let Some((_, radius)) = cs.interval()
            && used >= 8
            && radius <= margin
        {
            stopped_early = used < ensemble.members;
            break;
        }
    }
    let (p_hat, radius) = cs.interval().expect("observed at least one member");
    // Multi-fidelity MLMC: the sampler must be infallible, so a building
    // failure poisons the sample with NaN and is reported below.
    let mut failure: Option<SolidError> = None;
    let mut sampler = |level: usize, g: u64| -> f64 {
        let member = u32::try_from(g % u64::from(ensemble.members)).expect("small");
        let real = ensemble.realize(member).expect("ensemble realizes");
        let coarse = {
            let mut frame = StoryFrame::new(story);
            peak_drift(&frame.run(&real.values, dt), story.h)
        };
        if level == 0 {
            return coarse;
        }
        match base.clone().run(&real.values, dt) {
            Ok(r) => r.peak_drift_ratio - coarse,
            Err(e) => {
                failure.get_or_insert(e);
                f64::NAN
            }
        }
    };
    let mlmc = mlmc_estimate(&mut sampler, &[1.0, 40.0], 4, 1e-6);
    if let Some(e) = failure {
        return Err(e);
    }
    Ok(BuildingFragility {
        members_used: used,
        p_hat,
        radius,
        stopped_early,
        exceedances,
        peak_drifts,
        mlmc,
    })
}
