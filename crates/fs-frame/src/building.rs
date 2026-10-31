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
//! an implicit Newmark history with the floor weight still applied and an
//! energy ledger that includes its work.
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
use fs_solid::fiber::material::RcSection;
use fs_solid::fiber::rc_section;
use fs_solid::frame2d::{FrameMassFormulation, FrameModalAnalysis, FrameModeStatus};
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
    /// Uniform gravity in the fixed building axes `[horizontal, vertical]`
    /// (m/s²). Positive vertical points up. Zero is legal; the default is
    /// terrestrial downward gravity. This is a held physical load, not base
    /// excitation or a change of reference frame.
    pub gravity_m_s2: [f64; 2],
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
            gravity_m_s2: [0.0, -G],
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
    modes: FrameModalAnalysis,
    damping: Rayleigh,
    gravity: Vec<f64>,
    /// Accepted column and beam sections, including their immutable source bindings.
    material_sections: Option<[RcSection; 2]>,
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

fn validate_spec(spec: BuildingSpec) -> Result<(), SolidError> {
    let positive = [
        spec.story_height,
        spec.bay_width,
        spec.floor_mass,
        spec.scale,
        spec.column.depth,
        spec.column.width,
        spec.beam.depth,
        spec.beam.width,
    ];
    if spec.stories == 0
        || spec.bays == 0
        || positive.iter().any(|v| !v.is_finite() || *v <= 0.0)
        || spec.gravity_m_s2.iter().any(|v| !v.is_finite())
        || [spec.column.steel_area, spec.beam.steel_area]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        || !spec.zeta.is_finite()
        || !(0.0..1.0).contains(&spec.zeta)
    {
        return Err(SolidError::InvalidInput {
            what: "building spec needs finite positive geometry/mass/scale, nonnegative \
                   reinforcement, finite gravity, at least one story and bay, and zeta in [0, 1)"
                .into(),
        });
    }
    Ok(())
}

impl BuildingFrame {
    /// Build the frame, apply gravity, extract modes, fit damping.
    ///
    /// # Errors
    /// [`SolidError::InvalidInput`] for a degenerate spec; solver errors if
    /// the gravity preload or modal extraction fails.
    #[allow(clippy::needless_range_loop)] // floor/bay grid indexing
    pub fn new(spec: BuildingSpec) -> Result<BuildingFrame, SolidError> {
        Self::with_member_mass(spec, 0.0, 0.0, FrameMassFormulation::Lumped)
    }

    /// Include explicitly supplied column/beam line mass in gravity, modes and
    /// the earthquake history. Floor mass must exclude these member masses to
    /// avoid double counting. Consistent inertia remains in the original frame;
    /// this is not geometrically exact moving-mass dynamics.
    ///
    /// # Errors
    /// Invalid line masses/specification, gravity equilibrium or modal solve.
    #[allow(clippy::needless_range_loop)]
    pub fn with_member_mass(
        spec: BuildingSpec,
        column_mass_per_length: f64,
        beam_mass_per_length: f64,
        mass_formulation: FrameMassFormulation,
    ) -> Result<BuildingFrame, SolidError> {
        if [column_mass_per_length, beam_mass_per_length]
            .iter()
            .any(|m| !m.is_finite() || *m < 0.0)
        {
            return Err(SolidError::InvalidInput {
                what: "member line mass must be finite and nonnegative".into(),
            });
        }
        validate_spec(spec)?;
        Self::assemble(
            spec,
            scaled(spec.column, spec.scale),
            scaled(spec.beam, 1.0),
            [column_mass_per_length, beam_mass_per_length],
            mass_formulation,
            None,
        )
    }

    /// Build a nonlinear building from explicit material-bearing RC sections.
    ///
    /// The section geometry must match the corresponding member in spec.
    /// Supply unscaled sections: spec.scale scales column fiber areas and
    /// column line mass together, exactly once. Beam sections are unscaled.
    /// Floor mass must exclude member mass. The same material-derived mass
    /// enters gravity, modal extraction and base-acceleration forcing, while
    /// each integration point receives its own virgin constitutive history.
    /// Original source selections remain available through [`Self::material_sections`].
    /// The legacy [`Self::new`] constructor remains an authored RC fixture.
    ///
    /// # Errors
    /// Invalid or mismatched specifications, already-scaled sections, or a
    /// failed gravity equilibrium/modal solve. No partially built frame is returned.
    pub fn with_sections(
        spec: BuildingSpec,
        column: RcSection,
        beam: RcSection,
        mass_formulation: FrameMassFormulation,
    ) -> Result<Self, SolidError> {
        validate_spec(spec)?;
        for (section, member) in [(&column, spec.column), (&beam, spec.beam)] {
            let geometry = section.geometry();
            if section.area_scale() != 1.0
                || geometry.depth_m != member.depth
                || geometry.width_m != member.width
                || geometry.steel_area_per_face_m2 != member.steel_area
            {
                return Err(SolidError::InvalidInput {
                    what: "material sections must be unscaled and match building member geometry"
                        .into(),
                });
            }
        }
        let column = column
            .scaled(spec.scale)
            .map_err(|error| SolidError::InvalidInput {
                what: format!("column material section: {error}"),
            })?;
        let line_mass = [column.line_mass_kg_m(), beam.line_mass_kg_m()];
        Self::assemble(
            spec,
            SectionModel::Fiber(column.section().clone()),
            SectionModel::Fiber(beam.section().clone()),
            line_mass,
            mass_formulation,
            Some([column, beam]),
        )
    }

    #[allow(clippy::needless_range_loop)]
    fn assemble(
        spec: BuildingSpec,
        column_section: SectionModel,
        beam_section: SectionModel,
        member_line_mass: [f64; 2],
        mass_formulation: FrameMassFormulation,
        material_sections: Option<[RcSection; 2]>,
    ) -> Result<Self, SolidError> {
        let mut frame = Frame2d::with_mass_formulation(Geometry::Corotational, mass_formulation);
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
        let column = || column_section.clone();
        let beam = || beam_section.clone();
        for s in 1..=spec.stories {
            for b in 0..=spec.bays {
                frame.add_element(
                    floors[s - 1][b],
                    floors[s][b],
                    np,
                    member_line_mass[0],
                    &column,
                )?;
                frame.add_mass(floors[s][b], node_mass);
            }
            for b in 0..spec.bays {
                frame.add_element(
                    floors[s][b],
                    floors[s][b + 1],
                    np,
                    member_line_mass[1],
                    &beam,
                )?;
            }
        }
        // The same reference mass projects floor/member gravity and base forcing.
        let acceleration: Vec<f64> = (0..frame.ndof())
            .map(|d| match d % 3 {
                0 => spec.gravity_m_s2[0],
                1 => spec.gravity_m_s2[1],
                _ => 0.0,
            })
            .collect();
        let gravity = frame.inertial_forces(&acceleration)?;
        frame.static_load(&gravity, 4)?;
        let modes = frame.modal_analysis()?;
        if modes
            .base_influence_mass_kg
            .iter()
            .any(|mass| !mass.is_finite() || *mass <= 0.0)
        {
            return Err(SolidError::InvalidInput {
                what: "building reference physical mass must be finite and positive".into(),
            });
        }
        if let Some(mode) = modes
            .modes
            .iter()
            .find(|m| m.status() != FrameModeStatus::Oscillatory)
        {
            return Err(SolidError::InvalidInput {
                what: format!(
                    "gravity-loaded tangent has {:?} mode: signed eigenvalue {}, resolution {}",
                    mode.status(),
                    mode.eigenvalue_per_s2,
                    mode.eigenvalue_resolution_per_s2
                ),
            });
        }
        // Uniform horizontal base motion couples through Gamma_x = phi^T M r_x.
        // Do not fit its damping to a lower vertical or unexcited member mode.
        let participation_floor = 256.0 * f64::EPSILON * modes.base_influence_mass_kg[0];
        let lateral: Vec<f64> = modes
            .modes
            .iter()
            .filter(|m| m.effective_mass_kg[0] > participation_floor)
            .map(|m| m.eigenvalue_per_s2.sqrt())
            .collect();
        let w1 = *lateral
            .first()
            .ok_or_else(|| SolidError::InternalInvariant {
                what: "building has no resolved horizontally participating mode".into(),
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
            modes,
            damping,
            gravity,
            material_sections,
        })
    }

    /// Periods of the initial gravity-loaded modes that participate in uniform
    /// horizontal base motion (s), descending. All directions remain available
    /// through [`Self::modal_analysis`].
    #[must_use]
    pub fn periods(&self) -> &[f64] {
        &self.periods
    }

    /// Total reference physical mass [kg]: floor point masses and all member
    /// mass, including the portions attached to fixed supports.
    ///
    /// This is the native mass owner's full rigid-translation inventory. The
    /// sum of modal effective masses can be smaller because supported inertia
    /// also contributes directly to the base reaction. Time histories do not
    /// change this reference inventory.
    #[must_use]
    pub fn reference_mass_kg(&self) -> f64 {
        self.modes.base_influence_mass_kg[0]
    }

    /// Complete signed, mass-normalized initial gravity-loaded modal analysis.
    /// A later nonlinear time history does not refresh this initial analysis.
    #[must_use]
    pub fn modal_analysis(&self) -> &FrameModalAnalysis {
        &self.modes
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

    /// Full held gravity load used by both preload and every time-history
    /// step. Each node contributes `[Fx, Fy, Mz]` in N, N, N m, including
    /// supported nodes and consistent member end moments. Node order is
    /// bottom floor first, then left-to-right within each floor. These are
    /// applied loads, not support reactions.
    #[must_use]
    pub fn gravity_forces_n(&self) -> &[f64] {
        &self.gravity
    }

    /// Accepted column/beam material sections, with the actual column area
    /// scale and unchanged source identities. None identifies the authored
    /// legacy fixture whose line masses, if any, were supplied independently.
    /// These are immutable virgin blueprints, not integration-point histories.
    #[must_use]
    pub fn material_sections(&self) -> Option<&[RcSection; 2]> {
        self.material_sections.as_ref()
    }

    /// Run a horizontal base-acceleration record `ag` sampled every `dt`.
    /// Retains the configured floor/member gravity throughout the history.
    /// Base excitation remains horizontal regardless of gravity direction.
    /// Consumes the gravity-loaded state; clone the frame to rerun.
    ///
    /// # Errors
    /// Solver errors if a step fails to converge (reduce `dt`).
    pub fn run(&mut self, ag: &[f64], dt: f64) -> Result<BuildingResponse, SolidError> {
        let u_grav = self.frame.displacements().to_vec();
        let hist = self
            .frame
            .newmark_with_load(ag, dt, (1.0, 0.0), self.damping, &self.gravity)?;
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
    let base = BuildingFrame::new(spec)?;
    building_fragility_from_frame(ensemble, &base, story, drift_limit, alpha, margin)
}

/// Run the building fragility study from one explicitly admitted frame state.
///
/// Every ensemble member and every MLMC fine correction starts from a clone
/// of base, preserving its material selections, mass formulation and initial
/// constitutive state. The caller's frame is unchanged, including on a failed
/// history. Pass a freshly gravity-loaded [`BuildingFrame::with_sections`]
/// result for a study of the initial material-bearing building. The coarse
/// story model is still an explicit approximation; it need not share the
/// fine model's material laws for the paired MLMC correction to telescope.
///
/// # Errors
/// Propagates a failed fine-model time history without returning partial results.
///
/// # Panics
/// The ensemble and probability controls have the same admission contract as
/// [`building_fragility`].
pub fn building_fragility_from_frame(
    ensemble: &StochasticEnsemble,
    base: &BuildingFrame,
    story: StoryParams,
    drift_limit: f64,
    alpha: f64,
    margin: f64,
) -> Result<BuildingFragility, SolidError> {
    assert_ground_motion_ensemble(ensemble);
    let dt = ensemble.dt.value;
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
