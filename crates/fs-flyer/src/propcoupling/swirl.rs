//! Two-way wing/propeller coupling with native angular wake momentum.
//!
//! The original coupling loop consumes the actual swirling rotor thrust and
//! torque at each trial disk inflow. Wing and inter-disk wash retain the
//! existing axial actuator-disk field; no transported helical wake or swirl
//! velocity acting on the wings is implied. Shaft torque remains signed.

use super::spatial::DiskAxes;
use super::{
    CoupledStep, IncidentInflow, InfluenceOperator, Panel, PropCouplingSolverSpec, PropDisk,
    Refusal, Rotor, RotorLoads, SectionClosure, SectionFrame, StripSpec, coupled_step_with_rotors,
};
use fs_airscrew::swirl::{self as rotor_swirl, Controls, Induction, Power, Solution};
use fs_blake3::hash_domain;

/// Explicit numerical continuation of a previous complete accepted solve.
/// These are seeds, not additional physical wake inertia or shaft momentum.
#[derive(Clone, Debug, PartialEq)]
pub struct WarmStart {
    /// Axial actuator-disk wash used by the original wing coupling [m/s].
    pub w_slip: [f64; 2],
    /// One original axial/tangential induction per radial station, L/R.
    pub induction: [Vec<Induction>; 2],
}

/// Actual swirling rotor states and their accepted two-way wing solution.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledSwirlStep {
    /// Accepted wing/wash state and loads, with explicit swirl power.
    pub coupled: CoupledStep<Power>,
    /// L/R original annular linear/angular momentum solutions at the SAME
    /// disk inflows as the accepted loads, not a later diagnostic re-solve.
    pub rotors: [Solution; 2],
}

impl CoupledSwirlStep {
    /// Preserve the accepted numerical branches for the next physical stage.
    #[must_use]
    pub fn warm_start(&self) -> WarmStart {
        WarmStart {
            w_slip: self.coupled.w_slip,
            induction: core::array::from_fn(|k| {
                self.rotors[k]
                    .stations
                    .iter()
                    .map(|s| s.induction)
                    .collect()
            }),
        }
    }
}

/// Select native axial AND angular BEMT in the existing two-way iteration.
///
/// `axes=None` retains body-x wash; supplied axes select the existing spatial
/// disk-overlap/vector-wash geometry. Distributed incident flow supplies the
/// explicit strip and disk velocities; shared admission refuses reversed or
/// nonfinite disk inflow.
/// The original correction cap, trial rejection, relaxation arithmetic and
/// per-disk acceptance remain. A rejected outer trial never becomes a new
/// station-branch anchor: all rotor calls use the caller's fixed seed.
///
/// # Errors
/// Native wing, swirl and coupling refusals pass through, without partial
/// loads. Invalid warm wash receives `swirl-coupling-invalid`; station-seed
/// admission is owned by native BEMT. This adds no high-induction correction.
#[allow(clippy::too_many_arguments)]
pub fn coupled_step(
    axes: Option<DiskAxes>,
    op: &InfluenceOperator,
    panels: &[Panel],
    strips: &[StripSpec],
    closure: SectionClosure<'_>,
    rotor: &Rotor,
    disks: &[PropDisk; 2],
    freestream: [f64; 3],
    rho: f64,
    spec: &PropCouplingSolverSpec,
    controls: Controls,
    warm: Option<&WarmStart>,
    incident: IncidentInflow<'_>,
) -> Result<CoupledSwirlStep, Refusal> {
    coupled_impl(
        axes, op, panels, strips, closure, None, rotor, disks, freestream, rho, spec, controls,
        warm, incident,
    )
}

/// Native axial/angular rotor momentum coupled to explicit local wing frames.
///
/// Spatial axes are mandatory: the complete vector wash reaches each section
/// before incidence/induction projection, including vertical and canted
/// surfaces. The same accepted circulation drives the actual swirling rotor
/// solves and the returned force velocities. Station warm starts and signed
/// shaft torque retain the original owner in [`coupled_step`].
///
/// # Errors
/// The refusals of [`coupled_step`] and the framed nonlinear wing model.
#[allow(clippy::too_many_arguments)]
pub fn coupled_step_with_section_frames(
    axes: &DiskAxes,
    op: &InfluenceOperator,
    panels: &[Panel],
    strips: &[StripSpec],
    closure: SectionClosure<'_>,
    frames: &[SectionFrame],
    rotor: &Rotor,
    disks: &[PropDisk; 2],
    freestream: [f64; 3],
    rho: f64,
    spec: &PropCouplingSolverSpec,
    controls: Controls,
    warm: Option<&WarmStart>,
    incident: IncidentInflow<'_>,
) -> Result<CoupledSwirlStep, Refusal> {
    coupled_impl(
        Some(*axes),
        op,
        panels,
        strips,
        closure,
        Some(frames),
        rotor,
        disks,
        freestream,
        rho,
        spec,
        controls,
        warm,
        incident,
    )
}

#[allow(clippy::too_many_arguments)]
fn coupled_impl(
    axes: Option<DiskAxes>,
    op: &InfluenceOperator,
    panels: &[Panel],
    strips: &[StripSpec],
    closure: SectionClosure<'_>,
    frames: Option<&[SectionFrame]>,
    rotor: &Rotor,
    disks: &[PropDisk; 2],
    freestream: [f64; 3],
    rho: f64,
    spec: &PropCouplingSolverSpec,
    controls: Controls,
    warm: Option<&WarmStart>,
    incident: IncidentInflow<'_>,
) -> Result<CoupledSwirlStep, Refusal> {
    if warm.is_some_and(|w| !w.w_slip.into_iter().all(f64::is_finite)) {
        return Err(Refusal {
            code: "swirl-coupling-invalid",
            message: "nonfinite axial wash seed".into(),
            ranked_repairs: vec!["supply a previous accepted coupling state".into()],
        });
    }
    let mut rotors = [None, None];
    let mut coupled = coupled_step_with_rotors(
        op,
        panels,
        strips,
        closure,
        rotor,
        disks,
        freestream,
        rho,
        spec,
        warm.map(|w| w.w_slip),
        Some(incident),
        axes,
        frames,
        |k, axial, omega| {
            let solution = rotor_swirl::bemt_solve(
                rotor,
                rho,
                axial,
                omega,
                controls,
                warm.map(|w| w.induction[k].as_slice()),
            )
            .map_err(|e| Refusal {
                code: e.code,
                message: e.message,
                ranked_repairs: e.ranked_repairs,
            })?;
            let loads = RotorLoads {
                thrust_n: solution.thrust_n,
                torque_nm: solution.torque_nm,
                power: solution.power,
            };
            rotors[k] = Some(solution);
            Ok(loads)
        },
    )?;
    // Model selection and the native bounded solve controls are identity
    // ingredients; the default axial coupling digest remains unchanged.
    let mut identity = coupled.spec_digest.as_bytes().to_vec();
    identity.extend_from_slice(&controls.maximum_station_evaluations.to_le_bytes());
    identity.extend_from_slice(&controls.relative_tolerance.to_bits().to_le_bytes());
    coupled.spec_digest =
        hash_domain("org.frankensim.fs-flyer.prop-coupling-swirl.v1", &identity).to_hex();
    let [Some(left), Some(right)] = rotors else {
        return Err(Refusal {
            code: "swirl-coupling-invalid",
            message: "accepted coupling omitted a rotor solution".into(),
            ranked_repairs: vec!["require a complete two-rotor solve".into()],
        });
    };
    Ok(CoupledSwirlStep {
        coupled,
        rotors: [left, right],
    })
}
