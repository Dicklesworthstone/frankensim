//! Native mean/first-harmonic and skew wake storage on independent shafts.
//!
//! Fixed disk-frame modes feed the same prepared wing/blade evaluator as geared
//! flight. The original modal transport and its signed work stay separate from
//! physical shaft momentum, spring storage, bearing heat and carrier recoil.

use super::super::super::dynamic::BladeIncident;
use super::super::super::dynamic::first_harmonic::{HarmonicAerodynamics, HarmonicWake};
use super::{
    Refusal, RotorState, ShaftEvaluation, ShaftInputs, ShaftMotion, ShaftPropulsion, ShaftStep,
    reject,
};
use fs_airscrew::dynamic_inflow::first_harmonic::skew::{SkewInflowControls, SkewInflowStep};
use fs_airscrew::dynamic_inflow::first_harmonic::{
    FirstHarmonicInflow, HarmonicInduction, HarmonicInflowStep,
};
use fs_airscrew::oblique::DiskFrame;

/// Actual carrier/independent-shaft stage and all six physical wake modes.
/// Modes use the fixed body disk frames returned by `disk_frames`; they never
/// rotate into the instantaneous wind direction or follow blade phase.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicShaftState<const N: usize> {
    pub motion: RotorState<N>,
    pub induction: [HarmonicInduction; 2],
}

/// Accepted native mechanics and modal-fluid accounts for the same interval.
#[derive(Clone, Debug, PartialEq)]
pub struct HarmonicShaftStep<const N: usize> {
    pub state: HarmonicShaftState<N>,
    pub motion: ShaftStep<N>,
    pub inflow: [HarmonicInflowStep; 2],
    /// The same modal steps as `inflow`, with source skew transport and signed
    /// cross work exposed. These are never credited as additional wake updates.
    pub skew_inflow: Option<[SkewInflowStep; 2]>,
    pub initial: ShaftEvaluation<N, HarmonicAerodynamics>,
    pub midpoint: ShaftEvaluation<N, HarmonicAerodynamics>,
    /// Implicit wake input work minus actual midpoint blade induced work.
    /// This partitioned quadrature defect is not bearing heat or dissipation.
    pub coupling_work_residual_j: [f64; 2],
}

/// Physical asymmetric wakes over the retained compliant-shaft mechanical owner.
///
/// Each rotor receives its own signed absolute `p/J`; transverse hub moments
/// use the carrier port and axial aerodynamic torque uses the external shaft
/// port once. The original native modal apparent masses never enter aircraft
/// inertia. This is finite-state axial/modal transport, not resolved wake world
/// momentum, dynamic swirl, an unsteady wing or a whole-fluid passivity theorem.
pub struct HarmonicShaftPropulsion<'p, 'a, const N: usize> {
    propulsion: &'p ShaftPropulsion<'a, N>,
    wake: HarmonicWake<'a>,
}

impl<'p, 'a, const N: usize> HarmonicShaftPropulsion<'p, 'a, N> {
    /// Prepare native normal-flow first harmonics and oblique blade quadrature.
    pub fn new(
        propulsion: &'p ShaftPropulsion<'a, N>,
        density_kg_m3: f64,
        azimuth_points: usize,
    ) -> Result<Self, Refusal> {
        let axes = propulsion
            .geometry
            .mounts
            .map(|mount| propulsion.assembly.rotors[mount.rotor_index].axis_body);
        Ok(Self {
            propulsion,
            wake: HarmonicWake::new(propulsion.geometry, axes, density_kg_m3, azimuth_points)?,
        })
    }

    /// Select the original bounded Pitt--Peters skew transport. Its nonsymmetric
    /// gain and signed mean/harmonic cross work retain their native owner.
    pub fn new_skewed(
        propulsion: &'p ShaftPropulsion<'a, N>,
        density_kg_m3: f64,
        azimuth_points: usize,
        controls: SkewInflowControls,
    ) -> Result<Self, Refusal> {
        let axes = propulsion
            .geometry
            .mounts
            .map(|mount| propulsion.assembly.rotors[mount.rotor_index].axis_body);
        Ok(Self {
            propulsion,
            wake: HarmonicWake::new_skewed(
                propulsion.geometry,
                axes,
                density_kg_m3,
                azimuth_points,
                controls,
            )?,
        })
    }

    /// Use the existing pointwise material-wind or wind/wing field at actual
    /// blade positions, with the same accepted circulation and self-wash owner.
    pub fn with_blade_incident(mut self, incident: BladeIncident) -> Result<Self, Refusal> {
        self.wake = self.wake.with_blade_incident(incident)?;
        Ok(self)
    }

    #[must_use]
    pub fn blade_incident(&self) -> BladeIncident {
        self.wake.blade_incident()
    }

    #[must_use]
    pub fn disk_frames(&self) -> [DiskFrame; 2] {
        self.wake.disk_frames()
    }

    #[must_use]
    pub fn inflow_model(&self) -> FirstHarmonicInflow {
        self.wake.inflow_model()
    }

    /// Evaluate the actual independent spins and all stored modes without
    /// advancing memory or replacing it with a steady equilibrium solution.
    ///
    /// # Errors
    /// The native geometry, shaft, whole-disk throughflow, wing, blade and skew
    /// domains pass through unchanged. Density must remain fixed and modes finite.
    pub fn evaluate(
        &self,
        state: &HarmonicShaftState<N>,
        inputs: ShaftInputs<N>,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
    ) -> Result<ShaftEvaluation<N, HarmonicAerodynamics>, Refusal> {
        self.wake
            .admit_state(inputs.density_kg_m3, state.induction)?;
        self.propulsion
            .evaluate_with(&state.motion, inputs, wind_world, |stage| {
                self.wake.solve(stage, state.induction, None)
            })
    }

    /// Advance body, independent shafts and all six physical induction modes.
    /// Initial forces predict every state to half time. Actual midpoint forces
    /// then drive the original full native mechanical and modal wake maps.
    ///
    /// Full modal/skew proposals and finite coupling-work checks finish inside
    /// the second callback, before the canonical mechanical owner can publish.
    /// Predictor work, heat and transport are never accepted source accounts.
    ///
    /// # Errors
    /// Inventory mismatch refuses before callbacks. Any force, wake, accounting
    /// or late mechanical refusal preserves canonical and published mechanics
    /// and every caller-owned mode. Callback effects remain caller-owned.
    /// Smooth coupled trajectories require timestep refinement; no coupled
    /// passivity or total fluid/body energy closure follows from these receipts.
    #[allow(clippy::too_many_arguments)]
    pub fn step<C, W>(
        &self,
        motion: &mut ShaftMotion<N>,
        induction: &mut [HarmonicInduction; 2],
        t_s: f64,
        angular_origin_world: [f64; 3],
        mut inputs: C,
        wind_world: W,
    ) -> Result<HarmonicShaftStep<N>, Refusal>
    where
        C: FnMut(f64, &HarmonicShaftState<N>) -> Result<ShaftInputs<N>, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
    {
        if motion.assembly() != self.propulsion.assembly {
            return Err(reject(
                "propulsion-assembly-invalid",
                "prepared shaft motion has a different physical inventory",
            ));
        }
        let before = *induction;
        let dt_s = motion.duration_s();
        let mut stage_induction = before;
        let mut initial = None;
        let mut midpoint = None;
        let mut full_wake = None;
        let accepted = motion.advance_with(t_s, angular_origin_world, |time, mechanical| {
            let stage = HarmonicShaftState {
                motion: *mechanical,
                induction: stage_induction,
            };
            let evaluation = self.evaluate(&stage, inputs(time, &stage)?, &|point| {
                wind_world(time, point)
            })?;
            let loads = (evaluation.carrier, evaluation.torques);
            if initial.is_none() {
                stage_induction = self
                    .wake
                    .advance_wakes(before, &evaluation.coupled, 0.5 * dt_s)?
                    .0
                    .map(|wake| wake.induction);
                initial = Some(evaluation);
            } else {
                let (inflow, skew) = self.wake.advance_wakes(before, &evaluation.coupled, dt_s)?;
                let work =
                    HarmonicWake::coupling_work_residual(&inflow, &evaluation.coupled, dt_s)?;
                // Every fallible air/source operation precedes native mechanics.
                // A subsequent mechanical refusal discards this whole proposal.
                full_wake = Some((inflow, skew, work));
                midpoint = Some(evaluation);
            }
            Ok(loads)
        })?;
        // Native success guarantees exactly two completed stage callbacks.
        let (inflow, skew_inflow, coupling_work_residual_j) =
            full_wake.expect("successful native midpoint modal wake proposal");
        *induction = inflow.map(|wake| wake.induction);
        Ok(HarmonicShaftStep {
            state: HarmonicShaftState {
                motion: accepted.state,
                induction: *induction,
            },
            motion: accepted,
            inflow,
            skew_inflow,
            initial: initial.expect("successful native initial load evaluation"),
            midpoint: midpoint.expect("successful native midpoint load evaluation"),
            coupling_work_residual_j,
        })
    }
}
