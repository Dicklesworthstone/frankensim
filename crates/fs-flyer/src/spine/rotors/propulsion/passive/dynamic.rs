//! Independent compliant shafts and the existing physical dynamic-inflow owner.
//!
//! The native shaft and wake predictors form one joint midpoint. Wake trials
//! finish before the retained mechanical state can commit; only the accepted
//! full interval contributes shaft work, bearing heat and wake transport.

use super::super::dynamic::{BladeIncident, DynamicAerodynamics, DynamicWake};
use super::{
    Refusal, RotorState, ShaftEvaluation, ShaftInputs, ShaftMotion, ShaftPropulsion, ShaftStep,
    reject,
};
use fs_airscrew::dynamic_inflow::{InflowStep, UniformAxialInflow};

/// Native first-harmonic and skew wake storage on independent shafts.
#[path = "dynamic/first_harmonic.rs"]
pub mod first_harmonic;

/// Physical snapshot supplied to a joint shaft/body/wake load law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DynamicShaftState<const N: usize> {
    /// Complete independent-shaft and six-DOF carrier state.
    pub motion: RotorState<N>,
    /// Actual L/R uniform disk-induced velocities [m/s], not solver guesses.
    pub induced_mps: [f64; 2],
}

/// Accepted compliant-shaft motion, physical wake memory and source accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct DynamicShaftStep<const N: usize> {
    /// Joint endpoint; the corresponding canonical mechanics remain in `ShaftMotion`.
    pub state: DynamicShaftState<N>,
    /// Native spring/drag motion, storage, heat and world momentum balance.
    pub motion: ShaftStep<N>,
    /// Per-disk apparent-fluid storage, axial impulse and induced-work balances.
    pub inflow: [InflowStep; 2],
    /// Real start loads used only for the joint midpoint prediction.
    pub initial: ShaftEvaluation<N, DynamicAerodynamics>,
    /// Real loads at the predicted body, independent shafts and physical wakes.
    pub midpoint: ShaftEvaluation<N, DynamicAerodynamics>,
}

/// Two evolving uniform wakes driven by independent native compliant shafts.
///
/// This shares the geared consumer's wing, blade, mutual-wash and wake-flux
/// owner. Each propeller receives its own actual absolute spin `p/J`; springs,
/// bearing drag, recoil and canonical momentum stay in `ShaftMotion`.
/// Apparent wake mass is excluded from the physical aircraft mass/inertia.
/// The per-disk axial storage/work balances are not resolved-wake world
/// momentum or a whole-fluid/body passivity theorem. The wing remains steady;
/// no harmonic wake state, skew transport or ideal gear constraint is added.
pub struct DynamicShaftPropulsion<'p, 'a, const N: usize> {
    propulsion: &'p ShaftPropulsion<'a, N>,
    wake: DynamicWake<'a>,
}

impl<'p, 'a, const N: usize> DynamicShaftPropulsion<'p, 'a, N> {
    /// Bind the original uniform axial inflow law at fixed positive density.
    /// Actual three-dimensional disk overlap and thrust axes are retained.
    pub fn new(
        propulsion: &'p ShaftPropulsion<'a, N>,
        density_kg_m3: f64,
    ) -> Result<Self, Refusal> {
        Ok(Self {
            propulsion,
            wake: DynamicWake::new(propulsion.geometry, density_kg_m3)?,
        })
    }

    /// Select the original oblique blade quadrature and carrier tilt-rate
    /// kinematics, with the same uniform axial wake-transport approximation.
    /// The actual rotor axes and blade handedness determine spin direction.
    pub fn new_oblique(
        propulsion: &'p ShaftPropulsion<'a, N>,
        density_kg_m3: f64,
        azimuth_points: usize,
    ) -> Result<Self, Refusal> {
        let rotor_axes = propulsion
            .geometry
            .mounts
            .map(|mount| propulsion.assembly.rotors[mount.rotor_index].axis_body);
        Ok(Self {
            propulsion,
            wake: DynamicWake::new_oblique(
                propulsion.geometry,
                rotor_axes,
                density_kg_m3,
                azimuth_points,
            )?,
        })
    }

    /// Select the same pointwise material-wind/wing field as the geared
    /// oblique owner. A non-hub choice requires native blade quadrature.
    pub fn with_blade_incident(mut self, incident: BladeIncident) -> Result<Self, Refusal> {
        self.wake = self.wake.with_blade_incident(incident)?;
        Ok(self)
    }

    #[must_use]
    pub fn blade_incident(&self) -> BladeIncident {
        self.wake.blade_incident()
    }

    /// The explicit physical per-disk apparent-air-mass model.
    #[must_use]
    pub fn inflow_model(&self) -> UniformAxialInflow {
        self.wake.inflow_model()
    }

    /// Evaluate the wing and both blade loads at the stored physical induction
    /// and independent shaft momenta, without advancing any mechanical or wake state.
    /// Wind samples are world material-air velocities at actual world positions.
    ///
    /// # Errors
    /// The shared mechanical, material-air and native aerodynamic refusals;
    /// changing the admitted wake density or nonfinite induction also refuses.
    pub fn evaluate(
        &self,
        state: &DynamicShaftState<N>,
        inputs: ShaftInputs<N>,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
    ) -> Result<ShaftEvaluation<N, DynamicAerodynamics>, Refusal> {
        self.wake
            .admit_state(inputs.density_kg_m3, state.induced_mps)?;
        self.propulsion
            .evaluate_with(&state.motion, inputs, wind_world, |stage| {
                self.wake.solve(stage, state.induced_mps)
            })
    }

    /// Advance the actual body, independent shafts and both physical wakes.
    /// The initial loads predict all coordinates to half time; the midpoint
    /// loads drive the original full-interval native shaft and implicit wake
    /// maps. This is a second-order partitioned midpoint composition for smooth
    /// admitted cases, with load-feedback timestep refinement required.
    ///
    /// # Errors
    /// Inventory mismatch refuses before callbacks. Initial, midpoint, wake
    /// or final mechanical accounting failures preserve both `motion` and
    /// `induced_mps`, including retained canonical momentum. Both wake trials
    /// finish inside the load callbacks, before mechanics can publish. The
    /// exclusive mutable borrows are retained through the final infallible
    /// induction assignment; callback side effects remain caller-owned.
    #[allow(clippy::too_many_arguments)]
    pub fn step<C, W>(
        &self,
        motion: &mut ShaftMotion<N>,
        induced_mps: &mut [f64; 2],
        t_s: f64,
        angular_origin_world: [f64; 3],
        mut inputs: C,
        wind_world: W,
    ) -> Result<DynamicShaftStep<N>, Refusal>
    where
        C: FnMut(f64, &DynamicShaftState<N>) -> Result<ShaftInputs<N>, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
    {
        if motion.assembly() != self.propulsion.assembly {
            return Err(reject(
                "propulsion-assembly-invalid",
                "prepared shaft motion has a different physical inventory",
            ));
        }
        let before = *induced_mps;
        let dt_s = motion.duration_s();
        let mut stage_induction = before;
        let mut initial = None;
        let mut midpoint = None;
        let mut inflow = None;
        let accepted = motion.advance_with(t_s, angular_origin_world, |time, mechanical| {
            let stage = DynamicShaftState {
                motion: *mechanical,
                induced_mps: stage_induction,
            };
            let evaluation = self.evaluate(&stage, inputs(time, &stage)?, &|point| {
                wind_world(time, point)
            })?;
            let loads = (evaluation.carrier, evaluation.torques);
            if initial.is_none() {
                stage_induction = self
                    .wake
                    .advance_wakes(before, &evaluation.coupled, 0.5 * dt_s)?
                    .map(|wake| wake.induced_mps);
                initial = Some(evaluation);
            } else {
                // This must precede the mutable native owner's publication.
                // A late failure in either disk leaves all accepted states intact.
                inflow = Some(self.wake.advance_wakes(before, &evaluation.coupled, dt_s)?);
                midpoint = Some(evaluation);
            }
            Ok(loads)
        })?;
        // Successful advance_with guarantees exactly two completed callbacks.
        // No Result-returning operation follows its canonical state publication.
        let inflow = inflow.expect("successful native midpoint wake trial");
        *induced_mps = inflow.map(|wake| wake.induced_mps);
        Ok(DynamicShaftStep {
            state: DynamicShaftState {
                motion: accepted.state,
                induced_mps: *induced_mps,
            },
            motion: accepted,
            inflow,
            initial: initial.expect("successful native initial load evaluation"),
            midpoint: midpoint.expect("successful native midpoint load evaluation"),
        })
    }
}
