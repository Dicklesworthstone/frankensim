//! The existing wing/BEMT owner driven by independent compliant shafts.
//!
//! Absolute axial momentum determines each physical propeller's speed. No
//! common gear coordinate or prescribed-speed shortcut couples the rotors.

use super::super::passive::{ShaftMotion, ShaftStep};
use super::super::{RotorAssembly, RotorState, RotorTorques};
use super::{
    AerodynamicSolution, AerodynamicStage, CoupledStep, DiskAxes, OpenLoopDesign,
    PropellerMount, PropulsionGeometry, Refusal, StageEnvironment, TICK_HZ, reject,
    surfaces::LiftingSurfaces,
};
use crate::propcoupling::swirl::{CoupledSwirlStep, WarmStart};
use crate::spine::{Loads, admitted_state_and_rotation};
use fs_airscrew::swirl::Controls;

/// Actual environment and independent motor ports at an aerodynamic stage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShaftInputs<const N: usize> {
    /// Positive air density [kg/m³].
    pub density_kg_m3: f64,
    /// Prescribed canard setting [rad], with the original aerodynamic domain.
    pub canard_rad: f64,
    /// Internal torque on each physical shaft, conjugate to its relative phase.
    pub motor_nm: [f64; N],
    /// World-NED gravitational acceleration [m/s²], included exactly once.
    pub gravity_world_mps2: [f64; 3],
}

/// The shared aerodynamic solution with independently owned shaft ports.
#[derive(Clone, Debug, PartialEq)]
pub struct ShaftEvaluation<const N: usize, S = CoupledStep> {
    /// World force and body moment about the assembly COM, excluding axial
    /// propeller torque; that torque belongs only to `torques.aerodynamic_nm`.
    pub carrier: Loads,
    /// Internal motor drive and external aerodynamic torque on each rotor.
    pub torques: RotorTorques<N>,
    /// Aerodynamic body force before gravity is added [N].
    pub aerodynamic_force_body_n: [f64; 3],
    /// Signed absolute L/R physical propeller rates from axial momentum [rad/s].
    pub propeller_absolute_speed_rad_s: [f64; 2],
    /// The same accepted native two-way wing/rotor solution as the loads.
    pub coupled: S,
}

/// Accepted spring/drag motion and its two actual aerodynamic evaluations.
#[derive(Clone, Debug, PartialEq)]
pub struct ShaftPropulsionStep<const N: usize, S = CoupledStep> {
    pub motion: ShaftStep<N>,
    /// Prediction loads; their mechanical work and heat are not accepted.
    pub initial: ShaftEvaluation<N, S>,
    /// Real midpoint loads used by the accepted native full interval.
    pub midpoint: ShaftEvaluation<N, S>,
}

/// Two aerodynamic disks attached to independent passive shafts and one carrier.
///
/// Geometry, mass and mechanical inventory remain fixed. The same wing/BEMT,
/// local wind, spatial disk wash and optional lifting surfaces as geared
/// propulsion supply forces; `ShaftMotion` owns springs, bearing drag, recoil,
/// shaft phases and accepted canonical momentum. This lane is quasi-steady:
/// it does not add dynamic wake storage or an ideal gear constraint.
pub struct ShaftPropulsion<'a, const N: usize> {
    geometry: PropulsionGeometry<'a>,
    assembly: &'a RotorAssembly<N>,
}

impl<'a, const N: usize> ShaftPropulsion<'a, N> {
    /// Bind the legacy body-x disk-wash domain to physical independent shafts.
    pub fn new(
        design: &'a OpenLoopDesign,
        assembly: &'a RotorAssembly<N>,
        mounts: [PropellerMount; 2],
    ) -> Result<Self, Refusal> {
        Self::build(design, assembly, mounts, None)
    }

    /// Bind canted shafts with the existing physical 3-D disk-overlap owner.
    pub fn new_spatial(
        design: &'a OpenLoopDesign,
        assembly: &'a RotorAssembly<N>,
        mounts: [PropellerMount; 2],
    ) -> Result<Self, Refusal> {
        let axes = DiskAxes::new(mounts.map(|mount| mount.thrust_direction_body))?;
        Self::build(design, assembly, mounts, Some(axes))
    }

    fn build(
        design: &'a OpenLoopDesign,
        assembly: &'a RotorAssembly<N>,
        mounts: [PropellerMount; 2],
        spatial_axes: Option<DiskAxes>,
    ) -> Result<Self, Refusal> {
        assembly.admit()?;
        let geometry = PropulsionGeometry {
            design,
            mounts,
            spatial_axes,
            surfaces: None,
        };
        geometry.admit(
            assembly.reduced_body.mass_kg,
            mounts.map(|mount| {
                assembly
                    .rotors
                    .get(mount.rotor_index)
                    .map(|rotor| rotor.axis_body)
            }),
            "shaft",
        )?;
        Ok(Self { geometry, assembly })
    }

    /// Use admitted native lifting surfaces, with the same spatial-wash domain
    /// and normal-section load interpretation as the geared consumer.
    pub fn with_lifting_surfaces(
        mut self,
        surfaces: &'a LiftingSurfaces<'a>,
    ) -> Result<Self, Refusal> {
        if self.geometry.spatial_axes.is_none() {
            return Err(reject(
                "lifting-surfaces-invalid",
                "custom surfaces require spatial disk wash",
            ));
        }
        self.geometry.surfaces = Some(surfaces);
        Ok(self)
    }

    /// Solve the existing wing/BEMT coupling at the actual independent spins.
    /// `wind_world` samples material air velocity at WORLD positions. A stopped
    /// or reversed blade spin refuses through the existing positive-speed BEMT
    /// domain; no rotor speed is clamped or inferred from another shaft.
    ///
    /// # Errors
    /// Mechanical, environmental, wind and native aerodynamic refusals. Inputs
    /// are immutable; no mechanical or physical wake state is advanced here.
    pub fn evaluate(
        &self,
        state: &RotorState<N>,
        inputs: ShaftInputs<N>,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        warm: Option<[f64; 2]>,
    ) -> Result<ShaftEvaluation<N>, Refusal> {
        if warm.is_some_and(|w| !w.iter().all(|x| x.is_finite())) {
            return Err(reject(
                "propulsion-state-invalid",
                "nonfinite warm induction",
            ));
        }
        self.evaluate_with(state, inputs, wind_world, |stage| {
            self.geometry.equilibrium(stage, warm)
        })
    }

    fn evaluate_with<S: AerodynamicSolution>(
        &self,
        state: &RotorState<N>,
        inputs: ShaftInputs<N>,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        solve: impl FnOnce(AerodynamicStage<'_>) -> Result<S, Refusal>,
    ) -> Result<ShaftEvaluation<N, S>, Refusal> {
        let (rigid, rotation) =
            admitted_state_and_rotation(&self.assembly.reduced_body, &state.rigid, 1.0 / TICK_HZ)?;
        if !state
            .axial_momentum
            .iter()
            .chain(&state.relative_phase_rad)
            .chain(&inputs.motor_nm)
            .all(|x| x.is_finite())
        {
            return Err(reject(
                "propulsion-state-invalid",
                "nonfinite shaft state or motor torque",
            ));
        }
        let absolute = self.geometry.mounts.map(|mount| {
            state.axial_momentum[mount.rotor_index]
                / self.assembly.rotors[mount.rotor_index].axial_inertia_kgm2
        });
        let result = self.geometry.evaluate_with::<N, _>(
            &rigid,
            rotation,
            self.assembly.reduced_body.mass_kg,
            absolute,
            StageEnvironment {
                density_kg_m3: inputs.density_kg_m3,
                canard_rad: inputs.canard_rad,
                gravity_world_mps2: inputs.gravity_world_mps2,
            },
            wind_world,
            solve,
        )?;
        Ok(ShaftEvaluation {
            carrier: result.carrier,
            torques: RotorTorques {
                motor_nm: inputs.motor_nm,
                aerodynamic_nm: result.aerodynamic_nm,
            },
            aerodynamic_force_body_n: result.aerodynamic_force_body_n,
            propeller_absolute_speed_rad_s: result.propeller_absolute_speed_rad_s,
            coupled: result.aerodynamics,
        })
    }

    /// Advance the bound native spring/drag owner using two real aerodynamic
    /// solves. The first wash only warms the midpoint solve. Stage callbacks
    /// see current carrier motion and each shaft's actual momentum and phase,
    /// including after an impulse. Both stored mechanical states commit once.
    ///
    /// # Errors
    /// A different mechanical inventory refuses before callbacks. Any initial,
    /// midpoint or accepted-step refusal leaves `motion` unchanged; callback
    /// side effects remain the caller's responsibility. The mechanical source
    /// accounts do not claim total air/body energy closure for quasi-steady BEMT.
    #[allow(clippy::too_many_arguments)]
    pub fn step<C, W>(
        &self,
        motion: &mut ShaftMotion<N>,
        t_s: f64,
        angular_origin_world: [f64; 3],
        warm: Option<[f64; 2]>,
        inputs: C,
        wind_world: W,
    ) -> Result<ShaftPropulsionStep<N>, Refusal>
    where
        C: FnMut(f64, &RotorState<N>) -> Result<ShaftInputs<N>, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
    {
        let mut stage_warm = warm;
        self.step_with(
            motion,
            t_s,
            angular_origin_world,
            inputs,
            wind_world,
            |state, inputs, wind| {
                let evaluation = self.evaluate(state, inputs, wind, stage_warm)?;
                stage_warm = Some(evaluation.coupled.w_slip);
                Ok(evaluation)
            },
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_with<S, C, W, E>(
        &self,
        motion: &mut ShaftMotion<N>,
        t_s: f64,
        angular_origin_world: [f64; 3],
        mut inputs: C,
        wind_world: W,
        mut evaluate: E,
    ) -> Result<ShaftPropulsionStep<N, S>, Refusal>
    where
        C: FnMut(f64, &RotorState<N>) -> Result<ShaftInputs<N>, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
        E: FnMut(
            &RotorState<N>,
            ShaftInputs<N>,
            &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        ) -> Result<ShaftEvaluation<N, S>, Refusal>,
    {
        if motion.assembly() != self.assembly {
            return Err(reject(
                "propulsion-assembly-invalid",
                "prepared shaft motion has a different physical inventory",
            ));
        }
        let mut initial = None;
        let mut midpoint = None;
        let accepted = motion.advance_with(t_s, angular_origin_world, |time, state| {
            let evaluation = evaluate(
                state,
                inputs(time, state)?,
                &|point| wind_world(time, point),
            )?;
            let loads = (evaluation.carrier, evaluation.torques);
            if initial.is_none() {
                initial = Some(evaluation);
            } else {
                midpoint = Some(evaluation);
            }
            Ok(loads)
        })?;
        // Successful advance_with has completed exactly two accepted callbacks.
        // There is no new fallible operation after its atomic state publication.
        Ok(ShaftPropulsionStep {
            motion: accepted,
            initial: initial.expect("successful native initial load evaluation"),
            midpoint: midpoint.expect("successful native midpoint load evaluation"),
        })
    }
}

/// Native annular axial/angular wake momentum on independent passive shafts.
///
/// The selected rotor law changes the actual thrust, circulation and signed
/// aerodynamic torque. Springs, bearing drag, carrier recoil and absolute
/// shaft work remain owned by the same prepared `ShaftMotion`. Selecting this
/// model does not impose an ideal gear or change ordinary `ShaftPropulsion`.
pub struct SwirlShaftPropulsion<'p, 'a, const N: usize> {
    propulsion: &'p ShaftPropulsion<'a, N>,
    controls: Controls,
}

impl<'p, 'a, const N: usize> SwirlShaftPropulsion<'p, 'a, N> {
    /// Select the existing native bounded swirl law on an admitted assembly.
    /// The rotor solver admits its controls before publishing any load state.
    #[must_use]
    pub fn new(propulsion: &'p ShaftPropulsion<'a, N>, controls: Controls) -> Self {
        Self {
            propulsion,
            controls,
        }
    }

    /// Evaluate both original annular momentum equations at each shaft's
    /// actual signed absolute speed. Required torque reaches its external
    /// shaft port once, retaining counter-rotation and regenerative power.
    ///
    /// Wing and mutual-disk wash retain the original axial actuator-disk
    /// field. Numerical annulus seeds add no physical wake storage.
    ///
    /// # Errors
    /// The shared assembly/material-air admission and native wing, swirl and
    /// coupling refusals pass through without changing mechanical state.
    pub fn evaluate(
        &self,
        state: &RotorState<N>,
        inputs: ShaftInputs<N>,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        warm: Option<&WarmStart>,
    ) -> Result<ShaftEvaluation<N, CoupledSwirlStep>, Refusal> {
        self.propulsion
            .evaluate_with(state, inputs, wind_world, |stage| {
                self.propulsion
                    .geometry
                    .swirling(stage, self.controls, warm)
            })
    }

    /// Advance the independent shafts and carrier from real initial and
    /// predicted-midpoint wing/swirl loads. Only a complete accepted stage
    /// becomes the next numerical branch seed. Mechanical momentum, spring
    /// storage, motor work and bearing heat come from the native full step.
    ///
    /// # Errors
    /// Inventory mismatch refuses before callbacks. Initial, midpoint and
    /// accepted-step failures preserve both stored mechanical states; caller
    /// callback side effects remain outside that transaction. This steady
    /// aerodynamic model adds no full fluid/body energy closure claim.
    #[allow(clippy::too_many_arguments)]
    pub fn step<C, W>(
        &self,
        motion: &mut ShaftMotion<N>,
        t_s: f64,
        angular_origin_world: [f64; 3],
        warm: Option<&WarmStart>,
        inputs: C,
        wind_world: W,
    ) -> Result<ShaftPropulsionStep<N, CoupledSwirlStep>, Refusal>
    where
        C: FnMut(f64, &RotorState<N>) -> Result<ShaftInputs<N>, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
    {
        let mut stage_warm = warm.cloned();
        self.propulsion.step_with(
            motion,
            t_s,
            angular_origin_world,
            inputs,
            wind_world,
            |state, inputs, wind| {
                let evaluation = self.evaluate(state, inputs, wind, stage_warm.as_ref())?;
                stage_warm = Some(evaluation.coupled.warm_start());
                Ok(evaluation)
            },
        )
    }
}
