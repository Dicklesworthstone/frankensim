//! Native passive shaft dynamics coupled to the actual six-DOF Flyer carrier.
//!
//! Preparation converts the supplied absolute shaft momenta once, then retains
//! the native canonical state across accepted steps. Springs and viscous drag
//! act between shafts and carrier; their recoil, phase and physical energy
//! account come from the same fs-time owner. Translation uses the original
//! Flyer impulse admission and drift/kick/drift map.
//!
//! This lane requires total angular momentum to resolve the carrier velocity.
//! It never falls back to a different state representation. The existing
//! instantaneous and prescribed held-torque Flyer APIs keep their own domains.

use super::{
    CarrierImpulse, Loads, RotorAssembly, RotorState, RotorTorques, admit_derived_state,
    admit_loads, admitted_state_and_rotation, array, dd, dot_dd, finite, impulse_kinematics,
    reject, rotation_to_boundary, vec,
};
use crate::Refusal;
use fs_ga::Vec3;
use fs_time::lie::LieStepError;
use fs_time::lie::gyrostat::{
    HeldTorque,
    rotors::{
        AxialRotor as NativeRotor, RotorBank, RotorState as NativeState, TorsionalRotorStep,
        TorsionalRotorStepper,
    },
};

#[path = "passive/impulse.rs"]
mod impulse;
pub use impulse::ShaftImpulse;

/// Fixed, independent spring/drag laws on the existing relative shaft phases.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShaftProperties<const N: usize> {
    /// Nonnegative rotor-to-carrier torsional stiffness [N m/rad].
    pub stiffness_nm_per_rad: [f64; N],
    /// Finite, unwrapped relative spring rest phases [rad].
    pub rest_phase_rad: [f64; N],
    /// Nonnegative relative-shaft viscous drag [N m s/rad].
    pub damping_nms_per_rad: [f64; N],
}
impl<const N: usize> Default for ShaftProperties<N> {
    fn default() -> Self {
        Self {
            stiffness_nm_per_rad: [0.0; N],
            rest_phase_rad: [0.0; N],
            damping_nms_per_rad: [0.0; N],
        }
    }
}

/// Prepared native spring/drag evolution and its accepted six-DOF state.
///
/// Native canonical momentum is authoritative and is never reconstructed from
/// the published carrier velocity at later steps. A failed advance changes
/// neither this momentum nor the published Flyer state. Preparation may
/// allocate the native 3-by-3 factor; successful advances allocate no storage.
#[derive(Debug)]
pub struct ShaftMotion<const N: usize> {
    assembly: RotorAssembly<N>,
    bank: RotorBank<N>,
    stepper: TorsionalRotorStepper<N>,
    native_state: NativeState<N>,
    state: RotorState<N>,
}

/// Actual accepted carrier/shaft motion with source-native storage and losses.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShaftStep<const N: usize> {
    /// Accepted complete state in the original Flyer representation.
    pub state: RotorState<N>,
    /// Accepted relative shaft rates [rad/s].
    pub relative_speed_rad_s: [f64; N],
    /// Held external force integrated over the interval, in world axes [N s].
    pub external_linear_impulse_world: [f64; 3],
    /// Includes spin and orbital angular impulse about the supplied world origin.
    pub external_angular_impulse_world: [f64; 3],
    /// Authored internal motor torque paired with accepted relative phase change.
    pub motor_work_j: f64,
    /// Force work plus native external torque work and aerodynamic shaft work.
    /// The aerodynamic part pairs its torque with ABSOLUTE shaft motion:
    /// its relative-phase work here complements its carrier work in the native receipt.
    pub external_work_j: f64,
    /// Translation plus native rotational kinetic change, evaluated by increments.
    pub kinetic_energy_change_j: f64,
    /// Native spring potential at the accepted interval start [J].
    pub elastic_energy_before_j: f64,
    /// Native spring potential at the accepted interval end [J].
    pub elastic_energy_after_j: f64,
    /// Original native receipt, including actual shaft impulses, per-shaft heat,
    /// spring-energy changes, and the separate gyrostat drift energy defect.
    /// Its drive work combines the supplied motor and aerodynamic shaft drives.
    pub native: TorsionalRotorStep<N>,
    /// Delta kinetic + Delta elastic + heat - motor - external - native drift.
    /// Includes phase/torque quadrature rounding; it is never physical heat.
    pub work_energy_residual_j: f64,
    /// Actual linear momentum change minus applied force impulse [N s].
    pub linear_momentum_residual: [f64; 3],
    /// Actual spin-plus-orbital change minus external angular impulse [N m s].
    pub angular_momentum_residual: [f64; 3],
}

impl<const N: usize> RotorAssembly<N> {
    /// Bind actual Flyer state to the native passive-shaft time owner.
    ///
    /// All stiffness zero selects native viscous evolution; all stiffness and
    /// drag zero selects its held-torque map. This does not replace the distinct
    /// endpoint-prescribed step_held map. Axes and absolute momenta are retained
    /// without normalization or a conversion through relative shaft speed.
    ///
    /// Duration obeys the original Flyer timestep domain and the native
    /// spring/drag resolution bound. Gravity remains an explicit force supplied
    /// once by the caller. Reprepare after changing inertia, axes, coefficients
    /// or duration. Electrical, gear, gimbal and fluid states are not inferred.
    ///
    /// # Errors
    /// Original Flyer admission, unresolved canonical carrier momentum, invalid
    /// shaft coefficients, or an interval exceeding native spring/drag bounds.
    pub fn shaft_motion(
        &self,
        initial: &RotorState<N>,
        properties: ShaftProperties<N>,
        duration_s: f64,
    ) -> Result<ShaftMotion<N>, Refusal> {
        self.admit()?;
        let (mut rigid, rotation) =
            admitted_state_and_rotation(&self.reduced_body, &initial.rigid, duration_s)?;
        let placeholder =
            NativeRotor::from_unit_axis(Vec3::new(1.0, 0.0, 0.0), 1.0).map_err(native_refusal)?;
        let mut rotors = [placeholder; N];
        for (target, source) in rotors.iter_mut().zip(self.rotors) {
            *target = NativeRotor::from_unit_axis(vec(source.axis_body), source.axial_inertia_kgm2)
                .map_err(native_refusal)?;
        }
        let bank = RotorBank::from_reduced(self.reduced_body.inertia_tensor(), rotors)
            .map_err(native_refusal)?;
        let native_state = bank
            .state_from_axial_momenta(
                rotation,
                vec(rigid.omega_body),
                initial.axial_momentum,
                initial.relative_phase_rad,
            )
            .map_err(native_refusal)?;
        let stepper = bank
            .torsional_stepper(
                properties.stiffness_nm_per_rad,
                properties.rest_phase_rad,
                properties.damping_nms_per_rad,
                duration_s,
            )
            .map_err(native_refusal)?;
        stepper
            .potential_energy(&native_state)
            .map_err(native_refusal)?;
        rigid.omega_body = array(
            bank.angular_velocity(&native_state)
                .map_err(native_refusal)?,
        );
        Ok(ShaftMotion {
            assembly: *self,
            bank,
            stepper,
            native_state,
            state: RotorState { rigid, ..*initial },
        })
    }
}

impl<const N: usize> ShaftMotion<N> {
    /// The accepted state in the existing Flyer representation.
    pub const fn state(&self) -> &RotorState<N> {
        &self.state
    }

    /// Original canonical state retained without lossy velocity round-trips.
    pub const fn native_state(&self) -> &NativeState<N> {
        &self.native_state
    }

    /// Prepared positive full-step interval [s].
    pub fn duration_s(&self) -> f64 {
        self.stepper.duration_s()
    }

    /// Native trace-based frequency estimate, not a convergence certificate.
    pub fn angular_frequency_bound_rad_s(&self) -> f64 {
        self.stepper.angular_frequency_bound_rad_s()
    }

    /// Prepared native spring/drag half-stage resolution number, at most two.
    pub fn half_stage_resolution_number(&self) -> f64 {
        self.stepper.half_stage_resolution_number()
    }

    /// Advance held loads and publish the complete carrier/shaft state atomically.
    ///
    /// Force is world-frame and carrier torque is body-frame, excluding the
    /// separately supplied aerodynamic rotor torques. Internal motor drive
    /// changes shaft momentum; aerodynamic drive changes both shaft and total
    /// momentum exactly once. Native phases, spring storage, viscous heat and
    /// drift error remain distinct. State-dependent load laws require the
    /// caller's coupled stage and timestep refinement.
    ///
    /// # Errors
    /// Nonfinite loads, native shaft/phase/recoil refusal, unresolved translational
    /// impulse, or nonfinite endpoint/accounting. The accepted state is unchanged.
    #[allow(clippy::too_many_lines)] // Stage all physical accounts before atomic publication.
    pub fn advance(
        &mut self,
        loads: Loads,
        torques: RotorTorques<N>,
        angular_origin_world: [f64; 3],
    ) -> Result<ShaftStep<N>, Refusal> {
        let loads = admit_loads(loads, "native shaft")?;
        let origin = vec(angular_origin_world);
        if !finite(origin)
            || !torques
                .motor_nm
                .iter()
                .chain(&torques.aerodynamic_nm)
                .all(|value| value.is_finite())
        {
            return Err(reject("shaft torques and angular origin must be finite"));
        }
        let dt = self.duration_s();
        let drive = core::array::from_fn(|i| torques.motor_nm[i] + torques.aerodynamic_nm[i]);
        let mut external = loads.moment_nm.map(dd);
        for (rotor, aerodynamic) in self.assembly.rotors.iter().zip(torques.aerodynamic_nm) {
            for (component, axis) in external.iter_mut().zip(rotor.axis_body) {
                *component = *component + dd(axis) * dd(aerodynamic);
            }
        }
        let external = vec(external.map(|value| value.to_f64()));
        let native = self
            .stepper
            .step(&self.native_state, drive, HeldTorque::Body(external))
            .map_err(native_refusal)?;
        let impulse = vec(loads.force_n).scale(dt);
        for (force, applied) in loads.force_n.into_iter().zip(array(impulse)) {
            if !applied.is_finite() || (force != 0.0 && applied == 0.0) {
                return Err(super::impulse_refusal());
            }
        }
        let (mut rigid, actual_linear, _) = impulse_kinematics(
            self.assembly.reduced_body.mass_kg,
            self.assembly.reduced_body.inertia_tensor(),
            &self.state.rigid,
            CarrierImpulse {
                linear_world_ns: array(impulse),
                angular_body_nms: [0.0; 3],
            },
            Vec3::default(),
        )?;
        let v0 = vec(self.state.rigid.vel_mps);
        let v1 = vec(rigid.vel_mps);
        let mean_velocity = v0.scale(0.5) + v1.scale(0.5);
        let midpoint_position = vec(self.state.rigid.pos_m) + v0.scale(0.5 * dt);
        rigid.pos_m = array(midpoint_position + v1.scale(0.5 * dt));
        rigid.quat = rotation_to_boundary(native.motion.state_after.rotation);
        rigid.omega_body = array(
            self.bank
                .angular_velocity(&native.motion.state_after)
                .map_err(native_refusal)?,
        );
        admit_derived_state(&rigid, "native shaft endpoint")?;
        let next = RotorState {
            rigid,
            axial_momentum: native.motion.state_after.axial_momentum_nms,
            relative_phase_rad: native.motion.state_after.relative_angle_rad,
        };
        let rates = self
            .bank
            .relative_rates(&native.motion.state_after)
            .map_err(native_refusal)?;
        let before_elastic = self
            .stepper
            .potential_energy(&self.native_state)
            .map_err(native_refusal)?;
        let after_elastic = self
            .stepper
            .potential_energy(&native.motion.state_after)
            .map_err(native_refusal)?;
        let mut motor_work = dd(0.0);
        let mut external_work = dot_dd(impulse, mean_velocity) + dd(native.motion.external_work_j);
        for i in 0..N {
            let phase = dd(next.relative_phase_rad[i]) - dd(self.state.relative_phase_rad[i]);
            motor_work = motor_work + dd(torques.motor_nm[i]) * phase;
            external_work = external_work + dd(torques.aerodynamic_nm[i]) * phase;
        }
        let kinetic_change =
            dot_dd(actual_linear, mean_velocity) + dd(native.motion.kinetic_energy_change_j);
        let mut energy_residual =
            kinetic_change - motor_work - external_work - dd(native.motion.drift_energy_defect_j);
        for i in 0..N {
            energy_residual = energy_residual
                + dd(native.elastic_energy_change_j[i])
                + dd(native.dissipated_energy_j[i]);
        }

        let orbital_impulse = (midpoint_position - origin).cross(impulse);
        let angular_impulse = native.motion.external_angular_impulse_world + orbital_impulse;
        let linear_residual = actual_linear - impulse;
        // Incremental orbital balance avoids subtracting two large scene totals.
        let position_change = vec(core::array::from_fn(|i| {
            (dd(next.rigid.pos_m[i]) - dd(self.state.rigid.pos_m[i])).to_f64()
        }));
        let angular_residual = native.motion.angular_momentum_residual_world
            + (vec(self.state.rigid.pos_m) - origin).cross(actual_linear)
            + position_change.cross(v1.scale(self.assembly.reduced_body.mass_kg))
            - orbital_impulse;
        if ![
            motor_work.to_f64(),
            external_work.to_f64(),
            kinetic_change.to_f64(),
            energy_residual.to_f64(),
        ]
        .iter()
        .all(|value| value.is_finite())
            || !finite(angular_impulse)
            || !finite(linear_residual)
            || !finite(angular_residual)
        {
            return Err(reject("native shaft work or momentum accounting overflow"));
        }
        let receipt = ShaftStep {
            state: next,
            relative_speed_rad_s: rates,
            external_linear_impulse_world: array(impulse),
            external_angular_impulse_world: array(angular_impulse),
            motor_work_j: motor_work.to_f64(),
            external_work_j: external_work.to_f64(),
            kinetic_energy_change_j: kinetic_change.to_f64(),
            elastic_energy_before_j: before_elastic,
            elastic_energy_after_j: after_elastic,
            native,
            work_energy_residual_j: energy_residual.to_f64(),
            linear_momentum_residual: array(linear_residual),
            angular_momentum_residual: array(angular_residual),
        };
        self.native_state = native.motion.state_after;
        self.state = next;
        Ok(receipt)
    }
}

fn native_refusal(error: LieStepError) -> Refusal {
    Refusal {
        code: "rotor-shaft-refused",
        message: format!("native shaft dynamics refused: {error}"),
        ranked_repairs: vec![
            "resolve carrier momentum and shaft phase; reduce the interval for the supplied spring and drag"
                .into(),
        ],
    }
}
