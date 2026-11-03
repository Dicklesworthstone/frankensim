//! Instantaneous physical impulses on the retained native passive-shaft state.

use super::super::{
    CarrierImpulse, RotorImpulses, RotorState, admit_derived_state, admit_impulse_component, array,
    dd, dot_dd, finite, impulse_kinematics, reject, vec,
};
use super::{ShaftMotion, native_refusal};
use crate::Refusal;
use fs_ga::Vec3;
use fs_time::lie::gyrostat::rotors::RotorImpulse as NativeImpulse;

/// Accepted zero-duration event, with original native momentum and work accounts.
///
/// No spring phase, pose, physical time or damping stage advances. Spring
/// potential is unchanged. Any impact/contact loss belongs to the law supplying
/// the impulses; negative impulse work is not relabelled as shaft heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShaftImpulse<const N: usize> {
    /// Complete accepted state in the original Flyer representation.
    pub state: RotorState<N>,
    /// Accepted relative shaft rates [rad/s].
    pub relative_speed_rad_s: [f64; N],
    /// Declared external linear impulse in world axes [N s].
    pub external_linear_impulse_world: [f64; 3],
    /// Declared spin-plus-orbital impulse about the supplied world origin [N m s].
    pub external_angular_impulse_world: [f64; 3],
    /// Internal motor/brake impulse times its native mean relative shaft rate [J].
    pub motor_work_j: f64,
    /// Carrier force/moment and aerodynamic shaft work at mean physical rates [J].
    pub external_work_j: f64,
    /// Actual translation plus native rotational kinetic increment [J].
    pub kinetic_energy_change_j: f64,
    /// Unchanged native spring potential before and after the event [J].
    pub elastic_energy_j: f64,
    /// Original native accepted angular impulse, shaft recoil and work receipt.
    pub native: NativeImpulse<N>,
    /// Actual kinetic increment minus authored motor and external work [J].
    /// Rounded transfer/velocity defects remain explicit, never physical heat.
    pub work_energy_residual_j: f64,
    /// Actual linear momentum change minus authored impulse [N s].
    pub linear_momentum_residual: [f64; 3],
    /// Actual native spin-plus-orbital change minus authored impulse [N m s].
    pub angular_momentum_residual: [f64; 3],
}

impl<const N: usize> ShaftMotion<N> {
    /// Apply simultaneous carrier and shaft impulses without a fictitious timestep.
    ///
    /// The native bank consumes the retained canonical momentum directly.
    /// Internal motor impulses change shaft momentum with carrier recoil;
    /// aerodynamic shaft impulses also enter total external momentum exactly
    /// once. Carrier impulses are about the assembly COM; use
    /// CarrierImpulse::at_body_point for a force impulse at a material point.
    ///
    /// Accepted pose and relative phases stay fixed. The prepared spring/drag
    /// model and its interval remain intact for the next advance. This method
    /// does not solve a collision law, engage a shaft lock or allocate heat.
    ///
    /// # Errors
    /// Native canonical recoil/energy refusal, nonfinite data or a transfer
    /// unresolved under the original Flyer impulse admission. There is no
    /// fallback to the distinct velocity-based impulse lane. A late refusal
    /// changes neither canonical momentum nor the published six-DOF state.
    #[allow(clippy::too_many_lines)] // Admit the complete event before publishing either state.
    pub fn apply_impulse(
        &mut self,
        impulses: RotorImpulses<N>,
        angular_origin_world: [f64; 3],
    ) -> Result<ShaftImpulse<N>, Refusal> {
        let origin = vec(angular_origin_world);
        if !finite(origin)
            || !finite(vec(impulses.carrier.linear_world_ns))
            || !finite(vec(impulses.carrier.angular_body_nms))
            || !impulses
                .motor_nms
                .iter()
                .chain(&impulses.aerodynamic_nms)
                .all(|value| value.is_finite())
        {
            return Err(reject("shaft impulses and angular origin must be finite"));
        }
        let shaft = core::array::from_fn(|i| impulses.motor_nms[i] + impulses.aerodynamic_nms[i]);
        let mut external = impulses.carrier.angular_body_nms.map(dd);
        let mut external_terms = impulses
            .carrier
            .angular_body_nms
            .map(|value| dd(value.abs()));
        for (rotor, aerodynamic) in self.assembly.rotors.iter().zip(impulses.aerodynamic_nms) {
            for j in 0..3 {
                let axis = rotor.axis_body[j];
                external[j] = external[j] + dd(axis) * dd(aerodynamic);
                external_terms[j] = external_terms[j] + dd(axis.abs()) * dd(aerodynamic.abs());
            }
        }
        let external = vec(external.map(|value| value.to_f64()));
        let native = self
            .bank
            .apply_impulse(&self.native_state, shaft, external)
            .map_err(native_refusal)?;
        // Retain Flyer's explicit-impulse resolution budget. Native work is
        // computed from actual rounded changes, not the requested increments.
        for i in 0..N {
            admit_impulse_component(
                shaft[i],
                native.axial_momentum_change_nms[i],
                native.axial_momentum_change_nms[i].abs(),
            )?;
        }
        for i in 0..3 {
            admit_impulse_component(
                array(external)[i],
                array(native.external_angular_impulse_body)[i],
                external_terms[i].to_f64(),
            )?;
        }

        let linear = vec(impulses.carrier.linear_world_ns);
        let (mut rigid, actual_linear, _) = impulse_kinematics(
            self.assembly.reduced_body.mass_kg,
            self.assembly.reduced_body.inertia_tensor(),
            &self.state.rigid,
            CarrierImpulse {
                linear_world_ns: impulses.carrier.linear_world_ns,
                angular_body_nms: [0.0; 3],
            },
            Vec3::default(),
        )?;
        let omega_after = self
            .bank
            .angular_velocity(&native.state_after)
            .map_err(native_refusal)?;
        rigid.omega_body = array(omega_after);
        admit_derived_state(&rigid, "native shaft impulse endpoint")?;
        let next = RotorState {
            rigid,
            axial_momentum: native.state_after.axial_momentum_nms,
            relative_phase_rad: native.state_after.relative_angle_rad,
        };
        let rates = self
            .bank
            .relative_rates(&native.state_after)
            .map_err(native_refusal)?;
        let elastic = self
            .stepper
            .potential_energy(&native.state_after)
            .map_err(native_refusal)?;
        let mean_velocity =
            vec(self.state.rigid.vel_mps).scale(0.5) + vec(rigid.vel_mps).scale(0.5);
        let mean_omega = vec(self.state.rigid.omega_body).scale(0.5) + omega_after.scale(0.5);
        let mut motor_work = dd(0.0);
        let mut external_work = dot_dd(linear, mean_velocity)
            + dot_dd(vec(impulses.carrier.angular_body_nms), mean_omega);
        for i in 0..N {
            motor_work =
                motor_work + dd(impulses.motor_nms[i]) * dd(native.mean_relative_rate_rad_s[i]);
            let absolute_mean = (0.5 * self.state.axial_momentum[i] + 0.5 * next.axial_momentum[i])
                / self.assembly.rotors[i].axial_inertia_kgm2;
            external_work = external_work + dd(impulses.aerodynamic_nms[i]) * dd(absolute_mean);
        }
        let kinetic_change =
            dot_dd(actual_linear, mean_velocity) + dd(native.kinetic_energy_change_j);
        let work_residual = kinetic_change - motor_work - external_work;
        let arm = vec(self.state.rigid.pos_m) - origin;
        let rotation = self.native_state.rotation;
        let angular_impulse = rotation
            .rotate(external)
            .map_err(|error| native_refusal(error.into()))?
            + arm.cross(linear);
        let linear_residual = actual_linear - linear;
        let spin_residual = vec(core::array::from_fn(|i| {
            (dd(array(native.external_angular_impulse_body)[i]) - dd(array(external)[i])).to_f64()
        }));
        let angular_residual = rotation
            .rotate(spin_residual)
            .map_err(|error| native_refusal(error.into()))?
            + arm.cross(linear_residual);
        if ![
            motor_work.to_f64(),
            external_work.to_f64(),
            kinetic_change.to_f64(),
            work_residual.to_f64(),
        ]
        .iter()
        .all(|value| value.is_finite())
            || !finite(angular_impulse)
            || !finite(linear_residual)
            || !finite(angular_residual)
        {
            return Err(reject(
                "native shaft impulse work or momentum accounting overflow",
            ));
        }

        let receipt = ShaftImpulse {
            state: next,
            relative_speed_rad_s: rates,
            external_linear_impulse_world: impulses.carrier.linear_world_ns,
            external_angular_impulse_world: array(angular_impulse),
            motor_work_j: motor_work.to_f64(),
            external_work_j: external_work.to_f64(),
            kinetic_energy_change_j: kinetic_change.to_f64(),
            elastic_energy_j: elastic,
            native,
            work_energy_residual_j: work_residual.to_f64(),
            linear_momentum_residual: array(linear_residual),
            angular_momentum_residual: array(angular_residual),
        };
        self.native_state = native.state_after;
        self.state = next;
        Ok(receipt)
    }
}
