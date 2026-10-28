//! Six-DOF carrier motion with independent, body-fixed axial rotors.
//!
//! Motors exchange angular momentum internally; aerodynamic shaft torque
//! is external to the combined body/rotor system. Both affect rotor spin,
//! but only the latter belongs in the total external torque. The native
//! gyrostat map also includes gyroscopic precession without a lagged
//! `omega cross rotor_momentum` force or a separately invented recoil.
//!
//! This is a held-load integrator. It does not change registered SimLoop
//! presets or infer rotor inertia, shaft direction, gearing, or a load law.

use super::{
    Loads, RigidBody, SixDofState, TICK_HZ, admit_and_canonicalize_state, admit_derived_state,
    admit_loads, admitted_state_and_rotation, rotation_to_boundary,
};
use crate::Refusal;
use fs_ga::Vec3;
use fs_time::lie::gyrostat::{HeldTorque, driven_step};

/// Fixed per-tick rotor-count bound; successful steps allocate no heap storage.
pub const MAX_ROTORS: usize = 64;

/// Axisymmetric rotor mounted at a fixed orientation in the body.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxialRotor {
    /// Unit spin axis in the body's principal-axis frame; never normalized here.
    pub axis_body: [f64; 3],
    /// Positive axial spin inertia [kg m²].
    pub axial_inertia_kgm2: f64,
}

/// Mass/inertia ownership for a carrier and all attached rotors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotorAssembly<const N: usize> {
    /// Total assembly mass and REDUCED principal inertia: rotor transverse and
    /// offset mass inertia retained, rotor axial spin inertia excluded. Equivalently
    /// `I_reduced = I_locked - sum(J axis axis^T)`, in its principal frame.
    pub reduced_body: RigidBody,
    /// Independent rotor spin coordinates, in deterministic input order.
    pub rotors: [AxialRotor; N],
}

/// Complete axisymmetric body/rotor state, reusable for the next tick.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotorState<const N: usize> {
    /// Six-DOF carrier state in the native body-FRD/world-NED convention.
    pub rigid: SixDofState,
    /// Absolute axial angular momentum [kg m²/s], signed along each rotor axis.
    /// `p = J * (relative_spin + axis dot omega_body)`.
    pub axial_momentum: [f64; N],
    /// Unwrapped rotor angle relative to the carrier [rad].
    pub relative_phase_rad: [f64; N],
}

/// Torques held over one tick, signed along each rotor's declared axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotorTorques<const N: usize> {
    /// Internal motor torque on each rotor [N m]; opposite reaction is automatic.
    pub motor_nm: [f64; N],
    /// External aerodynamic torque on each rotor [N m]; drag usually opposes
    /// absolute rotor spin. It is included automatically in total external torque.
    pub aerodynamic_nm: [f64; N],
}

/// Accepted motion with separate internal-drive and external-work accounting.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotorStep<const N: usize> {
    /// Accepted complete endpoint state.
    pub state: RotorState<N>,
    /// Accepted relative spin rates [rad/s].
    pub relative_speed_rad_s: [f64; N],
    /// Net external linear impulse in world axes [N s].
    pub external_linear_impulse_world: [f64; 3],
    /// Net external angular impulse about the requested fixed world origin [N m s].
    pub external_angular_impulse_world: [f64; 3],
    /// Work supplied by the internal motor ports [J].
    pub motor_work_j: f64,
    /// External force, carrier-moment and rotor aerodynamic work [J].
    pub external_work_j: f64,
    /// Total translation, reduced-body rotation and rotor axial kinetic energy [J].
    pub kinetic_energy_before_j: f64,
    /// Actual accepted endpoint kinetic energy [J].
    pub kinetic_energy_after_j: f64,
    /// Actual kinetic change minus motor and external work [J]. Numerical
    /// splitting/quadrature error is reported, never called physical dissipation.
    pub kinetic_work_residual_j: f64,
    /// Actual global linear momentum change minus external impulse [N s].
    pub linear_momentum_residual: [f64; 3],
    /// Actual global spin-plus-orbital change minus external impulse [N m s].
    pub angular_momentum_residual: [f64; 3],
}

fn vec(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn array(a: Vec3) -> [f64; 3] {
    [a.x, a.y, a.z]
}
fn finite(a: Vec3) -> bool {
    a.x.is_finite() && a.y.is_finite() && a.z.is_finite()
}
fn reject(message: &str) -> Refusal {
    Refusal {
        code: "rotor-assembly-invalid",
        message: message.into(),
        ranked_repairs: vec![
            "supply finite admitted reduced inertia, unit axes, momenta and torques".into(),
        ],
    }
}

impl<const N: usize> RotorAssembly<N> {
    /// Validate fixed rotor count, axes, total mass and reduced inertia.
    ///
    /// # Errors
    /// Native rigid-body admission or `rotor-assembly-invalid`.
    pub fn admit(&self) -> Result<(), Refusal> {
        self.reduced_body.admit()?;
        if N > MAX_ROTORS {
            return Err(reject("rotor count exceeds 64"));
        }
        for rotor in &self.rotors {
            let axis = vec(rotor.axis_body);
            if !finite(axis)
                || (axis.dot(axis) - 1.0).abs() > 1e-12
                || !rotor.axial_inertia_kgm2.is_finite()
                || rotor.axial_inertia_kgm2 <= 0.0
            {
                return Err(reject(
                    "rotor axis must be unit and axial inertia positive finite",
                ));
            }
        }
        Ok(())
    }

    /// Build canonical absolute rotor momenta from measured relative spin rates.
    ///
    /// # Errors
    /// Model/state admission or non-finite rotor phase, speed or derived momentum.
    pub fn state_from_relative_spin(
        &self,
        rigid: SixDofState,
        phase_rad: [f64; N],
        speed_rad_s: [f64; N],
    ) -> Result<RotorState<N>, Refusal> {
        self.admit()?;
        let rigid = admit_and_canonicalize_state(&self.reduced_body, &rigid, 1.0 / TICK_HZ)?;
        if !phase_rad.iter().chain(&speed_rad_s).all(|x| x.is_finite()) {
            return Err(reject("rotor phase and relative speed must be finite"));
        }
        let axial_momentum = core::array::from_fn(|i| {
            self.rotors[i].axial_inertia_kgm2
                * (speed_rad_s[i] + vec(self.rotors[i].axis_body).dot(vec(rigid.omega_body)))
        });
        if !axial_momentum.iter().all(|x| x.is_finite()) {
            return Err(reject("rotor momentum overflow"));
        }
        Ok(RotorState {
            rigid,
            axial_momentum,
            relative_phase_rad: phase_rad,
        })
    }

    fn stored(&self, momenta: &[f64; N]) -> Vec3 {
        self.rotors
            .iter()
            .zip(momenta)
            .fold(Vec3::default(), |sum, (rotor, p)| {
                sum + vec(rotor.axis_body).scale(*p)
            })
    }

    fn body_momentum(&self, omega: [f64; 3]) -> Vec3 {
        vec(core::array::from_fn(|i| {
            self.reduced_body.inertia_kgm2[i] * omega[i]
        }))
    }

    fn energy(&self, state: &RotorState<N>) -> f64 {
        let velocity = vec(state.rigid.vel_mps);
        let omega = vec(state.rigid.omega_body);
        let mut energy = 0.5 * self.reduced_body.mass_kg * velocity.dot(velocity)
            + 0.5 * omega.dot(self.body_momentum(state.rigid.omega_body));
        for (rotor, momentum) in self.rotors.iter().zip(&state.axial_momentum) {
            energy += 0.5 * momentum * (momentum / rotor.axial_inertia_kgm2);
        }
        energy
    }
}

/// Advance an independent-rotor assembly under declared held loads.
///
/// `loads.force_n` is a net world force at the assembly COM; `loads.moment_nm`
/// is a body-frame carrier moment EXCLUDING the rotor aerodynamic torques
/// supplied in `torques`. Gravity belongs in the supplied force exactly once.
/// Translation uses drift/kick/drift; rotation uses the native driven gyrostat.
/// Every internal motor reaction and external rotor torque is counted once.
/// Phase uses the trapezoidal relative rate. The map is second order for these
/// held inputs, and successful steps use only fixed-size stack storage.
///
/// Geometry/mass properties are constant. Axes are body-fixed and spins are
/// independent; gear constraints, gimbals, imbalance, aero load evaluation,
/// moving masses and contact are not modeled. The caller must stage changing
/// load laws consistently and supply physical inertia rather than infer it
/// from an engine's equivalent shaft inertia. No historical calibration or
/// integration into the reduced longitudinal SimLoop is implied.
///
/// # Errors
/// Native spine/gyrostat admission and `rotor-assembly-invalid`. All inputs
/// are immutable and no physical state is published on any refusal.
pub fn step_held<const N: usize>(
    assembly: &RotorAssembly<N>,
    state: &RotorState<N>,
    loads: Loads,
    torques: RotorTorques<N>,
    angular_origin_world: [f64; 3],
    dt_s: f64,
) -> Result<RotorStep<N>, Refusal> {
    assembly.admit()?;
    let (canonical, rotation) =
        admitted_state_and_rotation(&assembly.reduced_body, &state.rigid, dt_s)?;
    let loads = admit_loads(loads, "rotor assembly held")?;
    let origin = vec(angular_origin_world);
    if !finite(origin)
        || !state
            .axial_momentum
            .iter()
            .chain(&state.relative_phase_rad)
            .chain(&torques.motor_nm)
            .chain(&torques.aerodynamic_nm)
            .all(|x| x.is_finite())
    {
        return Err(reject(
            "rotor state, torques and angular origin must be finite",
        ));
    }
    let mut next = *state;
    next.rigid = canonical;
    let stored_before = assembly.stored(&state.axial_momentum);
    let mut total_torque = vec(loads.moment_nm);
    for i in 0..N {
        next.axial_momentum[i] = dt_s.mul_add(
            torques.motor_nm[i] + torques.aerodynamic_nm[i],
            state.axial_momentum[i],
        );
        total_torque =
            total_torque + vec(assembly.rotors[i].axis_body).scale(torques.aerodynamic_nm[i]);
    }
    let stored_after = assembly.stored(&next.axial_momentum);
    if !next.axial_momentum.iter().all(|x| x.is_finite())
        || !finite(stored_before)
        || !finite(stored_after)
        || !finite(total_torque)
    {
        return Err(reject("rotor momentum or torque overflow"));
    }
    let (rotation_after, omega_after) = driven_step(
        rotation,
        vec(canonical.omega_body),
        vec(assembly.reduced_body.inertia_kgm2),
        stored_before,
        stored_after,
        HeldTorque::Body(total_torque),
        dt_s,
    )
    .map_err(|error| Refusal {
        code: "lie-step-refused",
        message: format!("native rotor/body step refused: {error}"),
        ranked_repairs: vec!["reduce the step or inspect rotor and carrier loads".into()],
    })?;
    next.rigid.quat = rotation_to_boundary(rotation_after);
    next.rigid.omega_body = array(omega_after);
    let force = vec(loads.force_n);
    let impulse = force.scale(dt_s);
    let velocity_before = vec(canonical.vel_mps);
    let velocity_after = velocity_before + impulse.scale(1.0 / assembly.reduced_body.mass_kg);
    let midpoint_position = vec(canonical.pos_m) + velocity_before.scale(dt_s * 0.5);
    next.rigid.vel_mps = array(velocity_after);
    next.rigid.pos_m = array(midpoint_position + velocity_after.scale(dt_s * 0.5));
    admit_derived_state(&next.rigid, "rotor assembly endpoint")?;
    let omega_mid = vec(canonical.omega_body).scale(0.5) + omega_after.scale(0.5);
    let mut speeds = [0.0; N];
    let mut motor_work = 0.0;
    let mut external_work = force
        .dot((velocity_before.scale(0.5) + velocity_after.scale(0.5)).scale(dt_s))
        + dt_s * vec(loads.moment_nm).dot(omega_mid);
    for i in 0..N {
        let rotor = assembly.rotors[i];
        let p_mid = 0.5 * state.axial_momentum[i] + 0.5 * next.axial_momentum[i];
        let absolute_mid_speed = p_mid / rotor.axial_inertia_kgm2;
        let phase_delta = dt_s * (absolute_mid_speed - vec(rotor.axis_body).dot(omega_mid));
        next.relative_phase_rad[i] += phase_delta;
        speeds[i] = next.axial_momentum[i] / rotor.axial_inertia_kgm2
            - vec(rotor.axis_body).dot(omega_after);
        motor_work += torques.motor_nm[i] * phase_delta;
        external_work += dt_s * torques.aerodynamic_nm[i] * absolute_mid_speed;
    }
    let p0 = velocity_before.scale(assembly.reduced_body.mass_kg);
    let p1 = velocity_after.scale(assembly.reduced_body.mass_kg);
    let l0 = rotation
        .rotate(assembly.body_momentum(canonical.omega_body) + stored_before)
        .map_err(|_| reject("initial angular momentum overflow"))?
        + (vec(canonical.pos_m) - origin).cross(p0);
    let l1 = rotation_after
        .rotate(assembly.body_momentum(next.rigid.omega_body) + stored_after)
        .map_err(|_| reject("endpoint angular momentum overflow"))?
        + (vec(next.rigid.pos_m) - origin).cross(p1);
    let angular_impulse = (midpoint_position - origin).cross(impulse)
        + (rotation
            .rotate(total_torque)
            .map_err(|_| reject("initial torque overflow"))?
            + rotation_after
                .rotate(total_torque)
                .map_err(|_| reject("endpoint torque overflow"))?)
        .scale(dt_s * 0.5);
    let kinetic_before = assembly.energy(state);
    let kinetic_after = assembly.energy(&next);
    let kinetic_residual = kinetic_after - kinetic_before - motor_work - external_work;
    let linear_residual = p1 - p0 - impulse;
    let angular_residual = l1 - l0 - angular_impulse;
    if !next
        .relative_phase_rad
        .iter()
        .chain(&speeds)
        .all(|x| x.is_finite())
        || ![
            kinetic_before,
            kinetic_after,
            kinetic_residual,
            motor_work,
            external_work,
        ]
        .iter()
        .all(|x| x.is_finite())
        || !finite(angular_impulse)
        || !finite(linear_residual)
        || !finite(angular_residual)
    {
        return Err(reject("rotor phase, work or momentum accounting overflow"));
    }
    Ok(RotorStep {
        state: next,
        relative_speed_rad_s: speeds,
        external_linear_impulse_world: array(impulse),
        external_angular_impulse_world: array(angular_impulse),
        motor_work_j: motor_work,
        external_work_j: external_work,
        kinetic_energy_before_j: kinetic_before,
        kinetic_energy_after_j: kinetic_after,
        kinetic_work_residual_j: kinetic_residual,
        linear_momentum_residual: array(linear_residual),
        angular_momentum_residual: array(angular_residual),
    })
}
