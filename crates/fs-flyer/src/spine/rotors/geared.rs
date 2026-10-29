//! One ideal, lossless gear coordinate coupled to full six-DOF carrier motion.
//!
//! With relative rotor rates `g_a s_dot`, the kinetic cross term is
//! `s_dot c.dot(omega)`, where `J_eq = sum(J_a g_a²)` and
//! `c = sum(J_a g_a e_a)`. The canonical shaft momentum is
//! `p = J_eq s_dot + c.dot(omega)`. Eliminating `s_dot` gives
//! `I_red = I_locked - c cᵀ/J_eq`, `h = c p/J_eq`, and total body
//! momentum `L = I_red omega + h`. In particular, counterrotation cancels
//! angular momentum according to SIGNED ratios, not reflected inertia.

use super::super::{
    Loads, SixDofState, TICK_HZ, admit_derived_state, admit_loads, admit_step, admit_time_interval,
    rotation_from_boundary, rotation_to_boundary,
};
use super::{MAX_ROTORS, array, finite, vec};
use crate::Refusal;
use fs_ga::{Mat3, So3, Vec3};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};

/// An axisymmetric rotor constrained to one common generalized shaft angle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GearedRotor {
    /// Unit spin axis in body coordinates. No implicit normalization.
    pub axis_body: [f64; 3],
    /// Physical axial spin inertia [kg m²], not a reflected engine inertia.
    pub axial_inertia_kgm2: f64,
    /// Signed, finite, nonzero ratio: relative rotor angle = ratio * s.
    pub ratio: f64,
}

/// Admitted, immutable mass properties for one ideal gear train.
#[derive(Clone, Debug)]
pub struct GearedAssembly<const N: usize> {
    mass_kg: f64,
    locked_inertia: Mat3,
    rotors: [GearedRotor; N],
    equivalent_inertia: f64,
    cross_inertia: Vec3,
    reduced_inertia: InertiaTensor3,
}

/// Carrier state and the single canonical gear coordinate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GearedState {
    /// Body-FRD/world-NED six-DOF state about the whole assembly COM.
    pub rigid: SixDofState,
    /// Generalized shaft momentum `p = J_eq s_dot + c.dot(omega)` [kg m²/s].
    pub momentum: f64,
    /// Unwrapped common relative gear angle `s` [rad].
    pub phase_rad: f64,
}

/// One evaluated external wrench and internal shaft drive.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GearedLoads<const N: usize> {
    /// World force at the assembly COM and body moment about it. EXCLUDES
    /// axial rotor aerodynamic torques below; includes gravity if desired.
    pub carrier: Loads,
    /// Generalized INTERNAL motor torque [N m], conjugate to `s`.
    pub motor_nm: f64,
    /// EXTERNAL axial aerodynamic torque on each rotor [N m], signed along
    /// its declared axis. The owner adds `sum(e_a tau_a)` exactly once.
    pub aerodynamic_nm: [f64; N],
}

/// Accepted gear/body state with work and global momentum accounting.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GearedStep<const N: usize> {
    /// Accepted endpoint, reusable without reconstructing shaft momentum.
    pub state: GearedState,
    /// Common relative shaft rate `s_dot` [rad/s].
    pub speed_rad_s: f64,
    /// Actual relative rotor rates `g_a s_dot` [rad/s].
    pub relative_speed_rad_s: [f64; N],
    /// Absolute axial rotor rates `g_a s_dot + e_a.dot(omega)` [rad/s].
    pub absolute_speed_rad_s: [f64; N],
    /// External linear impulse in world axes [N s].
    pub external_linear_impulse_world: [f64; 3],
    /// External angular impulse about the requested fixed world origin [N m s].
    pub external_angular_impulse_world: [f64; 3],
    /// Internal motor work `Q_motor * delta_s` [J].
    pub motor_work_j: f64,
    /// Work of world force, carrier moment and aerodynamic shaft ports [J].
    pub external_work_j: f64,
    /// Physical total kinetic energy at the input [J].
    pub kinetic_energy_before_j: f64,
    /// Physical total kinetic energy at the accepted endpoint [J].
    pub kinetic_energy_after_j: f64,
    /// Kinetic change minus declared work [J]; numerical error, not damping.
    pub kinetic_work_residual_j: f64,
    /// Actual world linear momentum change minus external impulse [N s].
    pub linear_momentum_residual: [f64; 3],
    /// Actual world spin-plus-orbital change minus external impulse [N m s].
    pub angular_momentum_residual: [f64; 3],
}

fn reject(message: &str) -> Refusal {
    Refusal {
        code: "geared-rotor-invalid",
        message: message.into(),
        ranked_repairs: vec![
            "supply physical locked inertia, unit axes, signed gear ratios and finite state".into(),
        ],
    }
}

fn map_lie(error: fs_time::lie::LieStepError) -> Refusal {
    Refusal {
        code: "lie-step-refused",
        message: format!("native ideal-gear/body step refused: {error}"),
        ranked_repairs: vec!["inspect reduced inertia and loads or reduce the step".into()],
    }
}

impl<const N: usize> GearedAssembly<N> {
    /// Assemble the full-inertia Schur reduction once. Locked inertia includes
    /// every rotor's axial/transverse and offset inertia about the assembly COM.
    /// Both locked and reduced tensors must be symmetric positive definite.
    ///
    /// # Errors
    /// `geared-rotor-invalid` for count, mass, rotor or derived-value admission;
    /// `lie-step-refused` for a refused physical inertia tensor.
    pub fn new(
        mass_kg: f64,
        locked_inertia_kgm2: Mat3,
        rotors: [GearedRotor; N],
    ) -> Result<Self, Refusal> {
        if N == 0 || N > MAX_ROTORS || !mass_kg.is_finite() || mass_kg <= 0.0 {
            return Err(reject(
                "require positive mass and between 1 and 64 geared rotors",
            ));
        }
        InertiaTensor3::new(locked_inertia_kgm2).map_err(map_lie)?;
        let mut equivalent_inertia = 0.0;
        let mut cross_inertia = Vec3::default();
        for rotor in &rotors {
            let axis = vec(rotor.axis_body);
            if !finite(axis)
                || (axis.dot(axis) - 1.0).abs() > 1e-12
                || !rotor.axial_inertia_kgm2.is_finite()
                || rotor.axial_inertia_kgm2 <= 0.0
                || !rotor.ratio.is_finite()
                || rotor.ratio == 0.0
            {
                return Err(reject(
                    "rotors require unit axes, positive physical axial inertia and nonzero signed ratios",
                ));
            }
            let reflected = rotor.axial_inertia_kgm2 * rotor.ratio;
            equivalent_inertia += reflected * rotor.ratio;
            cross_inertia = cross_inertia + axis.scale(reflected);
        }
        if !equivalent_inertia.is_finite() || equivalent_inertia <= 0.0 || !finite(cross_inertia) {
            return Err(reject("gear inertia reduction overflow or underflow"));
        }
        let c = array(cross_inertia);
        let mut reduced = locked_inertia_kgm2;
        for i in 0..3 {
            for j in 0..3 {
                reduced.m[3 * i + j] -= c[i] * (c[j] / equivalent_inertia);
            }
        }
        let reduced_inertia = InertiaTensor3::new(reduced).map_err(map_lie)?;
        Ok(Self {
            mass_kg,
            locked_inertia: locked_inertia_kgm2,
            rotors,
            equivalent_inertia,
            cross_inertia,
            reduced_inertia,
        })
    }

    /// Total physical mass [kg].
    #[must_use]
    pub fn mass_kg(&self) -> f64 {
        self.mass_kg
    }

    /// Caller-supplied physical locked inertia [kg m²].
    #[must_use]
    pub fn locked_inertia_kgm2(&self) -> Mat3 {
        self.locked_inertia
    }

    /// Fixed rotor definitions, in their admitted order.
    #[must_use]
    pub fn rotors(&self) -> &[GearedRotor; N] {
        &self.rotors
    }

    /// Reflected axial inertia `sum(J_a g_a²)` [kg m²].
    #[must_use]
    pub fn equivalent_inertia_kgm2(&self) -> f64 {
        self.equivalent_inertia
    }

    /// Signed kinetic cross coefficient `sum(J_a g_a e_a)` [kg m²].
    #[must_use]
    pub fn cross_inertia_kgm2(&self) -> [f64; 3] {
        array(self.cross_inertia)
    }

    /// Schur-reduced tensor; products of inertia are retained.
    #[must_use]
    pub fn reduced_inertia_kgm2(&self) -> Mat3 {
        self.reduced_inertia.matrix()
    }

    /// Construct canonical shaft momentum from a measured relative gear rate.
    ///
    /// # Errors
    /// Native spine admission or `geared-rotor-invalid` on nonfinite gear data.
    pub fn state_from_speed(
        &self,
        rigid: SixDofState,
        phase_rad: f64,
        speed_rad_s: f64,
    ) -> Result<GearedState, Refusal> {
        let mut state = GearedState {
            rigid,
            phase_rad,
            momentum: self.equivalent_inertia * speed_rad_s
                + self.cross_inertia.dot(vec(rigid.omega_body)),
        };
        let rotation = self.admit_state(&state, 1.0 / TICK_HZ)?;
        state.rigid.quat = rotation_to_boundary(rotation);
        Ok(state)
    }

    /// Recover the common relative gear rate [rad/s].
    ///
    /// # Errors
    /// `geared-rotor-invalid` if recovery is nonfinite.
    pub fn speed(&self, state: &GearedState) -> Result<f64, Refusal> {
        let speed = (state.momentum - self.cross_inertia.dot(vec(state.rigid.omega_body)))
            / self.equivalent_inertia;
        if speed.is_finite() {
            Ok(speed)
        } else {
            Err(reject("nonfinite generalized shaft speed"))
        }
    }

    fn stored(&self, momentum: f64) -> Vec3 {
        self.cross_inertia.scale(momentum / self.equivalent_inertia)
    }

    pub(super) fn admit_state(&self, state: &GearedState, dt_s: f64) -> Result<So3, Refusal> {
        admit_step(&state.rigid, dt_s)?;
        if !state.momentum.is_finite() || !state.phase_rad.is_finite() {
            return Err(reject("gear momentum and phase must be finite"));
        }
        self.speed(state)?;
        rotation_from_boundary(state.rigid.quat)
    }

    fn energy(&self, state: &GearedState) -> Result<f64, Refusal> {
        let v = vec(state.rigid.vel_mps);
        let w = vec(state.rigid.omega_body);
        Ok(0.5 * self.mass_kg * v.dot(v)
            + 0.5 * w.dot(self.reduced_inertia.momentum(w).map_err(map_lie)?)
            + 0.5 * state.momentum * (state.momentum / self.equivalent_inertia))
    }
}

/// Advance under one fixed world force and fixed body carrier/aerodynamic
/// torques. Internal motor torque changes `p` but never total external torque.
/// Aerodynamic torque changes `p` by `sum(g_a tau_a)` and total body torque
/// by `sum(e_a tau_a)` exactly once. Its shaft power uses ABSOLUTE rotor rate;
/// motor work uses the common RELATIVE phase change.
///
/// Translation is drift/kick/drift; native full-tensor gyrostat rotation is a
/// symmetric split. Phase/work use trapezoidal rates. This is second order
/// for the declared held loads. Use [`step`] for state-dependent load laws.
///
/// # Errors
/// Native state/step/group admission, load admission or gear/gyrostat refusal.
/// Inputs are immutable; a refusal publishes no physical endpoint.
pub fn step_held<const N: usize>(
    assembly: &GearedAssembly<N>,
    state: &GearedState,
    loads: GearedLoads<N>,
    angular_origin_world: [f64; 3],
    dt_s: f64,
) -> Result<GearedStep<N>, Refusal> {
    let rotation = assembly.admit_state(state, dt_s)?;
    advance_held(assembly, state, rotation, loads, angular_origin_world, dt_s)
}

/// Evaluate a smooth state-dependent load law at the start and a half-step
/// predictor, then advance with that midpoint load. Both evaluations receive
/// canonical full six-DOF state and the actual gear momentum. This explicit
/// midpoint staging retains order two for smooth admitted load laws. The load
/// function is called exactly twice in start/midpoint order; its own side
/// effects are the caller's responsibility if a later stage refuses.
///
/// # Errors
/// As [`step_held`], plus propagated load-law refusals and invalid clock data.
pub fn step<const N: usize, F>(
    assembly: &GearedAssembly<N>,
    state: &GearedState,
    t_s: f64,
    dt_s: f64,
    angular_origin_world: [f64; 3],
    mut loads: F,
) -> Result<GearedStep<N>, Refusal>
where
    F: FnMut(f64, &GearedState) -> Result<GearedLoads<N>, Refusal>,
{
    let rotation = assembly.admit_state(state, dt_s)?;
    admit_time_interval(t_s, dt_s)?;
    if !finite(vec(angular_origin_world)) {
        return Err(reject("angular momentum origin must be finite"));
    }
    let mut canonical = *state;
    canonical.rigid.quat = rotation_to_boundary(rotation);
    let initial = loads(t_s, &canonical)?;
    let predictor = advance_held(
        assembly,
        &canonical,
        rotation,
        initial,
        angular_origin_world,
        0.5 * dt_s,
    )?;
    let midpoint = loads(t_s + 0.5 * dt_s, &predictor.state)?;
    advance_held(
        assembly,
        &canonical,
        rotation,
        midpoint,
        angular_origin_world,
        dt_s,
    )
}

fn advance_held<const N: usize>(
    assembly: &GearedAssembly<N>,
    state: &GearedState,
    rotation: So3,
    loads: GearedLoads<N>,
    angular_origin_world: [f64; 3],
    dt_s: f64,
) -> Result<GearedStep<N>, Refusal> {
    admit_loads(loads.carrier, "geared rotor")?;
    let origin = vec(angular_origin_world);
    if !finite(origin)
        || !loads.motor_nm.is_finite()
        || !loads.aerodynamic_nm.iter().all(|x| x.is_finite())
    {
        return Err(reject("origin and shaft torques must be finite"));
    }
    let mut generalized_torque = loads.motor_nm;
    let mut total_torque = vec(loads.carrier.moment_nm);
    for (rotor, torque) in assembly.rotors.iter().zip(loads.aerodynamic_nm) {
        generalized_torque += rotor.ratio * torque;
        total_torque = total_torque + vec(rotor.axis_body).scale(torque);
    }
    let mut next = *state;
    next.momentum += dt_s * generalized_torque;
    let stored_before = assembly.stored(state.momentum);
    let stored_after = assembly.stored(next.momentum);
    if !next.momentum.is_finite()
        || !finite(total_torque)
        || !finite(stored_before)
        || !finite(stored_after)
    {
        return Err(reject("shaft or total momentum overflow"));
    }
    let omega0 = vec(state.rigid.omega_body);
    let (rotation_after, omega1) = assembly
        .reduced_inertia
        .driven_step(
            rotation,
            omega0,
            stored_before,
            stored_after,
            HeldTorque::Body(total_torque),
            dt_s,
        )
        .map_err(map_lie)?;
    next.rigid.quat = rotation_to_boundary(rotation_after);
    next.rigid.omega_body = array(omega1);
    let force = vec(loads.carrier.force_n);
    let impulse = force.scale(dt_s);
    let v0 = vec(state.rigid.vel_mps);
    let v1 = v0 + impulse.scale(1.0 / assembly.mass_kg);
    let kick_position = vec(state.rigid.pos_m) + v0.scale(0.5 * dt_s);
    next.rigid.pos_m = array(kick_position + v1.scale(0.5 * dt_s));
    next.rigid.vel_mps = array(v1);
    admit_derived_state(&next.rigid, "geared rotor endpoint")?;
    let speed0 = assembly.speed(state)?;
    let speed1 = assembly.speed(&next)?;
    let phase_delta = dt_s * (0.5 * speed0 + 0.5 * speed1);
    next.phase_rad += phase_delta;
    let omega_mid = omega0.scale(0.5) + omega1.scale(0.5);
    let motor_work = loads.motor_nm * phase_delta;
    let mut external_work = force.dot((v0.scale(0.5) + v1.scale(0.5)).scale(dt_s))
        + dt_s * vec(loads.carrier.moment_nm).dot(omega_mid);
    let mut relative_speed = [0.0; N];
    let mut absolute_speed = [0.0; N];
    for i in 0..N {
        let rotor = assembly.rotors[i];
        relative_speed[i] = rotor.ratio * speed1;
        absolute_speed[i] = relative_speed[i] + vec(rotor.axis_body).dot(omega1);
        external_work += loads.aerodynamic_nm[i]
            * (rotor.ratio * phase_delta + dt_s * vec(rotor.axis_body).dot(omega_mid));
    }
    let p0 = v0.scale(assembly.mass_kg);
    let p1 = v1.scale(assembly.mass_kg);
    let l0 = rotation
        .rotate(assembly.reduced_inertia.momentum(omega0).map_err(map_lie)? + stored_before)
        .map_err(|_| reject("initial angular momentum overflow"))?
        + (vec(state.rigid.pos_m) - origin).cross(p0);
    let l1 = rotation_after
        .rotate(assembly.reduced_inertia.momentum(omega1).map_err(map_lie)? + stored_after)
        .map_err(|_| reject("endpoint angular momentum overflow"))?
        + (vec(next.rigid.pos_m) - origin).cross(p1);
    let angular_impulse = (kick_position - origin).cross(impulse)
        + (rotation
            .rotate(total_torque)
            .map_err(|_| reject("initial torque overflow"))?
            + rotation_after
                .rotate(total_torque)
                .map_err(|_| reject("endpoint torque overflow"))?)
        .scale(0.5 * dt_s);
    let kinetic_before = assembly.energy(state)?;
    let kinetic_after = assembly.energy(&next)?;
    let kinetic_residual = kinetic_after - kinetic_before - motor_work - external_work;
    let linear_residual = p1 - p0 - impulse;
    let angular_residual = l1 - l0 - angular_impulse;
    if ![
        next.phase_rad,
        motor_work,
        external_work,
        kinetic_before,
        kinetic_after,
        kinetic_residual,
    ]
    .iter()
    .chain(&relative_speed)
    .chain(&absolute_speed)
    .all(|x| x.is_finite())
        || !finite(angular_impulse)
        || !finite(linear_residual)
        || !finite(angular_residual)
    {
        return Err(reject("gear phase, work or momentum accounting overflow"));
    }
    Ok(GearedStep {
        state: next,
        speed_rad_s: speed1,
        relative_speed_rad_s: relative_speed,
        absolute_speed_rad_s: absolute_speed,
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
