//! Finite shaft actuation on the native full-tensor gyrostat.
//!
//! Fixed body axes a_i and absolute axial momenta p_i give
//!
//!     L = I omega + sum_i a_i p_i
//!     T = omega.I.omega/2 + sum_i p_i^2/(2 J_i)
//!     relative_rate_i = p_i/J_i - a_i.omega
//!
//! I retains each rotor's transverse inertia and excludes its axial inertia.
//! A motor impulse changes p_i, with equal/opposite carrier recoil at fixed L.
//! The existing native gyrostat map owns all attitude/precession drift.
//! Rotors are balanced and axisymmetric with fixed axes and fixed inertia;
//! eccentric masses, moving gimbals, bearings and electrical states are not
//! inferred. Net shaft torque may drive or brake; no artificial loss is added.

use super::{HeldTorque, InertiaTensor3, LieStepError};
use fs_ga::{Mat3, So3, Vec3};
use fs_math::dd::Dd;

#[path = "rotors/viscous.rs"]
mod viscous;
pub use viscous::{ViscousRotorStep, ViscousRotorStepper};

#[path = "rotors/torsional.rs"]
mod torsional;
pub use torsional::{TorsionalRotorStep, TorsionalRotorStepper};

#[path = "rotors/locked.rs"]
mod locked;
pub use locked::{LockedRotorDrift, LockedRotorImpulse, RotorLock, RotorLockEngagement};

/// Maximum fixed-axis rotors in this allocation-free stepping lane.
pub const MAX_AXIAL_ROTORS: usize = 32;

/// Admitted spin axis and polar inertia of one balanced axisymmetric rotor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxialRotor {
    axis_body: Vec3,
    axial_inertia_kg_m2: f64,
}

impl AxialRotor {
    /// Normalize a finite nonzero body axis once and admit positive polar inertia.
    pub fn new(axis_body: Vec3, axial_inertia_kg_m2: f64) -> Result<Self, LieStepError> {
        vector(axis_body, "rotor axis")?;
        scalar(axial_inertia_kg_m2, "rotor axial inertia")?;
        let scale = axis_body
            .x
            .abs()
            .max(axis_body.y.abs())
            .max(axis_body.z.abs());
        if scale == 0.0 || axial_inertia_kg_m2 <= 0.0 {
            return Err(invalid("nonzero rotor axis and positive axial inertia"));
        }
        let scaled = axis_body.scale(1.0 / scale);
        // Divide components separately when 1/scale would overflow.
        let scaled = if finite(scaled) {
            scaled
        } else {
            Vec3::new(
                axis_body.x / scale,
                axis_body.y / scale,
                axis_body.z / scale,
            )
        };
        let axis_body = scaled.scale(1.0 / scaled.norm());
        vector(axis_body, "normalized rotor axis")?;
        Ok(Self {
            axis_body,
            axial_inertia_kg_m2,
        })
    }

    /// Unit spin direction in the original fixed body frame.
    pub const fn axis_body(self) -> Vec3 {
        self.axis_body
    }
    /// Polar inertia, kg m^2; it is not included in the reduced carrier tensor.
    pub const fn axial_inertia_kg_m2(self) -> f64 {
        self.axial_inertia_kg_m2
    }
}

/// Canonical rotational state, with unwrapped relative rotor angles.
///
/// Total body momentum includes all rotors. Absolute axial momentum is
/// J*(relative rotor rate + axis.dot(carrier omega)), not J*relative rate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotorState<const N: usize> {
    /// Carrier attitude, mapping body coordinates into the inertial frame.
    pub rotation: So3,
    /// Total carrier-plus-rotor angular momentum in body axes, N m s.
    pub total_momentum_body: Vec3,
    /// Absolute axial momentum of each rotor, N m s, in inventory order.
    pub axial_momentum_nms: [f64; N],
    /// Unwrapped shaft angles relative to the carrier, rad, in the same order.
    pub relative_angle_rad: [f64; N],
}

/// A fixed rotor inventory sharing one admitted full reduced inertia.
///
/// Preparation caches the existing principal-frame representation. A step is
/// O(N), allocation-free, and uses one native gyrostat drift. An empty inventory
/// reduces to the same rigid-body map; at most 32 rotors are admitted.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotorBank<const N: usize> {
    inertia: InertiaTensor3,
    rotors: [AxialRotor; N],
}

/// Actual simultaneous internal shaft and external angular impulse.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotorImpulse<const N: usize> {
    /// Accepted canonical state; attitude and relative angles are unchanged.
    pub state_after: RotorState<N>,
    /// Actual rounded changes of rotor axial momentum, N m s.
    pub axial_momentum_change_nms: [f64; N],
    /// Actual total-body-momentum change, N m s, excluding internal transfers.
    pub external_angular_impulse_body: Vec3,
    /// Mean relative shaft rates at the two impulse endpoints, rad/s.
    pub mean_relative_rate_rad_s: [f64; N],
    /// Individual motor/brake work: actual Delta p times mean relative rate, J.
    pub shaft_work_j: [f64; N],
    /// External torque-impulse work at mean carrier angular velocity, J.
    pub external_work_j: f64,
    /// Full rotor-plus-carrier kinetic change, using polarized increments, J.
    pub kinetic_energy_change_j: f64,
    /// Kinetic change minus external and shaft work; measured rounding defect.
    pub work_energy_residual_j: f64,
}

/// Accepted torque-driven attitude, spin and relative-phase trajectory.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotorStep<const N: usize> {
    /// Complete state for the next step.
    pub state_after: RotorState<N>,
    /// Actual rotor axial-momentum changes, N m s.
    pub axial_momentum_change_nms: [f64; N],
    /// Actual external angular impulse in inertial coordinates, N m s.
    pub external_angular_impulse_world: Vec3,
    /// Individual shaft work accumulated from the two finite transfer events, J.
    pub shaft_work_j: [f64; N],
    /// External torque work accumulated at its two endpoint-frame kicks, J.
    pub external_work_j: f64,
    /// Full kinetic change from the original to returned state, J.
    pub kinetic_energy_change_j: f64,
    /// Native free-gyrostat drift's measured kinetic defect, J; never motor loss.
    pub drift_energy_defect_j: f64,
    /// Kinetic change minus both work sources AND the reported drift defect, J.
    pub work_energy_residual_j: f64,
    /// Actual spatial angular-momentum change minus external impulse, N m s.
    pub angular_momentum_residual_world: Vec3,
}

impl<const N: usize> RotorBank<N> {
    /// Bind the native reduced tensor, already excluding rotor axial inertias.
    /// No geometry or physical mass distribution is inferred from this form.
    pub fn from_reduced(
        inertia: InertiaTensor3,
        rotors: [AxialRotor; N],
    ) -> Result<Self, LieStepError> {
        if N > MAX_AXIAL_ROTORS {
            return Err(invalid("rotor inventory exceeds 32"));
        }
        Ok(Self { inertia, rotors })
    }

    /// Reduce the complete LOCKED assembly tensor by sum_i J_i a_i a_i^T.
    ///
    /// The supplied tensor must already include every rotor's complete locked
    /// inertia, including transverse inertia and any fixed carrier mass.
    /// Native tensor admission owns symmetry, positivity and basis preparation.
    /// A nonpositive reduced form refuses; no axial inertia is counted twice.
    pub fn from_locked(
        locked_inertia_body: Mat3,
        rotors: [AxialRotor; N],
    ) -> Result<Self, LieStepError> {
        if N > MAX_AXIAL_ROTORS {
            return Err(invalid("rotor inventory exceeds 32"));
        }
        let locked = InertiaTensor3::new(locked_inertia_body)?.matrix();
        let mut reduced = [0.0; 9];
        for row in 0..3 {
            for col in row..3 {
                let mut value = dd(locked.m[3 * row + col]);
                for rotor in rotors {
                    let a = components(rotor.axis_body);
                    value = value - dd(rotor.axial_inertia_kg_m2) * dd(a[row]) * dd(a[col]);
                }
                reduced[3 * row + col] = value.to_f64();
                reduced[3 * col + row] = value.to_f64();
            }
        }
        Self::from_reduced(InertiaTensor3::new(Mat3 { m: reduced })?, rotors)
    }

    /// Native cached reduced inertia in the original body frame.
    pub const fn reduced_inertia(&self) -> InertiaTensor3 {
        self.inertia
    }
    /// Original admitted rotor inventory.
    pub const fn rotors(&self) -> &[AxialRotor; N] {
        &self.rotors
    }

    /// Initialize from carrier angular velocity and RELATIVE shaft speeds.
    ///
    /// Refuses when total momentum cannot resolve the requested carrier rate
    /// against stored rotor momentum. Canonical states may also be supplied
    /// directly; subsequent methods always interpret their actual stored values.
    pub fn state_from_velocities(
        &self,
        rotation: So3,
        omega_body: Vec3,
        relative_rate_rad_s: [f64; N],
        relative_angle_rad: [f64; N],
    ) -> Result<RotorState<N>, LieStepError> {
        vector(omega_body, "carrier angular velocity")?;
        let mut axial = [0.0; N];
        for i in 0..N {
            scalar(relative_rate_rad_s[i], "relative shaft rate")?;
            axial[i] = scalar(
                self.rotors[i].axial_inertia_kg_m2
                    * (relative_rate_rad_s[i] + self.rotors[i].axis_body.dot(omega_body)),
                "absolute rotor axial momentum",
            )?;
        }
        let state = RotorState {
            rotation,
            total_momentum_body: self.inertia.momentum(omega_body)?
                + self.rotor_momentum(&axial)?,
            axial_momentum_nms: axial,
            relative_angle_rad,
        };
        self.admit_state(&state)?;
        let recovered = self.angular_velocity(&state)?;
        let scale = components(omega_body)
            .into_iter()
            .map(f64::abs)
            .fold(0.0, f64::max);
        let error = components(recovered - omega_body)
            .into_iter()
            .map(f64::abs)
            .fold(0.0, f64::max);
        if error > 4096.0 * f64::EPSILON * scale {
            return Err(invalid(
                "unresolved carrier velocity against rotor momentum",
            ));
        }
        self.kinetic_energy(&state)?;
        Ok(state)
    }

    /// Recover carrier angular velocity from the actual canonical momenta.
    pub fn angular_velocity(&self, state: &RotorState<N>) -> Result<Vec3, LieStepError> {
        self.admit_state(state)?;
        self.inertia.velocity(self.carrier_momentum(state)?)
    }

    /// Recover relative shaft speeds, including the carrier's axial motion.
    pub fn relative_rates(&self, state: &RotorState<N>) -> Result<[f64; N], LieStepError> {
        let omega = self.angular_velocity(state)?;
        let mut rates = [0.0; N];
        for i in 0..N {
            rates[i] = scalar(
                state.axial_momentum_nms[i] / self.rotors[i].axial_inertia_kg_m2
                    - self.rotors[i].axis_body.dot(omega),
                "relative shaft rate",
            )?;
        }
        Ok(rates)
    }

    /// Complete kinetic energy, including absolute rotor axial energy.
    pub fn kinetic_energy(&self, state: &RotorState<N>) -> Result<f64, LieStepError> {
        self.admit_state(state)?;
        let carrier = self.carrier_momentum(state)?;
        let omega = self.inertia.velocity(carrier)?;
        let mut value = dot_dd(omega, carrier) * dd(0.5);
        for i in 0..N {
            value = value
                + dd(0.5 * state.axial_momentum_nms[i] / self.rotors[i].axial_inertia_kg_m2)
                    * dd(state.axial_momentum_nms[i]);
        }
        scalar(value.to_f64(), "rotor-bank kinetic energy")
    }

    /// Apply simultaneous finite shaft and external body angular impulses.
    ///
    /// Internal impulses change absolute rotor p while leaving TOTAL L fixed.
    /// All motors share the same original and returned carrier velocities:
    /// cross-axis recoil is not approximated by independent effective inertias.
    /// For fixed R,L the relative-speed mobility is diag(1/J)+A^T I^-1 A.
    /// Work uses actual rounded increments and the original kinetic form.
    /// Nonfinite or unresolved nonzero increments refuse atomically. The native
    /// aggregate rotor vector must also resolve the simultaneous carrier recoil:
    /// its increment is checked against actual shaft/external transfers using
    /// the carrier/transfer scale, never the potentially much larger stored spin.
    pub fn apply_impulse(
        &self,
        state: &RotorState<N>,
        shaft_impulse_nms: [f64; N],
        external_impulse_body: Vec3,
    ) -> Result<RotorImpulse<N>, LieStepError> {
        self.admit_state(state)?;
        let before_carrier = self.carrier_momentum(state)?;
        let before_omega = self.inertia.velocity(before_carrier)?;
        vector(external_impulse_body, "external angular impulse")?;
        let mut after = *state;
        after.total_momentum_body =
            super::super::torque_kick(state.total_momentum_body, external_impulse_body, 1.0)?;
        for ((old, new), impulse) in components(state.total_momentum_body)
            .into_iter()
            .zip(components(after.total_momentum_body))
            .zip(components(external_impulse_body))
        {
            if impulse != 0.0 && old == new {
                return Err(invalid("unresolved external angular impulse"));
            }
        }
        let mut actual = [0.0; N];
        for i in 0..N {
            after.axial_momentum_nms[i] = increment(
                state.axial_momentum_nms[i],
                shaft_impulse_nms[i],
                "rotor axial impulse",
            )?;
            actual[i] = after.axial_momentum_nms[i] - state.axial_momentum_nms[i];
        }
        let after_carrier = self.carrier_momentum(&after)?;
        let external = after.total_momentum_body - state.total_momentum_body;
        self.admit_recoil(before_carrier, after_carrier, external, &actual)?;
        let after_omega = self.inertia.velocity(after_carrier)?;
        let mean_omega = before_omega.scale(0.5) + after_omega.scale(0.5);
        let mut rates = [0.0; N];
        let mut shaft_work = [0.0; N];
        for i in 0..N {
            rates[i] = scalar(
                (0.5 * state.axial_momentum_nms[i] + 0.5 * after.axial_momentum_nms[i])
                    / self.rotors[i].axial_inertia_kg_m2
                    - self.rotors[i].axis_body.dot(mean_omega),
                "mean relative shaft rate",
            )?;
            shaft_work[i] = scalar(actual[i] * rates[i], "shaft impulse work")?;
        }
        let external_work = scalar(
            dot_dd(external, mean_omega).to_f64(),
            "external impulse work",
        )?;
        let kinetic_change = self.energy_change(state, &after)?;
        let residual = work_residual(kinetic_change, &shaft_work, external_work, 0.0)?;
        Ok(RotorImpulse {
            state_after: after,
            axial_momentum_change_nms: actual,
            external_angular_impulse_body: external,
            mean_relative_rate_rad_s: rates,
            shaft_work_j: shaft_work,
            external_work_j: external_work,
            kinetic_energy_change_j: kinetic_change,
            work_energy_residual_j: residual,
        })
    }

    /// Advance held motor torques and a body- or space-held external torque.
    ///
    /// Two simultaneous half-kicks surround ONE existing native gyrostat drift.
    /// Relative phases use the same kick-endpoint relative-speed quadrature as
    /// shaft work. With representable held torques, torque times this angle
    /// increment equals the reported shaft work to rounding precision.
    ///
    /// This is a symmetric second-order approximation with fixed axes/inertia.
    /// The native drift's measured energy defect remains separate from motor
    /// work; it is never physical damping. Resolve gyroscopic motion and phase
    /// by refinement. State-dependent motor laws need an appropriate coupled
    /// stage; this selector does not freeze or solve an electrical/controller
    /// state, impose speed limits, or silently clamp a torque.
    ///
    /// A finite nonzero signed duration supports reversibility checks. All work
    /// is bounded by 32 rotors and the native constant-work drift; there are no
    /// callbacks, allocations or iterative solves. Inputs are unchanged on error.
    pub fn step(
        &self,
        state: &RotorState<N>,
        shaft_torque_nm: [f64; N],
        external_torque: HeldTorque,
        duration_s: f64,
    ) -> Result<RotorStep<N>, LieStepError> {
        scalar(duration_s, "rotor-bank duration")?;
        if duration_s == 0.0 || 0.5 * duration_s == 0.0 {
            return Err(invalid("nonzero resolved rotor-bank duration"));
        }
        let half = 0.5 * duration_s;
        let mut shaft_impulse = [0.0; N];
        for i in 0..N {
            shaft_impulse[i] = product(shaft_torque_nm[i], half, "shaft torque impulse")?;
        }
        let torque_at = |rotation: So3| -> Result<Vec3, LieStepError> {
            match external_torque {
                HeldTorque::Body(value) => vector(value, "external body torque"),
                HeldTorque::Space(value) => Ok(rotation.inverse().rotate(value)?),
            }
        };
        let impulse_at = |rotation: So3| -> Result<Vec3, LieStepError> {
            let torque = torque_at(rotation)?;
            Ok(Vec3::new(
                product(torque.x, half, "external torque impulse")?,
                product(torque.y, half, "external torque impulse")?,
                product(torque.z, half, "external torque impulse")?,
            ))
        };
        let first = self.apply_impulse(state, shaft_impulse, impulse_at(state.rotation)?)?;
        let mut drifted = first.state_after;
        (drifted.rotation, drifted.total_momentum_body) = self.inertia.momentum_step(
            drifted.rotation,
            drifted.total_momentum_body,
            self.rotor_momentum(&drifted.axial_momentum_nms)?,
            duration_s,
        )?;
        let second = self.apply_impulse(&drifted, shaft_impulse, impulse_at(drifted.rotation)?)?;
        self.finish_step(state, &first, &drifted, &second, duration_s, None)
    }

    fn finish_step(
        &self,
        state: &RotorState<N>,
        first: &RotorImpulse<N>,
        drifted: &RotorState<N>,
        second: &RotorImpulse<N>,
        duration_s: f64,
        accepted_phase: Option<[f64; N]>,
    ) -> Result<RotorStep<N>, LieStepError> {
        let drift_defect = self.energy_change(&first.state_after, drifted)?;
        let mut after = second.state_after;
        let mut shaft_work = [0.0; N];
        let mut axial_change = [0.0; N];
        for i in 0..N {
            after.relative_angle_rad[i] = if let Some(phase) = accepted_phase.as_ref() {
                scalar(phase[i], "accepted relative rotor phase")?
            } else {
                let rate = 0.5 * first.mean_relative_rate_rad_s[i]
                    + 0.5 * second.mean_relative_rate_rad_s[i];
                let angle = product(duration_s, rate, "relative phase increment")?;
                increment(state.relative_angle_rad[i], angle, "relative rotor phase")?
            };
            shaft_work[i] = scalar(
                (dd(first.shaft_work_j[i]) + dd(second.shaft_work_j[i])).to_f64(),
                "shaft step work",
            )?;
            axial_change[i] = after.axial_momentum_nms[i] - state.axial_momentum_nms[i];
        }
        let external_work = scalar(
            (dd(first.external_work_j) + dd(second.external_work_j)).to_f64(),
            "external step work",
        )?;
        let external_impulse = state.rotation.rotate(first.external_angular_impulse_body)?
            + drifted
                .rotation
                .rotate(second.external_angular_impulse_body)?;
        vector(external_impulse, "external spatial angular impulse")?;
        let before_world = state.rotation.rotate(state.total_momentum_body)?;
        let after_world = after.rotation.rotate(after.total_momentum_body)?;
        let residual_world = subtract3(after_world, before_world, external_impulse);
        vector(residual_world, "rotor spatial momentum balance")?;
        let kinetic_change = self.energy_change(state, &after)?;
        Ok(RotorStep {
            state_after: after,
            axial_momentum_change_nms: axial_change,
            external_angular_impulse_world: external_impulse,
            shaft_work_j: shaft_work,
            external_work_j: external_work,
            kinetic_energy_change_j: kinetic_change,
            drift_energy_defect_j: drift_defect,
            work_energy_residual_j: work_residual(
                kinetic_change,
                &shaft_work,
                external_work,
                drift_defect,
            )?,
            angular_momentum_residual_world: residual_world,
        })
    }

    fn admit_state(&self, state: &RotorState<N>) -> Result<(), LieStepError> {
        vector(state.total_momentum_body, "total body momentum")?;
        for value in state
            .axial_momentum_nms
            .iter()
            .chain(&state.relative_angle_rad)
        {
            scalar(*value, "finite rotor momentum and phase")?;
        }
        Ok(())
    }

    fn rotor_momentum(&self, axial: &[f64; N]) -> Result<Vec3, LieStepError> {
        let mut sum = [dd(0.0); 3];
        for (rotor, p) in self.rotors.iter().zip(axial) {
            scalar(*p, "absolute rotor axial momentum")?;
            for (value, a) in sum.iter_mut().zip(components(rotor.axis_body)) {
                *value = *value + dd(a) * dd(*p);
            }
        }
        vector(
            Vec3::new(sum[0].to_f64(), sum[1].to_f64(), sum[2].to_f64()),
            "stored rotor vector",
        )
    }

    fn carrier_momentum(&self, state: &RotorState<N>) -> Result<Vec3, LieStepError> {
        vector(
            state.total_momentum_body - self.rotor_momentum(&state.axial_momentum_nms)?,
            "reduced carrier momentum",
        )
    }

    fn admit_recoil(
        &self,
        before_carrier: Vec3,
        after_carrier: Vec3,
        external: Vec3,
        axial_change: &[f64; N],
    ) -> Result<(), LieStepError> {
        // Each rotor can resolve its own impulse while a much larger rotor
        // erases that transfer in the aggregate h passed to the native drift.
        // Check the transfer independently of both endpoint sums and energy:
        // at zero carrier rate, lost recoil can falsely report zero work error.
        let mut required = components(external).map(dd);
        for (rotor, delta) in self.rotors.iter().zip(axial_change) {
            for (value, axis) in required.iter_mut().zip(components(rotor.axis_body)) {
                *value = *value - dd(axis) * dd(*delta);
            }
        }
        for ((before, after), expected) in components(before_carrier)
            .into_iter()
            .zip(components(after_carrier))
            .zip(required)
        {
            let expected = scalar(expected.to_f64(), "carrier impulse balance")?;
            let actual = scalar((dd(after) - dd(before)).to_f64(), "carrier impulse change")?;
            // Retain the lane's existing carrier-resolution allowance, with no
            // |h| term that could admit an erased (or materially partial) kick.
            let scale = before.abs().max(after.abs()).max(expected.abs());
            if (expected != 0.0 && actual == 0.0)
                || (actual - expected).abs() > 4096.0 * f64::EPSILON * scale
            {
                return Err(invalid("unresolved carrier recoil against rotor momentum"));
            }
        }
        Ok(())
    }

    fn energy_change(
        &self,
        before: &RotorState<N>,
        after: &RotorState<N>,
    ) -> Result<f64, LieStepError> {
        let p0 = self.carrier_momentum(before)?;
        let p1 = self.carrier_momentum(after)?;
        let mean_omega =
            self.inertia.velocity(p0)?.scale(0.5) + self.inertia.velocity(p1)?.scale(0.5);
        let mut change = dot_dd(mean_omega, p1 - p0);
        for i in 0..N {
            let average = 0.5 * before.axial_momentum_nms[i] + 0.5 * after.axial_momentum_nms[i];
            let delta = after.axial_momentum_nms[i] - before.axial_momentum_nms[i];
            change = change + dd(average / self.rotors[i].axial_inertia_kg_m2) * dd(delta);
        }
        scalar(change.to_f64(), "rotor kinetic increment")
    }
}

fn invalid(context: &'static str) -> LieStepError {
    LieStepError::InvalidParameter { context }
}
fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}
fn components(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
fn finite(value: Vec3) -> bool {
    components(value).iter().all(|v| v.is_finite())
}
fn scalar(value: f64, context: &'static str) -> Result<f64, LieStepError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(context))
    }
}
fn vector(value: Vec3, context: &'static str) -> Result<Vec3, LieStepError> {
    if finite(value) {
        Ok(value)
    } else {
        Err(invalid(context))
    }
}
fn product(a: f64, b: f64, context: &'static str) -> Result<f64, LieStepError> {
    scalar(a, context)?;
    scalar(b, context)?;
    let value = scalar(a * b, context)?;
    if a != 0.0 && b != 0.0 && value == 0.0 {
        return Err(invalid(context));
    }
    Ok(value)
}
fn increment(value: f64, delta: f64, context: &'static str) -> Result<f64, LieStepError> {
    scalar(delta, context)?;
    let next = scalar(value + delta, context)?;
    if delta != 0.0 && next == value {
        return Err(invalid(context));
    }
    Ok(next)
}
fn dot_dd(a: Vec3, b: Vec3) -> Dd {
    components(a)
        .into_iter()
        .zip(components(b))
        .fold(dd(0.0), |sum, (a, b)| sum + dd(a) * dd(b))
}
fn subtract3(a: Vec3, b: Vec3, c: Vec3) -> Vec3 {
    let a = components(a);
    let b = components(b);
    let c = components(c);
    let values: [f64; 3] = core::array::from_fn(|i| (dd(a[i]) - dd(b[i]) - dd(c[i])).to_f64());
    Vec3::new(values[0], values[1], values[2])
}
fn work_residual<const N: usize>(
    change: f64,
    shaft: &[f64; N],
    external: f64,
    drift: f64,
) -> Result<f64, LieStepError> {
    let mut value = dd(change) - dd(external) - dd(drift);
    for work in shaft {
        value = value - dd(*work);
    }
    scalar(value.to_f64(), "rotor work balance")
}
