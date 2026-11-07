//! Reciprocal finite-gap mesh compliance on original rotor coordinates.
//!
//! x = q_first - ratio*q_second - center, with half-gap g. The two
//! penetrations are z+ = max(x-g,0), z- = max(-x-g,0), and
//! V = k+*z+^2/2 + k-*z-^2/2. Compression-only dashpots oppose further
//! penetration; unloading remains elastic and never attracts the teeth.
//! A source-owned potential secant and path-averaged compression loss close
//! work even when a subflow enters/leaves the gap or changes flanks.

use super::{
    HeldTorque, LieStepError, RotorBank, RotorImpulse, RotorState, RotorStep, components, dd,
    dot_dd, increment, invalid, product, scalar, vector,
};
use fs_ga::{So3, Vec3};
use fs_math::dd::Dd;

const ROUNDING: f64 = 8192.0 * f64::EPSILON;
/// Maximum scalar bisections per shaft half-stage; affine cases need none.
pub const MAX_BACKLASH_ROOT_ITERATIONS: u8 = 96;

/// Authored reduced flank law; these are not tooth-geometry/material estimates.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct BacklashFlank {
    pub stiffness_nm_per_rad: f64,
    /// Viscous cushioning ONLY while penetration increases [N m s/rad].
    /// Elastic energy is returned without dashpot force while unloading.
    pub compression_damping_nms_per_rad: f64,
}

/// One fixed signed-ratio mesh, in unwrapped carrier-relative shaft angles.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BacklashCoupling {
    pub first_rotor: usize,
    pub second_rotor: usize,
    /// Port strain is q_first - ratio*q_second - center_angle_rad.
    /// Positive net effort gives shaft torques (-effort, ratio*effort).
    pub ratio: f64,
    pub center_angle_rad: f64,
    /// The unloaded interval is [-half_gap_rad, +half_gap_rad].
    pub half_gap_rad: f64,
    pub positive_flank: BacklashFlank,
    pub negative_flank: BacklashFlank,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BacklashContact {
    /// Includes the two zero-penetration boundaries.
    Gap,
    PositiveFlank,
    NegativeFlank,
}

/// Cancellation/refusal returns no partially advanced state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BacklashStepError {
    Native(LieStepError),
    Cancelled,
    /// The bounded monotone scalar solve did not resolve its original equation.
    ScalarSolve,
}

impl From<LieStepError> for BacklashStepError {
    fn from(value: LieStepError) -> Self {
        Self::Native(value)
    }
}
impl core::fmt::Display for BacklashStepError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Native(error) => error.fmt(f),
            Self::Cancelled => f.write_str("backlash step cancelled"),
            Self::ScalarSolve => f.write_str("unresolved bounded backlash scalar solve"),
        }
    }
}
impl std::error::Error for BacklashStepError {}

/// One prepared finite-compliance backlash port, with complete carrier recoil.
///
/// Two shaft half-stages surround the original gyrostat drift. Each solves one
/// strongly monotone scalar equation; free-gap and unchanged-flank affine cases
/// solve directly, crossings use at most 96 bisections. Preparation/steps are
/// O(N), allocation-free and bounded by the bank's 64 rotors. Cancellation is
/// polled at every bisection and between the bounded native stages.
///
/// This is a discrete path approximation, not a resolved contact-event
/// trajectory. In particular, zero average motion need not mean no physical
/// compression/unloading inside a subflow. Refine gap crossings and reversals;
/// no exact event time or nonsmooth convergence order is asserted.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BacklashStepper<const N: usize> {
    bank: RotorBank<N>,
    coupling: BacklashCoupling,
    duration_s: f64,
    mobility: Dd,
    alpha: Dd,
    resolution: f64,
    active: bool,
}

/// Actual mechanical motion and the two source-owned flank work accounts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BacklashStep<const N: usize> {
    pub motion: RotorStep<N>,
    pub drive_work_j: [f64; N],
    /// Exact polarized potential increments on actual stored phases, [+,-].
    pub elastic_energy_change_j: [f64; 2],
    /// Nonnegative compression-dashpot heat, in [+,-] flank order.
    pub dissipated_energy_j: [f64; 2],
    /// Nonnegative compressive impulse magnitudes of the discrete paths, [+,-].
    /// Signed net port impulse is positive minus negative; this is not a peak
    /// tooth load or a resolved time history of engagement.
    pub compressive_angular_impulse_nms: [f64; 2],
    pub contact_after: BacklashContact,
    /// Bisections used by each half-stage, zero for a direct affine solve.
    pub root_iterations: [u8; 2],
    /// Delta kinetic + Delta elastic + heat - drive - external - drift defect.
    pub work_energy_residual_j: f64,
    /// Sum of magnitudes of measured shaft-momentum addition, carrier momentum
    /// reconstruction and phase-storage work rounding in both stages. This
    /// numerical allowance is neither heat nor a source, and is NOT subtracted
    /// from the reported raw work residual.
    pub work_roundoff_allowance_j: f64,
}

impl BacklashCoupling {
    fn flanks(self) -> [BacklashFlank; 2] {
        [self.positive_flank, self.negative_flank]
    }
    fn ports(self) -> [(usize, f64); 2] {
        [(self.first_rotor, 1.0), (self.second_rotor, -self.ratio)]
    }
    fn mapped<const N: usize>(self, values: &[f64; N]) -> Dd {
        dd(values[self.first_rotor]) - dd(self.ratio) * dd(values[self.second_rotor])
    }
    fn coordinate<const N: usize>(self, phase: &[f64; N]) -> Dd {
        self.mapped(phase) - dd(self.center_angle_rad)
    }
    fn penetration(self, coordinate: Dd) -> [Dd; 2] {
        [
            positive(coordinate - dd(self.half_gap_rad)),
            positive(-coordinate - dd(self.half_gap_rad)),
        ]
    }
    fn potential_change(self, before: Dd, after: Dd) -> Result<[f64; 2], LieStepError> {
        let a = self.penetration(before);
        let b = self.penetration(after);
        let mut change = [0.0; 2];
        for (side, flank) in self.flanks().iter().enumerate() {
            change[side] = scalar(
                (dd(0.5)
                    * dd(flank.stiffness_nm_per_rad)
                    * (a[side] + b[side])
                    * (b[side] - a[side]))
                    .to_f64(),
                "backlash original elastic increment",
            )?;
        }
        Ok(change)
    }
    fn contact(self, x: Dd) -> BacklashContact {
        if dd(self.half_gap_rad).lt(x) {
            BacklashContact::PositiveFlank
        } else if x.lt(dd(-self.half_gap_rad)) {
            BacklashContact::NegativeFlank
        } else {
            BacklashContact::Gap
        }
    }

    // Integrate each elastic force along the declared straight phase path.
    // The secant remains signed correctly when crossing either/both flanks.
    fn effort(self, start: Dd, change: Dd, duration: f64) -> Result<Effort, LieStepError> {
        let mut result = Effort {
            elastic: [dd(0.0); 2],
            damping: [dd(0.0); 2],
        };
        for (side, flank) in self.flanks().iter().enumerate() {
            let sign = if side == 0 { 1.0 } else { -1.0 };
            let z0 = dd(sign) * start - dd(self.half_gap_rad);
            let dz = dd(sign) * change;
            let z1 = z0 + dz;
            let u0 = positive(z0);
            let u1 = positive(z1);
            let compression = if !z0.lt(dd(0.0)) && !z1.lt(dd(0.0)) {
                // Keep a small resolvable change beside a larger stored strain.
                result.elastic[side] = dd(flank.stiffness_nm_per_rad) * (z0 + dd(0.5) * dz);
                dz
            } else {
                let delta = u1 - u0;
                result.elastic[side] = if change == dd(0.0) {
                    dd(flank.stiffness_nm_per_rad) * u0
                } else {
                    dd(sign) * dd(0.5) * dd(flank.stiffness_nm_per_rad) * (u0 + u1) * delta / change
                };
                delta
            };
            result.damping[side] =
                dd(flank.compression_damping_nms_per_rad) * positive(compression) / dd(duration);
            for value in [result.elastic[side], result.damping[side]] {
                if scalar(value.to_f64(), "finite backlash flank effort")? < 0.0 {
                    return Err(invalid("nonnegative compressive backlash effort"));
                }
            }
        }
        Ok(result)
    }
}

impl<const N: usize> RotorBank<N> {
    /// Prepare fixed backlash geometry and two passive authored flank laws.
    /// With delta=duration/2, the local diagnostic is
    /// delta*H*c_max + delta^2*H*k_max/2 <= 2. This is a resolution limit,
    /// not a contact-event, accuracy or material-validity certificate.
    pub fn backlash_stepper(
        &self,
        coupling: BacklashCoupling,
        duration_s: f64,
    ) -> Result<BacklashStepper<N>, LieStepError> {
        scalar(duration_s, "backlash duration")?;
        let delta = 0.5 * duration_s;
        if delta <= 0.0
            || coupling.first_rotor >= N
            || coupling.second_rotor >= N
            || coupling.first_rotor == coupling.second_rotor
        {
            return Err(invalid(
                "positive backlash duration and distinct original shafts",
            ));
        }
        scalar(coupling.ratio, "backlash ratio")?;
        scalar(coupling.center_angle_rad, "backlash center phase")?;
        if scalar(coupling.half_gap_rad, "backlash half gap")? < 0.0 {
            return Err(invalid("nonnegative backlash half gap"));
        }
        let (mut k, mut c) = (0.0_f64, 0.0_f64);
        for flank in coupling.flanks() {
            let stiffness = scalar(flank.stiffness_nm_per_rad, "backlash stiffness")?;
            let damping = scalar(flank.compression_damping_nms_per_rad, "backlash damping")?;
            if stiffness < 0.0 || damping < 0.0 {
                return Err(invalid("nonnegative backlash flank coefficients"));
            }
            k = k.max(stiffness);
            c = c.max(damping);
        }
        let mut axis = [dd(0.0); 3];
        let mut mobility = dd(0.0);
        for (i, b) in coupling.ports() {
            mobility = mobility + dd(b) * dd(b) / dd(self.rotors[i].axial_inertia_kg_m2);
            for (component, value) in axis.iter_mut().zip(components(self.rotors[i].axis_body)) {
                *component = *component + dd(b) * dd(value);
            }
        }
        let axis = vector(
            Vec3::new(axis[0].to_f64(), axis[1].to_f64(), axis[2].to_f64()),
            "backlash recoil axis",
        )?;
        mobility = mobility + dot_dd(axis, self.inertia.velocity(axis)?);
        if scalar(mobility.to_f64(), "backlash port mobility")? <= 0.0 {
            return Err(invalid("positive backlash port mobility"));
        }
        let alpha = dd(0.5 * delta) * (dd(delta) * mobility);
        let resolution = scalar(
            ((dd(delta) * mobility * dd(c) + alpha * dd(k)) * dd(1.0 + 256.0 * f64::EPSILON))
                .to_f64(),
            "backlash half-stage resolution",
        )?;
        let active = k != 0.0 || c != 0.0;
        if resolution > 2.0 || (active && scalar(alpha.to_f64(), "backlash stage response")? == 0.0)
        {
            return Err(invalid("resolved backlash stiffness and damping"));
        }
        Ok(BacklashStepper {
            bank: *self,
            coupling,
            duration_s,
            mobility,
            alpha,
            resolution,
            active,
        })
    }
}

impl<const N: usize> BacklashStepper<N> {
    pub const fn coupling(&self) -> BacklashCoupling {
        self.coupling
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    pub const fn mobility(&self) -> f64 {
        self.mobility.to_f64()
    }
    pub const fn half_stage_resolution_number(&self) -> f64 {
        self.resolution
    }

    pub fn contact(&self, state: &RotorState<N>) -> Result<BacklashContact, LieStepError> {
        self.bank.admit_state(state)?;
        let x = self.coupling.coordinate(&state.relative_angle_rad);
        scalar(x.to_f64(), "backlash mesh coordinate")?;
        Ok(self.coupling.contact(x))
    }
    pub fn potential_energy(&self, state: &RotorState<N>) -> Result<f64, LieStepError> {
        self.bank.admit_state(state)?;
        let x = self.coupling.coordinate(&state.relative_angle_rad);
        let z = self.coupling.penetration(x);
        let mut energy = dd(0.0);
        for (side, flank) in self.coupling.flanks().iter().enumerate() {
            energy = energy + dd(0.5) * dd(flank.stiffness_nm_per_rad) * z[side] * z[side];
        }
        scalar(energy.to_f64(), "backlash original potential")
    }
    pub fn step(
        &self,
        state: &RotorState<N>,
        drive_torque_nm: [f64; N],
        external_torque: HeldTorque,
    ) -> Result<BacklashStep<N>, BacklashStepError> {
        self.step_controlled(state, drive_torque_nm, external_torque, || false)
    }
    pub fn step_controlled(
        &self,
        state: &RotorState<N>,
        drive_torque_nm: [f64; N],
        external_torque: HeldTorque,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BacklashStep<N>, BacklashStepError> {
        poll(&mut cancelled)?;
        if !self.active {
            let motion =
                self.bank
                    .step(state, drive_torque_nm, external_torque, self.duration_s)?;
            let contact_after = self.contact(&motion.state_after)?;
            poll(&mut cancelled)?;
            return Ok(BacklashStep {
                drive_work_j: motion.shaft_work_j,
                work_energy_residual_j: motion.work_energy_residual_j,
                motion,
                elastic_energy_change_j: [0.0; 2],
                dissipated_energy_j: [0.0; 2],
                compressive_angular_impulse_nms: [0.0; 2],
                contact_after,
                root_iterations: [0; 2],
                work_roundoff_allowance_j: 0.0,
            });
        }
        let first = self.stage(
            state,
            drive_torque_nm,
            self.external_impulse(state.rotation, external_torque)?,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        let mut drifted = first.impulse.state_after;
        drifted.relative_angle_rad = first.phase;
        (drifted.rotation, drifted.total_momentum_body) = self.bank.inertia.momentum_step(
            drifted.rotation,
            drifted.total_momentum_body,
            self.bank.rotor_momentum(&drifted.axial_momentum_nms)?,
            self.duration_s,
        )?;
        poll(&mut cancelled)?;
        let second = self.stage(
            &drifted,
            drive_torque_nm,
            self.external_impulse(drifted.rotation, external_torque)?,
            &mut cancelled,
        )?;
        let motion = self.bank.finish_step(
            state,
            &first.impulse,
            &drifted,
            &second.impulse,
            self.duration_s,
            Some(second.phase),
        )?;
        let elastic_energy_change_j = self.coupling.potential_change(
            self.coupling.coordinate(&state.relative_angle_rad),
            self.coupling.coordinate(&second.phase),
        )?;
        let mut drive_work_j = [0.0; N];
        let mut dissipated_energy_j = [0.0; 2];
        let mut compressive_angular_impulse_nms = [0.0; 2];
        let mut balance = dd(motion.kinetic_energy_change_j)
            - dd(motion.external_work_j)
            - dd(motion.drift_energy_defect_j);
        for i in 0..N {
            drive_work_j[i] = scalar(
                (dd(first.drive_work[i]) + dd(second.drive_work[i])).to_f64(),
                "backlash drive work",
            )?;
            balance = balance - dd(drive_work_j[i]);
        }
        for side in 0..2 {
            dissipated_energy_j[side] = scalar(
                (dd(first.heat[side]) + dd(second.heat[side])).to_f64(),
                "backlash heat",
            )?;
            compressive_angular_impulse_nms[side] = scalar(
                (dd(first.compressive_impulse[side]) + dd(second.compressive_impulse[side]))
                    .to_f64(),
                "backlash compressive impulse",
            )?;
            balance = balance + dd(elastic_energy_change_j[side]) + dd(dissipated_energy_j[side]);
        }
        let result = BacklashStep {
            contact_after: self.contact(&motion.state_after)?,
            motion,
            drive_work_j,
            elastic_energy_change_j,
            dissipated_energy_j,
            compressive_angular_impulse_nms,
            root_iterations: [first.iterations, second.iterations],
            work_energy_residual_j: scalar(balance.to_f64(), "backlash complete work balance")?,
            work_roundoff_allowance_j: scalar(
                (dd(first.roundoff_allowance) + dd(second.roundoff_allowance)).to_f64(),
                "backlash complete work roundoff allowance",
            )?,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    fn external_impulse(&self, rotation: So3, torque: HeldTorque) -> Result<Vec3, LieStepError> {
        let body = match torque {
            HeldTorque::Body(value) => vector(value, "backlash external body torque")?,
            HeldTorque::Space(value) => rotation.inverse().rotate(value)?,
        };
        let delta = 0.5 * self.duration_s;
        Ok(Vec3::new(
            product(delta, body.x, "backlash external impulse")?,
            product(delta, body.y, "backlash external impulse")?,
            product(delta, body.z, "backlash external impulse")?,
        ))
    }

    fn solved(&self, x: Dd, y: f64, free: Dd) -> Result<Option<Effort>, LieStepError> {
        let effort = self.coupling.effort(x, dd(y), 0.5 * self.duration_s)?;
        let residual = scalar(
            (dd(y) + self.alpha * effort.net() - free).to_f64(),
            "backlash scalar residual",
        )?;
        let scale = y
            .abs()
            .max(free.to_f64().abs())
            .max((self.alpha * effort.magnitude(0)).to_f64().abs())
            .max((self.alpha * effort.magnitude(1)).to_f64().abs());
        scalar(scale, "backlash scalar residual scale")?;
        Ok((residual.abs() <= 64.0 * f64::EPSILON * scale).then_some(effort))
    }

    fn solve(
        &self,
        x: Dd,
        free: Dd,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Effort, u8), BacklashStepError> {
        let delta = 0.5 * self.duration_s;
        let at_zero = self.coupling.effort(x, dd(0.0), delta)?;
        let base = free - self.alpha * at_zero.net();
        let bound = scalar(base.to_f64(), "backlash monotone root bracket")?;
        // Convex secants and compression damping are nondecreasing in y, so
        // F(y)=y+alpha*effort(y)-free has slope >=1. Its root lies between
        // zero and -F(0), without an unbounded bracket search.
        if let Some(effort) = self.solved(x, free.to_f64(), free)? {
            return Ok((effort, 0));
        }
        for (side, flank) in self.coupling.flanks().iter().enumerate() {
            let sign = if side == 0 { 1.0 } else { -1.0 };
            let z = dd(sign) * x - dd(self.coupling.half_gap_rad);
            if z.lt(dd(0.0)) {
                continue;
            }
            let damping = if sign * bound > 0.0 {
                flank.compression_damping_nms_per_rad
            } else {
                0.0
            };
            let denominator = dd(1.0)
                + self.alpha * (dd(0.5 * flank.stiffness_nm_per_rad) + dd(damping) / dd(delta));
            let candidate = scalar((base / denominator).to_f64(), "backlash affine root")?;
            if !(z + dd(sign * candidate)).lt(dd(0.0)) {
                if let Some(effort) = self.solved(x, candidate, free)? {
                    return Ok((effort, 0));
                }
            }
        }
        let (mut lo, mut hi) = (bound.min(0.0), bound.max(0.0));
        for iterations in 1..=MAX_BACKLASH_ROOT_ITERATIONS {
            poll(cancelled)?;
            let candidate = f64::midpoint(lo, hi);
            if let Some(effort) = self.solved(x, candidate, free)? {
                return Ok((effort, iterations));
            }
            if candidate == lo || candidate == hi {
                break;
            }
            let effort = self.coupling.effort(x, dd(candidate), delta)?;
            let residual = dd(candidate) + self.alpha * effort.net() - free;
            if residual.lt(dd(0.0)) {
                lo = candidate;
            } else {
                hi = candidate;
            }
        }
        Err(BacklashStepError::ScalarSolve)
    }

    fn stage(
        &self,
        state: &RotorState<N>,
        drive: [f64; N],
        external: Vec3,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Stage<N>, BacklashStepError> {
        poll(cancelled)?;
        self.bank.admit_state(state)?;
        let delta = 0.5 * self.duration_s;
        let initial_rates = self.bank.relative_rates(state)?;
        let next_total =
            super::super::super::torque_kick(state.total_momentum_body, external, 1.0)?;
        let mut free_carrier = components(next_total - state.total_momentum_body).map(dd);
        let mut drive_impulse = [dd(0.0); N];
        for i in 0..N {
            product(delta, drive[i], "backlash drive impulse")?;
            drive_impulse[i] = dd(delta) * dd(drive[i]);
            for (component, axis) in free_carrier
                .iter_mut()
                .zip(components(self.bank.rotors[i].axis_body))
            {
                *component = *component - dd(axis) * drive_impulse[i];
            }
        }
        let dv = self.bank.inertia.velocity(Vec3::new(
            free_carrier[0].to_f64(),
            free_carrier[1].to_f64(),
            free_carrier[2].to_f64(),
        ))?;
        let mut free_rate = self.coupling.mapped(&initial_rates);
        for (i, b) in self.coupling.ports() {
            free_rate = free_rate
                + dd(0.5 * b)
                    * (drive_impulse[i] / dd(self.bank.rotors[i].axial_inertia_kg_m2)
                        - dot_dd(self.bank.rotors[i].axis_body, dv));
        }
        let x = self.coupling.coordinate(&state.relative_angle_rad);
        scalar(x.to_f64(), "backlash initial coordinate")?;
        let free = dd(delta) * free_rate;
        scalar(free.to_f64(), "backlash free phase motion")?;
        let (effort, iterations) = self.solve(x, free, cancelled)?;
        let mut applied = drive_impulse;
        for (i, b) in self.coupling.ports() {
            applied[i] = applied[i] - dd(delta) * dd(b) * effort.net();
        }
        let requested = applied.map(|v| v.to_f64());
        let impulse = self.bank.apply_impulse(state, requested, external)?;
        let final_rates = self.bank.relative_rates(&impulse.state_after)?;
        let mut phase = state.relative_angle_rad;
        let mut actual_angle = [dd(0.0); N];
        for i in 0..N {
            let angle = product(
                delta,
                impulse.mean_relative_rate_rad_s[i],
                "backlash phase increment",
            )?;
            let reverses = (initial_rates[i] > 0.0 && final_rates[i] < 0.0)
                || (initial_rates[i] < 0.0 && final_rates[i] > 0.0);
            let scale = if reverses {
                (dd(0.5 * delta) * (dd(initial_rates[i].abs()) + dd(final_rates[i].abs()))).to_f64()
            } else {
                angle.abs()
            };
            phase[i] = if reverses {
                scalar(phase[i] + angle, "backlash stored phase")?
            } else {
                increment(phase[i], angle, "backlash stored phase")?
            };
            actual_angle[i] = dd(phase[i]) - dd(state.relative_angle_rad[i]);
            admit(
                (actual_angle[i] - dd(angle)).to_f64(),
                scale,
                "backlash original phase resolution",
            )?;
        }
        let actual_change = actual_angle[self.coupling.first_rotor]
            - dd(self.coupling.ratio) * actual_angle[self.coupling.second_rotor];
        let actual_effort = self.coupling.effort(x, actual_change, delta)?;
        let elastic = self
            .coupling
            .potential_change(x, self.coupling.coordinate(&phase))?;
        let mut expected = drive_impulse;
        let mut law_scale = drive_impulse.map(|value| value.to_f64().abs());
        for (i, b) in self.coupling.ports() {
            expected[i] = expected[i] - dd(delta) * dd(b) * actual_effort.net();
            law_scale[i] = law_scale[i]
                .max((dd(delta * b) * actual_effort.magnitude(0)).to_f64().abs())
                .max((dd(delta * b) * actual_effort.magnitude(1)).to_f64().abs());
        }
        let mut drive_work = [0.0; N];
        let mut heat = [0.0; 2];
        let mut compressive_impulse = [0.0; 2];
        let mut balance = dd(impulse.kinetic_energy_change_j) - dd(impulse.external_work_j);
        let mut work_scale = impulse
            .kinetic_energy_change_j
            .abs()
            .max(impulse.external_work_j.abs());
        let mut roundoff = dd(0.0);
        for i in 0..N {
            let actual = impulse.axial_momentum_change_nms[i];
            if requested[i] != 0.0 && actual == 0.0 {
                return Err(invalid("unresolved backlash shaft reaction").into());
            }
            admit(
                (dd(requested[i]) - expected[i]).to_f64(),
                law_scale[i].max(requested[i].abs()),
                "backlash original shaft law",
            )?;
            // Let P be the requested f64 impulse, P* the original law on the
            // stored phase, and e_q=Delta q_stored-delta*vbar. The raw work
            // residual contains (Delta p-P)*vbar - (P*/delta)*e_q, plus the
            // separately admitted law mismatch (P-P*)*vbar. Measure the two
            // storage terms componentwise; never use stored spin or energy as
            // a permissive scale. Native erased-transfer/recoil refusal above
            // remains authoritative even when this roundoff is measurable.
            let mean_rate = dd(impulse.mean_relative_rate_rad_s[i]);
            let phase_error = actual_angle[i] - dd(delta) * mean_rate;
            roundoff = roundoff
                + ((dd(actual) - dd(requested[i])) * mean_rate).abs()
                + ((expected[i] / dd(delta)) * phase_error).abs();
            drive_work[i] = scalar(
                (dd(drive[i]) * actual_angle[i]).to_f64(),
                "backlash original drive work",
            )?;
            balance = balance - dd(drive_work[i]);
            work_scale = work_scale
                .max(drive_work[i].abs())
                .max(impulse.shaft_work_j[i].abs())
                .max((expected[i] * mean_rate).to_f64().abs());
        }
        // The native carrier coordinate k=L-A*p is itself reconstructed in
        // finite storage. Its already-admitted recoil discrepancy contributes
        // omega_bar.(Delta k-Delta L+A*Delta p) to the work identity. Bound the
        // actual component errors, not |A*p| or the carrier's stored energy.
        let before_carrier = self.bank.carrier_momentum(state)?;
        let after_carrier = self.bank.carrier_momentum(&impulse.state_after)?;
        let mean_omega = self.bank.inertia.velocity(before_carrier)?.scale(0.5)
            + self.bank.inertia.velocity(after_carrier)?.scale(0.5);
        let mut required_carrier = components(impulse.external_angular_impulse_body).map(dd);
        for (rotor, transfer) in self
            .bank
            .rotors
            .iter()
            .zip(impulse.axial_momentum_change_nms)
        {
            for (component, axis) in required_carrier.iter_mut().zip(components(rotor.axis_body)) {
                *component = *component - dd(axis) * dd(transfer);
            }
        }
        for (((before, after), required), velocity) in components(before_carrier)
            .into_iter()
            .zip(components(after_carrier))
            .zip(required_carrier)
            .zip(components(mean_omega))
        {
            let actual_change = dd(after) - dd(before);
            roundoff = roundoff + ((actual_change - required) * dd(velocity)).abs();
            work_scale = work_scale
                .max((actual_change * dd(velocity)).to_f64().abs())
                .max((required * dd(velocity)).to_f64().abs());
        }
        for side in 0..2 {
            let sign = if side == 0 { 1.0 } else { -1.0 };
            heat[side] = scalar(
                (dd(sign) * actual_effort.damping[side] * actual_change).to_f64(),
                "backlash physical compression heat",
            )?;
            compressive_impulse[side] = scalar(
                (dd(delta) * actual_effort.magnitude(side)).to_f64(),
                "backlash physical compressive impulse",
            )?;
            if heat[side] < 0.0 || compressive_impulse[side] < 0.0 {
                return Err(invalid("passive non-attractive backlash flank").into());
            }
            balance = balance + dd(heat[side]) + dd(elastic[side]);
            work_scale = work_scale.max(heat[side]).max(elastic[side].abs());
        }
        let roundoff_allowance = scalar(roundoff.to_f64(), "backlash measured work roundoff")?;
        let work_limit = scalar(
            roundoff_allowance + ROUNDING * work_scale,
            "backlash original work admission limit",
        )?;
        if scalar(balance.to_f64(), "backlash original potential work")?.abs() > work_limit {
            return Err(invalid("backlash original potential work").into());
        }
        Ok(Stage {
            impulse,
            phase,
            drive_work,
            heat,
            compressive_impulse,
            iterations,
            roundoff_allowance,
        })
    }
}

#[derive(Clone, Copy)]
struct Effort {
    elastic: [Dd; 2],
    damping: [Dd; 2],
}
impl Effort {
    fn magnitude(self, side: usize) -> Dd {
        self.elastic[side] + self.damping[side]
    }
    fn net(self) -> Dd {
        self.magnitude(0) - self.magnitude(1)
    }
}
struct Stage<const N: usize> {
    impulse: RotorImpulse<N>,
    phase: [f64; N],
    drive_work: [f64; N],
    heat: [f64; 2],
    compressive_impulse: [f64; 2],
    iterations: u8,
    roundoff_allowance: f64,
}
fn positive(value: Dd) -> Dd {
    if value.lt(dd(0.0)) { dd(0.0) } else { value }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), BacklashStepError> {
    if cancelled() {
        Err(BacklashStepError::Cancelled)
    } else {
        Ok(())
    }
}
fn admit(residual: f64, scale: f64, context: &'static str) -> Result<(), LieStepError> {
    scalar(residual, context)?;
    scalar(scale, context)?;
    if residual.abs() > ROUNDING * scale {
        Err(invalid(context))
    } else {
        Ok(())
    }
}
