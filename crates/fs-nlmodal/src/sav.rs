//! Exact-linear scalar-auxiliary-variable stepping of the retained SOS model.
//!
//! Write `H = H_linear + V(q)`, where V is the complete quartic SOS potential,
//! and introduce `r = sqrt(V+C0)`, `g = grad(V)/sqrt(V+C0)`. The continuous
//! equations are `q'=p`, `p'=-Omega^2 q-Cp-g r+f`, `r'=g.p/2`.
//! One scalar represents the TOTAL potential; no stress channel is discarded.
//!
//! Freeze g at an exact-linear half-step predictor with the initial nonlinear
//! force held. Over the full step, hold `f-g*r_bar` and reuse fs-phs's exact
//! scalar oscillator trajectories, mean velocities and positive viscous heat.
//! With their free displacements d and force gains B, the only solve is
//! `r_bar = (r + g.d/4)/(1 + g.B.g/4)`. Then `delta_r = g.delta_q/2`.
//! Consequently `delta(H_linear+r^2-C0) = f.delta_q - heat`.
//!
//! This is a second-order exponential SAV discretization, not exact nonlinear
//! motion. The numerical invariant differs from the physical SOS energy;
//! every record reports that discrepancy. Time, mode and stress-basis refinement
//! remain necessary. No contact, full matrix viscosity, longitudinal-memory,
//! reciprocal termination or render-path migration is implied by this API.
//!
//! References: Jiang, Wang & Cai, JCP 419 (2020), 109690,
//! doi:10.1016/j.jcp.2020.109690 (exponential SAV); Russo, Bilbao & Ducceschi,
//! "Scalar Auxiliary Variable Techniques for Nonlinear Transverse String
//! Vibration", IFAC (2024), https://mdphys.org/PDF/IFAC2024_Strings.pdf.
//! The update and work identity above specify this held-force extension.

use crate::{NlModalError, SosModalStorage};
use fs_math::det;
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_phs::{PhsError, mass_spring_damper};

/// Maximum retained modal coordinates. Capacity never selects or drops modes.
pub const MAX_SAV_MODES: usize = 4096;

/// Explicit numerical inputs; the offset is an energy, not a stiffness change.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SavConfig {
    /// Positive finite C0 [J]. Its value affects finite-step auxiliary error.
    pub offset_j: f64,
    /// Caller-declared retained-mode capacity, at most [`MAX_SAV_MODES`].
    pub max_modes: usize,
    /// Maximum admitted omega*dt, positive and strictly below pi.
    /// This guards linear sampling only; it does not certify nonlinear bandwidth.
    pub max_omega_dt: f64,
}

/// A complete checkpoint, including the auxiliary degree of freedom.
/// Restart by cloning this value; recreating it from displacement/velocity
/// alone would reset the auxiliary energy defect and change the trajectory.
#[derive(Clone, Debug, PartialEq)]
pub struct SavState {
    coordinates: Vec<f64>,
    // Store r-sqrt(C0): tiny nonlinear energy must not be lost beside C0.
    auxiliary_shift: f64,
    offset_j: f64,
}

impl SavState {
    /// Original interleaved mass-normalized `[q0,p0,q1,p1,...]` coordinates.
    pub fn coordinates(&self) -> &[f64] {
        &self.coordinates
    }

    /// Auxiliary amplitude r [sqrt(J)], not sqrt(V+C0) after discrete steps.
    pub fn auxiliary_amplitude_sqrt_j(&self) -> f64 {
        det::sqrt(self.offset_j) + self.auxiliary_shift
    }
}

/// Physical heat/work and both distinct energies for the accepted endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SavRecord {
    /// H_linear + r^2 - C0 [J], evaluated without subtracting C0 from r^2.
    pub modified_energy_j: f64,
    /// The original H_linear + V(q) [J]. It is not the conserved quantity.
    pub physical_energy_j: f64,
    /// Modified minus physical energy [J]; a temporal-convergence observable.
    pub auxiliary_energy_defect_j: f64,
    /// Cancellation-stable change of the modified energy during this step [J].
    pub modified_energy_change_j: f64,
    /// Independent integral of C_i*p_i(t)^2 on the held-force trajectories [J].
    pub dissipated_j: f64,
    /// Held generalized force times its same-window displacement [J].
    pub port_work_j: f64,
    /// delta(modified energy) + heat - port work [J].
    pub balance_residual_j: f64,
    /// One scalar Sherman-Morrison solve if the model is nonlinear; zero otherwise.
    pub rank_one_solves: usize,
}

/// No state is published on any refusal, including cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum SavError {
    Model(NlModalError),
    Parameter {
        what: &'static str,
    },
    UnresolvedMode {
        mode: usize,
        omega_dt: f64,
        maximum: f64,
    },
    Linear(PhsError),
    EnergyBalance {
        residual_j: f64,
        scale_j: f64,
    },
    Cancelled,
}

impl core::fmt::Display for SavError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Model(error) => write!(f, "SAV storage: {error}"),
            Self::Parameter { what } => write!(f, "SAV requires {what}"),
            Self::UnresolvedMode {
                mode,
                omega_dt,
                maximum,
            } => write!(
                f,
                "SAV mode {mode} has omega*dt={omega_dt}, above declared {maximum}"
            ),
            Self::Linear(error) => write!(f, "SAV exact linear trajectory: {error}"),
            Self::EnergyBalance {
                residual_j,
                scale_j,
            } => write!(
                f,
                "SAV modified-energy residual {residual_j} J at scale {scale_j} J"
            ),
            Self::Cancelled => write!(f, "SAV step cancelled before publication"),
        }
    }
}
impl std::error::Error for SavError {}
impl From<NlModalError> for SavError {
    fn from(error: NlModalError) -> Self {
        Self::Model(error)
    }
}
impl From<PhsError> for SavError {
    fn from(error: PhsError) -> Self {
        Self::Linear(error)
    }
}
fn refuse(what: &'static str) -> SavError {
    SavError::Parameter { what }
}
fn finite(values: &[f64]) -> bool {
    values.iter().all(|v| v.is_finite())
}

/// Fixed-clock exact-linear SAV owner with allocation-free warm steps.
///
/// Diagonal stress matrices (including Kirchhoff-Carrier) cost O(N) per channel.
/// A dense channel retains its full O(N^2) multiplication. Each linear cache
/// is one original scalar oscillator; nonlinear coupling is retained through
/// the common SAV, without a dense 2N-state exponential or Newton iteration.
/// Damping coefficients are C_i [1/s] in `p_i' = ... - C_i*p_i`, not ratios.
pub struct SavStepper<'a> {
    storage: &'a SosModalStorage,
    config: SavConfig,
    dt: f64,
    root_offset: f64,
    nonlinear: bool,
    diagonal: Vec<bool>,
    linear: Vec<PreparedExactZoh<'static>>,
    // Exact half-interval mean velocity = a*q+b*p+beta*held_force.
    half_mean: Vec<[f64; 3]>,
    gains: Vec<f64>,
    candidate: Vec<f64>,
    q: Vec<f64>,
    gradient: Vec<f64>,
    product: Vec<f64>,
}

impl<'a> SavStepper<'a> {
    /// Admit every retained mode/channel and prepare exact scalar caches.
    /// There is no silent damping approximation or stiffness frequency matching.
    pub fn new(
        storage: &'a SosModalStorage,
        damping_per_s: &[f64],
        config: SavConfig,
        dt_s: f64,
    ) -> Result<Self, SavError> {
        let n = storage.omegas.len();
        if config.max_modes == 0
            || config.max_modes > MAX_SAV_MODES
            || n == 0
            || n > config.max_modes
        {
            return Err(refuse("a nonempty mode set within the explicit capacity"));
        }
        if !config.offset_j.is_finite()
            || config.offset_j <= 0.0
            || !dt_s.is_finite()
            || dt_s * 0.5 <= 0.0
            || !config.max_omega_dt.is_finite()
            || config.max_omega_dt <= 0.0
            || config.max_omega_dt >= core::f64::consts::PI
        {
            return Err(refuse(
                "positive finite offset/clock and 0 < max_omega_dt < pi",
            ));
        }
        storage.validate()?;
        if damping_per_s.len() != n || damping_per_s.iter().any(|c| !c.is_finite() || *c < 0.0) {
            return Err(refuse(
                "one finite nonnegative viscous coefficient per mode",
            ));
        }
        for (mode, &omega) in storage.omegas.iter().enumerate() {
            let omega_dt = omega * dt_s;
            if !omega_dt.is_finite() || omega_dt > config.max_omega_dt {
                return Err(SavError::UnresolvedMode {
                    mode,
                    omega_dt,
                    maximum: config.max_omega_dt,
                });
            }
        }
        let mut linear = Vec::with_capacity(n);
        let mut half_mean = Vec::with_capacity(n);
        let mut gains = Vec::with_capacity(n);
        let budget = LinearSplitBudget {
            max_modes: 2,
            max_omega_dt: config.max_omega_dt,
        };
        for (&omega, &damping) in storage.omegas.iter().zip(damping_per_s) {
            let full = PreparedExactZoh::new_owned(
                mass_spring_damper(1.0, omega * omega, damping)?,
                dt_s,
                budget,
            )?;
            let half = PreparedExactZoh::new_owned(
                mass_spring_damper(1.0, omega * omega, damping)?,
                dt_s * 0.5,
                budget,
            )?;
            let (a, beta) = half.average_relation(&[1.0, 0.0], &[0.0], 0)?;
            let (b, _) = half.average_relation(&[0.0, 1.0], &[0.0], 0)?;
            half_mean.push([a, b, beta]);
            let (_, beta) = full.average_relation(&[0.0, 0.0], &[0.0], 0)?;
            let gain = dt_s * beta;
            if !gain.is_finite() || gain < 0.0 {
                return Err(refuse(
                    "a finite nonnegative exact held-force displacement gain",
                ));
            }
            gains.push(gain);
            linear.push(full);
        }
        let diagonal = storage
            .channels
            .iter()
            .map(|channel| {
                (0..n).all(|i| (0..n).all(|j| i == j || channel.coupling[i * n + j] == 0.0))
            })
            .collect();
        let nonlinear = storage
            .channels
            .iter()
            .any(|channel| channel.coefficient > 0.0 && channel.coupling.iter().any(|e| *e != 0.0));
        Ok(Self {
            storage,
            config,
            dt: dt_s,
            root_offset: det::sqrt(config.offset_j),
            nonlinear,
            diagonal,
            linear,
            half_mean,
            gains,
            candidate: vec![0.0; 2 * n],
            q: vec![0.0; n],
            gradient: vec![0.0; n],
            product: vec![0.0; n],
        })
    }

    /// Initialize the auxiliary variable from the original physical potential.
    /// This is for a new trajectory; clone [`SavState`] for exact checkpointing.
    pub fn initial_state(&mut self, coordinates: &[f64]) -> Result<SavState, SavError> {
        if coordinates.len() != self.candidate.len() || !finite(coordinates) {
            return Err(refuse("complete finite interleaved initial coordinates"));
        }
        for (q, pair) in self.q.iter_mut().zip(coordinates.chunks_exact(2)) {
            *q = pair[0];
        }
        let potential = potential_gradient(
            self.storage,
            &self.diagonal,
            &self.q,
            &mut self.gradient,
            &mut self.product,
        )?;
        let root = det::sqrt(potential + self.config.offset_j);
        let shift = potential / (root + self.root_offset);
        if !root.is_finite()
            || !shift.is_finite()
            || !linear_energy(self.storage, coordinates).is_finite()
        {
            return Err(refuse(
                "representable initial physical and auxiliary energy",
            ));
        }
        Ok(SavState {
            coordinates: coordinates.to_vec(),
            auxiliary_shift: shift,
            offset_j: self.config.offset_j,
        })
    }

    /// Advance under one mass-normalized generalized force per mode.
    /// All accepted state, including its auxiliary variable, changes atomically.
    pub fn step(
        &mut self,
        state: &mut SavState,
        held_force: &[f64],
    ) -> Result<SavRecord, SavError> {
        self.step_controlled(state, held_force, || false)
    }

    /// Poll between complete modal operations and before committing the state.
    pub fn step_controlled(
        &mut self,
        state: &mut SavState,
        held_force: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SavRecord, SavError> {
        if cancelled() {
            return Err(SavError::Cancelled);
        }
        let n = self.storage.omegas.len();
        if state.coordinates.len() != 2 * n
            || held_force.len() != n
            || state.offset_j.to_bits() != self.config.offset_j.to_bits()
            || !finite(&state.coordinates)
            || !finite(held_force)
            || !state.auxiliary_shift.is_finite()
        {
            return Err(refuse(
                "a compatible complete finite checkpoint and force vector",
            ));
        }
        for (q, pair) in self.q.iter_mut().zip(state.coordinates.chunks_exact(2)) {
            *q = pair[0];
        }
        let old_potential = potential_gradient(
            self.storage,
            &self.diagonal,
            &self.q,
            &mut self.gradient,
            &mut self.product,
        )?;
        let r = self.root_offset + state.auxiliary_shift;
        if self.nonlinear {
            let root = det::sqrt(old_potential + self.config.offset_j);
            if !root.is_finite() {
                return Err(refuse("finite auxiliary predictor normalization"));
            }
            // Predictor uses the current EXTENDED state r, not a reset sqrt(V+C0).
            for i in 0..n {
                if cancelled() {
                    return Err(SavError::Cancelled);
                }
                let [a, b, beta] = self.half_mean[i];
                let q = state.coordinates[2 * i];
                let p = state.coordinates[2 * i + 1];
                let force = held_force[i] - (self.gradient[i] / root) * r;
                self.q[i] = q + 0.5 * self.dt * (a * q + b * p + beta * force);
            }
            let potential = potential_gradient(
                self.storage,
                &self.diagonal,
                &self.q,
                &mut self.gradient,
                &mut self.product,
            )?;
            let root = det::sqrt(potential + self.config.offset_j);
            if !root.is_finite() {
                return Err(refuse("finite midpoint auxiliary normalization"));
            }
            for g in &mut self.gradient {
                *g /= root;
            }
        } else {
            self.gradient.fill(0.0);
        }
        let mut numerator = r;
        let mut denominator = 1.0;
        for i in 0..n {
            if cancelled() {
                return Err(SavError::Cancelled);
            }
            let pair = &state.coordinates[2 * i..2 * i + 2];
            let (free, beta) = self.linear[i].average_relation(pair, &[0.0], 0)?;
            let displacement = self.dt * (free + beta * held_force[i]);
            let g = self.gradient[i];
            numerator += 0.25 * g * displacement;
            denominator += 0.25 * g * (self.gains[i] * g);
        }
        if !finite(&[numerator, denominator]) || denominator < 1.0 {
            return Err(refuse(
                "finite positive rank-one denominator and auxiliary numerator",
            ));
        }
        let r_bar = numerator / denominator;
        let mut auxiliary_increment = 0.0;
        let mut linear_change = 0.0;
        let mut heat = 0.0;
        let mut work = 0.0;
        for i in 0..n {
            if cancelled() {
                return Err(SavError::Cancelled);
            }
            let force = held_force[i] - self.gradient[i] * r_bar;
            let mut mean = [0.0];
            let record = self.linear[i].step_into(
                &state.coordinates[2 * i..2 * i + 2],
                &[force],
                &mut self.candidate[2 * i..2 * i + 2],
                &mut mean,
            )?;
            let displacement = self.dt * mean[0];
            auxiliary_increment += 0.5 * self.gradient[i] * displacement;
            linear_change += record.delta_h;
            heat += record.dissipated;
            work += held_force[i] * displacement;
        }
        let next_shift = state.auxiliary_shift + auxiliary_increment;
        let auxiliary_change = (next_shift - state.auxiliary_shift)
            * (2.0 * self.root_offset + next_shift + state.auxiliary_shift);
        let modified_change = linear_change + auxiliary_change;
        let auxiliary_energy = next_shift * (2.0 * self.root_offset + next_shift);
        for (q, pair) in self.q.iter_mut().zip(self.candidate.chunks_exact(2)) {
            *q = pair[0];
        }
        let potential = potential_gradient(
            self.storage,
            &self.diagonal,
            &self.q,
            &mut self.gradient,
            &mut self.product,
        )?;
        let linear = linear_energy(self.storage, &self.candidate);
        let physical = linear + potential;
        let modified = linear + auxiliary_energy;
        let residual = modified_change + heat - work;
        // C0 is deliberately excluded: a large arbitrary offset cannot hide
        // a failed balance at the physical or evolved auxiliary energy scale.
        let scale = (linear_energy(self.storage, &state.coordinates) + old_potential)
            .max(physical)
            .max(linear + auxiliary_energy.abs())
            .max(modified_change.abs())
            .max(heat.abs())
            .max(work.abs())
            .max(f64::MIN_POSITIVE);
        if !finite(&[next_shift, physical, modified, residual, scale]) || heat < 0.0 {
            return Err(refuse(
                "finite endpoint energies and independently nonnegative heat",
            ));
        }
        if residual.abs() > 1.0e-12 * scale {
            return Err(SavError::EnergyBalance {
                residual_j: residual,
                scale_j: scale,
            });
        }
        let record = SavRecord {
            modified_energy_j: modified,
            physical_energy_j: physical,
            auxiliary_energy_defect_j: auxiliary_energy - potential,
            modified_energy_change_j: modified_change,
            dissipated_j: heat,
            port_work_j: work,
            balance_residual_j: residual,
            rank_one_solves: usize::from(self.nonlinear),
        };
        if cancelled() {
            return Err(SavError::Cancelled);
        }
        state.coordinates.copy_from_slice(&self.candidate);
        state.auxiliary_shift = next_shift;
        Ok(record)
    }
}

fn linear_energy(storage: &SosModalStorage, x: &[f64]) -> f64 {
    storage
        .omegas
        .iter()
        .zip(x.chunks_exact(2))
        .map(|(&omega, pair)| {
            let scaled_q = omega * pair[0];
            0.5 * (scaled_q * scaled_q + pair[1] * pair[1])
        })
        .sum()
}

// Evaluate V and grad V directly. Subtracting the quadratic gradient would
// erase weak nonlinear restoring forces beside a high-frequency linear mode.
fn potential_gradient(
    storage: &SosModalStorage,
    diagonal: &[bool],
    q: &[f64],
    gradient: &mut [f64],
    product: &mut [f64],
) -> Result<f64, SavError> {
    if !finite(q) {
        return Err(refuse("finite nonlinear evaluation coordinates"));
    }
    let n = q.len();
    let mut potential = 0.0;
    gradient.fill(0.0);
    for (channel, &is_diagonal) in storage.channels.iter().zip(diagonal) {
        if channel.coefficient == 0.0 {
            continue;
        }
        for i in 0..n {
            product[i] = if is_diagonal {
                channel.coupling[i * n + i] * q[i]
            } else {
                // The antisymmetric part cancels from q^T E q. This is the
                // derivative of that same scalar even at admission roundoff.
                (0..n)
                    .map(|j| {
                        f64::midpoint(channel.coupling[i * n + j], channel.coupling[j * n + i])
                            * q[j]
                    })
                    .sum()
            };
        }
        let s: f64 = q.iter().zip(product.iter()).map(|(a, b)| a * b).sum();
        let weighted = channel.coefficient * s;
        potential += 0.25 * weighted * s;
        for (g, eq) in gradient.iter_mut().zip(product.iter()) {
            *g += weighted * eq;
        }
    }
    if !potential.is_finite() || potential < 0.0 || !finite(gradient) {
        return Err(refuse(
            "representable nonlinear potential and restoring force",
        ));
    }
    Ok(potential)
}
