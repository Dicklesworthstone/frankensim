//! Exact piecewise-linear stress control of the existing Maxwell history.
//!
//! With `x_j = sqrt(E_j/tau_j) q_j`, stress control gives
//! `xdot = -B x + v stress_rate/E0`, where
//! `B = diag(1/tau_j) - v v^T/E0` and `v_j = sqrt(E_j/tau_j)`.
//! This positive-definite symmetric matrix is diagonalized once. Each modal
//! trajectory is an affine-forcing exponential evaluated by the same scalar
//! owner as strain control. Heat independently integrates `sum(x_j^2)`;
//! external work integrates stress times the actual evolving strain rate.
//!
//! This is fixed-temperature, small-strain, linear viscoelasticity. Calibration,
//! finite-strain foam crushing, permanent set, damage and thermal evolution are
//! separate concerns. The returned interval sums are ordinary floating-point
//! bounds, not certified outward-rounded enclosures.

use super::{GeneralizedMaxwell, PronyState, ViscoError, prony_ramp};
use fs_la::eigen::jacobi_eigh;
use fs_math::det;

type Result<T> = core::result::Result<T, ViscoError>;
fn require(ok: bool, what: &'static str) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(ViscoError::Parameters { what })
    }
}
fn positive(value: f64, what: &'static str) -> Result<f64> {
    require(value.is_finite() && value > 0.0, what)?;
    Ok(value)
}
fn energy(law: &GeneralizedMaxwell, strain: f64, q: &[f64]) -> Result<f64> {
    let mut value = 0.5 * (law.e_inf * strain) * strain;
    for (&(e, _), &q) in law.terms.iter().zip(q) {
        value += 0.5 * (e * q) * q;
    }
    require(
        value.is_finite() && value >= 0.0,
        "unrepresentable stored Maxwell energy",
    )?;
    Ok(value)
}

/// Accepted stress ramp; energies are densities in J/m³.
#[derive(Debug, Clone, PartialEq)]
pub struct StressRampResponse {
    pub stress: f64,
    pub strain: f64,
    pub incremental_work_j_m3: f64,
    pub incremental_dissipated_j_m3: f64,
    pub stored_energy_j_m3: f64,
    /// Independent work minus change in storage minus independent dashpot heat.
    pub energy_balance_residual_j_m3: f64,
    /// Whole-ramp bounds from monotone modal endpoint intervals.
    pub strain_bounds: [f64; 2],
    /// Whole-ramp rate bounds [1/s]. A nontrivial instantaneous jump has
    /// `[-infinity,+infinity]`; a zero-duration no-op has `[0,0]`.
    pub strain_rate_bounds_s: [f64; 2],
}

/// Stress control sharing the original PronyState and its work/heat ledger.
/// Admission allows at most 64 strictly positive branches, or a pure elastic
/// law with no branches. Zero-weight branches remain valid in the original
/// GeneralizedMaxwell but are refused by this singularity-sensitive transform.
#[derive(Debug, Clone)]
pub struct StressControlledMaxwell {
    law: GeneralizedMaxwell,
    e0: f64,
    weights: Vec<f64>,
    rates: Vec<f64>,
    modes: Vec<f64>,
    forcing: Vec<f64>,
    branch_stress: Vec<f64>,
}
impl StressControlledMaxwell {
    pub fn new(law: &GeneralizedMaxwell) -> Result<Self> {
        positive(
            law.e_inf,
            "stress control needs positive finite equilibrium modulus",
        )?;
        let n = law.terms.len();
        require(n <= 64, "stress control admits at most 64 Maxwell branches")?;
        let mut e0 = law.e_inf;
        let mut weights = Vec::with_capacity(n);
        let mut branch_rates = Vec::with_capacity(n);
        for &(e, tau) in &law.terms {
            positive(e, "stress control requires strictly positive branch moduli")?;
            positive(
                tau,
                "stress control requires positive finite relaxation times",
            )?;
            e0 += e;
            weights.push(positive(
                e.sqrt() / tau.sqrt(),
                "unrepresentable Maxwell branch scaling",
            )?);
            branch_rates.push(positive(1.0 / tau, "unrepresentable Maxwell branch rate")?);
        }
        positive(e0, "unrepresentable instantaneous Maxwell modulus")?;
        if n == 0 {
            return Ok(Self {
                law: law.clone(),
                e0,
                weights,
                rates: vec![],
                modes: vec![],
                forcing: vec![],
                branch_stress: vec![],
            });
        }
        let scale = branch_rates.iter().copied().fold(0.0_f64, f64::max);
        let normalized: Vec<_> = branch_rates.iter().map(|r| r / scale).collect();
        let mut matrix = vec![0.0; n * n];
        for i in 0..n {
            positive(
                normalized[i],
                "Maxwell rate range exceeds the stress transform",
            )?;
            // Sum the other moduli instead of subtracting nearly equal values
            // in 1-E_i/E0 when one branch dominates.
            let other = law
                .terms
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .fold(law.e_inf, |sum, (_, &(e, _))| sum + e);
            matrix[i * n + i] = normalized[i] * (other / e0);
            for j in 0..i {
                let entry = -((law.terms[i].0 / e0).sqrt() * normalized[i].sqrt())
                    * ((law.terms[j].0 / e0).sqrt() * normalized[j].sqrt());
                matrix[i * n + j] = entry;
                matrix[j * n + i] = entry;
            }
        }
        let (eigenvalues, modes) = jacobi_eigh(&matrix, n);
        let tolerance = 4096.0 * f64::EPSILON * n as f64;
        for i in 0..n {
            positive(
                eigenvalues[i],
                "stress transform lost positive definiteness",
            )?;
            for row in 0..n {
                let actual = (0..n)
                    .map(|j| matrix[row * n + j] * modes[j * n + i])
                    .sum::<f64>();
                require(
                    (actual - eigenvalues[i] * modes[row * n + i]).abs() <= tolerance,
                    "stress transform eigensolver residual is unresolved",
                )?;
            }
            for j in 0..n {
                let inner = (0..n)
                    .map(|k| modes[k * n + i] * modes[k * n + j])
                    .sum::<f64>();
                require(
                    (inner - if i == j { 1.0 } else { 0.0 }).abs() <= tolerance,
                    "stress transform modes are not orthonormal",
                )?;
            }
        }
        let rates: Vec<_> = eigenvalues
            .iter()
            .map(|r| positive(r * scale, "unrepresentable stress relaxation rate"))
            .collect::<Result<_>>()?;
        let mut forcing = vec![0.0; n];
        let mut branch_stress = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                forcing[i] += modes[j * n + i] * (weights[j] / e0);
                branch_stress[i] += modes[j * n + i] * (law.terms[j].0 / weights[j]);
            }
            require(
                forcing[i].is_finite() && branch_stress[i].is_finite(),
                "unrepresentable stress transform coefficients",
            )?;
        }
        Ok(Self {
            law: law.clone(),
            e0,
            weights,
            rates,
            modes,
            forcing,
            branch_stress,
        })
    }

    pub fn law(&self) -> &GeneralizedMaxwell {
        &self.law
    }

    /// Advance a prescribed linear stress ramp, including holds, reversal and
    /// zero-load recovery. Any failure leaves every state field unchanged.
    /// `dt=0` is an elastic jump with zero dashpot heat. Positive-duration
    /// ramps refuse unrepresentable rates and unresolved energy balance.
    pub fn step(
        &self,
        state: &mut PronyState,
        stress_new: f64,
        dt: f64,
    ) -> Result<StressRampResponse> {
        let n = self.weights.len();
        require(
            dt.is_finite() && dt >= 0.0 && stress_new.is_finite(),
            "stress ramp needs finite stress and nonnegative finite duration",
        )?;
        require(
            state.q.len() == n
                && state.q.iter().all(|q| q.is_finite())
                && state.strain.is_finite()
                && state.stress.is_finite()
                && state.work.is_finite()
                && state.dissipated.is_finite()
                && state.dissipated >= 0.0,
            "stress ramp requires a finite matching Prony state",
        )?;
        let mut calculated = self.law.e_inf * state.strain;
        let mut magnitude = calculated.abs() + state.stress.abs();
        for (&(e, _), &q) in self.law.terms.iter().zip(&state.q) {
            calculated += e * q;
            magnitude += (e * q).abs();
        }
        require(
            calculated.is_finite()
                && magnitude.is_finite()
                && (calculated - state.stress).abs()
                    <= 128.0 * f64::EPSILON * (n + 1) as f64 * magnitude.max(f64::MIN_POSITIVE),
            "Prony state stress does not match its internal variables",
        )?;
        let before = energy(&self.law, state.strain, &state.q)?;
        let delta_stress = stress_new - state.stress;
        require(delta_stress.is_finite(), "unrepresentable stress increment")?;
        let mut q = state.q.clone();
        let mut work = (0.5 * state.stress + 0.5 * stress_new) * (delta_stress / self.e0);
        let mut heat = 0.0;
        let (strain, strain_bounds, rate_bounds);
        if dt == 0.0 || n == 0 {
            let delta_strain = delta_stress / self.e0;
            strain = state.strain + delta_strain;
            for q in &mut q {
                *q += delta_strain;
            }
            strain_bounds = [state.strain.min(strain), state.strain.max(strain)];
            rate_bounds = if dt > 0.0 {
                [delta_strain / dt; 2]
            } else if delta_stress == 0.0 {
                [0.0; 2]
            } else {
                [f64::NEG_INFINITY, f64::INFINITY]
            };
        } else {
            let mut next = vec![0.0; n];
            let mut branch_bounds = [0.0; 2];
            let mut rates = [(delta_stress / self.e0) / dt; 2];
            for i in 0..n {
                let initial = (0..n)
                    .map(|j| self.modes[j * n + i] * (self.weights[j] * state.q[j]))
                    .sum::<f64>();
                let delta = self.forcing[i] * delta_stress;
                let (end, mean, loss) = prony_ramp(initial, delta, self.rates[i] * dt);
                require(
                    end.is_finite() && mean.is_finite() && loss.is_finite() && loss >= 0.0,
                    "unrepresentable stress-controlled modal history",
                )?;
                next[i] = end;
                heat += loss / self.rates[i];
                let (zeroth, first) = modal_integrals(initial, delta, self.rates[i], dt, mean);
                work += (self.forcing[i] * zeroth) * state.stress
                    + (self.forcing[i] * first) * delta_stress;
                let a = self.branch_stress[i] * initial;
                let b = self.branch_stress[i] * end;
                branch_bounds[0] += a.min(b);
                branch_bounds[1] += a.max(b);
                let a = self.forcing[i] * initial;
                let b = self.forcing[i] * end;
                rates[0] += a.min(b);
                rates[1] += a.max(b);
            }
            for (j, q) in q.iter_mut().enumerate() {
                *q = (0..n).map(|i| self.modes[j * n + i] * next[i]).sum::<f64>() / self.weights[j];
            }
            let branch_sum = self
                .law
                .terms
                .iter()
                .zip(&q)
                .map(|(&(e, _), &q)| e * q)
                .sum::<f64>();
            strain = (stress_new - branch_sum) / self.law.e_inf;
            strain_bounds = [
                ((state.stress.min(stress_new) - branch_bounds[1]) / self.law.e_inf)
                    .min(state.strain)
                    .min(strain),
                ((state.stress.max(stress_new) - branch_bounds[0]) / self.law.e_inf)
                    .max(state.strain)
                    .max(strain),
            ];
            rate_bounds = rates;
        }
        require(
            strain.is_finite()
                && q.iter().all(|q| q.is_finite())
                && strain_bounds.iter().all(|x| x.is_finite())
                && (dt == 0.0 || rate_bounds.iter().all(|x| x.is_finite()))
                && work.is_finite()
                && heat.is_finite()
                && heat >= 0.0,
            "unrepresentable stress ramp response",
        )?;
        let stored = energy(&self.law, strain, &q)?;
        let residual = work - (stored - before) - heat;
        let balance_scale = work
            .abs()
            .max(before)
            .max(stored)
            .max(heat)
            .max(f64::MIN_POSITIVE);
        require(
            residual.is_finite() && residual.abs() <= 1e-8 * balance_scale,
            "stress ramp energy balance is numerically unresolved",
        )?;
        let cumulative_work = state.work + work;
        let cumulative_heat = state.dissipated + heat;
        require(
            cumulative_work.is_finite() && cumulative_heat.is_finite(),
            "unrepresentable cumulative Maxwell work or heat",
        )?;
        let response = StressRampResponse {
            stress: stress_new,
            strain,
            incremental_work_j_m3: work,
            incremental_dissipated_j_m3: heat,
            stored_energy_j_m3: stored,
            energy_balance_residual_j_m3: residual,
            strain_bounds,
            strain_rate_bounds_s: rate_bounds,
        };
        *state = PronyState {
            q,
            strain,
            stress: stress_new,
            work: cumulative_work,
            dissipated: cumulative_heat,
        };
        Ok(response)
    }
}

// Integrals of y(t) and (t/dt)y(t), respectively. Computing them directly for
// large lambda*dt retains finite creep work even when the mean underflows.
fn modal_integrals(y0: f64, delta: f64, rate: f64, dt: f64, mean: f64) -> (f64, f64) {
    let r = rate * dt;
    if r < 0.1 {
        // d=integral_0^1 u exp(-r*u)du; e=(1/2-d)/r.
        let d = 0.5
            + r * (-1.0 / 3.0
                + r * (1.0 / 8.0
                    + r * (-1.0 / 30.0
                        + r * (1.0 / 144.0
                            + r * (-1.0 / 840.0
                                + r * (1.0 / 5760.0
                                    + r * (-1.0 / 45360.0
                                        + r * (1.0 / 403200.0
                                            + r * (-1.0 / 3991680.0 + r / 43545600.0)))))))));
        let e = 1.0 / 3.0
            + r * (-1.0 / 8.0
                + r * (1.0 / 30.0
                    + r * (-1.0 / 144.0
                        + r * (1.0 / 840.0
                            + r * (-1.0 / 5760.0
                                + r * (1.0 / 45360.0
                                    + r * (-1.0 / 403200.0
                                        + r * (1.0 / 3991680.0 - r / 43545600.0))))))));
        (dt * mean, dt * (d * y0 + e * delta))
    } else {
        let a = det::exp(-r);
        let s = -det::expm1(-r);
        let b = s / r;
        let d = (b - a) / r;
        (
            (s * y0 + (1.0 - b) * delta) / rate,
            ((b - a) * y0 + (0.5 - d) * delta) / rate,
        )
    }
}
