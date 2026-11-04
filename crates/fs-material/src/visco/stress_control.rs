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

use super::{GeneralizedMaxwell, PronyRamp, PronyState, ViscoError};
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
    /// Whole-ramp bounds from integrated monotone modal response terms.
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
        for i in 0..n {
            for j in 0..n {
                forcing[i] += modes[j * n + i] * (weights[j] / e0);
            }
            require(
                forcing[i].is_finite(),
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
        if dt == 0.0 && delta_stress == 0.0 {
            // Entry coherence admits floating-point roundoff, but a no-op
            // must not use that allowance to alter an accepted history.
            return Ok(StressRampResponse {
                stress: state.stress,
                strain: state.strain,
                incremental_work_j_m3: 0.0,
                incremental_dissipated_j_m3: 0.0,
                stored_energy_j_m3: before,
                energy_balance_residual_j_m3: 0.0,
                strain_bounds: [state.strain; 2],
                strain_rate_bounds_s: [0.0; 2],
            });
        }
        let mut q = state.q.clone();
        let mut work = (0.5 * state.stress + 0.5 * stress_new) * (delta_stress / self.e0);
        let mut heat = 0.0;
        let mut reconstruction_scale = stress_new.abs();
        let (mut strain, mut strain_bounds, rate_bounds);
        if dt == 0.0 || n == 0 {
            let delta_strain = delta_stress / self.e0;
            strain = state.strain + delta_strain;
            reconstruction_scale +=
                self.law.e_inf * state.strain.abs() + self.law.e_inf * delta_strain.abs();
            for (&(e, _), q) in self.law.terms.iter().zip(&mut q) {
                reconstruction_scale += e * q.abs() + e * delta_strain.abs();
                *q += delta_strain;
            }
            strain_bounds = [state.strain.min(strain), state.strain.max(strain)];
            rate_bounds = if dt > 0.0 {
                [delta_strain / dt; 2]
            } else {
                [f64::NEG_INFINITY, f64::INFINITY]
            };
        } else {
            let mut next = vec![(0.0, 0.0); n];
            let elastic_increment = delta_stress / self.e0;
            let mut increment = elastic_increment;
            let mut increment_bounds = [elastic_increment.min(0.0), elastic_increment.max(0.0)];
            let mut rates = [elastic_increment / dt; 2];
            for i in 0..n {
                let mut initial = 0.0;
                let mut initial_magnitude = 0.0;
                for j in 0..n {
                    let term = self.modes[j * n + i] * (self.weights[j] * state.q[j]);
                    initial += term;
                    initial_magnitude += term.abs();
                }
                let delta = self.forcing[i] * delta_stress;
                let ramp = PronyRamp::new(self.rates[i] * dt);
                let (end, mean, loss) = ramp.apply(initial, delta);
                require(
                    end.is_finite() && mean.is_finite() && loss.is_finite() && loss >= 0.0,
                    "unrepresentable stress-controlled modal history",
                )?;
                next[i] = (
                    end,
                    ramp.decay * initial_magnitude + (ramp.ramp * delta).abs(),
                );
                heat += loss / self.rates[i];
                let (zeroth, first, parts) =
                    modal_integrals(initial, delta, self.rates[i], dt, mean, &ramp);
                increment += self.forcing[i] * zeroth;
                work += (self.forcing[i] * zeroth) * state.stress
                    + (self.forcing[i] * first) * delta_stress;
                // Each integrated initial/forcing term is monotone from zero
                // to its endpoint. Summing their intervals also bounds paths
                // that reverse direction, without dividing by E_inf or using
                // dt times a potentially enormous initial relaxation rate.
                for part in parts {
                    let contribution = self.forcing[i] * part;
                    increment_bounds[0] += contribution.min(0.0);
                    increment_bounds[1] += contribution.max(0.0);
                }
                let a = self.forcing[i] * initial;
                let b = self.forcing[i] * end;
                rates[0] += a.min(b);
                rates[1] += a.max(b);
            }
            for (j, q) in q.iter_mut().enumerate() {
                let mut value = 0.0;
                let mut magnitude = 0.0;
                for (i, &(end, scale)) in next.iter().enumerate() {
                    let mode = self.modes[j * n + i];
                    value += mode * end;
                    magnitude += mode.abs() * scale;
                }
                *q = value / self.weights[j];
                reconstruction_scale += (magnitude / self.weights[j]) * self.law.terms[j].0;
            }
            // Integrate the same strain rate used by external work. Recovering
            // strain from (stress - branch stress)/E_inf loses creep accuracy
            // when the equilibrium spring is weak compared with the branches.
            strain = state.strain + increment;
            reconstruction_scale += self.law.e_inf * state.strain.abs()
                + self.law.e_inf * increment_bounds[0].abs()
                + self.law.e_inf * increment_bounds[1].abs();
            strain_bounds = [
                (state.strain + increment_bounds[0])
                    .min(state.strain)
                    .min(strain),
                (state.strain + increment_bounds[1])
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
        // Roundoff in the modal transforms can leave a tiny constitutive
        // residual, especially when unloading cancels the branch history.
        // Close only that residual through the instantaneous modulus, never
        // through the weak equilibrium spring. The allowance uses the actual
        // reconstructed terms before cancellation, not just the target stress.
        // This is numerical closure, not an additional physical stress jump.
        let reconstructed = self.law.e_inf * strain
            + self
                .law
                .terms
                .iter()
                .zip(&q)
                .map(|(&(e, _), &q)| e * q)
                .sum::<f64>();
        let stress_residual = stress_new - reconstructed;
        let relative_roundoff = 128.0 * f64::EPSILON * (n + 1) as f64;
        require(
            reconstruction_scale.is_finite()
                && stress_residual.is_finite()
                && stress_residual.abs()
                    <= relative_roundoff * reconstruction_scale.max(f64::MIN_POSITIVE),
            "stress ramp constitutive reconstruction is numerically unresolved",
        )?;
        let correction = stress_residual / self.e0;
        if correction != 0.0 {
            strain += correction;
            for q in &mut q {
                *q += correction;
            }
            strain_bounds[0] = strain_bounds[0].min(strain);
            strain_bounds[1] = strain_bounds[1].max(strain);
        }
        let mut corrected_stress = self.law.e_inf * strain;
        let mut corrected_magnitude = corrected_stress.abs() + stress_new.abs();
        for (&(e, _), &q) in self.law.terms.iter().zip(&q) {
            corrected_stress += e * q;
            corrected_magnitude += (e * q).abs();
        }
        require(
            strain.is_finite()
                && q.iter().all(|q| q.is_finite())
                && corrected_stress.is_finite()
                && corrected_magnitude.is_finite()
                && (corrected_stress - stress_new).abs()
                    <= relative_roundoff * corrected_magnitude.max(f64::MIN_POSITIVE),
            "stress ramp numerical closure is unresolved",
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

// Integrals of y(t) and (t/dt)y(t), plus the separate initial and forcing
// contributions to the first integral. Both contributions are monotone over
// the ramp, so their endpoint intervals bound strain even through reversal.
// Computing them directly for large lambda*dt retains finite creep strain and
// work even when the mean underflows. Reuse the original Prony coefficients.
fn modal_integrals(
    y0: f64,
    delta: f64,
    rate: f64,
    dt: f64,
    mean: f64,
    ramp: &PronyRamp,
) -> (f64, f64, [f64; 2]) {
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
        (
            dt * mean,
            dt * (d * y0 + e * delta),
            [dt * (ramp.ramp * y0), dt * (ramp.mean_ramp * delta)],
        )
    } else {
        let a = ramp.decay;
        let s = -det::expm1(-r);
        let b = ramp.ramp;
        let d = (b - a) / r;
        (
            (s * y0 + (1.0 - b) * delta) / rate,
            ((b - a) * y0 + (0.5 - d) * delta) / rate,
            [(s * y0) / rate, ((1.0 - b) * delta) / rate],
        )
    }
}
