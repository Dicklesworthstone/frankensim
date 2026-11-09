//! Intrinsic series-creep felt: the original WoolFelt envelope and passive
//! Kelvin histories in strain/Pa units. No geometry, material defaults or
//! instrument identity is introduced here.
//!
//! Each prepared step holds the common creep stress over its interval. The
//! nonlinear spring follows its original quasistatic path to the trial elastic
//! strain; each Kelvin arm then has an exact held-stress trajectory. Under a
//! strain history this is an implicit first-order discretization and requires
//! refinement. Work integrates that declared path; heat is evaluated directly
//! from the arm loss law and irreversible compaction, not an energy remainder.
use super::{Uniaxial, WoolFelt, WoolFeltResponse, WoolFeltState};
use crate::MaterialError;
use fs_math::{dd::Dd, det};

mod creep;
mod prony;
use creep::PreparedBranch;
pub use creep::{CreepAdvance, CreepSpectrum, PreparedCreep};

/// The original source realization retains at most eight active material arms.
pub const MAX_FELT_RATE_BRANCHES: usize = 8;

/// A series Kelvin arm per unit strain: eta*z_dot + E*z = stress.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CreepBranch {
    pub modulus_pa: f64,
    pub viscosity_pa_s: f64,
}

/// A nonlinear unilateral spring in series with explicit passive creep arms.
/// Parameters are immutable after admission; preparing a clock caches all
/// branch exponentials. The supplied envelope retains its original meaning.
#[derive(Clone, Debug)]
pub struct RateFelt {
    envelope: WoolFelt,
    creep: Vec<CreepBranch>,
}

/// Complete physical checkpoint. Energy densities are J/m^3 (equivalently Pa).
/// Inactive entries of creep_strain must remain exactly zero.
#[derive(Clone, Debug, PartialEq)]
pub struct RateFeltState {
    pub hysteresis: WoolFeltState,
    pub strain: f64,
    pub elastic_strain: f64,
    pub creep_strain: [f64; MAX_FELT_RATE_BRANCHES],
    pub work: f64,
    pub stored: f64,
    pub dissipated: f64,
}

/// Trial response; commit its state only after the enclosing solve accepts it.
#[derive(Clone, Debug, PartialEq)]
pub struct RateFeltResponse {
    pub state: RateFeltState,
    pub stress_pa: f64,
    /// Consistent derivative of this strain-driven trial, history held fixed.
    pub tangent_pa: f64,
    pub work_density_j_m3: f64,
    pub stored_energy_change_j_m3: f64,
    pub creep_dissipation_j_m3: f64,
    pub crush_dissipation_j_m3: f64,
    pub balance_residual_j_m3: f64,
}

/// Prepared original law and fixed nonnegative interval. Warm trials allocate
/// no history buffers and evaluate no branch exponentials. Nonlinear spring
/// evaluations still use the original deterministic power law.
#[derive(Clone, Debug)]
pub struct PreparedRateFelt {
    law: RateFelt,
    duration: f64,
    branches: [PreparedBranch; MAX_FELT_RATE_BRANCHES],
    compliance: f64,
}

fn bad(what: &str) -> MaterialError {
    MaterialError::State {
        what: format!("rate felt: {what}"),
    }
}
fn finite(value: f64, what: &str) -> Result<f64, MaterialError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad(what))
    }
}
fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}

// Compensate 1+x rounding while retaining the deterministic logarithm owner.
fn log1p_positive(x: f64) -> f64 {
    if x.is_infinite() {
        return x;
    }
    let u = 1.0 + x;
    if u == 1.0 {
        x
    } else {
        det::ln(u) * (x / (u - 1.0))
    }
}

// Integral difference for the ORIGINAL nonnegative power law. A near-peak
// increment is retained instead of subtracting two large stored energies.
fn power_difference(high: f64, low: f64, power: f64) -> f64 {
    if high == low {
        return 0.0;
    }
    let value = det::pow(high, power);
    if low == 0.0 {
        value
    } else {
        value * -det::expm1(-power * log1p_positive((high - low) / low))
    }
}

impl RateFelt {
    pub fn new(envelope: WoolFelt, creep: Vec<CreepBranch>) -> Result<Self, MaterialError> {
        envelope.response(0.0, &envelope.initial_state())?;
        creep::validate_branches(&creep)?;
        Ok(Self { envelope, creep })
    }

    pub fn envelope(&self) -> &WoolFelt {
        &self.envelope
    }
    pub fn creep(&self) -> &[CreepBranch] {
        &self.creep
    }

    pub fn initial_state(&self) -> RateFeltState {
        RateFeltState {
            hysteresis: self.envelope.initial_state(),
            strain: 0.0,
            elastic_strain: 0.0,
            creep_strain: [0.0; MAX_FELT_RATE_BRANCHES],
            work: 0.0,
            stored: 0.0,
            dissipated: 0.0,
        }
    }

    pub fn prepare(&self, duration_s: f64) -> Result<PreparedRateFelt, MaterialError> {
        let (branches, compliance) = creep::prepare_branches(&self.creep, duration_s)?;
        Ok(PreparedRateFelt {
            law: self.clone(),
            duration: duration_s,
            branches,
            compliance,
        })
    }

    pub fn step_strain_driven(
        &self,
        state: &RateFeltState,
        strain: f64,
        duration_s: f64,
    ) -> Result<RateFeltResponse, MaterialError> {
        self.prepare(duration_s)?.step_strain_driven(state, strain)
    }

    pub fn step_stress_driven(
        &self,
        state: &RateFeltState,
        stress_pa: f64,
        duration_s: f64,
    ) -> Result<RateFeltResponse, MaterialError> {
        self.prepare(duration_s)?
            .step_stress_driven(state, stress_pa)
    }

    /// Original spring recoverable energy density at a committed history.
    fn spring_storage(&self, strain: f64, history: &WoolFeltState) -> Result<f64, MaterialError> {
        if history.eps_max == 0.0 {
            return Ok(0.0);
        }
        let residual = self.envelope.eps_residual(history);
        if strain <= residual {
            return Ok(0.0);
        }
        let span = history.eps_max - residual;
        finite(
            history.sig_max * span / (self.envelope.q + 1.0)
                * det::pow((strain - residual) / span, self.envelope.q + 1.0),
            "unrepresentable recovered spring energy",
        )
    }

    /// Original stress integral and independent positive compaction loss.
    /// An increased peak loses the explicit fraction
    /// ((q-p)+c*(p+1))/(q+1) of the added virgin-envelope work.
    fn spring_step(
        &self,
        before: f64,
        after: f64,
        history: &WoolFeltState,
    ) -> Result<(Dd, Dd, Dd), MaterialError> {
        let l = &self.envelope;
        let residual = l.eps_residual(history);
        let span = history.eps_max - residual;
        let return_change = |a: f64, b: f64| -> Result<Dd, MaterialError> {
            if history.eps_max == 0.0 || a == b {
                return Ok(Dd::ZERO);
            }
            let x = ((a - residual) / span).max(0.0);
            let y = ((b - residual) / span).max(0.0);
            let delta = power_difference(x.max(y), x.min(y), l.q + 1.0);
            let value = dd(history.sig_max) * dd(span) / dd(l.q + 1.0) * dd(delta);
            let value = if y >= x { value } else { -value };
            finite(value.to_f64(), "original spring return work overflow")?;
            Ok(value)
        };
        if after <= history.eps_max {
            let work = return_change(before, after)?;
            return Ok((work, work, Dd::ZERO));
        }
        let to_peak = return_change(before, history.eps_max)?;
        let added = dd(l.f_ref) * dd(l.eps_ref) / dd(l.p + 1.0)
            * dd(power_difference(
                after / l.eps_ref,
                history.eps_max / l.eps_ref,
                l.p + 1.0,
            ));
        let fraction = (dd(l.q - l.p) + dd(l.crush_fraction) * dd(l.p + 1.0)) / dd(l.q + 1.0);
        let crush = added * fraction;
        let work = to_peak + added;
        let storage = to_peak + added * (dd(1.0) - fraction);
        finite(work.to_f64(), "original spring loading work overflow")?;
        finite(crush.to_f64(), "original compaction heat overflow")?;
        Ok((work, storage, crush))
    }

    fn admit_state(&self, state: &RateFeltState) -> Result<WoolFeltResponse, MaterialError> {
        if !state.strain.is_finite()
            || state.strain > self.envelope.eps_densify
            || !state.elastic_strain.is_finite()
            || state.elastic_strain > state.hysteresis.eps_max
            || !state.work.is_finite()
            || !state.stored.is_finite()
            || state.stored < 0.0
            || !state.dissipated.is_finite()
            || state.dissipated < 0.0
            || state
                .creep_strain
                .iter()
                .any(|z| !z.is_finite() || *z < 0.0)
            || state.creep_strain[self.creep.len()..]
                .iter()
                .any(|z| *z != 0.0)
        {
            return Err(bad("inadmissible physical history or energy ledger"));
        }
        let response = self
            .envelope
            .response(state.elastic_strain, &state.hysteresis)?;
        let mut total = dd(state.elastic_strain);
        let mut stored = dd(self.spring_storage(state.elastic_strain, &state.hysteresis)?);
        for (branch, &z) in self.creep.iter().zip(&state.creep_strain) {
            total = total + dd(z);
            stored = stored + dd(0.5) * dd(branch.modulus_pa) * dd(z) * dd(z);
        }
        let strain_scale = state
            .strain
            .abs()
            .max(state.elastic_strain.abs())
            .max(state.creep_strain.iter().copied().fold(0.0, f64::max));
        if (total - dd(state.strain)).to_f64().abs() > 128.0 * f64::EPSILON * strain_scale {
            return Err(bad(
                "total and elastic/creep strains do not describe the same state",
            ));
        }
        let stored = finite(stored.to_f64(), "stored energy overflow")?;
        let scale = stored
            .abs()
            .max(state.stored.abs())
            .max(state.work.abs())
            .max(state.dissipated);
        if (stored - state.stored).abs() > 1e-12 * scale
            || (dd(state.work) - dd(state.stored) - dd(state.dissipated))
                .to_f64()
                .abs()
                > 1e-12 * scale
        {
            return Err(bad("physical history and original energy ledger disagree"));
        }
        Ok(response)
    }
}

impl PreparedRateFelt {
    pub fn duration_s(&self) -> f64 {
        self.duration
    }
    pub fn creep_compliance_pa_inverse(&self) -> f64 {
        self.compliance
    }

    /// Decayed free strain of the original history, for implicit contact solves.
    pub fn free_creep_strain(&self, state: &RateFeltState) -> Result<f64, MaterialError> {
        self.law.admit_state(state)?;
        self.free_strain(state)
    }
    fn free_strain(&self, state: &RateFeltState) -> Result<f64, MaterialError> {
        creep::free_strain(&self.branches[..self.law.creep.len()], &state.creep_strain)
    }

    /// Exact held nonnegative stress for the creep arms. A zero-stress command
    /// chooses the original spring's residual strain; an open gap has no unique
    /// strain under force control. Recovery continues after stress removal.
    pub fn step_stress_driven(
        &self,
        state: &RateFeltState,
        stress: f64,
    ) -> Result<RateFeltResponse, MaterialError> {
        let previous = self.law.admit_state(state)?;
        if !stress.is_finite() || stress < 0.0 {
            return Err(bad("finite compressive stress required"));
        }
        let l = &self.law.envelope;
        let h = &state.hysteresis;
        let elastic = if stress == 0.0 {
            l.eps_residual(h)
        } else if stress == previous.stress_pa {
            state.elastic_strain
        } else if stress == h.sig_max {
            h.eps_max
        } else if stress > h.sig_max {
            l.eps_ref * det::pow(stress / l.f_ref, 1.0 / l.p)
        } else {
            let residual = l.eps_residual(h);
            residual + (h.eps_max - residual) * det::pow(stress / h.sig_max, 1.0 / l.q)
        };
        let trial = l.response(elastic, h)?;
        let scale = stress.abs().max(trial.stress_pa.abs());
        if (trial.stress_pa - stress).abs() > 256.0 * f64::EPSILON * scale {
            return Err(bad("original felt stress inverse unresolved"));
        }
        self.finish(state, elastic, stress, trial.tangent_pa, None)
    }

    /// Implicit endpoint response for a total engineering strain. The monotone
    /// scalar solve retains every Kelvin arm on the same held-stress path.
    pub fn step_strain_driven(
        &self,
        state: &RateFeltState,
        strain: f64,
    ) -> Result<RateFeltResponse, MaterialError> {
        self.law.admit_state(state)?;
        if !strain.is_finite() || strain > self.law.envelope.eps_densify {
            return Err(bad(
                "finite total strain within densification domain required",
            ));
        }
        let free = self.free_strain(state)?;
        let available = finite(
            (dd(strain) - dd(free)).to_f64(),
            "available spring strain overflow",
        )?;
        let h = &state.hysteresis;
        let l = &self.law.envelope;
        let upper = l.response(available, h)?;
        if self.compliance == 0.0 || upper.stress_pa == 0.0 {
            return self.finish(
                state,
                available,
                upper.stress_pa,
                upper.tangent_pa,
                Some(strain),
            );
        }
        let mut lo = l.eps_residual(h).min(available);
        let mut hi = available;
        let mut elastic = hi;
        let mut resolved = None;
        for _ in 0..96 {
            let response = l.response(elastic, h)?;
            let residual = finite(
                (dd(elastic) + dd(self.compliance) * dd(response.stress_pa) - dd(available))
                    .to_f64(),
                "strain-driven creep residual overflow",
            )?;
            let scale = elastic
                .abs()
                .max(available.abs())
                .max((self.compliance * response.stress_pa).abs());
            if residual.abs() <= 32.0 * f64::EPSILON * scale {
                resolved = Some(response);
                break;
            }
            if residual > 0.0 {
                hi = elastic;
            } else {
                lo = elastic;
            }
            let slope = 1.0 + self.compliance * response.tangent_pa;
            let newton = elastic - residual / slope;
            let next = if newton.is_finite() && newton > lo && newton < hi {
                newton
            } else {
                0.5 * lo + 0.5 * hi
            };
            if next == elastic {
                return Err(bad("strain-driven creep solve stalled"));
            }
            elastic = next;
        }
        let response = resolved.ok_or_else(|| bad("strain-driven creep solve did not converge"))?;
        self.finish(
            state,
            elastic,
            response.stress_pa,
            response.tangent_pa,
            Some(strain),
        )
    }

    fn finish(
        &self,
        state: &RateFeltState,
        elastic: f64,
        stress: f64,
        spring_tangent: f64,
        requested_strain: Option<f64>,
    ) -> Result<RateFeltResponse, MaterialError> {
        let trial = self.law.envelope.response(elastic, &state.hysteresis)?;
        let after_spring = self.law.spring_storage(elastic, &trial.state)?;
        let (spring_work, spring_change, crush) =
            self.law
                .spring_step(state.elastic_strain, elastic, &state.hysteresis)?;
        if crush.to_f64() < 0.0 {
            return Err(bad("negative original compaction loss"));
        }
        let mut next = state.clone();
        next.hysteresis = trial.state;
        next.elastic_strain = elastic;
        let mut total = dd(elastic);
        let mut work = spring_work;
        let mut heat = Dd::ZERO;
        let mut storage_change = spring_change;
        let mut stored = dd(after_spring);
        for (i, b) in self.branches[..self.law.creep.len()].iter().enumerate() {
            let step = b.advance(state.creep_strain[i], stress)?;
            let after = step.after;
            let change = step.change;
            let actual = step.actual;
            let scale = after
                .abs()
                .max(state.creep_strain[i])
                .max(change.to_f64().abs());
            if (change.to_f64() != 0.0 && actual.to_f64() == 0.0)
                || (actual - change).to_f64().abs() > 8.0 * f64::EPSILON * scale
            {
                return Err(bad("unresolved physical creep increment"));
            }
            next.creep_strain[i] = after;
            total = total + dd(after);
            let lost = step.heat;
            let branch_work = step.work;
            let branch_change = step.stored_change;
            let branch_scale = branch_work
                .to_f64()
                .abs()
                .max(branch_change.to_f64().abs())
                .max(lost.to_f64().abs());
            if lost.to_f64() < 0.0
                || !lost.to_f64().is_finite()
                || (branch_work - branch_change - lost).to_f64().abs() > 1e-12 * branch_scale
            {
                return Err(bad("physical Kelvin work and heat do not close"));
            }
            work = work + branch_work;
            heat = heat + lost;
            storage_change = storage_change + branch_change;
            stored = stored + step.stored;
        }
        let represented = finite(total.to_f64(), "total felt strain overflow")?;
        next.strain = requested_strain.unwrap_or(represented);
        let strain_scale = represented.abs().max(next.strain.abs()).max(elastic.abs());
        if next.strain > self.law.envelope.eps_densify
            || (total - dd(next.strain)).to_f64().abs() > 128.0 * f64::EPSILON * strain_scale
        {
            return Err(bad(
                "accepted total strain unresolved or beyond densification",
            ));
        }
        let loss = heat + crush;
        next.work = finite((dd(state.work) + work).to_f64(), "cumulative work overflow")?;
        next.stored = finite(stored.to_f64(), "stored energy overflow")?;
        next.dissipated = finite(
            (dd(state.dissipated) + loss).to_f64(),
            "cumulative loss overflow",
        )?;
        let residual = work - storage_change - loss;
        let scale = work
            .to_f64()
            .abs()
            .max(storage_change.to_f64().abs())
            .max(loss.to_f64().abs());
        if residual.to_f64().abs() > 1e-12 * scale {
            return Err(bad("original felt work balance unresolved"));
        }
        self.law.admit_state(&next)?;
        let tangent = if self.compliance == 0.0 {
            spring_tangent
        } else if (self.compliance * spring_tangent).is_finite() {
            spring_tangent / (1.0 + self.compliance * spring_tangent)
        } else {
            1.0 / (self.compliance + 1.0 / spring_tangent)
        };
        Ok(RateFeltResponse {
            state: next,
            stress_pa: stress,
            tangent_pa: finite(tangent, "consistent tangent overflow")?,
            work_density_j_m3: finite(work.to_f64(), "work overflow")?,
            stored_energy_change_j_m3: finite(storage_change.to_f64(), "energy change overflow")?,
            creep_dissipation_j_m3: finite(heat.to_f64(), "creep heat overflow")?,
            crush_dissipation_j_m3: crush.to_f64(),
            balance_residual_j_m3: finite(residual.to_f64(), "work residual overflow")?,
        })
    }
}
