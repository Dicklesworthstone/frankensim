//! The original aqueous solve inside one monotone implicit formula extent.

use super::*;

impl MineralKineticsV1 {
    /// Advance only the selected mineral at fixed original source temperature
    /// and activity basis. Charge and acid/base families equilibrate at every
    /// candidate extent. Refusal or cancellation never changes the input.
    /// Prepared identities/buffers are reused for the entire window; the full
    /// dissolved-species report is materialized only for its final state.
    pub fn advance(
        &self,
        initial: &MineralKineticStateV1,
        duration: Time,
        controls: MineralKineticControlsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MineralKineticOutcomeV1> {
        poll(&mut cancelled)?;
        self.admit_owner(initial)?;
        admit(duration, controls)?;
        let mut workspace = Workspace::new(self, &initial.stocks, initial.basis, controls.aqueous, &mut cancelled)?;
        if duration.value() == 0.0 || self.rate().rate_mol_s() == 0.0 {
            poll(&mut cancelled)?;
            return Ok(MineralKineticOutcomeV1 {
                final_state: initial.clone(), duration, transfer: Amount::new(0.0),
                accepted: 0, rejected: 0, counts: Counts::default(), maximum_error: 0.0,
                blocked: false, exhausted: false,
            });
        }
        let mut current = initial.phase;
        let mut elapsed = 0.0;
        let mut interval = controls.initial_step.value();
        let (mut accepted, mut rejected) = (0, 0);
        let mut maximum_error: f64 = 0.0;
        let (mut blocked, mut exhausted) = (false, false);
        while elapsed < duration.value() {
            poll(&mut cancelled)?;
            if accepted + rejected == controls.maximum_attempts { return Err(AcidBaseError::Budget); }
            let remaining = duration.value() - elapsed;
            let h = interval.min(remaining);
            let end = if h == remaining { duration.value() } else { elapsed + h };
            let h = end - elapsed;
            if h <= 0.0 || 0.5 * h == 0.0 { return Err(AcidBaseError::Arithmetic("mineral time increment loses range")); }
            let full = workspace.stage(current, h, controls, &mut cancelled)?;
            let middle = workspace.stage(current, 0.5 * h, controls, &mut cancelled)?;
            let fine = workspace.stage(middle.phase, h - 0.5 * h, controls, &mut cancelled)?;
            let error = temporal_error(full.phase, fine.phase, controls);
            if error > 1.0 {
                rejected += 1;
                let next = 0.5 * h;
                if next < controls.minimum_step.value() || next <= 0.0 || next == h {
                    return Err(AcidBaseError::Arithmetic("mineral temporal tolerance cannot be met at minimum step"));
                }
                interval = next;
                continue;
            }
            current = fine.phase;
            blocked |= middle.blocked || fine.blocked;
            exhausted |= middle.exhausted || fine.exhausted;
            elapsed = end;
            accepted += 1;
            maximum_error = maximum_error.max(error);
            interval = if error <= 0.25 { (h + h).min(controls.maximum_step.value()) } else { h };
        }
        let transfer = if initial.phase.dissolved.max(current.dissolved) <= initial.phase.solid.max(current.solid) {
            current.dissolved - initial.phase.dissolved
        } else { initial.phase.solid - current.solid };
        let (final_state, counts) = workspace.finish(
            self, Arc::clone(&initial.stocks), current, initial.basis, controls.aqueous, &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(MineralKineticOutcomeV1 {
            final_state, duration, transfer: Amount::new(transfer), accepted, rejected, counts,
            maximum_error, blocked, exhausted,
        })
    }
}

struct Stage { phase: Phase, blocked: bool, exhausted: bool }
impl Workspace<'_> {
    fn stage(&mut self, initial: Phase, duration: f64, controls: MineralKineticControlsV1, cancelled: &mut impl FnMut() -> bool) -> Result<Stage> {
        let start = self.evaluate(initial, cancelled)?;
        let growing = start.log_omega.is_some_and(|log| log > 0.0);
        if initial.solid == 0.0 && self.model.rate().growth_policy() == MineralGrowthPolicyV1::ExistingSolidOnly {
            return Ok(Stage { phase: initial, blocked: growing, exhausted: false });
        }
        let direction = if growing { -1.0 } else { 1.0 };
        let donor = if growing { initial.dissolved } else { initial.solid };
        if donor == 0.0 { return Ok(Stage { phase: initial, blocked: false, exhausted: false }); }
        let log_time_rate = det::ln(duration) + det::ln(self.model.rate().rate_mol_s());
        let tolerance = controls.extent_absolute_tolerance.value().min(1e-8 * donor);
        if tolerance == 0.0 { return Err(AcidBaseError::Arithmetic("positive mineral extent tolerance loses range")); }
        let log_tolerance = det::ln(tolerance);
        let initial_force = force(0.0, direction, log_time_rate, &start)?;
        if start.log_omega == Some(0.0)
            || start.log_omega.is_some_and(|log| log.abs() <= start.spread)
        {
            // This is an explicitly tolerance-qualified unresolved increment,
            // never a claim that the mineral reached thermodynamic equilibrium.
            if initial_force.accepts(log_tolerance) {
                return Ok(Stage { phase: initial, blocked: false, exhausted: false });
            }
            return Err(AcidBaseError::Arithmetic("aqueous pH resolution does not qualify mineral driving"));
        }
        let endpoint = transferred(initial, donor, direction)?;
        let end_trial = self.evaluate(endpoint, cancelled)?;
        let end_force = force(donor, direction, log_time_rate, &end_trial)?;
        if !growing && end_force.upper() <= 0.0 {
            return Ok(Stage { phase: endpoint, blocked: false, exhausted: true });
        }
        if end_force.lower() <= 0.0 {
            return Err(AcidBaseError::Arithmetic("aqueous pH resolution does not resolve mineral exhaustion boundary"));
        }
        let minimum = f64::from_bits(1);
        if minimum >= donor { return Err(AcidBaseError::Arithmetic("nonzero interior mineral extent loses range")); }
        let mut low = 0.0;
        let mut high = donor;
        let predicted_log = log_time_rate + log_abs_drive(start.log_omega);
        let predicted = det::exp(predicted_log.min(det::ln(donor)));
        if predicted == 0.0 { return Err(AcidBaseError::Arithmetic("nonzero mineral transfer underflows")); }
        let mut candidate = if predicted > 0.0 && predicted < high { predicted } else { 0.5 * high };
        for _ in 0..controls.maximum_stage_iterations {
            poll(cancelled)?;
            if candidate <= low || candidate >= high {
                return Err(AcidBaseError::Arithmetic("implicit mineral extent reached floating-point resolution"));
            }
            let phase = transferred(initial, candidate, direction)?;
            let trial = self.evaluate(phase, cancelled)?;
            let residual = force(candidate, direction, log_time_rate, &trial)?;
            if residual.accepts(log_tolerance) {
                return Ok(Stage { phase, blocked: false, exhausted: false });
            }
            if residual.lower() > 0.0 {
                high = candidate;
            } else if residual.upper() < 0.0 {
                low = candidate;
            } else {
                return Err(AcidBaseError::Arithmetic("aqueous pH uncertainty exceeds implicit mineral extent tolerance"));
            }
            candidate = if low == 0.0 {
                // Logarithmic search retains very small finite transfers under
                // a large donor without an unbounded linear bisection count.
                let midpoint = 0.5 * det::ln(minimum) + 0.5 * det::ln(high);
                det::exp(midpoint)
            } else if high - low <= low {
                low + 0.5 * (high - low)
            } else {
                det::exp(0.5 * det::ln(low) + 0.5 * det::ln(high))
            };
        }
        Err(AcidBaseError::Budget)
    }
}

fn transferred(initial: Phase, magnitude: f64, direction: f64) -> Result<Phase> {
    let phase = if direction > 0.0 {
        Phase { dissolved: finite(initial.dissolved + magnitude, "dissolved mineral receiver overflow")?, solid: initial.solid - magnitude }
    } else {
        Phase { dissolved: initial.dissolved - magnitude, solid: finite(initial.solid + magnitude, "mineral solid receiver overflow")? }
    };
    nonnegative(phase.dissolved, "mineral extent exceeds dissolved donor")?;
    nonnegative(phase.solid, "mineral extent exceeds solid donor")?;
    Ok(phase)
}

// Signed physical extent residual divided by a positive logarithmic scale.
// This admits very large Omega without ever constructing exp(logOmega) or an
// overflowing h*K product. The numerical aqueous spread is scaled identically.
struct Force { value: f64, uncertainty: f64, log_scale: f64 }
impl Force {
    fn lower(&self) -> f64 { self.value - self.uncertainty }
    fn upper(&self) -> f64 { self.value + self.uncertainty }
    fn accepts(&self, log_tolerance: f64) -> bool {
        let bound = self.value.abs() + self.uncertainty;
        bound == 0.0 || det::ln(bound) + self.log_scale <= log_tolerance
    }
}
fn force(magnitude: f64, direction: f64, log_time_rate: f64, trial: &Trial) -> Result<Force> {
    let log_amount = if magnitude == 0.0 { f64::NEG_INFINITY } else { det::ln(magnitude) };
    let log_rate = log_time_rate + log_abs_drive(trial.log_omega);
    let log_spread = if trial.spread == 0.0 || trial.log_omega.is_none() { f64::NEG_INFINITY } else {
        log_time_rate + trial.log_omega.expect("positive IAP") + log_positive_expm1(trial.spread)
    };
    let scale = log_amount.max(log_rate).max(log_spread);
    if scale == f64::NEG_INFINITY { return Ok(Force { value: 0.0, uncertainty: 0.0, log_scale: 0.0 }); }
    finite(scale, "logarithmic mineral extent scale")?;
    let sign = if trial.log_omega.is_some_and(|log| log > 0.0) { 1.0 } else { -1.0 };
    let amount = if magnitude == 0.0 { 0.0 } else { det::exp(log_amount - scale) };
    let rate = if log_rate == f64::NEG_INFINITY { 0.0 } else { det::exp(log_rate - scale) };
    let uncertainty = if log_spread == f64::NEG_INFINITY { 0.0 } else { det::exp(log_spread - scale) };
    Ok(Force { value: amount + direction * sign * rate, uncertainty, log_scale: scale })
}
fn log_abs_drive(log_omega: Option<f64>) -> f64 {
    match log_omega {
        None => 0.0,
        Some(0.0) => f64::NEG_INFINITY,
        Some(log) if log > 0.0 => log_positive_expm1(log),
        Some(log) => det::ln(-det::expm1(log)),
    }
}
fn log_positive_expm1(value: f64) -> f64 {
    if value > 50.0 { value } else { det::ln(det::expm1(value)) }
}

fn admit(duration: Time, c: MineralKineticControlsV1) -> Result<()> {
    let relative = c.relative_tolerance.value();
    let absolute = c.amount_absolute_tolerance.value();
    if !duration.value().is_finite() || duration.value() < 0.0
        || !relative.is_finite() || !(0.0..=1.0).contains(&relative)
        || !absolute.is_finite() || absolute < 0.0 || (absolute == 0.0 && relative == 0.0)
        || !c.extent_absolute_tolerance.value().is_finite() || c.extent_absolute_tolerance.value() <= 0.0
        || [c.initial_step, c.minimum_step, c.maximum_step].iter().any(|t| !t.value().is_finite() || t.value() <= 0.0)
        || c.minimum_step.value() > c.initial_step.value() || c.initial_step.value() > c.maximum_step.value()
        || !(1..=256).contains(&c.maximum_stage_iterations) || !(1..=65536).contains(&c.maximum_attempts)
    { Err(AcidBaseError::Invalid("mineral kinetic duration, tolerances or bounded work controls")) } else { Ok(()) }
}
fn temporal_error(full: Phase, fine: Phase, controls: MineralKineticControlsV1) -> f64 {
    let mut error: f64 = 0.0;
    for (a, b) in [(full.dissolved, fine.dissolved), (full.solid, fine.solid)] {
        if a == b { continue; }
        let magnitude = a.max(b);
        let absolute = controls.amount_absolute_tolerance.value();
        let scale = magnitude.max(absolute);
        error = error.max(((a - b).abs() / scale)
            / (absolute / scale + controls.relative_tolerance.value() * (magnitude / scale)));
    }
    error
}
