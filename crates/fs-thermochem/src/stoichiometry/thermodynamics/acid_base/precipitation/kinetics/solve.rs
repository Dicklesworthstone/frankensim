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
            // Repeated nondecimal intervals can leave a sub-minimum final
            // time sliver. Qualify ALL remaining motion before omitting it;
            // this is never a per-stage licence to discard small transfers.
            if h == remaining && h < controls.minimum_step.value() {
                let trial = workspace.evaluate(current, &mut cancelled)?;
                let bound = force(0.0, 1.0, det::ln(h) + det::ln(self.rate().rate_mol_s()), &trial)?;
                if let Some(error) = omitted_tail_error_norm(&bound, current, controls) {
                    maximum_error = maximum_error.max(error);
                    break;
                }
            }
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
        // Affinity decreases along the physical direction, so the initial
        // rate (including aqueous uncertainty) bounds the implicit transfer.
        // A short supported step need not speciate the distant depleted bath.
        let mut high = forward_extent_upper(&initial_force, donor)?;
        let endpoint = transferred(initial, high, direction)?;
        let end_trial = self.evaluate(endpoint, cancelled)?;
        let end_force = force(high, direction, log_time_rate, &end_trial)?;
        if high == donor && !growing && end_force.upper() <= 0.0 {
            return Ok(Stage { phase: endpoint, blocked: false, exhausted: true });
        }
        if high < donor && end_force.accepts(log_tolerance) {
            let phase = accept_transfer(initial, endpoint, high)?;
            return Ok(Stage { phase, blocked: false, exhausted: false });
        }
        if end_force.lower() <= 0.0 {
            return Err(AcidBaseError::Arithmetic("aqueous pH resolution does not resolve mineral extent bound"));
        }
        let minimum = f64::from_bits(1);
        if minimum >= high { return Err(AcidBaseError::Arithmetic("nonzero interior mineral extent loses range")); }
        let mut low = 0.0;
        let predicted_log = log_time_rate + log_abs_drive(start.log_omega);
        let predicted = det::exp(predicted_log.min(det::ln(high)));
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
                let phase = accept_transfer(initial, phase, candidate)?;
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

// Tiny trial extents are valid during logarithmic bracketing. An accepted
// nonzero extent must change at least one independently retained phase stock;
// otherwise repeated windows would silently erase a finite accumulated transfer.
fn accept_transfer(initial: Phase, phase: Phase, magnitude: f64) -> Result<Phase> {
    if magnitude > 0.0
        && phase.dissolved == initial.dissolved
        && phase.solid == initial.solid
    {
        Err(AcidBaseError::Arithmetic(
            "nonzero mineral transfer is below phase-stock resolution",
        ))
    } else {
        Ok(phase)
    }
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

fn forward_extent_upper(initial_force: &Force, donor: f64) -> Result<f64> {
    let bound = initial_force.value.abs() + initial_force.uncertainty;
    let log_bound = det::ln(bound) + initial_force.log_scale;
    let log_donor = det::ln(donor);
    if log_bound >= log_donor { return Ok(donor); }
    if det::exp(log_bound) == 0.0 {
        return Err(AcidBaseError::Arithmetic("nonzero mineral transfer underflows"));
    }
    // Include the logarithm/exponential rounding scale, not merely a single
    // ULP of the reconstructed amount; the inventory cap remains exact.
    let padded_log = log_bound + 32.0 * f64::EPSILON * (1.0 + log_bound.abs());
    if padded_log >= log_donor { return Ok(donor); }
    Ok(det::exp(padded_log).next_up().min(donor))
}

// In this autonomous monotone finite bath, the initial affinity bounds the
// magnitude of the subsequent rate. Force already includes the aqueous
// uncertainty, so its zero-extent bound covers the entire omitted tail.
fn omitted_tail_error_norm(bound: &Force, phase: Phase, controls: MineralKineticControlsV1) -> Option<f64> {
    let magnitude = bound.value.abs() + bound.uncertainty;
    if magnitude == 0.0 { return Some(0.0); }
    let log_bound = det::ln(magnitude) + bound.log_scale;
    if log_bound > det::ln(controls.extent_absolute_tolerance.value()) { return None; }
    let mut log_error = f64::NEG_INFINITY;
    for stock in [phase.dissolved, phase.solid] {
        let absolute = controls.amount_absolute_tolerance.value();
        let scale = stock.max(absolute);
        if scale == 0.0 { return None; }
        let tolerance = absolute / scale + controls.relative_tolerance.value() * (stock / scale);
        if tolerance == 0.0 { return None; }
        log_error = log_error.max(log_bound - det::ln(scale) - det::ln(tolerance));
    }
    if log_error > 0.0 { None } else { Some(det::exp(log_error)) }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tiny_bracket_trials_are_allowed_but_unresolved_acceptance_is_refused() {
        let initial = Phase { dissolved: 0.001, solid: 0.109 };
        for direction in [1.0, -1.0] {
            let phase = transferred(initial, 1.0e-21, direction).unwrap();
            assert_eq!(phase.dissolved, initial.dissolved);
            assert_eq!(phase.solid, initial.solid);
            assert!(matches!(
                accept_transfer(initial, phase, 1.0e-21),
                Err(AcidBaseError::Arithmetic(_))
            ));
        }
    }

    #[test]
    fn trace_donor_exhaustion_remains_resolved_with_an_unchanged_receiver() {
        let dissolving = Phase { dissolved: 1.0, solid: 1.0e-30 };
        let phase = transferred(dissolving, 1.0e-30, 1.0).unwrap();
        assert_eq!(phase.dissolved, dissolving.dissolved);
        assert_eq!(phase.solid, 0.0);
        assert!(accept_transfer(dissolving, phase, 1.0e-30).is_ok());

        let growing = Phase { dissolved: 1.0e-30, solid: 1.0 };
        let phase = transferred(growing, 1.0e-30, -1.0).unwrap();
        assert_eq!(phase.dissolved, 0.0);
        assert_eq!(phase.solid, growing.solid);
        assert!(accept_transfer(growing, phase, 1.0e-30).is_ok());
    }

    #[test]
    fn zero_extent_can_retain_both_phase_stocks() {
        let initial = Phase { dissolved: 0.001, solid: 0.109 };
        assert!(accept_transfer(initial, initial, 0.0).is_ok());
    }

    #[test]
    fn nondecimal_time_tail_has_a_reported_whole_tail_error_bound() {
        let elapsed = (0..80).fold(0.0, |elapsed, _| elapsed + 0.00625);
        let remaining: f64 = 0.5 - elapsed;
        let controls = MineralKineticControlsV1::default();
        assert!(remaining > 0.0 && remaining < controls.minimum_step.value());
        let bound = Force { value: 1.0, uncertainty: 1.0e-10, log_scale: det::ln(remaining) + det::ln(1.0e-4) };
        let phase = Phase { dissolved: 0.01, solid: 0.1 };
        let error = omitted_tail_error_norm(&bound, phase, controls).unwrap();
        assert!(error > 0.0 && error <= 1.0);
    }

    #[test]
    fn final_tail_is_not_omitted_when_its_motion_exceeds_stock_accuracy() {
        let mut controls = MineralKineticControlsV1::default();
        controls.relative_tolerance = Dimensionless::new(1.0e-16);
        controls.amount_absolute_tolerance = Amount::new(0.0);
        controls.extent_absolute_tolerance = Amount::new(1.0e-12);
        let bound = Force { value: 1.0, uncertainty: 0.0, log_scale: det::ln(1.0e-18) };
        let phase = Phase { dissolved: 0.001, solid: 0.109 };
        assert!(omitted_tail_error_norm(&bound, phase, controls).is_none());
    }

    #[test]
    fn forward_extent_bound_is_local_uncertainty_aware_and_inventory_capped() {
        let bound = Force { value: -1.0, uncertainty: 0.125, log_scale: det::ln(0.25) };
        let upper = forward_extent_upper(&bound, 100.0).unwrap();
        assert!(upper >= 0.28125 && upper < 0.2812500000001);
        assert_eq!(forward_extent_upper(&bound, 0.1).unwrap(), 0.1);
        let large = Force { value: 1.0, uncertainty: 1.0, log_scale: 1000.0 };
        assert_eq!(forward_extent_upper(&large, 0.1).unwrap(), 0.1);
    }
}
