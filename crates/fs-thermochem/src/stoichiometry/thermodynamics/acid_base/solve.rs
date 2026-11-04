//! Logarithmic family fractions and one monotone electroneutrality solve.
use super::*;
use fs_math::det;

const LN_10: f64 = core::f64::consts::LN_10;

/// Prepared identities, source owner and bounded inventory buffers for nested
/// phase trials. Only declared amounts may vary after the original admission.
pub(super) struct PreparedAcidBase<'a> {
    model: &'a IdealAcidBaseSystem,
    ordered: Vec<f64>,
    spectators: Vec<SpectatorIon>,
    basis: f64,
    options: SpeciationOptions,
}

/// Compact accepted charge root. The complete species/source report is built
/// only after the enclosing phase solve accepts its component inventories.
#[derive(Clone, Copy)]
pub(super) struct AcidBaseTrial {
    point: Point,
    bracket: [f64; 2],
    iterations: usize,
}
impl AcidBaseTrial {
    pub(super) fn ph(&self) -> f64 {
        self.point.ph
    }
    pub(super) fn ph_bracket(&self) -> [f64; 2] {
        self.bracket
    }
    pub(super) fn iterations(&self) -> usize {
        self.iterations
    }
}

impl IdealAcidBaseSystem {
    /// Solve all declared families simultaneously on one explicit ideal basis.
    /// Supply every analytical-family total exactly once, including zero pools.
    /// Spectator ions are fixed, caller-declared stocks; missing counterions
    /// are never fabricated. The answer is specific to those actual declarations.
    pub fn speciate(
        &self,
        totals: &[(String, Amount)],
        spectators: &[SpectatorIon],
        amount_at_unit_activity: Amount,
        options: SpeciationOptions,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AcidBaseSpeciation> {
        self.speciate_partitioned(
            totals,
            spectators,
            amount_at_unit_activity,
            options,
            &[],
            cancelled,
        )
    }

    // A finite neutral gas adds a zero-charge statistical weight to its
    // analytical family. Charge remains monotone in pH, so the original ONE
    // scalar solve closes gas partition and every aqueous family together.
    pub(super) fn speciate_partitioned(
        &self,
        totals: &[(String, Amount)],
        spectators: &[SpectatorIon],
        amount_at_unit_activity: Amount,
        options: SpeciationOptions,
        gas_weights: &[Option<gas_liquid::GasWeight>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AcidBaseSpeciation> {
        let prepared = self.prepare_speciation(
            totals,
            spectators,
            amount_at_unit_activity,
            options,
            &mut cancelled,
        )?;
        let trial = prepared.solve_partitioned(gas_weights, &mut cancelled)?;
        prepared.finish_partitioned(trial, gas_weights, &mut cancelled)
    }

    /// Admit and order the original species once for repeated phase trials.
    /// This is internal to the existing chemistry owners, not a new activity law.
    pub(super) fn prepare_speciation(
        &self,
        totals: &[(String, Amount)],
        spectators: &[SpectatorIon],
        amount_at_unit_activity: Amount,
        options: SpeciationOptions,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedAcidBase<'_>> {
        poll(cancelled)?;
        let basis = amount_at_unit_activity.value();
        if !basis.is_finite() || basis <= 0.0 {
            return Err(AcidBaseError::Invalid(
                "common amount at unit activity must be positive",
            ));
        }
        let bracket = options.ph_bracket.map(|x| x.value());
        if bracket.iter().any(|x| !x.is_finite())
            || bracket[0] >= bracket[1]
            || !(bracket[1] - bracket[0]).is_finite()
            || !options.ph_tolerance.value().is_finite()
            || options.ph_tolerance.value() <= 0.0
            || options.ph_tolerance.value() > 1.0
            || !options.absolute_charge_tolerance.value().is_finite()
            || options.absolute_charge_tolerance.value() <= 0.0
            || !options.relative_charge_tolerance.value().is_finite()
            || !(0.0..=0.01).contains(&options.relative_charge_tolerance.value())
            || !(1..=1024).contains(&options.max_iterations)
        {
            return Err(AcidBaseError::Invalid(
                "invalid pH bracket charge tolerances or work budget",
            ));
        }
        if totals.len() != self.inner.families.len() {
            return Err(AcidBaseError::Invalid(
                "supply every analytical-family total exactly once",
            ));
        }
        let mut indexed = BTreeMap::new();
        for (id, amount) in totals {
            nonnegative(
                amount.value(),
                "analytical family totals must be finite and nonnegative",
            )?;
            if indexed.insert(id.as_str(), amount.value()).is_some()
                || !self.inner.families.iter().any(|family| family.id == *id)
            {
                return Err(AcidBaseError::Invalid(
                    "duplicate or unknown analytical-family total",
                ));
            }
        }
        let ordered: Vec<_> = self
            .inner
            .families
            .iter()
            .map(|f| indexed[f.id.as_str()])
            .collect();
        if spectators.len() > 256 {
            return Err(AcidBaseError::Invalid("too many spectator ions"));
        }
        let mut spectators = spectators.to_vec();
        spectators.sort_by(|a, b| a.species.cmp(&b.species));
        if spectators.windows(2).any(|p| p[0].species == p[1].species) {
            return Err(AcidBaseError::Invalid("duplicate spectator identity"));
        }
        for ion in &spectators {
            nonnegative(
                ion.amount.value(),
                "spectator amounts must be finite and nonnegative",
            )?;
            if ion.charge == 0
                || !(-32..=32).contains(&ion.charge)
                || self
                    .inner
                    .network
                    .species()
                    .binary_search(&ion.species)
                    .is_ok()
            {
                return Err(AcidBaseError::Invalid(
                    "spectators need nonzero bounded charge and cannot duplicate reactive species",
                ));
            }
        }
        Ok(PreparedAcidBase {
            model: self,
            ordered,
            spectators,
            basis,
            options,
        })
    }

    fn report(
        &self,
        totals: &[f64],
        spectators: Vec<SpectatorIon>,
        basis: f64,
        answer: Point,
        bracket: [f64; 2],
        iterations: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcidBaseSpeciation> {
        let mut families = Vec::with_capacity(totals.len());
        let ln_h = det::ln(answer.h);
        let mut maximum = (ln_h + det::ln(answer.oh)
            - self.inner.water.standard.log_equilibrium_constant().value())
        .abs();
        for (family, &total) in self.inner.families.iter().zip(totals) {
            poll(cancelled)?;
            let fractions = fractions(family, ln_h)?;
            let mut forms = Vec::with_capacity(family.species.len());
            let mut sum = Sum::default();
            for (j, id) in family.species.iter().enumerate() {
                let amount = member_amount(total, fractions[j])?;
                let activity = if amount == 0.0 {
                    0.0
                } else {
                    positive(amount / basis, "family activity loses range")?
                };
                sum.add(amount)?;
                forms.push(AcidBaseSpecies {
                    species: id.clone(),
                    amount: Amount::new(amount),
                    fraction: Dimensionless::new(fractions[j]),
                    activity: Dimensionless::new(activity),
                });
            }
            if (sum.finish()? - total).abs() > 64.0 * f64::EPSILON * total {
                return Err(AcidBaseError::Arithmetic(
                    "speciated family no longer closes its analytical total",
                ));
            }
            if total > 0.0 {
                for (j, step) in family.steps.iter().enumerate() {
                    let residual = ln_h + det::ln(forms[j + 1].activity.value())
                        - det::ln(forms[j].activity.value())
                        - step.standard.log_equilibrium_constant().value();
                    if !residual.is_finite() {
                        return Err(AcidBaseError::Arithmetic(
                            "mass-action residual loses range",
                        ));
                    }
                    maximum = maximum.max(residual.abs());
                }
            }
            families.push(FamilySpeciation {
                id: family.id.clone(),
                total: Amount::new(total),
                forms,
            });
        }
        poll(cancelled)?;
        Ok(AcidBaseSpeciation {
            model: self.clone(),
            amount_at_unit_activity: Amount::new(basis),
            ph: Dimensionless::new(answer.ph),
            proton_activity: Dimensionless::new(answer.h),
            hydroxide_activity: Dimensionless::new(answer.oh),
            proton_amount: Amount::new(answer.n_h),
            hydroxide_amount: Amount::new(answer.n_oh),
            families,
            spectators,
            charge_balance_residual: Amount::new(answer.charge),
            max_log_mass_action_residual: Dimensionless::new(maximum),
            ph_bracket: bracket.map(Dimensionless::new),
            iterations,
        })
    }
}

impl PreparedAcidBase<'_> {
    /// Canonical slot lookup is paid once when preparing the enclosing owner.
    pub(super) fn spectator_index(&self, species: &SpeciesId) -> Option<usize> {
        self.spectators
            .binary_search_by(|ion| ion.species.cmp(species))
            .ok()
    }

    pub(super) fn set_family_total(&mut self, index: usize, amount: Amount) -> Result<()> {
        nonnegative(amount.value(), "invalid prepared family amount")?;
        let target = self
            .ordered
            .get_mut(index)
            .ok_or(AcidBaseError::Invalid("prepared family index"))?;
        *target = amount.value();
        Ok(())
    }

    pub(super) fn set_spectator_amount(&mut self, index: usize, amount: Amount) -> Result<()> {
        nonnegative(amount.value(), "invalid prepared spectator amount")?;
        let target = self
            .spectators
            .get_mut(index)
            .ok_or(AcidBaseError::Invalid("prepared spectator index"))?;
        target.amount = amount;
        Ok(())
    }

    /// Same charge equation and acceptance bounds, with no trial allocation.
    pub(super) fn solve(&self, cancelled: &mut impl FnMut() -> bool) -> Result<AcidBaseTrial> {
        self.solve_partitioned(&[], cancelled)
    }

    fn solve_partitioned(
        &self,
        gas_weights: &[Option<gas_liquid::GasWeight>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcidBaseTrial> {
        if !gas_weights.is_empty() && gas_weights.len() != self.model.inner.families.len() {
            return Err(AcidBaseError::Invalid(
                "gas weights must match the family axis",
            ));
        }
        let options = self.options;
        let mut bracket = options.ph_bracket.map(|x| x.value());
        let evaluate = |ph, cancelled: &mut _| {
            point(
                &self.model.inner,
                &self.ordered,
                &self.spectators,
                self.basis,
                ph,
                gas_weights,
                cancelled,
            )
        };
        let lower = evaluate(bracket[0], cancelled)?;
        let upper = evaluate(bracket[1], cancelled)?;
        // Charge decreases strictly with pH: free H+ decreases, OH- increases,
        // and every family's mean charge decreases as it loses protons.
        if lower.charge < 0.0 || upper.charge > 0.0 {
            return Err(AcidBaseError::OutsideBracket);
        }
        let (answer, iterations) = if lower.charge == 0.0 {
            bracket = [bracket[0]; 2];
            (lower, 0)
        } else if upper.charge == 0.0 {
            bracket = [bracket[1]; 2];
            (upper, 0)
        } else {
            let mut accepted = None;
            for iteration in 1..=options.max_iterations {
                let ph = bracket[0] + 0.5 * (bracket[1] - bracket[0]);
                if ph == bracket[0] || ph == bracket[1] {
                    return Err(AcidBaseError::Arithmetic(
                        "pH inversion reached floating resolution",
                    ));
                }
                let value = evaluate(ph, cancelled)?;
                let tolerance = options.absolute_charge_tolerance.value()
                    + options.relative_charge_tolerance.value() * value.unsigned_charge;
                positive(tolerance, "charge tolerance loses range")?;
                if bracket[1] - bracket[0] <= options.ph_tolerance.value()
                    && value.charge.abs() <= tolerance
                {
                    accepted = Some((value, iteration));
                    break;
                }
                if value.charge > 0.0 {
                    bracket[0] = ph;
                } else {
                    bracket[1] = ph;
                }
            }
            accepted.ok_or(AcidBaseError::Budget)?
        };
        Ok(AcidBaseTrial {
            point: answer,
            bracket,
            iterations,
        })
    }

    /// Identical terminal-member arithmetic to the complete aqueous report.
    pub(super) fn terminal_activity(&self, index: usize, trial: &AcidBaseTrial) -> Result<f64> {
        let family = self
            .model
            .inner
            .families
            .get(index)
            .ok_or(AcidBaseError::Invalid("prepared terminal family index"))?;
        let fractions = fractions(family, det::ln(trial.point.h))?;
        let amount = member_amount(self.ordered[index], fractions[family.species.len() - 1])?;
        if amount == 0.0 {
            Ok(0.0)
        } else {
            positive(amount / self.basis, "family activity loses range")
        }
    }

    /// Restore the inventories of the accepted trial before consuming these
    /// buffers. A stale trial cannot publish a report for different charge data.
    pub(super) fn finish(
        self,
        trial: AcidBaseTrial,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcidBaseSpeciation> {
        self.finish_partitioned(trial, &[], cancelled)
    }

    fn finish_partitioned(
        mut self,
        trial: AcidBaseTrial,
        gas_weights: &[Option<gas_liquid::GasWeight>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcidBaseSpeciation> {
        let restored = point(
            &self.model.inner,
            &self.ordered,
            &self.spectators,
            self.basis,
            trial.ph(),
            gas_weights,
            cancelled,
        )?;
        if restored != trial.point {
            return Err(AcidBaseError::Invalid(
                "accepted aqueous trial no longer matches the prepared inventories",
            ));
        }
        // The original aqueous report owns only dissolved stocks, while the
        // enclosing gas owner retains the finite cross-phase totals.
        for (index, total) in self.ordered.iter_mut().enumerate() {
            if let Some(Some(weight)) = gas_weights.get(index) {
                let fractions = fractions(&self.model.inner.families[index], -trial.ph() * LN_10)?;
                *total = weight.partition(*total, fractions[weight.neutral])?[0];
            }
        }
        self.model.report(
            &self.ordered,
            self.spectators,
            self.basis,
            trial.point,
            trial.bracket,
            trial.iterations,
            cancelled,
        )
    }
}

#[derive(Clone, Copy, PartialEq)]
struct Point {
    ph: f64,
    h: f64,
    oh: f64,
    n_h: f64,
    n_oh: f64,
    charge: f64,
    unsigned_charge: f64,
}
#[allow(clippy::too_many_arguments)]
fn point(
    layout: &Layout,
    totals: &[f64],
    spectators: &[SpectatorIon],
    basis: f64,
    ph: f64,
    gas_weights: &[Option<gas_liquid::GasWeight>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Point> {
    poll(cancelled)?;
    let ln_h = -ph * LN_10;
    let h = positive(det::exp(ln_h), "proton activity loses range")?;
    let oh = positive(
        det::exp(layout.water.standard.log_equilibrium_constant().value() - ln_h),
        "hydroxide activity loses range",
    )?;
    let n_h = positive(basis * h, "proton amount loses range")?;
    let n_oh = positive(basis * oh, "hydroxide amount loses range")?;
    let mut net = Sum::default();
    let mut absolute = Sum::default();
    add_charge(&mut net, &mut absolute, n_h, 1)?;
    add_charge(&mut net, &mut absolute, n_oh, -1)?;
    for (index, (family, &total)) in layout.families.iter().zip(totals).enumerate() {
        poll(cancelled)?;
        let fractions = fractions(family, ln_h)?;
        let total = match gas_weights.get(index) {
            Some(Some(weight)) => weight.partition(total, fractions[weight.neutral])?[0],
            _ => total,
        };
        for (j, &fraction) in fractions.iter().enumerate().take(family.species.len()) {
            let amount = member_amount(total, fraction)?;
            add_charge(
                &mut net,
                &mut absolute,
                amount,
                family.steps[0].acid_charge - j as i32,
            )?;
        }
    }
    for ion in spectators {
        add_charge(&mut net, &mut absolute, ion.amount.value(), ion.charge)?;
    }
    Ok(Point {
        ph,
        h,
        oh,
        n_h,
        n_oh,
        charge: net.finish()?,
        unsigned_charge: absolute.finish()?,
    })
}

pub(super) fn fractions(family: &AcidFamily, ln_h: f64) -> Result<[f64; MAX_DISSOCIATIONS + 1]> {
    let mut logs = [0.0; MAX_DISSOCIATIONS + 1];
    let mut maximum: f64 = 0.0;
    for (j, step) in family.steps.iter().enumerate() {
        logs[j + 1] = logs[j] + step.standard.log_equilibrium_constant().value() - ln_h;
        if !logs[j + 1].is_finite() {
            return Err(AcidBaseError::Arithmetic(
                "acid-family log weights lose range",
            ));
        }
        maximum = maximum.max(logs[j + 1]);
    }
    let mut normalization = Sum::default();
    for &value in logs.iter().take(family.species.len()) {
        // No numerical seed or zeroed physical member: every relative weight
        // must remain representable under this bounded ideal closure.
        normalization.add(positive(
            det::exp(value - maximum),
            "positive acid-family weight underflows",
        )?)?;
    }
    let log_normalization = det::ln(positive(
        normalization.finish()?,
        "family normalization loses range",
    )?);
    let mut result = [0.0; MAX_DISSOCIATIONS + 1];
    for j in 0..family.species.len() {
        result[j] = positive(
            det::exp((logs[j] - maximum) - log_normalization),
            "positive acid-family fraction underflows",
        )?;
        if result[j] > 1.0 {
            return Err(AcidBaseError::Arithmetic("family fraction exceeds one"));
        }
    }
    Ok(result)
}
fn member_amount(total: f64, fraction: f64) -> Result<f64> {
    if total == 0.0 {
        Ok(0.0)
    } else {
        positive(total * fraction, "positive family member amount underflows")
    }
}
fn add_charge(net: &mut Sum, absolute: &mut Sum, amount: f64, charge: i32) -> Result<()> {
    if charge == 0 || amount == 0.0 {
        return Ok(());
    }
    let z = f64::from(charge);
    let value = z * amount;
    if !value.is_finite() || value == 0.0 {
        return Err(AcidBaseError::Arithmetic(
            "ionic charge account loses range",
        ));
    }
    net.add(value)?;
    net.add(z.mul_add(amount, -value))?;
    absolute.add(value.abs())?;
    Ok(())
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(AcidBaseError::Cancelled)
    } else {
        Ok(())
    }
}
