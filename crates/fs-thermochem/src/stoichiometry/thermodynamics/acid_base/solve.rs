//! Logarithmic family fractions and one monotone electroneutrality solve.
use super::*;
use fs_math::det;

const LN_10: f64 = core::f64::consts::LN_10;

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
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AcidBaseSpeciation> {
        poll(&mut cancelled)?;
        let basis = amount_at_unit_activity.value();
        if !basis.is_finite() || basis <= 0.0 {
            return Err(AcidBaseError::Invalid(
                "common amount at unit activity must be positive",
            ));
        }
        let mut bracket = options.ph_bracket.map(|x| x.value());
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
        let evaluate =
            |ph, cancelled: &mut _| point(&self.inner, &ordered, &spectators, basis, ph, cancelled);
        let lower = evaluate(bracket[0], &mut cancelled)?;
        let upper = evaluate(bracket[1], &mut cancelled)?;
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
                let value = evaluate(ph, &mut cancelled)?;
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
        self.report(
            &ordered,
            spectators,
            basis,
            answer,
            bracket,
            iterations,
            &mut cancelled,
        )
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

struct Point {
    ph: f64,
    h: f64,
    oh: f64,
    n_h: f64,
    n_oh: f64,
    charge: f64,
    unsigned_charge: f64,
}
fn point(
    layout: &Layout,
    totals: &[f64],
    spectators: &[SpectatorIon],
    basis: f64,
    ph: f64,
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
    for (family, &total) in layout.families.iter().zip(totals) {
        poll(cancelled)?;
        let fractions = fractions(family, ln_h)?;
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

fn fractions(family: &AcidFamily, ln_h: f64) -> Result<[f64; MAX_DISSOCIATIONS + 1]> {
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
