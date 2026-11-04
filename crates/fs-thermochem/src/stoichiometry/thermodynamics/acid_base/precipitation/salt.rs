//! One finite stoichiometric salt coupled to an acid/base family and charge balance.
//!
//! The admitted reaction is `p M(z_M+) + q A(z_A-) + w H2O -> solid`, where
//! A is the terminal member of one original acid family, p and q are positive
//! integers, p*z_M = q*abs(z_A), and w is a nonnegative integer. Pure solid and
//! fixed solvent have unit activity. Both analytical component totals include
//! the corresponding p or q units per mole of any initial solid. Other families
//! and spectator inventories are independently fixed.
//!
//! With u the limiting dissolved formula-equivalent amount, dissolved metal is
//! p*u plus its stoichiometric excess and dissolved family is q*u plus its
//! excess. Every trial recomputes the shared pH and every acid/base family.
//! This retains positive traces without subtracting a nearly complete reaction
//! extent from a large stock. Saturation uses the original activity product
//! a_M^p*a_A^q, not an apparent one-particle salt concentration.
//!
//! The terminal-anion saturation residual is strictly increasing in ln(u),
//! with derivative at least the limiting component's integer coefficient.
//! The family's mean charge is at least z_A, so adding p cations and q family
//! units adds nonnegative charge at fixed pH. Charge balance raises pH and the
//! terminal fraction. This supplies the original bounded scalar solve.
//!
//! This is one explicitly selected solid under the existing ideal-activity,
//! fixed-solvent model. It does not select among competing solids, include
//! omitted metal complexes or redox reactions, or assert a closed H/O solvent
//! inventory. Source uncertainty is retained as evidence, not propagated into
//! the numerical pH or saturation brackets.

use super::super::{
    AcidBaseError, AcidBaseSpeciation, ConservedReaction, FormationGibbsSource,
    IdealAcidBaseSystem, ReactionNetwork, Result, SpeciesId, SpectatorIon, StandardReactionGibbs,
    Sum, charge, consistent_sources, nonnegative, positive, require_phase, solve::AcidBaseTrial,
};
use super::{PrecipitationOptions, PrecipitationStatus};
use fs_math::det;
use fs_qty::{Amount, Dimensionless};

/// Original source-bound integer salt reaction and the complete aqueous owner.
#[derive(Clone, Debug)]
pub struct SaltPrecipitation {
    aqueous: IdealAcidBaseSystem,
    reaction: ConservedReaction,
    standard: StandardReactionGibbs,
    cation: SpeciesId,
    solid: SpeciesId,
    anion: SpeciesId,
    charge: i32,
    cation_stoichiometry: u32,
    anion_stoichiometry: u32,
    family_index: usize,
    network: ReactionNetwork,
}

/// One coupled charge, acid/base and finite salt partition.
#[derive(Clone, Debug)]
pub struct SaltSpeciation {
    model: SaltPrecipitation,
    /// Analytical cation-component amount, including any initial solid [mol].
    pub cation_total: Amount,
    /// Analytical family-component amount, including any initial solid [mol].
    pub family_total: Amount,
    /// Dissolved free cation [mol], retained independently of solid subtraction.
    pub dissolved_cation: Amount,
    /// Total dissolved amount of all forms in the selected acid family [mol].
    pub dissolved_family: Amount,
    /// Selected pure solid [mol of the original one-formula-unit reaction].
    pub solid_amount: Amount,
    /// Independently reduced dissolved cation + p*solid - cation total [mol].
    pub cation_balance_residual: Amount,
    /// Independently reduced dissolved family + q*solid - family total [mol].
    pub family_balance_residual: Amount,
    /// Complete aqueous result at the simultaneously solved partition and pH.
    pub aqueous: AcidBaseSpeciation,
    /// Natural log of IAP/K for dissolution; positive means supersaturation.
    /// Undefined if either analytical component is empty. This is not log10 SI.
    pub log_saturation_ratio: Option<Dimensionless>,
    /// Bound from the retained inner pH bracket on natural-log saturation
    /// variation. A numerical resolution measure, not an interval certificate.
    pub log_saturation_numerical_spread: Option<Dimensionless>,
    /// Numerical bracket on limiting dissolved formula-equivalent amount [mol].
    /// Divide dissolved cation by p or dissolved family by q, whichever limits.
    pub dissolved_limiting_bracket: [Amount; 2],
    /// Present, absent, or at least one empty analytical component.
    pub status: PrecipitationStatus,
    /// Outer salt bisections, excluding endpoint evaluations.
    pub iterations: usize,
    /// Number of complete bounded aqueous solves used.
    pub aqueous_solves: usize,
    /// Total inner pH bisections, including endpoint solves.
    pub aqueous_iterations: usize,
}
impl SaltSpeciation {
    /// Original phase equation, source receipts and complete aqueous model.
    pub fn model(&self) -> &SaltPrecipitation {
        &self.model
    }
}

impl SaltPrecipitation {
    /// Bind exactly p M(z_M+) + q A(z_A-) + w H2O -> neutral solid.
    ///
    /// Cation and terminal family anion have negative integer coefficients
    /// -p and -q with 1 <= p,q <= u32::MAX; solid has coefficient +1.
    /// Optional water has a negative integer coefficient. The original exact
    /// equation must satisfy p*z_M + q*z_A = 0. The family anion is inferred
    /// from its original complete protonation chain, never from a formula label.
    /// All shared species must retain identical original selected sources.
    pub fn from_sources(
        aqueous: IdealAcidBaseSystem,
        reaction: ConservedReaction,
        cation: SpeciesId,
        solid: SpeciesId,
        family_id: &str,
        sources: &[FormationGibbsSource<'_>],
    ) -> Result<Self> {
        let family_index = aqueous
            .families()
            .iter()
            .position(|family| family.id() == family_id)
            .ok_or(AcidBaseError::Invalid("missing salt acid/base family"))?;
        let anion = aqueous.families()[family_index]
            .species()
            .last()
            .ok_or(AcidBaseError::Invalid("empty salt acid/base family"))?
            .clone();
        let water = aqueous.water().water();
        let coefficient = |id: &SpeciesId| {
            reaction
                .terms()
                .iter()
                .find(|term| &term.species == id)
                .map(|term| term.coefficient)
        };
        let z = charge(&reaction, &cation)?;
        let anion_charge = charge(&reaction, &anion)?;
        let cation_stoichiometry = coefficient(&cation)
            .and_then(|value| value.checked_neg())
            .and_then(|value| u32::try_from(value).ok())
            .filter(|value| *value > 0);
        let anion_stoichiometry = coefficient(&anion)
            .and_then(|value| value.checked_neg())
            .and_then(|value| u32::try_from(value).ok())
            .filter(|value| *value > 0);
        let water_coefficient = coefficient(water);
        if !(1..=32).contains(&z)
            || !(-32..=-1).contains(&anion_charge)
            || charge(&reaction, &solid)? != 0
            || cation_stoichiometry.is_none()
            || anion_stoichiometry.is_none()
            || cation_stoichiometry
                .zip(anion_stoichiometry)
                .is_some_and(|(p, q)| {
                    i128::from(p) * i128::from(z) + i128::from(q) * i128::from(anion_charge) != 0
                })
            || coefficient(&solid) != Some(1)
            || !matches!(water_coefficient, None | Some(..=-1))
            || reaction.terms().len() != 3 + usize::from(water_coefficient.is_some())
            || aqueous
                .reaction_network()
                .species()
                .binary_search(&cation)
                .is_ok()
            || aqueous
                .reaction_network()
                .species()
                .binary_search(&solid)
                .is_ok()
        {
            return Err(AcidBaseError::Invalid(
                "requires positive integer cation and terminal-anion coefficients forming one exactly neutral salt, with optional hydrate water",
            ));
        }
        let standard = reaction.standard_gibbs_change(sources)?;
        for source in sources {
            let (phase, pure) = if *source.species == solid {
                ("source-phase-solid", true)
            } else if source.species == water {
                ("source-phase-liquid", true)
            } else {
                ("source-phase-aqueous", false)
            };
            require_phase(source, &standard, phase, pure)?;
        }
        let mut standards = vec![&standard, aqueous.water().standard()];
        for family in aqueous.families() {
            for step in family.steps() {
                standards.push(step.standard());
            }
        }
        consistent_sources(standards)?;
        let mut reactions = aqueous.reaction_network().reactions().to_vec();
        reactions.push(reaction.clone());
        let network = ReactionNetwork::new(reactions)?;
        Ok(Self {
            aqueous,
            reaction,
            standard,
            cation,
            solid,
            anion,
            charge: z,
            cation_stoichiometry: cation_stoichiometry.expect("admitted cation coefficient"),
            anion_stoichiometry: anion_stoichiometry.expect("admitted anion coefficient"),
            family_index,
            network,
        })
    }

    /// Original aqueous owner with all complete analytical families.
    pub fn aqueous_model(&self) -> &IdealAcidBaseSystem {
        &self.aqueous
    }
    /// Original exactly conserved precipitation equation.
    pub fn reaction(&self) -> &ConservedReaction {
        &self.reaction
    }
    /// Original source evidence and standard precipitation equilibrium constant.
    pub fn standard(&self) -> &StandardReactionGibbs {
        &self.standard
    }
    /// All original equations with shared atom/charge definitions checked.
    pub fn reaction_network(&self) -> &ReactionNetwork {
        &self.network
    }
    /// Selected free cation identity and its original equation-defined charge.
    pub fn cation(&self) -> &SpeciesId {
        &self.cation
    }
    /// Explicitly selected pure solid.
    pub fn solid(&self) -> &SpeciesId {
        &self.solid
    }
    /// Terminal member of the selected original acid/base family.
    pub fn anion(&self) -> &SpeciesId {
        &self.anion
    }
    /// Metal-component moles per mole of the original solid formula unit.
    pub const fn cation_stoichiometry(&self) -> u32 {
        self.cation_stoichiometry
    }
    /// Analytical family units per mole of the original solid formula unit.
    pub const fn anion_stoichiometry(&self) -> u32 {
        self.anion_stoichiometry
    }
    /// Analytical family consumed with its original coefficient by the solid.
    pub fn family_id(&self) -> &str {
        self.aqueous.families()[self.family_index].id()
    }

    /// Solve both component balances, pH and every acid family simultaneously.
    ///
    /// The cation total and selected family total include p and q component
    /// units, respectively, per mole of any initial solid formula unit.
    /// Other analytical totals and spectator amounts stay fixed. The caller
    /// must not duplicate owned cation, solid, or any aqueous family member
    /// among spectators. All aqueous activities share the supplied basis.
    ///
    /// The minimum dissolved amount applies to the limiting FORMULA-EQUIVALENT
    /// amount, dissolved cation/p or dissolved family/q. A root below that
    /// caller-declared floor refuses rather than erasing a trace.
    /// Both logarithmic bracket width and saturation residual, including the
    /// inner pH uncertainty, must meet their respective tolerances.
    #[allow(clippy::too_many_arguments)]
    pub fn speciate(
        &self,
        cation_total: Amount,
        family_totals: &[(String, Amount)],
        spectators: &[SpectatorIon],
        amount_at_unit_activity: Amount,
        options: PrecipitationOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SaltSpeciation> {
        poll(&mut cancelled)?;
        nonnegative(cation_total.value(), "invalid salt cation total")?;
        let basis = positive(
            amount_at_unit_activity.value(),
            "positive common aqueous basis required",
        )?;
        let floor = positive(
            options.minimum_dissolved_amount.value(),
            "positive dissolved search floor required",
        )?;
        if [
            options.log_activity_tolerance,
            options.log_saturation_tolerance,
        ]
        .iter()
        .any(|v| !v.value().is_finite() || v.value() <= 0.0)
            || !(1..=256).contains(&options.max_iterations)
            || spectators
                .iter()
                .any(|ion| ion.species == self.cation || ion.species == self.solid)
        {
            return Err(AcidBaseError::Invalid(
                "invalid salt controls or duplicated owned participant",
            ));
        }
        let mut selected = family_totals
            .iter()
            .filter(|(id, _)| id == self.family_id());
        let (_, family_total) = selected
            .next()
            .ok_or(AcidBaseError::Invalid("missing salt analytical total"))?;
        if selected.next().is_some() {
            return Err(AcidBaseError::Invalid("duplicate salt analytical total"));
        }
        let family_total = *family_total;
        nonnegative(family_total.value(), "invalid salt family total")?;
        let (limit, excess, cation_limits) = partition_basis(
            cation_total.value(),
            family_total.value(),
            self.cation_stoichiometry,
            self.anion_stoichiometry,
        )?;
        let p = f64::from(self.cation_stoichiometry);
        let q = f64::from(self.anion_stoichiometry);
        let minimum_slope = if cation_limits { p } else { q };
        let mut ions = spectators.to_vec();
        ions.push(SpectatorIon {
            species: self.cation.clone(),
            charge: self.charge,
            amount: cation_total,
        });
        // Identity/source admission and canonical inventory buffers are shared
        // by every nested pH solve. Only two scalar amounts vary per trial;
        // the complete dissolved-species report is materialized once, below.
        let mut prepared = self.aqueous.prepare_speciation(
            family_totals,
            &ions,
            amount_at_unit_activity,
            options.aqueous,
            &mut cancelled,
        )?;
        let cation_index = prepared
            .spectator_index(&self.cation)
            .ok_or(AcidBaseError::Invalid("missing owned salt cation"))?;
        let mut solves = 0;
        let mut inner_iterations = 0;
        let mut evaluate = |u: f64, cancelled: &mut _| -> Result<Trial> {
            poll(cancelled)?;
            // Preserve the supplied upper endpoint bit-for-bit: reconstructing
            // it as excess + limit can round differently from the caller total.
            let (metal, family) = if u == limit {
                (cation_total.value(), family_total.value())
            } else if cation_limits {
                (p * u, excess + q * u)
            } else {
                (excess + p * u, q * u)
            };
            if !metal.is_finite()
                || !family.is_finite()
                || metal < 0.0
                || family < 0.0
                || metal > cation_total.value()
                || family > family_total.value()
            {
                return Err(AcidBaseError::Arithmetic(
                    "salt component partition loses inventory range",
                ));
            }
            if u > 0.0 && (metal == 0.0 || family == 0.0) {
                return Err(AcidBaseError::Arithmetic(
                    "positive dissolved salt component underflows",
                ));
            }
            prepared.set_family_total(self.family_index, Amount::new(family))?;
            prepared.set_spectator_amount(cation_index, Amount::new(metal))?;
            let aqueous = prepared.solve(cancelled)?;
            solves += 1;
            inner_iterations += aqueous.iterations();
            let residual = if u == 0.0 {
                None
            } else {
                let anion_activity = prepared.terminal_activity(self.family_index, &aqueous)?;
                let mut sum = Sum::default();
                sum.add(p * det::ln(metal))?;
                sum.add(-p * det::ln(basis))?;
                sum.add(q * det::ln(anion_activity))?;
                sum.add(self.standard.log_equilibrium_constant().value())?;
                Some(sum.finish()?)
            };
            // d ln(alpha_terminal)/d pH = ln(10) * (N - mean(step)),
            // The anion enters Q to power q, so its saturation contribution
            // is bounded by q*N*ln(10). Keep this uncertainty in every
            // phase/sign decision and in the final saturation criterion.
            let spread = q
                * self.aqueous.families()[self.family_index].steps().len() as f64
                * core::f64::consts::LN_10
                * aqueous
                    .ph_bracket()
                    .iter()
                    .map(|p| (p - aqueous.ph()).abs())
                    .fold(0.0, f64::max);
            if !spread.is_finite() {
                return Err(AcidBaseError::Arithmetic(
                    "salt aqueous pH resolution loses range",
                ));
            }
            Ok(Trial {
                limiting: u,
                metal,
                family,
                residual,
                spread,
                aqueous,
            })
        };
        let upper = evaluate(limit, &mut cancelled)?;
        let mut bracket = [limit; 2];
        let (answer, status, iterations) = if limit == 0.0 {
            (upper, PrecipitationStatus::EmptyComponent, 0)
        } else if upper.value()? + upper.spread <= 0.0 {
            (upper, PrecipitationStatus::SolidAbsent, 0)
        } else {
            if upper.value()? - upper.spread <= 0.0 {
                return Err(AcidBaseError::Arithmetic(
                    "aqueous pH bracket does not resolve salt phase absence",
                ));
            }
            if floor >= limit {
                return Err(AcidBaseError::Invalid(
                    "dissolved search floor excludes a precipitating salt partition",
                ));
            }
            let lower = evaluate(floor, &mut cancelled)?;
            if lower.value()? - lower.spread > 0.0 {
                return Err(AcidBaseError::OutsideBracket);
            }
            if lower.value()? + lower.spread > 0.0 {
                return Err(AcidBaseError::Arithmetic(
                    "aqueous pH bracket does not resolve salt search floor",
                ));
            }
            let mut logs = [det::ln(floor), det::ln(limit)];
            let mut accepted = if lower.value()? == 0.0 && lower.spread == 0.0 {
                Some((lower, 0))
            } else {
                None
            };
            if accepted.is_some() {
                bracket = [floor; 2];
            }
            for iteration in 1..=options.max_iterations {
                if accepted.is_some() {
                    break;
                }
                poll(&mut cancelled)?;
                let midpoint = logs[0].midpoint(logs[1]);
                if midpoint <= logs[0] || midpoint >= logs[1] {
                    return Err(AcidBaseError::Arithmetic(
                        "salt partition reached logarithmic resolution",
                    ));
                }
                let u = positive(det::exp(midpoint), "limiting salt component loses range")?;
                if u <= floor || u >= limit {
                    return Err(AcidBaseError::Arithmetic(
                        "salt partition reached inventory resolution",
                    ));
                }
                let trial = evaluate(u, &mut cancelled)?;
                let residual = trial.value()?;
                if residual == 0.0 && trial.spread == 0.0 {
                    bracket = [u; 2];
                    accepted = Some((trial, iteration));
                    break;
                }
                let residual_low = residual - trial.spread;
                let residual_high = residual + trial.spread;
                if residual_low > 0.0 {
                    logs[1] = midpoint;
                } else if residual_high < 0.0 {
                    logs[0] = midpoint;
                } else {
                    // dr/d ln(u) >= the limiting component's coefficient.
                    // Both possible signs remain inside these
                    // bounds even when the inner pH bracket straddles zero.
                    logs[0] = logs[0].max(midpoint - residual_high.max(0.0) / minimum_slope);
                    logs[1] = logs[1].min(midpoint - residual_low.min(0.0) / minimum_slope);
                }
                if logs[1] - logs[0] <= options.log_activity_tolerance.value()
                    && residual.abs() + trial.spread <= options.log_saturation_tolerance.value()
                {
                    bracket = logs.map(det::exp);
                    accepted = Some((trial, iteration));
                }
            }
            let (answer, iterations) = accepted.ok_or(AcidBaseError::Budget)?;
            (answer, PrecipitationStatus::SolidPresent, iterations)
        };
        poll(&mut cancelled)?;
        let solid = limit - answer.limiting;
        if !solid.is_finite()
            || solid < 0.0
            || (status == PrecipitationStatus::SolidPresent && solid == 0.0)
        {
            return Err(AcidBaseError::Arithmetic(
                "salt solid partition loses range",
            ));
        }
        let cation_residual = component_residual(answer.metal, p * solid, cation_total.value())?;
        let family_residual = component_residual(answer.family, q * solid, family_total.value())?;
        prepared.set_family_total(self.family_index, Amount::new(answer.family))?;
        prepared.set_spectator_amount(cation_index, Amount::new(answer.metal))?;
        let aqueous = prepared.finish(answer.aqueous, &mut cancelled)?;
        Ok(SaltSpeciation {
            model: self.clone(),
            cation_total,
            family_total,
            dissolved_cation: Amount::new(answer.metal),
            dissolved_family: Amount::new(answer.family),
            solid_amount: Amount::new(solid),
            cation_balance_residual: Amount::new(cation_residual),
            family_balance_residual: Amount::new(family_residual),
            aqueous,
            log_saturation_ratio: answer.residual.map(Dimensionless::new),
            log_saturation_numerical_spread: answer
                .residual
                .map(|_| Dimensionless::new(answer.spread)),
            dissolved_limiting_bracket: bracket.map(Amount::new),
            status,
            iterations,
            aqueous_solves: solves,
            aqueous_iterations: inner_iterations,
        })
    }
}

// Choose the limiting analytical component using the original integer
// stoichiometry. Comparing rounded M/p and F/q alone can tie even when
// reconstructing the other stock would produce a negative excess. The
// two-product residuals retain the sign of p*F-q*M without an invented floor.
fn partition_basis(metal: f64, family: f64, p: u32, q: u32) -> Result<(f64, f64, bool)> {
    if p == 1 && q == 1 {
        return Ok((metal.min(family), (metal - family).abs(), metal <= family));
    }
    let (mut pf, mut qf) = (f64::from(p), f64::from(q));
    if !(pf * family).is_finite() || !(qf * metal).is_finite() {
        // A common exact binary scale prevents overflow of the comparison,
        // retaining the same stoichiometric ratio and excess component units.
        let power = 32 - p.max(q).leading_zeros();
        let scale = (1_u64 << power) as f64;
        pf /= scale;
        qf /= scale;
    }
    let left = pf * family;
    let right = qf * metal;
    let mut difference = Sum::default();
    difference.add(left)?;
    difference.add(-right)?;
    difference.add(pf.mul_add(family, -left))?;
    difference.add(-qf.mul_add(metal, -right))?;
    let difference = difference.finish()?;
    let cation_limits = difference >= 0.0;
    let (original, coefficient, excess) = if cation_limits {
        (metal, f64::from(p), difference / pf)
    } else {
        (family, f64::from(q), -difference / qf)
    };
    let limit = original / coefficient;
    if !limit.is_finite()
        || (original > 0.0 && limit == 0.0)
        || !excess.is_finite()
        || excess < 0.0
        || (difference != 0.0 && excess == 0.0)
    {
        return Err(AcidBaseError::Arithmetic(
            "positive stoichiometric salt inventory or excess loses range",
        ));
    }
    Ok((limit, excess, cation_limits))
}

struct Trial {
    limiting: f64,
    metal: f64,
    family: f64,
    residual: Option<f64>,
    spread: f64,
    aqueous: AcidBaseTrial,
}
impl Trial {
    fn value(&self) -> Result<f64> {
        self.residual
            .ok_or(AcidBaseError::Arithmetic("empty salt component activity"))
    }
}
fn component_residual(dissolved: f64, solid: f64, total: f64) -> Result<f64> {
    let mut sum = Sum::default();
    sum.add(dissolved)?;
    sum.add(solid)?;
    sum.add(-total)?;
    let residual = sum.finish()?;
    if residual.abs() > 64.0 * f64::EPSILON * total {
        return Err(AcidBaseError::Arithmetic(
            "salt analytical component loses closure",
        ));
    }
    Ok(residual)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(AcidBaseError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
