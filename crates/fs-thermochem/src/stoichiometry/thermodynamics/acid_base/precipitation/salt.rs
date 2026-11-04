//! One finite 1:1 salt coupled to an entire acid/base family and charge balance.
//!
//! The admitted reaction is `M(z+) + A(z-) + w H2O -> solid`, where A is
//! the terminal member of one original acid family and w is a nonnegative
//! integer. Pure solid and fixed solvent have unit activity. Both analytical
//! component totals include any initial solid. The other families and
//! spectator inventories are independently fixed.
//!
//! With u the smaller dissolved component, the other dissolved component is
//! u plus the excess analytical inventory. Every trial recomputes the shared
//! pH and every acid/base family. This retains positive trace concentrations
//! without subtracting a nearly complete reaction extent from a large stock.
//!
//! The terminal-anion saturation residual is strictly increasing in ln(u),
//! with derivative at least one. For an N-step family ending at charge -z,
//! its mean charge is at least -z; adding one cation and one family unit
//! therefore adds nonnegative charge at fixed pH. Charge balance raises pH,
//! increasing the terminal fraction. This supplies a bounded scalar solve.
//!
//! This is one explicitly selected solid under the existing ideal-activity,
//! fixed-solvent model. It does not select among competing solids, include
//! omitted metal complexes or redox reactions, or assert a closed H/O solvent
//! inventory. Source uncertainty is retained as evidence, not propagated into
//! the numerical pH or saturation brackets.

use super::{PrecipitationOptions, PrecipitationStatus};
use super::super::{
    AcidBaseError, AcidBaseSpeciation, ConservedReaction, FormationGibbsSource,
    IdealAcidBaseSystem, ReactionNetwork, Result, SpeciesId, SpectatorIon,
    StandardReactionGibbs, Sum, charge, consistent_sources, nonnegative, positive,
    require_phase,
};
use fs_math::det;
use fs_qty::{Amount, Dimensionless};

/// Original source-bound 1:1 salt reaction and the complete aqueous owner.
#[derive(Clone, Debug)]
pub struct SaltPrecipitation {
    aqueous: IdealAcidBaseSystem,
    reaction: ConservedReaction,
    standard: StandardReactionGibbs,
    cation: SpeciesId,
    solid: SpeciesId,
    anion: SpeciesId,
    charge: i32,
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
    /// Selected pure solid, on the same one-formula-unit basis [mol].
    pub solid_amount: Amount,
    /// Independently reduced dissolved cation + solid - cation total [mol].
    pub cation_balance_residual: Amount,
    /// Independently reduced dissolved family + solid - family total [mol].
    pub family_balance_residual: Amount,
    /// Complete aqueous result at the simultaneously solved partition and pH.
    pub aqueous: AcidBaseSpeciation,
    /// Natural log of IAP/K for dissolution; positive means supersaturation.
    /// Undefined if either analytical component is empty. This is not log10 SI.
    pub log_saturation_ratio: Option<Dimensionless>,
    /// Bound from the retained inner pH bracket on natural-log saturation
    /// variation. A numerical resolution measure, not an interval certificate.
    pub log_saturation_numerical_spread: Option<Dimensionless>,
    /// Numerical bracket on the smaller dissolved component [mol].
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
    /// Bind exactly M(z+) + A(z-) + w H2O -> neutral solid.
    ///
    /// Cation and terminal family anion have coefficients -1, solid +1,
    /// and optional water has a negative integer coefficient. The family
    /// anion is inferred from the original complete protonation chain.
    /// All shared species must retain identical original selected sources.
    pub fn from_sources(
        aqueous: IdealAcidBaseSystem,
        reaction: ConservedReaction,
        cation: SpeciesId,
        solid: SpeciesId,
        family_id: &str,
        sources: &[FormationGibbsSource<'_>],
    ) -> Result<Self> {
        let family_index = aqueous.families().iter()
            .position(|family| family.id() == family_id)
            .ok_or(AcidBaseError::Invalid("missing salt acid/base family"))?;
        let anion = aqueous.families()[family_index].species().last()
            .ok_or(AcidBaseError::Invalid("empty salt acid/base family"))?
            .clone();
        let water = aqueous.water().water();
        let coefficient = |id: &SpeciesId| {
            reaction.terms().iter().find(|term| &term.species == id)
                .map(|term| term.coefficient)
        };
        let z = charge(&reaction, &cation)?;
        let water_coefficient = coefficient(water);
        if !(1..=32).contains(&z)
            || charge(&reaction, &anion)? != -z
            || charge(&reaction, &solid)? != 0
            || coefficient(&cation) != Some(-1)
            || coefficient(&anion) != Some(-1)
            || coefficient(&solid) != Some(1)
            || !matches!(water_coefficient, None | Some(..=-1))
            || reaction.terms().len() != 3 + usize::from(water_coefficient.is_some())
            || aqueous.reaction_network().species().binary_search(&cation).is_ok()
            || aqueous.reaction_network().species().binary_search(&solid).is_ok()
        {
            return Err(AcidBaseError::Invalid(
                "requires one cation and its charge-matched terminal family anion to one neutral salt, with optional hydrate water",
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
            aqueous, reaction, standard, cation, solid, anion,
            charge: z, family_index, network,
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
    /// Analytical family that shares one component with each solid formula unit.
    pub fn family_id(&self) -> &str {
        self.aqueous.families()[self.family_index].id()
    }

    /// Solve both component balances, pH and every acid family simultaneously.
    ///
    /// The cation total and selected family total include any initial solid.
    /// Other analytical totals and spectator amounts stay fixed. The caller
    /// must not duplicate owned cation, solid, or any aqueous family member
    /// among spectators. All aqueous activities share the supplied basis.
    ///
    /// The minimum dissolved amount applies to the limiting component. A root
    /// below that caller-declared floor refuses rather than erasing a trace.
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
            amount_at_unit_activity.value(), "positive common aqueous basis required",
        )?;
        let floor = positive(
            options.minimum_dissolved_amount.value(), "positive dissolved search floor required",
        )?;
        if [options.log_activity_tolerance, options.log_saturation_tolerance]
            .iter().any(|v| !v.value().is_finite() || v.value() <= 0.0)
            || !(1..=256).contains(&options.max_iterations)
            || spectators.iter().any(|ion| ion.species == self.cation || ion.species == self.solid)
        {
            return Err(AcidBaseError::Invalid(
                "invalid salt controls or duplicated owned participant",
            ));
        }
        let mut selected = family_totals.iter().enumerate()
            .filter(|(_, (id, _))| id == self.family_id());
        let (input_index, (_, family_total)) = selected.next()
            .ok_or(AcidBaseError::Invalid("missing salt analytical total"))?;
        if selected.next().is_some() {
            return Err(AcidBaseError::Invalid("duplicate salt analytical total"));
        }
        let family_total = *family_total;
        nonnegative(family_total.value(), "invalid salt family total")?;
        let limit = cation_total.value().min(family_total.value());
        let excess = (cation_total.value() - family_total.value()).abs();
        let cation_limits = cation_total.value() <= family_total.value();
        let mut families = family_totals.to_vec();
        let mut ions = spectators.to_vec();
        ions.push(SpectatorIon {
            species: self.cation.clone(), charge: self.charge, amount: cation_total,
        });
        let mut solves = 0;
        let mut inner_iterations = 0;
        let mut evaluate = |u: f64, cancelled: &mut _| -> Result<Trial> {
            poll(cancelled)?;
            // Preserve the supplied upper endpoint bit-for-bit: reconstructing
            // it as excess + limit can round differently from the caller total.
            let (metal, family) = if u == limit {
                (cation_total.value(), family_total.value())
            } else if cation_limits {
                (u, excess + u)
            } else {
                (excess + u, u)
            };
            if !metal.is_finite() || !family.is_finite()
                || metal < 0.0 || family < 0.0
                || metal > cation_total.value() || family > family_total.value()
            {
                return Err(AcidBaseError::Arithmetic(
                    "salt component partition loses inventory range",
                ));
            }
            families[input_index].1 = Amount::new(family);
            ions.last_mut().ok_or(AcidBaseError::Invalid("missing owned salt cation"))?
                .amount = Amount::new(metal);
            let aqueous = self.aqueous.speciate(
                &families, &ions, amount_at_unit_activity, options.aqueous, cancelled,
            )?;
            solves += 1;
            inner_iterations += aqueous.iterations;
            let residual = if u == 0.0 {
                None
            } else {
                let anion = aqueous.families[self.family_index].forms.last()
                    .ok_or(AcidBaseError::Invalid("missing terminal salt anion"))?;
                let mut sum = Sum::default();
                sum.add(det::ln(metal))?;
                sum.add(-det::ln(basis))?;
                sum.add(det::ln(anion.activity.value()))?;
                sum.add(self.standard.log_equilibrium_constant().value())?;
                Some(sum.finish()?)
            };
            // d ln(alpha_terminal)/d pH = ln(10) * (N - mean(step)),
            // bounded by N*ln(10). Keep this inner-solve uncertainty in every
            // phase/sign decision and in the final saturation criterion.
            let spread = self.aqueous.families()[self.family_index].steps().len() as f64
                * core::f64::consts::LN_10
                * aqueous.ph_bracket.iter()
                    .map(|p| (p.value() - aqueous.ph.value()).abs())
                    .fold(0.0, f64::max);
            if !spread.is_finite() {
                return Err(AcidBaseError::Arithmetic(
                    "salt aqueous pH resolution loses range",
                ));
            }
            Ok(Trial { limiting: u, metal, family, residual, spread, aqueous })
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
                    // dr/d ln(u) >= 1. Both possible signs remain inside these
                    // bounds even when the inner pH bracket straddles zero.
                    logs[0] = logs[0].max(midpoint - residual_high.max(0.0));
                    logs[1] = logs[1].min(midpoint - residual_low.min(0.0));
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
        if !solid.is_finite() || solid < 0.0
            || (status == PrecipitationStatus::SolidPresent && solid == 0.0)
        {
            return Err(AcidBaseError::Arithmetic("salt solid partition loses range"));
        }
        let cation_residual = component_residual(answer.metal, solid, cation_total.value())?;
        let family_residual = component_residual(answer.family, solid, family_total.value())?;
        Ok(SaltSpeciation {
            model: self.clone(),
            cation_total, family_total,
            dissolved_cation: Amount::new(answer.metal),
            dissolved_family: Amount::new(answer.family),
            solid_amount: Amount::new(solid),
            cation_balance_residual: Amount::new(cation_residual),
            family_balance_residual: Amount::new(family_residual),
            aqueous: answer.aqueous,
            log_saturation_ratio: answer.residual.map(Dimensionless::new),
            log_saturation_numerical_spread: answer.residual.map(|_| Dimensionless::new(answer.spread)),
            dissolved_limiting_bracket: bracket.map(Amount::new),
            status, iterations,
            aqueous_solves: solves, aqueous_iterations: inner_iterations,
        })
    }
}

struct Trial {
    limiting: f64,
    metal: f64,
    family: f64,
    residual: Option<f64>,
    spread: f64,
    aqueous: AcidBaseSpeciation,
}
impl Trial {
    fn value(&self) -> Result<f64> {
        self.residual.ok_or(AcidBaseError::Arithmetic("empty salt component activity"))
    }
}
fn component_residual(dissolved: f64, solid: f64, total: f64) -> Result<f64> {
    let mut sum = Sum::default();
    sum.add(dissolved)?;
    sum.add(solid)?;
    sum.add(-total)?;
    let residual = sum.finish()?;
    if residual.abs() > 64.0 * f64::EPSILON * total {
        return Err(AcidBaseError::Arithmetic("salt analytical component loses closure"));
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
