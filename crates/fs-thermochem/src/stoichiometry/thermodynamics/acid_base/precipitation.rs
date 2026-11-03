//! One hydrolyzing mineral phase coupled to aqueous charge and proton balances.
//!
//! For `M(z+) + w H2O -> solid + z H+`, the pure solid and solvent have unit
//! activity. The ideal aqueous partition satisfies `a_M = a_H^z / K_precip`
//! when solid is present, and the corresponding undersaturation inequality
//! when the entire finite component stays dissolved. All acid/base families
//! share the proton activity determined by electroneutrality. Proton activity
//! is never held at the value calculated before precipitation.
//!
//! This is the fixed-solvent, ideal-activity model of the original acid/base
//! owner, extended by ONE explicitly selected solid. Other metal complexes,
//! competing solids, nonideal activities, nucleation, rates and solvent H/O
//! inventories require other models. USGS describes pure-phase equilibrium
//! and absence inequalities in
//! <https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-13.htm>.

use super::*;
use fs_math::det;

/// Source-bound, exactly conserved hydrolysis phase and its aqueous owner.
#[derive(Clone, Debug)]
pub struct HydrolysisPrecipitation {
    aqueous: IdealAcidBaseSystem,
    reaction: ConservedReaction,
    standard: StandardReactionGibbs,
    cation: SpeciesId,
    solid: SpeciesId,
    charge: i32,
    network: ReactionNetwork,
}

/// Numerical limits for the outer mineral partition and inner charge solve.
#[derive(Clone, Copy, Debug)]
pub struct PrecipitationOptions {
    /// Original bounded pH/charge solve. Its accuracy must resolve the desired
    /// mineral residual; an inadequate inner tolerance does not earn success.
    pub aqueous: SpeciationOptions,
    /// Positive lower search bound on dissolved component [mol]. A solution
    /// below this caller-declared bound refuses instead of erasing trace metal.
    pub minimum_dissolved_amount: Amount,
    /// Maximum final natural-log dissolved-activity bracket width.
    pub log_activity_tolerance: Dimensionless,
    /// Maximum absolute natural-log saturation residual for a present solid.
    pub log_saturation_tolerance: Dimensionless,
    /// At most 256 outer bisections, each with the bounded aqueous solve above.
    pub max_iterations: usize,
}
impl Default for PrecipitationOptions {
    fn default() -> Self {
        Self {
            aqueous: SpeciationOptions {
                ph_tolerance: Dimensionless::new(1e-12),
                absolute_charge_tolerance: Amount::new(1e-15),
                relative_charge_tolerance: Dimensionless::new(1e-12),
                ..SpeciationOptions::default()
            },
            minimum_dissolved_amount: Amount::new(1e-100),
            log_activity_tolerance: Dimensionless::new(1e-9),
            log_saturation_tolerance: Dimensionless::new(1e-9),
            max_iterations: 256,
        }
    }
}

/// Equilibrium along the declared mineral/acid-base reaction set only.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PrecipitationStatus {
    /// The finite component is all dissolved; saturation ratio is at most one.
    SolidAbsent,
    /// Positive solid and dissolved amounts coexist at the requested tolerance.
    SolidPresent,
    /// No analytical component was supplied; its saturation is undefined.
    EmptyComponent,
}

/// Coupled ideal solution and finite solid partition, with original evidence.
#[derive(Clone, Debug)]
pub struct PrecipitationSpeciation {
    model: HydrolysisPrecipitation,
    /// Total of dissolved cation and this solid, each on a one-formula-unit
    /// basis. It is not an inferred total across omitted complexes or phases.
    pub component_total: Amount,
    /// Positive trace concentrations remain representable independently of the
    /// much larger solid stock; they are not obtained by subtracting an extent.
    pub dissolved_amount: Amount,
    /// Nonnegative amount of the selected solid [mol].
    pub solid_amount: Amount,
    /// Independently reduced `dissolved + solid - component_total` [mol].
    pub component_balance_residual: Amount,
    /// Full simultaneous aqueous result; its spectator list includes the
    /// outer-solved cation in addition to the caller's fixed spectator stocks.
    pub aqueous: AcidBaseSpeciation,
    /// Natural log of IAP/K for the DISSOLUTION equation. Positive means
    /// supersaturation; absent-phase results must be <= 0. Not log10(SI).
    pub log_saturation_ratio: Option<Dimensionless>,
    /// Maximum natural-log saturation variation over the retained inner pH
    /// bracket. Included in the convergence check, not a rigorous enclosure.
    pub log_saturation_numerical_spread: Option<Dimensionless>,
    /// Final numerical dissolved-amount bracket; no interval certificate.
    pub dissolved_bracket: [Amount; 2],
    /// Present, absent or an empty component.
    pub status: PrecipitationStatus,
    /// Outer mineral bisections, excluding endpoint evaluations.
    pub iterations: usize,
    /// Number of complete bounded aqueous solves used.
    pub aqueous_solves: usize,
    /// Sum of the inner pH bisection counts, including endpoint solves.
    pub aqueous_iterations: usize,
}
impl PrecipitationSpeciation {
    /// Original phase equation, formation sources and aqueous-family owners.
    pub fn model(&self) -> &HydrolysisPrecipitation {
        &self.model
    }
}

impl HydrolysisPrecipitation {
    /// Bind `M(z+) + w H2O -> solid + z H+`, with exactly four species,
    /// cation/solid coefficients -1/+1, positive integer w and 1 <= z <= 32.
    /// Atom/charge identities stay in the original reaction/network owner.
    /// All shared species must retain IDENTICAL original selected sources.
    pub fn from_sources(
        aqueous: IdealAcidBaseSystem,
        reaction: ConservedReaction,
        cation: SpeciesId,
        solid: SpeciesId,
        sources: &[FormationGibbsSource<'_>],
    ) -> Result<Self> {
        let water = aqueous.water();
        let coefficient = |id: &SpeciesId| {
            reaction
                .terms()
                .iter()
                .find(|term| &term.species == id)
                .map(|term| term.coefficient)
        };
        let z = charge(&reaction, &cation)?;
        if !(1..=32).contains(&z)
            || charge(&reaction, &solid)? != 0
            || reaction.terms().len() != 4
            || coefficient(&cation) != Some(-1)
            || coefficient(&solid) != Some(1)
            || coefficient(water.proton()) != Some(i128::from(z))
            || !matches!(coefficient(water.water()), Some(n) if n < 0)
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
                "requires distinct cation/water to one solid plus charge-matched protons",
            ));
        }
        let standard = reaction.standard_gibbs_change(sources)?;
        for source in sources {
            let (phase, pure) = if *source.species == solid {
                ("source-phase-solid", true)
            } else if source.species == water.water() {
                ("source-phase-liquid", true)
            } else {
                ("source-phase-aqueous", false)
            };
            require_phase(source, &standard, phase, pure)?;
        }
        let mut standards = vec![&standard, water.standard()];
        let mut reactions = vec![reaction.clone(), water.reaction().clone()];
        for family in aqueous.families() {
            for step in family.steps() {
                standards.push(step.standard());
                reactions.push(step.reaction().clone());
            }
        }
        consistent_sources(standards)?;
        let network = ReactionNetwork::new(reactions)?;
        Ok(Self {
            aqueous,
            reaction,
            standard,
            cation,
            solid,
            charge: z,
            network,
        })
    }

    /// Original fixed-solvent acid/base owner and its source references.
    pub fn aqueous_model(&self) -> &IdealAcidBaseSystem {
        &self.aqueous
    }
    /// Original exactly conserved precipitation equation.
    pub fn reaction(&self) -> &ConservedReaction {
        &self.reaction
    }
    /// Original formation Gibbs evidence for the precipitation direction.
    pub fn standard(&self) -> &StandardReactionGibbs {
        &self.standard
    }
    /// All original equations, with shared atom/charge definitions checked.
    pub fn reaction_network(&self) -> &ReactionNetwork {
        &self.network
    }
    /// Explicit dissolved cation identity; its charge comes from the equation.
    pub fn cation(&self) -> &SpeciesId {
        &self.cation
    }
    /// Selected pure solid identity, not an automatic stable-phase choice.
    pub fn solid(&self) -> &SpeciesId {
        &self.solid
    }

    /// Solve pH, every acid/base family and one finite solid partition together.
    /// `component_total` includes all initial dissolved cation plus this solid.
    /// Family totals and other ions are independently supplied and fixed; the
    /// cation must NOT also appear among caller spectators or family members.
    /// The fixed solvent basis is common to all aqueous activities.
    ///
    /// For increasing dissolved cation, charge balance increases pH and reduces
    /// proton activity. Thus ln(Q_precip/K_precip) decreases strictly. The outer
    /// solve uses log dissolved amount to retain trace dissolved metal even
    /// when a reaction-extent subtraction could not resolve it.
    #[allow(clippy::too_many_arguments)]
    pub fn speciate(
        &self,
        component_total: Amount,
        family_totals: &[(String, Amount)],
        spectators: &[SpectatorIon],
        amount_at_unit_activity: Amount,
        options: PrecipitationOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PrecipitationSpeciation> {
        poll(&mut cancelled)?;
        nonnegative(component_total.value(), "invalid mineral component total")?;
        let basis = positive(
            amount_at_unit_activity.value(),
            "positive common aqueous basis required",
        )?;
        let floor = positive(
            options.minimum_dissolved_amount.value(),
            "positive dissolved search floor required",
        )?;
        for tolerance in [
            options.log_activity_tolerance,
            options.log_saturation_tolerance,
        ] {
            if !tolerance.value().is_finite() || tolerance.value() <= 0.0 {
                return Err(AcidBaseError::Invalid(
                    "invalid mineral logarithmic tolerance",
                ));
            }
        }
        if !(1..=256).contains(&options.max_iterations)
            || spectators
                .iter()
                .any(|ion| ion.species == self.cation || ion.species == self.solid)
        {
            return Err(AcidBaseError::Invalid(
                "invalid work budget or duplicated phase participant",
            ));
        }
        let mut inputs = spectators.to_vec();
        inputs.push(SpectatorIon {
            species: self.cation.clone(),
            charge: self.charge,
            amount: component_total,
        });
        let mut solves = 0;
        let mut inner_iterations = 0;
        let mut evaluate = |amount: f64, cancelled: &mut _| -> Result<Trial> {
            poll(cancelled)?;
            inputs
                .last_mut()
                .ok_or(AcidBaseError::Invalid("missing owned cation"))?
                .amount = Amount::new(amount);
            let aqueous = self.aqueous.speciate(
                family_totals,
                &inputs,
                amount_at_unit_activity,
                options.aqueous,
                cancelled,
            )?;
            solves += 1;
            inner_iterations += aqueous.iterations;
            let residual = if amount == 0.0 {
                None
            } else {
                let log_activity = det::ln(amount) - det::ln(basis);
                let residual = f64::from(self.charge) * det::ln(aqueous.proton_activity.value())
                    - log_activity
                    - self.standard.log_equilibrium_constant().value();
                if !residual.is_finite() {
                    return Err(AcidBaseError::Arithmetic("mineral residual loses range"));
                }
                Some(residual)
            };
            let spread = f64::from(self.charge)
                * core::f64::consts::LN_10
                * aqueous
                    .ph_bracket
                    .iter()
                    .map(|p| (p.value() - aqueous.ph.value()).abs())
                    .fold(0.0, f64::max);
            if !spread.is_finite() {
                return Err(AcidBaseError::Arithmetic(
                    "aqueous pH resolution loses range",
                ));
            }
            Ok(Trial {
                amount,
                residual,
                spread,
                aqueous,
            })
        };
        let upper = evaluate(component_total.value(), &mut cancelled)?;
        let mut bracket = [component_total.value(); 2];
        let (answer, status, iterations) = if component_total.value() == 0.0 {
            (upper, PrecipitationStatus::EmptyComponent, 0)
        } else if upper.value()? - upper.spread >= 0.0 {
            (upper, PrecipitationStatus::SolidAbsent, 0)
        } else {
            if upper.value()? + upper.spread >= 0.0 {
                return Err(AcidBaseError::Arithmetic(
                    "aqueous pH bracket does not resolve phase absence",
                ));
            }
            if floor >= component_total.value() {
                return Err(AcidBaseError::Invalid(
                    "dissolved search floor excludes a precipitating partition",
                ));
            }
            let lower = evaluate(floor, &mut cancelled)?;
            if lower.value()? + lower.spread < 0.0 {
                return Err(AcidBaseError::OutsideBracket);
            }
            if lower.value()? - lower.spread < 0.0 {
                return Err(AcidBaseError::Arithmetic(
                    "aqueous pH bracket does not resolve dissolved search floor",
                ));
            }
            let mut logs = [det::ln(floor), det::ln(component_total.value())];
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
                        "mineral partition reached logarithmic resolution",
                    ));
                }
                let amount = positive(det::exp(midpoint), "dissolved amount loses range")?;
                if amount <= floor || amount >= component_total.value() {
                    return Err(AcidBaseError::Arithmetic(
                        "mineral partition reached inventory resolution",
                    ));
                }
                let trial = evaluate(amount, &mut cancelled)?;
                let residual = trial.value()?;
                if residual == 0.0 && trial.spread == 0.0 {
                    bracket = [amount; 2];
                    accepted = Some((trial, iteration));
                    break;
                }
                let residual_low = residual - trial.spread;
                let residual_high = residual + trial.spread;
                if residual_low > 0.0 {
                    logs[0] = midpoint;
                } else if residual_high < 0.0 {
                    logs[1] = midpoint;
                } else {
                    // For u=ln(n_M), dr/du = -1 - z*ln(10)*d(pH)/du <= -1.
                    // An uncertain sign can still narrow the root interval by
                    // its residual range without discarding either sign.
                    logs[0] = logs[0].max(midpoint + residual_low.min(0.0));
                    logs[1] = logs[1].min(midpoint + residual_high.max(0.0));
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
        let solid = component_total.value() - answer.amount;
        if !solid.is_finite()
            || solid < 0.0
            || (status == PrecipitationStatus::SolidPresent && solid == 0.0)
        {
            return Err(AcidBaseError::Arithmetic("solid partition loses range"));
        }
        let mut balance = Sum::default();
        balance.add(answer.amount)?;
        balance.add(solid)?;
        balance.add(-component_total.value())?;
        let residual = balance.finish()?;
        if residual.abs() > 64.0 * f64::EPSILON * component_total.value() {
            return Err(AcidBaseError::Arithmetic(
                "mineral component balance loses closure",
            ));
        }
        Ok(PrecipitationSpeciation {
            model: self.clone(),
            component_total,
            dissolved_amount: Amount::new(answer.amount),
            solid_amount: Amount::new(solid),
            component_balance_residual: Amount::new(residual),
            aqueous: answer.aqueous,
            log_saturation_ratio: answer.residual.map(|v| Dimensionless::new(-v)),
            log_saturation_numerical_spread: answer
                .residual
                .map(|_| Dimensionless::new(answer.spread)),
            dissolved_bracket: bracket.map(Amount::new),
            status,
            iterations,
            aqueous_solves: solves,
            aqueous_iterations: inner_iterations,
        })
    }
}

struct Trial {
    amount: f64,
    residual: Option<f64>,
    spread: f64,
    aqueous: AcidBaseSpeciation,
}
impl Trial {
    fn value(&self) -> Result<f64> {
        self.residual
            .ok_or(AcidBaseError::Arithmetic("empty mineral activity"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(AcidBaseError::Cancelled)
    } else {
        Ok(())
    }
}
