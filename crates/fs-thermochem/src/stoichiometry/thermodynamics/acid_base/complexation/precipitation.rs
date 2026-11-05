//! One selected Cu2S phase in the original finite cyanide/pH system.
//!
//! The exact phase equation is `2 Cu+ + HS- -> Cu2S + H+`. Pure solid has
//! unit activity; its dissolution saturation residual is
//! `r = 2 ln(aCu+) + ln(aHS-) - ln(aH+) + ln(K_precip)`.
//! Positive solid requires r=0; an absent solid requires r<=0. The original
//! `IdealLigandSystem` recomputes every copper complex, free cyanide and pH
//! for every trial partition. No sequential frozen-pH approximation is used.
//!
//! Let u=min(Cu_dissolved/2, S_dissolved). The excess analytical component is
//! retained separately, so traces never come from total minus precipitate.
//! The reaction direction adds 2Cu and 1S to the dissolved component vector.
//! The ideal aqueous Gibbs minimum is convex in its analytical components;
//! its derivative in this direction is the phase chemical-potential sum.
//! Thus r increases with u. Bisection retains sign-resolved endpoints and
//! includes the inner ligand/pH numerical spread in every sign decision.
//!
//! This is a selected H2S/HS- and copper-cyanide assemblage, with a conservative
//! admitted final pH interval [2,10]. This domain is an explicit approximation
//! policy, not an experimental accuracy bound. S(2-), copper hydrosulfide
//! complexes, CuCN and other solids, redox, gas loss, nucleation and rates are
//! absent. The source-bound ideal activity and fixed solvent conventions of
//! the original owner remain unchanged. There is no full solvent H/O ledger.

use super::super::precipitation::PrecipitationStatus;
use super::*;

/// Conservative selected-assemblage domain, not a measured validity interval.
pub const COPPER_SULFIDE_PH_DOMAIN: [f64; 2] = [2.0, 10.0];

/// Original copper/ligand owner and exactly one proton-coupled Cu2S phase.
#[derive(Clone, Debug)]
pub struct CopperSulfidePrecipitation {
    ligand: IdealLigandSystem,
    reaction: ConservedReaction,
    standard: StandardReactionGibbs,
    copper_index: usize,
    sulfide_index: usize,
    solid: SpeciesId,
    network: ReactionNetwork,
}

/// Bounded phase search and original bounded inner speciation controls.
#[derive(Clone, Copy, Debug)]
pub struct CopperSulfideOptions {
    /// Original solver controls, including its charge and ligand residuals.
    pub ligand: LigandSpeciationOptions,
    /// Positive lower search bound on min(Cu_dissolved/2,S_dissolved), in mol.
    /// A root below this floor refuses instead of replacing a trace by zero.
    pub minimum_dissolved_amount: Amount,
    /// Maximum final natural-log width of the limiting-component bracket.
    pub log_dissolved_tolerance: Dimensionless,
    /// Maximum absolute ln(IAP/K), including inner numerical spread.
    pub log_saturation_tolerance: Dimensionless,
    /// At most256 outer rounds; ambiguous rounds may evaluate both quarters.
    pub max_iterations: usize,
    /// At most1024 complete original ligand-system solves, including endpoints.
    pub max_ligand_solves: usize,
}
impl Default for CopperSulfideOptions {
    fn default() -> Self {
        Self {
            ligand: LigandSpeciationOptions {
                log_ligand_tolerance: Dimensionless::new(1e-10),
                ..LigandSpeciationOptions::default()
            },
            minimum_dissolved_amount: Amount::new(1e-100),
            log_dissolved_tolerance: Dimensionless::new(1e-8),
            log_saturation_tolerance: Dimensionless::new(1e-7),
            max_iterations: 192,
            max_ligand_solves: 512,
        }
    }
}

/// Complete finite Cu/S/CN partition, with the original source owners retained.
#[derive(Clone, Debug)]
pub struct CopperSulfideSpeciation {
    model: CopperSulfidePrecipitation,
    /// Analytical Cu atoms, including two per supplied solid formula unit [mol].
    pub copper_total: Amount,
    /// Analytical sulfide S, including one per supplied solid formula unit [mol].
    pub sulfide_total: Amount,
    /// Sum of all original dissolved Cu(I) family members [mol].
    pub dissolved_copper: Amount,
    /// Sum of the original dissolved H2S and HS- amounts [mol].
    pub dissolved_sulfide: Amount,
    /// Selected pure Cu2S formula units [mol].
    pub solid_amount: Amount,
    /// Independently reduced sum(Cu members)+2*solid-total [mol].
    pub copper_balance_residual: Amount,
    /// Independently reduced H2S+HS-+solid-total [mol].
    pub sulfide_balance_residual: Amount,
    /// Complete original simultaneous aqueous, metal and ligand state.
    /// Its ligand balance includes all dissolved cyanide; Cu2S contains none.
    pub ligand: LigandSpeciation,
    /// Natural log of dissolution IAP/K; undefined if either component is zero.
    pub log_saturation_ratio: Option<Dimensionless>,
    /// Retained inner-resolution spread, not a rigorous interval certificate.
    pub log_saturation_numerical_spread: Option<Dimensionless>,
    /// Numerical bracket on min(Cu_dissolved/2,S_dissolved) [mol].
    pub dissolved_limiting_bracket: [Amount; 2],
    /// Absent, present, or no supplied Cu/S component.
    pub status: PrecipitationStatus,
    /// Outer rounds, excluding endpoint evaluations.
    pub iterations: usize,
    /// Complete original ligand solves used.
    pub ligand_solves: usize,
    /// All original inner charge solves used, including discarded trials.
    pub aqueous_solves: usize,
}
impl CopperSulfideSpeciation {
    /// Original source-bound phase and all original aqueous owners.
    pub fn model(&self) -> &CopperSulfidePrecipitation {
        &self.model
    }
}

impl CopperSulfidePrecipitation {
    /// Bind exactly `2 Cu+ + HS- -> Cu2S + H+` to the named Cu(I) and H2S/HS-
    /// families. The ligand must be CN-, and the copper family must start with
    /// uncomplexed Cu+. Exact atoms, charges, phase/source conventions and
    /// shared original Gibbs receipts are checked; names infer no chemistry.
    pub fn from_sources(
        ligand: IdealLigandSystem,
        reaction: ConservedReaction,
        copper_family: &str,
        sulfide_family: &str,
        solid: SpeciesId,
        sources: &[FormationGibbsSource<'_>],
    ) -> Result<Self> {
        let copper_index = ligand
            .families()
            .iter()
            .position(|f| f.id() == copper_family)
            .ok_or(LigandSpeciationError::Invalid(
                "missing original copper family",
            ))?;
        let sulfide_index = ligand
            .aqueous_model()
            .families()
            .iter()
            .position(|f| f.id() == sulfide_family)
            .ok_or(LigandSpeciationError::Invalid(
                "missing original sulfide family",
            ))?;
        let copper = &ligand.families()[copper_index];
        let sulfide = &ligand.aqueous_model().families()[sulfide_index];
        if sulfide.steps().len() != 1
            || sulfide_family == ligand.ligand_family_id()
            || ligand
                .reaction_network()
                .species()
                .binary_search(&solid)
                .is_ok()
        {
            return Err(LigandSpeciationError::Invalid(
                "require one disjoint H2S/HS- family and a new solid",
            ));
        }
        let proton = ligand.aqueous_model().water().proton();
        let cu = &copper.species()[0];
        let hs = &sulfide.species()[1];
        require_terms(&reaction, &[(cu, -2), (hs, -1), (&solid, 1), (proton, 1)])?;
        require_formula(&reaction, cu, &[("Cu", 1)], 1)?;
        require_formula(&reaction, hs, &[("H", 1), ("S", 1)], -1)?;
        require_formula(&reaction, &solid, &[("Cu", 2), ("S", 1)], 0)?;
        require_formula(&reaction, proton, &[("H", 1)], 1)?;
        require_formula(
            sulfide.steps()[0].reaction(),
            &sulfide.species()[0],
            &[("H", 2), ("S", 1)],
            0,
        )?;
        require_formula(
            copper.steps()[0].reaction(),
            copper.ligand(),
            &[("C", 1), ("N", 1)],
            -1,
        )?;
        let standard = reaction.standard_gibbs_change(sources)?;
        for source in sources {
            require_phase(
                source,
                &standard,
                if *source.species == solid {
                    "source-phase-solid"
                } else {
                    "source-phase-aqueous"
                },
                *source.species == solid,
            )?;
        }
        let mut standards = vec![&standard, ligand.aqueous_model().water().standard()];
        for family in ligand.aqueous_model().families() {
            standards.extend(family.steps().iter().map(|s| s.standard()));
        }
        for family in ligand.families() {
            standards.extend(family.steps().iter().map(|s| s.standard()));
        }
        consistent_sources(standards)?;
        let mut reactions = ligand.reaction_network().reactions().to_vec();
        reactions.push(reaction.clone());
        let network = ReactionNetwork::new(reactions)?;
        Ok(Self {
            ligand,
            reaction,
            standard,
            copper_index,
            sulfide_index,
            solid,
            network,
        })
    }
    /// The existing complete finite metal/ligand and acid/base model.
    pub fn ligand_model(&self) -> &IdealLigandSystem {
        &self.ligand
    }
    /// Original one-formula-unit precipitation equation.
    pub fn reaction(&self) -> &ConservedReaction {
        &self.reaction
    }
    /// Selected direct-G sources and the original precipitation ln K.
    pub fn standard(&self) -> &StandardReactionGibbs {
        &self.standard
    }
    /// Original shared-species network, including the selected phase.
    pub fn reaction_network(&self) -> &ReactionNetwork {
        &self.network
    }
    /// Original selected solid identity; its phase reference is not renamed.
    pub fn solid(&self) -> &SpeciesId {
        &self.solid
    }
    /// Original Cu(I) analytical family identifier.
    pub fn copper_family_id(&self) -> &str {
        self.ligand.families()[self.copper_index].id()
    }
    /// Original truncated H2S/HS- analytical family identifier.
    pub fn sulfide_family_id(&self) -> &str {
        self.ligand.aqueous_model().families()[self.sulfide_index].id()
    }

    /// Solve finite Cu/S partition, copper complexes, cyanide and charge together.
    ///
    /// The named Cu and sulfide totals include any initial Cu2S solid in the
    /// ratio 2:1. Other metal/acid totals and spectators stay fixed. Supply all
    /// original totals exactly once. The selected solid cannot be a spectator.
    /// The entire final pH bracket must lie in [2,10]; no partial result survives
    /// any source, arithmetic, budget, domain or cancellation refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn speciate(
        &self,
        ligand_total: Amount,
        metal_totals: &[(String, Amount)],
        other_acid_totals: &[(String, Amount)],
        spectators: &[SpectatorIon],
        amount_at_unit_activity: Amount,
        options: CopperSulfideOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CopperSulfideSpeciation> {
        poll(&mut cancelled)?;
        let floor = require_positive(
            options.minimum_dissolved_amount.value(),
            "positive dissolved floor required",
        )?;
        require_positive(
            options.log_dissolved_tolerance.value(),
            "positive phase bracket tolerance required",
        )?;
        require_positive(
            options.log_saturation_tolerance.value(),
            "positive saturation tolerance required",
        )?;
        if !(1..=256).contains(&options.max_iterations)
            || !(1..=1024).contains(&options.max_ligand_solves)
            || spectators.iter().any(|s| s.species == self.solid)
        {
            return Err(LigandSpeciationError::Invalid(
                "invalid phase work limits or duplicated solid",
            ));
        }
        let copper_slot = unique_total(metal_totals, self.copper_family_id())?;
        let sulfide_slot = unique_total(other_acid_totals, self.sulfide_family_id())?;
        let copper_total = metal_totals[copper_slot].1.value();
        let sulfide_total = other_acid_totals[sulfide_slot].1.value();
        let (limit, excess, copper_limits) = partition(copper_total, sulfide_total)?;
        let mut metals = metal_totals.to_vec();
        let mut acids = other_acid_totals.to_vec();
        let mut ligand_solves = 0;
        let mut aqueous_solves = 0;
        let mut evaluate = |u: f64, cancelled: &mut _| -> Result<Trial> {
            poll(cancelled)?;
            if ligand_solves >= options.max_ligand_solves {
                return Err(LigandSpeciationError::Budget);
            }
            let (cu, sulfur) = if u == limit {
                (copper_total, sulfide_total)
            } else if copper_limits {
                (2.0 * u, excess + u)
            } else {
                (excess + 2.0 * u, u)
            };
            if !cu.is_finite()
                || !sulfur.is_finite()
                || cu < 0.0
                || sulfur < 0.0
                || cu > copper_total
                || sulfur > sulfide_total
                || (u > 0.0 && (cu == 0.0 || sulfur == 0.0))
            {
                return Err(LigandSpeciationError::Arithmetic(
                    "dissolved phase components lose range",
                ));
            }
            metals[copper_slot].1 = Amount::new(cu);
            acids[sulfide_slot].1 = Amount::new(sulfur);
            let state = self.ligand.speciate(
                ligand_total,
                &metals,
                &acids,
                spectators,
                amount_at_unit_activity,
                options.ligand,
                &mut *cancelled,
            )?;
            ligand_solves += 1;
            aqueous_solves += state.aqueous_solves;
            let (residual, spread) = if u == 0.0 {
                (None, 0.0)
            } else {
                let cu_activity = computed_positive(
                    state.metals[self.copper_index].members[0].activity.value(),
                    "free copper trace loses range",
                )?;
                let hs_activity = computed_positive(
                    state.aqueous.families[self.sulfide_index].forms[1]
                        .activity
                        .value(),
                    "hydrosulfide trace loses range",
                )?;
                let h_activity = computed_positive(
                    state.aqueous.proton_activity.value(),
                    "proton trace loses range",
                )?;
                let r = add([
                    2.0 * det::ln(cu_activity),
                    det::ln(hs_activity),
                    -det::ln(h_activity),
                    self.standard.log_equilibrium_constant().value(),
                ])?;
                let free = state.free_ligand.value();
                let ligand_spread = if free == 0.0 {
                    0.0
                } else {
                    let log_free = det::ln(free);
                    state
                        .free_ligand_bracket
                        .iter()
                        .map(|n| (det::ln(n.value()) - log_free).abs())
                        .fold(0.0, f64::max)
                };
                let ph_spread = state
                    .aqueous
                    .ph_bracket
                    .iter()
                    .map(|p| (p.value() - state.aqueous.ph.value()).abs())
                    .fold(0.0, f64::max);
                let max_ligands = *self.ligand.families()[self.copper_index]
                    .ligand_counts()
                    .last()
                    .unwrap() as f64;
                // |d ln(free Cu)/d ln(free CN)| <= maximum ligand count.
                // d ln(HS/H)/d pH lies between ln(10) and 2ln(10).
                let spread = add([
                    2.0 * max_ligands * ligand_spread,
                    2.0 * LN_10 * ph_spread,
                    2.0 * state.max_log_mass_action_residual.value(),
                ])?;
                (Some(r), spread)
            };
            Ok(Trial {
                u,
                state,
                residual,
                spread,
            })
        };
        // Upper endpoint also admits every original total, identity and control.
        let upper = evaluate(limit, &mut cancelled)?;
        let mut bracket = [limit; 2];
        let (answer, status, iterations) = if limit == 0.0 {
            (upper, PrecipitationStatus::EmptyComponent, 0)
        } else if upper.value()? + upper.spread <= 0.0 {
            (upper, PrecipitationStatus::SolidAbsent, 0)
        } else {
            if upper.value()? - upper.spread <= 0.0 {
                return Err(LigandSpeciationError::Arithmetic(
                    "inner resolution cannot decide phase absence",
                ));
            }
            if floor >= limit {
                return Err(LigandSpeciationError::OutsideBracket);
            }
            let lower = evaluate(floor, &mut cancelled)?;
            if lower.value()? - lower.spread > 0.0 {
                return Err(LigandSpeciationError::OutsideBracket);
            }
            if lower.value()? + lower.spread > 0.0 {
                return Err(LigandSpeciationError::Arithmetic(
                    "inner resolution cannot resolve phase search floor",
                ));
            }
            let mut logs = [det::ln(floor), det::ln(limit)];
            let mut accepted = None;
            for iteration in 1..=options.max_iterations {
                poll(&mut cancelled)?;
                let old = logs;
                let midpoint = logs[0].midpoint(logs[1]);
                if midpoint <= logs[0] || midpoint >= logs[1] {
                    return Err(LigandSpeciationError::Arithmetic(
                        "phase bracket reaches floating resolution",
                    ));
                }
                let u = exp_positive(midpoint, "phase trace loses range")?;
                let trial = evaluate(u, &mut cancelled)?;
                let residual = trial.value()?;
                if logs[1] - logs[0] <= options.log_dissolved_tolerance.value()
                    && residual.abs() + trial.spread <= options.log_saturation_tolerance.value()
                {
                    bracket = logs.map(det::exp);
                    accepted = Some((trial, iteration));
                    break;
                }
                if residual - trial.spread > 0.0 {
                    logs[1] = midpoint;
                } else if residual + trial.spread < 0.0 {
                    logs[0] = midpoint;
                } else {
                    // No unproved slope bound: probe each half and retain only
                    // sign-resolved endpoints. An unresolved inner solve cannot
                    // silently shrink the outer root out of its bracket.
                    for probe in [logs[0].midpoint(midpoint), midpoint.midpoint(logs[1])] {
                        let trial = evaluate(
                            exp_positive(probe, "phase quarter loses range")?,
                            &mut cancelled,
                        )?;
                        let r = trial.value()?;
                        if r - trial.spread > 0.0 {
                            logs[1] = logs[1].min(probe);
                        } else if r + trial.spread < 0.0 {
                            logs[0] = logs[0].max(probe);
                        }
                    }
                }
                if logs == old {
                    return Err(LigandSpeciationError::Arithmetic(
                        "inner speciation cannot resolve requested phase bracket",
                    ));
                }
            }
            let (answer, iterations) = accepted.ok_or(LigandSpeciationError::Budget)?;
            (answer, PrecipitationStatus::SolidPresent, iterations)
        };
        poll(&mut cancelled)?;
        if answer.state.aqueous.ph_bracket[0].value() < COPPER_SULFIDE_PH_DOMAIN[0]
            || answer.state.aqueous.ph_bracket[1].value() > COPPER_SULFIDE_PH_DOMAIN[1]
        {
            return Err(LigandSpeciationError::Invalid(
                "selected H2S/HS- copper-sulfide model requires final pH in [2,10]",
            ));
        }
        let solid = limit - answer.u;
        if !solid.is_finite()
            || solid < 0.0
            || (status == PrecipitationStatus::SolidPresent && solid == 0.0)
        {
            return Err(LigandSpeciationError::Arithmetic(
                "finite solid partition loses range",
            ));
        }
        let dissolved_copper = add_iter(
            answer.state.metals[self.copper_index]
                .members
                .iter()
                .map(|s| s.amount.value()),
        )?;
        let dissolved_sulfide = add_iter(
            answer.state.aqueous.families[self.sulfide_index]
                .forms
                .iter()
                .map(|s| s.amount.value()),
        )?;
        let copper_residual = component_residual(dissolved_copper, 2.0 * solid, copper_total)?;
        let sulfide_residual = component_residual(dissolved_sulfide, solid, sulfide_total)?;
        Ok(CopperSulfideSpeciation {
            model: self.clone(),
            copper_total: Amount::new(copper_total),
            sulfide_total: Amount::new(sulfide_total),
            dissolved_copper: Amount::new(dissolved_copper),
            dissolved_sulfide: Amount::new(dissolved_sulfide),
            solid_amount: Amount::new(solid),
            copper_balance_residual: Amount::new(copper_residual),
            sulfide_balance_residual: Amount::new(sulfide_residual),
            log_saturation_ratio: answer.residual.map(Dimensionless::new),
            log_saturation_numerical_spread: answer
                .residual
                .map(|_| Dimensionless::new(answer.spread)),
            ligand: answer.state,
            dissolved_limiting_bracket: bracket.map(Amount::new),
            status,
            iterations,
            ligand_solves,
            aqueous_solves,
        })
    }
}

fn unique_total(totals: &[(String, Amount)], id: &str) -> Result<usize> {
    let mut matches = totals
        .iter()
        .enumerate()
        .filter(|(_, (name, _))| name == id);
    let (index, (_, amount)) = matches.next().ok_or(LigandSpeciationError::Invalid(
        "missing precipitating component total",
    ))?;
    if matches.next().is_some() {
        return Err(LigandSpeciationError::Invalid(
            "duplicate precipitating component total",
        ));
    }
    require_nonnegative(amount.value(), "invalid precipitating component total")?;
    Ok(index)
}

fn partition(copper: f64, sulfide: f64) -> Result<(f64, f64, bool)> {
    // Binary division by two is exact unless a positive subnormal loses range.
    let cu_equivalent = copper * 0.5;
    if copper > 0.0 && (cu_equivalent == 0.0 || 2.0 * cu_equivalent != copper) {
        return Err(LigandSpeciationError::Arithmetic(
            "copper formula-equivalent amount loses range",
        ));
    }
    let copper_limits = cu_equivalent <= sulfide;
    let (limit, excess) = if copper_limits {
        (cu_equivalent, sulfide - cu_equivalent)
    } else {
        (sulfide, copper - 2.0 * sulfide)
    };
    if !excess.is_finite() || excess < 0.0 {
        return Err(LigandSpeciationError::Arithmetic(
            "Cu/S stoichiometric excess loses range",
        ));
    }
    Ok((limit, excess, copper_limits))
}
struct Trial {
    u: f64,
    state: LigandSpeciation,
    residual: Option<f64>,
    spread: f64,
}
impl Trial {
    fn value(&self) -> Result<f64> {
        self.residual
            .ok_or(LigandSpeciationError::Arithmetic("empty phase saturation"))
    }
}
fn component_residual(dissolved: f64, solid: f64, total: f64) -> Result<f64> {
    let residual = add([dissolved, solid, -total])?;
    if residual.abs() > 128.0 * f64::EPSILON * total {
        return Err(LigandSpeciationError::Arithmetic(
            "independent Cu/S component balance fails",
        ));
    }
    Ok(residual)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(LigandSpeciationError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
