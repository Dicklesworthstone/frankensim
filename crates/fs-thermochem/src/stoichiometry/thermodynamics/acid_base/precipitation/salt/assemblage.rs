//! Finite salt competition for one shared cation under the original aqueous law.
//!
//! Each selected phase consumes the common cation and a DISTINCT analytical
//! acid family. A coordinate update holds the other solid stocks fixed and
//! invokes the original `SaltPrecipitation` solve, including its shared pH.
//! Cyclic updates minimize the selected ideal Gibbs model along each phase
//! coordinate. Their stopping condition is joint phase complementarity at ONE
//! final aqueous state, never the collection of stale coordinate residuals.
//!
//! The complete cation and family ledgers are audited independently. Dissolved
//! amounts are carried directly, not reconstructed as total minus all solids.
//! Canonical phase order, bounded aqueous work and cancellation are explicit.
//! Distinct family ownership excludes polymorph/phase-amount degeneracy; this
//! is not a general Gibbs minimizer or arbitrary solid-assemblage engine.
//! Pure solvent activity, ideal aqueous activities and source uncertainty keep
//! the original meanings. No solvent H/O ledger, kinetics or gas loss is added.

use super::{
    AcidBaseError, AcidBaseSpeciation, ConservedReaction, IdealAcidBaseSystem,
    PrecipitationOptions, PrecipitationStatus, ReactionNetwork, Result, SaltPrecipitation,
    SpeciesId, SpectatorIon, Sum, consistent_sources, nonnegative, poll, positive,
};
use fs_math::det;
use fs_qty::{Amount, Dimensionless};
use std::collections::BTreeSet;

/// Bounded selected assemblage; each phase owns a distinct original acid family.
pub const MAX_SALT_PHASES: usize = 8;

/// Canonically ordered source-bound phases sharing one aqueous model and cation.
#[derive(Clone, Debug)]
pub struct CommonCationSaltSystem {
    phases: Vec<SaltPrecipitation>,
    network: ReactionNetwork,
}

/// Joint accuracy and total work limits, in addition to the original phase limits.
#[derive(Clone, Copy, Debug)]
pub struct SaltAssemblageOptions {
    /// Original coordinate and inner aqueous controls, unchanged in meaning.
    pub phase: PrecipitationOptions,
    /// Positive absolute allowance for each complete analytical ledger [mol].
    pub absolute_balance_tolerance: Amount,
    /// Nonnegative relative allowance against each original component total.
    pub relative_balance_tolerance: Dimensionless,
    /// At most 256 complete canonical passes through the selected phases.
    pub max_sweeps: usize,
    /// Hard total bound, including initial and all nested aqueous solves.
    /// At most 65536; each uses the original bounded charge-balance solver.
    pub max_aqueous_solves: usize,
}
impl Default for SaltAssemblageOptions {
    fn default() -> Self {
        Self {
            phase: PrecipitationOptions::default(),
            absolute_balance_tolerance: Amount::new(1e-14),
            relative_balance_tolerance: Dimensionless::new(1e-11),
            max_sweeps: 128,
            max_aqueous_solves: 16_384,
        }
    }
}

/// One phase evaluated at the final COMMON solution, not its last local update.
#[derive(Clone, Debug)]
pub struct SaltPhaseEquilibrium {
    /// Original phase-specific identity and formula-unit basis.
    pub solid: SpeciesId,
    /// Original acid family whose terminal ion enters this solid.
    pub family: String,
    /// Nonnegative pure-solid formula units [mol].
    pub amount: Amount,
    /// Present, absent, or an empty original analytical component.
    pub status: PrecipitationStatus,
    /// Natural log dissolution IAP/K, undefined for an empty component.
    pub log_saturation_ratio: Option<Dimensionless>,
    /// Inner pH resolution contribution, not an interval certificate.
    pub log_saturation_numerical_spread: Option<Dimensionless>,
}

/// Independently reduced dissolved plus solid analytical-family account.
#[derive(Clone, Debug)]
pub struct SaltFamilyBalance {
    /// Original complete aqueous family identity.
    pub family: String,
    /// Input total, including its component in any supplied solids [mol].
    pub total: Amount,
    /// Independently summed original aqueous member amounts [mol].
    pub dissolved: Amount,
    /// Family-component units held in the selected solid, if any [mol].
    pub in_solid: Amount,
    /// Dissolved plus solid minus original total [mol].
    pub residual: Amount,
}

/// Complete admitted finite assemblage with its original source owners and limits.
#[derive(Clone, Debug)]
pub struct SaltAssemblageSpeciation {
    model: CommonCationSaltSystem,
    /// Original cation-component amount in all dissolved and solid forms [mol].
    pub cation_total: Amount,
    /// Free cation stock, carried directly to retain positive traces [mol].
    pub dissolved_cation: Amount,
    /// Independent free cation plus all solid contributions minus total [mol].
    pub cation_balance_residual: Amount,
    /// Every original analytical family, including those without a selected phase.
    pub families: Vec<SaltFamilyBalance>,
    /// Final complementarity diagnostics in canonical solid-identity order.
    pub phases: Vec<SaltPhaseEquilibrium>,
    /// One final shared pH, original family speciation and fixed spectator stocks.
    pub aqueous: AcidBaseSpeciation,
    /// Original numerical controls; they are not source uncertainty estimates.
    pub options: SaltAssemblageOptions,
    /// Completed canonical phase sweeps; zero for initially undersaturated inputs.
    pub sweeps: usize,
    /// Original finite single-salt calls completed.
    pub coordinate_solves: usize,
    /// Total original aqueous solves, including initial and discarded trials.
    pub aqueous_solves: usize,
    /// Total inner pH bisections across those aqueous solves.
    pub aqueous_iterations: usize,
}
impl SaltAssemblageSpeciation {
    /// Original phase, aqueous and formation-source owners.
    pub fn model(&self) -> &CommonCationSaltSystem {
        &self.model
    }
}

impl CommonCationSaltSystem {
    /// Join 1..=8 original salts with the SAME cation and complete aqueous model.
    /// Each phase must own a distinct acid family and distinct solid identity.
    /// Shared equations, species definitions and selected Gibbs receipts must
    /// agree exactly. Polymorph competition and overlapping phase families refuse.
    pub fn new(mut phases: Vec<SaltPrecipitation>) -> Result<Self> {
        if phases.is_empty() || phases.len() > MAX_SALT_PHASES {
            return Err(AcidBaseError::Invalid("require 1..=8 selected salt phases"));
        }
        phases.sort_by(|a, b| a.solid.cmp(&b.solid));
        let first = &phases[0];
        let mut families = BTreeSet::new();
        let mut solids = BTreeSet::new();
        let mut standards = Vec::new();
        let mut reactions = first.aqueous.reaction_network().reactions().to_vec();
        for phase in &phases {
            if phase.cation != first.cation
                || phase.charge != first.charge
                || !same_aqueous(&first.aqueous, &phase.aqueous)
                || !families.insert(phase.family_id())
                || !solids.insert(&phase.solid)
            {
                return Err(AcidBaseError::Invalid(
                    "salts require identical aqueous owners and cation, with disjoint family/phase identities",
                ));
            }
            standards.push(&phase.standard);
            standards.push(phase.aqueous.water().standard());
            for family in phase.aqueous.families() {
                standards.extend(family.steps().iter().map(|step| step.standard()));
            }
            reactions.push(phase.reaction.clone());
        }
        consistent_sources(standards)?;
        let network = ReactionNetwork::new(reactions)?;
        Ok(Self { phases, network })
    }

    /// Original phases in canonical solid-identity order.
    pub fn phases(&self) -> &[SaltPrecipitation] {
        &self.phases
    }
    /// Complete shared original aqueous owner.
    pub fn aqueous_model(&self) -> &IdealAcidBaseSystem {
        &self.phases[0].aqueous
    }
    /// Shared source-defined free-cation identity.
    pub fn cation(&self) -> &SpeciesId {
        &self.phases[0].cation
    }
    /// All original exact equations, with common species definitions admitted.
    pub fn reaction_network(&self) -> &ReactionNetwork {
        &self.network
    }

    /// Equilibrate finite totals with all selected phases and one common pH.
    ///
    /// Inputs are analytical totals INCLUDING their solid contents, supplied once.
    /// Initial solid partition is immaterial to this selected equilibrium model;
    /// it supplies no rate or metastable persistence. Each coordinate can dissolve
    /// only its currently assigned solid; other solid stocks remain unavailable.
    /// All complete aqueous family totals and spectator ions are explicit.
    /// Cancellation, range loss, unresolved complementarity or exhausted work
    /// returns no partial state and never mutates the original inputs.
    pub fn speciate(
        &self,
        cation_total: Amount,
        family_totals: &[(String, Amount)],
        spectators: &[SpectatorIon],
        amount_at_unit_activity: Amount,
        options: SaltAssemblageOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SaltAssemblageSpeciation> {
        poll(&mut cancelled)?;
        admit_options(options)?;
        nonnegative(cation_total.value(), "invalid assemblage cation total")?;
        let basis = positive(
            amount_at_unit_activity.value(),
            "invalid common aqueous basis",
        )?;
        if spectators.iter().any(|ion| {
            ion.species == *self.cation()
                || self.phases.iter().any(|phase| phase.solid == ion.species)
        }) {
            return Err(AcidBaseError::Invalid(
                "spectator duplicates an owned cation or phase",
            ));
        }
        let mut ions = spectators.to_vec();
        ions.push(SpectatorIon {
            species: self.cation().clone(),
            charge: self.phases[0].charge,
            amount: cation_total,
        });
        // Original owner admits/canonicalizes all identities, totals and controls.
        let mut aqueous = self.aqueous_model().speciate(
            family_totals,
            &ions,
            amount_at_unit_activity,
            options.phase.aqueous,
            &mut cancelled,
        )?;
        let totals: Vec<_> = aqueous
            .families
            .iter()
            .map(|f| (f.id.clone(), f.total))
            .collect();
        let mut dissolved = totals.clone();
        let mut metal = cation_total.value();
        let mut solids = vec![0.0; self.phases.len()];
        let mut aqueous_solves = 1;
        let mut aqueous_iterations = aqueous.iterations;
        let mut coordinate_solves = 0;
        for sweep in 0..=options.max_sweeps {
            poll(&mut cancelled)?;
            let (phase_reports, converged) = self.phase_reports(
                &aqueous,
                metal,
                &solids,
                &totals,
                cation_total,
                basis,
                options.phase.log_saturation_tolerance.value(),
            )?;
            let (cation_residual, family_balances) =
                self.balances(&aqueous, metal, &solids, &totals, cation_total, options)?;
            if converged {
                return Ok(SaltAssemblageSpeciation {
                    model: self.clone(),
                    cation_total,
                    dissolved_cation: Amount::new(metal),
                    cation_balance_residual: Amount::new(cation_residual),
                    families: family_balances,
                    phases: phase_reports,
                    aqueous,
                    options,
                    sweeps: sweep,
                    coordinate_solves,
                    aqueous_solves,
                    aqueous_iterations,
                });
            }
            if sweep == options.max_sweeps {
                return Err(AcidBaseError::Budget);
            }
            let old_metal = metal;
            let old_solids = solids.clone();
            for (i, phase) in self.phases.iter().enumerate() {
                poll(&mut cancelled)?;
                let remaining = options.max_aqueous_solves - aqueous_solves;
                // Each original phase call evaluates at most two endpoints plus
                // max_iterations trials. Reserve that bound BEFORE starting it.
                let rounds = remaining
                    .saturating_sub(2)
                    .min(options.phase.max_iterations);
                if rounds == 0 {
                    return Err(AcidBaseError::Budget);
                }
                let local_metal = sum([metal, f64::from(phase.cation_stoichiometry) * solids[i]])?;
                let local_family = sum([
                    dissolved[phase.family_index].1.value(),
                    f64::from(phase.anion_stoichiometry) * solids[i],
                ])?;
                let mut local_totals = dissolved.clone();
                local_totals[phase.family_index].1 = Amount::new(local_family);
                let answer = phase.speciate(
                    Amount::new(local_metal),
                    &local_totals,
                    spectators,
                    amount_at_unit_activity,
                    PrecipitationOptions {
                        max_iterations: rounds,
                        ..options.phase
                    },
                    &mut cancelled,
                )?;
                aqueous_solves += answer.aqueous_solves;
                aqueous_iterations += answer.aqueous_iterations;
                coordinate_solves += 1;
                metal = answer.dissolved_cation.value();
                dissolved[phase.family_index].1 = answer.dissolved_family;
                solids[i] = answer.solid_amount.value();
                aqueous = answer.aqueous;
            }
            if metal == old_metal && solids == old_solids {
                return Err(AcidBaseError::Arithmetic(
                    "coordinate updates cannot resolve joint salt complementarity",
                ));
            }
        }
        Err(AcidBaseError::Budget)
    }

    #[allow(clippy::too_many_arguments)]
    fn phase_reports(
        &self,
        aqueous: &AcidBaseSpeciation,
        metal: f64,
        solids: &[f64],
        totals: &[(String, Amount)],
        cation_total: Amount,
        basis: f64,
        tolerance: f64,
    ) -> Result<(Vec<SaltPhaseEquilibrium>, bool)> {
        let mut reports = Vec::with_capacity(self.phases.len());
        let mut all_satisfied = true;
        for (i, phase) in self.phases.iter().enumerate() {
            let empty = cation_total.value() == 0.0 || totals[phase.family_index].1.value() == 0.0;
            let (ratio, spread) = if empty {
                if solids[i] != 0.0 {
                    return Err(AcidBaseError::Arithmetic(
                        "solid without its finite component",
                    ));
                }
                (None, None)
            } else {
                let family = &aqueous.families[phase.family_index];
                let anion = family
                    .forms
                    .last()
                    .ok_or(AcidBaseError::Invalid("empty original family"))?;
                positive(metal, "positive shared cation trace loses range")?;
                let activity = positive(
                    anion.activity.value(),
                    "positive salt anion trace loses range",
                )?;
                let p = f64::from(phase.cation_stoichiometry);
                let q = f64::from(phase.anion_stoichiometry);
                let r = sum([
                    p * (det::ln(metal) - det::ln(basis)),
                    q * det::ln(activity),
                    phase.standard.log_equilibrium_constant().value(),
                ])?;
                let spread = q
                    * phase.aqueous.families()[phase.family_index].steps().len() as f64
                    * core::f64::consts::LN_10
                    * aqueous
                        .ph_bracket
                        .iter()
                        .map(|v| (v.value() - aqueous.ph.value()).abs())
                        .fold(0.0, f64::max);
                nonnegative(spread, "invalid final saturation resolution")?;
                all_satisfied &= if solids[i] > 0.0 {
                    r.abs() + spread <= tolerance
                } else {
                    r + spread <= 0.0
                };
                (
                    Some(Dimensionless::new(r)),
                    Some(Dimensionless::new(spread)),
                )
            };
            reports.push(SaltPhaseEquilibrium {
                solid: phase.solid.clone(),
                family: phase.family_id().into(),
                amount: Amount::new(solids[i]),
                status: if empty {
                    PrecipitationStatus::EmptyComponent
                } else if solids[i] > 0.0 {
                    PrecipitationStatus::SolidPresent
                } else {
                    PrecipitationStatus::SolidAbsent
                },
                log_saturation_ratio: ratio,
                log_saturation_numerical_spread: spread,
            });
        }
        Ok((reports, all_satisfied))
    }

    fn balances(
        &self,
        aqueous: &AcidBaseSpeciation,
        metal: f64,
        solids: &[f64],
        totals: &[(String, Amount)],
        cation_total: Amount,
        options: SaltAssemblageOptions,
    ) -> Result<(f64, Vec<SaltFamilyBalance>)> {
        let mut cation = Sum::default();
        cation.add(metal)?;
        cation.add(-cation_total.value())?;
        for (phase, &solid) in self.phases.iter().zip(solids) {
            cation.add(f64::from(phase.cation_stoichiometry) * solid)?;
        }
        let cation_residual = cation.finish()?;
        admit_balance(cation_residual, cation_total.value(), options)?;
        let mut balances = Vec::with_capacity(totals.len());
        for (i, (id, total)) in totals.iter().enumerate() {
            let dissolved = sum(aqueous.families[i]
                .forms
                .iter()
                .map(|form| form.amount.value()))?;
            let solid = sum(self
                .phases
                .iter()
                .zip(solids)
                .filter(|(phase, _)| phase.family_index == i)
                .map(|(phase, &amount)| f64::from(phase.anion_stoichiometry) * amount))?;
            let residual = sum([dissolved, solid, -total.value()])?;
            admit_balance(residual, total.value(), options)?;
            balances.push(SaltFamilyBalance {
                family: id.clone(),
                total: *total,
                dissolved: Amount::new(dissolved),
                in_solid: Amount::new(solid),
                residual: Amount::new(residual),
            });
        }
        Ok((cation_residual, balances))
    }
}

fn same_reaction(a: &ConservedReaction, b: &ConservedReaction) -> bool {
    a.id() == b.id()
        && a.terms() == b.terms()
        && a.terms().iter().all(|term| {
            a.species_charge(&term.species) == b.species_charge(&term.species)
                && a.species_elements(&term.species) == b.species_elements(&term.species)
        })
}
fn same_aqueous(a: &IdealAcidBaseSystem, b: &IdealAcidBaseSystem) -> bool {
    same_reaction(a.water().reaction(), b.water().reaction())
        && a.families().len() == b.families().len()
        && a.families().iter().zip(b.families()).all(|(x, y)| {
            x.id() == y.id()
                && x.species() == y.species()
                && x.steps().len() == y.steps().len()
                && x.steps()
                    .iter()
                    .zip(y.steps())
                    .all(|(u, v)| same_reaction(u.reaction(), v.reaction()))
        })
}
fn admit_options(options: SaltAssemblageOptions) -> Result<()> {
    positive(
        options.absolute_balance_tolerance.value(),
        "invalid absolute component tolerance",
    )?;
    nonnegative(
        options.relative_balance_tolerance.value(),
        "invalid relative component tolerance",
    )?;
    positive(
        options.phase.minimum_dissolved_amount.value(),
        "invalid phase trace floor",
    )?;
    positive(
        options.phase.log_activity_tolerance.value(),
        "invalid phase bracket tolerance",
    )?;
    positive(
        options.phase.log_saturation_tolerance.value(),
        "invalid phase saturation tolerance",
    )?;
    if !(1..=256).contains(&options.max_sweeps)
        || !(1..=65_536).contains(&options.max_aqueous_solves)
        || !(1..=256).contains(&options.phase.max_iterations)
    {
        return Err(AcidBaseError::Invalid(
            "invalid finite assemblage work limits",
        ));
    }
    Ok(())
}
fn admit_balance(residual: f64, total: f64, options: SaltAssemblageOptions) -> Result<()> {
    let allowance = options.absolute_balance_tolerance.value()
        + options.relative_balance_tolerance.value() * total;
    if !allowance.is_finite() || residual.abs() > allowance {
        Err(AcidBaseError::Arithmetic(
            "complete salt assemblage component ledger does not close",
        ))
    } else {
        Ok(())
    }
}
fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64> {
    let mut sum = Sum::default();
    for value in values {
        sum.add(value)?;
    }
    sum.finish().map_err(Into::into)
}

#[cfg(test)]
mod tests;
