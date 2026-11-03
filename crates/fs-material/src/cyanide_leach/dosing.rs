//! Finite additions of already-dissolved NaCN and NaOH equivalents.
//!
//! A NaCN equivalent adds one mol to the uncomplexed HCN+CN analytical family
//! and one mol of the explicitly identified Na[+](aq) spectator. A NaOH
//! equivalent adds one mol of Na[+](aq) alkalinity. The original acid/base
//! solver determines free CN, HCN, H+ and OH- simultaneously; neither reagent
//! directly overwrites pH or an equilibrium species amount. This strong-ion
//! input convention is the usual analytical-total/charge-balance description
//! of solution preparation, as in USGS PHREEQC solution definitions:
//! <https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-48.htm>.
//!
//! Amounts are actual delivered dissolved equivalents, not an assertion that a
//! weighed solid dissolved completely. Purity, solubility, dissolution kinetics,
//! carrier volume and mixing heat require separate evidence. The caller declares
//! fixed solution volume, temperature and activity-reference amount as in the
//! original leach owner. No finite solvent H/O balance is claimed.
//!
//! The original initial charge is immutable. Cumulative NaCN, NaOH and sodium
//! inputs remain separate from cumulative gold-reaction progress. Each addition
//! re-speciates the actual solution and any finite competing-metal pools at the
//! retained absolute clock. Bound metal-cyanide ligands remain in the balance;
//! their original equilibrium residual is reported without normalization. No oxygen,
//! gas-transfer or numerical-error history is reset. Source, work-budget,
//! cancellation and arithmetic refusals publish neither a dose nor a state.

use super::*;

/// Explicit neutral-formula-unit delivery into the aqueous solution.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DissolvedSodiumReagent {
    /// One mol Na+ plus one mol of the uncomplexed HCN+CN analytical pool.
    SodiumCyanide,
    /// One mol Na+ alkalinity; actual OH- is owned by simultaneous speciation.
    SodiumHydroxide,
}

/// Cumulative accepted delivered amounts, independent of the original charge.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReagentInputs {
    /// Delivered dissolved NaCN formula units, mol.
    pub sodium_cyanide: Amount,
    /// Delivered dissolved NaOH formula units, mol.
    pub sodium_hydroxide: Amount,
    /// Added Na+, mol; the checked stoichiometric sum of both reagent inputs.
    pub sodium: Amount,
}
impl Default for ReagentInputs {
    fn default() -> Self {
        Self {
            sodium_cyanide: Amount::new(0.0),
            sodium_hydroxide: Amount::new(0.0),
            sodium: Amount::new(0.0),
        }
    }
}

/// Accepted addition and actual post-addition state, retaining original sources.
#[derive(Clone, Debug)]
pub struct ReagentDoseReceipt {
    layout: Arc<Layout>,
    /// Explicit chemical delivery convention.
    pub reagent: DissolvedSodiumReagent,
    /// Delivered dissolved equivalents in this call, mol.
    pub delivered: Amount,
    /// Cumulative accepted inputs after this call.
    pub inputs: ReagentInputs,
    /// Re-speciated solution at unchanged elapsed time and gold progress.
    pub state: CyanidationState,
    /// Uncomplexed + gold/competing-metal ligands - original - delivered NaCN, mol.
    /// Retains the original equilibrium owner's bounded ligand residual.
    pub cyanide_balance: Amount,
    /// Actual Na+ spectator - original Na+ spectator - delivered sodium, mol.
    pub sodium_balance: Amount,
    /// Shared work, including all nested speciation and final cancellation checks.
    pub work: usize,
}
impl ReagentDoseReceipt {
    /// Original empirical kinetics and source-backed acid/base system.
    pub fn model(&self) -> &GoldCyanidation {
        &self.layout.model
    }
    /// Original initial charge, never relabeled as a freshly dosed initial state.
    pub fn initial_charge(&self) -> &CyanidationCharge {
        &self.layout.charge
    }
}

pub(super) const SODIUM: &str = "Na[+](aq)";

pub(super) fn work<F: FnMut() -> bool>(maximum: usize, cancelled: &mut F) -> Result<Work<'_, F>> {
    require(
        (1..=10_000_000).contains(&maximum),
        "reagent addition work must be in 1..=10000000",
    )?;
    Ok(Work {
        used: 0,
        maximum,
        cancelled,
        stopped: None,
    })
}

/// Prepare a candidate layout only. Its caller retains ownership of publication.
/// The returned ceiling excludes oxygen; the closed owner reapplies its finite DO.
pub(super) fn prepare(
    original: &Layout,
    reagent: DissolvedSodiumReagent,
    delivered: Amount,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<Layout> {
    work.tick()?;
    nonnegative(
        delivered.value(),
        "delivered reagent must be finite and nonnegative",
    )?;
    require(
        original.model.chemistry.families().iter().all(|family| {
            family
                .species()
                .iter()
                .all(|species| species.as_str() != SODIUM)
        }),
        "dosed sodium cannot also be owned by an acid/base family",
    )?;
    require(
        original.competition.as_ref().is_none_or(|competition| {
            competition.model.families().iter().all(|family| {
                family
                    .species()
                    .iter()
                    .all(|species| species.as_str() != SODIUM)
            })
        }),
        "dosed sodium cannot also be owned by a competing-metal family",
    )?;
    for ion in &original.charge.spectators {
        if ion.species.as_str() == SODIUM {
            require(
                ion.charge == 1,
                "the original canonical Na[+](aq) spectator must have charge +1",
            )?;
            nonnegative(
                ion.amount.value(),
                "original sodium must be finite and nonnegative",
            )?;
        }
    }
    let mut next = original.clone();
    let input = match reagent {
        DissolvedSodiumReagent::SodiumCyanide => &mut next.reagents.sodium_cyanide,
        DissolvedSodiumReagent::SodiumHydroxide => &mut next.reagents.sodium_hydroxide,
    };
    *input = Amount::new(add(input.value(), delivered.value())?);
    next.reagents.sodium = Amount::new(sum([
        next.reagents.sodium_cyanide.value(),
        next.reagents.sodium_hydroxide.value(),
    ]
    .into_iter())?);
    let cyanide = cyanide_total(&next)?;
    let sodium = sodium_total(&next)?;
    if delivered.value() > 0.0 {
        if next.reagents.sodium.value() <= original.reagents.sodium.value()
            || sodium <= sodium_total(original)?
            || (reagent == DissolvedSodiumReagent::SodiumCyanide
                && cyanide <= cyanide_total(original)?)
        {
            return Err(CyanidationError::Arithmetic(
                "positive reagent dose is below an inventory's resolution",
            ));
        }
    }
    next.limit = next
        .charge
        .accessible_gold
        .value()
        .min(product(0.5, cyanide)?);
    Ok(next)
}

fn add(original: f64, increment: f64) -> Result<f64> {
    let value = sum([original, increment].into_iter())?;
    if increment > 0.0 && value <= original {
        Err(CyanidationError::Arithmetic(
            "positive reagent dose is below cumulative-input resolution",
        ))
    } else {
        Ok(value)
    }
}

pub(super) fn cyanide_total(layout: &Layout) -> Result<f64> {
    sum([
        layout.charge.uncomplexed_cyanide.value(),
        layout.reagents.sodium_cyanide.value(),
    ]
    .into_iter())
}
fn original_sodium(layout: &Layout) -> f64 {
    layout
        .charge
        .spectators
        .iter()
        .find(|ion| ion.species.as_str() == SODIUM)
        .map_or(0.0, |ion| ion.amount.value())
}
fn sodium_total(layout: &Layout) -> Result<f64> {
    sum([original_sodium(layout), layout.reagents.sodium.value()].into_iter())
}
pub(super) fn spectators(layout: &Layout) -> Result<Vec<SpectatorIon>> {
    let mut ions = layout.charge.spectators.clone();
    if layout.reagents.sodium.value() > 0.0 {
        let amount = Amount::new(sodium_total(layout)?);
        if let Some(sodium) = ions.iter_mut().find(|ion| ion.species.as_str() == SODIUM) {
            sodium.amount = amount;
        } else {
            ions.push(SpectatorIon {
                species: SpeciesId::new(SODIUM).expect("fixed sodium identity"),
                charge: 1,
                amount,
            });
        }
    }
    Ok(ions)
}

pub(super) fn balances(layout: &Layout, state: &CyanidationState) -> Result<(f64, f64)> {
    let c = &layout.charge;
    let cyanide = sum([
        state.uncomplexed_cyanide.value(),
        bound_cyanide(state),
        product(2.0, state.gold_complex.value())?,
        -c.uncomplexed_cyanide.value(),
        -product(2.0, c.gold_complex.value())?,
        -layout.reagents.sodium_cyanide.value(),
    ]
    .into_iter())?;
    // Count the fixed counterion pool only. An undosed model may instead own
    // Na through a metal family; such a family is explicitly refused on dosing.
    let tracked_sodium = if layout.reagents.sodium.value() > 0.0
        || c.spectators
            .iter()
            .any(|ion| ion.species.as_str() == SODIUM)
    {
        state
            .speciation
            .spectators
            .iter()
            .find(|ion| ion.species.as_str() == SODIUM)
            .map_or(0.0, |ion| ion.amount.value())
    } else {
        0.0
    };
    let sodium = sum([
        tracked_sodium,
        -original_sodium(layout),
        -layout.reagents.sodium.value(),
    ]
    .into_iter())?;
    for (residual, scale) in [
        (
            cyanide - ligand_balance_residual(state),
            cyanide_total(layout)?.max(product(2.0, c.gold_complex.value())?),
        ),
        (sodium, sodium_total(layout)?),
    ] {
        if residual.abs() > 128.0 * f64::EPSILON * scale {
            return Err(CyanidationError::Arithmetic(
                "dosed cyanide/sodium analytical balance does not close",
            ));
        }
    }
    Ok((cyanide, sodium))
}

pub(super) fn receipt(
    layout: Arc<Layout>,
    reagent: DissolvedSodiumReagent,
    delivered: Amount,
    state: &CyanidationState,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<ReagentDoseReceipt> {
    let (cyanide, sodium) = balances(&layout, state)?;
    work.tick()?;
    Ok(ReagentDoseReceipt {
        inputs: layout.reagents,
        layout,
        reagent,
        delivered,
        state: state.clone(),
        cyanide_balance: Amount::new(cyanide),
        sodium_balance: Amount::new(sodium),
        work: work.used,
    })
}
