//! Atomic dissolved-reagent additions to an operating carbon-in-leach bath.
//!
//! The existing delivery owner supplies NaCN/NaOH equivalents. This module
//! retains the complete CIL inventory, including cyanide and counterions on
//! every carbon population. A dose changes no elapsed time, leach progress,
//! sorbent loading, oxygen stock, or previously accumulated numerical error.
//! Fixed solution volume and temperature retain the original delivery model.

use super::*;
use crate::cyanide_leach::dosing::DissolvedSodiumReagent;

/// Descriptive alias for the existing public CIL dose receipt.
pub type CilReagentDoseReceipt = CilDoseReceipt;

impl CarbonInLeachBatch {
    /// Add delivered dissolved NaCN/NaOH and re-admit the complete current bath.
    ///
    /// Carbon inventories and bound ligands remain present. The finite initial
    /// dissolved-oxygen ceiling stays in force. A failure publishes no dose,
    /// state, or history. No target pH or instantaneous adsorption is imposed.
    pub fn dose_reagent(
        &mut self,
        reagent: DissolvedSodiumReagent,
        delivered: Amount,
        maximum_work: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CilReagentDoseReceipt> {
        let mut work = dosing::work(maximum_work, &mut cancelled)?;
        let next = prepare_addition(self, reagent, delivered, false, &mut work)?;
        let conservation = conservation_cil(&next.layout, &next.state)?;
        let receipt = receipt(&next, reagent, delivered, conservation, &mut work)?;
        *self = next;
        Ok(receipt)
    }
}

pub(super) fn prepare_addition(
    original: &CarbonInLeachBatch,
    reagent: DissolvedSodiumReagent,
    delivered: Amount,
    open_oxygen: bool,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<CarbonInLeachBatch> {
    let mut leach = dosing::prepare(&original.layout.leach, reagent, delivered, work)?;
    if !open_oxygen {
        leach.limit = leach
            .limit
            .min(product(4.0, leach.charge.dissolved_oxygen.value())?);
    }
    let mut layout = (*original.layout).clone();
    let mut stocks = Stocks::from_state(&original.state);
    if reagent == DissolvedSodiumReagent::SodiumCyanide {
        layout.total_cyanide = added(layout.total_cyanide, delivered.value())?;
    }
    // A potassium (or other admitted) pair is not converted to a sodium pair
    // merely because the delivered reagent contains sodium.
    if layout.pair.counterion().as_str() == dosing::SODIUM {
        stocks.counterion = added(stocks.counterion, delivered.value())?;
        layout.total_counterion = added(layout.total_counterion, delivered.value())?;
    }
    let remaining = positive_or_zero(
        leach.limit - stocks.gold,
        "dosed CIL progress exceeds its admitted leach ceiling",
    )?;
    ensure(
        stocks.counterion >= sum([stocks.dissolved, remaining].into_iter())?,
        "reagent addition would exceed the ion-pair counterion reserve",
    )?;
    layout.leach = Arc::new(leach);
    let layout = Arc::new(layout);
    let state = evaluate_cil(&layout, original.state.elapsed.value(), &stocks, work)?;
    Ok(CarbonInLeachBatch {
        layout,
        state,
        accumulated_error: original.accumulated_error,
    })
}

pub(super) fn receipt(
    next: &CarbonInLeachBatch,
    reagent: DissolvedSodiumReagent,
    delivered: Amount,
    conservation: CilConservation,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<CilReagentDoseReceipt> {
    let sodium = if next.layout.pair.counterion().as_str() == dosing::SODIUM {
        conservation.counterion.value()
    } else {
        let actual = next
            .state
            .solution
            .speciation
            .spectators
            .iter()
            .find(|ion| ion.species.as_str() == dosing::SODIUM)
            .map_or(0.0, |ion| ion.amount.value());
        let initial = next
            .layout
            .leach
            .charge
            .spectators
            .iter()
            .find(|ion| ion.species.as_str() == dosing::SODIUM)
            .map_or(0.0, |ion| ion.amount.value());
        let expected = sum([initial, next.layout.leach.reagents.sodium.value()].into_iter())?;
        let residual = sum([actual, -expected].into_iter())?;
        check_balance(residual, expected)?;
        residual
    };
    work.tick()?;
    Ok(CilReagentDoseReceipt {
        layout: next.layout.clone(),
        reagent,
        delivered,
        inputs: next.layout.leach.reagents,
        state: next.state.clone(),
        conservation,
        sodium_balance: Amount::new(sodium),
        work: work.used,
    })
}

fn added(original: f64, delivered: f64) -> Result<f64> {
    let next = sum([original, delivered].into_iter())?;
    if delivered > 0.0 && next <= original {
        Err(CilError::Arithmetic(
            "reagent addition is below CIL inventory resolution",
        ))
    } else {
        Ok(next)
    }
}

#[cfg(test)]
mod tests;
