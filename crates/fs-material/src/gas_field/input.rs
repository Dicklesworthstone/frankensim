//! Explicit species columns shared by prescribed and generated gas-field inputs.
//!
//! This adapter binds the original thermochemical pack owner and reorders only
//! authored primitive coordinates. It never averages derived gas properties or
//! inserts a physical trace concentration for an absent species.

use super::PrimitiveGasState;
use crate::gas_species::input::{GasInputError, ThermochemicalPhase, ThermochemicalSources};
use crate::thermochemical_gas::ThermochemicalGasModel;
use crate::thermochemical_gas::packs::{SpeciesPackRef, from_packs};
use fs_qty::semantic::{Composition, CompositionBasis};
use std::str::SplitAsciiWhitespace;

/// Bound retained per-state component receipts before reading the source packs.
pub const MAX_FIELD_STATE_COMPONENTS: usize = 16_384;
const MAX_FIELD_SPECIES: usize = 128;

/// One original source catalogue and one authored-to-canonical column mapping.
pub struct PrimitiveGasColumns {
    order: Vec<usize>,
    model: ThermochemicalGasModel,
    phase: ThermochemicalPhase,
}

fn word<'a>(fields: &mut SplitAsciiWhitespace<'a>) -> Result<&'a str, GasInputError> {
    fields.next().ok_or(GasInputError::Input("missing gas-field source value"))
}
fn scalar(fields: &mut SplitAsciiWhitespace<'_>) -> Result<f64, GasInputError> {
    let value: f64 = word(fields)?.parse()
        .map_err(|_| GasInputError::Input("invalid gas-field scalar"))?;
    if !value.is_finite() {
        return Err(GasInputError::Input("gas-field scalars must be finite"));
    }
    Ok(value)
}
fn material(what: String) -> GasInputError {
    GasInputError::Material(crate::MaterialError::Parameters { what })
}

impl PrimitiveGasColumns {
    /// Consume mole COUNT ID..., retaining suffix fields for the caller.
    /// The states argument is the complete station/inlet/injection count.
    /// Input is consumed only after every source and budget check succeeds.
    pub fn read(
        fields: &mut SplitAsciiWhitespace<'_>,
        sources: Option<&ThermochemicalSources<'_>>,
        states: usize,
    ) -> Result<Self, GasInputError> {
        let mut remaining = fields.clone();
        if word(&mut remaining)? != "mole" {
            return Err(GasInputError::Input(
                "gas-field primitive compositions require explicit mole fractions",
            ));
        }
        let count: usize = word(&mut remaining)?.parse()
            .map_err(|_| GasInputError::Input("invalid gas-field species count"))?;
        if count == 0 || count > MAX_FIELD_SPECIES {
            return Err(GasInputError::Input("gas-field species count is outside 1..=128"));
        }
        if states == 0 || states.checked_mul(count)
            .is_none_or(|terms| terms > MAX_FIELD_STATE_COMPONENTS)
        {
            return Err(GasInputError::Input(
                "gas-field state/species product exceeds the retained source budget",
            ));
        }
        let sources = sources.ok_or(GasInputError::Input(
            "gas field requires an explicit thermochemical source registry",
        ))?;
        let mut selected: Vec<(usize, SpeciesPackRef<'_>)> = Vec::with_capacity(count);
        for column in 0..count {
            let name = word(&mut remaining)?;
            let mut candidates = sources.species.iter()
                .filter(|source| source.species.association().species().as_str() == name);
            let source = candidates.next().ok_or(GasInputError::Input(
                "gas-field species has no source registry entry",
            ))?;
            if candidates.next().is_some() {
                return Err(GasInputError::Input("ambiguous gas-field source registry entry"));
            }
            if selected.iter().any(|(_, source)|
                source.species.association().species().as_str() == name)
            {
                return Err(GasInputError::Input("duplicate gas-field species column"));
            }
            selected.push((column, (*source).clone()));
        }
        selected.sort_by(|(_, a), (_, b)| {
            a.species.association().species().cmp(b.species.association().species())
        });
        // A strictly positive dyadic seed retains the complete source catalogue.
        // It is never evaluated as a physical mixture. Every actual evaluation
        // uses the original station's fractions, including exact zeros.
        let weight = 1.0 / count.next_power_of_two() as f64;
        let mut seed = vec![weight; count];
        seed[0] = 1.0 - weight * (count - 1) as f64;
        let composition = Composition::new(CompositionBasis::MoleFraction, seed)
            .map_err(|error| material(format!("gas-field catalogue composition: {error}")))?;
        let order = selected.iter().map(|(column, _)| *column).collect();
        let selected: Vec<_> = selected.into_iter().map(|(_, source)| source).collect();
        let water = sources.water_species.filter(|id| {
            selected.iter().any(|source| source.species.association().species() == *id)
        });
        let model = from_packs(&selected, composition, water).map_err(GasInputError::Material)?;
        *fields = remaining;
        Ok(Self { order, model, phase: sources.phase })
    }

    /// Consume T_K FRACTION... in the declared column order.
    /// Positive temperature and finite scalars are checked here; composition
    /// admission and every subsequent EOS/transport evaluation retain their
    /// original physical owners. Input is unchanged on refusal.
    pub fn read_state(
        &self, fields: &mut SplitAsciiWhitespace<'_>,
    ) -> Result<PrimitiveGasState, GasInputError> {
        let mut remaining = fields.clone();
        let temperature_k = scalar(&mut remaining)?;
        if temperature_k <= 0.0 {
            return Err(GasInputError::Input("gas-field temperature must be positive"));
        }
        let supplied = (0..self.order.len())
            .map(|_| scalar(&mut remaining))
            .collect::<Result<Vec<_>, _>>()?;
        let fractions = self.order.iter().map(|column| supplied[*column]).collect();
        let composition = Composition::new(CompositionBasis::MoleFraction, fractions)
            .map_err(|error| material(format!("gas-field primitive composition: {error}")))?;
        *fields = remaining;
        Ok(PrimitiveGasState { temperature_k, composition })
    }

    /// Original bound model and explicit phase policy, without reevaluation.
    #[must_use]
    pub fn into_parts(self) -> (ThermochemicalGasModel, ThermochemicalPhase) {
        (self.model, self.phase)
    }
}
