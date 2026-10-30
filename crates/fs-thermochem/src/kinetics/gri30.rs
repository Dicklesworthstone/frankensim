//! Complete 53-species, 325-channel GRI-Mech 3.0 from Cantera 3.2.0.
//!
//! Source: <https://github.com/Cantera/cantera/blob/4a8358eb80cfeb50474386b5f9ec0b3a83519889/data/gri30.yaml>.
//! Original GRI-Mech version 3.0 is dated 1999-07-30. All original duplicate
//! channels, twelve generic third-body laws, nine explicit collision partners,
//! three Lindemann and 26 Troe falloff channels are retained. The 16 originally
//! irreversible channels stay irreversible: no reverse rate is manufactured.
//! Consequently this complete source does not have a universal detailed-balance
//! or nonnegative chemical entropy-production guarantee. Consumers retain their
//! finite-step physical-domain and entropy acceptance.
//!
//! Original NASA7 coefficients embed in the existing NASA9 owner. Each original
//! species interval and seam is retained; the common interval is 300..3000 K,
//! limited above by CH3O. The optional explicit continuity treatment adjusts
//! only constant high-region enthalpy and entropy offsets, never Cp or the
//! source cards. Original source h/s is also available.
//!
//! Source A values use cm, mol, s; activation energies use cal/mol. Construction
//! converts them once to the existing reference-rate convention. All subsequent
//! thermodynamics, reaction progress, derivatives, colliders and falloff use the
//! original reusable owners. No transport law or fuel-specific experimental
//! validation follows from loading this methane/natural-gas mechanism. In
//! particular, its small hydrocarbon chemistry is not a soot model or a wood or
//! PMMA volatile mechanism. Original attribution and redistribution notice are
//! retained in `gri30/LICENSE-Cantera.txt`.

use super::{
    ArrheniusRateV1, ElementaryReactionV1, FalloffBroadeningV1, IdealGasKineticsV1,
    KineticsErrorV1, PressureDependenceV1, ReactionProgressRateQuantityV1, ThirdBodyEfficienciesV1,
    TroeFalloffV1,
};
use crate::equilibrium::composition::{CompositionEquilibriumErrorV1, IdealGasEquilibriumV1};
use crate::{
    ChargeVector, ChemistryError, ElementId, ElementalMatrix, ElementalReferenceIdV1,
    MolarEnergyQuantityV1, MolarThermalQuantityV1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReactionId, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    required_nasa9_parameter_dims,
};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{MolarMass, Pressure, Temperature};
use std::collections::BTreeMap;

mod data;
pub use data::{CANONICAL_SPECIES, MOLAR_MASSES_KG_PER_MOL};

/// Pinned Cantera release commit, not a mutable source branch.
pub const SOURCE_COMMIT: &str = "4a8358eb80cfeb50474386b5f9ec0b3a83519889";
/// Exact original `data/gri30.yaml` Git blob.
pub const SOURCE_BLOB: &str = "9a62004ec17bbec8ac11e97785932d426dd58610";
/// Original `src/thermo/Elements.cpp` Git blob supplying atomic weights.
pub const ATOMIC_WEIGHTS_BLOB: &str = "ebe9da8b1a3e7944d1de070851b36fdb858ad8d4";
/// Common lower source bound in kelvin.
pub const MINIMUM_TEMPERATURE_K: f64 = 300.0;
/// Common upper source bound in kelvin, set by the original CH3O card.
pub const MAXIMUM_TEMPERATURE_K: f64 = 3000.0;

/// Explicit choice of original versus continuous piecewise caloric references.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Gri30Calorics {
    /// Original NASA7 h/s, including independently fitted seam discontinuities.
    Original,
    /// Retain the entire low-temperature region; reconcile upper-region h/s
    /// through the existing continuity owner. Bounds are 0.2 J/mol and
    /// 2e-4 J/(mol K), numerical source adjustments, not empirical uncertainty.
    ContinuousAt300K,
}

/// Original owner refusing construction; no partial mechanism is returned.
#[derive(Clone, Debug, PartialEq)]
pub enum Gri30MechanismError {
    /// Species or exact elemental/charge assembly refused.
    Chemistry(ChemistryError),
    /// Original standard-state or explicit continuity owner refused.
    Thermodynamics(ThermochemErrorV1),
    /// Original neutral species bank refused.
    Equilibrium(CompositionEquilibriumErrorV1),
    /// Original elementary kinetic owner refused.
    Kinetics(KineticsErrorV1),
}
impl core::fmt::Display for Gri30MechanismError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Cantera 3.2.0 GRI-Mech 3.0 source: {self:?}")
    }
}
impl std::error::Error for Gri30MechanismError {}
impl From<ChemistryError> for Gri30MechanismError {
    fn from(error: ChemistryError) -> Self {
        Self::Chemistry(error)
    }
}
impl From<ThermochemErrorV1> for Gri30MechanismError {
    fn from(error: ThermochemErrorV1) -> Self {
        Self::Thermodynamics(error)
    }
}
impl From<CompositionEquilibriumErrorV1> for Gri30MechanismError {
    fn from(error: CompositionEquilibriumErrorV1) -> Self {
        Self::Equilibrium(error)
    }
}
impl From<KineticsErrorV1> for Gri30MechanismError {
    fn from(error: KineticsErrorV1) -> Self {
        Self::Kinetics(error)
    }
}

/// Load all source channels once into the existing immutable kinetic mechanism.
///
/// Construction is bounded by 53 species, 106 NASA regions and 325 reactions.
/// Prepared rate evaluations allocate no per-query buffers. Source molecular
/// mass balance is admitted at 1e-12 relative without adjusting any mass.
pub fn cantera_gri30(calorics: Gri30Calorics) -> Result<IdealGasKineticsV1, Gri30MechanismError> {
    let models = (0..CANONICAL_SPECIES.len())
        .map(|i| species_model(i, calorics))
        .collect::<Result<Vec<_>, _>>()?;
    let ids = models
        .iter()
        .map(|model| model.species().clone())
        .collect::<Vec<_>>();
    let elements = ["Ar", "C", "H", "N", "O"]
        .iter()
        .map(|symbol| ElementId::new(*symbol))
        .collect::<Result<Vec<_>, _>>()?;
    let atoms = (0..elements.len())
        .map(|element| {
            data::SPECIES
                .iter()
                .map(|species| species.atoms[element])
                .collect()
        })
        .collect();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(elements, ids.clone(), atoms)?,
        ChargeVector::new(ids, vec![0; CANONICAL_SPECIES.len()])?,
        models,
    )?;
    let mut reactions = Vec::with_capacity(data::REACTIONS.len());
    for (i, source) in data::REACTIONS.iter().enumerate() {
        let order = source
            .reactants
            .iter()
            .map(|(_, exponent)| *exponent)
            .sum::<u32>();
        let generic_collider = matches!(source.pressure, SourcePressure::ThirdBody(_));
        let forward = arrhenius(source.rate, order + u32::from(generic_collider))?;
        let pressure = match source.pressure {
            SourcePressure::Independent => PressureDependenceV1::Independent,
            SourcePressure::ThirdBody(efficiencies) => {
                PressureDependenceV1::ThirdBody(colliders(efficiencies)?)
            }
            SourcePressure::Falloff {
                low,
                efficiencies,
                troe,
            } => PressureDependenceV1::Falloff {
                low_pressure_rate: arrhenius(low, order + 1)?,
                colliders: colliders(efficiencies)?,
                broadening: match troe {
                    None => FalloffBroadeningV1::Lindemann,
                    Some([alpha, t3, t1, t2]) => FalloffBroadeningV1::Troe(TroeFalloffV1::new(
                        alpha,
                        Temperature::new(t3),
                        Temperature::new(t1),
                        Some(Temperature::new(t2)),
                    )?),
                },
            },
        };
        let reaction = ElementaryReactionV1::new(
            ReactionId::new(format!("cantera-gri30-r{:03}", i + 1))?,
            side(source.reactants)?,
            side(source.products)?,
            forward,
        )?
        .with_pressure_dependence(pressure);
        reactions.push(if source.reversible {
            reaction
        } else {
            reaction.irreversible()
        });
    }
    let identity = match calorics {
        Gri30Calorics::Original => "cantera-v3.2.0-gri30-original",
        Gri30Calorics::ContinuousAt300K => "cantera-v3.2.0-gri30-continuous-300k",
    };
    Ok(IdealGasKineticsV1::new(
        bank, identity, 1, reactions, 1e-12,
    )?)
}

type Side = &'static [(&'static str, u32)];
type Efficiencies = &'static [(&'static str, f64)];
struct SourceSpecies {
    /// Ar, C, H, N, O in the canonical element order.
    atoms: [i128; 5],
    temperatures: [f64; 3],
    nasa7: [[f64; 7]; 2],
    note: &'static str,
}
enum SourcePressure {
    Independent,
    ThirdBody(Efficiencies),
    Falloff {
        low: [f64; 3],
        efficiencies: Efficiencies,
        troe: Option<[f64; 4]>,
    },
}
struct SourceReaction {
    reactants: Side,
    products: Side,
    /// Original A, b and Ea in cm/mol/s and cal/mol.
    rate: [f64; 3],
    pressure: SourcePressure,
    reversible: bool,
}

fn side(source: Side) -> Result<Vec<(SpeciesId, u32)>, ChemistryError> {
    source
        .iter()
        .map(|&(name, n)| Ok((SpeciesId::new(name)?, n)))
        .collect()
}
fn colliders(source: Efficiencies) -> Result<ThirdBodyEfficienciesV1, Gri30MechanismError> {
    Ok(ThirdBodyEfficienciesV1::new(
        1.0,
        source
            .iter()
            .map(|&(name, value)| Ok((SpeciesId::new(name)?, value)))
            .collect::<Result<Vec<_>, ChemistryError>>()?,
    )?)
}
fn arrhenius(
    [a, b, ea_cal_mol]: [f64; 3],
    total_order: u32,
) -> Result<ArrheniusRateV1, KineticsErrorV1> {
    // SI progress is mol/(m3 s). Source concentrations are mol/cm3, so the
    // coefficient multiplier is (1e-6)^(total_order-1). Generic M contributes
    // one order; explicit collision partners are already on both source sides.
    let ea = 4.184 * ea_cal_mol;
    let reference = 1000.0;
    let log_scale = fs_math::det::ln(a)
        + (total_order as f64 - 1.0) * fs_math::det::ln(1e-6)
        + b * fs_math::det::ln(reference)
        - ea / (R * reference);
    ArrheniusRateV1::new(
        Temperature::new(reference),
        ReactionProgressRateQuantityV1::new(fs_math::det::exp(log_scale)),
        b,
        MolarEnergyQuantityV1::new(ea),
    )
}

fn species_model(
    index: usize,
    treatment: Gri30Calorics,
) -> Result<Nasa9StandardStateModelV1, ThermochemErrorV1> {
    let source = &data::SPECIES[index];
    let mut regions = Vec::with_capacity(2);
    for region in 0..2 {
        let original = source.nasa7[region];
        let coefficients = [
            0.0,
            0.0,
            original[0],
            original[1],
            original[2],
            original[3],
            original[4],
            original[5],
            original[6],
        ];
        let mut parameters = BTreeMap::new();
        for ((key, dims), value) in required_nasa9_parameter_dims()
            .into_iter()
            .take(9)
            .zip(coefficients)
        {
            parameters.insert(key.to_owned(), LawParameter { value, dims });
        }
        parameters.insert(
            "reference_pressure".into(),
            LawParameter {
                value: 101325.0,
                dims: Pressure::DIMS,
            },
        );
        let lo = source.temperatures[region];
        let hi = source.temperatures[region + 1];
        let mut card = ConstitutiveModelCard {
            law: LawId("nasa9-standard-state".into()), law_version: 1, parameters,
            state_schema_version: 0, initial_state: InitialStatePolicy::ZeroInternalState,
            validity: Default::default(), sources: vec![],
            provenance: Provenance {
                source: format!("GRI-Mech 3.0 (1999-07-30); Cantera v3.2.0 {SOURCE_COMMIT}; data/gri30.yaml Git blob {SOURCE_BLOB}; original {} NASA7 {}--{} K; source note {}", CANONICAL_SPECIES[index], lo, hi, source.note),
                license: "Original GRI-Mech attribution retained; Cantera BSD-3-Clause notice in kinetics/gri30/LICENSE-Cantera.txt".into(), artifact: None,
            },
        };
        card.validity = card.validity.with("T", lo, hi);
        regions.push(Nasa9RegionV1::from_card(card)?);
    }
    let model = Nasa9StandardStateModelV1::new(
        SpeciesId::new(CANONICAL_SPECIES[index]).expect("fixed source species identifier"),
        MolarMass::new(MOLAR_MASSES_KG_PER_MOL[index]),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(101325.0),
            ElementalReferenceIdV1::new("cantera-gri30-original-reference-elements")?,
        )?,
        regions,
    )?;
    match treatment {
        Gri30Calorics::Original => Ok(model),
        Gri30Calorics::ContinuousAt300K => model.with_caloric_continuity(
            Temperature::new(300.0),
            MolarEnergyQuantityV1::new(0.2),
            MolarThermalQuantityV1::new(2e-4),
        ),
    }
}

#[cfg(test)]
mod tests;
