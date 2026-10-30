//! Complete ideal-gas `ohmech` from Cantera v3.2.0 `data/h2o2.yaml`.
//!
//! The source describes this as the hydrogen/oxygen submechanism extracted from
//! GRI-Mech 3.0, with N2 added. All ten species and 29 reversible channels are
//! retained, including duplicate channels, explicit collider species, five
//! generic third-body laws and one Troe falloff law. This is not the complete
//! GRI methane mechanism. No combustion rate or radical seed is inferred.
//!
//! Source revision: <https://github.com/Cantera/cantera/blob/4a8358eb80cfeb50474386b5f9ec0b3a83519889/data/h2o2.yaml>.
//! Original NASA7 coefficients embed exactly in the original NASA9 owner with
//! its first two coefficients zero. The common source interval is 300..3500 K.
//! Original thermodynamics remains available; crossing its independently fitted
//! region seam requires the explicitly selected bounded continuity treatment.
//! Forward laws retain the original cm/mol/s and cal/mol conversion. Reverse
//! laws, affinities, Jacobians and pressure dependence remain in `kinetics`.
//! No time integrator, transport coefficients, real-gas phase or experimental
//! accuracy claim is supplied by loading this mechanism.

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

/// Pinned Cantera release commit, not a mutable "latest" source.
pub const SOURCE_COMMIT: &str = "4a8358eb80cfeb50474386b5f9ec0b3a83519889";
/// Exact original `data/h2o2.yaml` Git blob.
pub const SOURCE_BLOB: &str = "65e6321a60bafb1feb4f3b1be5347e6d1e89c937";
/// Exact original `src/thermo/Elements.cpp` Git blob supplying atomic weights.
pub const ATOMIC_WEIGHTS_BLOB: &str = "ebe9da8b1a3e7944d1de070851b36fdb858ad8d4";
/// Complete axis returned by the original species-bank owner.
pub const CANONICAL_SPECIES: [&str; 10] =
    ["AR", "H", "H2", "H2O", "H2O2", "HO2", "N2", "O", "O2", "OH"];
/// Source atomic weights (H=1.008, O=15.999, N=14.007, Ar=39.95 g/mol) summed
/// on `CANONICAL_SPECIES`, in kg/mol.
pub const MOLAR_MASSES_KG_PER_MOL: [f64; 10] = [
    0.03995, 0.001008, 0.002016, 0.018015, 0.034014, 0.033006, 0.028014, 0.015999, 0.031998,
    0.017007,
];
/// Common lower bound; the source AR and N2 cards begin at 300 K.
pub const MINIMUM_TEMPERATURE_K: f64 = 300.0;
/// Common upper bound; the reactive species cards end at 3500 K.
pub const MAXIMUM_TEMPERATURE_K: f64 = 3500.0;

/// Explicit choice of original versus reconciled piecewise caloric references.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HydrogenCalorics {
    /// Exact source NASA7 h/s, including the original small seam jumps.
    Original,
    /// Preserve the entire low-temperature region and add constant h/s offsets
    /// to the upper region using the original continuity owner. Budgets are
    /// 0.01 J/mol and 2e-5 J/(mol K); these bound numerical source adjustment,
    /// not empirical uncertainty. Cp and source cards are unchanged.
    ContinuousAt300K,
}

/// Original owner refusing construction; no partial mechanism is returned.
#[derive(Clone, Debug, PartialEq)]
pub enum HydrogenMechanismError {
    /// Species identities or exact elemental/charge assembly refused.
    Chemistry(ChemistryError),
    /// Original standard-state or explicit continuity owner refused.
    Thermodynamics(ThermochemErrorV1),
    /// Original neutral species bank refused.
    Equilibrium(CompositionEquilibriumErrorV1),
    /// Original elementary kinetic owner refused.
    Kinetics(KineticsErrorV1),
}
impl core::fmt::Display for HydrogenMechanismError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Cantera 3.2.0 hydrogen source: {self:?}")
    }
}
impl std::error::Error for HydrogenMechanismError {}
impl From<ChemistryError> for HydrogenMechanismError {
    fn from(e: ChemistryError) -> Self {
        Self::Chemistry(e)
    }
}
impl From<ThermochemErrorV1> for HydrogenMechanismError {
    fn from(e: ThermochemErrorV1) -> Self {
        Self::Thermodynamics(e)
    }
}
impl From<CompositionEquilibriumErrorV1> for HydrogenMechanismError {
    fn from(e: CompositionEquilibriumErrorV1) -> Self {
        Self::Equilibrium(e)
    }
}
impl From<KineticsErrorV1> for HydrogenMechanismError {
    fn from(e: KineticsErrorV1) -> Self {
        Self::Kinetics(e)
    }
}

/// Load all source channels once into the existing immutable mechanism.
///
/// Construction is bounded by ten species, twenty NASA regions and 29 channels;
/// repeated evaluations use the original prepared kinetic workspace. Molecular
/// mass admission uses 1e-12 relative, without adjusting any species mass.
pub fn cantera_h2o2(
    calorics: HydrogenCalorics,
) -> Result<IdealGasKineticsV1, HydrogenMechanismError> {
    let models = (0..10)
        .map(|i| species_model(i, calorics))
        .collect::<Result<Vec<_>, _>>()?;
    let ids = models
        .iter()
        .map(|m| m.species().clone())
        .collect::<Vec<_>>();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            ["Ar", "H", "N", "O"]
                .iter()
                .map(|s| ElementId::new(*s))
                .collect::<Result<Vec<_>, _>>()?,
            ids.clone(),
            vec![
                vec![1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                vec![0, 1, 2, 2, 2, 1, 0, 0, 0, 1],
                vec![0, 0, 0, 0, 0, 0, 2, 0, 0, 0],
                vec![0, 0, 0, 1, 2, 2, 0, 1, 2, 1],
            ],
        )?,
        ChargeVector::new(ids, vec![0; 10])?,
        models,
    )?;
    let mut reactions = Vec::with_capacity(REACTIONS.len());
    for (i, &(reactants, products, source_rate)) in REACTIONS.iter().enumerate() {
        let number = i + 1;
        let order: u32 = reactants.iter().map(|(_, exponent)| exponent).sum();
        let generic_collider = matches!(number, 1 | 2 | 6 | 12 | 15);
        let forward = arrhenius(source_rate, order + u32::from(generic_collider))?;
        let pressure = match number {
            1 => PressureDependenceV1::ThirdBody(colliders(&[
                ("H2", 2.4),
                ("H2O", 15.4),
                ("AR", 0.83),
            ])?),
            2 => PressureDependenceV1::ThirdBody(colliders(&[
                ("H2", 2.0),
                ("H2O", 6.0),
                ("AR", 0.7),
            ])?),
            6 => PressureDependenceV1::ThirdBody(colliders(&[
                ("O2", 0.0),
                ("H2O", 0.0),
                ("N2", 0.0),
                ("AR", 0.0),
            ])?),
            12 => PressureDependenceV1::ThirdBody(colliders(&[
                ("H2", 0.0),
                ("H2O", 0.0),
                ("AR", 0.63),
            ])?),
            15 => PressureDependenceV1::ThirdBody(colliders(&[
                ("H2", 0.73),
                ("H2O", 3.65),
                ("AR", 0.38),
            ])?),
            22 => PressureDependenceV1::Falloff {
                low_pressure_rate: arrhenius([2.3e18, -0.9, -1700.0], order + 1)?,
                colliders: colliders(&[("H2", 2.0), ("H2O", 6.0), ("AR", 0.7)])?,
                broadening: FalloffBroadeningV1::Troe(TroeFalloffV1::new(
                    0.7346,
                    Temperature::new(94.0),
                    Temperature::new(1756.0),
                    Some(Temperature::new(5182.0)),
                )?),
            },
            _ => PressureDependenceV1::Independent,
        };
        reactions.push(
            ElementaryReactionV1::new(
                ReactionId::new(format!("cantera-h2o2-r{number:02}"))?,
                side(reactants)?,
                side(products)?,
                forward,
            )?
            .with_pressure_dependence(pressure),
        );
    }
    let identity = match calorics {
        HydrogenCalorics::Original => "cantera-v3.2.0-h2o2-ohmech-original",
        HydrogenCalorics::ContinuousAt300K => "cantera-v3.2.0-h2o2-ohmech-continuous-300k",
    };
    Ok(IdealGasKineticsV1::new(
        bank, identity, 1, reactions, 1e-12,
    )?)
}

fn side(source: Side) -> Result<Vec<(SpeciesId, u32)>, ChemistryError> {
    source
        .iter()
        .map(|&(s, n)| Ok((SpeciesId::new(s)?, n)))
        .collect()
}
fn colliders(source: &[(&str, f64)]) -> Result<ThirdBodyEfficienciesV1, HydrogenMechanismError> {
    Ok(ThirdBodyEfficienciesV1::new(
        1.0,
        source
            .iter()
            .map(|&(s, e)| Ok((SpeciesId::new(s)?, e)))
            .collect::<Result<Vec<_>, ChemistryError>>()?,
    )?)
}
fn arrhenius(
    [a, b, ea_cal_mol]: [f64; 3],
    total_order: u32,
) -> Result<ArrheniusRateV1, KineticsErrorV1> {
    // q_SI = 1e6 q_cm = A_cm (1e-6)^(order-1) product(c_SI^order).
    // Generic third-body and low-pressure orders include their extra collider.
    // Explicit collider species remain on BOTH stoichiometric sides and are
    // already counted in total_order. A source calorie is exactly 4.184 J.
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
    i: usize,
    treatment: HydrogenCalorics,
) -> Result<Nasa9StandardStateModelV1, ThermochemErrorV1> {
    let (minimum, maximum) = if matches!(i, 0 | 6) {
        (300.0, 5000.0)
    } else {
        (200.0, 3500.0)
    };
    let mut regions = Vec::with_capacity(2);
    for (region, (lo, hi)) in [(minimum, 1000.0), (1000.0, maximum)]
        .into_iter()
        .enumerate()
    {
        let original = NASA7[i][region];
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
        let mut card = ConstitutiveModelCard {
            law: LawId("nasa9-standard-state".into()), law_version: 1, parameters,
            state_schema_version: 0, initial_state: InitialStatePolicy::ZeroInternalState,
            validity: Default::default(), sources: vec![],
            provenance: Provenance {
                source: format!("GRI-Mech 3.0 hydrogen/oxygen submechanism; Cantera v3.2.0 {SOURCE_COMMIT}; data/h2o2.yaml Git blob {SOURCE_BLOB}; original {} NASA7 {}--{} K; source note {}", CANONICAL_SPECIES[i], lo, hi, NOTES[i]),
                license: "Original GRI-Mech attribution retained; Cantera distribution under BSD-3-Clause".into(), artifact: None,
            },
        };
        card.validity = card.validity.with("T", lo, hi);
        regions.push(Nasa9RegionV1::from_card(card)?);
    }
    // Fixed source identities have been checked when constructing the bank.
    let species = SpeciesId::new(CANONICAL_SPECIES[i]).expect("fixed source species identifier");
    let reference = ElementalReferenceIdV1::new("cantera-h2o2-original-reference-elements")?;
    let model = Nasa9StandardStateModelV1::new(
        species,
        MolarMass::new(MOLAR_MASSES_KG_PER_MOL[i]),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(101325.0),
            reference,
        )?,
        regions,
    )?;
    match treatment {
        HydrogenCalorics::Original => Ok(model),
        HydrogenCalorics::ContinuousAt300K => model.with_caloric_continuity(
            Temperature::new(300.0),
            MolarEnergyQuantityV1::new(0.01),
            MolarThermalQuantityV1::new(2e-5),
        ),
    }
}

const NOTES: [&str; 10] = [
    "120186", "L7/88", "TPIS78", "L8/89", "L7/88", "L5/89", "121286", "L1/90", "TPIS89", "RUS78",
];

// Original source arrays, canonically reordered without refitting.
#[rustfmt::skip]
const NASA7: [[[f64; 7]; 2]; 10] = [
    // AR
    [[2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.366],
     [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.366]],
    // H
    [[2.5, 7.05332819e-13, -1.99591964e-15, 2.30081632e-18, -9.27732332e-22, 2.54736599e4, -0.446682853],
     [2.50000001, -2.30842973e-11, 1.61561948e-14, -4.73515235e-18, 4.98197357e-22, 2.54736599e4, -0.446682914]],
    // H2
    [[2.34433112, 7.98052075e-3, -1.9478151e-5, 2.01572094e-8, -7.37611761e-12, -917.935173, 0.683010238],
     [3.3372792, -4.94024731e-5, 4.99456778e-7, -1.79566394e-10, 2.00255376e-14, -950.158922, -3.20502331]],
    // H2O
    [[4.19864056, -2.0364341e-3, 6.52040211e-6, -5.48797062e-9, 1.77197817e-12, -3.02937267e4, -0.849032208],
     [3.03399249, 2.17691804e-3, -1.64072518e-7, -9.7041987e-11, 1.68200992e-14, -3.00042971e4, 4.9667701]],
    // H2O2
    [[4.27611269, -5.42822417e-4, 1.67335701e-5, -2.15770813e-8, 8.62454363e-12, -1.77025821e4, 3.43505074],
     [4.16500285, 4.90831694e-3, -1.90139225e-6, 3.71185986e-10, -2.87908305e-14, -1.78617877e4, 2.91615662]],
    // HO2
    [[4.30179801, -4.74912051e-3, 2.11582891e-5, -2.42763894e-8, 9.29225124e-12, 294.80804, 3.71666245],
     [4.0172109, 2.23982013e-3, -6.3365815e-7, 1.1424637e-10, -1.07908535e-14, 111.856713, 3.78510215]],
    // N2
    [[3.298677, 1.4082404e-3, -3.963222e-6, 5.641515e-9, -2.444854e-12, -1020.8999, 3.950372],
     [2.92664, 1.4879768e-3, -5.68476e-7, 1.0097038e-10, -6.753351e-15, -922.7977, 5.980528]],
    // O
    [[3.1682671, -3.27931884e-3, 6.64306396e-6, -6.12806624e-9, 2.11265971e-12, 2.91222592e4, 2.05193346],
     [2.56942078, -8.59741137e-5, 4.19484589e-8, -1.00177799e-11, 1.22833691e-15, 2.92175791e4, 4.78433864]],
    // O2
    [[3.78245636, -2.99673416e-3, 9.84730201e-6, -9.68129509e-9, 3.24372837e-12, -1063.94356, 3.65767573],
     [3.28253784, 1.48308754e-3, -7.57966669e-7, 2.09470555e-10, -2.16717794e-14, -1088.45772, 5.45323129]],
    // OH
    [[3.99201543, -2.40131752e-3, 4.61793841e-6, -3.88113333e-9, 1.3641147e-12, 3615.08056, -0.103925458],
     [3.09288767, 5.48429716e-4, 1.26505228e-7, -8.79461556e-11, 1.17412376e-14, 3858.657, 4.4766961]],
];

type Side = &'static [(&'static str, u32)];
// Source reaction order and A [cm/mol/s], b, Ea [cal/mol]. Explicit colliders
// stay on each side; generic M and (+M) are supplied by the pressure-law match.
#[rustfmt::skip]
const REACTIONS: [(Side, Side, [f64; 3]); 29] = [
    (&[("O", 2)], &[("O2", 1)], [1.2e17, -1.0, 0.0]),
    (&[("O", 1), ("H", 1)], &[("OH", 1)], [5.0e17, -1.0, 0.0]),
    (&[("O", 1), ("H2", 1)], &[("H", 1), ("OH", 1)], [3.87e4, 2.7, 6260.0]),
    (&[("O", 1), ("HO2", 1)], &[("OH", 1), ("O2", 1)], [2.0e13, 0.0, 0.0]),
    (&[("O", 1), ("H2O2", 1)], &[("OH", 1), ("HO2", 1)], [9.63e6, 2.0, 4000.0]),
    (&[("H", 1), ("O2", 1)], &[("HO2", 1)], [2.8e18, -0.86, 0.0]),
    (&[("H", 1), ("O2", 2)], &[("HO2", 1), ("O2", 1)], [2.08e19, -1.24, 0.0]),
    (&[("H", 1), ("O2", 1), ("H2O", 1)], &[("HO2", 1), ("H2O", 1)], [1.126e19, -0.76, 0.0]),
    (&[("H", 1), ("O2", 1), ("N2", 1)], &[("HO2", 1), ("N2", 1)], [2.6e19, -1.24, 0.0]),
    (&[("H", 1), ("O2", 1), ("AR", 1)], &[("HO2", 1), ("AR", 1)], [7.0e17, -0.8, 0.0]),
    (&[("H", 1), ("O2", 1)], &[("O", 1), ("OH", 1)], [2.65e16, -0.6707, 1.7041e4]),
    (&[("H", 2)], &[("H2", 1)], [1.0e18, -1.0, 0.0]),
    (&[("H", 2), ("H2", 1)], &[("H2", 2)], [9.0e16, -0.6, 0.0]),
    (&[("H", 2), ("H2O", 1)], &[("H2", 1), ("H2O", 1)], [6.0e19, -1.25, 0.0]),
    (&[("H", 1), ("OH", 1)], &[("H2O", 1)], [2.2e22, -2.0, 0.0]),
    (&[("H", 1), ("HO2", 1)], &[("O", 1), ("H2O", 1)], [3.97e12, 0.0, 671.0]),
    (&[("H", 1), ("HO2", 1)], &[("O2", 1), ("H2", 1)], [4.48e13, 0.0, 1068.0]),
    (&[("H", 1), ("HO2", 1)], &[("OH", 2)], [8.4e13, 0.0, 635.0]),
    (&[("H", 1), ("H2O2", 1)], &[("HO2", 1), ("H2", 1)], [1.21e7, 2.0, 5200.0]),
    (&[("H", 1), ("H2O2", 1)], &[("OH", 1), ("H2O", 1)], [1.0e13, 0.0, 3600.0]),
    (&[("OH", 1), ("H2", 1)], &[("H", 1), ("H2O", 1)], [2.16e8, 1.51, 3430.0]),
    (&[("OH", 2)], &[("H2O2", 1)], [7.4e13, -0.37, 0.0]),
    (&[("OH", 2)], &[("O", 1), ("H2O", 1)], [3.57e4, 2.4, -2110.0]),
    (&[("OH", 1), ("HO2", 1)], &[("O2", 1), ("H2O", 1)], [1.45e13, 0.0, -500.0]),
    (&[("OH", 1), ("H2O2", 1)], &[("HO2", 1), ("H2O", 1)], [2.0e12, 0.0, 427.0]),
    (&[("OH", 1), ("H2O2", 1)], &[("HO2", 1), ("H2O", 1)], [1.7e18, 0.0, 2.941e4]),
    (&[("HO2", 2)], &[("O2", 1), ("H2O2", 1)], [1.3e11, 0.0, -1630.0]),
    (&[("HO2", 2)], &[("O2", 1), ("H2O2", 1)], [4.2e14, 0.0, 1.2e4]),
    (&[("OH", 1), ("HO2", 1)], &[("O2", 1), ("H2O", 1)], [5.0e15, 0.0, 1.733e4]),
];

#[cfg(test)]
mod tests;
