//! Reciprocal free-molecular gas exchange with its molecular energy flux.
//!
//! An incident Maxwellian beam carries n_dot=A*tau*p_i/sqrt(2*pi*M_i*R*T)
//! mol/s. Both ORIGINAL directional beams survive, including when their net
//! mass vanishes. Each crossing molecule carries h_i(T)-R*T/(2*M_i) J/kg:
//! translational selection raises the bulk internal energy by R*T/(2*M_i).
//! Polyatomic use assumes equilibrated internal modes independent of velocity
//! and transmission independent of molecular energy AND internal state.
//!
//! Primary monatomic kinetic theory and reciprocity: Gaspard and Andrieux,
//! "Nonlinear transport effects in mass separation by effusion", J. Stat.
//! Mech. (2011) P03024, https://gaspard.pierre.web.ulb.be/GA.JSM.11.pdf.
//! This is a thin, ballistic aperture without wall heat, work, accommodation,
//! adsorption or chemical reaction. It is not a long-pore thermal-creep law,
//! continuum nozzle, DGM correction, or inferred intrinsic permeability.
//!
//! Reservoirs are stationary, well mixed and locally Maxwellian with negligible
//! bulk drift. Crossing momentum is thermalized in reservoirs or absorbed by
//! fixed supports; no moving-interface work or momentum state is claimed.
//! The caller supplies a free-molecular qualification covering EVERY mixture
//! of the retained gas catalog over the entire authored T/p support. No mean
//! free path, thinness, species selectivity or empirical accuracy is inferred.

use super::{
    ClosedGasParcel, GasEnergyError, SpeciesId, Sum, catalog, invalid, poll, positive,
    source_enthalpy, source_error,
};
use crate::state_point::ResolvedMaterialStatePoint;
use crate::thermochemical_gas::ThermochemicalGasModel;
use fs_matdb::EvaluationDecision;
use fs_math::det;
use fs_qty::{
    Dims, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use std::{collections::BTreeMap, sync::Arc};

pub const MAX_EFFUSIVE_SPECIES: usize = 32;
pub const EFFUSIVE_TRANSMISSION_DIMS: Dims = Dims([0, 0, 0, 0, 0, 0]);
pub const EFFUSIVE_TRANSMISSION_PROPERTY: &str = "effusive_transmission_probability";

/// Caller-certified support for all mixtures of the bound original catalog.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffusiveSupport {
    pub temperature_k: (f64, f64),
    pub pressure_pa: (f64, f64),
}
impl EffusiveSupport {
    fn validate(self) -> Result<(), GasEnergyError> {
        if ![self.temperature_k, self.pressure_pa]
            .iter()
            .all(|&(lo, hi)| positive(lo) && positive(hi) && lo <= hi)
        {
            return Err(invalid(
                "effusive support needs positive ordered absolute T/p bounds",
            ));
        }
        Ok(())
    }
    fn admit(self, temperature: f64, pressure: f64) -> Result<(), GasEnergyError> {
        if !positive(temperature)
            || !positive(pressure)
            || temperature < self.temperature_k.0
            || temperature > self.temperature_k.1
            || pressure < self.pressure_pa.0
            || pressure > self.pressure_pa.1
        {
            return Err(invalid(
                "gas endpoint is outside the authored effusive T/p support",
            ));
        }
        Ok(())
    }
}

/// Authored geometric and collisional qualification, independent of exposed area.
///
/// The minimum mean free path must cover every mixture/subset of the interface's
/// exact original gas catalog throughout EffusiveSupport. Diameter/lambda <=0.1
/// and thickness/diameter <=0.01 are conservative admission gates, not measured
/// accuracy or rigorous error bounds. The caller may require tighter limits.
/// Query area is total exposed area of apertures satisfying this qualification.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusiveRegime {
    aperture_diameter_m: f64,
    wall_thickness_m: f64,
    minimum_mean_free_path_m: f64,
    provenance: String,
}
impl EffusiveRegime {
    pub fn try_new(
        aperture_diameter_m: f64,
        wall_thickness_m: f64,
        minimum_mean_free_path_m: f64,
        provenance: impl Into<String>,
    ) -> Result<Self, GasEnergyError> {
        let provenance = provenance.into();
        if !positive(aperture_diameter_m)
            || !positive(minimum_mean_free_path_m)
            || !wall_thickness_m.is_finite()
            || wall_thickness_m < 0.0
            || aperture_diameter_m / minimum_mean_free_path_m > 0.1
            || wall_thickness_m / aperture_diameter_m > 0.01
            || !valid_provenance(&provenance)
        {
            return Err(invalid(
                "effusion needs a thin ballistic aperture and bounded explicit qualification provenance",
            ));
        }
        Ok(Self {
            aperture_diameter_m,
            wall_thickness_m,
            minimum_mean_free_path_m,
            provenance,
        })
    }
    pub const fn aperture_diameter_m(&self) -> f64 {
        self.aperture_diameter_m
    }
    pub const fn wall_thickness_m(&self) -> f64 {
        self.wall_thickness_m
    }
    pub const fn minimum_mean_free_path_m(&self) -> f64 {
        self.minimum_mean_free_path_m
    }
    pub fn provenance(&self) -> &str {
        &self.provenance
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum SourceAxis {
    Temperature,
    Pressure,
}
fn source_axis(axis: &str, quantity: Option<QuantitySpec>) -> Result<SourceAxis, GasEnergyError> {
    let (kind, dims, owned) = match axis {
        "T" | "temperature" => (
            QuantityKind::AbsoluteTemperature,
            Temperature::DIMS,
            SourceAxis::Temperature,
        ),
        "pressure" => (QuantityKind::Pressure, Pressure::DIMS, SourceAxis::Pressure),
        _ => {
            return Err(invalid(
                "effusive transmission source has an unowned condition axis",
            ));
        }
    };
    let semantic = QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    if !quantity.is_none_or(|q| q == QuantitySpec::dimensional(dims) || q == semantic) {
        return Err(invalid(
            "effusive transmission source has incompatible absolute T/p semantics",
        ));
    }
    Ok(owned)
}
fn valid_provenance(value: &str) -> bool {
    !value.trim().is_empty() && value.len() <= 4096
}

/// Reciprocal transmission probability for one explicitly named source species.
///
/// A single probability applies in both directions and is independent of speed,
/// internal molecular state, temperature, pressure and mixture composition over
/// its admitted support. These are authored physical assumptions; a dimensionless
/// source value or a species label alone cannot establish them.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusiveSpeciesTransmission {
    species: SpeciesId,
    probability: f64,
    provenance: String,
    source: Option<Arc<ResolvedMaterialStatePoint>>,
    source_property: Option<String>,
    source_bounds: Vec<(SourceAxis, (f64, f64))>,
}
impl EffusiveSpeciesTransmission {
    /// Explicit model input; zero closes this species and one is an open aperture.
    pub fn try_new(
        species: SpeciesId,
        probability: f64,
        provenance: impl Into<String>,
    ) -> Result<Self, GasEnergyError> {
        let provenance = provenance.into();
        if !probability.is_finite()
            || !(0.0..=1.0).contains(&probability)
            || !valid_provenance(&provenance)
        {
            return Err(invalid(
                "effusive transmission needs a probability in [0,1] and bounded provenance",
            ));
        }
        Ok(Self {
            species,
            probability,
            provenance,
            source: None,
            source_property: None,
            source_bounds: Vec::new(),
        })
    }

    /// Bind an original resolved, validity-wide scalar constant in dimensionless
    /// probability units. Retain its source receipt, uncertainty and identity.
    /// Only named absolute T/p applicability axes are owned here; composition,
    /// humidity, stress, time and other conditions need an additional owner.
    /// The explicit provenance must qualify reciprocal, energy-independent
    /// transmission for this species, not merely repeat a material name.
    pub fn from_material(
        species: SpeciesId,
        source: ResolvedMaterialStatePoint,
        property_name: &str,
        provenance: impl Into<String>,
    ) -> Result<Self, GasEnergyError> {
        let property = source
            .property(property_name)
            .ok_or(invalid("missing selected effusive transmission property"))?;
        if property.requirement().quantity()
            != QuantitySpec::dimensional(EFFUSIVE_TRANSMISSION_DIMS)
            || !property.source_axes().is_empty()
            || property.answer().receipt.decision != EvaluationDecision::ConstantWithinValidity
        {
            return Err(invalid(
                "effusive transmission requires a validity-wide dimensionless scalar constant",
            ));
        }
        let mut result = Self::try_new(species, property.value_si(), provenance)?;
        let mut seen = [false; 2];
        for (axis, value) in source.query_point() {
            let owned = source_axis(axis, source.axis_quantities().get(axis).copied())?;
            let index = usize::from(matches!(owned, SourceAxis::Pressure));
            if seen[index] || !positive(*value) {
                return Err(invalid(
                    "effusive source needs unique positive absolute T/p coordinates",
                ));
            }
            seen[index] = true;
        }
        for (axis, &(lo, hi)) in property.answer().evidence.model.validity.bounds() {
            let owned = source_axis(
                axis,
                property.answer().receipt.axis_quantities.get(axis).copied(),
            )?;
            if !positive(lo) || !positive(hi) || lo > hi {
                return Err(invalid(
                    "effusive material validity requires positive ordered absolute T/p bounds",
                ));
            }
            result.source_bounds.push((owned, (lo, hi)));
        }
        result.source_property = Some(property_name.to_owned());
        result.source = Some(Arc::new(source));
        Ok(result)
    }
    pub const fn species(&self) -> &SpeciesId {
        &self.species
    }
    pub const fn probability(&self) -> f64 {
        self.probability
    }
    pub fn provenance(&self) -> &str {
        &self.provenance
    }
    pub fn source(&self) -> Option<&ResolvedMaterialStatePoint> {
        self.source.as_deref()
    }
    pub fn source_property(&self) -> Option<&str> {
        self.source_property.as_deref()
    }
    fn admit(&self, temperature: f64, pressure: f64) -> Result<(), GasEnergyError> {
        for &(axis, (lo, hi)) in &self.source_bounds {
            let value = match axis {
                SourceAxis::Temperature => temperature,
                SourceAxis::Pressure => pressure,
            };
            if value < lo || value > hi {
                return Err(invalid(
                    "actual gas endpoint is outside original effusive transmission source support",
                ));
            }
        }
        Ok(())
    }
}

/// Immutable gas-qualified interface; composition changes preserve source laws.
/// The catalog's initial fractions identify sources, not a privileged mixture:
/// qualification covers all mixtures/subsets over the entire authored T/p box.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusiveInterface {
    gas_catalog: Arc<ThermochemicalGasModel>,
    transmissions: Arc<[EffusiveSpeciesTransmission]>,
    support: EffusiveSupport,
    regime: EffusiveRegime,
}
impl EffusiveInterface {
    pub fn try_new(
        gas_catalog: ThermochemicalGasModel,
        mut transmissions: Vec<EffusiveSpeciesTransmission>,
        support: EffusiveSupport,
        regime: EffusiveRegime,
    ) -> Result<Self, GasEnergyError> {
        support.validate()?;
        let components = gas_catalog.thermodynamics().components();
        if components.is_empty()
            || components.len() > MAX_EFFUSIVE_SPECIES
            || transmissions.len() != components.len()
        {
            return Err(invalid(
                "effusive transmission must cover the bounded original gas catalog exactly once",
            ));
        }
        transmissions.sort_by(|a, b| a.species.cmp(&b.species));
        if !transmissions
            .iter()
            .map(|t| &t.species)
            .eq(components.iter().map(|s| s.species()))
        {
            return Err(invalid(
                "effusive transmission must cover each original gas species exactly once",
            ));
        }
        Ok(Self {
            gas_catalog: Arc::new(gas_catalog),
            transmissions: transmissions.into(),
            support,
            regime,
        })
    }
    pub fn gas_catalog(&self) -> &ThermochemicalGasModel {
        &self.gas_catalog
    }
    pub fn transmissions(&self) -> &[EffusiveSpeciesTransmission] {
        &self.transmissions
    }
    pub const fn support(&self) -> EffusiveSupport {
        self.support
    }
    pub const fn regime(&self) -> &EffusiveRegime {
        &self.regime
    }
}

/// One nonnegative gross beam from its unchanged original reservoir.
/// Inactive/absent species have exact zero rates and no caloric query.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EffusiveBeam {
    pub molar_rate_mol_s: f64,
    pub mass_rate_kg_s: f64,
    pub source_enthalpy_j_kg: Option<f64>,
    pub crossing_energy_j_kg: Option<f64>,
    /// Kept separately so a large formation reference cannot erase this term.
    pub additional_specific_energy_j_kg: Option<f64>,
    /// Total molecular energy, which may be negative in a valid source datum.
    pub energy_rate_w: f64,
    /// Correction relative to this beam's original bulk species enthalpy.
    pub additional_energy_rate_w: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct EffusiveSpeciesFlux {
    pub species: SpeciesId,
    pub left_to_right: EffusiveBeam,
    pub right_to_left: EffusiveBeam,
}
#[derive(Clone, Debug, PartialEq)]
pub struct EffusiveFlux {
    pub interface: EffusiveInterface,
    pub species: Vec<EffusiveSpeciesFlux>,
    pub net_molar_rate_left_to_right_mol_s: f64,
    pub net_mass_rate_left_to_right_kg_s: f64,
    /// Stable net energy uses the original NASA relative-enthalpy owner.
    pub energy_rate_left_to_right_w: f64,
    /// Diagnostic rate; a finite update multiplies each actual rounded dm by
    /// its beam's additional_specific_energy_j_kg before the single gas solve.
    pub additional_energy_rate_left_to_right_w: f64,
}

impl ClosedGasParcel {
    /// Sample finite-reservoir rates without modifying either gas state.
    ///
    /// Send BOTH gross beams through aggregate original-species withdrawal
    /// admission. A finite update must combine their original enthalpy and
    /// additional energy before reconstructing each parcel ONCE. Do not apply
    /// nozzle pressure-order or isothermal partial-pressure gates to effusion.
    /// Its zero-mass-flow condition is p_i,left/sqrt(T_left)=p_i,right/sqrt(T_right);
    /// heat generally remains nonzero there. Finite endpoint entropy admission
    /// and physical timestep refinement still belong to the transport owner.
    ///
    /// Zero geometric area is exactly inactive without querying gas/calibration
    /// support or joining inactive sources. A calibrated zero transmission still
    /// requires its source support; closure is not extrapolated beyond that box.
    pub fn effusive_flux_to(
        &self,
        right: &Self,
        interface: &EffusiveInterface,
        area_m2: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EffusiveFlux, GasEnergyError> {
        poll(&mut cancelled)?;
        if !area_m2.is_finite() || area_m2 < 0.0 {
            return Err(invalid("effusive area must be finite and nonnegative"));
        }
        let mut result = EffusiveFlux {
            interface: interface.clone(),
            species: Vec::new(),
            net_molar_rate_left_to_right_mol_s: 0.0,
            net_mass_rate_left_to_right_kg_s: 0.0,
            energy_rate_left_to_right_w: 0.0,
            additional_energy_rate_left_to_right_w: 0.0,
        };
        if area_m2 == 0.0 {
            poll(&mut cancelled)?;
            return Ok(result);
        }
        let (left_state, right_state) = (self.state().state(), right.state().state());
        let (tl, tr) = (left_state.temperature, right_state.temperature);
        let (pl, pr) = (left_state.pressure, right_state.pressure);
        for (t, p) in [(tl, pl), (tr, pr)] {
            interface.support.admit(t, p)?;
        }
        if [self, right]
            .iter()
            .any(|p| p.model().thermodynamics().components().len() > MAX_EFFUSIVE_SPECIES)
        {
            return Err(invalid(
                "effusive original gas exceeds the species work bound",
            ));
        }
        let joined = catalog(self, right)?;
        let components = joined.thermodynamics().components();
        if components.len() > MAX_EFFUSIVE_SPECIES
            || components.iter().any(|s| {
                !interface
                    .transmissions
                    .iter()
                    .any(|t| &t.species == s.species())
            })
        {
            return Err(invalid(
                "effusive interface qualification does not cover these gas species",
            ));
        }
        let fractions = |gas: &ThermochemicalGasModel| -> BTreeMap<SpeciesId, f64> {
            gas.thermodynamics()
                .components()
                .iter()
                .zip(gas.thermodynamics().mole_composition().fractions())
                .map(|(s, &x)| (s.species().clone(), x))
                .collect()
        };
        admit_calibration(&joined, interface.gas_catalog())?;
        let (left_x, right_x) = (fractions(self.model()), fractions(right.model()));
        let (mut moles, mut mass, mut energy, mut extra) = (
            Sum::default(),
            Sum::default(),
            Sum::default(),
            Sum::default(),
        );
        for transmission in interface.transmissions.iter() {
            poll(&mut cancelled)?;
            let id = &transmission.species;
            let Some(species) = components.iter().find(|s| s.species() == id) else {
                continue;
            };
            let molar_mass = species.molar_mass().value();
            let lx = left_x.get(id).copied().unwrap_or(0.0);
            let rx = right_x.get(id).copied().unwrap_or(0.0);
            for (t, p) in [(tl, pl), (tr, pr)] {
                transmission.admit(t, p)?;
            }
            let l = beam(self, id, lx, molar_mass, transmission.probability, area_m2)?;
            poll(&mut cancelled)?;
            let r = beam(right, id, rx, molar_mass, transmission.probability, area_m2)?;
            let delta_mass = l.mass_rate_kg_s - r.mass_rate_kg_s;
            let paired_mass = l.mass_rate_kg_s.min(r.mass_rate_kg_s);
            moles.add(l.molar_rate_mol_s - r.molar_rate_mol_s)?;
            mass.add(delta_mass)?;
            let mut net_energy = Sum::default();
            let mut net_extra = Sum::default();
            if delta_mass != 0.0 {
                let donor = if delta_mass > 0.0 { l } else { r };
                net_energy.add(signed_product(
                    delta_mass,
                    donor
                        .source_enthalpy_j_kg
                        .ok_or(invalid("active effusive beam is missing original enthalpy"))?,
                )?)?;
                net_extra.add(signed_product(
                    delta_mass,
                    donor.additional_specific_energy_j_kg.ok_or(invalid(
                        "active effusive beam is missing molecular energy correction",
                    ))?,
                )?)?;
            }
            if paired_mass > 0.0 && tl != tr {
                let dh = species
                    .molar_enthalpy_change(Temperature::new(tr), Temperature::new(tl))
                    .map_err(source_error)?
                    .value()
                    / molar_mass;
                if !dh.is_finite() {
                    return Err(invalid(
                        "effusive relative source enthalpy is unrepresentable",
                    ));
                }
                net_energy.add(signed_product(paired_mass, dh)?)?;
                let correction = signed_product(-0.5 * R, tl - tr)?;
                let paired_extra =
                    positive_product_ratio(&[paired_mass, correction.abs()], &[molar_mass])?;
                net_extra.add(paired_extra.copysign(correction))?;
            }
            net_energy.add(net_extra.value)?;
            energy.add(net_energy.value)?;
            extra.add(net_extra.value)?;
            result.species.push(EffusiveSpeciesFlux {
                species: id.clone(),
                left_to_right: l,
                right_to_left: r,
            });
        }
        poll(&mut cancelled)?;
        result.net_molar_rate_left_to_right_mol_s = moles.value;
        result.net_mass_rate_left_to_right_kg_s = mass.value;
        result.energy_rate_left_to_right_w = energy.value;
        result.additional_energy_rate_left_to_right_w = extra.value;
        Ok(result)
    }
}

// Compare the existing immutable physical sources directly. A calibration
// covering every mixture must not reconstruct/renormalize a second composition.
fn admit_calibration(
    gas: &ThermochemicalGasModel,
    bound: &ThermochemicalGasModel,
) -> Result<(), GasEnergyError> {
    for (species, transport) in gas
        .thermodynamics()
        .components()
        .iter()
        .zip(gas.transport_materials())
    {
        let index = bound
            .thermodynamics()
            .components()
            .binary_search_by(|candidate| candidate.species().cmp(species.species()))
            .map_err(|_| invalid("effusive qualification is missing an active gas species"))?;
        if &bound.thermodynamics().components()[index] != species
            || &bound.transport_materials()[index] != transport
        {
            return Err(invalid(
                "effusive qualification substitutes an original caloric or transport source",
            ));
        }
        let original = gas
            .source_packs()
            .iter()
            .find(|receipt| &receipt.species == species.species());
        let retained = bound
            .source_packs()
            .iter()
            .find(|receipt| &receipt.species == species.species());
        if original != retained {
            return Err(invalid(
                "effusive qualification substitutes an original source-pack binding",
            ));
        }
    }
    let active_water = bound.water_species().filter(|id| {
        gas.thermodynamics()
            .components()
            .iter()
            .any(|species| species.species() == *id)
    });
    if gas.water_species() != active_water {
        return Err(invalid(
            "effusive qualification changes the explicit active water identity",
        ));
    }
    Ok(())
}

fn beam(
    parcel: &ClosedGasParcel,
    species: &SpeciesId,
    fraction: f64,
    molar_mass: f64,
    transmission: f64,
    area: f64,
) -> Result<EffusiveBeam, GasEnergyError> {
    if fraction == 0.0 || transmission == 0.0 {
        return Ok(EffusiveBeam::default());
    }
    let state = parcel.state().state();
    let (temperature, pressure) = (state.temperature, state.pressure);
    let molar = positive_product_ratio(
        &[area, transmission, fraction, pressure],
        &[
            det::sqrt(2.0 * core::f64::consts::PI * R),
            det::sqrt(molar_mass),
            det::sqrt(temperature),
        ],
    )?;
    let mass = positive_product_ratio(&[molar, molar_mass], &[])?;
    let h = source_enthalpy(parcel, species)?;
    let correction = -positive_product_ratio(&[R, temperature], &[2.0, molar_mass])?;
    let crossing = h + correction;
    if !crossing.is_finite() {
        return Err(invalid("effusive crossing energy is unrepresentable"));
    }
    Ok(EffusiveBeam {
        molar_rate_mol_s: molar,
        mass_rate_kg_s: mass,
        source_enthalpy_j_kg: Some(h),
        crossing_energy_j_kg: Some(crossing),
        additional_specific_energy_j_kg: Some(correction),
        energy_rate_w: signed_product(mass, crossing)?,
        additional_energy_rate_w: signed_product(mass, correction)?,
    })
}

fn signed_product(a: f64, b: f64) -> Result<f64, GasEnergyError> {
    if !a.is_finite() || !b.is_finite() {
        return Err(invalid("nonfinite effusive energy factor"));
    }
    if a == 0.0 || b == 0.0 {
        return Ok(0.0);
    }
    Ok(positive_product_ratio(&[a.abs(), b.abs()], &[])?.copysign(
        if a.is_sign_negative() == b.is_sign_negative() {
            1.0
        } else {
            -1.0
        },
    ))
}

/// Positive product/quotient with normalized binary factors. Intermediate
/// overflow/underflow does not close a real path; an unrepresentable final
/// positive result refuses explicitly. At most seven factors occur per beam.
fn positive_product_ratio(numerator: &[f64], denominator: &[f64]) -> Result<f64, GasEnergyError> {
    fn split(mut value: f64) -> (f64, i32) {
        let mut adjustment = 0;
        if value < f64::MIN_POSITIVE {
            value *= 18_014_398_509_481_984.0;
            adjustment = -54;
        }
        let bits = value.to_bits();
        (
            f64::from_bits((bits & ((1u64 << 52) - 1)) | (1023u64 << 52)),
            ((bits >> 52) & 0x7ff) as i32 - 1023 + adjustment,
        )
    }
    let mut mantissa = 1.0;
    let mut exponent = 0i32;
    for (&value, divide) in numerator
        .iter()
        .map(|v| (v, false))
        .chain(denominator.iter().map(|v| (v, true)))
    {
        if !positive(value) {
            return Err(invalid(
                "effusive rate needs positive representable factors",
            ));
        }
        let (factor, power) = split(value);
        if divide {
            mantissa /= factor;
            exponent -= power;
        } else {
            mantissa *= factor;
            exponent += power;
        }
        if mantissa >= 2.0 {
            mantissa *= 0.5;
            exponent += 1;
        } else if mantissa < 1.0 {
            mantissa *= 2.0;
            exponent -= 1;
        }
    }
    let result = if (-1022..=1023).contains(&exponent) {
        mantissa * f64::from_bits(((exponent + 1023) as u64) << 52)
    } else if (-1075..-1022).contains(&exponent) {
        (mantissa * f64::MIN_POSITIVE) * f64::from_bits(((exponent + 1022 + 1023) as u64) << 52)
    } else {
        return Err(invalid(
            "positive effusive rate or energy is unrepresentable",
        ));
    };
    if !positive(result) {
        return Err(invalid(
            "positive effusive rate or energy vanished or overflowed",
        ));
    }
    Ok(result)
}

#[cfg(test)]
#[path = "effusion/tests.rs"]
mod tests;
