//! Finite Langmuir site storage with the original gas and sourced solid calorics.
//!
//! One adsorbate occupies one homogeneous population of noninteracting sites.
//! Capacity is mol/kg DRY skeleton; affinity is 1/Pa, not a dimensionless K.
//! The source reference affinity fixes the bound-phase entropy datum. Its
//! temperature dependence follows the ORIGINAL gas and bound-phase caloric
//! charts, so adsorption heat is their explicit datum-aligned difference.
//!
//! Frozen forcing integrates both gross site streams analytically. This is
//! exact for fixed gas partial pressure and surface temperature, not for a
//! finite nonisothermal bath. Original gas withdrawal, site turnover, thermal
//! changes, final affinity and combined entropy bound the accepted interval.
//! Both finite temperatures evolve; there is no implicit thermostat.
//!
//! Streams carry original bulk gas partial enthalpy. This is a continuum
//! well-mixed pore/film exchange convention, NOT a free-molecular impingement
//! model. The latter needs a different rate/energy closure (h-RT/2 alone does
//! not establish nonisothermal detailed balance). Solid pressure work,
//! competitive sites, surface diffusion, hysteresis, swelling and heat
//! conduction are outside this homogeneous rigid-solid rung.
//!
//! Langmuir kinetics and standard-state distinctions: Knopf & Ammann (2021),
//! https://doi.org/10.5194/acp-21-15725-2021 . Source applicability is required
//! independently for the actual adsorbent, gas and caloric charts.

/// Shared-site multicomponent uptake with the original source and caloric owners.
pub mod competitive;

use crate::phase::{PhaseStateError, ResolvedHeatCapacityCurve};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{
        GasEnergyError,
        parcel::{
            ClosedGasParcel,
            transfer::species::boundary::{
                GasSpeciesBoundaryCandidate, GasSpeciesBoundarySolve, SpeciesBoundaryRelativeTransfer,
                SpeciesBoundaryTransfer,
            },
        },
    },
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_math::det;
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    Nasa9StandardStateModelV1, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};
use std::sync::Arc;

pub const LANGMUIR_CAPACITY_PROPERTY: &str = "langmuir_capacity";
pub const LANGMUIR_REFERENCE_AFFINITY_PROPERTY: &str = "langmuir_reference_affinity";
pub const LANGMUIR_DESORPTION_RATE_PROPERTY: &str = "langmuir_desorption_rate";
pub const LANGMUIR_CAPACITY_DIMS: Dims = Dims([0, -1, 0, 0, 0, 1]);
pub const LANGMUIR_AFFINITY_DIMS: Dims = Dims([1, -1, 2, 0, 0, 0]);
pub const LANGMUIR_RATE_DIMS: Dims = Dims([0, 0, -1, 0, 0, 0]);

/// Capacity and reference affinity are selected validity-wide constants.
/// Desorption rate may depend on canonical temperature. Other source
/// conditions are fixed except actual total gas pressure applicability.
#[derive(Clone, Debug, PartialEq)]
pub struct AdsorptionSource {
    pub material: MaterialCard,
    pub point: QueryPoint,
    pub selection: MaterialPropertySelection,
    pub reference_temperature_k: f64,
}

/// Existing source-qualified single-phase charts. The bound-adsorbate offset
/// aligns its h datum with the ORIGINAL gas species; it is not supplied heat.
/// Negligible solid pV and rigid geometry identify solid h storage with u.
#[derive(Clone, Debug, PartialEq)]
pub struct AdsorptionCaloricBinding {
    pub dry_skeleton: ResolvedHeatCapacityCurve,
    pub bound_adsorbate: ResolvedHeatCapacityCurve,
    pub bound_enthalpy_offset_j_kg: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AdsorptionError {
    Invalid(&'static str),
    Material(MaterialStatePointError),
    Thermochemistry(String),
    Phase(PhaseStateError),
    Gas(GasEnergyError),
    InventoryLimit {
        gross_gas_fraction: f64,
        site_turnover: f64,
    },
    ThermalRange,
    NotConverged {
        residual_j: f64,
    },
    EquilibriumCrossed,
    EntropyDecreased {
        change_j_k: f64,
    },
    Cancelled,
}
impl core::fmt::Display for AdsorptionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite Langmuir adsorption: {self:?}")
    }
}
impl std::error::Error for AdsorptionError {}
impl From<MaterialStatePointError> for AdsorptionError {
    fn from(e: MaterialStatePointError) -> Self {
        Self::Material(e)
    }
}
impl From<PhaseStateError> for AdsorptionError {
    fn from(e: PhaseStateError) -> Self {
        Self::Phase(e)
    }
}
impl From<GasEnergyError> for AdsorptionError {
    fn from(e: GasEnergyError) -> Self {
        if e == GasEnergyError::Cancelled {
            Self::Cancelled
        } else {
            Self::Gas(e)
        }
    }
}
type Result<T> = core::result::Result<T, AdsorptionError>;
fn invalid(s: &'static str) -> AdsorptionError {
    AdsorptionError::Invalid(s)
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(AdsorptionError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(x: f64) -> Result<f64> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(invalid("adsorption arithmetic is unrepresentable"))
    }
}
fn temp_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn gas_failure(e: impl core::fmt::Display) -> AdsorptionError {
    AdsorptionError::Thermochemistry(e.to_string())
}

fn source_point(context: &QueryPoint, temperature: f64, pressure: f64) -> Result<QueryPoint> {
    context
        .clone()
        .with_quantity("temperature", temp_quantity(), temperature)
        .and_then(|p| {
            p.with_quantity(
                "pressure",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::Pressure,
                    ValueForm::Static,
                )),
                pressure,
            )
        })
        .map_err(|_| invalid("adsorption source needs finite typed temperature and pressure"))
}

/// Immutable physical and source identity shared by successive solid states.
#[derive(Clone, Debug, PartialEq)]
pub struct LangmuirAdsorbent {
    species: SpeciesId,
    catalog: ThermochemicalGasModel,
    gas_index: usize,
    source: AdsorptionSource,
    resolved: ResolvedMaterialStatePoint,
    calorics: AdsorptionCaloricBinding,
    capacity_mol_kg_dry: f64,
    reference_affinity_pa_inverse: f64,
    reference_heat_j_mol: f64,
    temperature_range_k: [f64; 2],
}
#[derive(Clone, Debug, PartialEq)]
pub struct LangmuirProperties {
    pub temperature_k: f64,
    pub affinity_pa_inverse: f64,
    pub desorption_rate_s_inverse: f64,
    /// h_g(T)-h_bound(T), positive for exothermic adsorption.
    pub adsorption_heat_j_mol: f64,
    pub source: ResolvedMaterialStatePoint,
}
impl LangmuirAdsorbent {
    pub fn try_new(
        catalog: &ThermochemicalGasModel,
        species: SpeciesId,
        source: AdsorptionSource,
        calorics: AdsorptionCaloricBinding,
    ) -> Result<Self> {
        let index = catalog
            .thermodynamics()
            .components()
            .iter()
            .position(|s| s.species() == &species)
            .ok_or(invalid(
                "adsorption requires a retained original gas species",
            ))?;
        let reference = source.reference_temperature_k;
        if !positive(reference) || !calorics.bound_enthalpy_offset_j_kg.is_finite() {
            return Err(invalid(
                "positive reference temperature and explicit finite bound-phase datum required",
            ));
        }
        let requirements = [
            ScalarPropertyRequirement::try_new(
                LANGMUIR_CAPACITY_PROPERTY,
                LANGMUIR_CAPACITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?,
            ScalarPropertyRequirement::try_new(
                LANGMUIR_REFERENCE_AFFINITY_PROPERTY,
                LANGMUIR_AFFINITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?,
            ScalarPropertyRequirement::try_new(
                LANGMUIR_DESORPTION_RATE_PROPERTY,
                LANGMUIR_RATE_DIMS,
                ScalarAdmissibility::NonNegative,
            )?,
        ];
        let pressure = source.point.axes().get("pressure").copied().ok_or(invalid(
            "adsorption source needs canonical typed gas pressure",
        ))?;
        let point = source_point(&source.point, reference, pressure)?;
        let resolved = resolve_material_state_point(
            &source.material,
            &point,
            &requirements,
            source.selection.clone(),
        )?;
        for name in [
            LANGMUIR_CAPACITY_PROPERTY,
            LANGMUIR_REFERENCE_AFFINITY_PROPERTY,
            LANGMUIR_DESORPTION_RATE_PROPERTY,
        ] {
            let property = resolved
                .property(name)
                .ok_or(invalid("missing selected adsorption source"))?;
            if !property
                .has_only_owned_and_fixed_external_conditions("temperature", temp_quantity())
                || property
                    .source_axes()
                    .iter()
                    .any(|(axis, _)| axis != "temperature")
                || property
                    .answer()
                    .evidence
                    .model
                    .validity
                    .bounds()
                    .keys()
                    .any(|axis| {
                        let quantity = property.answer().receipt.axis_quantities.get(axis);
                        axis == "T"
                            || axis == "p"
                            || axis == "P"
                            || (axis != "temperature"
                                && quantity.is_some_and(|q| q.dims() == Temperature::DIMS))
                            || (axis != "pressure"
                                && quantity.is_some_and(|q| q.dims() == fs_qty::Pressure::DIMS))
                    })
                || (name != LANGMUIR_DESORPTION_RATE_PROPERTY && !property.source_axes().is_empty())
            {
                return Err(invalid(
                    "capacity/reference affinity must be constant and kinetic temperature must be canonical; no stale aliases",
                ));
            }
        }
        let capacity = resolved
            .property(LANGMUIR_CAPACITY_PROPERTY)
            .unwrap()
            .value_si();
        let affinity = resolved
            .property(LANGMUIR_REFERENCE_AFFINITY_PROPERTY)
            .unwrap()
            .value_si();
        let dry_knots = calorics.dry_skeleton.curve().knots();
        let bound_knots = calorics.bound_adsorbate.curve().knots();
        if dry_knots
            .iter()
            .chain(bound_knots)
            .any(|k| k.liquid_mass_fraction != 0.0)
        {
            return Err(invalid(
                "this rigid adsorbent requires solid single-phase caloric charts",
            ));
        }
        let range = [
            dry_knots[0].temperature_k.max(bound_knots[0].temperature_k),
            dry_knots[dry_knots.len() - 1]
                .temperature_k
                .min(bound_knots[bound_knots.len() - 1].temperature_k),
        ];
        if range[0] >= range[1] || reference < range[0] || reference > range[1] {
            return Err(invalid(
                "dry and bound caloric charts need an overlapping reference-temperature span",
            ));
        }
        let gas = &catalog.thermodynamics().components()[index];
        for t in range {
            gas.evaluate(Temperature::new(t)).map_err(gas_failure)?;
        }
        let gas_h = gas
            .evaluate(Temperature::new(reference))
            .map_err(gas_failure)?
            .properties()
            .h()
            .value();
        let bound_h = calorics
            .bound_adsorbate
            .state_at_temperature(reference)?
            .specific_enthalpy_j_kg();
        let heat = finite(
            gas_h - (bound_h + calorics.bound_enthalpy_offset_j_kg) * gas.molar_mass().value(),
        )?;
        Ok(Self {
            species,
            catalog: catalog.clone(),
            gas_index: index,
            source,
            resolved,
            calorics,
            capacity_mol_kg_dry: capacity,
            reference_affinity_pa_inverse: affinity,
            reference_heat_j_mol: heat,
            temperature_range_k: range,
        })
    }
    pub fn species(&self) -> &SpeciesId {
        &self.species
    }
    pub const fn capacity_mol_kg_dry(&self) -> f64 {
        self.capacity_mol_kg_dry
    }
    pub const fn calorics(&self) -> &AdsorptionCaloricBinding {
        &self.calorics
    }
    pub const fn temperature_range_k(&self) -> [f64; 2] {
        self.temperature_range_k
    }
    fn gas_source(&self) -> &Nasa9StandardStateModelV1 {
        &self.catalog.thermodynamics().components()[self.gas_index]
    }
    fn molar_mass(&self) -> f64 {
        self.gas_source().molar_mass().value()
    }
    fn admit_catalog(&self, catalog: &ThermochemicalGasModel) -> Result<()> {
        let i = catalog
            .thermodynamics()
            .components()
            .iter()
            .position(|s| s.species() == &self.species)
            .ok_or(invalid("retained catalog is missing the adsorbate"))?;
        if &catalog.thermodynamics().components()[i] != self.gas_source()
            || catalog.transport_materials()[i]
                != self.catalog.transport_materials()[self.gas_index]
            || catalog
                .source_packs()
                .iter()
                .find(|r| r.species == self.species)
                != self
                    .catalog
                    .source_packs()
                    .iter()
                    .find(|r| r.species == self.species)
        {
            return Err(invalid(
                "adsorbate caloric, transport or source-pack identity changed",
            ));
        }
        Ok(())
    }
    fn gas_dh(&self, a: f64, b: f64) -> Result<f64> {
        Ok(self
            .gas_source()
            .molar_enthalpy_change(Temperature::new(a), Temperature::new(b))
            .map_err(gas_failure)?
            .value())
    }
    fn gas_ds(&self, a: f64, b: f64) -> Result<f64> {
        Ok(self
            .gas_source()
            .molar_entropy_change(Temperature::new(a), Temperature::new(b))
            .map_err(gas_failure)?
            .value())
    }
    fn bound_ds(&self, a: f64, b: f64) -> Result<f64> {
        Ok(self
            .calorics
            .bound_adsorbate
            .specific_entropy_change(a, b)?
            * self.molar_mass())
    }
    /// Temperature-dependent affinity from the SAME caloric and reference data.
    pub fn properties(
        &self,
        temperature_k: f64,
        total_pressure_pa: f64,
    ) -> Result<LangmuirProperties> {
        if !positive(total_pressure_pa)
            || temperature_k < self.temperature_range_k[0]
            || temperature_k > self.temperature_range_k[1]
            || !positive(temperature_k)
        {
            return Err(AdsorptionError::ThermalRange);
        }
        let t0 = self.source.reference_temperature_k;
        let dh = self.gas_dh(t0, temperature_k)?
            - self
                .calorics
                .bound_adsorbate
                .specific_enthalpy_change(t0, temperature_k)?
                * self.molar_mass();
        let heat = finite(self.reference_heat_j_mol + dh)?;
        // Stable Gibbs/van't Hoff relation; constant enthalpy datums cancel first.
        let log_ratio = (self.reference_heat_j_mol / R)
            * ((t0 - temperature_k) / t0 / temperature_k)
            + dh / (R * temperature_k)
            + (self.bound_ds(t0, temperature_k)? - self.gas_ds(t0, temperature_k)?) / R;
        let affinity = if temperature_k == t0 {
            self.reference_affinity_pa_inverse
        } else {
            det::exp(det::ln(self.reference_affinity_pa_inverse) + log_ratio)
        };
        if !positive(affinity) {
            return Err(invalid("Langmuir affinity vanished or overflowed"));
        }
        let point = source_point(&self.source.point, temperature_k, total_pressure_pa)?;
        let resolved = self
            .resolved
            .at_query_point(&self.source.material, &point)?;
        let rate = resolved
            .property(LANGMUIR_DESORPTION_RATE_PROPERTY)
            .unwrap()
            .value_si();
        Ok(LangmuirProperties {
            temperature_k,
            affinity_pa_inverse: affinity,
            desorption_rate_s_inverse: rate,
            adsorption_heat_j_mol: heat,
            source: resolved,
        })
    }
}

/// Dry skeleton mass is fixed. Bound adsorbate is a separate, finite inventory;
/// its sensible and binding storage are both part of the finite solid energy.
#[derive(Clone, Debug, PartialEq)]
pub struct AdsorptionSolid {
    model: Arc<LangmuirAdsorbent>,
    dry_mass_kg: f64,
    adsorbed_moles: f64,
    temperature_k: f64,
}
impl AdsorptionSolid {
    pub fn try_new(
        model: Arc<LangmuirAdsorbent>,
        dry_mass_kg: f64,
        adsorbed_moles: f64,
        temperature_k: f64,
    ) -> Result<Self> {
        let capacity = dry_mass_kg * model.capacity_mol_kg_dry;
        if !positive(dry_mass_kg)
            || !positive(capacity)
            || !nonnegative(adsorbed_moles)
            || adsorbed_moles > capacity
            || !positive(temperature_k)
        {
            return Err(invalid(
                "finite dry mass, representable capacity and loading within the site inventory required",
            ));
        }
        let coverage = adsorbed_moles / capacity;
        let bound_mass = adsorbed_moles * model.molar_mass();
        if (adsorbed_moles > 0.0 && (coverage == 0.0 || !positive(bound_mass)))
            || (adsorbed_moles < capacity && coverage == 1.0)
        {
            return Err(invalid(
                "nonzero occupied/vacant sites and bound mass must remain representable",
            ));
        }
        model
            .calorics
            .dry_skeleton
            .state_at_temperature(temperature_k)?;
        model
            .calorics
            .bound_adsorbate
            .state_at_temperature(temperature_k)?;
        Ok(Self {
            model,
            dry_mass_kg,
            adsorbed_moles,
            temperature_k,
        })
    }
    pub fn model(&self) -> &LangmuirAdsorbent {
        &self.model
    }
    pub const fn dry_mass_kg(&self) -> f64 {
        self.dry_mass_kg
    }
    pub const fn adsorbed_moles(&self) -> f64 {
        self.adsorbed_moles
    }
    pub fn adsorbed_mass_kg(&self) -> f64 {
        self.adsorbed_moles * self.model.molar_mass()
    }
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub fn capacity_moles(&self) -> f64 {
        self.dry_mass_kg * self.model.capacity_mol_kg_dry
    }
    pub fn coverage(&self) -> f64 {
        self.adsorbed_moles / self.capacity_moles()
    }
}

/// Physical freezing error and numerical balance controls remain distinct.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdsorptionSolve {
    pub gas: GasSpeciesBoundarySolve,
    pub maximum_gas_gross_fraction: f64,
    /// (gross adsorption + gross desorption)/total site moles, not net loading.
    pub maximum_site_turnover: f64,
    pub maximum_gas_surface_temperature_difference_k: f64,
    pub maximum_temperature_change_k: f64,
    pub equilibrium_coverage_tolerance: f64,
    pub entropy_tolerance_j_k: f64,
    pub maximum_iterations: usize,
}
impl AdsorptionSolve {
    fn validate(self, duration: f64) -> Result<()> {
        self.gas.gas.validate()?;
        if !nonnegative(duration)
            || !positive(self.gas.mass_tolerance_kg)
            || !positive(0.5 * self.gas.energy_tolerance_j)
            || !positive(self.maximum_gas_gross_fraction)
            || self.maximum_gas_gross_fraction >= 1.0
            || !positive(self.maximum_site_turnover)
            || !nonnegative(self.maximum_gas_surface_temperature_difference_k)
            || !nonnegative(self.maximum_temperature_change_k)
            || !nonnegative(self.equilibrium_coverage_tolerance)
            || self.equilibrium_coverage_tolerance >= 1.0
            || !positive(self.entropy_tolerance_j_k)
            || self.maximum_iterations == 0
            || self.maximum_iterations > 256
        {
            return Err(invalid(
                "finite physical interval, bounded inventory/thermal controls and positive numerical tolerances required",
            ));
        }
        Ok(())
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct AdsorptionStep {
    pub duration_s: f64,
    pub coverage: [f64; 2],
    pub equilibrium_coverage: [f64; 2],
    pub gas_temperature_k: [f64; 2],
    pub solid_temperature_k: [f64; 2],
    pub adsorbed_moles: [f64; 2],
    pub gross_adsorption_mol: f64,
    pub gross_desorption_mol: f64,
    pub incoming_gas_enthalpy_j_mol: f64,
    pub outgoing_gas_enthalpy_j_mol: f64,
    pub energy_into_gas_j: f64,
    pub dry_sensible_energy_change_j: f64,
    pub bound_sensible_energy_change_j: f64,
    pub solid_energy_change_j: f64,
    pub gas_energy_change_j: f64,
    pub mass_residual_kg: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_k: f64,
    pub source: [LangmuirProperties; 2],
    pub iterations: usize,
}
#[derive(Clone, Debug, PartialEq)]
pub struct AdsorptionCandidate<'catalog> {
    solid: AdsorptionSolid,
    gas: GasSpeciesBoundaryCandidate<'catalog>,
    step: AdsorptionStep,
}
impl<'catalog> AdsorptionCandidate<'catalog> {
    pub const fn solid(&self) -> &AdsorptionSolid {
        &self.solid
    }
    pub fn gas(&self) -> &ClosedGasParcel {
        self.gas.gas()
    }
    pub const fn gas_exchange(&self) -> &GasSpeciesBoundaryCandidate<'catalog> {
        &self.gas
    }
    pub const fn step(&self) -> &AdsorptionStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        AdsorptionSolid,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        AdsorptionStep,
    ) {
        let (gas, catalog, _) = self.gas.into_parts();
        (self.solid, gas, catalog, self.step)
    }
}

fn equilibrium(affinity: f64, pressure: f64) -> Result<f64> {
    let product = affinity * pressure;
    let denominator = 1.0 + product;
    if !nonnegative(product)
        || !positive(denominator)
        || (affinity > 0.0 && pressure > 0.0 && product == 0.0)
    {
        return Err(invalid("Langmuir pressure activity is unrepresentable"));
    }
    let theta = product / denominator;
    if product > 0.0 && (theta == 0.0 || theta == 1.0) {
        return Err(invalid(
            "finite Langmuir equilibrium cannot lose occupied or vacant sites",
        ));
    }
    Ok(theta)
}

/// Returns 1-exp(-u) and 1-(1-exp(-u))/u without small-u cancellation.
fn relaxation(u: f64) -> Result<(f64, f64)> {
    if !nonnegative(u) {
        return Err(invalid("Langmuir kinetic interval is unrepresentable"));
    }
    if u == 0.0 {
        return Ok((0.0, 0.0));
    }
    if u < 0.125 {
        let mut term = u;
        let mut f = term;
        let mut g = 0.5 * u;
        // g = u/2! - u²/3! + ...; fixed finite work and deterministic arithmetic.
        let mut gt = g;
        for k in 2..=18 {
            term *= -u / f64::from(k);
            f += term;
            gt *= -u / f64::from(k + 1);
            g += gt;
        }
        Ok((f, g))
    } else {
        let f = 1.0 - det::exp(-u);
        Ok((f, 1.0 - f / u))
    }
}
fn integrated_fraction(
    initial: f64,
    target: f64,
    duration: f64,
    rate: f64,
    relaxation: f64,
    remainder: f64,
) -> f64 {
    if initial >= target {
        target * duration + (initial - target) * (relaxation / rate)
    } else {
        initial * duration + (target - initial) * (duration * remainder)
    }
}
#[derive(Clone, Copy, Debug)]
struct Kinetics {
    final_moles: f64,
    adsorption: f64,
    desorption: f64,
    equilibrium: f64,
}
fn integrate_sites(
    initial: f64,
    capacity: f64,
    affinity: f64,
    desorption_rate: f64,
    pressure: f64,
    duration: f64,
) -> Result<Kinetics> {
    let theta = initial / capacity;
    let target = equilibrium(affinity, pressure)?;
    if duration == 0.0 || desorption_rate == 0.0 {
        return Ok(Kinetics {
            final_moles: initial,
            adsorption: 0.0,
            desorption: 0.0,
            equilibrium: target,
        });
    }
    let adsorption_rate = desorption_rate * (affinity * pressure);
    if pressure > 0.0 && adsorption_rate == 0.0 {
        return Err(invalid("nonzero adsorption rate underflow"));
    }
    let rate = finite(adsorption_rate + desorption_rate)?;
    if !positive(rate) {
        return Err(invalid("nonzero adsorption kinetics lost their rate"));
    }
    let u = rate * duration;
    if u == 0.0 {
        return Err(invalid("nonzero kinetic interval underflow"));
    }
    let (relax, remainder) = relaxation(u)?;
    let delta = capacity * ((target - theta) * relax);
    let next = initial + delta;
    if !nonnegative(next)
        || next > capacity
        || (target != theta && (delta == 0.0 || next == initial))
        || (initial > 0.0 && next == 0.0)
        || (initial < capacity && next == capacity)
    {
        return Err(invalid(
            "nonzero Langmuir loading change vanished or left the finite sites",
        ));
    }
    let occupied = integrated_fraction(theta, target, duration, rate, relax, remainder);
    let vacant = integrated_fraction(1.0 - theta, 1.0 - target, duration, rate, relax, remainder);
    if ((theta > 0.0 || target > 0.0) && occupied == 0.0)
        || ((theta < 1.0 || target < 1.0) && vacant == 0.0)
    {
        return Err(invalid(
            "nonzero integrated occupied/vacant site time underflow",
        ));
    }
    let adsorption = capacity * (adsorption_rate * vacant);
    let desorption = capacity * (desorption_rate * occupied);
    if !nonnegative(adsorption)
        || !nonnegative(desorption)
        || (adsorption_rate > 0.0 && vacant > 0.0 && adsorption == 0.0)
        || (occupied > 0.0 && desorption == 0.0)
    {
        return Err(invalid("gross Langmuir streams are unrepresentable"));
    }
    Ok(Kinetics {
        final_moles: next,
        adsorption,
        desorption,
        equilibrium: target,
    })
}

fn sensible(solid: &AdsorptionSolid, final_bound_mass: f64, temperature: f64) -> Result<[f64; 2]> {
    let charts = &solid.model.calorics;
    Ok([
        solid.dry_mass_kg
            * charts
                .dry_skeleton
                .specific_enthalpy_change(solid.temperature_k, temperature)?,
        final_bound_mass
            * charts
                .bound_adsorbate
                .specific_enthalpy_change(solid.temperature_k, temperature)?,
    ])
}
fn thermal_endpoint(
    solid: &AdsorptionSolid,
    final_bound_mass: f64,
    heat: f64,
    controls: AdsorptionSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, [f64; 2], usize)> {
    if heat == 0.0 {
        return Ok((solid.temperature_k, [0.0; 2], 0));
    }
    let [mut lo, mut hi] = solid.model.temperature_range_k;
    let a = sensible(solid, final_bound_mass, lo)?;
    let b = sensible(solid, final_bound_mass, hi)?;
    if heat < a[0] + a[1] || heat > b[0] + b[1] {
        return Err(AdsorptionError::ThermalRange);
    }
    let mut temperature = solid.temperature_k;
    let mut residual = f64::INFINITY;
    for iteration in 1..=controls.maximum_iterations {
        poll(cancelled)?;
        let change = sensible(solid, final_bound_mass, temperature)?;
        residual = finite((change[0] + change[1]) - heat)?;
        // A nonzero supplied heat must never be admitted as exact stasis,
        // even when it is below the caller's absolute residual allowance.
        if residual.abs() <= 0.5 * controls.gas.energy_tolerance_j
            && temperature != solid.temperature_k
        {
            return Ok((temperature, change, iteration));
        }
        if residual > 0.0 {
            hi = temperature;
        } else {
            lo = temperature;
        }
        if lo == hi {
            break;
        }
        // The existing chart owns the exact local chord slope, including its
        // knot-side convention. Newton is exact within a shared linear span;
        // the original bracket and bisection remain the convergence safeguard.
        let charts = &solid.model.calorics;
        let capacity = solid.dry_mass_kg
            * charts
                .dry_skeleton
                .secant_heat_capacity(temperature, temperature)?
            + final_bound_mass
                * charts
                    .bound_adsorbate
                    .secant_heat_capacity(temperature, temperature)?;
        if !positive(capacity) {
            return Err(invalid("finite positive combined caloric slope required"));
        }
        let newton = temperature - residual / capacity;
        let trial = if newton.is_finite() && newton >= lo && newton <= hi && newton != temperature {
            newton
        } else {
            0.5 * lo + 0.5 * hi
        };
        if trial == temperature {
            break;
        }
        temperature = trial;
    }
    Err(AdsorptionError::NotConverged {
        residual_j: residual,
    })
}

fn log_ratio(a: f64, b: f64) -> f64 {
    if a == b {
        return 0.0;
    }
    let (lo, hi, sign) = if a > b { (b, a, 1.0) } else { (a, b, -1.0) };
    let delta = (hi - lo) / lo;
    let rounded = 1.0 + delta;
    sign * if rounded == 1.0 {
        delta
    } else if rounded.is_finite() {
        det::ln(rounded) * (delta / (rounded - 1.0))
    } else {
        det::ln(hi) - det::ln(lo)
    }
}
fn mixing_remainder(a: f64, b: f64) -> f64 {
    if a == b {
        return 0.0;
    }
    let z = (b - a) / a;
    if z.abs() <= 0.125 {
        let mut polynomial = 1.0 / (22.0 * 21.0);
        for k in (2..22).rev() {
            polynomial = 1.0 / f64::from(k * (k - 1)) - z * polynomial;
        }
        ((a * z) * z) * polynomial
    } else {
        b * log_ratio(b, a) - (b - a)
    }
}
fn site_entropy_change(initial: f64, next: f64, capacity: f64) -> f64 {
    if initial == next {
        return 0.0;
    }
    if initial > 0.0 && initial < capacity && next > 0.0 && next < capacity {
        -R * ((next - initial) * log_ratio(initial, capacity - initial)
            + mixing_remainder(initial, next)
            + mixing_remainder(capacity - initial, capacity - next))
    } else {
        let xlogx = |n: f64| {
            if n == 0.0 {
                0.0
            } else {
                n * det::ln(n / capacity)
            }
        };
        -R * ((xlogx(next) + xlogx(capacity - next)) - (xlogx(initial) + xlogx(capacity - initial)))
    }
}
// Use arbitrary 1 Pa reference explicitly; the same reference enters ln(b*p_ref).
// Original species entropy constants cancel by conserved gas+site inventory.
fn combined_entropy(
    solid: &AdsorptionSolid,
    next: &AdsorptionSolid,
    gas: &ClosedGasParcel,
    next_gas: &ClosedGasParcel,
    catalog: &ThermochemicalGasModel,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let law = &solid.model;
    let t0 = law.source.reference_temperature_k;
    let dn = next.adsorbed_moles - solid.adsorbed_moles;
    let dry = &law.calorics.dry_skeleton;
    let mut entropy =
        solid.dry_mass_kg * dry.specific_entropy_change(solid.temperature_k, next.temperature_k)?;
    entropy += next.adsorbed_moles * law.bound_ds(solid.temperature_k, next.temperature_k)?
        + dn * law.bound_ds(t0, solid.temperature_k)?
        + dn * (-law.reference_heat_j_mol / t0 + R * det::ln(law.reference_affinity_pa_inverse))
        + site_entropy_change(
            solid.adsorbed_moles,
            next.adsorbed_moles,
            solid.capacity_moles(),
        );
    let (tg0, tg1) = (
        gas.state().state().temperature,
        next_gas.state().state().temperature,
    );
    for source in catalog.thermodynamics().components() {
        poll(cancelled)?;
        let n0 = gas.species_mass_kg(source.species())? / source.molar_mass().value();
        let n1 = next_gas.species_mass_kg(source.species())? / source.molar_mass().value();
        if n0 == 0.0 && n1 == 0.0 {
            continue;
        }
        let delta = if source.species() == &law.species {
            -dn
        } else {
            0.0
        };
        let ds = |a, b| {
            source
                .molar_entropy_change(Temperature::new(a), Temperature::new(b))
                .map(|s| s.value())
                .map_err(gas_failure)
        };
        entropy += n1 * ds(tg0, tg1)? + delta * ds(t0, tg0)?;
        let log_p =
            |n: f64, t: f64| det::ln(n) + det::ln(R) + det::ln(t) - det::ln(gas.volume_m3());
        entropy -= R * if n0 == 0.0 {
            n1 * log_p(n1, tg1)
        } else if n1 == 0.0 {
            -n0 * log_p(n0, tg0)
        } else {
            delta * log_p(n0, tg0) + n1 * (log_ratio(n1, n0) + log_ratio(tg1, tg0))
        };
    }
    finite(entropy)
}

/// One conservative finite gas/adsorbent proposal. Inputs never mutate.
/// The frozen rate law is deliberately admitted only near equal gas/surface T.
pub fn advance_candidate<'catalog>(
    solid: &AdsorptionSolid,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: AdsorptionSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AdsorptionCandidate<'catalog>> {
    controls.validate(duration_s)?;
    poll(&mut cancelled)?;
    let law = &solid.model;
    law.admit_catalog(catalog)?;
    let tg = gas.state().state().temperature;
    let ts = solid.temperature_k;
    if (tg - ts).abs() > controls.maximum_gas_surface_temperature_difference_k {
        return Err(invalid(
            "frozen Langmuir rates require the declared near-isothermal gas/surface regime",
        ));
    }
    let before = law.properties(ts, gas.state().state().pressure)?;
    let pressure = gas.species_partial_pressure_pa(&law.species)?;
    let site = integrate_sites(
        solid.adsorbed_moles,
        solid.capacity_moles(),
        before.affinity_pa_inverse,
        before.desorption_rate_s_inverse,
        pressure,
        duration_s,
    )?;
    let mass = law.molar_mass();
    let initial_gas_mass = gas.species_mass_kg(&law.species)?;
    let gross_fraction = if site.adsorption == 0.0 {
        0.0
    } else {
        site.adsorption * mass / initial_gas_mass
    };
    let turnover = (site.adsorption + site.desorption) / solid.capacity_moles();
    if !nonnegative(gross_fraction)
        || !nonnegative(turnover)
        || gross_fraction > controls.maximum_gas_gross_fraction
        || turnover > controls.maximum_site_turnover
    {
        return Err(AdsorptionError::InventoryLimit {
            gross_gas_fraction: gross_fraction,
            site_turnover: turnover,
        });
    }
    let dn = site.final_moles - solid.adsorbed_moles;
    let delta_mass = finite(dn * mass)?;
    if dn != 0.0 && delta_mass == 0.0 {
        return Err(invalid("adsorbed species mass change underflow"));
    }
    let site_residual = ((site.adsorption - site.desorption) - dn) * mass;
    if !site_residual.is_finite() || site_residual.abs() > controls.gas.mass_tolerance_kg {
        return Err(invalid(
            "integrated gross streams failed the original site inventory balance",
        ));
    }
    let incoming_h = catalog.species_specific_enthalpy_j_kg(&law.species, tg)? * mass;
    let outgoing_h = catalog.species_specific_enthalpy_j_kg(&law.species, ts)? * mass;
    let stream_difference = -law.gas_dh(ts, tg)?;
    let relative_energy = finite(site.desorption * stream_difference)?;
    let gas_energy = finite(-dn * incoming_h + relative_energy)?;
    let thermal_heat = finite(
        dn * (before.adsorption_heat_j_mol - stream_difference) - relative_energy,
    )?;
    let (temperature, sensible_changes, iterations) = thermal_endpoint(
        solid,
        site.final_moles * mass,
        thermal_heat,
        controls,
        &mut cancelled,
    )?;
    if (temperature - ts).abs() > controls.maximum_temperature_change_k {
        return Err(invalid(
            "solid temperature change exceeds the frozen-kinetics allowance",
        ));
    }
    let next = AdsorptionSolid::try_new(
        Arc::clone(&solid.model),
        solid.dry_mass_kg,
        site.final_moles,
        temperature,
    )?;
    poll(&mut cancelled)?;
    // Keep the source's small nonisothermal desorption heat separate from the
    // original gas partial-enthalpy datum throughout the finite gas inversion.
    let transfer = SpeciesBoundaryRelativeTransfer {
        species: law.species.clone(),
        mass_into_gas_kg: -delta_mass,
        relative_energy_into_gas_j: relative_energy,
    };
    // Reserve half of the common balance budget for each finite energy owner.
    let gas_controls = GasSpeciesBoundarySolve {
        energy_tolerance_j: 0.5 * controls.gas.energy_tolerance_j,
        ..controls.gas
    };
    let candidate =
        gas.exchange_species_boundary_relative(catalog, &[transfer], gas_controls, &mut cancelled)?;
    let gas_next = candidate.gas();
    let tg_next = gas_next.state().state().temperature;
    if (tg_next - tg).abs() > controls.maximum_temperature_change_k
        || (tg_next - temperature).abs() > controls.maximum_gas_surface_temperature_difference_k
    {
        return Err(invalid(
            "final gas/surface temperatures exceed the frozen-kinetics support",
        ));
    }
    let after = law.properties(temperature, gas_next.state().state().pressure)?;
    let target = equilibrium(
        after.affinity_pa_inverse,
        gas_next.species_partial_pressure_pa(&law.species)?,
    )?;
    if (dn > 0.0 && next.coverage() - target > controls.equilibrium_coverage_tolerance)
        || (dn < 0.0 && target - next.coverage() > controls.equilibrium_coverage_tolerance)
    {
        return Err(AdsorptionError::EquilibriumCrossed);
    }
    let entropy = combined_entropy(solid, &next, gas, gas_next, catalog, &mut cancelled)?;
    if entropy < -controls.entropy_tolerance_j_k {
        return Err(AdsorptionError::EntropyDecreased {
            change_j_k: entropy,
        });
    }
    // Cancel conserved thermochemical reference energies BEFORE rounding.
    let residual = finite(
        (sensible_changes[0] + sensible_changes[1] - thermal_heat)
            + candidate.step().energy_residual_j,
    )?;
    if residual.abs() > controls.gas.energy_tolerance_j {
        return Err(invalid(
            "combined gas/adsorbent energy balance exceeds tolerance",
        ));
    }
    let bound_h0 = law
        .calorics
        .bound_adsorbate
        .state_at_temperature(ts)?
        .specific_enthalpy_j_kg()
        + law.calorics.bound_enthalpy_offset_j_kg;
    let solid_energy = finite(sensible_changes[0] + sensible_changes[1] + delta_mass * bound_h0)?;
    let mass_residual = finite(candidate.step().mass_residual_kg + site_residual)?;
    if mass_residual.abs() > controls.gas.mass_tolerance_kg {
        return Err(invalid(
            "combined gas/site species inventory residual exceeds tolerance",
        ));
    }
    let step = AdsorptionStep {
        duration_s,
        coverage: [solid.coverage(), next.coverage()],
        equilibrium_coverage: [site.equilibrium, target],
        gas_temperature_k: [tg, tg_next],
        solid_temperature_k: [ts, temperature],
        adsorbed_moles: [solid.adsorbed_moles, next.adsorbed_moles],
        gross_adsorption_mol: site.adsorption,
        gross_desorption_mol: site.desorption,
        incoming_gas_enthalpy_j_mol: incoming_h,
        outgoing_gas_enthalpy_j_mol: outgoing_h,
        energy_into_gas_j: gas_energy,
        dry_sensible_energy_change_j: sensible_changes[0],
        bound_sensible_energy_change_j: sensible_changes[1],
        solid_energy_change_j: solid_energy,
        gas_energy_change_j: candidate.step().stored_energy_change_j,
        mass_residual_kg: mass_residual,
        energy_residual_j: residual,
        entropy_change_j_k: entropy,
        source: [before, after],
        iterations,
    };
    poll(&mut cancelled)?;
    Ok(AdsorptionCandidate {
        solid: next,
        gas: candidate,
        step,
    })
}

#[cfg(test)]
mod tests;
