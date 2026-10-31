//! Linear electric and magnetic storage derived from material data and geometry.
//!
//! These are the uniform-field constitutive reductions C = epsilon A / d and
//! L = mu N^2 A / l. They lower into the existing circuit owner; they do not
//! integrate a circuit or guess a material from its name. Absolute SI material
//! coefficients are required, so relative permittivity/permeability cannot be
//! mistaken for their dimensioned counterparts.
//!
//! The uniform reductions below are isotropic, linear, nondispersive and lossless.
//! The planar module additionally admits explicitly oriented dielectric sources;
//! [`dielectric`] supplies a separate passive Debye polarization-memory model.
//! Fringing, leakage,
//! saturation, hysteresis, dielectric loss and geometry evolution require other
//! laws. Source curves and surfaces may be resolved at fixed temperature,
//! pressure, composition, moisture and concentration, with typed external roles.
//! Frequency, field, strain or stress dependence cannot be frozen into linear
//! storage. These are snapshots; evolving a coefficient requires its parameter
//! work and, where relevant, thermal entropy to be accounted by the caller.

/// Source-bound reversible permanent-magnet recoil at a fixed material state.
pub mod permanent;
/// Source-resolved reversible nonlinear isotropic magnetic storage.
pub mod nonlinear;
/// Opt-in scalar magnetic history with explicit stored energy and pinning heat.
pub mod hysteresis;
/// Source-bound reversible thermal field energy and entropy at fixed geometry.
pub mod thermal;
/// Source-backed scalar or oriented planar dielectric response.
pub mod planar;
/// Passive source-backed dielectric polarization memory and its RC realization.
pub mod dielectric;

use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Dims, QuantitySpec};

use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};

/// Absolute electric permittivity, in farads per metre.
pub const ELECTRIC_PERMITTIVITY_PROPERTY: &str = "electric_permittivity";
/// SI dimensions of absolute permittivity, not relative permittivity.
pub const ELECTRIC_PERMITTIVITY_DIMS: Dims = Dims([-3, -1, 4, 0, 2, 0]);
/// Absolute magnetic permeability, in henries per metre.
pub const MAGNETIC_PERMEABILITY_PROPERTY: &str = "magnetic_permeability";
/// SI dimensions of absolute permeability, not porous permeability.
pub const MAGNETIC_PERMEABILITY_DIMS: Dims = Dims([1, 1, -2, 0, -2, 0]);

/// Source, model-domain or arithmetic refusal; no fallback coefficient is used.
#[derive(Debug)]
pub enum FieldStorageError {
    /// Shared source selection, quantity or validity refusal.
    Material(MaterialStatePointError),
    /// Unsupported model, geometry or unrepresentable storage.
    Invalid(&'static str),
}

impl From<MaterialStatePointError> for FieldStorageError {
    fn from(value: MaterialStatePointError) -> Self {
        Self::Material(value)
    }
}
impl core::fmt::Display for FieldStorageError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material field storage: {self:?}")
    }
}
impl std::error::Error for FieldStorageError {}

fn positive(value: f64, what: &'static str) -> Result<f64, FieldStorageError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(FieldStorageError::Invalid(what))
    }
}

fn coefficient(
    material: &ResolvedMaterialStatePoint,
    name: &str,
    dims: Dims,
) -> Result<f64, FieldStorageError> {
    let property = material.property(name).ok_or(FieldStorageError::Invalid(
        "missing or ambiguous linear field-storage property",
    ))?;
    // Reuse the material-wide source-domain policy: every value axis AND every
    // applicability bound must describe a fixed physical context. In particular,
    // pressure is not solid stress and composition/moisture are not strain just
    // because their dimensions happen to agree. Frequency applicability is not
    // admitted for this nondispersive time-domain storage law, even for a scalar.
    if property.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !property.is_static_at_fixed_external_conditions()
    {
        return Err(FieldStorageError::Invalid(
            "linear field storage requires an absolute SI scalar at fixed external conditions, without frequency, field, strain or stress dependence or bounds",
        ));
    }
    positive(
        property.value_si(),
        "field-storage coefficient must be finite and positive",
    )
}

/// Absolute scalar permittivity for a linear isotropic lossless field owner.
/// Uses the same source/model admission as the uniform circuit capacitor.
/// Fixed external conditions remain attached to the original material bundle;
/// frequency, field, strain and stress dependence or bounds require other laws.
pub fn linear_electric_permittivity(
    material: &ResolvedMaterialStatePoint,
) -> Result<f64, FieldStorageError> {
    coefficient(
        material,
        ELECTRIC_PERMITTIVITY_PROPERTY,
        ELECTRIC_PERMITTIVITY_DIMS,
    )
}

/// Absolute permittivity explicitly constant across its whole source domain.
/// Retained charge owners without dielectric parameter work/thermal entropy
/// must use this admission rather than infer constancy from two snapshots.
/// Source applicability still follows the original nondispersive static law;
/// runtime refreshes must retain the claim and remain inside its validity.
pub fn constant_electric_permittivity(
    material: &ResolvedMaterialStatePoint,
) -> Result<f64, FieldStorageError> {
    constant_coefficient(material, ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS)
}

fn constant_coefficient(
    material: &ResolvedMaterialStatePoint,
    name: &'static str,
    dims: Dims,
) -> Result<f64, FieldStorageError> {
    let value = coefficient(material, name, dims)?;
    let source = material.property(name).expect("linear admission retains its source");
    if !source.source_axes().is_empty()
        || source.answer().receipt.decision != fs_matdb::EvaluationDecision::ConstantWithinValidity
    {
        return Err(FieldStorageError::Invalid(
            "live dielectric storage without parameter work requires validity-wide constant permittivity",
        ));
    }
    Ok(value)
}

/// Absolute scalar permeability for a linear isotropic lossless field owner.
/// Shares the exact source/model admission used by the circuit inductance:
/// temperature and typed pressure/composition/moisture/concentration may be
/// fixed at the resolved state. Frequency, field, strain and stress dependence
/// or bounds require a different model. The original bundle remains owned by
/// the caller; changing its coefficient is not a zero-work storage update.
pub fn linear_magnetic_permeability(
    material: &ResolvedMaterialStatePoint,
) -> Result<f64, FieldStorageError> {
    coefficient(
        material,
        MAGNETIC_PERMEABILITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS,
    )
}

/// Admit a temperature refresh without changing prepared linear magnetic storage.
///
/// The complete original card, coordinate schema, requirements and selected
/// claims stay pinned. Only one owned absolute-temperature coordinate may move;
/// all other physical coordinates remain fixed. A permeability law with actual
/// temperature dependence refuses any temperature change, even if its endpoint
/// values happen to coincide: fixed-field storage cannot pay its parameter work.
/// This is shared by spatial induction and circuit-derived winding histories.
pub fn admit_fixed_magnetic_refresh(
    prepared: &ResolvedMaterialStatePoint,
    refreshed: &ResolvedMaterialStatePoint,
) -> Result<(), FieldStorageError> {
    admit_magnetic_temperature_refresh(prepared, refreshed)?;
    let mu = linear_magnetic_permeability(refreshed)?;
    let original_mu = linear_magnetic_permeability(prepared)?;
    let property = refreshed.property(MAGNETIC_PERMEABILITY_PROPERTY).expect("admitted permeability");
    if mu.to_bits() != original_mu.to_bits()
        || prepared.query_point().iter().zip(refreshed.query_point()).any(|((axis, a), (_, b))| {
            a.to_bits() != b.to_bits() && property.source_axes().iter().any(|(name, _)| name == axis)
        })
    {
        return Err(FieldStorageError::Invalid("heated magnetic region changes prepared storage; new extraction and magnetic parameter work required"));
    }
    Ok(())
}

/// Admit a source-pinned temperature change for an owner that REASSEMBLES its
/// magnetic operator and explicitly accounts for the resulting parameter work.
/// This only checks source/context continuity, not a thermodynamic derivative
/// or energy balance; thermal owners must additionally bind a supported law.
pub fn admit_magnetic_temperature_refresh(
    prepared: &ResolvedMaterialStatePoint,
    refreshed: &ResolvedMaterialStatePoint,
) -> Result<(), FieldStorageError> {
    admit_source_temperature_refresh(prepared, refreshed)?;
    linear_magnetic_permeability(prepared)?;
    linear_magnetic_permeability(refreshed)?;
    Ok(())
}

/// Admit a temperature-only refresh of one complete original source bundle.
/// Card identity, quantity conventions, requirements and selected claims remain
/// pinned; all nonthermal coordinates stay fixed. This shared check also covers
/// electrical/thermal properties carried beside a nonlinear magnetic law.
/// It does not admit a changed constitutive law or account for parameter work;
/// the physical owner must separately enforce its storage and heat semantics.
pub fn admit_source_temperature_refresh(
    prepared: &ResolvedMaterialStatePoint,
    refreshed: &ResolvedMaterialStatePoint,
) -> Result<(), FieldStorageError> {
    use fs_qty::{Temperature, semantic::{QuantityKind, ValueForm}};
    use FieldStorageError::Invalid;
    if prepared.card_identity() != refreshed.card_identity()
        || prepared.material() != refreshed.material()
        || prepared.axis_quantities() != refreshed.axis_quantities()
        || prepared.query_point().len() != refreshed.query_point().len()
        || prepared.properties().len() != refreshed.properties().len()
        || prepared.properties().iter().zip(refreshed.properties()).any(|(a, b)| {
            a.requirement() != b.requirement() || a.source_axes() != b.source_axes()
                || a.answer().receipt.selected != b.answer().receipt.selected
        })
    {
        return Err(Invalid("magnetic refresh must retain the complete original card, requirements, coordinate schema and selected claims"));
    }
    let mut temperatures = 0;
    for ((name, value), (next_name, next_value)) in prepared.query_point().iter().zip(refreshed.query_point()) {
        if name != next_name { return Err(Invalid("magnetic source coordinates changed correspondence")); }
        let quantity = prepared.axis_quantities().get(name).copied();
        let is_temperature = quantity.and_then(QuantitySpec::semantic_type).is_some_and(|ty| {
            ty.kind() == QuantityKind::AbsoluteTemperature && ty.form() == ValueForm::Static
        }) || name == "T" && quantity.is_none_or(|q| q == QuantitySpec::dimensional(Temperature::DIMS));
        if is_temperature {
            temperatures += 1;
            if !value.is_finite() || *value <= 0.0 || !next_value.is_finite() || *next_value <= 0.0 {
                return Err(Invalid("magnetic temperature must be positive and finite"));
            }
        }
        if value.to_bits() != next_value.to_bits() {
            if !is_temperature { return Err(Invalid("only the owned magnetic temperature may change during a fixed-geometry step")); }
        }
    }
    if temperatures > 1 { return Err(Invalid("magnetic temperature coordinate is ambiguous")); }
    Ok(())
}

fn resolve(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    name: &'static str,
    dims: Dims,
) -> Result<ResolvedMaterialStatePoint, FieldStorageError> {
    let requirement =
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)?;
    Ok(resolve_material_state_point(
        card,
        point,
        &[requirement],
        selection,
    )?)
}

fn energy(coordinate: f64, storage: f64) -> Result<f64, FieldStorageError> {
    if !coordinate.is_finite() {
        return Err(FieldStorageError::Invalid(
            "storage coordinate must be finite",
        ));
    }
    let result = 0.5 * coordinate * (coordinate / storage);
    if !result.is_finite() || (coordinate != 0.0 && result == 0.0) {
        return Err(FieldStorageError::Invalid(
            "field energy is not representable",
        ));
    }
    Ok(result)
}

/// A uniform dielectric gap with its complete material snapshot retained.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedCapacitor {
    material: ResolvedMaterialStatePoint,
    area_m2: f64,
    gap_m: f64,
    capacitance_f: f64,
}
impl ResolvedCapacitor {
    /// Adopt the SAME bundle used by other physics owners; preserve all receipts.
    /// Area and separation are current physical geometry, not source properties.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        area_m2: f64,
        gap_m: f64,
    ) -> Result<Self, FieldStorageError> {
        positive(area_m2, "capacitor area must be finite and positive")?;
        positive(gap_m, "capacitor gap must be finite and positive")?;
        let epsilon = linear_electric_permittivity(&material)?;
        let capacitance_f = positive(
            epsilon * (area_m2 / gap_m),
            "capacitance must be representable and positive",
        )?;
        positive(
            capacitance_f.recip(),
            "inverse capacitance must be representable and positive",
        )?;
        Ok(Self {
            material,
            area_m2,
            gap_m,
            capacitance_f,
        })
    }
    /// Resolve the dielectric at the complete physical query, then lower geometry.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
        area_m2: f64,
        gap_m: f64,
    ) -> Result<Self, FieldStorageError> {
        Self::try_from_material(
            resolve(
                card,
                point,
                selection,
                ELECTRIC_PERMITTIVITY_PROPERTY,
                ELECTRIC_PERMITTIVITY_DIMS,
            )?,
            area_m2,
            gap_m,
        )
    }
    /// Original source bundle; no reciprocal or derived claim has been fabricated.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Current facing area [m2].
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    /// Current uniform dielectric separation [m].
    pub const fn gap_m(&self) -> f64 {
        self.gap_m
    }
    /// Derived C = epsilon A / d [F].
    pub const fn capacitance_f(&self) -> f64 {
        self.capacitance_f
    }
    /// Stored field energy q^2 / (2 C), with signed charge in coulombs.
    pub fn stored_energy_j(&self, charge_c: f64) -> Result<f64, FieldStorageError> {
        energy(charge_c, self.capacitance_f)
    }
}

/// A uniform linear magnetic path and ideal winding, without leakage or fringe.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedInductor {
    material: ResolvedMaterialStatePoint,
    area_m2: f64,
    path_length_m: f64,
    turns: u32,
    inductance_h: f64,
}
impl ResolvedInductor {
    /// Adopt a source snapshot; N is the actual winding count, not a fitted gain.
    /// The conserved circuit coordinate is flux LINKAGE N*Phi, in webers.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        area_m2: f64,
        path_length_m: f64,
        turns: u32,
    ) -> Result<Self, FieldStorageError> {
        positive(area_m2, "magnetic area must be finite and positive")?;
        positive(
            path_length_m,
            "magnetic path length must be finite and positive",
        )?;
        if turns == 0 {
            return Err(FieldStorageError::Invalid("winding needs positive turns"));
        }
        let mu = linear_magnetic_permeability(&material)?;
        let n = f64::from(turns);
        let inductance_h = positive(
            mu * (area_m2 / path_length_m) * n * n,
            "inductance must be representable and positive",
        )?;
        positive(
            inductance_h.recip(),
            "inverse inductance must be representable and positive",
        )?;
        Ok(Self {
            material,
            area_m2,
            path_length_m,
            turns,
            inductance_h,
        })
    }
    /// Resolve the magnetic material at the complete physical query.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
        area_m2: f64,
        path_length_m: f64,
        turns: u32,
    ) -> Result<Self, FieldStorageError> {
        Self::try_from_material(
            resolve(
                card,
                point,
                selection,
                MAGNETIC_PERMEABILITY_PROPERTY,
                MAGNETIC_PERMEABILITY_DIMS,
            )?,
            area_m2,
            path_length_m,
            turns,
        )
    }
    /// Complete retained source material and state.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Current uniform magnetic-path area [m2].
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    /// Current uniform path length [m].
    pub const fn path_length_m(&self) -> f64 {
        self.path_length_m
    }
    /// Actual winding count.
    pub const fn turns(&self) -> u32 {
        self.turns
    }
    /// Derived L = mu N^2 A / l [H].
    pub const fn inductance_h(&self) -> f64 {
        self.inductance_h
    }
    /// Stored field energy lambda^2 / (2 L), with flux LINKAGE in webers.
    pub fn stored_energy_j(&self, flux_linkage_wb: f64) -> Result<f64, FieldStorageError> {
        energy(flux_linkage_wb, self.inductance_h)
    }
}

