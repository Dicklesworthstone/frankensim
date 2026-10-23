//! Linear electric and magnetic storage derived from material data and geometry.
//!
//! These are the uniform-field constitutive reductions C = epsilon A / d and
//! L = mu N^2 A / l. They lower into the existing circuit owner; they do not
//! integrate a circuit or guess a material from its name. Absolute SI material
//! coefficients are required, so relative permittivity/permeability cannot be
//! mistaken for their dimensioned counterparts.
//!
//! The model is isotropic, linear, nondispersive and lossless. Fringing, leakage,
//! saturation, hysteresis, dielectric loss and geometry evolution require other
//! laws. Temperature-only source curves may be resolved at a fixed temperature;
//! frequency/field-dependent curves cannot be frozen into linear storage.

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
    name: &'static str,
    dims: Dims,
) -> Result<f64, FieldStorageError> {
    let property = material.property(name).ok_or(FieldStorageError::Invalid(
        "missing or ambiguous linear field-storage property",
    ))?;
    if property.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !property.is_constant_at_fixed_temperature()
    {
        return Err(FieldStorageError::Invalid(
            "linear field storage needs an absolute SI scalar constant or temperature-only coefficient",
        ));
    }
    // A time-domain RLC reduction has no frequency-domain admission or field
    // amplitude tracker. Even a constant value with such a restricted support
    // must not turn a point measurement into an unrestricted storage law.
    if property.angular_frequency_band_rad_s()?.is_some() {
        return Err(FieldStorageError::Invalid(
            "frequency-bounded field storage needs a dispersion-aware circuit model",
        ));
    }
    for (axis, _) in &property.answer().receipt.query_point {
        let quantity = property.answer().receipt.axis_quantities.get(axis);
        if property
            .answer()
            .evidence
            .model
            .validity
            .bound(axis)
            .is_some()
            && axis != "T"
            && quantity.is_none()
        {
            return Err(FieldStorageError::Invalid(
                "bounded non-temperature storage coordinates need explicit quantity descriptors",
            ));
        }
        if quantity.is_some_and(|q| q.dims() == fs_qty::Frequency::DIMS)
            && property
                .answer()
                .evidence
                .model
                .validity
                .bound(axis)
                .is_some()
        {
            return Err(FieldStorageError::Invalid(
                "rate-bounded field storage needs an operating-frequency model",
            ));
        }
        let field_axis = quantity.is_some_and(|q| {
            [
                Dims([1, 1, -3, 0, -1, 0]), // electric field [V/m]
                Dims([-2, 0, 1, 0, 1, 0]),  // electric displacement [C/m2]
                Dims([0, 1, -2, 0, -1, 0]), // magnetic flux density [T]
                Dims([-1, 0, 0, 0, 1, 0]),  // magnetic field strength [A/m]
            ]
            .contains(&q.dims())
        }) || quantity.is_none() && ["E", "D", "B", "H"].contains(&axis.as_str());
        if field_axis
            && property
                .answer()
                .evidence
                .model
                .validity
                .bound(axis)
                .is_some()
        {
            return Err(FieldStorageError::Invalid(
                "field-bounded source needs an operating-field admission model",
            ));
        }
    }
    positive(
        property.value_si(),
        "field-storage coefficient must be finite and positive",
    )
}

/// Absolute scalar permeability for a linear isotropic lossless field owner.
/// Shares the exact source/model admission used by the circuit inductance:
/// temperature may be fixed at the resolved state, while frequency/field
/// dependence and operating-field restrictions require a different model.
/// The original material bundle remains owned by the caller.
pub fn linear_magnetic_permeability(
    material: &ResolvedMaterialStatePoint,
) -> Result<f64, FieldStorageError> {
    coefficient(
        material,
        MAGNETIC_PERMEABILITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS,
    )
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
        let epsilon = coefficient(
            &material,
            ELECTRIC_PERMITTIVITY_PROPERTY,
            ELECTRIC_PERMITTIVITY_DIMS,
        )?;
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
