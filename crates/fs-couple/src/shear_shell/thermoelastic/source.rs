//! Thermal coefficients from the exact material bundle of each shell region.
//!
//! These roles describe instantaneous expansion, specific heat at constant
//! stress and through-thickness conductivity at the original held environment.
//! Existing Kelvin and interface losses remain independent; their source data
//! must exclude the thermal contribution represented by this memory model.

use super::super::{
    MaterialShellImage, ResolvedShellRegion, invalid, owner, poll, restrict_band,
};
use crate::acoustic_realize::AcousticRealizeError;
use crate::thin_plate::{PlateMaterialModel, plate_property};
use fs_exec::CancelGate;
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
};
use fs_material::visco::slab_ted::SlabThermoelasticBending;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec};

/// Explicit source roles for instantaneous linear thermal expansion [1/K].
/// Names address properties already present in the original material bundle.
/// Negative and zero coefficients are admitted; no coefficient is fitted.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ShellThermalExpansionProperties {
    /// One coefficient in a source explicitly selected as isotropic.
    /// The resulting engineering expansion is [alpha, alpha, 0].
    Isotropic {
        /// Original instantaneous linear-expansion property.
        linear_per_k: String,
    },
    /// Full in-plane engineering expansion in the ORIGINAL card axes.
    /// Orthotropic21 exchanges the two normal entries and retains shear sign.
    Engineering {
        /// Original axis-1 linear expansion.
        normal_1_per_k: String,
        /// Original axis-2 linear expansion.
        normal_2_per_k: String,
        /// Original engineering shear expansion, 2 alpha_12, including zero.
        shear_12_per_k: String,
    },
}

/// Thermal source roles for one shell region, in original region order.
/// All properties must be finite dimension-only SI scalars at the complete
/// held material environment. The bundle must include a positive, explicitly
/// typed absolute-temperature `T` coordinate; no reference temperature is guessed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ShellThermalProperties {
    /// Explicit isotropic or engineering expansion selection.
    pub expansion: ShellThermalExpansionProperties,
    /// Specific heat at constant stress [J/(kg K)].
    pub specific_heat: String,
    /// Conductivity normal to the original shell midsurface [W/(m K)].
    /// This is an explicit material direction, not an inferred isotropic value.
    pub normal_conductivity: String,
}

#[derive(Clone, Debug)]
pub(super) struct BoundShellThermalRegion {
    pub(super) law: SlabThermoelasticBending,
    pub(super) source_band_rad_s: (f64, f64),
}

pub(super) fn bind_regions(
    image: &MaterialShellImage,
    properties: &[ShellThermalProperties],
    band_rad_s: (f64, f64),
    gate: &CancelGate,
) -> Result<Vec<BoundShellThermalRegion>, AcousticRealizeError> {
    poll(gate)?;
    let source = image.source();
    if properties.len() != source.regions().len() || properties.is_empty()
        || !band_rad_s.0.is_finite() || !band_rad_s.1.is_finite()
        || band_rad_s.0 < 0.0 || band_rad_s.1 <= band_rad_s.0
    {
        return Err(invalid("shell thermal binding requires one selection per region and a finite frequency band"));
    }
    if !source.coatings.is_empty() || source.coating_relaxation.is_some()
        || source.substrate_relaxation.is_some() || image.substrate_relaxation().is_some()
        || image.coating_relaxation().is_some()
    {
        return Err(invalid("shell slab thermal memory requires an uncoated homogeneous substrate without other material memory"));
    }
    let mut bound = Vec::with_capacity(properties.len());
    for (region, thermal) in source.regions().iter().zip(properties) {
        poll(gate)?;
        bound.push(bind_region(region, thermal, band_rad_s)?);
    }
    poll(gate)?;
    Ok(bound)
}

fn bind_region(
    region: &ResolvedShellRegion,
    properties: &ShellThermalProperties,
    requested_band: (f64, f64),
) -> Result<BoundShellThermalRegion, AcousticRealizeError> {
    if region.constant_loss_reference_omega_rad_s().is_some() {
        return Err(invalid("shell thermal memory requires independent Kelvin viscosity, not frequency-anchored substrate loss"));
    }
    let input = &region.input().region;
    let state = &input.material;
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static,
    ));
    let temperature = state.state_coordinate("T")
        .filter(|t| t.is_finite() && *t > 0.0)
        .ok_or_else(|| invalid("shell thermal memory requires positive absolute temperature at the original state"))?;
    if state.axis_quantities().get("T") != Some(&absolute) {
        return Err(invalid("shell thermal temperature must be explicitly typed absolute and static"));
    }
    let support = region.omega_band_rad_s();
    let mut band = (requested_band.0.max(support.0), requested_band.1.min(support.1));
    let mut read = |name: &str, dims| {
        let value = plate_property(state, name, dims)?;
        restrict_band(state, name, &mut band, false)?;
        Ok::<f64, AcousticRealizeError>(value)
    };
    let mut alpha = match &properties.expansion {
        ShellThermalExpansionProperties::Isotropic { linear_per_k } => {
            if !input.model.is_isotropic() {
                return Err(invalid("isotropic thermal expansion requires an explicitly isotropic shell source model"));
            }
            let a = read(linear_per_k, INVERSE_TEMPERATURE_DIMS)?;
            [a, a, 0.0]
        }
        ShellThermalExpansionProperties::Engineering {
            normal_1_per_k, normal_2_per_k, shear_12_per_k,
        } => [
            read(normal_1_per_k, INVERSE_TEMPERATURE_DIMS)?,
            read(normal_2_per_k, INVERSE_TEMPERATURE_DIMS)?,
            read(shear_12_per_k, INVERSE_TEMPERATURE_DIMS)?,
        ],
    };
    let cp = read(&properties.specific_heat, SPECIFIC_HEAT_CAPACITY_DIMS)?;
    let conductivity = read(&properties.normal_conductivity, THERMAL_CONDUCTIVITY_DIMS)?;
    if cp <= 0.0 || conductivity <= 0.0 {
        return Err(invalid("shell thermal specific heat and normal conductivity must be positive"));
    }
    if !band.0.is_finite() || !band.1.is_finite() || band.0 < 0.0 || band.1 <= band.0 {
        return Err(invalid("shell thermal and mechanical applicability have no nonempty requested frequency interval"));
    }
    // Q is already in the original section's selected principal frame. Axes21
    // swaps x/y and reverses z; in-plane engineering shear does not change sign.
    if input.model == PlateMaterialModel::Orthotropic21 {
        alpha.swap(0, 1);
    }
    let h = region.thickness_m();
    let rho = plate_property(state, "density", Density::DIMS)?;
    let section = region.section();
    let (a, d) = section.principal_matrices();
    let a = *a.coefficients();
    let d = d.coefficients();
    // These checks admit only the roundoff from the existing homogeneous-solid
    // constructor. They do not turn a laminate or an effective sheet into a slab.
    if section.coupling_pa_m2().iter().any(|v| *v != 0.0)
        || section.mass_first_moment_kg_m() != 0.0
        || !same_coefficient(section.areal_mass_kg_m2(), rho * h)
        || !same_coefficient(section.rotary_inertia_kg(), (rho * h * h) * (h / 12.0))
        || a.iter().zip(d).any(|(&a, &d)| !same_coefficient(d, (a * h) * (h / 12.0)))
    {
        return Err(invalid("shell thermal memory requires the original symmetric homogeneous solid section"));
    }
    let law = SlabThermoelasticBending::new_coupled(
        a.map(|value| value / h), alpha, temperature, rho * cp, conductivity, h,
    ).map_err(owner)?;
    Ok(BoundShellThermalRegion { law, source_band_rad_s: band })
}

fn same_coefficient(actual: f64, expected: f64) -> bool {
    actual.is_finite() && expected.is_finite()
        && (actual == expected
            || (actual - expected).abs()
                <= 64.0 * f64::EPSILON * actual.abs().max(expected.abs()))
}

#[cfg(test)]
pub(super) mod tests;
