//! Explicit material-card coefficients for the original Hill48 return map.
//!
//! The reference yield radius, the work-conjugate hardening modulus and all
//! six dimensionless Hill coefficients must be independently selected. A
//! uniaxial yield observation does not define the reference radius unless its
//! direction and coefficient normalization have already been resolved by the
//! source. A J2 hardening slope, hardness or elongation cannot supply missing
//! directional data. The caller explicitly binds fixed material axes.
//! Mixed hardening additionally selects its own Prager coefficient C in
//! X=C*eps_p; the J2 uniaxial H_k convention is never substituted. Optional
//! Voce Q and b are also Hill-specific: b is per Hill work-conjugate alpha,
//! which generally differs from the J2 equivalent plastic strain.
//!
//! Original claims, conditions, units, uncertainty and usage receipts remain
//! in the complete material bundle. This adapter freezes supported external
//! conditions; it does not evolve temperature, texture, rate or damage, and
//! makes no empirical calibration claim for authored coefficient fixtures.

use core::fmt;

/// Calibration from measured planar plastic strain ratios and explicit shear data.
pub mod planar;

use fs_evidence::ModelCard;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Dims, Pressure, QuantitySpec};

use super::oriented::OrientedHill48;
use super::{Hill48, Hill48Coefficients, HillError, HillResponse, HillState};
use crate::elastic::IsotropicElastic;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, POISSON_RATIO_PROPERTY,
    ResolvedMaterialStatePoint, ResolvedScalarProperty, ScalarAdmissibility,
    ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};

/// Reference radius Y [Pa] of the declared dimensionless Hill quadratic.
/// This is not a direction-free uniaxial tensile yield observation.
pub const HILL_REFERENCE_YIELD_STRESS_PROPERTY: &str = "hill_reference_yield_stress";
/// Energetic linear hardening H [Pa] for Hill's work-conjugate accumulated flow.
/// An isotropic J2 equivalent-strain modulus is not substituted.
pub const HILL_HARDENING_MODULUS_PROPERTY: &str = "hill_hardening_modulus";
/// Energetic Prager coefficient C [Pa] in `X=C*eps_p`, with storage
/// `C*eps_p:eps_p/2`. This is not the J2 uniaxial hardening modulus H_k.
pub const HILL_KINEMATIC_MODULUS_PROPERTY: &str = "hill_kinematic_hardening_modulus";
/// Added isotropic saturation strength Q [Pa] of the Hill yield radius.
pub const HILL_VOCE_SATURATION_INCREMENT_PROPERTY: &str = "hill_voce_saturation_increment";
/// Saturation rate b per Hill work-conjugate alpha [dimensionless].
/// A J2-equivalent-strain rate is not interchangeable for an anisotropic metric.
pub const HILL_VOCE_SATURATION_RATE_PROPERTY: &str = "hill_voce_saturation_rate";
/// Dimensionless coefficients in the exact order F, G, H, L, M, N.
/// Normal stress differences use F/G/H; shear yz/zx/xy uses 2L/2M/2N.
pub const HILL_COEFFICIENT_PROPERTIES: [&str; 6] = [
    "hill_coefficient_f",
    "hill_coefficient_g",
    "hill_coefficient_h",
    "hill_coefficient_l",
    "hill_coefficient_m",
    "hill_coefficient_n",
];

/// All ten explicit scalar requirements. Valid negative individual F/G/H
/// remain possible; the native law checks the complete deviatoric quadratic.
pub fn hill48_requirements() -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    let mut requirements = vec![
        ScalarPropertyRequirement::try_new(
            YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            POISSON_RATIO_PROPERTY,
            Dims([0; 6]),
            ScalarAdmissibility::OpenInterval {
                lower: -1.0,
                upper: 0.5,
            },
        )?,
        ScalarPropertyRequirement::try_new(
            HILL_REFERENCE_YIELD_STRESS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            HILL_HARDENING_MODULUS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::NonNegative,
        )?,
    ];
    for (index, name) in HILL_COEFFICIENT_PROPERTIES.iter().enumerate() {
        requirements.push(ScalarPropertyRequirement::try_new(
            *name,
            Dims([0; 6]),
            if index < 3 {
                ScalarAdmissibility::Finite
            } else {
                ScalarAdmissibility::StrictlyPositive
            },
        )?);
    }
    Ok(requirements)
}

/// All eleven independent coefficients of mixed isotropic/Prager Hill flow.
/// An explicitly sourced zero C preserves the isotropic-hardening model.
pub fn hill48_mixed_requirements() -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError>
{
    let mut requirements = hill48_requirements()?;
    requirements.push(ScalarPropertyRequirement::try_new(
        HILL_KINEMATIC_MODULUS_PROPERTY,
        Pressure::DIMS,
        ScalarAdmissibility::NonNegative,
    )?);
    Ok(requirements)
}

/// Elasticity, Hill metric and linear hardening plus independently sourced Q/b.
/// Q=0 is an explicit measured/model selection; absent Q or b is not zero.
/// Add the separate C requirement when selecting Prager hardening as well.
pub fn hill48_voce_requirements() -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError>
{
    let mut requirements = hill48_requirements()?;
    requirements.extend(voce_requirements()?);
    Ok(requirements)
}

pub(super) fn voce_requirements() -> Result<[ScalarPropertyRequirement; 2], MaterialStatePointError>
{
    Ok([
        ScalarPropertyRequirement::try_new(
            HILL_VOCE_SATURATION_INCREMENT_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::NonNegative,
        )?,
        ScalarPropertyRequirement::try_new(
            HILL_VOCE_SATURATION_RATE_PROPERTY,
            Dims([0; 6]),
            ScalarAdmissibility::StrictlyPositive,
        )?,
    ])
}

/// Original source bundle, immutable orientation and one native constitutive law.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedHill48 {
    material: ResolvedMaterialStatePoint,
    law: OrientedHill48,
}

impl ResolvedHill48 {
    /// Bind coefficients at one original material state. Existing extra physics
    /// receipts are preserved. The strain bound must lie in (0, 0.05], and
    /// columns of `material_to_spatial` explicitly identify the material axes.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        strain_limit: f64,
        material_to_spatial: [[f64; 3]; 3],
    ) -> Result<Self, HillSourceError> {
        if !(strain_limit.is_finite() && strain_limit > 0.0 && strain_limit <= 0.05) {
            return Err(HillSourceError::Invalid(
                "small-strain bound must be in (0, 0.05]",
            ));
        }
        let youngs = coefficient(&material, YOUNG_MODULUS_PROPERTY, Pressure::DIMS)?.value_si();
        let poisson = coefficient(&material, POISSON_RATIO_PROPERTY, Dims([0; 6]))?.value_si();
        let yield_stress = coefficient(
            &material,
            HILL_REFERENCE_YIELD_STRESS_PROPERTY,
            Pressure::DIMS,
        )?
        .value_si();
        let hardening =
            coefficient(&material, HILL_HARDENING_MODULUS_PROPERTY, Pressure::DIMS)?.value_si();
        let mut values = [0.0; 6];
        for (value, name) in values.iter_mut().zip(HILL_COEFFICIENT_PROPERTIES) {
            *value = coefficient(&material, name, Dims([0; 6]))?.value_si();
        }
        let [f, g, h, l, m, n] = values;
        let local = Hill48::new(
            IsotropicElastic::new(youngs, poisson, strain_limit)?,
            yield_stress,
            hardening,
            Hill48Coefficients { f, g, h, l, m, n },
        )?;
        let law = OrientedHill48::new(local, material_to_spatial)?;
        Ok(Self { material, law })
    }

    /// Opt into the original native Prager law using the independently selected
    /// C already retained in this material bundle. Every original receipt and
    /// the immutable orientation survive. The ordinary constructor retains C=0
    /// even if an extra C receipt is present, until this explicit opt-in.
    ///
    /// Select this before retaining specimen history; this constructor does
    /// not rebind an existing accepted state or infer C from J2, tensile or
    /// isotropic-hardening data. Missing, negative or unowned-axis C refuses.
    pub fn with_linear_kinematic_hardening(mut self) -> Result<Self, HillSourceError> {
        let modulus = coefficient(
            &self.material,
            HILL_KINEMATIC_MODULUS_PROPERTY,
            Pressure::DIMS,
        )?
        .value_si();
        let local = self
            .law
            .law()
            .clone()
            .with_linear_kinematic_hardening(modulus)?;
        self.law = OrientedHill48::new(local, *self.law.material_to_spatial())?;
        Ok(self)
    }

    /// Select the independently sourced Hill Q/b already retained in this
    /// bundle. The original axes, Prager coefficient and all receipts survive.
    /// Select before loading; resolving a new external temperature alone does
    /// not transport history or account for temperature-dependent free energy.
    pub fn with_voce_isotropic_hardening(mut self) -> Result<Self, HillSourceError> {
        let saturation = coefficient(
            &self.material,
            HILL_VOCE_SATURATION_INCREMENT_PROPERTY,
            Pressure::DIMS,
        )?
        .value_si();
        let rate = coefficient(
            &self.material,
            HILL_VOCE_SATURATION_RATE_PROPERTY,
            Dims([0; 6]),
        )?
        .value_si();
        let local = self
            .law
            .law()
            .clone()
            .with_voce_isotropic_hardening(saturation, rate)?;
        self.law = OrientedHill48::new(local, *self.law.material_to_spatial())?;
        Ok(self)
    }

    /// Complete original material identity, query, source selection and receipts.
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }

    /// Native oriented return, stored energy, and physical/numerical work accounts.
    pub const fn law(&self) -> &OrientedHill48 {
        &self.law
    }

    /// Original local coefficients with their original numerical applicability.
    pub const fn local_law(&self) -> &Hill48 {
        self.law.law()
    }

    /// Spatial response with retained local-axis material history.
    pub fn response(&self, strain: &Voigt, state: &HillState) -> Result<HillResponse, HillError> {
        self.law.response(strain, state)
    }

    /// Explicitly selected reference-radius source, not a substituted yield datum.
    pub fn reference_yield_source(&self) -> &ResolvedScalarProperty {
        self.material
            .property(HILL_REFERENCE_YIELD_STRESS_PROPERTY)
            .expect("immutable Hill bundle retains its admitted reference radius")
    }

    /// Original modulus against Hill accumulated flow, retaining its usage receipt.
    pub fn hardening_source(&self) -> &ResolvedScalarProperty {
        self.material
            .property(HILL_HARDENING_MODULUS_PROPERTY)
            .expect("immutable Hill bundle retains its admitted hardening source")
    }
}

impl SmallStrainLaw for ResolvedHill48 {
    type State = HillState;
    fn initial_state(&self) -> HillState {
        self.law.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &HillState) -> Voigt {
        self.law.stress(strain, state)
    }
    fn tangent(&self, strain: &Voigt, state: &HillState) -> Tangent6 {
        self.law.tangent(strain, state)
    }
    fn update_state(&self, strain: &Voigt, state: &HillState) -> HillState {
        self.law.update_state(strain, state)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.law.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.law.card()
    }
}

/// Resolve the exact coefficients at one explicit source point and orientation.
pub fn resolve_hill48(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    strain_limit: f64,
    material_to_spatial: [[f64; 3]; 3],
) -> Result<ResolvedHill48, HillSourceError> {
    let material = resolve_material_state_point(card, point, &hill48_requirements()?, selection)?;
    ResolvedHill48::try_from_material(material, strain_limit, material_to_spatial)
}

/// Resolve all independent mixed-hardening coefficients at one fixed physical
/// point, retaining their source selection and explicit material orientation.
pub fn resolve_hill48_mixed(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    strain_limit: f64,
    material_to_spatial: [[f64; 3]; 3],
) -> Result<ResolvedHill48, HillSourceError> {
    let material =
        resolve_material_state_point(card, point, &hill48_mixed_requirements()?, selection)?;
    ResolvedHill48::try_from_material(material, strain_limit, material_to_spatial)?
        .with_linear_kinematic_hardening()
}

/// Resolve an explicitly selected Hill saturation law at one fixed source point.
pub fn resolve_hill48_voce(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    strain_limit: f64,
    material_to_spatial: [[f64; 3]; 3],
) -> Result<ResolvedHill48, HillSourceError> {
    let material =
        resolve_material_state_point(card, point, &hill48_voce_requirements()?, selection)?;
    ResolvedHill48::try_from_material(material, strain_limit, material_to_spatial)?
        .with_voce_isotropic_hardening()
}

/// Original source, elasticity, orientation or native Hill refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum HillSourceError {
    MaterialState(MaterialStatePointError),
    Elasticity(MaterialError),
    Constitutive(HillError),
    Invalid(&'static str),
}
impl fmt::Display for HillSourceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MaterialState(source) => write!(f, "Hill48 source refusal: {source}"),
            Self::Elasticity(source) => write!(f, "Hill48 elasticity refusal: {source}"),
            Self::Constitutive(source) => write!(f, "Hill48 constitutive refusal: {source}"),
            Self::Invalid(what) => write!(f, "Hill48 reduction refusal: {what}"),
        }
    }
}
impl std::error::Error for HillSourceError {}
impl From<MaterialStatePointError> for HillSourceError {
    fn from(source: MaterialStatePointError) -> Self {
        Self::MaterialState(source)
    }
}
impl From<MaterialError> for HillSourceError {
    fn from(source: MaterialError) -> Self {
        Self::Elasticity(source)
    }
}
impl From<HillError> for HillSourceError {
    fn from(source: HillError) -> Self {
        Self::Constitutive(source)
    }
}

fn coefficient<'a>(
    material: &'a ResolvedMaterialStatePoint,
    name: &'static str,
    dims: Dims,
) -> Result<&'a ResolvedScalarProperty, HillSourceError> {
    let property = material.property(name).ok_or(HillSourceError::Invalid(
        "missing or ambiguous Hill coefficient",
    ))?;
    if property.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !property.value_si().is_finite()
        || !property.has_only_fixed_external_conditions(None, false)
        || property.angular_frequency_band_rad_s()?.is_some()
    {
        return Err(HillSourceError::Invalid(
            "Hill48 needs absolute SI coefficients at fixed external conditions, without evolving mechanical, history, rate or frequency coordinates",
        ));
    }
    Ok(property)
}
