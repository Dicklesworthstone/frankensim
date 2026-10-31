//! Fixed-condition material sources for the native Voce J2 law.
//!
//! The saturation contribution is `Q * (1 - exp(-b * alpha))`, following the
//! [Voce convention in the Sandia LAME manual](https://www.sandia.gov/files/sierra/SM_LAME_5_28/main/models/j2_plasticity.html#voce).
//! `Q` is an additional yield-stress increment, not the absolute saturation
//! stress. The independent linear isotropic and optional Prager moduli keep
//! their original meanings. A single post-yield tangent cannot identify these
//! coefficients and is not converted into a Voce hardening law.
//!
//! The complete original source bundle survives this reduction, including
//! other physics receipts. Coefficients are fixed for a load path. Resolving
//! them at another temperature constructs a different law; it does not provide
//! history transfer, thermo-plastic energy exchange, creep or phase evolution.

use fs_evidence::ModelCard;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Dims, Pressure, QuantitySpec};

use super::VoceJ2;
use crate::plastic::J2State;
use crate::plastic::source::{
    HardeningPropertyBasis, J2SourceError, ResolvedJ2Plasticity, j2_requirements,
};
use crate::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ResolvedScalarProperty,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};

/// Additional asymptotic isotropic hardening Q [Pa], excluding initial yield
/// stress and the independently sourced linear hardening term.
pub const VOCE_SATURATION_INCREMENT_PROPERTY: &str = "voce_saturation_increment";
/// Saturation rate b per unit equivalent plastic strain [dimensionless]. This
/// is not a temporal rate or an independently evolving strain coordinate.
pub const VOCE_SATURATION_RATE_PROPERTY: &str = "voce_saturation_rate";

/// Exact native elastic, linear hardening and independent Voce requirements.
/// An explicit Q=0 retains the original J2 limit; an absent datum is not zero.
///
/// # Errors
/// A total-strain tangent basis cannot identify a saturation law.
pub fn voce_requirements(
    basis: HardeningPropertyBasis,
) -> Result<Vec<ScalarPropertyRequirement>, J2SourceError> {
    admit_basis(basis)?;
    let mut requirements = j2_requirements(basis)?;
    requirements.push(ScalarPropertyRequirement::try_new(
        VOCE_SATURATION_INCREMENT_PROPERTY,
        Pressure::DIMS,
        ScalarAdmissibility::NonNegative,
    )?);
    requirements.push(ScalarPropertyRequirement::try_new(
        VOCE_SATURATION_RATE_PROPERTY,
        Dims([0; 6]),
        ScalarAdmissibility::StrictlyPositive,
    )?);
    Ok(requirements)
}

/// Native nonlinear J2 response with its complete immutable material snapshot.
/// The original J2 source owner admits the backbone and retains every receipt;
/// this adapter adds only the two independently selected saturation parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedVoceJ2 {
    backbone: ResolvedJ2Plasticity,
    law: VoceJ2,
}

impl ResolvedVoceJ2 {
    /// Adopt a shared material bundle without dropping unrelated source roles.
    /// Static curves or surfaces at supported fixed external conditions are
    /// admitted; evolving strain, stress, rate, history and frequency conditions
    /// cannot be frozen into constitutive coefficients.
    ///
    /// # Errors
    /// Missing or wrong-quantity coefficients, unsupported source conditions,
    /// nonrepresentable parameters or a tangent-only hardening basis refuse.
    /// The original J2 owner enforces the strain ceiling in `(0, 0.05]`.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        basis: HardeningPropertyBasis,
        strain_limit: f64,
    ) -> Result<Self, J2SourceError> {
        admit_basis(basis)?;
        let saturation_pa = coefficient(
            &material,
            VOCE_SATURATION_INCREMENT_PROPERTY,
            Pressure::DIMS,
        )?
        .value_si();
        let rate_per_strain =
            coefficient(&material, VOCE_SATURATION_RATE_PROPERTY, Dims([0; 6]))?.value_si();
        let backbone = ResolvedJ2Plasticity::try_from_material(material, basis, strain_limit)?;
        let law = VoceJ2::new(backbone.law().clone(), saturation_pa, rate_per_strain)?;
        Ok(Self { backbone, law })
    }

    /// Complete original snapshot, with original quantities, claims and usage
    /// receipts for both saturation parameters and every additional property.
    #[must_use]
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        self.backbone.material()
    }

    /// Original source-bound elastic/linear-hardening owner. Its bundle also
    /// retains the Voce parameters; no second resolution or attribution occurs.
    #[must_use]
    pub const fn backbone(&self) -> &ResolvedJ2Plasticity {
        &self.backbone
    }

    /// Native return map, consistent tangent and energetic saturation model.
    #[must_use]
    pub const fn law(&self) -> &VoceJ2 {
        &self.law
    }

    /// Explicit selection of linear isotropic or mixed isotropic/Prager data.
    #[must_use]
    pub const fn hardening_basis(&self) -> HardeningPropertyBasis {
        self.backbone.hardening_basis()
    }

    /// Original Q claim, retaining its source uncertainty and usage receipt.
    #[must_use]
    pub fn saturation_source(&self) -> &ResolvedScalarProperty {
        self.material()
            .property(VOCE_SATURATION_INCREMENT_PROPERTY)
            .expect("the immutable Voce bundle retains its admitted saturation increment")
    }

    /// Original b claim, retaining its source uncertainty and usage receipt.
    #[must_use]
    pub fn rate_source(&self) -> &ResolvedScalarProperty {
        self.material()
            .property(VOCE_SATURATION_RATE_PROPERTY)
            .expect("the immutable Voce bundle retains its admitted saturation rate")
    }
}

impl SmallStrainLaw for ResolvedVoceJ2 {
    type State = J2State;

    fn initial_state(&self) -> J2State {
        self.law.initial_state()
    }

    fn stress(&self, strain: &Voigt, state: &J2State) -> Voigt {
        self.law.stress(strain, state)
    }

    fn tangent(&self, strain: &Voigt, state: &J2State) -> Tangent6 {
        self.law.tangent(strain, state)
    }

    fn update_state(&self, strain: &Voigt, state: &J2State) -> J2State {
        self.law.update_state(strain, state)
    }

    fn admissibility(&self) -> MaterialAdmissibility {
        self.law.admissibility()
    }

    fn card(&self) -> ModelCard {
        self.law.card()
    }
}

/// Resolve six isotropic or seven mixed coefficients at the same source point.
/// Constitutive updates delegate to the native Voce owner; no material name,
/// missing measurement or monotonic tangent supplies a saturation parameter.
///
/// # Errors
/// Source selection, quantity, validity and native parameter failures refuse.
pub fn resolve_voce_j2(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    basis: HardeningPropertyBasis,
    strain_limit: f64,
) -> Result<ResolvedVoceJ2, J2SourceError> {
    let material =
        resolve_material_state_point(card, point, &voce_requirements(basis)?, selection)?;
    ResolvedVoceJ2::try_from_material(material, basis, strain_limit)
}

fn admit_basis(basis: HardeningPropertyBasis) -> Result<(), J2SourceError> {
    if basis == HardeningPropertyBasis::UniaxialTangentModulus {
        return Err(J2SourceError::Invalid(
            "Voce hardening requires independent plastic-strain moduli and Q/b; a total-strain tangent cannot identify saturation",
        ));
    }
    Ok(())
}

fn coefficient<'a>(
    material: &'a ResolvedMaterialStatePoint,
    name: &'static str,
    dims: Dims,
) -> Result<&'a ResolvedScalarProperty, J2SourceError> {
    let property = material.property(name).ok_or(J2SourceError::Invalid(
        "missing or ambiguous Voce coefficient",
    ))?;
    if property.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !property.value_si().is_finite()
        || !property.has_only_fixed_external_conditions(None, false)
        || property.angular_frequency_band_rad_s()?.is_some()
    {
        return Err(J2SourceError::Invalid(
            "Voce needs absolute SI coefficients at fixed external conditions without evolving mechanical, history, rate or frequency coordinates",
        ));
    }
    Ok(property)
}
