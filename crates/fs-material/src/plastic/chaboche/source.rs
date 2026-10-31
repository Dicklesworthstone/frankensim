//! Fixed-condition sources for one independently fitted Chaboche backstress.
//!
//! The branch follows `dX = (2/3) C d(epsilon_p) - gamma X d(alpha)`, with
//! the dynamic recovery convention described by the
//! [MOOSE combined nonlinear hardening model](https://mooseframework.inl.gov/source/materials/CombinedNonlinearHardeningPlasticity.html).
//! C is a plastic-strain modulus [Pa]; gamma is a recovery coefficient per unit
//! equivalent plastic strain. Neither is a temporal rate or a saturation stress.
//! Both are independent of the Voce Q/b and optional linear Prager modulus.
//! A monotonic total-strain tangent cannot identify this nonlinear cyclic law.
//!
//! The complete source bundle, including the original Voce and other physics
//! receipts, survives adoption. Conditions remain fixed for the entire load
//! path. Resolving another temperature does not transfer history or establish
//! a thermo-plastic energy coupling, static recovery or creep model.

use fs_evidence::ModelCard;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Dims, Pressure, QuantitySpec};

use super::{BackstressBranch, ChabocheJ2, ChabocheState};
use crate::plastic::source::{HardeningPropertyBasis, J2SourceError};
use crate::plastic::voce::source::{ResolvedVoceJ2, voce_requirements};
use crate::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ResolvedScalarProperty,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};

/// Independently named, ordered sources for one to eight retained backstresses.
pub mod branches;
pub use branches::{
    BackstressProperties, ResolvedChabocheBranches, chaboche_branch_requirements,
    resolve_chaboche_branches,
};

/// Independent branch modulus C [Pa], in `dX = (2/3) C d(epsilon_p) - gamma X d(alpha)`.
pub const CHABOCHE_BACKSTRESS_MODULUS_PROPERTY: &str = "chaboche_backstress_modulus";
/// Dynamic recovery gamma per unit equivalent plastic strain [dimensionless].
/// A source with a time-rate, stress or evolving strain coordinate is not this datum.
pub const CHABOCHE_DYNAMIC_RECOVERY_PROPERTY: &str = "chaboche_dynamic_recovery";

/// Complete Voce requirements plus independently sourced branch C and gamma.
/// Explicit C=0 gives an inactive branch; explicit gamma=0 gives linear
/// kinematic hardening. Missing data never select either limit implicitly.
///
/// # Errors
/// A total-strain tangent cannot identify independent nonlinear cyclic data.
pub fn chaboche_requirements(
    basis: HardeningPropertyBasis,
) -> Result<Vec<ScalarPropertyRequirement>, J2SourceError> {
    admit_basis(basis)?;
    let mut requirements = voce_requirements(basis)?;
    requirements.push(ScalarPropertyRequirement::try_new(
        CHABOCHE_BACKSTRESS_MODULUS_PROPERTY,
        Pressure::DIMS,
        ScalarAdmissibility::NonNegative,
    )?);
    requirements.push(ScalarPropertyRequirement::try_new(
        CHABOCHE_DYNAMIC_RECOVERY_PROPERTY,
        Dims([0; 6]),
        ScalarAdmissibility::NonNegative,
    )?);
    Ok(requirements)
}

/// One native nonlinear backstress and its complete immutable material snapshot.
/// The existing Voce owner retains all original receipts; this owner adds the
/// two independently selected branch coefficients without fitting or inference.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedChabocheJ2 {
    backbone: ResolvedVoceJ2,
    law: ChabocheJ2<1>,
}

impl ResolvedChabocheJ2 {
    /// Adopt the original shared bundle, including all additional physics roles.
    /// Supported curves and surfaces at fixed external conditions are admitted;
    /// mechanical history, stress, strain, frequency and rate conditions cannot
    /// be frozen while the native plastic and backstress histories evolve.
    ///
    /// # Errors
    /// Missing or wrong-quantity coefficients, evolving source coordinates,
    /// invalid native parameters and a tangent-only hardening basis refuse.
    /// The original J2 owner enforces the strain ceiling in `(0, 0.05]`.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        basis: HardeningPropertyBasis,
        strain_limit: f64,
    ) -> Result<Self, J2SourceError> {
        admit_basis(basis)?;
        let modulus_pa = coefficient(
            &material,
            CHABOCHE_BACKSTRESS_MODULUS_PROPERTY,
            Pressure::DIMS,
        )?
        .value_si();
        let recovery_per_strain =
            coefficient(&material, CHABOCHE_DYNAMIC_RECOVERY_PROPERTY, Dims([0; 6]))?.value_si();
        let branch = BackstressBranch::new(modulus_pa, recovery_per_strain)?;
        let backbone = ResolvedVoceJ2::try_from_material(material, basis, strain_limit)?;
        let law = ChabocheJ2::new(backbone.law().clone(), [branch])?;
        Ok(Self { backbone, law })
    }

    /// Complete original snapshot, with original quantities, claims and usage
    /// receipts for the cyclic coefficients, the Voce law and every extra role.
    #[must_use]
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        self.backbone.material()
    }

    /// Original source-bound Voce owner, preserving the complete shared bundle.
    #[must_use]
    pub const fn backbone(&self) -> &ResolvedVoceJ2 {
        &self.backbone
    }

    /// Native implicit cyclic update, retained backstress and consistent tangent.
    #[must_use]
    pub const fn law(&self) -> &ChabocheJ2<1> {
        &self.law
    }

    /// Explicit selection of linear isotropic or mixed isotropic/Prager data.
    #[must_use]
    pub const fn hardening_basis(&self) -> HardeningPropertyBasis {
        self.backbone.hardening_basis()
    }

    /// Original branch C claim, including uncertainty and usage receipt.
    #[must_use]
    pub fn modulus_source(&self) -> &ResolvedScalarProperty {
        self.material()
            .property(CHABOCHE_BACKSTRESS_MODULUS_PROPERTY)
            .expect("the immutable Chaboche bundle retains its admitted branch modulus")
    }

    /// Original gamma claim, including uncertainty and usage receipt.
    #[must_use]
    pub fn recovery_source(&self) -> &ResolvedScalarProperty {
        self.material()
            .property(CHABOCHE_DYNAMIC_RECOVERY_PROPERTY)
            .expect("the immutable Chaboche bundle retains its admitted dynamic recovery")
    }
}

impl SmallStrainLaw for ResolvedChabocheJ2 {
    type State = ChabocheState<1>;

    fn initial_state(&self) -> ChabocheState<1> {
        self.law.initial_state()
    }

    fn stress(&self, strain: &Voigt, state: &ChabocheState<1>) -> Voigt {
        self.law.stress(strain, state)
    }

    fn tangent(&self, strain: &Voigt, state: &ChabocheState<1>) -> Tangent6 {
        self.law.tangent(strain, state)
    }

    fn update_state(&self, strain: &Voigt, state: &ChabocheState<1>) -> ChabocheState<1> {
        self.law.update_state(strain, state)
    }

    fn admissibility(&self) -> MaterialAdmissibility {
        self.law.admissibility()
    }

    fn card(&self) -> ModelCard {
        self.law.card()
    }
}

/// Resolve eight isotropic or nine mixed coefficients at one fixed source point.
/// Every native constitutive operation delegates to the Chaboche owner; a
/// material name or monotonic tangent never supplies missing cyclic parameters.
///
/// # Errors
/// Source selection, quantity, validity and native parameter failures refuse.
pub fn resolve_chaboche_j2(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    basis: HardeningPropertyBasis,
    strain_limit: f64,
) -> Result<ResolvedChabocheJ2, J2SourceError> {
    let material =
        resolve_material_state_point(card, point, &chaboche_requirements(basis)?, selection)?;
    ResolvedChabocheJ2::try_from_material(material, basis, strain_limit)
}

fn admit_basis(basis: HardeningPropertyBasis) -> Result<(), J2SourceError> {
    if basis == HardeningPropertyBasis::UniaxialTangentModulus {
        return Err(J2SourceError::Invalid(
            "Chaboche hardening requires independent plastic-strain moduli, Voce Q/b and branch C/gamma; a total-strain tangent cannot identify nonlinear cyclic recovery",
        ));
    }
    Ok(())
}

fn coefficient<'a>(
    material: &'a ResolvedMaterialStatePoint,
    name: &str,
    dims: Dims,
) -> Result<&'a ResolvedScalarProperty, J2SourceError> {
    let property = material.property(name).ok_or(J2SourceError::Invalid(
        "missing or ambiguous Chaboche branch coefficient",
    ))?;
    if property.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !property.value_si().is_finite()
        || !property.has_only_fixed_external_conditions(None, false)
        || property.angular_frequency_band_rad_s()?.is_some()
    {
        return Err(J2SourceError::Invalid(
            "Chaboche needs absolute SI coefficients at fixed external conditions without evolving mechanical, history, rate or frequency coordinates",
        ));
    }
    Ok(property)
}
