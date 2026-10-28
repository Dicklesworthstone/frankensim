//! Material-card coefficients for the existing small-strain J2 law.
//!
//! All coefficients are resolved at one physical state. The native radial
//! return, plastic history and consistent tangent remain owned by
//! [`super::J2Plasticity`]. This adapter retains the complete source bundle;
//! it does not turn a derived hardening modulus into a measured property.
//!
//! Conditions are frozen for one load path. Re-resolving at a new temperature
//! gives a different law, not a history-preserving thermo-plastic update.
//! Creep, rate effects, softening, phase changes and finite strain are outside
//! this reduction. The caller must keep the load path within the declared
//! strain bound and the source conditions; the raw law does not police them.

use core::fmt;

use fs_evidence::ModelCard;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Dims, Pressure, QuantitySpec};

use super::{J2Plasticity, J2State};
use crate::elastic::IsotropicElastic;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, POISSON_RATIO_PROPERTY,
    ResolvedMaterialStatePoint, ResolvedScalarProperty, ScalarAdmissibility,
    ScalarPropertyRequirement, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};

/// Slope of the yield stress with respect to equivalent plastic strain [Pa].
pub const ISOTROPIC_HARDENING_MODULUS_PROPERTY: &str = "isotropic_hardening_modulus";
/// Prager modulus H_k [Pa], with backstress a=(2/3)H_k dev(epsilon_p).
pub const KINEMATIC_HARDENING_MODULUS_PROPERTY: &str = "kinematic_hardening_modulus";
/// Monotonic uniaxial post-yield stress slope with respect to total strain [Pa].
pub const POST_YIELD_TANGENT_MODULUS_PROPERTY: &str = "post_yield_tangent_modulus";

/// The explicitly chosen independent hardening datum.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HardeningPropertyBasis {
    /// The native J2 hardening modulus `H = d(sigma_y)/d(alpha)`.
    PlasticStrainModulus,
    /// Independent H_iso and H_k sources. The isotropic property remains
    /// d(sigma_y)/d(alpha); H_k translates the surface through retained plastic
    /// strain. A monotonic total-strain tangent cannot identify this split.
    MixedPlasticStrainModuli,
    /// A uniaxial total-strain tangent `Et`, converted by `H = Et E/(E-Et)`.
    ///
    /// This conversion assumes monotonic uniaxial stress in the same bilinear
    /// rate-independent law. It is not a general multiaxial tangent conversion.
    UniaxialTangentModulus,
}

impl HardeningPropertyBasis {
    /// Exact selected source property name; the other basis is never a fallback.
    #[must_use]
    pub const fn property_name(self) -> &'static str {
        match self {
            Self::PlasticStrainModulus | Self::MixedPlasticStrainModuli => ISOTROPIC_HARDENING_MODULUS_PROPERTY,
            Self::UniaxialTangentModulus => POST_YIELD_TANGENT_MODULUS_PROPERTY,
        }
    }

    /// Include the primary hardening source in a shared material bundle.
    /// For mixed hardening, `j2_requirements` also requires its independent H_k.
    pub fn requirement(self) -> Result<ScalarPropertyRequirement, MaterialStatePointError> {
        ScalarPropertyRequirement::try_new(
            self.property_name(),
            Pressure::DIMS,
            ScalarAdmissibility::NonNegative,
        )
    }
}

/// Exact requirements of this reduction, for composition with other consumers.
/// Perfect plasticity requires an explicit zero hardening datum.
pub fn j2_requirements(
    basis: HardeningPropertyBasis,
) -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
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
            YIELD_STRESS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        basis.requirement()?,
    ];
    if basis == HardeningPropertyBasis::MixedPlasticStrainModuli {
        requirements.push(ScalarPropertyRequirement::try_new(
            KINEMATIC_HARDENING_MODULUS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::NonNegative,
        )?);
    }
    Ok(requirements)
}

/// An immutable native J2 law together with its complete source-use evidence.
///
/// Implements [`SmallStrainLaw`] directly so native uniform or per-element
/// solvers can borrow the complete resolved inputs. Constitutive operations
/// delegate to the same J2 owner; the material receipts remain available through
/// [`Self::material`] alongside the native model assumptions.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedJ2Plasticity {
    material: ResolvedMaterialStatePoint,
    basis: HardeningPropertyBasis,
    law: J2Plasticity,
}

impl ResolvedJ2Plasticity {
    /// Adopt an already resolved bundle without losing other physics receipts.
    ///
    /// `basis` remains explicit even if the bundle contains both independent
    /// hardening sources. They are not reconciled or averaged. Curves and
    /// surfaces at supported fixed external conditions are evaluated by the
    /// common material owner before conversion. Mechanical source coordinates
    /// and rates cannot be frozen into these coefficients.
    ///
    /// # Errors
    /// Missing/ambiguous sources, wrong quantities, unsupported coefficient
    /// forms, invalid parameters or nonrepresentable native elastic/plastic
    /// coefficients refuse. `strain_limit` must be finite and in `(0, 0.05]`.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        basis: HardeningPropertyBasis,
        strain_limit: f64,
    ) -> Result<Self, J2SourceError> {
        if !(strain_limit.is_finite() && strain_limit > 0.0 && strain_limit <= 0.05) {
            return Err(J2SourceError::Invalid(
                "small-strain bound must be in (0, 0.05]",
            ));
        }
        let youngs = coefficient(&material, YOUNG_MODULUS_PROPERTY, Pressure::DIMS)?.value_si();
        let poisson = coefficient(&material, POISSON_RATIO_PROPERTY, Dims([0; 6]))?.value_si();
        let yield_stress =
            coefficient(&material, YIELD_STRESS_PROPERTY, Pressure::DIMS)?.value_si();
        let source_hardening =
            coefficient(&material, basis.property_name(), Pressure::DIMS)?.value_si();
        if !(source_hardening >= 0.0 && source_hardening.is_finite()) {
            return Err(J2SourceError::Invalid(
                "hardening source must be finite and nonnegative",
            ));
        }
        let elastic = IsotropicElastic::new(youngs, poisson, strain_limit)?;
        let hardening = match basis {
            HardeningPropertyBasis::PlasticStrainModulus
                | HardeningPropertyBasis::MixedPlasticStrainModuli => source_hardening,
            HardeningPropertyBasis::UniaxialTangentModulus => {
                if source_hardening >= youngs {
                    return Err(J2SourceError::Invalid(
                        "uniaxial tangent must be smaller than Young's modulus",
                    ));
                }
                // Evaluate the ratio first: E*Et can overflow even when H fits.
                source_hardening * (youngs / (youngs - source_hardening))
            }
        };
        if !(hardening.is_finite()
            && hardening >= 0.0
            && (source_hardening == 0.0 || hardening > 0.0))
        {
            return Err(J2SourceError::Invalid(
                "hardening conversion is not representable",
            ));
        }
        let kinematic = if basis == HardeningPropertyBasis::MixedPlasticStrainModuli {
            coefficient(&material, KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS)?.value_si()
        } else { 0.0 };
        let law = J2Plasticity::new_mixed(elastic, yield_stress, hardening, kinematic)?;
        let (lambda, mu) = law.elastic.lame();
        let bulk = lambda + 2.0 * mu / 3.0;
        let denominator = 3.0 * mu + law.combined_hardening();
        if !(mu.is_finite()
            && mu > 0.0
            && bulk.is_finite()
            && bulk > 0.0
            && denominator.is_finite()
            && denominator > 0.0)
            || !law
                .elastic
                .stiffness()
                .iter()
                .flatten()
                .all(|v| v.is_finite())
        {
            return Err(J2SourceError::Invalid(
                "native elastic or return-map coefficients are not representable",
            ));
        }
        Ok(Self {
            material,
            basis,
            law,
        })
    }

    /// Complete original material snapshot, including additional properties.
    #[must_use]
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }

    /// The existing law for native small-strain solvers, including 3-D solids.
    #[must_use]
    pub const fn law(&self) -> &J2Plasticity {
        &self.law
    }

    /// The independent source quantity used in this reduction.
    #[must_use]
    pub const fn hardening_basis(&self) -> HardeningPropertyBasis {
        self.basis
    }

    /// Native yield-surface hardening modulus `H` [Pa].
    #[must_use]
    pub const fn hardening_modulus_pa(&self) -> f64 {
        self.law.hardening
    }

    /// Separately sourced Prager modulus, zero in the original isotropic modes.
    #[must_use]
    pub const fn kinematic_hardening_modulus_pa(&self) -> f64 {
        self.law.kinematic_hardening
    }

    /// Selected independent H_k source, only when mixed hardening is active.
    #[must_use]
    pub fn kinematic_hardening_source(&self) -> Option<&ResolvedScalarProperty> {
        if self.basis == HardeningPropertyBasis::MixedPlasticStrainModuli {
            self.material.property(KINEMATIC_HARDENING_MODULUS_PROPERTY)
        } else { None }
    }

    /// Original H or Et source, retaining its selection receipt and uncertainty.
    #[must_use]
    pub fn hardening_source(&self) -> &ResolvedScalarProperty {
        self.material
            .property(self.basis.property_name())
            .expect("the immutable J2 bundle retains its admitted hardening source")
    }
}

impl SmallStrainLaw for ResolvedJ2Plasticity {
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

/// Resolve all four isotropic or five mixed coefficients at one material state.
/// No material name, absent hardening source or unsupported state supplies a
/// default metal or an implicit perfectly plastic law.
///
/// # Errors
/// Source selection, dimensions, validity, coefficient or conversion refusals.
pub fn resolve_j2_plasticity(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    basis: HardeningPropertyBasis,
    strain_limit: f64,
) -> Result<ResolvedJ2Plasticity, J2SourceError> {
    let material = resolve_material_state_point(card, point, &j2_requirements(basis)?, selection)?;
    ResolvedJ2Plasticity::try_from_material(material, basis, strain_limit)
}

/// Source admission and native constitutive-parameter failures.
#[derive(Clone, Debug, PartialEq)]
pub enum J2SourceError {
    /// The shared source resolver refused the requested state or claim.
    MaterialState(MaterialStatePointError),
    /// The existing elastic/J2 owner refused a parameter.
    Parameters(MaterialError),
    /// The requested reduction or a derived coefficient is inadmissible.
    Invalid(&'static str),
}

impl fmt::Display for J2SourceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MaterialState(source) => write!(f, "J2 source refusal: {source}"),
            Self::Parameters(source) => write!(f, "J2 parameter refusal: {source}"),
            Self::Invalid(what) => write!(f, "J2 reduction refusal: {what}"),
        }
    }
}

impl std::error::Error for J2SourceError {}

impl From<MaterialStatePointError> for J2SourceError {
    fn from(source: MaterialStatePointError) -> Self {
        Self::MaterialState(source)
    }
}

impl From<MaterialError> for J2SourceError {
    fn from(source: MaterialError) -> Self {
        Self::Parameters(source)
    }
}

fn coefficient<'a>(
    material: &'a ResolvedMaterialStatePoint,
    name: &'static str,
    dims: Dims,
) -> Result<&'a ResolvedScalarProperty, J2SourceError> {
    let property = material.property(name).ok_or(J2SourceError::Invalid(
        "missing or ambiguous J2 coefficient",
    ))?;
    if property.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !property.value_si().is_finite()
        || !property.has_only_fixed_external_conditions(None, false)
    {
        return Err(J2SourceError::Invalid(
            "J2 needs absolute SI coefficients at fixed external conditions, without evolving mechanical or frequency coordinates",
        ));
    }
    if property.angular_frequency_band_rad_s()?.is_some() {
        return Err(J2SourceError::Invalid(
            "frequency-bounded coefficients cannot define this rate-independent law",
        ));
    }
    Ok(property)
}
