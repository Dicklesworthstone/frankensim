//! Original calibration of post-initiation equivalent plastic work per area.
//!
//! The exact source property supplies degraded equivalent-work capacity for
//! attenuation driven by `h * integral (sigma_y + H_iso * alpha) d(alpha)`.
//! Resolution retains the default exponential shape. An optional explicit
//! builder authors a different normalized profile while preserving that capacity
//! receipt; the receipt alone does not validate the new shape. This capacity is
//! neither total fracture dissipation nor a cohesive traction-opening work.
//! Damage releases stored energy separately; it is not added to plastic heat
//! through this property. Toughness, engineering elongation and other fracture
//! energy properties are never substituted for this model-specific datum.
//!
//! Native J2 coefficients, the work calibration and the initiation locus must
//! retain the same original material card and fixed external conditions. The
//! receiving geometry supplies the characteristic length explicitly. No width,
//! crack plane, source value or history transfer is inferred here.

use super::{
    DuctileEvolutionError, DuctileJ2, DuctileJ2Conditions, DuctileJ2Law, DuctileJ2Settings,
    DuctileWorkEvolution,
};
use crate::ductile::source::ResolvedDuctileLocus;
use crate::plastic::energy::EnergyTolerances;
use crate::plastic::source::{
    HardeningPropertyBasis, J2SourceError, ResolvedJ2Plasticity, j2_requirements,
};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Dims, QuantitySpec};

/// Calibration of the post-initiation equivalent plastic-work attenuation.
/// This distinct property does not denote total fracture loss or cohesive Gc.
pub const POST_INITIATION_EQUIVALENT_PLASTIC_WORK_PROPERTY: &str =
    "ductile_post_initiation_equivalent_plastic_work_per_area";
/// Work per area [J/m2]; equal dimensions do not authorize other property names.
pub const PLASTIC_WORK_PER_AREA_DIMS: Dims = Dims([0, 1, -2, 0, 0, 0]);

/// Receiving geometry and explicit numerical gates, with no calibration override.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctileJ2SourceSettings {
    /// Fixed original crack-band width [m], supplied by the geometry owner.
    pub characteristic_length_m: f64,
    /// Maximum permitted damage change in one accepted physical increment.
    pub maximum_damage_increment: f64,
    /// Native yield and independent work-account tolerances.
    pub energy_tolerances: EnergyTolerances,
}

/// Compose the exact work calibration with the original native J2 requirements.
/// Additional multiphysics properties may remain in the resolved bundle.
pub fn ductile_j2_requirements(
    basis: HardeningPropertyBasis,
) -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    let mut requirements = j2_requirements(basis)?;
    requirements.push(ScalarPropertyRequirement::try_new(
        POST_INITIATION_EQUIVALENT_PLASTIC_WORK_PROPERTY,
        PLASTIC_WORK_PER_AREA_DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )?);
    Ok(requirements)
}

/// The original J2/material bundle and its source-bound ductile evolution law.
/// Runtime updates use the prepared numeric law without querying source data.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedDuctileJ2 {
    plastic: ResolvedJ2Plasticity,
    law: DuctileJ2<ResolvedDuctileLocus>,
}

impl ResolvedDuctileJ2 {
    /// Adopt native J2 and initiation sources without selecting either again.
    ///
    /// The J2 bundle must also retain the explicitly named post-initiation work
    /// property. Each calibration axis must describe a supported fixed external
    /// condition; evolving rate, damage, strain and stress restrictions refuse.
    /// Source identity and every held coordinate agree in both directions with
    /// the initiation locus, including typed temperature and extra conditions.
    pub fn try_from_j2(
        plastic: ResolvedJ2Plasticity,
        initiation: ResolvedDuctileLocus,
        settings: DuctileJ2SourceSettings,
    ) -> Result<Self, DuctileJ2SourceError> {
        if !initiation.source().matches_material(plastic.material()) {
            return Err(DuctileJ2SourceError::Invalid(
                "ductile evolution requires the original J2 card, material state and complete fixed initiation context",
            ));
        }
        let property = plastic
            .material()
            .property(POST_INITIATION_EQUIVALENT_PLASTIC_WORK_PROPERTY)
            .ok_or(DuctileJ2SourceError::Invalid(
                "missing or ambiguous post-initiation equivalent plastic-work source",
            ))?;
        let requirement = property.requirement();
        if requirement.quantity() != QuantitySpec::dimensional(PLASTIC_WORK_PER_AREA_DIMS)
            || requirement.hardness_test().is_some()
            || requirement.elastic_component().is_some()
            || requirement.strain_component().is_some()
            || requirement.stress_component().is_some()
            || !property.is_static_at_fixed_external_conditions()
            || !property.value_si().is_finite()
            || property.value_si() <= 0.0
        {
            return Err(DuctileJ2SourceError::Invalid(
                "post-initiation equivalent plastic work requires positive scalar J/m2 data at fixed external conditions",
            ));
        }
        let law = DuctileJ2::new(
            plastic.law().clone(),
            initiation,
            DuctileJ2Settings {
                characteristic_length_m: settings.characteristic_length_m,
                post_initiation_plastic_work_j_m2: property.value_si(),
                maximum_damage_increment: settings.maximum_damage_increment,
                energy_tolerances: settings.energy_tolerances,
            },
        )?;
        Ok(Self { plastic, law })
    }

    /// Prepare native J2 from one complete bundle, preserving extra receipts.
    /// The selected hardening basis and small-strain bound remain explicit.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        basis: HardeningPropertyBasis,
        strain_limit: f64,
        initiation: ResolvedDuctileLocus,
        settings: DuctileJ2SourceSettings,
    ) -> Result<Self, DuctileJ2SourceError> {
        Self::try_from_j2(
            ResolvedJ2Plasticity::try_from_material(material, basis, strain_limit)?,
            initiation,
            settings,
        )
    }

    /// Explicitly author a normalized post-initiation shape before borrowing
    /// histories. The original equivalent-work CAPACITY, native J2, initiation,
    /// typed context and source receipts remain unchanged. A scalar capacity
    /// claim does not establish measured support for this shape or authorize
    /// material removal, loss of compression support, or erased inventories.
    pub fn with_evolution(mut self, evolution: DuctileWorkEvolution) -> Self {
        self.law = self.law.with_evolution(evolution);
        self
    }

    /// Original mechanical coefficients, work calibration and all source receipts.
    #[must_use]
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        self.plastic.material()
    }

    /// Original source-bound J2 model, including its chosen hardening basis.
    #[must_use]
    pub const fn plastic(&self) -> &ResolvedJ2Plasticity {
        &self.plastic
    }

    /// Prepared native law; its initiation locus retains the original observations.
    #[must_use]
    pub const fn law(&self) -> &DuctileJ2<ResolvedDuctileLocus> {
        &self.law
    }

    /// Borrow the same law for a physical increment at actual runtime conditions.
    /// The numeric owner checks the calibration temperature and positive duration.
    pub fn small_strain_law(
        &self,
        conditions: DuctileJ2Conditions,
    ) -> Result<DuctileJ2Law<'_, ResolvedDuctileLocus>, DuctileEvolutionError> {
        self.law.small_strain_law(conditions)
    }
}

/// Resolve the work calibration and mechanical coefficients from one original card.
/// The independently resolved initiation locus must come from that same snapshot.
pub fn resolve_ductile_j2(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    basis: HardeningPropertyBasis,
    strain_limit: f64,
    initiation: ResolvedDuctileLocus,
    settings: DuctileJ2SourceSettings,
) -> Result<ResolvedDuctileJ2, DuctileJ2SourceError> {
    let material =
        resolve_material_state_point(card, point, &ductile_j2_requirements(basis)?, selection)?;
    ResolvedDuctileJ2::try_from_material(material, basis, strain_limit, initiation, settings)
}

/// Original source, mechanical-law or ductile-evolution refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum DuctileJ2SourceError {
    MaterialState(MaterialStatePointError),
    Plastic(J2SourceError),
    Native(DuctileEvolutionError),
    Invalid(&'static str),
}

impl core::fmt::Display for DuctileJ2SourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MaterialState(error) => write!(f, "ductile evolution material source: {error}"),
            Self::Plastic(error) => write!(f, "ductile evolution J2 source: {error}"),
            Self::Native(error) => write!(f, "ductile evolution law: {error}"),
            Self::Invalid(reason) => write!(f, "ductile evolution source: {reason}"),
        }
    }
}
impl std::error::Error for DuctileJ2SourceError {}
impl From<MaterialStatePointError> for DuctileJ2SourceError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::MaterialState(error)
    }
}
impl From<J2SourceError> for DuctileJ2SourceError {
    fn from(error: J2SourceError) -> Self {
        Self::Plastic(error)
    }
}
impl From<DuctileEvolutionError> for DuctileJ2SourceError {
    fn from(error: DuctileEvolutionError) -> Self {
        Self::Native(error)
    }
}
