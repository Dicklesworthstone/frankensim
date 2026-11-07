//! Original material-card hardening curves on the native J2 owner.
//!
//! The exact abscissa is accumulated equivalent plastic strain. The two exact
//! property names distinguish Cauchy and Kirchhoff yield stress; an engineering
//! tensile curve, total strain or pressure claim cannot supply an alias.
//! Preparation selects one immutable claim at zero plastic strain, pins it,
//! and queries every original knot at the same fixed external conditions.
//! Piecewise-linear constitutive interpolation is explicit model behavior;
//! original observations, uncertainty and source interpolation policy survive.
//! No source query or receipt allocation occurs during a constitutive return.
//!
//! These are fixed-condition laws. A newly resolved temperature is a different
//! law, not a history transfer or a thermal work/entropy account.

use super::{HardeningPoint, TabulatedError, TabulatedHardening, TabulatedJ2, TabulatedJ2Response};
use crate::elastic::IsotropicElastic;
use crate::plastic::energy::{
    EnergyError, EnergyIncrement, EnergyTolerances, IncrementEnergy, StoredEnergy,
};
use crate::plastic::section::energy::NativeSectionEnergyLaw;
use crate::plastic::source::KINEMATIC_HARDENING_MODULUS_PROPERTY;
use crate::plastic::{J2Plasticity, J2State};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, POISSON_RATIO_PROPERTY,
    ResolvedMaterialStatePoint, ResolvedScalarProperty, ScalarAdmissibility,
    ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY, is_fixed_external_coordinate,
    resolve_material_state_point,
};
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
use fs_blake3::ContentHash;
use fs_evidence::ModelCard;
use fs_matdb::{
    ClaimSelection, MatDbError, MaterialAnswer, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, QueryPoint,
};
use fs_qty::{Dims, Pressure, QuantitySpec};

/// Absolute Cauchy yield radius versus accumulated equivalent plastic strain.
pub const CAUCHY_YIELD_CURVE_PROPERTY: &str = "j2_cauchy_yield_stress";
/// Absolute Kirchhoff yield radius in the explicit original-volume finite law.
pub const KIRCHHOFF_YIELD_CURVE_PROPERTY: &str = "j2_kirchhoff_yield_stress";
/// Dimensionless accumulated path, not a signed or total strain component.
pub const ACCUMULATED_PLASTIC_STRAIN_AXIS: &str = "accumulated-equivalent-plastic-strain";

/// Exact source interpretation; no current/reference-volume conversion occurs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum YieldStressMeasure {
    Cauchy,
    Kirchhoff,
}
impl YieldStressMeasure {
    pub const fn property_name(self) -> &'static str {
        match self {
            Self::Cauchy => CAUCHY_YIELD_CURVE_PROPERTY,
            Self::Kirchhoff => KIRCHHOFF_YIELD_CURVE_PROPERTY,
        }
    }
}

/// Original stress convention and a hard bound on source preparation work.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HardeningCurveRequest {
    pub stress_measure: YieldStressMeasure,
    pub maximum_points: usize,
}

/// The original complete curve, fixed context and every queried observation.
#[derive(Clone, Debug, PartialEq)]
pub struct HardeningCurveSource {
    card_identity: ContentHash,
    material_state: MaterialStateId,
    context: QueryPoint,
    request: HardeningCurveRequest,
    selection: ClaimSelection,
    selection_answer: MaterialAnswer,
    claim: PropertyClaim,
    observations: Vec<(QueryPoint, MaterialAnswer)>,
}
impl HardeningCurveSource {
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }
    pub const fn material_state(&self) -> &MaterialStateId {
        &self.material_state
    }
    /// Complete fixed query, excluding the runtime accumulated-plastic-strain axis.
    pub const fn context(&self) -> &QueryPoint {
        &self.context
    }
    pub const fn request(&self) -> &HardeningCurveRequest {
        &self.request
    }
    pub const fn selection(&self) -> ClaimSelection {
        self.selection
    }
    /// Original policy/pin decision at zero plastic strain, before knot pinning.
    pub const fn selection_answer(&self) -> &MaterialAnswer {
        &self.selection_answer
    }
    pub const fn claim(&self) -> &PropertyClaim {
        &self.claim
    }
    /// Original typed coordinates and exact knot answers, in increasing strain.
    pub fn observations(&self) -> &[(QueryPoint, MaterialAnswer)] {
        &self.observations
    }

    pub fn matches_material(&self, material: &ResolvedMaterialStatePoint) -> bool {
        self.card_identity == material.card_identity()
            && &self.material_state == material.material()
            && self.context == material.query()
    }
}

/// Prepared native curve retaining its actual stress convention and evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedTabulatedHardening {
    hardening: TabulatedHardening,
    source: HardeningCurveSource,
}
impl ResolvedTabulatedHardening {
    pub const fn hardening(&self) -> &TabulatedHardening {
        &self.hardening
    }
    pub const fn source(&self) -> &HardeningCurveSource {
        &self.source
    }
    pub const fn stress_measure(&self) -> YieldStressMeasure {
        self.source.request.stress_measure
    }
}

/// Source selection, coordinate, native-law or preparation-domain refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum TabulatedSourceError {
    Source(MatDbError),
    MaterialState(MaterialStatePointError),
    Material(MaterialError),
    Curve(TabulatedError),
    Invalid(&'static str),
}
impl core::fmt::Display for TabulatedSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(e) => write!(f, "tabulated hardening source: {e}"),
            Self::MaterialState(e) => write!(f, "tabulated material state: {e}"),
            Self::Material(e) => write!(f, "tabulated native material: {e}"),
            Self::Curve(e) => e.fmt(f),
            Self::Invalid(e) => write!(f, "tabulated source admission: {e}"),
        }
    }
}
impl std::error::Error for TabulatedSourceError {}
impl From<MatDbError> for TabulatedSourceError {
    fn from(e: MatDbError) -> Self {
        Self::Source(e)
    }
}
impl From<MaterialStatePointError> for TabulatedSourceError {
    fn from(e: MaterialStatePointError) -> Self {
        Self::MaterialState(e)
    }
}
impl From<MaterialError> for TabulatedSourceError {
    fn from(e: MaterialError) -> Self {
        Self::Material(e)
    }
}
impl From<TabulatedError> for TabulatedSourceError {
    fn from(e: TabulatedError) -> Self {
        Self::Curve(e)
    }
}

fn require(accepted: bool, reason: &'static str) -> Result<(), TabulatedSourceError> {
    if accepted {
        Ok(())
    } else {
        Err(TabulatedSourceError::Invalid(reason))
    }
}

fn query(
    card: &MaterialCard,
    key: &PropertyKey,
    point: &QueryPoint,
    selection: ClaimSelection,
) -> Result<MaterialAnswer, TabulatedSourceError> {
    let answer = match selection {
        ClaimSelection::Policy(policy) => card.claims().query_typed(key, point, policy)?,
        ClaimSelection::Pinned(pin) => card.claims().query_pinned_typed(key, point, pin)?,
    };
    card.claims().verify_receipt(&answer.receipt)?;
    Ok(answer)
}

/// Prepare all original knots from one selected immutable source curve.
///
/// `fixed_context` contains only supported static external conditions. The
/// resolver introduces the exact typed accumulated-plastic-strain coordinate
/// itself, selects at zero and pins that claim for the full preparation. Every
/// knot must be covered by the original validity domain. Tables require a zero
/// origin, positive nondecreasing stress and at least two increasing knots.
/// Piecewise-linear interpolation does not change a TabulatedOnly source claim.
pub fn resolve_tabulated_hardening(
    card: &MaterialCard,
    fixed_context: &QueryPoint,
    request: &HardeningCurveRequest,
    selection: ClaimSelection,
) -> Result<ResolvedTabulatedHardening, TabulatedSourceError> {
    require(
        request.maximum_points >= 2,
        "declare at least two supported source points",
    )?;
    for (axis, &value) in fixed_context.axes() {
        require(
            axis != ACCUMULATED_PLASTIC_STRAIN_AXIS
                && is_fixed_external_coordinate(
                    axis,
                    fixed_context.axis_quantities().get(axis).copied(),
                    Some(value),
                ),
            "curve preparation needs only fixed external conditions; strain, stress, rates and history cannot be frozen",
        )?;
    }
    let quantity = QuantitySpec::dimensional(Dims::NONE);
    let origin =
        fixed_context
            .clone()
            .with_quantity(ACCUMULATED_PLASTIC_STRAIN_AXIS, quantity, 0.0)?;
    let key = PropertyKey::new(request.stress_measure.property_name(), Pressure::DIMS);
    let selected = query(card, &key, &origin, selection)?;
    let claim =
        card.claims()
            .claim(selected.receipt.selected)
            .ok_or(TabulatedSourceError::Invalid(
                "selected hardening claim is absent",
            ))?;
    let knots = match &claim.value {
        PropertyValue::Curve {
            abscissa,
            abscissa_dims,
            knots,
            dims,
        } if abscissa == ACCUMULATED_PLASTIC_STRAIN_AXIS
            && *abscissa_dims == Dims::NONE
            && *dims == Pressure::DIMS =>
        {
            knots
        }
        _ => {
            return Err(TabulatedSourceError::Invalid(
                "an original yield-stress curve on accumulated equivalent plastic strain is required",
            ));
        }
    };
    require(
        knots.len() >= 2 && knots.len() <= request.maximum_points,
        "original hardening curve exceeds its preparation budget or has fewer than two knots",
    )?;
    require(
        claim
            .validity
            .axis_quantities()
            .get(ACCUMULATED_PLASTIC_STRAIN_AXIS)
            == Some(&quantity)
            && claim
                .validity
                .bounds()
                .contains_key(ACCUMULATED_PLASTIC_STRAIN_AXIS),
        "the accumulated equivalent plastic-strain axis needs exact dimensionless source bounds",
    )?;
    for axis in claim.validity.bounds().keys() {
        if axis == ACCUMULATED_PLASTIC_STRAIN_AXIS {
            continue;
        }
        require(
            is_fixed_external_coordinate(
                axis,
                fixed_context.axis_quantities().get(axis).copied(),
                fixed_context.axes().get(axis).copied(),
            ),
            "a one-dimensional hardening curve cannot freeze an evolving response or rate coordinate",
        )?;
    }
    let hardening = TabulatedHardening::new(
        knots
            .iter()
            .map(|&(plastic_strain, yield_stress_pa)| HardeningPoint {
                plastic_strain,
                yield_stress_pa,
            })
            .collect(),
    )?;
    let pin = ClaimSelection::Pinned(selected.receipt.selected);
    let mut observations = Vec::with_capacity(knots.len());
    for &(alpha, stress) in knots {
        let at = origin
            .clone()
            .with(ACCUMULATED_PLASTIC_STRAIN_AXIS, alpha)?;
        let answer = query(card, &key, &at, pin)?;
        require(
            answer.receipt.selected == selected.receipt.selected
                && answer.evidence.value.value == stress,
            "every original knot must retain its exact selected claim and yield stress",
        )?;
        observations.push((at, answer));
    }
    Ok(ResolvedTabulatedHardening {
        hardening,
        source: HardeningCurveSource {
            card_identity: card.content_hash(),
            material_state: card.id().clone(),
            context: fixed_context.clone(),
            request: *request,
            selection,
            selection_answer: selected,
            claim: claim.clone(),
            observations,
        },
    })
}

/// The optional independent translation law is an explicit model selection.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TabulatedKinematicHardening {
    /// Explicit isotropic-only selection; no backstress datum is inferred.
    None,
    /// Requires the original independent nonnegative Prager modulus [Pa].
    Prager,
}

/// Elastic/Prager source roles; the curve supplies the complete yield radius.
/// A total-strain tangent or a separate linear isotropic modulus is not used.
pub fn tabulated_j2_requirements(
    kinematic: TabulatedKinematicHardening,
) -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    let mut requirements = vec![
        ScalarPropertyRequirement::try_new(
            YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            POISSON_RATIO_PROPERTY,
            Dims::NONE,
            ScalarAdmissibility::OpenInterval {
                lower: -1.0,
                upper: 0.5,
            },
        )?,
    ];
    if kinematic == TabulatedKinematicHardening::Prager {
        requirements.push(ScalarPropertyRequirement::try_new(
            KINEMATIC_HARDENING_MODULUS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::NonNegative,
        )?);
    }
    Ok(requirements)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TabulatedJ2SourceRequest {
    pub kinematic_hardening: TabulatedKinematicHardening,
    /// Positive Frobenius small-strain limit, no greater than 0.05.
    pub strain_limit: f64,
    pub maximum_points: usize,
}

/// Native small-strain law retaining the complete elastic bundle and curve.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedTabulatedJ2 {
    material: ResolvedMaterialStatePoint,
    curve: ResolvedTabulatedHardening,
    kinematic: TabulatedKinematicHardening,
    law: TabulatedJ2,
}
impl ResolvedTabulatedJ2 {
    /// Bind an already resolved bundle without dropping additional receipts.
    /// The same immutable card, complete fixed context and Cauchy curve are
    /// mandatory. This constructs a virgin law; it never transfers history.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        curve: ResolvedTabulatedHardening,
        kinematic: TabulatedKinematicHardening,
        strain_limit: f64,
    ) -> Result<Self, TabulatedSourceError> {
        require(
            curve.stress_measure() == YieldStressMeasure::Cauchy,
            "small-strain J2 requires an explicit Cauchy yield curve",
        )?;
        require(
            curve.source.matches_material(&material),
            "elastic coefficients and yield curve must retain the same original card and fixed context",
        )?;
        require(
            strain_limit.is_finite() && strain_limit > 0.0 && strain_limit <= 0.05,
            "small-strain norm limit must lie in (0, 0.05]",
        )?;
        let youngs = coefficient(&material, YOUNG_MODULUS_PROPERTY, Pressure::DIMS)?.value_si();
        let poisson = coefficient(&material, POISSON_RATIO_PROPERTY, Dims::NONE)?.value_si();
        let prager = if kinematic == TabulatedKinematicHardening::Prager {
            coefficient(
                &material,
                KINEMATIC_HARDENING_MODULUS_PROPERTY,
                Pressure::DIMS,
            )?
            .value_si()
        } else {
            0.0
        };
        let backbone = J2Plasticity::new_mixed(
            IsotropicElastic::new(youngs, poisson, strain_limit)?,
            curve.hardening.initial_yield_stress_pa(),
            0.0,
            prager,
        )?;
        let law = TabulatedJ2::new(backbone, curve.hardening.clone())?;
        Ok(Self {
            material,
            curve,
            kinematic,
            law,
        })
    }
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    pub const fn curve(&self) -> &ResolvedTabulatedHardening {
        &self.curve
    }
    pub const fn law(&self) -> &TabulatedJ2 {
        &self.law
    }
    pub const fn kinematic_hardening(&self) -> TabulatedKinematicHardening {
        self.kinematic
    }
    pub fn response(
        &self,
        strain: &Voigt,
        state: &J2State,
    ) -> Result<TabulatedJ2Response, TabulatedError> {
        self.law.response(strain, state)
    }
    pub fn validate_committed_history(
        &self,
        strain: &Voigt,
        state: &J2State,
        tolerance_pa: f64,
    ) -> Result<(), EnergyError> {
        self.law
            .validate_committed_history(strain, state, tolerance_pa)
    }
    pub fn stored_energy(
        &self,
        strain: &Voigt,
        state: &J2State,
    ) -> Result<StoredEnergy, EnergyError> {
        self.law.stored_energy(strain, state)
    }
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<EnergyIncrement, EnergyError> {
        self.law
            .advance_with_energy(old_strain, old_state, strain, tolerances)
    }
}
impl SmallStrainLaw for ResolvedTabulatedJ2 {
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
impl NativeSectionEnergyLaw for ResolvedTabulatedJ2 {
    fn native_energy_increment(
        &self,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<(J2State, Voigt, IncrementEnergy), MaterialError> {
        let increment = self
            .advance_with_energy(old_strain, old_state, strain, tolerances)
            .map_err(|error| MaterialError::State {
                what: error.to_string(),
            })?;
        Ok((increment.state, increment.stress_pa, increment.energy))
    }
}

pub fn resolve_tabulated_j2(
    card: &MaterialCard,
    fixed_context: &QueryPoint,
    request: &TabulatedJ2SourceRequest,
    material_selection: MaterialPropertySelection,
    curve_selection: ClaimSelection,
) -> Result<ResolvedTabulatedJ2, TabulatedSourceError> {
    let curve = resolve_tabulated_hardening(
        card,
        fixed_context,
        &HardeningCurveRequest {
            stress_measure: YieldStressMeasure::Cauchy,
            maximum_points: request.maximum_points,
        },
        curve_selection,
    )?;
    let material = resolve_material_state_point(
        card,
        fixed_context,
        &tabulated_j2_requirements(request.kinematic_hardening)?,
        material_selection,
    )?;
    ResolvedTabulatedJ2::try_from_material(
        material,
        curve,
        request.kinematic_hardening,
        request.strain_limit,
    )
}

fn coefficient<'a>(
    material: &'a ResolvedMaterialStatePoint,
    name: &str,
    dims: Dims,
) -> Result<&'a ResolvedScalarProperty, TabulatedSourceError> {
    let property = material
        .property(name)
        .ok_or(TabulatedSourceError::Invalid(
            "missing or ambiguous elastic/Prager coefficient",
        ))?;
    require(
        property.requirement().quantity() == QuantitySpec::dimensional(dims)
            && property.value_si().is_finite()
            && property.has_only_fixed_external_conditions(None, false)
            && property.angular_frequency_band_rad_s()?.is_none(),
        "elastic/Prager coefficients need exact SI quantities at fixed external conditions",
    )?;
    Ok(property)
}
