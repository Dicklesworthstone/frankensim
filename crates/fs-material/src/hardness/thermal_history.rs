//! Separately sourced indentation hardness on accepted native thermoplastic history.
//!
//! The original `CaloricPoint` owns temperature, accumulated plastic strain,
//! complete cyclic backstresses, coefficient transport, physical heat and work.
//! This observer pins one independent indentation test and evaluates H(T, alpha)
//! at that owner's accepted endpoints. It neither converts yield stress to
//! hardness nor adds hardness changes to any mechanical or thermal account.
//!
//! `advance` returns the original feedback report unchanged, together with its
//! typed hardness endpoint, only after both owners accept. Transport callbacks
//! retain the native candidate-only contract; callers commit their finite
//! reservoir report and the returned point together. Rejection leaves the input
//! point and source untouched. There is no additional return map or heat solve.
//!
//! These are endpoint observations, not an integration of hardness along the
//! physical loading path: source support between accepted endpoints is not
//! certified. In particular, H(T, alpha) alone supplies no annealing kinetics,
//! dwell-history evolution, phase transformation, indentation size correction or
//! wear law. Native J2, Voce, viscous and N-branch Chaboche histories retain their
//! existing physical scope. Synthetic numerical checks do not qualify a material.

use core::{convert::Infallible, marker::PhantomData};

use fs_blake3::ContentHash;
use fs_matdb::{MatDbError, PropertyKey, QueryPoint};
use fs_qty::{Dims, QuantitySpec};

use super::IndentationHardnessStatePoint;
use super::plastic_history::ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA;
use crate::phase::{UniformEnthalpyStep, UniformEnthalpyStepInput};
use crate::plastic::J2State;
use crate::plastic::chaboche::ChabocheState;
use crate::plastic::source::ResolvedJ2Plasticity;
use crate::plastic::thermal::family::ThermoplasticSourceFamily;
use crate::plastic::thermal::{
    CaloricPoint, FeedbackConfig, FeedbackError, FeedbackStep, ThermoplasticSource,
};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use crate::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
use crate::{SmallStrainLaw, Voigt};

mod sealed {
    pub trait NativeHistory {}
    impl NativeHistory for crate::plastic::J2State {}
    impl<const N: usize> NativeHistory for crate::plastic::chaboche::ChabocheState<N> {}
}

/// Read the accumulated equivalent strain from a complete native history.
/// This sealed view cannot substitute total strain or fabricate a history.
pub trait AccumulatedPlasticHistory: sealed::NativeHistory {
    fn accumulated_equivalent_plastic_strain(&self) -> f64;
}
impl AccumulatedPlasticHistory for J2State {
    fn accumulated_equivalent_plastic_strain(&self) -> f64 {
        self.alpha
    }
}
impl<const N: usize> AccumulatedPlasticHistory for ChabocheState<N> {
    fn accumulated_equivalent_plastic_strain(&self) -> f64 {
        self.plastic.alpha
    }
}

/// One original indentation test pinned to a native thermal source family.
/// The same binding can observe accepted points from that original source bundle;
/// it does not own a second plastic history or a second thermal clock.
#[derive(Clone, Debug)]
pub struct ThermoplasticHardnessSource<
    'a,
    S: ThermoplasticSourceFamily<'a> = ThermoplasticSource<'a>,
> {
    source: IndentationHardnessStatePoint,
    reference: ResolvedJ2Plasticity,
    caloric_identity: ContentHash,
    family: PhantomData<fn() -> CaloricPoint<'a, S>>,
}

/// Unmodified native feedback report and its separately sourced hardness.
#[derive(Debug)]
pub struct ThermoplasticHardnessStep<
    'a,
    R,
    S: ThermoplasticSourceFamily<'a> = ThermoplasticSource<'a>,
> {
    pub feedback: FeedbackStep<'a, R, S>,
    pub hardness: IndentationHardnessStatePoint,
    /// Signed response change [Pa], never an additional heat or storage term.
    pub hardness_change_pa: f64,
}

impl<'a, S> ThermoplasticHardnessSource<'a, S>
where
    S: ThermoplasticSourceFamily<'a>,
    <S::Material as SmallStrainLaw>::State: AccumulatedPlasticHistory,
{
    /// Pin an exact original indentation test at an accepted native point.
    /// The original card, caloric chart, selected mechanical source bundle,
    /// fixed environment and hardening convention remain binding thereafter.
    /// This is a source observer shared by compatible specimens, not a unique
    /// specimen-history identity or a registration of cyclic branch ordering.
    pub fn try_new(
        point: &CaloricPoint<'a, S>,
        hardness_key: PropertyKey,
        selection: MaterialPropertySelection,
    ) -> Result<Self, ThermoplasticHardnessError> {
        let query = hardness_query(point)?;
        let thermal = point.thermal_source().thermal_material_point();
        let requirement = ScalarPropertyRequirement::try_with_key(
            &hardness_key,
            ScalarAdmissibility::StrictlyPositive,
        )?;
        let source = IndentationHardnessStatePoint::try_from_resolved(
            resolve_material_state_point(thermal.material, &query, &[requirement], selection)?,
            hardness_key,
        )?;
        admit(thermal, &query, &source)?;
        Ok(Self {
            source,
            reference: point.thermal_source().reference().clone(),
            caloric_identity: point.caloric_state().phase_curve_identity(),
            family: PhantomData,
        })
    }

    /// Original selected hardness, evidence, test convention and claim receipt.
    #[must_use]
    pub const fn hardness_source(&self) -> &IndentationHardnessStatePoint {
        &self.source
    }

    /// Observe H at actual accepted caloric T and native accumulated alpha.
    /// This reuses the original claim and refuses changed material/source
    /// owners, missing endpoint coverage, unowned axes or nonpositive values.
    pub fn at_accepted(
        &self,
        point: &CaloricPoint<'a, S>,
    ) -> Result<IndentationHardnessStatePoint, ThermoplasticHardnessError> {
        if point.thermal_source().reference() != &self.reference
            || point.caloric_state().phase_curve_identity() != self.caloric_identity
        {
            return Err(ThermoplasticHardnessError::Invalid(
                "accepted point belongs to a different original mechanical or caloric source",
            ));
        }
        let query = hardness_query(point)?;
        let thermal = point.thermal_source().thermal_material_point();
        let source = self.source.at_query_point(thermal.material, &query)?;
        admit(thermal, &query, &source)?;
        Ok(source)
    }

    /// Propose one native interval, then publish the original report and typed
    /// hardness together. Every native temperature/work/yield gate is unchanged.
    /// Hardness is queried only at old and accepted endpoints, not within native
    /// fixed-point iterations. No source objects or native histories are copied
    /// beyond the original feedback candidate and the two endpoint receipts.
    pub fn advance<R, E>(
        &self,
        point: &CaloricPoint<'a, S>,
        strain: Voigt,
        duration_s: f64,
        config: FeedbackConfig,
        transport: impl FnMut(UniformEnthalpyStepInput<'_>) -> Result<UniformEnthalpyStep<R>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ThermoplasticHardnessStep<'a, R, S>, ThermoplasticHardnessStepError<E>> {
        if cancelled() {
            return Err(ThermoplasticHardnessStepError::Feedback(
                FeedbackError::Cancelled,
            ));
        }
        let previous = self
            .at_accepted(point)
            .map_err(ThermoplasticHardnessStepError::Hardness)?;
        let feedback = point
            .advance(strain, duration_s, config, transport, &mut cancelled)
            .map_err(ThermoplasticHardnessStepError::Feedback)?;
        let hardness = self
            .at_accepted(&feedback.next)
            .map_err(ThermoplasticHardnessStepError::Hardness)?;
        if cancelled() {
            return Err(ThermoplasticHardnessStepError::Feedback(
                FeedbackError::Cancelled,
            ));
        }
        let hardness_change_pa = hardness.hardness_pa() - previous.hardness_pa();
        Ok(ThermoplasticHardnessStep {
            feedback,
            hardness,
            hardness_change_pa,
        })
    }
}

fn hardness_query<'a, S>(
    point: &CaloricPoint<'a, S>,
) -> Result<QueryPoint, ThermoplasticHardnessError>
where
    S: ThermoplasticSourceFamily<'a>,
    <S::Material as SmallStrainLaw>::State: AccumulatedPlasticHistory,
{
    let alpha = point.history().accumulated_equivalent_plastic_strain();
    if !alpha.is_finite() || alpha < 0.0 {
        return Err(ThermoplasticHardnessError::Invalid(
            "native accumulated equivalent plastic strain must be finite and nonnegative",
        ));
    }
    let thermal = point.thermal_source().thermal_material_point();
    let query = thermal.at(*point.caloric_state())?;
    if point.resolved().material().query() != query {
        return Err(ThermoplasticHardnessError::Invalid(
            "accepted mechanics and caloric state must have the same original context",
        ));
    }
    Ok(query.with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), alpha)?)
}

fn admit(
    thermal: ThermalMaterialPoint<'_>,
    query: &QueryPoint,
    source: &IndentationHardnessStatePoint,
) -> Result<(), ThermoplasticHardnessError> {
    if !source
        .property()
        .has_only_owned_and_fixed_external_conditions(ALPHA, QuantitySpec::dimensional(Dims::NONE))
    {
        return Err(ThermoplasticHardnessError::Invalid(
            "hardness response and applicability must use only native history, temperature and fixed external conditions",
        ));
    }
    ThermalMaterialPoint {
        context: query,
        ..thermal
    }
    .admit_static_with_axes(
        source.property(),
        &[(ALPHA, QuantitySpec::dimensional(Dims::NONE))],
    )?;
    Ok(())
}

/// Original source/test, query, actual-state or endpoint-support refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermoplasticHardnessError {
    Material(MaterialStatePointError),
    ThermalMaterial(ThermalMaterialError),
    Query(MatDbError),
    Invalid(&'static str),
}
impl core::fmt::Display for ThermoplasticHardnessError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => error.fmt(f),
            Self::ThermalMaterial(error) => error.fmt(f),
            Self::Query(error) => error.fmt(f),
            Self::Invalid(reason) => write!(f, "thermoplastic hardness: {reason}"),
        }
    }
}
impl std::error::Error for ThermoplasticHardnessError {}
impl From<MaterialStatePointError> for ThermoplasticHardnessError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
impl From<ThermalMaterialError> for ThermoplasticHardnessError {
    fn from(error: ThermalMaterialError) -> Self {
        Self::ThermalMaterial(error)
    }
}
impl From<MatDbError> for ThermoplasticHardnessError {
    fn from(error: MatDbError) -> Self {
        Self::Query(error)
    }
}

/// Preserve the native feedback refusal separately from a hardness refusal.
#[derive(Debug)]
pub enum ThermoplasticHardnessStepError<E = Infallible> {
    Feedback(FeedbackError<E>),
    Hardness(ThermoplasticHardnessError),
}
impl<E: core::fmt::Display> core::fmt::Display for ThermoplasticHardnessStepError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Feedback(error) => error.fmt(f),
            Self::Hardness(error) => error.fmt(f),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ThermoplasticHardnessStepError<E> {}
