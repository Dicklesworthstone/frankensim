//! Original indentation hardness on complete native cyclic-plastic histories.
//!
//! Voce saturation and each Chaboche backstress remain owned by their existing
//! constitutive laws. This observer admits their actual endpoint history and
//! nonlinear yield surface, then samples the independent original H(alpha)
//! source prepared by `plastic_history`. It neither projects cyclic history
//! onto a linear J2 yield surface nor repeats a return map. Dynamic recovery of
//! a backstress does not erase accumulated plastic strain or imply annealing.
//!
//! Temperature and other external source coordinates stay fixed. The complete
//! declared alpha interval has the original source's continuous support; every
//! scalar observation is allocation-free. A typed endpoint is resolved only
//! when a consumer requests its original test and claim receipts. No hardness
//! change supplies stress, stored energy, heat, damage or a fatigue-life law.

use fs_matdb::{MaterialCard, PropertyKey, PropertyLineSlice};

use super::IndentationHardnessStatePoint;
use super::plastic_history::{
    ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS, J2HardnessError, ResolvedJ2Hardness,
};
use crate::plastic::chaboche::source::{ResolvedChabocheBranches, ResolvedChabocheJ2};
use crate::plastic::chaboche::{ChabocheError, ChabocheJ2, ChabocheState};
use crate::plastic::energy::EnergyError;
use crate::plastic::source::ResolvedJ2Plasticity;
use crate::plastic::voce::source::ResolvedVoceJ2;
use crate::state_point::MaterialPropertySelection;
use crate::{SmallStrainLaw, Voigt};

mod sealed {
    pub trait Sealed {}
}

/// Existing source-bound rate-independent cyclic laws with complete histories.
/// Sealing prevents an adapter from dropping branches or silently substituting
/// a different accumulated-strain convention or instantaneous yield surface.
pub trait SourceBoundHardnessLaw: SmallStrainLaw + sealed::Sealed {
    /// Full original source bundle carried by the native elastic/J2 backbone.
    /// This is used only for fixed-condition hardness-source preparation.
    fn hardness_backbone(&self) -> &ResolvedJ2Plasticity;

    /// Exact source, native coefficients, hardening basis and branch-role order.
    fn same_hardness_source(&self, other: &Self) -> bool;

    /// Validate the entire native endpoint and return its accumulated strain.
    /// The explicit positive numerical yield allowance is in pascals.
    fn accepted_hardness_alpha(
        &self,
        strain: &Voigt,
        state: &Self::State,
        yield_tolerance_pa: f64,
    ) -> Result<f64, CyclicHardnessError>;
}

impl sealed::Sealed for ResolvedVoceJ2 {}
impl SourceBoundHardnessLaw for ResolvedVoceJ2 {
    fn hardness_backbone(&self) -> &ResolvedJ2Plasticity {
        self.backbone()
    }

    fn same_hardness_source(&self, other: &Self) -> bool {
        same_backbone(self.backbone(), other.backbone()) && self.law() == other.law()
    }

    fn accepted_hardness_alpha(
        &self,
        strain: &Voigt,
        state: &Self::State,
        yield_tolerance_pa: f64,
    ) -> Result<f64, CyclicHardnessError> {
        self.law()
            .validate_committed_history(strain, state, yield_tolerance_pa)?;
        Ok(state.alpha)
    }
}

impl sealed::Sealed for ResolvedChabocheJ2 {}
impl SourceBoundHardnessLaw for ResolvedChabocheJ2 {
    fn hardness_backbone(&self) -> &ResolvedJ2Plasticity {
        self.backbone().backbone()
    }

    fn same_hardness_source(&self, other: &Self) -> bool {
        self.backbone().same_hardness_source(other.backbone()) && self.law() == other.law()
    }

    fn accepted_hardness_alpha(
        &self,
        strain: &Voigt,
        state: &Self::State,
        yield_tolerance_pa: f64,
    ) -> Result<f64, CyclicHardnessError> {
        accepted_chaboche_alpha(self.law(), strain, state, yield_tolerance_pa)
    }
}

impl<const N: usize> sealed::Sealed for ResolvedChabocheBranches<N> {}
impl<const N: usize> SourceBoundHardnessLaw for ResolvedChabocheBranches<N> {
    fn hardness_backbone(&self) -> &ResolvedJ2Plasticity {
        self.backbone().backbone()
    }

    fn same_hardness_source(&self, other: &Self) -> bool {
        self.backbone().same_hardness_source(other.backbone())
            && self.law() == other.law()
            && self.properties() == other.properties()
    }

    fn accepted_hardness_alpha(
        &self,
        strain: &Voigt,
        state: &Self::State,
        yield_tolerance_pa: f64,
    ) -> Result<f64, CyclicHardnessError> {
        accepted_chaboche_alpha(self.law(), strain, state, yield_tolerance_pa)
    }
}

fn same_backbone(a: &ResolvedJ2Plasticity, b: &ResolvedJ2Plasticity) -> bool {
    a.material().identity() == b.material().identity()
        && a.hardening_basis() == b.hardening_basis()
        && a.law() == b.law()
}

fn accepted_chaboche_alpha<const N: usize>(
    law: &ChabocheJ2<N>,
    strain: &Voigt,
    state: &ChabocheState<N>,
    yield_tolerance_pa: f64,
) -> Result<f64, CyclicHardnessError> {
    if !yield_tolerance_pa.is_finite() || yield_tolerance_pa <= 0.0 {
        return Err(EnergyError::Invalid("yield allowance must be finite and positive").into());
    }
    // The native owner checks every branch's deviatoric, accumulated-strain and
    // saturation bounds. Its incoming-history gate deliberately allows trial
    // overstress, so an accepted endpoint also needs its actual cyclic yield.
    law.validate_history(strain, state)?;
    let elastic = core::array::from_fn(|i| strain[i] - state.plastic.plastic_strain[i]);
    let stress = law.backbone().backbone().elastic.stress(&elastic, &());
    let residual = law.yield_function_with_state(&stress, state);
    if !residual.is_finite() {
        return Err(EnergyError::Invalid("committed Chaboche stress is unrepresentable").into());
    }
    if residual > yield_tolerance_pa {
        return Err(EnergyError::YieldResidual {
            residual_pa: residual,
            tolerance_pa: yield_tolerance_pa,
        }
        .into());
    }
    Ok(state.plastic.alpha)
}

/// Complete native cyclic law and its separately sourced indentation response.
/// Mechanical and hardness source preparation is amortized across all points
/// sharing this immutable owner. No source receipt is cloned during observation.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedCyclicHardness<L: SourceBoundHardnessLaw> {
    plastic: L,
    prepared: ResolvedJ2Hardness,
}

impl<L: SourceBoundHardnessLaw> ResolvedCyclicHardness<L> {
    /// Prepare the full virgin-to-maximum interval at fixed native conditions.
    /// Independent measured hardness may increase or decrease with alpha; the
    /// native saturation, backstress and hardening moduli never replace it.
    pub fn new(
        plastic: L,
        card: &MaterialCard,
        hardness_key: PropertyKey,
        selection: MaterialPropertySelection,
        maximum_alpha: f64,
    ) -> Result<Self, CyclicHardnessError> {
        let prepared = ResolvedJ2Hardness::new(
            plastic.hardness_backbone().clone(),
            card,
            hardness_key,
            selection,
            maximum_alpha,
        )?;
        Ok(Self { plastic, prepared })
    }

    /// Original full cyclic law, source conditions and ordered branch roles.
    #[must_use]
    pub const fn plastic(&self) -> &L {
        &self.plastic
    }

    /// Original alpha=0 hardness with the exact indentation test and claim.
    #[must_use]
    pub fn hardness_source(&self) -> &IndentationHardnessStatePoint {
        self.prepared.hardness_source()
    }

    /// Every retained original-source breakpoint and its canonical receipt.
    #[must_use]
    pub fn source_line(&self) -> &PropertyLineSlice {
        self.prepared.source_line()
    }

    /// Upper bound of the complete supported accumulated-strain interval.
    #[must_use]
    pub fn maximum_accumulated_plastic_strain(&self) -> f64 {
        self.prepared.maximum_accumulated_plastic_strain()
    }

    /// Observe an accepted native endpoint without a second return map.
    /// Full source/branch identity, native history, instantaneous cyclic yield
    /// and the prepared support interval must all pass before a scalar returns.
    pub fn hardness_pa_at_accepted(
        &self,
        plastic: &L,
        strain: &Voigt,
        state: &L::State,
        yield_tolerance_pa: f64,
    ) -> Result<f64, CyclicHardnessError> {
        let alpha = self.accepted_alpha(plastic, strain, state, yield_tolerance_pa)?;
        Ok(self.prepared.hardness_pa_at(alpha)?)
    }

    /// Resolve one original typed endpoint for a contact/wear consumer.
    /// This preserves the prepared support restriction even if the card has a
    /// wider interval. The original immutable card and pinned claim are checked.
    pub fn resolve_hardness_at_accepted(
        &self,
        card: &MaterialCard,
        plastic: &L,
        strain: &Voigt,
        state: &L::State,
        yield_tolerance_pa: f64,
    ) -> Result<IndentationHardnessStatePoint, CyclicHardnessError> {
        let alpha = self.accepted_alpha(plastic, strain, state, yield_tolerance_pa)?;
        self.prepared.hardness_pa_at(alpha)?;
        let point = self
            .hardness_source()
            .resolved()
            .query()
            .with(ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS, alpha)
            .map_err(J2HardnessError::from)?;
        self.hardness_source()
            .at_query_point(card, &point)
            .map_err(J2HardnessError::from)
            .map_err(Into::into)
    }

    fn accepted_alpha(
        &self,
        plastic: &L,
        strain: &Voigt,
        state: &L::State,
        yield_tolerance_pa: f64,
    ) -> Result<f64, CyclicHardnessError> {
        if !self.plastic.same_hardness_source(plastic) {
            return Err(J2HardnessError::Invalid(
                "accepted history must retain its complete native cyclic source, law and branch order",
            )
            .into());
        }
        plastic.accepted_hardness_alpha(strain, state, yield_tolerance_pa)
    }
}

/// Original hardness/source refusal or complete native cyclic-history refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum CyclicHardnessError {
    Hardness(J2HardnessError),
    Cyclic(ChabocheError),
}

impl core::fmt::Display for CyclicHardnessError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Hardness(error) => write!(f, "cyclic hardness: {error}"),
            Self::Cyclic(error) => write!(f, "cyclic hardness history: {error}"),
        }
    }
}
impl std::error::Error for CyclicHardnessError {}
impl From<J2HardnessError> for CyclicHardnessError {
    fn from(error: J2HardnessError) -> Self {
        Self::Hardness(error)
    }
}
impl From<EnergyError> for CyclicHardnessError {
    fn from(error: EnergyError) -> Self {
        Self::Hardness(J2HardnessError::Mechanical(error))
    }
}
impl From<ChabocheError> for CyclicHardnessError {
    fn from(error: ChabocheError) -> Self {
        Self::Cyclic(error)
    }
}
