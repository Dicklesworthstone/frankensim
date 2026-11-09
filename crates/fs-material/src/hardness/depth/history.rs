//! Original depth profiles on independently accepted native plastic history.
//!
//! A complete J2, Voce or Chaboche endpoint supplies its actual accumulated
//! equivalent plastic strain. The separately selected indentation source then
//! supplies H(z, alpha) over a continuously supported exposed-depth interval.
//! The accepted alpha is frozen throughout this profile: it must represent the
//! preconditioned material column, not the surface contact's indentation depth.
//! This transfers a property only, without residual stress, plastic work, heat,
//! a second return map, or an assumed relation between hardness and yield.

use fs_blake3::ContentHash;
use fs_matdb::{MaterialCard, PropertyKey};
use fs_qty::{Dims, QuantitySpec};

use super::{DepthHardnessError, EXPOSED_MATERIAL_DEPTH_AXIS, ResolvedDepthHardness, require};
use crate::Voigt;
use crate::hardness::IndentationHardnessStatePoint;
use crate::hardness::cyclic_history::{
    CyclicHardnessError, SourceBoundHardnessLaw, native_source_identity,
};
use crate::hardness::plastic_history::ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA;
use crate::plastic::J2State;
use crate::plastic::energy::EnergyError;
use crate::plastic::source::ResolvedJ2Plasticity;
use crate::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, is_fixed_external_coordinate, resolve_material_state_point,
};

/// Immutable provenance from complete native endpoint admission. This receipt
/// retains the original constitutive source and actual alpha without copying
/// tensor histories into each later contact/wear proposal. It does not certify
/// an unobserved loading trajectory or spatial uniformity of a supplied column.
#[derive(Clone, Debug, PartialEq)]
pub struct AcceptedDepthHardnessHistory {
    mechanical: ResolvedMaterialStatePoint,
    native_identity: ContentHash,
    alpha: f64,
}

impl AcceptedDepthHardnessHistory {
    /// Complete original mechanical claims, fixed context and source receipts.
    #[must_use]
    pub const fn mechanical_source(&self) -> &ResolvedMaterialStatePoint {
        &self.mechanical
    }

    /// Original law, hardening convention and complete cyclic branch ordering.
    #[must_use]
    pub const fn native_source_identity(&self) -> ContentHash {
        self.native_identity
    }

    /// Actual admitted accumulated equivalent plastic strain, dimensionless.
    #[must_use]
    pub const fn alpha(&self) -> f64 {
        self.alpha
    }
}

/// No profile or partial state is published when native history or original
/// continuous indentation-source coverage refuses.
#[derive(Clone, Debug, PartialEq)]
pub enum DepthHardnessHistoryError {
    /// Complete native J2 tensor history or accepted yield surface was refused.
    NativeJ2(EnergyError),
    /// Complete native cyclic history, including every branch, was refused.
    NativeCyclic(CyclicHardnessError),
    /// Original source selection, continuous depth support or positivity failed.
    Profile(DepthHardnessError),
}
impl core::fmt::Display for DepthHardnessHistoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "accepted-history depth hardness: {self:?}")
    }
}
impl std::error::Error for DepthHardnessHistoryError {}
impl From<DepthHardnessError> for DepthHardnessHistoryError {
    fn from(error: DepthHardnessError) -> Self {
        Self::Profile(error)
    }
}

impl ResolvedDepthHardness {
    /// Prepare original H(z, alpha) from one complete accepted native J2 state.
    ///
    /// The original owner validates tensor history and its instantaneous yield
    /// surface before alpha is read. Every depth retains the same actual alpha,
    /// exact original indentation test and fixed mechanical environment. The
    /// caller owns the assumption of uniform preconditioning through the column.
    /// No contact indentation, removal depth or caller-provided scalar alpha is
    /// substituted for native accumulated strain. All inputs remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn from_j2_accepted(
        plastic: &ResolvedJ2Plasticity,
        strain: &Voigt,
        state: &J2State,
        yield_tolerance_pa: f64,
        card: &MaterialCard,
        hardness_key: PropertyKey,
        selection: MaterialPropertySelection,
        maximum_depth_m: f64,
    ) -> Result<Self, DepthHardnessHistoryError> {
        plastic
            .law()
            .validate_committed_history(strain, state, yield_tolerance_pa)
            .map_err(DepthHardnessHistoryError::NativeJ2)?;
        Ok(prepare(
            plastic,
            native_source_identity(b"j2", plastic, &[]),
            state.alpha,
            card,
            hardness_key,
            selection,
            maximum_depth_m,
        )?)
    }

    /// Prepare the same original depth profile from full native Voce/Chaboche
    /// history. Every ordered tensor backstress and the nonlinear yield surface
    /// pass their existing owner; no bare-J2 projection is used for admission.
    #[allow(clippy::too_many_arguments)]
    pub fn from_cyclic_accepted<L: SourceBoundHardnessLaw>(
        plastic: &L,
        strain: &Voigt,
        state: &L::State,
        yield_tolerance_pa: f64,
        card: &MaterialCard,
        hardness_key: PropertyKey,
        selection: MaterialPropertySelection,
        maximum_depth_m: f64,
    ) -> Result<Self, DepthHardnessHistoryError> {
        let alpha = plastic
            .accepted_hardness_alpha(strain, state, yield_tolerance_pa)
            .map_err(DepthHardnessHistoryError::NativeCyclic)?;
        Ok(prepare(
            plastic.hardness_backbone(),
            plastic.hardness_source_identity(),
            alpha,
            card,
            hardness_key,
            selection,
            maximum_depth_m,
        )?)
    }
}

#[allow(clippy::too_many_arguments)]
fn prepare(
    plastic: &ResolvedJ2Plasticity,
    native_identity: ContentHash,
    alpha: f64,
    card: &MaterialCard,
    hardness_key: PropertyKey,
    selection: MaterialPropertySelection,
    maximum_depth_m: f64,
) -> Result<ResolvedDepthHardness, DepthHardnessError> {
    require(
        maximum_depth_m.is_finite() && maximum_depth_m > 0.0,
        "finite positive exposed-material-depth limit required",
    )?;
    let mechanical = plastic.material();
    require(
        card.content_hash() == mechanical.card_identity() && card.id() == mechanical.material(),
        "native history and depth hardness require the same original material card",
    )?;
    for (axis, value) in mechanical.query_point() {
        require(
            is_fixed_external_coordinate(
                axis,
                mechanical.axis_quantities().get(axis).copied(),
                Some(*value),
            ),
            "native depth hardness requires a complete fixed external mechanical context",
        )?;
    }
    let query = mechanical
        .query()
        .with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), alpha)?
        .with_quantity(
            EXPOSED_MATERIAL_DEPTH_AXIS,
            QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
            0.0,
        )?;
    let requirement = ScalarPropertyRequirement::try_with_key(
        &hardness_key,
        ScalarAdmissibility::StrictlyPositive,
    )?;
    let source = IndentationHardnessStatePoint::try_from_resolved(
        resolve_material_state_point(card, &query, &[requirement], selection)?,
        hardness_key,
    )?;
    ResolvedDepthHardness::prepare_from_source(
        card,
        source,
        maximum_depth_m,
        Some(AcceptedDepthHardnessHistory {
            mechanical: mechanical.clone(),
            native_identity,
            alpha,
        }),
    )
}
