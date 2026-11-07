//! Prepared original indentation-hardness profiles below a reference surface.
//!
//! The depth coordinate is exposed original MATERIAL depth after removal. It is
//! not indenter penetration, accumulated plastic indentation, or a coating law.
//! An original test-qualified source supplies the profile. Each prepared owner
//! holds temperature and every other admitted external condition fixed;
//! `at_temperature` re-prepares the same original test over its complete depth
//! interval at a new prescribed temperature. Continuous support and positivity
//! are checked including interior breakpoints and holes, with no extrapolation.
//! Runtime depth evaluation is a binary search without source allocation or
//! querying. Temperature refresh adds no heat, annealing kinetics or depth change.

mod integral;
pub use integral::DepthHardnessAdvance;

use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{
    ClaimSelection, MatDbError, MaterialCard, PropertyKey, PropertyLineSlice, PropertySupportError,
    QueryPoint,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};

use super::IndentationHardnessStatePoint;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ScalarAdmissibility,
    ScalarPropertyRequirement, is_fixed_external_coordinate, resolve_material_state_point,
};

/// Original material depth exposed by removal [m]; distinct from penetration.
pub const EXPOSED_MATERIAL_DEPTH_AXIS: &str = "exposed-material-depth";

/// One exact original hardness test, fixed environment and continuous depth line.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedDepthHardness {
    source: IndentationHardnessStatePoint,
    source_line: PropertyLineSlice,
    knots: Vec<(f64, f64)>,
    maximum_depth_m: f64,
    identity: ContentHash,
}

impl ResolvedDepthHardness {
    /// Prepare a positive original hardness profile from zero exposed depth.
    /// Source line preparation is bounded by fs-matdb's native knot budget.
    /// This is a measured/modelled response profile, not inferred case depth,
    /// layered elasticity, a yield conversion or a changing indentation test.
    pub fn new(
        card: &MaterialCard,
        fixed: &QueryPoint,
        hardness_key: PropertyKey,
        selection: MaterialPropertySelection,
        maximum_depth_m: f64,
    ) -> Result<Self, DepthHardnessError> {
        require(
            maximum_depth_m.is_finite() && maximum_depth_m > 0.0,
            "finite positive exposed-material-depth limit required",
        )?;
        for (axis, value) in fixed.axes() {
            require(
                is_fixed_external_coordinate(
                    axis,
                    fixed.axis_quantities().get(axis).copied(),
                    Some(*value),
                ),
                "depth hardness requires otherwise fixed external conditions",
            )?;
        }
        let depth_quantity = QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0]));
        let lower =
            fixed
                .clone()
                .with_quantity(EXPOSED_MATERIAL_DEPTH_AXIS, depth_quantity, 0.0)?;
        let requirement = ScalarPropertyRequirement::try_with_key(
            &hardness_key,
            ScalarAdmissibility::StrictlyPositive,
        )?;
        let source = IndentationHardnessStatePoint::try_from_resolved(
            resolve_material_state_point(card, &lower, &[requirement], selection)?,
            hardness_key,
        )?;
        Self::prepare_from_source(card, source, maximum_depth_m)
    }

    // Both initial preparation and temperature continuation retain the same
    // full-depth support, positivity, original source receipts and identity.
    fn prepare_from_source(
        card: &MaterialCard,
        source: IndentationHardnessStatePoint,
        maximum_depth_m: f64,
    ) -> Result<Self, DepthHardnessError> {
        let depth_quantity = QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0]));
        require(
            source
                .property()
                .has_only_owned_and_fixed_external_conditions(
                    EXPOSED_MATERIAL_DEPTH_AXIS,
                    depth_quantity,
                ),
            "hardness may depend only on exposed material depth and fixed external conditions",
        )?;
        let lower = source.resolved().query();
        let upper = lower
            .clone()
            .with(EXPOSED_MATERIAL_DEPTH_AXIS, maximum_depth_m)?;
        let source_line = card.claims().query_line_typed(
            source.key(),
            EXPOSED_MATERIAL_DEPTH_AXIS,
            &lower,
            &upper,
            ClaimSelection::Pinned(source.property().answer().receipt.selected),
        )?;
        let knots: Vec<_> = source_line
            .knots
            .iter()
            .map(|(depth, answer)| (*depth, answer.evidence.value.value))
            .collect();
        require(
            knots.len() >= 2
                && knots.first().is_some_and(|knot| knot.0 == 0.0)
                && knots.last().is_some_and(|knot| knot.0 == maximum_depth_m)
                && knots.iter().all(|(depth, hardness)| {
                    depth.is_finite() && hardness.is_finite() && *hardness > 0.0
                })
                && knots.windows(2).all(|pair| pair[0].0 < pair[1].0),
            "every supported depth breakpoint needs finite positive indentation hardness",
        )?;
        let mut hash = DomainHasher::new("org.frankensim.fs-material.depth-hardness.v1");
        hash.update(source.identity().as_bytes());
        hash.update(&maximum_depth_m.to_bits().to_le_bytes());
        Ok(Self {
            source,
            source_line,
            knots,
            maximum_depth_m,
            identity: hash.finalize(),
        })
    }

    /// Re-prepare the SAME original depth profile at a prescribed temperature.
    ///
    /// Only the existing canonical absolute `T` coordinate changes. Its exact
    /// untyped, dimensional-temperature or semantic-absolute convention remains
    /// unchanged, as do every other condition, the complete test, pinned claim
    /// and maximum depth. The original card is required even for a no-op.
    ///
    /// Full continuous depth support and positive hardness at every retained
    /// breakpoint must pass again. Positive surface/deep endpoints alone are
    /// insufficient. This checks the new temperature's full depth line, not
    /// support along the intervening thermal path. No thermal history, heat,
    /// annealing, expansion, wear recession or indentation state is advanced.
    pub fn at_temperature(
        &self,
        card: &MaterialCard,
        temperature_k: f64,
    ) -> Result<Self, DepthHardnessError> {
        let original = self.source.resolved().query();
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        require(
            temperature_k.is_finite()
                && temperature_k > 0.0
                && original
                    .axes()
                    .get("T")
                    .is_some_and(|t| t.is_finite() && *t > 0.0)
                && original.axis_quantities().get("T").is_none_or(|quantity| {
                    *quantity == absolute
                        || *quantity == QuantitySpec::dimensional(Temperature::DIMS)
                }),
            "depth hardness refresh requires its original absolute-temperature T coordinate",
        )?;
        require(
            card.content_hash() == self.source.resolved().card_identity()
                && card.id() == self.source.resolved().material(),
            "depth hardness refresh requires the original immutable material card",
        )?;
        // Preserve the exact initial selection receipt and prepared identity.
        // A changed temperature correctly records a pinned-claim requery, but
        // an exact no-op needs neither another source query nor another proof.
        if original.axes()["T"].to_bits() == temperature_k.to_bits() {
            return Ok(self.clone());
        }
        let query = original.with("T", temperature_k)?;
        let source = self.source.at_query_point(card, &query)?;
        Self::prepare_from_source(card, source, self.maximum_depth_m)
    }

    /// Original surface endpoint with the exact selected claim and indentation test.
    #[must_use]
    pub const fn reference_source(&self) -> &IndentationHardnessStatePoint {
        &self.source
    }

    /// Complete original continuous-support breakpoints and receipts.
    #[must_use]
    pub const fn source_line(&self) -> &PropertyLineSlice {
        &self.source_line
    }

    /// Largest supported exposed material depth [m].
    #[must_use]
    pub const fn maximum_depth_m(&self) -> f64 {
        self.maximum_depth_m
    }

    /// Identity of the selected source, fixed environment and admitted depth range.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Original affine line value at an exposed depth, with exact knot endpoints.
    /// No runtime source query, allocation, extrapolation or physical evolution.
    pub fn hardness_pa_at(&self, depth_m: f64) -> Result<f64, DepthHardnessError> {
        require(
            depth_m.is_finite() && depth_m >= 0.0 && depth_m <= self.maximum_depth_m,
            "exposed depth lies outside the prepared original hardness profile",
        )?;
        let right = self.knots.partition_point(|knot| knot.0 < depth_m);
        let (x1, y1) = self.knots[right];
        if x1 == depth_m {
            return Ok(y1);
        }
        let (x0, y0) = self.knots[right - 1];
        let value = y0 + ((depth_m - x0) / (x1 - x0)) * (y1 - y0);
        require(
            value.is_finite() && value > 0.0,
            "positive representable depth hardness required",
        )?;
        Ok(value)
    }

    /// Export a typed original source endpoint only when a consumer needs it.
    /// The original card, selected claim and indentation test remain authoritative.
    pub fn resolve_at_depth(
        &self,
        card: &MaterialCard,
        depth_m: f64,
    ) -> Result<IndentationHardnessStatePoint, DepthHardnessError> {
        self.hardness_pa_at(depth_m)?;
        let query = self
            .source
            .resolved()
            .query()
            .with(EXPOSED_MATERIAL_DEPTH_AXIS, depth_m)?;
        Ok(self.source.at_query_point(card, &query)?)
    }
}

/// Original source, continuous-support, physical-coordinate or arithmetic refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum DepthHardnessError {
    /// Original typed source selection or endpoint admission failed.
    Material(MaterialStatePointError),
    /// Original query coordinates or conventions were refused.
    Query(MatDbError),
    /// Continuous original depth support is absent or exceeds its native budget.
    Support(PropertySupportError),
    /// A physical or arithmetic requirement was not met.
    Invalid(&'static str),
}
impl core::fmt::Display for DepthHardnessError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "exposed-depth indentation hardness: {self:?}")
    }
}
impl std::error::Error for DepthHardnessError {}
impl From<MaterialStatePointError> for DepthHardnessError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
impl From<MatDbError> for DepthHardnessError {
    fn from(error: MatDbError) -> Self {
        Self::Query(error)
    }
}
impl From<PropertySupportError> for DepthHardnessError {
    fn from(error: PropertySupportError) -> Self {
        Self::Support(error)
    }
}
fn require(condition: bool, reason: &'static str) -> Result<(), DepthHardnessError> {
    if condition {
        Ok(())
    } else {
        Err(DepthHardnessError::Invalid(reason))
    }
}
