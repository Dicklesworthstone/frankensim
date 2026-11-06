//! Original indentation hardness evaluated on native accumulated plastic history.
//!
//! The canonical coordinate below is J2's accumulated equivalent plastic strain:
//! alpha = integral sqrt(2/3 d(eps_p):d(eps_p)). It is neither a signed strain
//! component, total strain, strain rate, indentation depth nor contact recession.
//! The independent pressure-valued hardness source retains its original test.
//! Yield stress and either J2 hardening modulus never supply a missing hardness.
//!
//! Preparation uses the original material database's continuous line-support
//! owner. Scalars, curves and connected cuts through sampled surfaces retain
//! their pinned claim, interpolation policy, all crossed source-cell boundaries,
//! and fixed environmental conditions. Holes and exact-sample-only spans refuse.
//! Positive hardness is checked at every retained affine breakpoint.
//!
//! Runtime updates reuse the native checked J2 energy update and one binary
//! search over prepared numeric knots. They allocate no source objects and
//! perform no database query. Hardness is a response of the retained history;
//! changing it adds no mechanical storage, plastic heat, damage or wear.
//! Each candidate is atomic: a source-range or native work refusal publishes
//! neither advanced mechanical history nor a new hardness.
//!
//! Temperature and all other external conditions stay fixed at the original
//! mechanical state. This does not implement annealing, dynamic strain aging,
//! thermal history transfer, size effects beyond the declared test, or a
//! universal hardness-to-yield conversion. No measured-material accuracy or
//! service-life qualification follows from the numerical G0/G1/G3 regressions.

use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{
    ClaimSelection, MatDbError, MaterialCard, PropertyKey, PropertyLineSlice,
    PropertySupportError,
};
use fs_qty::{Dims, QuantitySpec};

use super::IndentationHardnessStatePoint;
use crate::plastic::J2State;
use crate::plastic::energy::{EnergyError, EnergyTolerances, IncrementEnergy};
use crate::plastic::source::{HardeningPropertyBasis, ResolvedJ2Plasticity};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ScalarAdmissibility,
    ScalarPropertyRequirement, is_fixed_external_coordinate, resolve_material_state_point,
};
use crate::{SmallStrainLaw, Voigt};

/// Accumulated J2 equivalent plastic strain; its exact name is part of the
/// source convention. Generic strain coordinates and strain components are
/// not aliases even though they are also dimensionless.
pub const ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS: &str =
    "accumulated-equivalent-plastic-strain";

/// Original native J2 law and independently sourced indentation-hardness path.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedJ2Hardness {
    plastic: ResolvedJ2Plasticity,
    source: IndentationHardnessStatePoint,
    source_line: PropertyLineSlice,
    // Keep runtime searches independent of the size/layout of source receipts.
    knots: Vec<(f64, f64)>,
    maximum_alpha: f64,
    identity: ContentHash,
}

/// Accepted native strain/history and the corresponding indentation hardness.
/// Private fields prevent changing hardness independently of the actual history.
#[derive(Clone, Debug, PartialEq)]
pub struct J2HardnessState {
    owner: ContentHash,
    strain: Voigt,
    plastic: J2State,
    hardness_pa: f64,
}

impl J2HardnessState {
    /// Accepted total small strain in the native tensor-Voigt convention.
    #[must_use]
    pub const fn strain(&self) -> &Voigt {
        &self.strain
    }

    /// Original native J2 history, including irreversible accumulated strain.
    #[must_use]
    pub const fn plastic(&self) -> &J2State {
        &self.plastic
    }

    /// Positive force/area hardness from the separately selected original test.
    #[must_use]
    pub const fn hardness_pa(&self) -> f64 {
        self.hardness_pa
    }
}

/// Checked mechanical increment with its consistent hardness endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct J2HardnessAdvance {
    /// Candidate accepted state; the input state remains untouched.
    pub state: J2HardnessState,
    /// Native endpoint stress, with no hardness-derived correction [Pa].
    pub stress_pa: Voigt,
    /// Original mechanical storage, physical dissipation and numerical loss.
    /// A hardness change is not an additional term in this work account.
    pub energy: IncrementEnergy,
    /// Signed indentation-hardness change [Pa]; data may harden or soften.
    pub hardness_change_pa: f64,
}

impl ResolvedJ2Hardness {
    /// Prepare a virgin-to-maximum-alpha path at the original J2 conditions.
    ///
    /// The original material card must own both mechanical and hardness data.
    /// The caller selects an exact pressure-valued hardness key, including its
    /// test apparatus/protocol/specimen, and an explicit accumulated-strain
    /// limit. The complete interval must have continuous original-source support.
    /// All retained coordinates outside the canonical alpha axis are held fixed.
    ///
    /// Source preparation is bounded by fs-matdb's MAX_PROPERTY_LINE_KNOTS.
    /// Neither a missing virgin datum nor coverage beyond the requested/source
    /// interval is inferred. No history has advanced when preparation refuses.
    pub fn new(
        plastic: ResolvedJ2Plasticity,
        card: &MaterialCard,
        hardness_key: PropertyKey,
        selection: MaterialPropertySelection,
        maximum_alpha: f64,
    ) -> Result<Self, J2HardnessError> {
        require(
            maximum_alpha.is_finite() && maximum_alpha > 0.0,
            "declare a finite positive accumulated-plastic-strain limit",
        )?;
        require(
            card.content_hash() == plastic.material().card_identity()
                && card.id() == plastic.material().material(),
            "hardness and native J2 must retain the same original material card",
        )?;
        for (axis, value) in plastic.material().query_point() {
            require(
                is_fixed_external_coordinate(
                    axis,
                    plastic.material().axis_quantities().get(axis).copied(),
                    Some(*value),
                ),
                "J2 hardness preparation requires a complete fixed external context",
            )?;
        }
        let alpha_quantity = QuantitySpec::dimensional(Dims::NONE);
        let lower = plastic.material().query().with_quantity(
            ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS,
            alpha_quantity,
            0.0,
        )?;
        let requirement = ScalarPropertyRequirement::try_with_key(
            &hardness_key,
            ScalarAdmissibility::StrictlyPositive,
        )?;
        let source = IndentationHardnessStatePoint::try_from_resolved(
            resolve_material_state_point(card, &lower, &[requirement], selection)?,
            hardness_key,
        )?;
        require(
            source.property().has_only_owned_and_fixed_external_conditions(
                ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS,
                alpha_quantity,
            ),
            "hardness must depend only on accumulated J2 strain and fixed external conditions",
        )?;
        let upper = lower.clone().with(
            ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS,
            maximum_alpha,
        )?;
        let source_line = card.claims().query_line_typed(
            source.key(),
            ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS,
            &lower,
            &upper,
            ClaimSelection::Pinned(source.property().answer().receipt.selected),
        )?;
        let knots: Vec<_> = source_line
            .knots
            .iter()
            .map(|(alpha, answer)| (*alpha, answer.evidence.value.value))
            .collect();
        require(
            knots.len() >= 2
                && knots.first().is_some_and(|knot| knot.0 == 0.0)
                && knots.last().is_some_and(|knot| knot.0 == maximum_alpha)
                && knots.iter().all(|(alpha, hardness)| {
                    alpha.is_finite() && hardness.is_finite() && *hardness > 0.0
                })
                && knots.windows(2).all(|pair| pair[0].0 < pair[1].0),
            "every supported source breakpoint needs finite positive indentation hardness",
        )?;
        let mut identity = DomainHasher::new("org.frankensim.fs-material.j2-hardness.v1");
        identity.update(plastic.material().identity().as_bytes());
        identity.update(source.identity().as_bytes());
        identity.update(&maximum_alpha.to_bits().to_le_bytes());
        identity.update(&[match plastic.hardening_basis() {
            HardeningPropertyBasis::PlasticStrainModulus => 0,
            HardeningPropertyBasis::MixedPlasticStrainModuli => 1,
            HardeningPropertyBasis::UniaxialTangentModulus => 2,
        }]);
        let law = plastic.law();
        for value in [
            law.elastic.youngs,
            law.elastic.poisson,
            law.elastic.strain_limit,
            law.yield_stress,
            law.hardening,
            law.kinematic_hardening,
        ] {
            identity.update(&value.to_bits().to_le_bytes());
        }
        Ok(Self {
            plastic,
            source,
            source_line,
            knots,
            maximum_alpha,
            identity: identity.finalize(),
        })
    }

    /// Original mechanical source and native constitutive owner.
    #[must_use]
    pub const fn plastic(&self) -> &ResolvedJ2Plasticity {
        &self.plastic
    }

    /// Original hardness state at alpha=0, including the complete test context.
    #[must_use]
    pub const fn hardness_source(&self) -> &IndentationHardnessStatePoint {
        &self.source
    }

    /// Every original continuous-support breakpoint and its canonical receipt.
    #[must_use]
    pub const fn source_line(&self) -> &PropertyLineSlice {
        &self.source_line
    }

    /// Maximum admitted accumulated equivalent plastic strain.
    #[must_use]
    pub const fn maximum_accumulated_plastic_strain(&self) -> f64 {
        self.maximum_alpha
    }

    /// Native virgin history with the original source's alpha=0 response.
    #[must_use]
    pub fn initial_state(&self) -> J2HardnessState {
        J2HardnessState {
            owner: self.identity,
            strain: [0.0; 6],
            plastic: self.plastic.initial_state(),
            hardness_pa: self.knots[0].1,
        }
    }

    /// Evaluate the prepared source path without evolving physical history.
    /// Original knots are exact; interiors use the original affine source law.
    /// This operation is bounded, allocation-free and refuses extrapolation.
    pub fn hardness_pa_at(&self, alpha: f64) -> Result<f64, J2HardnessError> {
        require(
            alpha.is_finite() && alpha >= 0.0 && alpha <= self.maximum_alpha,
            "accumulated plastic strain exceeds the prepared hardness-source interval",
        )?;
        let right = self.knots.partition_point(|knot| knot.0 < alpha);
        let (x1, y1) = self.knots[right];
        if x1 == alpha {
            return Ok(y1);
        }
        let (x0, y0) = self.knots[right - 1];
        let value = y0 + ((alpha - x0) / (x1 - x0)) * (y1 - y0);
        require(
            value.is_finite() && value > 0.0,
            "interpolated indentation hardness is not representable and positive",
        )?;
        Ok(value)
    }

    /// Propose one physical loading/unloading/reversal step through the native
    /// J2 return map and its existing energy/yield gates. Native accumulated
    /// strain alone drives hardness; absolute total strain is never substituted.
    ///
    /// The original stress and energetic-hardening interpretation are retained.
    /// In particular, stored hardening and backward-Euler numerical loss do not
    /// become heat when hardness changes. No source objects are cloned here.
    pub fn advance_candidate(
        &self,
        state: &J2HardnessState,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<J2HardnessAdvance, J2HardnessError> {
        self.require_state(state)?;
        let native = self.plastic.law().advance_with_energy(
            &state.strain,
            &state.plastic,
            strain,
            tolerances,
        )?;
        let hardness_pa = self.hardness_pa_at(native.state.alpha)?;
        Ok(J2HardnessAdvance {
            state: J2HardnessState {
                owner: self.identity,
                strain: native.strain,
                plastic: native.state,
                hardness_pa,
            },
            stress_pa: native.stress_pa,
            energy: native.energy,
            hardness_change_pa: hardness_pa - state.hardness_pa,
        })
    }

    /// Export a complete typed source endpoint only when a consumer needs it.
    /// The ordinary immutable-card/claim/quantity checks remain authoritative;
    /// this requery does not change history or allocate during ordinary stepping.
    pub fn resolve_hardness_state(
        &self,
        card: &MaterialCard,
        state: &J2HardnessState,
    ) -> Result<IndentationHardnessStatePoint, J2HardnessError> {
        self.require_state(state)?;
        let point = self.source.resolved().query().with(
            ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS,
            state.plastic.alpha,
        )?;
        Ok(self.source.at_query_point(card, &point)?)
    }

    fn require_state(&self, state: &J2HardnessState) -> Result<(), J2HardnessError> {
        require(
            state.owner == self.identity,
            "hardness history belongs to a different source, native J2 law or support interval",
        )
    }
}

/// Original-source, continuous-support or native mechanical refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum J2HardnessError {
    /// Original material-property/key admission.
    Material(MaterialStatePointError),
    /// Original physical query construction.
    Query(MatDbError),
    /// Original source's full continuous support or interpolation policy.
    Support(PropertySupportError),
    /// Native return-map, strain, yield or work-account failure.
    Mechanical(EnergyError),
    /// Incompatible source/history or unrepresentable hardness.
    Invalid(&'static str),
}

impl core::fmt::Display for J2HardnessError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "J2 hardness material: {error}"),
            Self::Query(error) => write!(f, "J2 hardness query: {error}"),
            Self::Support(error) => write!(f, "J2 hardness source support: {error}"),
            Self::Mechanical(error) => write!(f, "J2 hardness mechanics: {error}"),
            Self::Invalid(reason) => write!(f, "J2 hardness refused: {reason}"),
        }
    }
}
impl std::error::Error for J2HardnessError {}

impl From<MaterialStatePointError> for J2HardnessError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
impl From<MatDbError> for J2HardnessError {
    fn from(error: MatDbError) -> Self {
        Self::Query(error)
    }
}
impl From<PropertySupportError> for J2HardnessError {
    fn from(error: PropertySupportError) -> Self {
        Self::Support(error)
    }
}
impl From<EnergyError> for J2HardnessError {
    fn from(error: EnergyError) -> Self {
        Self::Mechanical(error)
    }
}

fn require(accepted: bool, reason: &'static str) -> Result<(), J2HardnessError> {
    if accepted {
        Ok(())
    } else {
        Err(J2HardnessError::Invalid(reason))
    }
}
