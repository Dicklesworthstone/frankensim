//! Source-paired sensible-enthalpy storage and temperature-dependent transport.
//!
//! One immutable card supplies Cp(T), rho(T) and k(T) at the same fixed external
//! conditions. The existing caloric owner integrates Cp; the existing conduction
//! table evaluates k. Reference density is resolved ONCE at an explicit reference
//! state and never replaced by the chart's changing equilibrium density.
//! This is stationary single-phase transport, without expansion or pressure work.

use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MatDbError, MaterialAnswer, MaterialCard, PropertyKey, QueryPoint};
use fs_material::phase::{
    HeatCapacitySourceQuery, PhaseStateError, ResolvedHeatCapacityCurve, SolidLiquidPhase,
};

use super::heterogeneous::ReferenceEnthalpyMaterial;
use super::{EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyError};
use crate::{ConductionError, ConductionMesh, ConductivityModel, ConductivityTable};
use fs_exec::Cx;

/// Explicit source selection and caloric datum for one stationary material.
#[derive(Clone, Copy, Debug)]
pub struct MaterialEnthalpySource<'a> {
    /// Density/Cp keys and complete lower/upper queries; all external axes stay fixed.
    pub caloric: HeatCapacitySourceQuery<'a>,
    /// Scalar conductivity in W/(m K), from the SAME material card.
    /// The constructor explicitly selects an isotropic conduction model.
    pub conductivity: &'a PropertyKey,
    /// Density is frozen at this point; only its temperature may differ from the endpoints.
    pub reference: &'a QueryPoint,
    /// Explicit single phase. Mixed-phase behavior is not inferred from Cp data.
    pub phase: SolidLiquidPhase,
    /// Declared specific enthalpy at the lower temperature endpoint, J/kg.
    pub lower_specific_enthalpy_j_kg: f64,
    /// Maximum caloric-chart interpolation error, J/kg; not source uncertainty.
    pub maximum_interpolation_error_j_kg: f64,
}

/// A source refusal produces no paired thermal material.
#[derive(Debug)]
pub enum MaterialEnthalpyError {
    /// Invalid or inconsistent reference state.
    Invalid(&'static str),
    /// Original density/Cp source or caloric integration refusal.
    Caloric(PhaseStateError),
    /// Original reference-density query refusal.
    Reference(MatDbError),
    /// Original conductivity source refusal.
    Conductivity(ConductionError),
}

impl std::fmt::Display for MaterialEnthalpyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "material enthalpy transport: {what}"),
            Self::Caloric(error) => error.fmt(f),
            Self::Reference(error) => error.fmt(f),
            Self::Conductivity(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for MaterialEnthalpyError {}

/// Source-bound material consumed by the original uniform or heterogeneous solver.
///
/// Every coefficient retains its original source usage; no nominal Cp or density
/// fills a missing property. Different keys/conditions/cards change the retained
/// identity. The identity does not assert propagated material uncertainty.
#[derive(Clone, Debug)]
pub struct MaterialEnthalpyTransport {
    card_identity: ContentHash,
    caloric: ResolvedHeatCapacityCurve,
    conductivity: ConductivityModel,
    reference_query: QueryPoint,
    reference_density: MaterialAnswer,
    identity: ContentHash,
}

impl MaterialEnthalpyTransport {
    /// Resolve the full interval for an explicitly isotropic conduction model.
    ///
    /// Selecting this constructor declares that the scalar source k(T) describes
    /// isotropic transport. Directional samples require a directional/tensor owner.
    ///
    /// Density, Cp and k must have explicit source temperature support. An
    /// unknown test-temperature flag never creates a heating interval. Existing
    /// owners preserve source knots, typed quantities, fixed external conditions
    /// and continuous support; all source selection stays on this one card.
    pub fn from_isotropic_card(
        card: &MaterialCard,
        source: MaterialEnthalpySource<'_>,
    ) -> Result<Self, MaterialEnthalpyError> {
        let axis = source.caloric.temperature_axis;
        let reference_k =
            source
                .reference
                .axes()
                .get(axis)
                .copied()
                .ok_or(MaterialEnthalpyError::Invalid(
                    "reference state lacks the temperature coordinate",
                ))?;
        let low_k = source.caloric.lower.axes().get(axis).copied().ok_or(
            MaterialEnthalpyError::Invalid("lower state lacks the temperature coordinate"),
        )?;
        let high_k = source.caloric.upper.axes().get(axis).copied().ok_or(
            MaterialEnthalpyError::Invalid("upper state lacks the temperature coordinate"),
        )?;
        let expected_reference = source
            .caloric
            .lower
            .clone()
            .with(axis, reference_k)
            .map_err(MaterialEnthalpyError::Reference)?;
        if !reference_k.is_finite()
            || reference_k <= 0.0
            || reference_k < low_k
            || reference_k > high_k
            || source.reference.axes() != expected_reference.axes()
            || source.reference.axis_quantities() != expected_reference.axis_quantities()
        {
            return Err(MaterialEnthalpyError::Invalid(
                "reference temperature must be supported and all other typed conditions unchanged",
            ));
        }
        let caloric = ResolvedHeatCapacityCurve::try_new(
            card,
            source.caloric,
            source.phase,
            source.lower_specific_enthalpy_j_kg,
            source.maximum_interpolation_error_j_kg,
        )
        .map_err(MaterialEnthalpyError::Caloric)?;
        let table = ConductivityTable::from_claims_at_query_points(
            card.claims(),
            source.conductivity,
            axis,
            &[source.caloric.lower.clone(), source.caloric.upper.clone()],
            source.caloric.selection,
        )
        .map_err(MaterialEnthalpyError::Conductivity)?;
        for receipt in caloric.receipts().iter().chain(table.receipts()) {
            if !card
                .claims()
                .claim(receipt.selected)
                .is_some_and(|claim| claim.validity.bound(axis).is_some())
            {
                return Err(MaterialEnthalpyError::Invalid(
                    "each thermal property requires explicit source temperature support",
                ));
            }
        }
        let reference_density = card
            .claims()
            .query_typed(
                source.caloric.density,
                source.reference,
                source.caloric.selection,
            )
            .map_err(MaterialEnthalpyError::Reference)?;
        let rho0 = reference_density.evidence.value.value;
        // The caloric owner records Cp then density at each retained source knot.
        if !rho0.is_finite()
            || rho0 <= 0.0
            || caloric
                .receipts()
                .iter()
                .skip(1)
                .step_by(2)
                .any(|receipt| receipt.selected != reference_density.receipt.selected)
        {
            return Err(MaterialEnthalpyError::Invalid(
                "reference density must be positive and use the caloric chart's exact density claim",
            ));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-conduction.material-enthalpy.v1");
        hash.update(card.content_hash().as_bytes());
        hash.update(caloric.curve().identity().as_bytes());
        hash.update(reference_density.receipt.content_hash().as_bytes());
        for receipt in table.receipts() {
            hash.update(receipt.content_hash().as_bytes());
        }
        Ok(Self {
            card_identity: card.content_hash(),
            caloric,
            conductivity: ConductivityModel::isotropic(table),
            reference_query: source.reference.clone(),
            reference_density,
            identity: hash.finalize(),
        })
    }

    /// Exact source-card identity shared by storage and transport.
    #[must_use]
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }

    /// Paired caloric, conductivity and reference-density identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Integrated Cp chart with all original Cp/density receipts.
    #[must_use]
    pub const fn caloric(&self) -> &ResolvedHeatCapacityCurve {
        &self.caloric
    }

    /// Source conductivity for the existing nonlinear conduction operator.
    #[must_use]
    pub const fn conductivity(&self) -> &ConductivityModel {
        &self.conductivity
    }

    /// Exact reference conditions, distinct from the evolving temperature.
    #[must_use]
    pub const fn reference_query(&self) -> &QueryPoint {
        &self.reference_query
    }

    /// Original reference-density answer, including uncertainty and usage receipt.
    #[must_use]
    pub const fn reference_density(&self) -> &MaterialAnswer {
        &self.reference_density
    }

    /// Frozen reference mass per unit reference volume, kg/m3.
    #[must_use]
    pub fn reference_density_kg_m3(&self) -> f64 {
        self.reference_density.evidence.value.value
    }

    /// Borrow the original chart and frozen density for heterogeneous assembly.
    #[must_use]
    pub fn reference_material(&self) -> ReferenceEnthalpyMaterial<'_> {
        ReferenceEnthalpyMaterial {
            curve: self.caloric.curve(),
            reference_density_kg_m3: self.reference_density_kg_m3(),
        }
    }

    /// Prepare the existing uniform spatial solver using this paired source.
    pub fn prepare_uniform<'m, 's>(
        &'s self,
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        budget: EnthalpyBudget,
    ) -> Result<EnthalpyBackwardEuler<'m, 's>, EnthalpyError> {
        EnthalpyBackwardEuler::uniform(
            cx,
            mesh,
            self.caloric.curve(),
            self.reference_density_kg_m3(),
            budget,
        )
    }
}
