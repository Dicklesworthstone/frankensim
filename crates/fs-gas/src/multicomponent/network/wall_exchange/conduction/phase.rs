//! Phase-aware conductivity from the same immutable chart as finite wall h.

use super::WallExchangeError;
use crate::multicomponent::wall_exchange::MaterialWallCaloric;
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point,
};
use fs_material::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use std::sync::Arc;

/// Original conductivity source with owned temperature and liquid MASS fraction.
///
/// The retained equilibrium chart determines both coordinates from specific h,
/// including distinct inventories on an isothermal latent plateau. The source
/// supplies k(T,phi); no universal phase interpolation or invented liquid k is
/// inserted. Pressure, composition and other admitted external conditions stay
/// fixed. Exactly one claim is selected at preparation and pinned on refresh.
///
/// This is a stationary material heat-transfer binding. It neither supplies
/// liquid momentum nor determines wall support, shape or moving-front geometry.
/// Source evaluations retain their original point receipts; they establish
/// support at the sampled states, not a continuous-path certificate.
#[derive(Clone, Debug, PartialEq)]
pub struct PhaseWallConductionSource {
    curve: Arc<EquilibriumEnthalpyPhaseCurve>,
    card: Arc<MaterialCard>,
    context: QueryPoint,
    temperature_axis: String,
    liquid_fraction_axis: String,
    reference: ResolvedMaterialStatePoint,
}

impl PhaseWallConductionSource {
    /// Resolve the original static conductivity once at an explicit h anchor.
    /// The anchor need not be the first knot of the caloric curve; source
    /// coverage may be a bounded subdomain. Every actual runtime endpoint is
    /// checked again with the same pinned claim and coordinate meanings.
    pub fn new(
        curve: Arc<EquilibriumEnthalpyPhaseCurve>,
        point: ThermalMaterialPoint<'_>,
        conductivity: &PropertyKey,
        selection: MaterialPropertySelection,
        reference_specific_h_j_kg: f64,
    ) -> Result<Self, WallExchangeError> {
        let quantity = conductivity.quantity();
        if quantity != QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS)
            && quantity
                != QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::ThermalConductivity,
                    ValueForm::Static,
                ))
        {
            return Err(WallExchangeError::Invalid(
                "phase wall requires a static scalar conductivity source",
            ));
        }
        let fraction = point
            .liquid_fraction_axis
            .ok_or(WallExchangeError::Invalid(
                "phase wall conductivity requires an explicit liquid mass fraction coordinate",
            ))?;
        let phase = curve
            .state_at_specific_enthalpy(reference_specific_h_j_kg)
            .map_err(WallExchangeError::PhaseConductionState)?;
        let query = point
            .at(phase)
            .map_err(WallExchangeError::PhaseConduction)?;
        let requirement = ScalarPropertyRequirement::try_with_key(
            conductivity,
            ScalarAdmissibility::StrictlyPositive,
        )
        .map_err(material_error)?;
        let reference =
            resolve_material_state_point(point.material, &query, &[requirement], selection)
                .map_err(material_error)?;
        point
            .admit_static(&reference.properties()[0])
            .map_err(WallExchangeError::PhaseConduction)?;
        Ok(Self {
            curve,
            card: Arc::new(point.material.clone()),
            context: point.context.clone(),
            temperature_axis: point.temperature_axis.to_owned(),
            liquid_fraction_axis: fraction.to_owned(),
            reference,
        })
    }

    /// Original immutable material and all of its measured/declared claims.
    pub fn material(&self) -> &MaterialCard {
        &self.card
    }

    /// Exact equilibrium chart used to derive both T and liquid fraction.
    pub fn curve(&self) -> &Arc<EquilibriumEnthalpyPhaseCurve> {
        &self.curve
    }

    /// Preparation-state source receipt; refresh retains its selected claim.
    pub fn reference(&self) -> &ResolvedMaterialStatePoint {
        &self.reference
    }

    /// Resolve the same conductivity at the actual original enthalpy inventory.
    /// Every returned value retains a complete source receipt. No accepted gas
    /// or wall state is changed by this query.
    pub fn at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
    ) -> Result<ResolvedMaterialStatePoint, WallExchangeError> {
        let phase = self
            .curve
            .state_at_specific_enthalpy(specific_h_j_kg)
            .map_err(WallExchangeError::PhaseConductionState)?;
        let point = ThermalMaterialPoint {
            material: &self.card,
            context: &self.context,
            temperature_axis: &self.temperature_axis,
            liquid_fraction_axis: Some(&self.liquid_fraction_axis),
        };
        let query = point
            .at(phase)
            .map_err(WallExchangeError::PhaseConduction)?;
        self.reference
            .at_query_point(&self.card, &query)
            .map_err(material_error)
    }

    pub(super) fn admit_caloric(
        &self,
        caloric: Option<&MaterialWallCaloric<'_>>,
    ) -> Result<(), WallExchangeError> {
        let chart = caloric.and_then(MaterialWallCaloric::phase_curve).ok_or(
            WallExchangeError::Invalid(
                "phase conduction requires the original phase-chart heat storage binding",
            ),
        )?;
        if chart.identity() != self.curve.identity() {
            return Err(WallExchangeError::Invalid(
                "phase conduction and caloric storage must retain the same original chart",
            ));
        }
        Ok(())
    }
}

fn material_error(error: fs_material::state_point::MaterialStatePointError) -> WallExchangeError {
    WallExchangeError::PhaseConduction(ThermalMaterialError::Material(error))
}
