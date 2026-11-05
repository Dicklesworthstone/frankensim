//! Original bulk conductivity on the same kinetic states as interface exchange.
//!
//! Each scalar conductivity is an admitted isotropic, static material law at
//! actual T and explicitly owned crystal coordinates. Fixed positive distances
//! describe two homogeneous layers; no geometric motion, homogenization rule,
//! bulk temperature profile or conductivity calibration is inferred here.

use super::{CaloricBinding, CrystallizationInterfaceError, invalid};
use crate::phase::crystallization::CrystallizationCaloricPair;
use crate::phase::crystallization::properties::{
    CrystallizationMaterialPoint, CrystallizationPropertySource, CrystallizationPropertyState,
};
use crate::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{Pressure, QuantitySpec};
use fs_thermal::crystallization::CrystallizationSolve;
use fs_thermal::crystallization::pressure::parcel::IsobaricCrystallizationParcel;

/// Both original bulk source answers at the ACTUAL accepted endpoint histories.
#[derive(Clone, Debug, PartialEq)]
pub struct CrystallizationBulkConductionResponse {
    pub sources: [CrystallizationPropertyState; 2],
    pub conductivity_w_m_k: [f64; 2],
}

/// Immutable two-layer preparation; original claims are selected once.
#[derive(Clone)]
pub struct CrystallizationBulkConduction<'a> {
    sources: [CrystallizationPropertySource<'a>; 2],
    pub(super) bindings: [CaloricBinding; 2],
    distances_m: [f64; 2],
}
impl core::fmt::Debug for CrystallizationBulkConduction<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("CrystallizationBulkConduction")
            .field("sources", &self.original_sources())
            .field("bindings", &self.bindings)
            .field("distances_m", &self.distances_m)
            .finish()
    }
}

impl<'a> CrystallizationBulkConduction<'a> {
    /// Prepare the SAME original material cards, caloric pairs and kinetic laws.
    ///
    /// Conductivity uses its canonical quantity and nonnegative admissibility:
    /// an explicitly supported zero conductivity insulates. A source may use
    /// absolute crystalline mass fraction, original-density volume fraction,
    /// both, or neither; omitting a coordinate does not establish independence.
    /// Any supplied bulk pressure coordinate must equal the parcel's pressure.
    /// Other admitted fixed external conditions retain their original values.
    /// Later domain failures cannot select a different conductivity claim.
    pub fn new(
        points: [CrystallizationMaterialPoint<'a>; 2],
        pairs: [CrystallizationCaloricPair<'a>; 2],
        parcels: [&IsobaricCrystallizationParcel<'_, 'a>; 2],
        distances_m: [f64; 2],
        selections: [MaterialPropertySelection; 2],
        solves: [CrystallizationSolve; 2],
    ) -> Result<Self, CrystallizationInterfaceError> {
        if distances_m.iter().any(|d| !d.is_finite() || *d <= 0.0) {
            return Err(invalid(
                "bulk contact needs two explicit positive finite distances",
            ));
        }
        for side in 0..2 {
            admit_pressure(points[side], parcels[side].pressure_pa())?;
        }
        let bindings = [
            CaloricBinding::new(&pairs[0], parcels[0])?,
            CaloricBinding::new(&pairs[1], parcels[1])?,
        ];
        let requirement = ScalarPropertyRequirement::try_new(
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            ScalarAdmissibility::NonNegative,
        )?;
        let [left_selection, right_selection] = selections;
        let sources = [
            CrystallizationPropertySource::new(
                pairs[0],
                *parcels[0].law(),
                parcels[0].state(),
                points[0],
                core::slice::from_ref(&requirement),
                left_selection,
                solves[0],
            )?,
            CrystallizationPropertySource::new(
                pairs[1],
                *parcels[1].law(),
                parcels[1].state(),
                points[1],
                core::slice::from_ref(&requirement),
                right_selection,
                solves[1],
            )?,
        ];
        Ok(Self {
            sources,
            bindings,
            distances_m,
        })
    }

    /// Fixed layer distances in the same A/B order as the original interface.
    pub const fn distances_m(&self) -> [f64; 2] {
        self.distances_m
    }

    pub fn original_sources(&self) -> [&ResolvedMaterialStatePoint; 2] {
        self.sources
            .each_ref()
            .map(|source| source.original_source())
    }

    /// Refresh both pinned bulk claims without advancing physical history.
    /// The caller can use this same answer in the original network conductance
    /// owner. Either both complete source snapshots are returned or neither is.
    pub fn evaluate(
        &self,
        parcels: [&IsobaricCrystallizationParcel<'_, '_>; 2],
    ) -> Result<CrystallizationBulkConductionResponse, CrystallizationInterfaceError> {
        for (binding, parcel) in self.bindings.iter().zip(parcels) {
            binding.admit(parcel)?;
        }
        let sources = [
            self.sources[0].evaluate(parcels[0].state())?,
            self.sources[1].evaluate(parcels[1].state())?,
        ];
        let conductivity_w_m_k = sources
            .each_ref()
            .map(|source| source.material().properties()[0].value_si());
        Ok(CrystallizationBulkConductionResponse {
            sources,
            conductivity_w_m_k,
        })
    }
}

fn admit_pressure(
    point: CrystallizationMaterialPoint<'_>,
    pressure_pa: f64,
) -> Result<(), CrystallizationInterfaceError> {
    for (axis, quantity) in point.context.axis_quantities() {
        if quantity.dims() != Pressure::DIMS {
            continue;
        }
        let is_pressure = *quantity == QuantitySpec::dimensional(Pressure::DIMS)
            || quantity.semantic_type().is_some_and(|semantic| {
                semantic.kind() == QuantityKind::Pressure && semantic.form() == ValueForm::Static
            });
        if !is_pressure || point.context.axes().get(axis).copied() != Some(pressure_pa) {
            return Err(invalid(
                "bulk conductivity pressure must be the actual fixed parcel pressure",
            ));
        }
    }
    Ok(())
}
