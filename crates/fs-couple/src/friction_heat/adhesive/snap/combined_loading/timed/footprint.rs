//! One owned physical contact footprint on the original caloric timeline.
//!
//! Footprint changes preserve the admitted material/geometry resistance per
//! unit area. They change actual native coverage as well as conductance. A
//! detached owner retains its link identity and refuses any scheduled revival.

use super::super::{JkrCombinedHeat, JkrCombinedHeatError};
use crate::friction_heat::{FrictionHeatRoute, adhesive::snap::JkrSnapHeatError};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{
    CellId, ContactPatchUpdate, DepositionError, DepositionObserver, DepositionSimulation, Event,
};
use std::collections::BTreeMap;

/// Original contact identity, ordered finite receivers and accepted footprint.
/// Only a source-owned timed proposal may change its area or detach it.
#[derive(Clone, Debug)]
pub struct JkrCombinedHeatFootprint {
    link_id: u64,
    cells: [CellId; 2],
    // Every area uses the same reference G/A, avoiding accumulated rescaling
    // roundoff. The original per-side resistances and geometry stay intact.
    reference: Option<FrictionHeatRoute>,
    route: Option<FrictionHeatRoute>,
}

impl JkrCombinedHeatFootprint {
    /// Bind the current attached source to its already admitted physical route.
    /// The exact original receivers and source area must agree; this operation
    /// creates no contact or caloric energy.
    pub fn new(
        history: &JkrCombinedHeat<'_>,
        thermal: &DepositionSimulation<'_>,
        route: FrictionHeatRoute,
    ) -> Result<Self, JkrCombinedHeatError> {
        let result = Self {
            link_id: route.link_id(),
            cells: route.cells(),
            reference: Some(route.clone()),
            route: Some(route),
        };
        result.admit(history, thermal)?;
        Ok(result)
    }

    /// Continue already separated source history while guarding an absent link.
    /// No attached preload, area, conductivity or new bond is manufactured.
    pub fn from_detached(
        history: &JkrCombinedHeat<'_>,
        thermal: &DepositionSimulation<'_>,
        link_id: u64,
    ) -> Result<Self, JkrCombinedHeatError> {
        let result = Self {
            link_id,
            cells: history.receivers.map(|receiver| receiver.cell),
            reference: None,
            route: None,
        };
        result.admit(history, thermal)?;
        Ok(result)
    }

    pub const fn link_id(&self) -> u64 {
        self.link_id
    }

    pub const fn cells(&self) -> [CellId; 2] {
        self.cells
    }

    pub fn area_m2(&self) -> Option<f64> {
        self.route.as_ref().map(FrictionHeatRoute::contact_area_m2)
    }

    pub const fn route(&self) -> Option<&FrictionHeatRoute> {
        self.route.as_ref()
    }

    pub(super) fn admit(
        &self,
        history: &JkrCombinedHeat<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), JkrCombinedHeatError> {
        self.admit_transport(history, thermal)?;
        match (history.state().contact(), self.area_m2()) {
            (Some(point), Some(area))
                if (area - point.contact_area_m2).abs()
                    <= 128.0 * f64::EPSILON * area.max(point.contact_area_m2) =>
            {
                Ok(())
            }
            (None, None) => Ok(()),
            _ => Err(invalid(
                "footprint differs from the accepted combined contact",
            )),
        }
    }

    /// Geometry samples within a ramp retain its incoming accepted mechanical
    /// history until the original sealed endpoint is replayed. Their source
    /// area is supplied by that ramp, while these bindings remain unchanged.
    pub(super) fn admit_transport(
        &self,
        history: &JkrCombinedHeat<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), JkrCombinedHeatError> {
        history.admit_receivers(thermal)?;
        if self.cells != history.receivers.map(|receiver| receiver.cell) {
            return Err(invalid(
                "footprint needs the original ordered caloric receivers",
            ));
        }
        self.admit_topology(thermal)
    }

    fn admit_topology(
        &self,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), JkrCombinedHeatError> {
        if let Some(route) = &self.route {
            route.validate(thermal)?;
        } else if thermal.contact_patch(self.link_id).is_some() {
            return Err(invalid("separated combined contact link was recreated"));
        }
        Ok(())
    }

    /// Admit the complete source radius envelope, including a reversal's peak,
    /// against its original physical faces and currently competing coverage.
    /// One native sweep resolves the free areas of both receivers. No thermal
    /// candidate clone is needed just to reserve this numerical area envelope.
    pub(super) fn reserve_maximum_area(
        &self,
        thermal: &DepositionSimulation<'_>,
        maximum_area_m2: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), JkrCombinedHeatError> {
        checkpoint()?;
        self.admit_topology(thermal)?;
        let Some(current) = &self.route else {
            return if maximum_area_m2 == 0.0 {
                Ok(())
            } else {
                Err(invalid("a separated footprint has zero covered area"))
            };
        };
        self.at_area(maximum_area_m2)?;
        let exposed = thermal.exposed_areas_m2(&self.cells, &mut checkpoint)?;
        for cell in self.cells {
            let available = Dd::from_f64(exposed[&cell]) + Dd::from_f64(current.contact_area_m2());
            let excess = (Dd::from_f64(maximum_area_m2) - available).to_f64();
            let tolerance = 128.0 * f64::EPSILON * maximum_area_m2.max(available.to_f64());
            if !excess.is_finite() || excess > tolerance {
                return Err(invalid(
                    "combined contact path exceeds available surface coverage",
                ));
            }
        }
        checkpoint()?;
        Ok(())
    }

    fn at_area(&self, area_m2: f64) -> Result<FrictionHeatRoute, JkrCombinedHeatError> {
        if !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(invalid("an attached footprint needs positive finite area"));
        }
        let reference = self
            .reference
            .as_ref()
            .ok_or_else(|| invalid("a detached-only owner cannot create a contact"))?;
        for (geometry, distance) in reference.geometry.iter().zip(reference.distances_m) {
            let face = geometry.volume_m3 / (2.0 * distance);
            if !face.is_finite() || area_m2 - face > 64.0 * f64::EPSILON * area_m2.max(face) {
                return Err(invalid(
                    "combined contact exceeds its original prismatic face",
                ));
            }
        }
        let scale = area_m2 / reference.contact_area_m2();
        let conductance = reference.patch.link.conductance_w_k * scale;
        if !scale.is_finite()
            || !conductance.is_finite()
            || (reference.patch.link.conductance_w_k > 0.0 && conductance == 0.0)
        {
            return Err(DepositionError::Arithmetic.into());
        }
        let mut route = reference.clone();
        route.patch.link.area_m2 = area_m2;
        route.patch.covered_area_m2 = [area_m2; 2];
        route.patch.link.conductance_w_k = conductance;
        Ok(route)
    }

    pub(super) fn set_area<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        area_m2: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, Self), JkrCombinedHeatError> {
        checkpoint()?;
        self.admit_topology(thermal)?;
        let current = self
            .route
            .as_ref()
            .ok_or_else(|| invalid("a separated footprint cannot reattach"))?;
        let route = self.at_area(area_m2)?;
        let thermal = if current.patch == route.patch {
            thermal.clone()
        } else {
            thermal.replace_contact_patches(
                &[ContactPatchUpdate {
                    id: self.link_id,
                    expected: Some(current.patch),
                    replacement: Some(route.patch),
                }],
                &mut checkpoint,
            )?
        };
        let mut next = self.clone();
        next.route = Some(route);
        checkpoint()?;
        Ok((thermal, next))
    }

    pub(super) fn detach<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, Self), JkrCombinedHeatError> {
        checkpoint()?;
        self.admit_topology(thermal)?;
        let current = self
            .route
            .as_ref()
            .ok_or_else(|| invalid("the owned contact is already detached"))?;
        let thermal = thermal.replace_contact_patches(
            &[ContactPatchUpdate {
                id: self.link_id,
                expected: Some(current.patch),
                replacement: None,
            }],
            &mut checkpoint,
        )?;
        let mut next = self.clone();
        next.route = None;
        checkpoint()?;
        Ok((thermal, next))
    }

    pub(super) fn guard<'g, 's>(
        &'g self,
        history: &'g JkrCombinedHeat<'s>,
        maximum_area_m2: f64,
        gate: &'g CancelGate,
    ) -> FootprintGuard<'g, 's> {
        FootprintGuard {
            footprint: self,
            history,
            maximum_area_m2,
            gate,
        }
    }
}

#[derive(Clone)]
pub(super) struct FootprintGuard<'g, 's> {
    footprint: &'g JkrCombinedHeatFootprint,
    history: &'g JkrCombinedHeat<'s>,
    maximum_area_m2: f64,
    gate: &'g CancelGate,
}

impl FootprintGuard<'_, '_> {
    fn check(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), JkrCombinedHeatError> {
        let mut check = || {
            if self.gate.is_requested() {
                return Err(DepositionError::Cancelled);
            }
            checkpoint()?;
            if self.gate.is_requested() {
                return Err(DepositionError::Cancelled);
            }
            Ok(())
        };
        check()?;
        self.history.admit_receivers(thermal)?;
        self.footprint
            .reserve_maximum_area(thermal, self.maximum_area_m2, &mut check)?;
        check()?;
        Ok(())
    }
}

impl DepositionObserver for FootprintGuard<'_, '_> {
    type Error = JkrCombinedHeatError;

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.check(thermal, checkpoint)?;
        Ok(BTreeMap::new())
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.check(thermal, &mut || Ok(()))
    }

    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.check(thermal, checkpoint)
    }
}

fn invalid(message: &'static str) -> JkrCombinedHeatError {
    JkrSnapHeatError::Invalid(message).into()
}
