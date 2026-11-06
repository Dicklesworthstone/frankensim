//! Time-integrated contact coverage on the original finite thermal cells.
//!
//! A native mechanical owner supplies the actual mean and endpoint real areas.
//! During backward Euler, temperatures and conductivities are endpoint values;
//! replacing coverage by its time mean integrates that supplied area path in
//! both contact conduction and ambient exposure. Cell geometry, caloric laws,
//! material identity, mass and source histories remain with their original owners.
//! This does not describe moving face centres, separation or new contact discovery.

use crate::friction_heat::{FrictionHeatRoute, finite, invalid};
use fs_thermal::enthalpy::{ContactPatchUpdate, DepositionError, DepositionSimulation};
use std::borrow::Cow;

#[derive(Clone, Copy)]
pub(super) struct ContactAreaPath<'a> {
    pub route: &'a FrictionHeatRoute,
    pub mean_area_m2: f64,
    pub endpoint_area_m2: f64,
}

/// One batched geometry stage, followed by the same owner's caloric solve.
pub(super) struct ContactAreaEvolution<'s, 'm> {
    thermal: Cow<'s, DepositionSimulation<'m>>,
    routes: Vec<FrictionHeatRoute>,
    endpoint_routes: Vec<FrictionHeatRoute>,
}

impl<'s, 'm> ContactAreaEvolution<'s, 'm> {
    pub fn new(
        initial: &'s DepositionSimulation<'m>,
        paths: &[ContactAreaPath<'_>],
        maximum_interfaces: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if paths.is_empty() || paths.len() > maximum_interfaces {
            return Err(invalid("bounded nonempty evolving contact areas required"));
        }
        let mut ordered: Vec<_> = paths.iter().collect();
        ordered.sort_by_key(|path| path.route.link_id());
        let mut routes = Vec::with_capacity(paths.len());
        let mut endpoint_routes = Vec::with_capacity(paths.len());
        let mut updates = Vec::with_capacity(paths.len());
        for (i, path) in ordered.iter().enumerate() {
            checkpoint()?;
            if i > 0 && ordered[i - 1].route.link_id() == path.route.link_id() {
                return Err(invalid(
                    "one evolving area owner per original contact required",
                ));
            }
            path.route.validate(initial)?;
            let incoming = path.route.contact_area_m2();
            let high = incoming.max(path.endpoint_area_m2);
            let low = incoming.min(path.endpoint_area_m2);
            let allowance = 128.0 * f64::EPSILON * high;
            if !path.mean_area_m2.is_finite()
                || !path.endpoint_area_m2.is_finite()
                || path.mean_area_m2 <= 0.0
                || path.endpoint_area_m2 <= 0.0
                || path.mean_area_m2 < low - allowance
                || path.mean_area_m2 > high + allowance
            {
                return Err(invalid(
                    "positive mean area within the original path extrema required",
                ));
            }
            let mean = with_area(path.route, path.mean_area_m2)?;
            let endpoint = with_area(path.route, path.endpoint_area_m2)?;
            updates.push(ContactPatchUpdate {
                id: path.route.link_id(),
                expected: Some(path.route.patch),
                replacement: Some(mean.patch),
            });
            routes.push(mean);
            endpoint_routes.push(endpoint);
        }
        let thermal = if updates.iter().all(|u| u.expected == u.replacement) {
            Cow::Borrowed(initial)
        } else {
            // One complete batch retains shared-cell coverage and the native
            // surface-capacity admission; no order-dependent intermediate patch.
            Cow::Owned(initial.replace_contact_patches(&updates, &mut checkpoint)?)
        };
        checkpoint()?;
        Ok(Self {
            thermal,
            routes,
            endpoint_routes,
        })
    }

    pub fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }

    /// Canonical link order, independently of mechanical construction order.
    pub fn routes(&self) -> &[FrictionHeatRoute] {
        &self.routes
    }

    /// Install endpoint coverage only after the original caloric solve accepts.
    /// Native pending-event, coverage and cancellation refusals preserve every
    /// caller-owned state. Endpoint routes retain the original reference
    /// conductance per unit area; each next heat trial still resolves actual k.
    pub fn finish(
        self,
        accepted: DepositionSimulation<'m>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, Vec<FrictionHeatRoute>), DepositionError> {
        checkpoint()?;
        if accepted.time_s() <= self.thermal.time_s() {
            return Err(invalid(
                "evolving coverage requires the accepted later thermal endpoint",
            ));
        }
        let mut updates = Vec::with_capacity(self.routes.len());
        for (mean, endpoint) in self.routes.iter().zip(&self.endpoint_routes) {
            checkpoint()?;
            mean.validate(&accepted)?;
            updates.push(ContactPatchUpdate {
                id: mean.link_id(),
                expected: Some(mean.patch),
                replacement: Some(endpoint.patch),
            });
        }
        let accepted = if updates.iter().all(|u| u.expected == u.replacement) {
            accepted
        } else {
            accepted.replace_contact_patches(&updates, &mut checkpoint)?
        };
        checkpoint()?;
        Ok((accepted, self.endpoint_routes))
    }
}

fn with_area(route: &FrictionHeatRoute, area: f64) -> Result<FrictionHeatRoute, DepositionError> {
    for (geometry, distance) in route.geometry.iter().zip(route.distances_m) {
        let face = finite(geometry.volume_m3 / (2.0 * distance))?;
        if area - face > 64.0 * f64::EPSILON * area.max(face) {
            return Err(invalid(
                "evolving real contact exceeds its original prismatic face",
            ));
        }
    }
    let scale = finite(area / route.contact_area_m2())?;
    let conductance = finite(route.patch.link.conductance_w_k * scale)?;
    if route.patch.link.conductance_w_k > 0.0 && conductance == 0.0 {
        return Err(DepositionError::Arithmetic);
    }
    let mut result = route.clone();
    result.patch.link.area_m2 = area;
    result.patch.covered_area_m2 = [area; 2];
    result.patch.link.conductance_w_k = conductance;
    Ok(result)
}
