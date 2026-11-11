//! Current transport measures on the original finite-cell heat owner.
//!
//! Geometry is supplied by an independently solved mechanical configuration.
//! The existing series-resistance law still describes resolved material paths.
//! Planar distances are physical normal lengths. An explicitly radial material
//! owner may lower its retained physical coordinates to A_face*integral(ds/A),
//! an exact resistance-equivalent distance that is NEVER a cell volume/length
//! measure. These metrics do not turn a two-point graph into a general Fourier
//! operator on nonorthogonal distorted tetrahedra. The geometric owner must
//! justify its reduction, current footprints and exposed-area coverage.
//! The reference cell geometry, inventory and material source never change.

use super::*;

/// Both CURRENT material resistance coordinates and their common flux area.
/// Planar paths use physical normal half-distances. An explicit radial owner
/// uses A_contact*integral(ds/A), retaining its physical coordinates separately.
/// Endpoint order is the original link's ordered `(a, b)` pair.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurrentConductiveLink {
    pub area_m2: f64,
    pub half_distances_m: [f64; 2],
}

/// One CURRENT boundary footprint and its material resistance coordinate.
/// An explicit radial owner supplies A_face*integral(ds/A), not radial thickness.
/// Surface resistance retains its original areal units and material meaning.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurrentConductiveBoundary {
    pub area_m2: f64,
    pub half_distance_m: f64,
}

/// Complete current metrics for precisely the original observer's channels.
///
/// The ordered maps preserve physical channel identity. Empty maps are valid
/// only when that original channel is absent. Positive measures keep the same
/// topology; opening contact, vanishing layers and newly exposed surfaces need
/// an explicit topology/coverage transition, not a zero or clipped metric.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CurrentConductionGeometry {
    pub links: BTreeMap<u64, CurrentConductiveLink>,
    pub reservoirs: BTreeMap<CellId, CurrentConductiveBoundary>,
    pub ambient_faces: BTreeMap<CellId, CurrentConductiveBoundary>,
}

fn admit_measure(value: f64) -> Result<(), DepositionError> {
    if positive(value) {
        Ok(())
    } else {
        Err(DepositionError::InvalidInput(
            "current conduction areas and normal distances must be finite and positive",
        ))
    }
}

impl CurrentConductionGeometry {
    fn admit(
        &self,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        for face in self.links.values() {
            checkpoint()?;
            admit_measure(face.area_m2)?;
            for distance in face.half_distances_m {
                admit_measure(distance)?;
            }
        }
        for face in self.reservoirs.values().chain(self.ambient_faces.values()) {
            checkpoint()?;
            admit_measure(face.area_m2)?;
            admit_measure(face.half_distance_m)?;
        }
        checkpoint()
    }

    /// Largest relative difference of ANY supplied area or resistance distance.
    /// Opposing geometric changes cannot cancel. The denominator is the larger
    /// positive endpoint measure, so this comparison has no dimensional floor.
    /// This bounds a heat-metric iterate change, not physical-coordinate or
    /// geometric discretization error; the geometric owner retains those.
    pub fn relative_difference(
        &self,
        other: &Self,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, DepositionError> {
        checkpoint()?;
        if self.links.keys().ne(other.links.keys())
            || self.reservoirs.keys().ne(other.reservoirs.keys())
            || self.ambient_faces.keys().ne(other.ambient_faces.keys())
        {
            return Err(DepositionError::InvalidInput(
                "current conduction geometry changed physical channel identity",
            ));
        }
        self.admit(&mut checkpoint)?;
        other.admit(&mut checkpoint)?;
        let relative = |a: f64, b: f64| (a - b).abs() / a.max(b);
        let mut worst = 0.0_f64;
        for (a, b) in self.links.values().zip(other.links.values()) {
            checkpoint()?;
            worst = worst.max(relative(a.area_m2, b.area_m2));
            for (a, b) in a.half_distances_m.into_iter().zip(b.half_distances_m) {
                worst = worst.max(relative(a, b));
            }
        }
        for (a, b) in self
            .reservoirs
            .values()
            .chain(self.ambient_faces.values())
            .zip(
                other
                    .reservoirs
                    .values()
                    .chain(other.ambient_faces.values()),
            )
        {
            checkpoint()?;
            worst = worst.max(relative(a.area_m2, b.area_m2));
            worst = worst.max(relative(a.half_distance_m, b.half_distance_m));
        }
        checkpoint()?;
        Ok(worst)
    }
}

impl<'a> ConductionObserver<'a> {
    /// Prepare complete current metrics without modifying the accepted owner.
    ///
    /// The original observer still owns every conductivity and resistance law.
    /// Its original cell, mass, caloric material, link, bath and reference-area
    /// identities are rechecked. Only areas and transport resistance distances used
    /// by that SAME series law change. Source conductivity refreshes at every
    /// native heat substep; current geometry is held for this proposed interval.
    ///
    /// This first-order endpoint-geometry composition needs independent physical
    /// interval and thermal-substep refinement. Explicit timeline fluxes and CFL,
    /// and implicit chain/graph fluxes and Jacobians, consume the same metrics.
    /// Implicit Newton trials and SDIRK stages hold this configuration fixed;
    /// removing the thermal CFL does not resolve mechanical geometry evolution.
    pub fn with_current_geometry(
        &self,
        thermal: &DepositionSimulation<'_>,
        geometry: CurrentConductionGeometry,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if !self.outside.is_empty() {
            return Err(DepositionError::InvalidInput(
                "current conduction metrics are not admitted beside links carried outside",
            ));
        }
        if self.has_phase_layers
            || geometry.links.keys().ne(self.links.keys())
            || geometry.reservoirs.keys().ne(self.reservoirs.keys())
            || geometry.ambient_faces.keys().ne(self.ambient_faces.keys())
        {
            return Err(DepositionError::InvalidInput(
                "complete original conduction channels and a single geometry owner required",
            ));
        }
        geometry.admit(&mut checkpoint)?;
        let mut next = self.clone();
        next.current_geometry = Some(geometry);
        next.prepare(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(next)
    }

    /// Actual metrics used by an explicitly selected deformed interval.
    /// `None` retains the original reference/planar-layer behavior.
    pub fn current_geometry(&self) -> Option<&CurrentConductionGeometry> {
        self.current_geometry.as_ref()
    }

    // Transport resistance measures before any separately owned planar-density
    // response. Current mechanical metrics and planar layers are exclusive.
    pub(super) fn link_metrics(&self, binding: &BoundLink) -> CurrentConductiveLink {
        self.current_geometry.as_ref().map_or(
            CurrentConductiveLink {
                area_m2: binding.area_m2,
                half_distances_m: binding.geometry.half_distances_m,
            },
            |current| current.links[&binding.geometry.link_id],
        )
    }

    pub(super) fn reservoir_metrics(&self, binding: &BoundReservoir) -> CurrentConductiveBoundary {
        self.current_geometry.as_ref().map_or(
            CurrentConductiveBoundary {
                area_m2: binding.area_m2,
                half_distance_m: binding.geometry.half_distance_m,
            },
            |current| current.reservoirs[&binding.geometry.cell_id],
        )
    }

    pub(super) fn ambient_metrics(&self, binding: &BoundAmbient) -> CurrentConductiveBoundary {
        self.current_geometry.as_ref().map_or(
            CurrentConductiveBoundary {
                area_m2: binding.area_m2,
                half_distance_m: binding.geometry.half_distance_m,
            },
            |current| current.ambient_faces[&binding.geometry.cell_id],
        )
    }
}
