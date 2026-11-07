//! Export native Archard removal from the same finite thermal receivers.
//!
//! Heat/sliding is advanced first. This optional endpoint operation then removes
//! each body's newly worn volume, its original sourced density times that volume,
//! and the material's current specific enthalpy through the existing thermal
//! inventory owner. It changes neither time, temperature, nor mechanical history.
//! The next heat solve uses the actual remaining masses and shortened geometry.
//!
//! The caller supplies the remaining prismatic geometry and every affected
//! contact patch. Cross-sectional area and material/grid registration are held:
//! only shortening and normal translation are admitted. Local recession is not
//! silently converted into a spatial thermal mesh. The existing cumulative wear
//! and half-space bounds remain; this heat-then-export split requires refinement.
//! Debris motion, subsequent debris heating, surface creation energy, and
//! mechanical work needed to move the thermal control volumes are not modeled.

use fs_exec::Cx;
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_thermal::enthalpy::removal::{
    MaterialRemoval, MaterialRemovalAudit, MaterialRemovalTolerances,
};
use fs_thermal::enthalpy::{CellGeometry, ContactPatchUpdate, DepositionSimulation};

use super::{
    FrictionHeatRoute, PlasticSlidingHeatContact, PlasticSlidingHeatError, checkpoint, finite,
};

/// Original contact history and remaining thermal inventory for joint adoption.
pub struct PlasticSlidingWearExportCandidate<'a, 't> {
    contact: PlasticSlidingHeatContact<'a>,
    thermal: DepositionSimulation<'t>,
    report: MaterialRemovalAudit,
}

impl<'a, 't> PlasticSlidingWearExportCandidate<'a, 't> {
    pub const fn contact(&self) -> &PlasticSlidingHeatContact<'a> {
        &self.contact
    }

    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }

    /// Native requested/represented volumes, masses, and signed enthalpy outflow.
    pub const fn report(&self) -> &MaterialRemovalAudit {
        &self.report
    }

    pub fn into_parts(
        self,
    ) -> (
        PlasticSlidingHeatContact<'a>,
        DepositionSimulation<'t>,
        MaterialRemovalAudit,
    ) {
        (self.contact, self.thermal, self.report)
    }
}

impl<'a> PlasticSlidingHeatContact<'a> {
    /// Original A/B wear accumulated since binding or the last accepted export.
    ///
    /// Wear already present when the contact is bound is the inventory baseline;
    /// it is never removed retroactively. Volumes are requests in the original
    /// Archard convention; the native export report retains arithmetic roundoff
    /// between each request and the actually represented inventory decrement.
    pub fn unexported_wear_volume_m3(&self) -> [f64; 2] {
        let cumulative = self.state.wear_state().volumes_m3();
        std::array::from_fn(|side| cumulative[side] - self.thermal_wear_checkpoint_m3[side])
    }

    /// Export all newly worn material at this accepted heat/contact endpoint.
    ///
    /// Supply complete remaining geometry and explicit updates for every old
    /// thermal link incident on a changed receiver. The original contact link,
    /// ordered cells, loaded area, prismatic face areas and tangent registration
    /// must survive. Conductivities belong to the original materials at their
    /// unchanged current temperatures; their source evidence remains caller-owned.
    ///
    /// Density comes from each original sourced caloric chart at current h. The
    /// native owner refuses if it disagrees with stored mass/volume: no expansion
    /// model or inferred m/V density is substituted. Remaining h and temperature
    /// are unchanged, including for a negative enthalpy datum. Outgoing enthalpy
    /// is material transport, never an extra heat source or cooling term.
    ///
    /// Pending events exactly at this endpoint remain pending. A source, geometry,
    /// inventory, patch, or cancellation refusal returns no partial update.
    /// Adopt both candidates together; later heat steps retain the export
    /// checkpoint so the same worn volume cannot be removed twice.
    #[allow(clippy::too_many_arguments)]
    pub fn export_wear_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        replacement_geometry: [CellGeometry; 2],
        patches: &[ContactPatchUpdate],
        conductivity_w_m_k: [f64; 2],
        tolerances: MaterialRemovalTolerances,
        cx: &Cx<'_>,
    ) -> Result<PlasticSlidingWearExportCandidate<'a, 't>, PlasticSlidingHeatError> {
        checkpoint(cx)?;
        self.admit_thermal(thermal)?;
        let area = self.admit_footprint(&self.state, cx)?;
        let volumes = self.unexported_wear_volume_m3();
        let (thermal, route, report) = export_original_receivers(
            thermal,
            &self.route,
            self.caloric,
            volumes,
            area,
            replacement_geometry,
            patches,
            conductivity_w_m_k,
            tolerances,
            cx,
        )?;
        let mut contact = self.clone();
        contact.route = route;
        contact.thermal_wear_checkpoint_m3 = self.state.wear_state().volumes_m3();
        contact.admit_thermal(&thermal)?;
        checkpoint(cx)?;
        Ok(PlasticSlidingWearExportCandidate {
            contact,
            thermal,
            report,
        })
    }
}

/// Both concrete contact owners admit their original thermal checkpoint first.
/// Inventory removal, caloric transport and prismatic registration then have one
/// native path, independent of how the contact obtained its hardness field.
#[allow(clippy::too_many_arguments)]
pub(super) fn export_original_receivers<'t>(
    thermal: &DepositionSimulation<'t>,
    original_route: &FrictionHeatRoute,
    caloric: [&ResolvedHeatCapacityCurve; 2],
    volumes: [f64; 2],
    area: f64,
    replacement_geometry: [CellGeometry; 2],
    patches: &[ContactPatchUpdate],
    conductivity_w_m_k: [f64; 2],
    tolerances: MaterialRemovalTolerances,
    cx: &Cx<'_>,
) -> Result<
    (
        DepositionSimulation<'t>,
        FrictionHeatRoute,
        MaterialRemovalAudit,
    ),
    PlasticSlidingHeatError,
> {
    if volumes
        .iter()
        .any(|volume| !volume.is_finite() || *volume < 0.0)
    {
        return Err(PlasticSlidingHeatError::Invalid(
            "thermal wear checkpoint exceeds the accepted original wear history",
        ));
    }
    let cells = original_route.cells();
    let mut removals = Vec::with_capacity(2);
    for side in 0..2 {
        let expected = thermal.cells()[&cells[side]].clone();
        let curve = caloric[side].curve();
        let values = curve.values_at_specific_enthalpy(expected.specific_h_j_kg)?;
        let caloric_identity = thermal
            .material_caloric_identity(expected.material_id)?
            .ok_or(PlasticSlidingHeatError::Invalid(
                "wear export requires the original caloric identity",
            ))?;
        removals.push(MaterialRemoval {
            expected,
            removed_volume_m3: volumes[side],
            density_kg_m3: values.bulk_density_kg_m3,
            replacement_geometry: replacement_geometry[side],
            caloric_identity,
        });
    }
    let (thermal, report) =
        thermal.remove_material_before_pending_events(&removals, patches, tolerances, || {
            checkpoint(cx)
        })?;
    let route = FrictionHeatRoute::from_link_with_resistance(
        &thermal,
        original_route.link_id(),
        conductivity_w_m_k,
        original_route.area_resistance_m2_k_w(),
    )?;
    if route.cells() != cells
        || (route.contact_area_m2() - area).abs()
            > 64.0 * f64::EPSILON * route.contact_area_m2().max(area)
    {
        return Err(PlasticSlidingHeatError::Invalid(
            "wear export must retain ordered original cells and actual loaded area",
        ));
    }
    admit_registration(original_route, &route)?;
    Ok((thermal, route, report))
}

// Both routes have already passed the native opposing-prism admission.
fn shared_face(route: &FrictionHeatRoute) -> Result<[f64; 3], PlasticSlidingHeatError> {
    let a = route.geometry[0];
    let b = route.geometry[1];
    for face in [a.start_m, a.end_m] {
        if face == b.start_m || face == b.end_m {
            return Ok(face);
        }
    }
    Err(PlasticSlidingHeatError::Invalid(
        "wear export requires the original shared thermal face",
    ))
}

fn inward_direction(
    route: &FrictionHeatRoute,
    side: usize,
    face: [f64; 3],
) -> Result<[f64; 3], PlasticSlidingHeatError> {
    let geometry = route.geometry[side];
    let other = if geometry.start_m == face {
        geometry.end_m
    } else {
        geometry.start_m
    };
    let length = finite(2.0 * route.distances_m()[side])?;
    let mut direction = [0.0; 3];
    for axis in 0..3 {
        direction[axis] = finite((other[axis] - face[axis]) / length)?;
    }
    Ok(direction)
}

fn admit_registration(
    before: &FrictionHeatRoute,
    after: &FrictionHeatRoute,
) -> Result<(), PlasticSlidingHeatError> {
    let old_face = shared_face(before)?;
    let new_face = shared_face(after)?;
    let normal = inward_direction(before, 0, old_face)?;
    let mut shift = [0.0; 3];
    let mut projection = 0.0;
    let scale = before.distances_m()[0]
        .max(before.distances_m()[1])
        .max(after.distances_m()[0])
        .max(after.distances_m()[1]);
    for axis in 0..3 {
        shift[axis] = finite(new_face[axis] - old_face[axis])?;
        projection = finite(projection + shift[axis] * normal[axis])?;
    }
    for axis in 0..3 {
        if finite(shift[axis] - projection * normal[axis])?.abs() > 256.0 * f64::EPSILON * scale {
            return Err(PlasticSlidingHeatError::Invalid(
                "wear export cannot translate the original grid tangentially",
            ));
        }
    }
    for side in 0..2 {
        let old_direction = inward_direction(before, side, old_face)?;
        let new_direction = inward_direction(after, side, new_face)?;
        if (0..3)
            .any(|axis| (old_direction[axis] - new_direction[axis]).abs() > 256.0 * f64::EPSILON)
        {
            return Err(PlasticSlidingHeatError::Invalid(
                "wear export cannot rotate the original material/grid registration",
            ));
        }
        let old_area =
            finite(before.geometry[side].volume_m3 / (2.0 * before.distances_m()[side]))?;
        let new_area = finite(after.geometry[side].volume_m3 / (2.0 * after.distances_m()[side]))?;
        if (old_area - new_area).abs() > 128.0 * f64::EPSILON * old_area.max(new_area) {
            return Err(PlasticSlidingHeatError::Invalid(
                "wear export must retain the original prismatic face areas",
            ));
        }
    }
    Ok(())
}
