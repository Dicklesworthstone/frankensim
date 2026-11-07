//! Actual worn mass and enthalpy exported from the original depth-hardness bodies.
//!
//! The shared native receiver export owns material inventory, current carried
//! enthalpy and prismatic registration. Depth hardness, normal indentation and
//! genuine shear history remain with the accepted contact owner.

use fs_exec::Cx;
use fs_thermal::enthalpy::removal::{MaterialRemovalAudit, MaterialRemovalTolerances};
use fs_thermal::enthalpy::{CellGeometry, ContactPatchUpdate, DepositionSimulation};

use super::super::{PlasticSlidingHeatError, checkpoint, removal::export_original_receivers};
use super::DepthSlidingHeatContact;

/// Complete depth-hardness contact and remaining thermal inventory for adoption.
pub struct DepthSlidingWearExportCandidate<'a, 't> {
    contact: DepthSlidingHeatContact<'a>,
    thermal: DepositionSimulation<'t>,
    report: MaterialRemovalAudit,
}

impl<'a, 't> DepthSlidingWearExportCandidate<'a, 't> {
    pub const fn contact(&self) -> &DepthSlidingHeatContact<'a> {
        &self.contact
    }

    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }

    /// Original requested/represented removal and signed enthalpy outflow.
    pub const fn report(&self) -> &MaterialRemovalAudit {
        &self.report
    }

    pub fn into_parts(
        self,
    ) -> (
        DepthSlidingHeatContact<'a>,
        DepositionSimulation<'t>,
        MaterialRemovalAudit,
    ) {
        (self.contact, self.thermal, self.report)
    }
}

impl<'a> DepthSlidingHeatContact<'a> {
    /// Actual original A/B wear since binding or the last accepted export.
    /// Existing wear at binding is the thermal inventory baseline and is never
    /// removed retroactively. Export does not reset cumulative wear limits.
    pub fn unexported_wear_volume_m3(&self) -> [f64; 2] {
        let cumulative = self.source.wear_state().volumes_m3();
        std::array::from_fn(|side| cumulative[side] - self.thermal_wear_checkpoint_m3[side])
    }

    /// Export newly worn material after accepting its heat/contact endpoint.
    ///
    /// The original sourced density supplies mass; the native thermal owner
    /// checks mass/volume consistency and exports the current signed specific
    /// enthalpy. Remaining h, temperature, time, depth profiles, indentation and
    /// actual shear histories stay unchanged. Outgoing material enthalpy is not
    /// a second heat sink, even when its chosen datum is negative.
    ///
    /// Supply shortened prismatic geometry and all affected contact patches.
    /// Original cross sections, ordered cells and tangent registration survive;
    /// no local depth field is silently remeshed into thermal cells. The next
    /// heat solve uses actual remaining masses and conduction distances. Pending
    /// endpoint events remain pending, and any refusal leaves both owners intact.
    #[allow(clippy::too_many_arguments)]
    pub fn export_wear_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        replacement_geometry: [CellGeometry; 2],
        patches: &[ContactPatchUpdate],
        conductivity_w_m_k: [f64; 2],
        tolerances: MaterialRemovalTolerances,
        cx: &Cx<'_>,
    ) -> Result<DepthSlidingWearExportCandidate<'a, 't>, PlasticSlidingHeatError> {
        checkpoint(cx)?;
        self.admit_thermal(thermal)?;
        let area = self.admit_footprint(&self.source, cx)?;
        let (thermal, route, report) = export_original_receivers(
            thermal,
            &self.route,
            self.caloric,
            self.unexported_wear_volume_m3(),
            area,
            replacement_geometry,
            patches,
            conductivity_w_m_k,
            tolerances,
            cx,
        )?;
        let mut contact = self.clone();
        contact.route = route;
        contact.thermal_wear_checkpoint_m3 = self.source.wear_state().volumes_m3();
        contact.admit_thermal(&thermal)?;
        checkpoint(cx)?;
        Ok(DepthSlidingWearExportCandidate {
            contact,
            thermal,
            report,
        })
    }
}
