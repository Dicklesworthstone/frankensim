//! Complete material heat transport on the retained three-dimensional fields.
//! The existing interface/layer/radiosity owners evaluate the same endpoint
//! temperatures as conductivity and magnetic storage; no second heat path is
//! stepped and no closed-enclosure transfer becomes external energy supply.

use super::nonlinear::SpatialNonlinearInductionThermalCandidate;
use super::*;
use crate::circuit_enthalpy::{
    contacts::heat::{HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork, closed_balance},
    magnetic::{MagneticThermalError, caloric::MagneticConstitutiveThermalLaw},
};
use fs_em::spatial::{
    nonlinear::SpatialNonlinearMaterial,
    vector_potential::induction::nonlinear::{
        SpatialNonlinearInduction, SpatialNonlinearInductionControls,
    },
};
use fs_exec::Cx;

fn heat_error(error: HeatNetworkError) -> Error {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(MagneticThermalError::Heat(other)),
    }
}

fn close_heat_balance<D>(
    report: &mut SpatialInductionThermalReport<D, HeatNetworkSample>,
    dt: f64,
    tolerance: f64,
) -> Result<(), Error> {
    let sample = &report.coupled.source.network;
    // Both the original native/free-energy account and the physical magnetic
    // internal-energy account must exclude the same internal radiation net.
    closed_balance(
        dt,
        report.coupled.total_energy_residual_j,
        sample,
        tolerance,
    )
    .map_err(heat_error)?;
    closed_balance(dt, report.total_energy_residual_j, sample, tolerance).map_err(heat_error)?;
    let internal_net = finite(dt * sample.radiation_net_power_w())?;
    report.external_heat_j = finite(report.external_heat_j - internal_net)?;
    report.total_energy_residual_j = finite(report.total_energy_residual_j + internal_net)?;
    Ok(())
}

impl PreparedHeatNetwork<'_, '_> {
    /// Continue three-dimensional induction with sourced thermal interfaces,
    /// integrated layer conduction and closed-enclosure radiation.
    ///
    /// Every trial uses this assembly's original material cards and caloric
    /// cells. The existing native A–phi field supplies local conductor-frame
    /// Joule heat and signed Lorentz work; the existing magnetic caloric owner
    /// supplies affine-mu(T) entropy/internal energy and reversible heat. Their
    /// complete field/heat history is published together after both balances
    /// pass. Contact resistance and emissivity are re-evaluated at those same
    /// endpoint temperatures, with their original selected claims retained.
    ///
    /// `regions` and `bindings` retain the original whole-volume, solid-chart
    /// and magnetic-capacity admissions. Geometry, contact areas and enclosure
    /// view factors remain fixed. A thermal contact is not an electrical
    /// contact: graph electrical bindings require their own current/work owner
    /// and refuse here. Thermal boundaries cannot depend on unowned E/J/B/H.
    ///
    /// The top-level external heat and physical residual exclude the closed
    /// enclosure's numerical net; the raw native account remains inspectable.
    /// The existing first-order field/thermal approximation is unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_spatial_induction_candidate<'t, 'm>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialInduction<'m>,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        regions: &[SpatialMagneticRegion],
        volume_tolerance_m3: f64,
        bindings: &[SpatialMagneticThermalLaw],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialInductionControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
    ) -> Result<SpatialInductionThermalCandidate<'t, 'm, HeatNetworkSample>, Error> {
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        check(&cancel)?;
        let bulk = self.thermal_material_network().map_err(heat_error)?;
        self.admit_field_independent_boundaries(cx, &[])
            .map_err(heat_error)?;
        let layout = SpatialMagneticThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            magnetic_regions: regions,
            volume_tolerance_m3,
        };
        let mut next = layout.step_resolved_candidate(
            thermal,
            target_s,
            field,
            endpoint_current_a_m2,
            velocity_m_s,
            bindings,
            controls,
            field_controls,
            assembly,
            initial_energy_increments_j,
            |trial| {
                let coefficients = self.resolve(cx, trial).map_err(heat_error)?;
                if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                    return Err(invalid(
                        "spatial induction cannot discard graph electrical work",
                    ));
                }
                Ok(coefficients.bulk)
            },
            |sample| &sample.material,
            &cancel,
        )?;
        close_heat_balance(
            &mut next.report,
            target_s - thermal.time_s(),
            controls.total_energy_tolerance_j,
        )?;
        check(&cancel)?;
        Ok(next)
    }

    /// Continue the native saturating/fixed-recoil A–phi history with the SAME
    /// sourced interface/layer/radiation and finite-enthalpy assembly.
    ///
    /// Conductivity follows each actual caloric temperature while the complete
    /// original magnetic law remains fixed over its source field range. Only
    /// physical local Joule heat enters the finite stores; backward-Euler
    /// magnetic loss remains a separate numerical export. No restart, fitted
    /// core loss, separate surface-temperature history or heat duplication is
    /// introduced. All source, geometry, receiver and cancellation refusals
    /// retain both accepted owners. The fixed-volume fully-solid limits of
    /// `SpatialMagneticThermalNetwork::step_nonlinear_candidate` still apply.
    #[allow(clippy::too_many_arguments)]
    pub fn step_spatial_nonlinear_induction_candidate<'t, 'm>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialNonlinearInduction<'m>,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        regions: &[SpatialMagneticRegion],
        volume_tolerance_m3: f64,
        controls: CircuitEnthalpySolve,
        field_controls: SpatialNonlinearInductionControls,
        initial_energy_increments_j: &[f64],
    ) -> Result<SpatialNonlinearInductionThermalCandidate<'t, 'm, HeatNetworkSample>, Error> {
        self.step_spatial_nonlinear_induction_with_laws_candidate(
            cx,
            thermal,
            target_s,
            field,
            endpoint_current_a_m2,
            velocity_m_s,
            regions,
            volume_tolerance_m3,
            None,
            controls,
            field_controls,
            initial_energy_increments_j,
        )
    }

    /// Continue the original three-dimensional nonlinear field with sourced
    /// magnetic entropy and the complete contact/layer/radiation heat assembly.
    ///
    /// Linear and source-separable nonlinear thermal laws share the same
    /// material cells, source selections and trial temperatures as conductivity
    /// and boundary heat. The native field owns coefficient work and its
    /// retained A/current history; the common caloric owner supplies reversible
    /// heat and U=psi+T*S. Numerical magnetic loss never enters finite heat.
    /// The geometry, source support and physical scope of
    /// `SpatialMagneticThermalNetwork::step_nonlinear_thermal_candidate` apply.
    #[allow(clippy::too_many_arguments)]
    pub fn step_spatial_nonlinear_induction_thermal_candidate<'t, 'm>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialNonlinearInduction<'m>,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        regions: &[SpatialMagneticRegion],
        volume_tolerance_m3: f64,
        bindings: &[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialNonlinearInductionControls,
        initial_energy_increments_j: &[f64],
    ) -> Result<SpatialNonlinearInductionThermalCandidate<'t, 'm, HeatNetworkSample>, Error> {
        self.step_spatial_nonlinear_induction_with_laws_candidate(
            cx,
            thermal,
            target_s,
            field,
            endpoint_current_a_m2,
            velocity_m_s,
            regions,
            volume_tolerance_m3,
            Some(bindings),
            controls,
            field_controls,
            initial_energy_increments_j,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_spatial_nonlinear_induction_with_laws_candidate<'t, 'm>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialNonlinearInduction<'m>,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        regions: &[SpatialMagneticRegion],
        volume_tolerance_m3: f64,
        bindings: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
        controls: CircuitEnthalpySolve,
        field_controls: SpatialNonlinearInductionControls,
        initial_energy_increments_j: &[f64],
    ) -> Result<SpatialNonlinearInductionThermalCandidate<'t, 'm, HeatNetworkSample>, Error> {
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        check(&cancel)?;
        let bulk = self.thermal_material_network().map_err(heat_error)?;
        let mut magnetic_axes = Vec::new();
        for region in regions {
            check(&cancel)?;
            if let SpatialMagneticRegion::Cell { tetrahedron, cell } = *region {
                if let Some(SpatialNonlinearMaterial::Reversible(law)) =
                    field.problem().materials().get(tetrahedron)
                {
                    magnetic_axes.push((cell, law.source_path().axis.as_str()));
                }
            }
        }
        self.admit_field_independent_boundaries(cx, &magnetic_axes)
            .map_err(heat_error)?;
        let layout = SpatialMagneticThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            magnetic_regions: regions,
            volume_tolerance_m3,
        };
        let mut next = layout.step_constitutive_resolved_candidate(
            thermal,
            target_s,
            field,
            endpoint_current_a_m2,
            velocity_m_s,
            bindings,
            controls,
            field_controls,
            initial_energy_increments_j,
            |trial| {
                let coefficients = self.resolve(cx, trial).map_err(heat_error)?;
                if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                    return Err(invalid(
                        "spatial induction cannot discard graph electrical work",
                    ));
                }
                Ok(coefficients.bulk)
            },
            |sample| &sample.material,
            &cancel,
        )?;
        close_heat_balance(
            &mut next.report,
            target_s - thermal.time_s(),
            controls.total_energy_tolerance_j,
        )?;
        check(&cancel)?;
        Ok(next)
    }
}
