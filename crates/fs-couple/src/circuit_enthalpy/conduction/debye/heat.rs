//! The original complete heat assembly receives the original dielectric losses.

use super::*;
use crate::circuit_enthalpy::contacts::heat::{
    HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork,
};
use fs_exec::Cx;

fn heat(error: HeatNetworkError) -> Error {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ConductionThermalError::Heat(other)),
    }
}

impl PreparedHeatNetwork<'_, '_> {
    /// Planar convenience entry point: solve the field once, then delegate to
    /// [`Self::step_dielectric_field_candidate`] for the complete original heat
    /// assembly and joint publication of both physical histories.
    #[allow(clippy::too_many_arguments)]
    pub fn step_debye_candidate<'t>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        dielectric: &PlanarDebyeDielectric,
        prepared: &PreparedDebyeStep,
        endpoint_voltage_v: &[(usize, f64)],
        regions: &[ConductionRegion],
        volume_tolerance_m3: f64,
        field_controls: ElectrostaticControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<DebyeEnthalpyCandidate<'t, HeatNetworkSample>, Error> {
        if dielectric.time_s().to_bits() != thermal.time_s().to_bits()
            || (dielectric.time_s() + prepared.time_step_s()).to_bits() != target_s.to_bits()
        {
            return Err(invalid(
                "matching Debye/thermal clocks and prepared interval required",
            ));
        }
        let proposal = dielectric
            .step_candidate(prepared, endpoint_voltage_v, field_controls, &mut || {
                cx.checkpoint().is_err()
            })
            .map_err(field)?;
        self.step_dielectric_field_candidate(
            cx,
            thermal,
            &proposal,
            regions,
            volume_tolerance_m3,
            controls,
            initial_energy_increments_j,
        )
    }

    /// Advance the original finite enthalpy timeline with a native planar,
    /// axisymmetric or tetrahedral dielectric candidate. Sourced Fourier links,
    /// thermal interfaces, integrated layers and closed gray-diffuse radiation
    /// all resolve at the same actual thermal endpoint.
    ///
    /// Total physical field loss is transferred once per original electrical
    /// cell; its report keeps polarization and Ohmic loss separate. Stored
    /// electric energy is never deposited. Enclosure radiation is internal
    /// transfer: its computed net is excluded from external supply in the joint
    /// balance. The source, fixed-volume and fully solid restrictions of
    /// [`ConductionThermalNetwork::step_dielectric_field_candidate`] apply.
    /// Graph electrical branches and unowned field-dependent thermal boundary
    /// laws refuse. No field or thermal history publishes on refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn step_dielectric_field_candidate<'t, F: DielectricFieldState, const D: usize>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        proposal: &DielectricFieldCandidate<F, D>,
        regions: &[ConductionRegion],
        volume_tolerance_m3: f64,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<DebyeEnthalpyCandidate<'t, HeatNetworkSample, F, D>, Error> {
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        check(&cancel)?;
        let bulk = self.thermal_material_network().map_err(heat)?;
        self.admit_field_independent_boundaries(cx, &[])
            .map_err(heat)?;
        let layout = ConductionThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            regions,
            volume_tolerance_m3,
        };
        layout.step_dielectric_field_resolved_candidate(
            thermal,
            proposal,
            controls,
            initial_energy_increments_j,
            |trial| {
                let resolved = self.resolve(cx, trial).map_err(heat)?;
                if !resolved.contacts.is_empty() || !resolved.bulk.resistors.is_empty() {
                    return Err(invalid(
                        "Debye field cannot silently omit graph electrical work",
                    ));
                }
                Ok(resolved.bulk)
            },
            |sample| &sample.material,
            HeatNetworkSample::radiation_net_power_w,
            &cancel,
        )
    }
}
