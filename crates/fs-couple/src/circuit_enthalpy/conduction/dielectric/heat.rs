//! Complete material heat transport shares the retained spatial charge owner.
//!
//! No temperature is prescribed after an electrical step: the original
//! interface, integrated-layer and radiosity owners resolve inside the same
//! nonlinear enthalpy residual as the original charge/displacement-current law.

use super::thermal::DielectricCaloricCandidate;
use super::*;
use crate::circuit_enthalpy::contacts::heat::{
    HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork, closed_balance,
};
use fs_exec::Cx;

fn heat_error(error: HeatNetworkError) -> Error {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ConductionThermalError::Heat(other)),
    }
}

impl PreparedHeatNetwork<'_, '_> {
    /// Continue spatial charge relaxation with the complete material heat
    /// assembly: sourced thermal contacts, integrated conduction layers and
    /// closed-enclosure radiation. Scalar or oriented conductivity is resolved
    /// at the SAME endpoint as every thermal source; the native midpoint E.J
    /// heats its original cells once. Potential, nodal charge, enthalpy and all
    /// source receipts publish together after the original energy checks.
    ///
    /// `regions` maps every original electrical triangle to its complete fixed
    /// caloric volume. This constant-permittivity path retains the original
    /// equal-density confined phase-change support. Geometry, contact areas,
    /// dielectric frames and radiation view factors remain fixed. Nonconstant
    /// permittivity requires [`Self::step_dielectric_thermal_candidate`].
    ///
    /// Graph electrical branches and unowned field-dependent thermal boundary
    /// laws refuse: neither their electrical work nor their local interface
    /// field can be inferred by a spatial charge owner. The enclosure's computed
    /// net is internal exchange, independently checked against the same energy
    /// allowance. External heat is `report.non_electrical_heat_j - dt *
    /// report.source.radiation_net_power_w()`. Endpoint transport and thermal
    /// backward Euler retain their first-order coupling approximation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_dielectric_candidate<'t, C: PlanarConductorMaterial>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &PlanarLeakyDielectric<C>,
        endpoint_voltage_v: &[(usize, f64)],
        regions: &[ConductionRegion],
        volume_tolerance_m3: f64,
        assembly: WeightedAssemblyLimits,
        field_controls: DielectricRelaxationControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<DielectricEnthalpyCandidate<'t, C, HeatNetworkSample>, Error> {
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        check(&cancel)?;
        let bulk = self.thermal_material_network().map_err(heat_error)?;
        self.admit_field_independent_boundaries(cx, &[])
            .map_err(heat_error)?;
        let layout = ConductionThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            regions,
            volume_tolerance_m3,
        };
        let next = layout.step_dielectric_resolved_candidate(
            thermal,
            target_s,
            field,
            endpoint_voltage_v,
            assembly,
            field_controls,
            controls,
            initial_energy_increments_j,
            |trial| {
                let resolved = self.resolve(cx, trial).map_err(heat_error)?;
                if !resolved.contacts.is_empty() || !resolved.bulk.resistors.is_empty() {
                    return Err(invalid(
                        "spatial charge cannot silently omit graph electrical work",
                    ));
                }
                Ok(resolved.bulk)
            },
            |sample| &sample.material,
            &cancel,
        )?;
        closed_balance(
            target_s - thermal.time_s(),
            next.report().total_energy_residual_j,
            &next.report().source,
            controls.total_energy_tolerance_j,
        )
        .map_err(heat_error)?;
        check(&cancel)?;
        Ok(next)
    }

    /// Couple source-derived scalar/oriented epsilon(T), retained spatial
    /// charge and complete material heat transport in one native residual.
    ///
    /// The same finite temperatures determine conductivity, dielectric
    /// Helmholtz storage, field entropy, contact resistance, layer conduction
    /// and emissivity. Native Ohmic heat and signed reversible field transfer
    /// `Qrev = -P - delta(T*S)` enter their original cells separately and once.
    /// Total physical internal energy closes using dielectric U=psi+T*S and
    /// caloric enthalpy; the radiosity net is excluded from external heat and
    /// retained explicitly in `report().closed_radiation_net_j`.
    ///
    /// `maximum_displacement_c_m2` retains the original per-triangle field and
    /// combined-capacity admission. This affine-source lane requires fully
    /// solid constant-density caloric charts covering the thermal source
    /// intervals. No moving interface, phase polarization, dielectric dispersion
    /// or surface-charge law is inferred. Refusal or cancellation returns no
    /// partial charge, heat, source or observer history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_dielectric_thermal_candidate<'t, C: PlanarConductorMaterial>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &PlanarLeakyDielectric<C>,
        endpoint_voltage_v: &[(usize, f64)],
        maximum_displacement_c_m2: &[f64],
        regions: &[ConductionRegion],
        volume_tolerance_m3: f64,
        assembly: WeightedAssemblyLimits,
        field_controls: DielectricRelaxationControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<DielectricCaloricCandidate<'t, C, HeatNetworkSample>, Error> {
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        check(&cancel)?;
        let bulk = self.thermal_material_network().map_err(heat_error)?;
        self.admit_field_independent_boundaries(cx, &[])
            .map_err(heat_error)?;
        let layout = ConductionThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            regions,
            volume_tolerance_m3,
        };
        let next = layout.step_dielectric_thermal_resolved_candidate(
            thermal,
            target_s,
            field,
            endpoint_voltage_v,
            maximum_displacement_c_m2,
            assembly,
            field_controls,
            controls,
            initial_energy_increments_j,
            |trial| {
                let resolved = self.resolve(cx, trial).map_err(heat_error)?;
                if !resolved.contacts.is_empty() || !resolved.bulk.resistors.is_empty() {
                    return Err(invalid(
                        "spatial charge cannot silently omit graph electrical work",
                    ));
                }
                Ok(resolved.bulk)
            },
            |sample| &sample.material,
            HeatNetworkSample::radiation_net_power_w,
            &cancel,
        )?;
        // The native free-energy account also closes when internal enclosure
        // transfer is removed. The original driver already checked physical U.
        closed_balance(
            target_s - thermal.time_s(),
            next.report().coupled.total_energy_residual_j,
            &next.report().coupled.source.material,
            controls.total_energy_tolerance_j,
        )
        .map_err(heat_error)?;
        check(&cancel)?;
        Ok(next)
    }
}
