//! Spatial current and complete sourced heat transport on one physical endpoint.
//! The prepared heat assembly supplies the SAME cells/cards used by the field
//! material refresh. Its original interface and radiosity owners supply heat;
//! the original spatial electrical solve alone supplies Joule dissipation.

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
    /// Joint spatial current, material-temperature/phase feedback, thermal
    /// contact resistance and closed-enclosure radiation.
    ///
    /// The field's actual E.J integral heats the explicitly mapped original
    /// cells once. At every nonlinear trial the existing prepared heat network
    /// re-evaluates pinned bulk/interface/emissivity sources, and the existing
    /// field owner re-evaluates its pinned scalar/oriented conductivity at that
    /// SAME thermal state. No graph resistor, second timeline or heat law is
    /// introduced. Current-controlled electrodes still solve their voltage
    /// from the full spatial operator on every trial.
    ///
    /// `regions` and `volume_tolerance_m3` have the original
    /// [`ConductionThermalNetwork`] whole-volume semantics. Material/caloric
    /// identities, fixed coordinates, constant phase-chart density, and source
    /// support retain the original strict checks. The heat assembly must have
    /// no electrical graph bindings: those require their own circuit owner.
    /// A thermal-interface resistance is not an electrical contact resistance.
    ///
    /// Geometry, enclosure view factors and contact/radiating areas remain
    /// fixed. This does not model material motion, volume-changing melting or
    /// evolving interface topology. In particular a non-solid phase chart must
    /// meet the same explicitly confined constant-density assumptions as the
    /// original spatial conduction lane; no changing surface is inferred.
    ///
    /// The returned report retains the complete material/contact/radiation
    /// samples. Closed-enclosure radiation is internal exchange: its signed
    /// numerical net is separately checked against the ORIGINAL total-energy
    /// allowance. The physically external non-electrical heat is
    /// `report.non_electrical_heat_j - dt * report.source.radiation_net_power_w()`.
    /// First-order backward-Euler coupling requires time/spatial refinement;
    /// this composition does not establish experimental material validation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_conduction_candidate<'t, C: PlanarConductorMaterial>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &PlanarConduction<C>,
        drive: ConductionDrive<'_>,
        regions: &[ConductionRegion],
        volume_tolerance_m3: f64,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<ConductionEnthalpyCandidate<'t, C, HeatNetworkSample>, Error> {
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        check(&cancel)?;
        let bulk = self.thermal_material_network().map_err(heat_error)?;
        let layout = ConductionThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            regions,
            volume_tolerance_m3,
        };
        let next = layout.step_resolved_candidate(
            thermal,
            target_s,
            problem,
            drive,
            assembly,
            field_controls,
            controls,
            initial_energy_increments_j,
            |trial| {
                let coefficients = self.resolve(cx, trial).map_err(heat_error)?;
                // The original assembly has already excluded graph bindings;
                // preserve a fail-closed boundary if another owner is extended.
                if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                    return Err(invalid(
                        "spatial current cannot silently omit graph electrical work",
                    ));
                }
                Ok(coefficients.bulk)
            },
            |sample| &sample.material,
            &cancel,
        )?;
        let report = next.report();
        closed_balance(
            target_s - thermal.time_s(),
            report.total_energy_residual_j,
            &report.source,
            controls.total_energy_tolerance_j,
        )
        .map_err(heat_error)?;
        check(&cancel)?;
        Ok(next)
    }
}
