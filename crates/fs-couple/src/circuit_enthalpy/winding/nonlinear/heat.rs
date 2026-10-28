//! Sourced thermal interfaces and radiation share the native winding residual.
//! Source selection, current/field history and finite heat remain with their
//! original owners. The only additional work here is their physical composition.

use super::*;
use crate::circuit_enthalpy::{
    ResistorCell,
    contacts::heat::{HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork, closed_balance},
    finite,
    winding::MagneticRegion,
};
use fs_exec::Cx;

fn heat_error(error: HeatNetworkError) -> Error {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(WindingThermalError::Heat(other)),
    }
}

impl PreparedHeatNetwork<'_, '_> {
    /// Evolve native winding/core currents and finite material heat together
    /// with sourced thermal-contact resistance and closed-enclosure radiation.
    ///
    /// Every nonlinear trial uses this prepared assembly's original cards,
    /// pinned claims and actual endpoint temperatures. The original field owns
    /// saturation and eddy currents; its physical winding/core/graph losses
    /// enter the same enthalpy residual exactly once. Numerical dissipation is
    /// retained in the balance and never deposited as material heat.
    ///
    /// `magnetic_thermal` enables the existing reversible magnetic caloric
    /// owner: its source-bound affine permeability, field entropy, internal
    /// energy and positive combined heat-capacity checks are unchanged. None
    /// retains the original fixed magnetic law, including its saturation law;
    /// it does not freeze or silently discard an authored temperature response.
    ///
    /// `graph_resistors` routes EVERY native ordinary resistor's actual heat.
    /// If this assembly binds their materials, all graph resistors must be
    /// covered, with exactly matching destinations and separate whole-cell
    /// volumes. Otherwise their native authored resistance remains fixed.
    /// Electrical series contacts are refused: thermal interface resistance
    /// alone does not supply their voltage jump, work or electrical heat split.
    ///
    /// Geometry, magnetic-region/cell correspondence, contact areas and enclosure
    /// view factors remain fixed. Mapped magnetic cells retain their original
    /// wholly solid constant-density charts. Time, field and thermal spatial
    /// refinement remain necessary; this is the existing first-order coupled
    /// approximation, not a measured-device or unconditional-stability claim.
    ///
    /// Both free-energy and physical magnetic/caloric balances exclude the
    /// closed enclosure's roundoff net from external supply. The retained
    /// caloric balance reports that corrected physical external heat directly.
    /// The ordinary report's external non-electrical heat is
    /// `non_electrical_heat_j - reversible_heat_j - dt*net_radiation_power_w`.
    /// All state remains unchanged unless the complete candidate is accepted.
    #[allow(clippy::too_many_arguments)]
    pub fn step_nonlinear_winding_candidate<'t, 'm>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &NonlinearWindingCircuit<'m>,
        external_inputs: &[f64],
        magnetic_regions: &[MagneticRegion],
        graph_resistors: &[ResistorCell],
        volume_tolerance_m3: f64,
        magnetic_thermal: Option<(&[MagneticThermalLaw], WeightedAssemblyLimits)>,
        controls: CircuitEnthalpySolve,
        winding_controls: NonlinearWindingControls,
        initial_energy_increments_j: &[f64],
    ) -> Result<NonlinearWindingEnthalpyCandidate<'t, 'm, HeatNetworkSample>, Error> {
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        check(&cancel)?;
        let bulk = self.field_material_network().map_err(heat_error)?;
        let layout = WindingThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            magnetic_regions,
            graph_resistors,
            volume_tolerance_m3,
        };
        let mut next = layout.step_nonlinear_resolved(
            thermal,
            target_s,
            circuit,
            external_inputs,
            (!bulk.resistors.is_empty()).then_some(bulk.resistors),
            magnetic_thermal,
            controls,
            winding_controls,
            initial_energy_increments_j,
            |state| {
                let coefficients = self.resolve(cx, state).map_err(heat_error)?;
                if !coefficients.contacts.is_empty() {
                    return Err(invalid(
                        "native winding cannot discard electrical contact work",
                    ));
                }
                Ok(coefficients.bulk)
            },
            |sample| &sample.material,
            &cancel,
        )?;
        let dt = target_s - thermal.time_s();
        let report = &mut next.report;
        let sample = &report.source.network;
        closed_balance(
            dt,
            report.total_energy_residual_j,
            sample,
            controls.total_energy_tolerance_j,
        )
        .map_err(heat_error)?;
        if let Some(balance) = &mut report.source.magnetic_balance {
            closed_balance(
                dt,
                balance.total_energy_residual_j,
                sample,
                controls.total_energy_tolerance_j,
            )
            .map_err(heat_error)?;
            let internal_net = finite(dt * sample.radiation_net_power_w())?;
            balance.external_heat_j = finite(balance.external_heat_j - internal_net)?;
            balance.total_energy_residual_j =
                finite(balance.total_energy_residual_j + internal_net)?;
        }
        check(&cancel)?;
        Ok(next)
    }
}
