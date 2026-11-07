//! Original material R(T) on the same solid-electrode/field/enthalpy solve.

use super::super::terminals::circuit::{SpatialTerminalCircuitThermalNetwork, TerminalError};
use super::*;
use crate::circuit_enthalpy::material::CellResistor;
use fs_material::conductor::ElectricalPropertyBasis;

impl SpatialTerminalCircuitThermalNetwork<'_> {
    /// Advance the linear A–phi field, original RLC graph and finite heat while
    /// each authored graph resistor follows its own original material R(T).
    ///
    /// The explicit graph destinations, source selections, electrical bases
    /// and uniform path geometries must describe complete, separate caloric
    /// resistor bodies. Accepted source state and resistance must agree before
    /// stepping. Every thermal trial then refreshes those same source choices,
    /// conductor sigma(T), and the admitted magnetic mu(T)/entropy law without
    /// resetting magnetic potential, native flux/charge or caloric history.
    ///
    /// Resistance-only native reassembly preserves circuit storage and creates
    /// no stored-energy parameter work. Original resistor loss and conductor
    /// Joule heat enter their finite receivers once; reversible magnetic heat
    /// retains its own original caloric account. Geometry, topology, L and C
    /// remain fixed. Prescribed conductor motion exchanges external work.
    #[allow(clippy::too_many_arguments)]
    pub fn step_linear_material_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialTerminalCircuit<'m, C>,
        external_inputs: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        resistors: &[CellResistor],
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        magnetic_laws: &[SpatialMagneticThermalLaw],
        controls: CircuitEnthalpySolve,
        circuit_controls: SpatialTerminalCircuitControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialCircuitThermalCandidate<'t, 'm, C>, TerminalError> {
        let cancel = RefCell::new(cancelled);
        check(&cancel).map_err(native_error)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        let network = MaterialNetwork {
            cells: self.magnetic.cells,
            resistors,
            links: self.magnetic.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, electrical_bases, || (*cancel.borrow_mut())())
            .map_err(material)
            .map_err(native_error)?;
        self.step_linear_resolved_candidate(
            thermal,
            target_s,
            circuit,
            external_inputs,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            Some(resistors),
            magnetic_laws,
            controls,
            circuit_controls,
            assembly,
            initial_energy_increments_j,
            |state| {
                prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)
                    .map_err(native_error)
            },
            |sample| sample,
            &cancel,
        )
    }
}
