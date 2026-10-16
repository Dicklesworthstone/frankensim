//! Both gas contacts use the same material-network storage and transport sources.
use super::*;
use crate::circuit_enthalpy::material::{MaterialNetwork, MaterialNetworkError, MaterialNetworkSample};
use crate::cavity_phs::piston::thermal::boundary::network::material::{GasConductionGeometry, resolve_contact};
use fs_material::state_point::ResolvedMaterialStatePoint;

/// Complete source-bound network and paired-gas proposal.
pub type MaterialDoubleNetworkCandidate<'m> = DoubleNetworkCandidate<'m, MaterialNetworkSample, ResolvedMaterialStatePoint>;
/// Original material, boundary and numerical failure with no partial endpoint.
pub type MaterialDoubleNetworkError = CircuitEnthalpyError<NetworkGasError<MaterialNetworkError, MaterialNetworkError>>;

impl DoubleActingGasPiston {
    /// Compose the original electrical, caloric and conductivity material bindings
    /// with both gas ports. Only material data and contact geometry are inputs;
    /// there is no handwritten constitutive callback, copied coefficient or
    /// independently evolved opposing wall. Current resistance and gas-contact
    /// conductivity see the same cell temperature/phase/source evaluation.
    ///
    /// Patch areas are summed per cell, including when both faces name that cell.
    /// They must fit its exposed surface after native contact/reservoir coverage.
    /// Ambient convection on a gas-contact cell must be zero, so the occupied
    /// patch is never also counted as an independent ambient boundary. Missing
    /// conductivity still refuses an isolated patch. Other cells retain their
    /// original reservoirs, sources and thermal links.
    #[allow(clippy::too_many_arguments)]
    pub fn material_network_candidate<'m>(
        &self, binding: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
        circuit: &MaterialCircuit, ports: &[f64], drive: DoubleNetworkDrive, controls: NetworkGasSolve,
        geometry: [GasConductionGeometry; 2], guess_j: &[f64], cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialDoubleNetworkCandidate<'m>, MaterialDoubleNetworkError> {
        let cancel = RefCell::new(cancelled);
        if (*cancel.borrow_mut())() { return Err(CircuitEnthalpyError::Cancelled); }
        if binding.resistors.len() > crate::material_circuit::MAX_MATERIAL_CIRCUIT_BRANCHES {
            return Err(CircuitEnthalpyError::Invalid("material circuit binding exceeds capacity"));
        }
        let mut coverage = BTreeMap::new();
        for (cell, patch) in drive.cells.into_iter().zip(geometry) {
            if (*cancel.borrow_mut())() { return Err(CircuitEnthalpyError::Cancelled); }
            if !patch.area_m2.is_finite() || patch.area_m2 < 0.0
                || !positive(patch.gas_distance_m) || !positive(patch.cell_distance_m) {
                return Err(CircuitEnthalpyError::Invalid("invalid paired gas-contact geometry"));
            }
            let area = coverage.entry(cell).or_insert(0.0);
            *area += patch.area_m2;
            if !area.is_finite() || *area > thermal.exposed_area_m2(cell)? {
                return Err(CircuitEnthalpyError::Invalid("paired gas patches exceed shared cell exposure"));
            }
        }
        let destinations: Vec<_> = binding.resistors.iter().map(|r| r.destination).collect();
        let result = self.network_candidate(thermal, target_s,
            CircuitEnthalpyInput { circuit, ports, destinations: &destinations }, drive, controls, guess_j,
            |state| binding.resolve(state, || (*cancel.borrow_mut())()),
            |state, source, gases| Ok([
                resolve_contact(binding, drive.cells[0], state, source, &gases[0], geometry[0])?,
                resolve_contact(binding, drive.cells[1], state, source, &gases[1], geometry[1])?,
            ]), || (*cancel.borrow_mut())());
        match result {
            Err(CircuitEnthalpyError::Source(NetworkGasError::Material(MaterialNetworkError::Thermal(DepositionError::Cancelled)))) =>
                Err(CircuitEnthalpyError::Cancelled),
            other => other,
        }
    }
}
