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

use crate::cavity_phs::piston::material::MaterialDoubleActingGasPiston;

/// Circuit, finite thermal network and original material-bound moving specimen.
/// No numeric-only mechanism escapes while its density/geometry source is lost.
pub struct MaterialPistonNetworkCandidate<'m> {
    thermal: DepositionSimulation<'m>,
    circuit: MaterialCircuit,
    piston: MaterialDoubleActingGasPiston,
    report: DoubleNetworkReport<MaterialNetworkSample, ResolvedMaterialStatePoint>,
}
impl<'m> MaterialPistonNetworkCandidate<'m> {
    /// Proposed original thermal cells, laws and pending-event history.
    pub const fn thermal(&self) -> &DepositionSimulation<'m> { &self.thermal }
    /// Proposed original finite electrical history.
    pub const fn circuit(&self) -> &MaterialCircuit { &self.circuit }
    /// Both gases and momentum with the original canonical body inertia/source.
    pub const fn piston(&self) -> &MaterialDoubleActingGasPiston { &self.piston }
    /// Unchanged actual circuit, thermal, gas and mechanical work accounts.
    pub const fn report(&self) -> &DoubleNetworkReport<MaterialNetworkSample, ResolvedMaterialStatePoint> { &self.report }
    /// Publish all histories together; no returned heat is an additional source.
    pub fn into_parts(self) -> (DepositionSimulation<'m>, MaterialCircuit, MaterialDoubleActingGasPiston,
        DoubleNetworkReport<MaterialNetworkSample, ResolvedMaterialStatePoint>) {
        (self.thermal, self.circuit, self.piston, self.report)
    }
}

impl MaterialDoubleActingGasPiston {
    /// Advance the material-bound specimen through the SAME circuit/thermal/gas
    /// solve. Density and geometry remain owned by its original ResolvedBody;
    /// there is no competing mass input and no source resampling on publication.
    ///
    /// Both gas endpoints, momentum, exposed areas and material records are
    /// carried into the existing material mechanism constructor unchanged. This
    /// rebinds ownership only: it performs no second timestep or constitutive
    /// update. Thermal cell temperatures do not silently heat, expand, melt or
    /// soften the rigid mechanical specimen, which may be a different body.
    ///
    /// Changing supported cards before constructing the material body or thermal
    /// network changes the ensuing coupled response through data. Replacing a
    /// live body's mass or chemical energy without accounting is not supported.
    #[allow(clippy::too_many_arguments)]
    pub fn material_network_candidate<'m>(
        &self, binding: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
        circuit: &MaterialCircuit, ports: &[f64], drive: DoubleNetworkDrive, controls: NetworkGasSolve,
        geometry: [GasConductionGeometry; 2], guess_j: &[f64], mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialPistonNetworkCandidate<'m>, MaterialDoubleNetworkError> {
        let candidate = self.piston().material_network_candidate(binding, thermal, target_s,
            circuit, ports, drive, controls, geometry, guess_j, &mut cancelled)?;
        let (thermal, circuit, next, report) = candidate.into_parts();
        let piston = Self::new(next.primary().gas().clone(), next.opposing().clone(), self.body().clone(),
            [next.primary().area_m2(), next.opposing_area_m2()], next.primary().velocity_m_s())
            .map_err(|e| -> MaterialDoubleNetworkError {
                CircuitEnthalpyError::Source(NetworkGasError::Boundary(ThermalPistonError::from(e)))
            })?;
        if cancelled() { return Err(CircuitEnthalpyError::Cancelled); }
        Ok(MaterialPistonNetworkCandidate { thermal, circuit, piston, report })
    }
}
