//! One original material binding for an open gas and finite circuit/thermal network.

use super::*;
use crate::cavity_phs::piston::thermal::boundary::network::material::GasConductionGeometry;
use crate::circuit_enthalpy::material::{MaterialNetwork, MaterialNetworkError, MaterialNetworkSample};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::state_point::ResolvedMaterialStatePoint;

impl GasPiston {
    /// Data-only nozzle, electrical/thermal network and material wall coupling.
    ///
    /// Each circuit branch declares whether its source property is conductivity
    /// or resistivity. The existing prepared network pins those source claims
    /// while reevaluating values at trial temperatures. The wall interface uses
    /// that SAME endpoint bundle, and gas conductivity comes from the transferred
    /// mixture's actual midpoint. No second material resolution or nominal gas
    /// coefficient is supplied here.
    ///
    /// Original patch coverage and ambient-convection admission still apply.
    /// The thermal cells can change temperature/phase; the translating piston
    /// remains a rigid body. See `nozzle_network_candidate` for split accuracy,
    /// source accounting, atomic publication and the admitted physical domain.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn nozzle_material_network_candidate<'m>(
        &self, reservoir: &ClosedGasParcel, nozzle: GasNozzle,
        binding: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
        circuit: &MaterialCircuit, ports: &[f64],
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        drive: NetworkGasDrive, flow: GasNozzleStepSolve, controls: NetworkGasSolve,
        geometry: GasConductionGeometry, guess_j: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NozzleNetworkCandidate<'m, MaterialNetworkSample, ResolvedMaterialStatePoint>,
        CircuitEnthalpyError<NetworkGasError<MaterialNetworkError, MaterialNetworkError>>>
    {
        let next = self.nozzle_network_transfer(reservoir, nozzle, thermal, target_s,
            drive, flow, controls, &mut cancelled)?;
        let network = next.piston.material_network_candidate_with_bases(binding, thermal,
            target_s, circuit, ports, electrical_bases, drive, controls, geometry, guess_j,
            &mut cancelled)?;
        self.finish_nozzle_network(reservoir, next, network, drive, flow, controls, &mut cancelled)
    }
}
