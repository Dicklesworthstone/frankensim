//! Data-only gas surface bindings on the existing material/circuit network.
use super::*;
use crate::circuit_enthalpy::material::{MaterialNetwork, MaterialNetworkError, MaterialNetworkSample};
use crate::cavity_phs::piston::thermal::boundary::network::material::{GasConductionGeometry, resolve_contact};
use fs_material::state_point::ResolvedMaterialStatePoint;

/// One non-overlapping patch of a material-bound thermal cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialGasSurface {
    /// Stable surface identity and its existing thermal cell.
    pub port: NetworkThermalPort,
    /// Current effective thermal path, distinct from electrical geometry.
    pub geometry: GasConductionGeometry,
}

/// One common gas/momentum input plus all material surfaces in contact with it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialSurfacesDrive<'a> {
    /// Complete bounded list of material patches. Non-overlap is caller supplied;
    /// the adapter additionally checks TOTAL covered area per receiving cell.
    pub surfaces: &'a [MaterialGasSurface],
    /// Existing gas sources and mechanical loads, held for this interval.
    pub input: PistonInput,
    /// Original pressure-work/momentum controls.
    pub piston: PistonSolve,
    /// Net gas heat feedback controls, not individual surface heat limits.
    pub exchange: HeatExchangeSolve,
}

impl GasPiston {
    /// Bind multiple finite material surfaces to the SAME gas/pressure-work step.
    /// Each uses the existing two-layer conduction law, at its own cell endpoint
    /// and the one common gas midpoint. Material and original source receipts
    /// must agree with the network. No independent wall store or material preset.
    ///
    /// Sum of patch areas on EACH cell must fit its remaining exposed surface,
    /// after existing links/reservoirs. Duplicate ambient convection on contacted
    /// cells refuses. This scalar area check does not prove spatial non-overlap.
    /// Missing conductivity refuses even at zero area or zero net gas heating.
    /// Other network cells retain their original authored boundaries and sources.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn material_network_surfaces_candidate<'m>(
        &self, binding: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
        circuit: &MaterialCircuit, electrical_ports: &[f64], drive: MaterialSurfacesDrive<'_>,
        controls: NetworkGasSolve, guess_j: &[f64], cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkSurfacesCandidate<'m, MaterialNetworkSample, ResolvedMaterialStatePoint>,
        CircuitEnthalpyError<NetworkGasError<MaterialNetworkError, MaterialNetworkError>>> {
        let cancel = RefCell::new(cancelled);
        if (*cancel.borrow_mut())() { return Err(CircuitEnthalpyError::Cancelled); }
        if drive.surfaces.is_empty() || drive.surfaces.len() > MAX_NETWORK_SURFACES
            || binding.resistors.len() > crate::material_circuit::MAX_MATERIAL_CIRCUIT_BRANCHES {
            return Err(CircuitEnthalpyError::Invalid("bounded complete material surfaces and resistor bindings required"));
        }
        let mut areas = BTreeMap::new(); let mut surfaces = BTreeMap::new();
        for surface in drive.surfaces {
            if (*cancel.borrow_mut())() { return Err(CircuitEnthalpyError::Cancelled); }
            let geometry = surface.geometry;
            if surfaces.insert(surface.port.id, surface).is_some() || !geometry.area_m2.is_finite()
                || geometry.area_m2 < 0.0 || !positive(geometry.gas_distance_m) || !positive(geometry.cell_distance_m) {
                return Err(CircuitEnthalpyError::Invalid("unique surface identities and admitted conduction geometry required"));
            }
            let cell = thermal.cells().get(&surface.port.cell)
                .ok_or(CircuitEnthalpyError::Invalid("inactive gas-contact cell"))?;
            let area = areas.entry(surface.port.cell).or_insert(0.0);
            *area += geometry.area_m2;
            if !area.is_finite() || *area > thermal.exposed_area_m2(surface.port.cell)?
                || cell.boundary.htc_w_m2_k != 0.0 {
                return Err(CircuitEnthalpyError::Invalid("summed gas patches exceed exposed area or duplicate ambient convection"));
            }
        }
        let ports: Vec<_> = surfaces.values().map(|s| s.port).collect();
        let destinations: Vec<_> = binding.resistors.iter().map(|r| r.destination).collect();
        let result = self.network_surfaces_candidate(thermal, target_s,
            CircuitEnthalpyInput { circuit, ports: electrical_ports, destinations: &destinations },
            NetworkSurfacesDrive { ports: &ports, input: drive.input, piston: drive.piston, exchange: drive.exchange },
            controls, guess_j, |state| binding.resolve(state, || (*cancel.borrow_mut())()),
            |port, state, source, gas| {
                resolve_contact(binding, port.cell, state, source, gas, surfaces[&port.id].geometry)
            }, || (*cancel.borrow_mut())());
        match result {
            Err(CircuitEnthalpyError::Source(NetworkGasError::Material(MaterialNetworkError::Thermal(DepositionError::Cancelled)))) =>
                Err(CircuitEnthalpyError::Cancelled),
            other => other,
        }
    }
}
