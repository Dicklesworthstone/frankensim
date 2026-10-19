//! Data-only gas contact on the original material/circuit/enthalpy binding.
use super::*;
use crate::circuit_enthalpy::material::{MaterialNetwork,MaterialNetworkError,MaterialNetworkSample};
use crate::cavity_phs::piston::thermal::layers::MaterialConductionLayers;
use fs_material::state_point::ResolvedMaterialStatePoint;

/// Explicit uniform conduction geometry, distinct from electrical current paths.
#[derive(Clone,Copy,Debug,PartialEq)]
pub struct GasConductionGeometry {
    /// Patch on the selected cell's exposed surface [m2]; zero isolates it.
    pub area_m2:f64,
    /// Positive effective gas layer distance [m].
    pub gas_distance_m:f64,
    /// Positive effective cell layer distance [m].
    pub cell_distance_m:f64,
}

impl GasPiston {
    /// Couple a finite source-bound circuit/thermal network through one exposed
    /// conductive cell patch. Resistance, storage and both sides' conductivity
    /// come from their original material owners at the actual trial states.
    /// No nominal gas, material-name switch or source callback is needed.
    ///
    /// The cell's ambient convection must be zero: this API does not silently
    /// count the gas patch again as exposed ambient. The patch must fit the
    /// remaining exposed area after existing links and reservoir coverage.
    /// Other cells may retain arbitrary admitted external boundaries. A spatial
    /// convection or gas-wall boundary-layer calculation is not implied.
    #[allow(clippy::too_many_arguments,clippy::type_complexity)]
    pub fn material_network_candidate<'m>(
        &self,binding:&MaterialNetwork<'_>,thermal:&DepositionSimulation<'m>,target_s:f64,
        circuit:&MaterialCircuit,ports:&[f64],drive:NetworkGasDrive,controls:NetworkGasSolve,
        geometry:GasConductionGeometry,guess_j:&[f64],cancelled:impl FnMut()->bool,
    )->Result<NetworkGasCandidate<'m,MaterialNetworkSample,ResolvedMaterialStatePoint>,
        CircuitEnthalpyError<NetworkGasError<MaterialNetworkError,MaterialNetworkError>>> {
        let cancel=RefCell::new(cancelled);
        if (*cancel.borrow_mut())() {return Err(CircuitEnthalpyError::Cancelled);}
        if binding.resistors.len()>crate::material_circuit::MAX_MATERIAL_CIRCUIT_BRANCHES {
            return Err(CircuitEnthalpyError::Invalid("material circuit binding exceeds capacity"));
        }
        let destinations:Vec<_>=binding.resistors.iter().map(|r|r.destination).collect();
        let result=self.network_candidate(thermal,target_s,CircuitEnthalpyInput {circuit,ports,destinations:&destinations},
            drive,controls,guess_j,|state|binding.resolve(state,||(*cancel.borrow_mut())()),
            |state,source,gas| resolve_contact(binding,drive.cell,state,source,gas,geometry),
            ||(*cancel.borrow_mut())());
        match result {
            Err(CircuitEnthalpyError::Source(NetworkGasError::Material(MaterialNetworkError::Thermal(DepositionError::Cancelled)))) =>
                Err(CircuitEnthalpyError::Cancelled),
            other=>other,
        }
    }
}

// The same source and conduction admission for single- and multi-chamber users.
// Multi-contact callers ALSO check aggregate coverage before entering a solve.
pub(crate) fn resolve_contact(
    binding:&MaterialNetwork<'_>,id:CellId,state:&DepositionSimulation<'_>,
    source:&MaterialNetworkSample,gas:&ResolvedThermochemicalGasState,geometry:GasConductionGeometry,
) -> Result<ConductanceSample<ResolvedMaterialStatePoint>,MaterialNetworkError> {
    let invalid=MaterialNetworkError::Invalid;
    let cell=state.cells().get(&id).ok_or(invalid("inactive gas contact cell"))?;
    if cell.boundary.htc_w_m2_k!=0.0 || !geometry.area_m2.is_finite() || geometry.area_m2<0.0
        || geometry.area_m2>state.exposed_area_m2(id)? {
        return Err(invalid("gas contact needs an exposed patch without duplicate ambient convection"));
    }
    let material=binding.cells.iter().find(|c|c.cell==id)
        .ok_or(invalid("gas contact lacks original cell material binding"))?;
    let sample=source.cells.iter().find(|c|c.cell==id)
        .ok_or(invalid("gas contact lacks resolved endpoint phase"))?;
    let (key,selection)=material.conductivity.as_ref()
        .ok_or(invalid("gas contact requires explicit cell conductivity"))?;
    let layers=MaterialConductionLayers {material:material.point.material,context:material.point.context,
        temperature_axis:material.point.temperature_axis,liquid_fraction_axis:material.point.liquid_fraction_axis,
        conductivity:key,selection:selection.clone(),area_m2:geometry.area_m2,
        gas_distance_m:geometry.gas_distance_m,body_distance_m:geometry.cell_distance_m};
    // Reuse the complete endpoint bundle, including electrical properties.
    // Re-querying just k creates a different bundle and rejects valid coupled
    // states. The layer owner still checks card, coordinates, key and source pin.
    let snapshot=sample.conductivity.as_ref()
        .ok_or(invalid("gas contact lacks resolved endpoint conductivity"))?;
    layers.resolve_from_material(gas,sample.phase,snapshot.clone())
        .map_err(MaterialNetworkError::Material)
}
