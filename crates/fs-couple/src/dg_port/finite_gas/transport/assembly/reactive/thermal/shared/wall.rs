//! Finite thermal storage pays the original fixed-geometry wall material work.
//!
//! The native wall owner alone transfers displacement, momentum and acoustic
//! histories. This adapter resolves its same material at the shared body's
//! endpoint and admits the finite thermal inventory that supplies that work.

use super::{
    EquilibriumPhaseState, Error, FiniteGasAcousticNode, ResolvedMaterialStatePoint,
    SharedFiniteGasAcousticThermalPort, Sum, UniformEnthalpyStepInput, poll, positive,
};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::wall::material::{
    FiniteGasDgBinding, FiniteGasDgError, FiniteGasDgTransportError, MaterialWallDg,
};
use fs_duct::dg::DgWallMaterialUpdate;
use fs_material::state_point::DENSITY_PROPERTY;
use std::collections::BTreeMap;

/// Bind an original uniform isotropic elastic wall to this thermal body.
/// Geometry, density, mass, restraint and source claims remain fixed. The native
/// displacement/momentum transfer supplies the signed work charged to the body.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedFiniteGasAcousticWallSolve {
    /// Absolute allowance for this original wall's signed material work [J].
    pub maximum_absolute_parameter_work_j: f64,
}

/// The actual endpoint material and native zero-duration wall intervention.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedFiniteGasAcousticWallStep {
    pub material: ResolvedMaterialStatePoint,
    pub native: DgWallMaterialUpdate,
}

fn native_error(parcel: u64, error: AcousticRealizeError) -> Error {
    match error {
        AcousticRealizeError::Cancelled => Error::Cancelled,
        error => Error::Native(FiniteGasDgTransportError::Acoustic {
            parcel,
            error: FiniteGasDgError::Acoustic(error),
        }),
    }
}

fn uniform(
    nodes: &BTreeMap<u64, FiniteGasAcousticNode>,
    id: u64,
) -> Result<&MaterialWallDg, Error> {
    let Some(FiniteGasAcousticNode::Wall(node)) = nodes.get(&id) else {
        return Err(Error::Invalid(
            "thermal wall work requires an original wall-DG owner",
        ));
    };
    let FiniteGasDgBinding::Uniform(wall) = &node.binding else {
        return Err(Error::Invalid(
            "shared-body wall feedback currently requires one uniform wall source",
        ));
    };
    Ok(wall)
}

fn same_density(a: f64, b: f64) -> bool {
    positive(a) && positive(b) && (a / b - 1.0).abs() <= 64.0 * f64::EPSILON
}

fn within_inventory(amount: f64, available: f64) -> bool {
    amount.is_finite()
        && amount >= 0.0
        && (amount <= available || amount / available - 1.0 <= 64.0 * f64::EPSILON)
}

fn admit_phase(state: EquilibriumPhaseState, density: f64) -> Result<(), Error> {
    if state.liquid_mass_fraction() != 0.0 || !same_density(state.bulk_density_kg_m3(), density) {
        return Err(Error::Invalid(
            "fixed-geometry thermal wall feedback requires an all-solid constant-density body",
        ));
    }
    Ok(())
}

/// Check every selected original wall before the first transport/source sample.
/// Returns false without changing the pre-existing no-wall thermal path.
pub(super) fn admit(
    nodes: &BTreeMap<u64, FiniteGasAcousticNode>,
    ports: &BTreeMap<u64, SharedFiniteGasAcousticThermalPort<'_>>,
    body: UniformEnthalpyStepInput<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<bool, Error> {
    let mut selected = false;
    let mut volume = Sum::default();
    let mut mass = Sum::default();
    let body_density = body.mass_kg / body.volume_m3;
    for (&id, port) in ports {
        let Some(controls) = port.wall else { continue };
        poll(cancelled)?;
        selected = true;
        admit_phase(body.initial, body_density)?;
        if !controls.maximum_absolute_parameter_work_j.is_finite()
            || controls.maximum_absolute_parameter_work_j < 0.0
        {
            return Err(Error::Invalid(
                "thermal wall feedback requires a finite nonnegative native work allowance",
            ));
        }
        let wall = uniform(nodes, id)?;
        // The same helper checks initial temperature, card, isotropic/elastic
        // model, complete pinned query bundle and fixed material density.
        let original = wall
            .wall()
            .at_temperature_pinned(
                port.layers.material,
                port.layers.temperature_axis,
                body.initial.temperature_k(),
                body.initial.temperature_k(),
                cancelled,
            )
            .map_err(|error| native_error(id, error))?;
        let density = original
            .material()
            .property(DENSITY_PROPERTY)
            .ok_or(Error::Invalid(
                "original thermal wall has no density source",
            ))?
            .value_si();
        if !same_density(density, body_density) {
            return Err(Error::Invalid(
                "thermal inventory density must match its original modeled wall material",
            ));
        }
        let radius = wall
            .duct()
            .uniform_circular_radius_m()
            .ok_or(Error::Invalid(
                "thermal wall work requires its original circular fixed geometry",
            ))?;
        // This is the same thin-membrane measure used by native surface inertia,
        // not a second solid geometry or a changed acoustic volume.
        let wall_volume =
            core::f64::consts::TAU * radius * wall.duct().total_length() * original.thickness_m();
        let wall_mass = wall_volume * density;
        if !positive(wall_volume) || !positive(wall_mass) {
            return Err(Error::Invalid("thermal wall inventory is unrepresentable"));
        }
        volume.add(wall_volume)?;
        mass.add(wall_mass)?;
        if !within_inventory(volume.value, body.volume_m3)
            || !within_inventory(mass.value, body.mass_kg)
        {
            return Err(Error::Invalid(
                "shared thermal body cannot contain less material than its selected native walls",
            ));
        }
    }
    poll(cancelled)?;
    Ok(selected)
}

pub(super) struct Updates {
    pub steps: BTreeMap<u64, SharedFiniteGasAcousticWallStep>,
    pub parameter_work_j: f64,
    pub final_acoustic_storage_change_j: f64,
}

/// Apply material work only on private post-transport/post-gas-refresh owners.
/// Every trial resolves from the original source and accepted body temperature.
pub(super) fn apply(
    original_nodes: &BTreeMap<u64, FiniteGasAcousticNode>,
    next_nodes: &mut BTreeMap<u64, FiniteGasAcousticNode>,
    ports: &BTreeMap<u64, SharedFiniteGasAcousticThermalPort<'_>>,
    body: UniformEnthalpyStepInput<'_>,
    state: EquilibriumPhaseState,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Updates, Error> {
    poll(cancelled)?;
    admit_phase(state, body.mass_kg / body.volume_m3)?;
    let mut steps = BTreeMap::new();
    let mut work = Sum::default();
    for (&id, port) in ports {
        let Some(controls) = port.wall else { continue };
        poll(cancelled)?;
        let original = uniform(original_nodes, id)?;
        let replacement = original
            .wall()
            .at_temperature_pinned(
                port.layers.material,
                port.layers.temperature_axis,
                body.initial.temperature_k(),
                state.temperature_k(),
                cancelled,
            )
            .map_err(|error| native_error(id, error))?;
        let material = replacement.material().clone();
        let Some(FiniteGasAcousticNode::Wall(node)) = next_nodes.get_mut(&id) else {
            return Err(Error::Invalid(
                "thermal trial lost its original wall-DG owner",
            ));
        };
        let FiniteGasDgBinding::Uniform(wall) = &mut node.binding else {
            return Err(Error::Invalid(
                "thermal trial changed its original wall source layout",
            ));
        };
        if wall.wall() != original.wall() {
            return Err(Error::Invalid(
                "gas trial changed the thermal wall's original material source",
            ));
        }
        let native = wall
            .rebind_controlled(
                replacement,
                &mut node.state,
                &mut node.propagator,
                controls.maximum_absolute_parameter_work_j,
                cancelled,
            )
            .map_err(|error| native_error(id, error))?;
        work.add(native.parameter_work_j)?;
        steps.insert(id, SharedFiniteGasAcousticWallStep { material, native });
    }
    // Measure the published endpoint independently of the parameter-work sum.
    let mut storage = Sum::default();
    for (&id, node) in next_nodes.iter() {
        poll(cancelled)?;
        let previous = original_nodes.get(&id).ok_or(Error::Invalid(
            "thermal wall endpoint lost an original acoustic identity",
        ))?;
        storage.add(node.stored_energy_j() - previous.stored_energy_j())?;
    }
    poll(cancelled)?;
    Ok(Updates {
        steps,
        parameter_work_j: work.value,
        final_acoustic_storage_change_j: storage.value,
    })
}
