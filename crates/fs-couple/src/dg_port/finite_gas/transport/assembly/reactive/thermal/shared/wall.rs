//! Finite thermal storage receives original wall loss and pays material work.
//!
//! The native wall owner alone transfers displacement, momentum and acoustic
//! histories. This adapter resolves its same material at the shared body's
//! endpoint and admits the finite thermal inventory that supplies that work.

use super::{
    EquilibriumPhaseState, Error, FiniteGasAcousticInput, FiniteGasAcousticNode,
    FiniteGasTransportNode, ResolvedMaterialStatePoint, SharedFiniteGasAcousticThermalPort, Sum,
    UniformEnthalpyStepInput, poll, positive,
};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::wall::material::{
    FiniteGasDgBinding, FiniteGasDgError, FiniteGasDgNode, FiniteGasDgTransportError,
    LocalMaterialWall,
};
use fs_duct::dg::DgWallMaterialUpdate;
use fs_material::state_point::DENSITY_PROPERTY;
use fs_phs::linear_exact::DenseExponentialBudget;
use std::collections::BTreeMap;

/// Bind an original uniform or regional material wall to this thermal body.
/// Geometry, density, mass, restraint and source claims remain fixed. The native
/// displacement/momentum/history transfer supplies signed material work. Elastic,
/// Kelvin viscosity and source-pair Maxwell memory retain their original laws.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedFiniteGasAcousticWallSolve {
    /// Absolute allowance for this original wall's signed material work [J].
    pub maximum_absolute_parameter_work_j: f64,
    /// Opt in to physical Kelvin/Maxwell heat from the original native sample.
    /// The explicit complete-state source-integration budget is additional to
    /// propagation admission. None keeps the loss exported from this boundary.
    pub dissipation_heat: Option<DenseExponentialBudget>,
}

/// Original source layout at the common thermal endpoint. Regional receipts
/// retain original DG element order; no averaged or representative card is made.
#[derive(Clone, Debug, PartialEq)]
pub enum SharedFiniteGasAcousticWallMaterials {
    Uniform(ResolvedMaterialStatePoint),
    Regional(Vec<ResolvedMaterialStatePoint>),
}

/// Actual endpoint sources and the native zero-duration wall intervention.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedFiniteGasAcousticWallStep {
    pub materials: SharedFiniteGasAcousticWallMaterials,
    pub native: DgWallMaterialUpdate,
    /// Original solid-only source integral deposited in the finite body [J].
    /// Gas boundary-layer, molecular, terminal and numerical loss are excluded.
    pub dissipated_heat_j: f64,
}

impl SharedFiniteGasAcousticWallStep {
    /// The unchanged uniform source receipt, when the original binding is uniform.
    pub fn uniform_material(&self) -> Option<&ResolvedMaterialStatePoint> {
        match &self.materials {
            SharedFiniteGasAcousticWallMaterials::Uniform(material) => Some(material),
            SharedFiniteGasAcousticWallMaterials::Regional(_) => None,
        }
    }

    /// Every original element's source receipt, when the binding is regional.
    pub fn regional_materials(&self) -> Option<&[ResolvedMaterialStatePoint]> {
        match &self.materials {
            SharedFiniteGasAcousticWallMaterials::Uniform(_) => None,
            SharedFiniteGasAcousticWallMaterials::Regional(materials) => Some(materials),
        }
    }
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

fn wall_node(
    nodes: &BTreeMap<u64, FiniteGasAcousticNode>,
    id: u64,
) -> Result<&FiniteGasDgNode, Error> {
    let Some(FiniteGasAcousticNode::Wall(node)) = nodes.get(&id) else {
        return Err(Error::Invalid(
            "thermal wall work requires an original wall-DG owner",
        ));
    };
    Ok(node)
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

#[allow(clippy::too_many_arguments)]
fn admit_inventory(
    id: u64,
    source: &LocalMaterialWall,
    radius: f64,
    length: f64,
    port: &SharedFiniteGasAcousticThermalPort<'_>,
    body: UniformEnthalpyStepInput<'_>,
    volume: &mut Sum,
    mass: &mut Sum,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    // Each element uses its own original query and source claims. The one
    // thermal body still requires its original card and fixed density.
    let original = source
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
    if !same_density(density, body.mass_kg / body.volume_m3) {
        return Err(Error::Invalid(
            "thermal inventory density must match its original modeled wall material",
        ));
    }
    // This is the original native thin-membrane measure. The uniform branch
    // retains its former total-length arithmetic; regions use their own length.
    let wall_volume = core::f64::consts::TAU * radius * length * original.thickness_m();
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
        let node = wall_node(nodes, id)?;
        let duct = node.binding.duct();
        let radius = duct.uniform_circular_radius_m().ok_or(Error::Invalid(
            "thermal wall work requires its original circular fixed geometry",
        ))?;
        match &node.binding {
            FiniteGasDgBinding::Uniform(wall) => admit_inventory(
                id,
                wall.wall(),
                radius,
                duct.total_length(),
                port,
                body,
                &mut volume,
                &mut mass,
                cancelled,
            )?,
            FiniteGasDgBinding::Regional(wall) => {
                let positions = duct.node_positions();
                let local = duct.mesh().order + 1;
                if wall.walls().len() != duct.mesh().elements || wall.walls().is_empty() {
                    return Err(Error::Invalid(
                        "thermal wall regions must retain every original element",
                    ));
                }
                for (element, source) in wall.walls().iter().enumerate() {
                    poll(cancelled)?;
                    let base = element * local;
                    let length = positions[base + local - 1] - positions[base];
                    admit_inventory(
                        id,
                        source,
                        radius,
                        length,
                        port,
                        body,
                        &mut volume,
                        &mut mass,
                        cancelled,
                    )?;
                }
            }
        }
    }
    poll(cancelled)?;
    Ok(selected)
}

/// Each source is integrated once, before discarded body/contact/work trials.
/// No public caller can supply an arbitrary deposited wall heat.
pub(super) struct PreparedHeat {
    pub by_wall: BTreeMap<u64, f64>,
    pub total_j: f64,
}

pub(super) fn prepare_heat(
    nodes: &BTreeMap<u64, FiniteGasAcousticNode>,
    inputs: &BTreeMap<u64, FiniteGasAcousticInput>,
    ports: &BTreeMap<u64, SharedFiniteGasAcousticThermalPort<'_>>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedHeat, Error> {
    let mut by_wall = BTreeMap::new();
    let mut heat = Sum::default();
    for (&id, port) in ports {
        let Some(budget) = port.wall.and_then(|controls| controls.dissipation_heat) else {
            continue;
        };
        poll(cancelled)?;
        let node = wall_node(nodes, id)?;
        let input = inputs.get(&id).ok_or(Error::Invalid(
            "wall heat requires the original acoustic input",
        ))?;
        let input = FiniteGasAcousticNode::transport_input(input, true).map_err(Error::Invalid)?;
        let convert = |error: fs_duct::DuctError| {
            native_error(
                id,
                match error {
                    fs_duct::DuctError::Cancelled { .. } => AcousticRealizeError::Cancelled,
                    error => AcousticRealizeError::Nonlinear(error.to_string()),
                },
            )
        };
        let channel = node
            .binding
            .duct()
            .prepare_compliant_wall_heat(&node.propagator, budget, cancelled)
            .map_err(convert)?;
        let value = channel
            .integrated_heat_j(&node.state, &input.ports, cancelled)
            .map_err(convert)?;
        if !value.is_finite() || value < 0.0 {
            return Err(Error::Invalid(
                "original wall source heat must be finite and nonnegative",
            ));
        }
        heat.add(value)?;
        by_wall.insert(id, value);
    }
    poll(cancelled)?;
    Ok(PreparedHeat {
        by_wall,
        total_j: heat.value,
    })
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
    heat: &PreparedHeat,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Updates, Error> {
    poll(cancelled)?;
    admit_phase(state, body.mass_kg / body.volume_m3)?;
    let mut steps = BTreeMap::new();
    let mut work = Sum::default();
    for (&id, port) in ports {
        let Some(controls) = port.wall else { continue };
        poll(cancelled)?;
        let original = wall_node(original_nodes, id)?;
        let Some(FiniteGasAcousticNode::Wall(node)) = next_nodes.get_mut(&id) else {
            return Err(Error::Invalid(
                "thermal trial lost its original wall-DG owner",
            ));
        };
        let mut refresh = |source: &LocalMaterialWall| {
            source
                .at_temperature_pinned(
                    port.layers.material,
                    port.layers.temperature_axis,
                    body.initial.temperature_k(),
                    state.temperature_k(),
                    cancelled,
                )
                .map_err(|error| native_error(id, error))
        };
        let (materials, native) = match (&original.binding, &mut node.binding) {
            (FiniteGasDgBinding::Uniform(original), FiniteGasDgBinding::Uniform(wall)) => {
                if wall.wall() != original.wall() {
                    return Err(Error::Invalid(
                        "gas trial changed the thermal wall's original material source",
                    ));
                }
                let replacement = refresh(original.wall())?;
                let materials =
                    SharedFiniteGasAcousticWallMaterials::Uniform(replacement.material().clone());
                let native = wall
                    .rebind_controlled(
                        replacement,
                        &mut node.state,
                        &mut node.propagator,
                        controls.maximum_absolute_parameter_work_j,
                        cancelled,
                    )
                    .map_err(|error| native_error(id, error))?;
                (materials, native)
            }
            (FiniteGasDgBinding::Regional(original), FiniteGasDgBinding::Regional(wall)) => {
                if wall.walls() != original.walls() {
                    return Err(Error::Invalid(
                        "gas trial changed an original regional wall source",
                    ));
                }
                let replacements = original
                    .walls()
                    .iter()
                    .map(&mut refresh)
                    .collect::<Result<Vec<_>, Error>>()?;
                let materials = SharedFiniteGasAcousticWallMaterials::Regional(
                    replacements
                        .iter()
                        .map(|source| source.material().clone())
                        .collect(),
                );
                let native = wall
                    .rebind_controlled(
                        replacements,
                        &mut node.state,
                        &mut node.propagator,
                        controls.maximum_absolute_parameter_work_j,
                        cancelled,
                    )
                    .map_err(|error| native_error(id, error))?;
                (materials, native)
            }
            _ => {
                return Err(Error::Invalid(
                    "thermal trial changed its original wall source layout",
                ));
            }
        };
        work.add(native.parameter_work_j)?;
        steps.insert(
            id,
            SharedFiniteGasAcousticWallStep {
                materials,
                native,
                dissipated_heat_j: heat.by_wall.get(&id).copied().unwrap_or(0.0),
            },
        );
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
