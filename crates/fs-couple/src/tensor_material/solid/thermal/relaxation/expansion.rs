//! Source-integrated eigenstrain and stable reciprocal tensor thermodynamics.
use super::*;
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::QueryPoint;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::visco_tensor::TensorThermalExpansion;
use fs_solid::linear3::source::{ThermalExpansionModel, ThermalExpansionSource};
use fs_solid::linear3::{TetElasticMaterial, TetThermalStrainState};

/// Original stress-free source state and fixed expansion axes for one element.
/// Elastic and expansion symmetry/orientation are independent explicit data.
/// The source card is the SAME card already bound by RelaxationThermalCell.
#[derive(Clone)]
pub struct ThermalExpansionCell<'m> {
    pub reference: &'m QueryPoint,
    pub model: ThermalExpansionModel,
    pub selection: MaterialPropertySelection,
    pub mesh_frame: ContentHash,
}

pub(super) struct PreparedExpansion {
    pub derivatives: Vec<TensorThermalExpansion>,
    pub identities: Vec<ContentHash>,
}

fn resolve<E>(
    source: &super::super::super::source::SolidMaterialSource,
    cell: &RelaxationThermalCell<'_>,
    binding: &ThermalExpansionCell<'_>,
    point: &QueryPoint,
) -> Result<ThermalExpansionSource, Error<E>> {
    let material = source
        .material
        .at_query_point(cell.point.material, point)
        .map_err(|e| solid(SolidError::MaterialSource(e)))?;
    let law = material
        .rotated_law(source.material_to_mesh)
        .map_err(|e| solid(SolidError::MaterialSource(e)))?;
    let elastic = TetElasticMaterial::try_new_mandel(
        material.density_kg_m3(),
        *law.equilibrium_mandel(),
        material.state().identity().0,
    )
    .map_err(|e| CircuitEnthalpyError::Source(RelaxationThermalError::Expansion(e.into())))?;
    binding
        .model
        .resolve(
            cell.point.material,
            binding.reference,
            point,
            binding.selection.clone(),
            &elastic,
            binding.mesh_frame,
            law.maximum_strain_norm(),
        )
        .map_err(|e| CircuitEnthalpyError::Source(RelaxationThermalError::Expansion(e)))
}

pub(super) fn prepare<E>(
    history: &TensorSolidHistory<'_, '_>,
    cells: &[RelaxationThermalCell<'_>],
    bindings: &[ThermalExpansionCell<'_>],
    gate: &CancelGate,
) -> Result<PreparedExpansion, Error<E>> {
    if bindings.len() != cells.len() {
        return Err(invalid(
            "complete expansion sources per original element required",
        ));
    }
    let sources = history
        .solid
        .material_sources()
        .ok_or_else(|| invalid("expansion needs original tensor sources"))?;
    let mut prepared = PreparedExpansion {
        derivatives: Vec::with_capacity(cells.len()),
        identities: Vec::with_capacity(cells.len()),
    };
    let mut expected = Vec::with_capacity(cells.len());
    for ((source, cell), binding) in sources.iter().zip(cells).zip(bindings) {
        check(gate)?;
        // The common caller already admitted the original stored
        // tensors over every selected source knot and retained their derivatives.
        // Expansion still requires its own complete affine-alpha source proof.
        if cell.point.temperature_axis != "T" || binding.mesh_frame.0 == [0; 32] {
            return Err(invalid(
                "native integrated expansion requires its T coordinate and a named mesh frame",
            ));
        }
        let query = source.material.state().query();
        let [low, high] = cell.temperature_interval_k;
        let lower = query
            .clone()
            .with("T", low)
            .map_err(|_| invalid("invalid expansion lower source query"))?;
        let upper = query
            .clone()
            .with("T", high)
            .map_err(|_| invalid("invalid expansion upper source query"))?;
        let accepted = resolve(source, cell, binding, &query)?;
        // Preserve the original integrated-strain admission at both ends;
        // derivatives come from the pinned coefficient path, never its secant.
        resolve(source, cell, binding, &lower)?;
        resolve(source, cell, binding, &upper)?;
        let derivative = accepted
            .affine_temperature_derivatives(cell.point.material, cell.temperature_interval_k)
            .map_err(|e| CircuitEnthalpyError::Source(RelaxationThermalError::Expansion(e)))?;
        let mut identity = DomainHasher::new("org.frankensim.tensor-solid.thermal-expansion.v1");
        identity.update(accepted.reference().identity().as_bytes());
        identity.update(binding.mesh_frame.as_bytes());
        match binding.model {
            ThermalExpansionModel::Isotropic => identity.update(&[0]),
            ThermalExpansionModel::OrientedOrthotropic {
                principal_to_mesh,
                principal_frame,
            } => {
                identity.update(&[1]);
                identity.update(principal_frame.as_bytes());
                for value in principal_to_mesh.iter().flatten() {
                    identity.update(&value.to_bits().to_le_bytes());
                }
            }
        }
        prepared.identities.push(identity.finalize());
        prepared.derivatives.push(derivative);
        expected.push(accepted.thermal_strain().free_strain_mandel());
    }
    if history
        .solid
        .thermal_expansion_identity
        .as_ref()
        .is_some_and(|old| *old != prepared.identities)
    {
        return Err(invalid(
            "continuing thermal memory must preserve expansion sources, stress-free reference and frames",
        ));
    }
    for (point, &element) in history
        .solid
        .projection
        .material_points()
        .iter()
        .zip(history.solid.projection.point_elements())
    {
        check(gate)?;
        if point.free_strain_mandel() != &expected[element] {
            return Err(invalid(
                "accepted tensor eigenstrain and integrated thermal source must agree",
            ));
        }
    }
    check(gate)?;
    Ok(prepared)
}

pub(super) fn derivatives<E>(
    model: &PreparedExpansion,
    temperatures_k: &[f64],
    gate: &CancelGate,
) -> Result<Vec<TensorThermalExpansion>, Error<E>> {
    if temperatures_k.len() != model.derivatives.len() {
        return Err(invalid("complete expansion derivative temperatures required"));
    }
    model.derivatives.iter().zip(temperatures_k)
        .map(|(template, &temperature)| {
            check(gate)?;
            template.at_temperature(temperature)
                .map_err(|e| solid(SolidError::Dynamics(
                    crate::tensor_material::TensorMaterialError::Material(e),
                )))
        }).collect()
}

pub(super) fn strains<E>(
    history: &TensorSolidHistory<'_, '_>,
    cells: &[RelaxationThermalCell<'_>],
    bindings: &[ThermalExpansionCell<'_>],
    points: &[QueryPoint],
    gate: &CancelGate,
) -> Result<Vec<TetThermalStrainState>, Error<E>> {
    let sources = history
        .solid
        .material_sources()
        .ok_or_else(|| invalid("expansion needs original tensor sources"))?;
    let mut result = Vec::with_capacity(cells.len());
    for (((source, cell), binding), point) in sources.iter().zip(cells).zip(bindings).zip(points) {
        check(gate)?;
        result.push(*resolve(source, cell, binding, point)?.thermal_strain());
    }
    Ok(result)
}
