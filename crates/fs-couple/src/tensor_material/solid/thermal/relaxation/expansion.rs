//! Source-integrated eigenstrain and stable reciprocal tensor thermodynamics.
use super::*;
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, PropertyKey, QueryPoint};
use fs_material::state_point::MaterialPropertySelection;
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
    pub alpha_mandel: Vec<[f64; 6]>,
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
    thermal: &DepositionSimulation<'_>,
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
        alpha_mandel: Vec::with_capacity(cells.len()),
        identities: Vec::with_capacity(cells.len()),
    };
    let mut expected = Vec::with_capacity(cells.len());
    for ((source, cell), binding) in sources.iter().zip(cells).zip(bindings) {
        check(gate)?;
        // Reuse the native constant-storage admission over ALL original knots.
        source
            .material
            .admit_thermal_relaxation(
                cell.point.material,
                cell.point.temperature_axis,
                cell.temperature_interval_k,
            )
            .map_err(|e| solid(SolidError::MaterialSource(e)))?;
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
        let cold = resolve(source, cell, binding, &lower)?;
        let warm = resolve(source, cell, binding, &upper)?;
        for property in accepted.reference().properties() {
            check(gate)?;
            if !property.has_only_fixed_external_conditions(Some("T"), true) {
                return Err(invalid("expansion cannot freeze unowned source conditions"));
            }
            let requirement = property.requirement();
            let path = cell
                .point
                .material
                .claims()
                .query_line_typed(
                    &PropertyKey::new(requirement.name(), requirement.dims()),
                    "T",
                    &lower,
                    &upper,
                    ClaimSelection::Pinned(property.answer().receipt.selected),
                )
                .map_err(|_| {
                    invalid(
                        "expansion requires continuous support from its original selected source",
                    )
                })?;
            if path.knots.len() < 2
                || path
                    .knots
                    .iter()
                    .any(|(_, answer)| answer.evidence.value.value != property.value_si())
            {
                return Err(invalid(
                    "tensor reciprocal expansion requires constant alpha over every original source knot",
                ));
            }
        }
        let a = cold.thermal_strain().free_strain_mandel();
        let b = warm.thermal_strain().free_strain_mandel();
        let alpha = core::array::from_fn(|i| (b[i] - a[i]) / (high - low));
        if alpha.iter().any(|v| !v.is_finite()) {
            return Err(invalid("source-integrated expansion derivative overflow"));
        }
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
        prepared.alpha_mandel.push(alpha);
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
    let mut coefficients = vec![0.0; cells.len()];
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
        coefficients[element] = finite(
            coefficients[element]
                + point
                    .expansion_capacity_coefficient_j_k2(prepared.alpha_mandel[element])
                    .map_err(|e| {
                        solid(SolidError::Dynamics(
                            crate::tensor_material::TensorMaterialError::Material(e),
                        ))
                    })?,
        )?;
    }
    // The actual caloric owner is piecewise linear in h(T). Check its positive
    // fixed-strain total capacity on every segment, including the hottest end.
    // C_eff = mass*dh/dT - T*V*alpha:(C_eq + sum C_j):alpha.
    for (cell, coefficient) in cells.iter().zip(coefficients) {
        check(gate)?;
        let [low, high] = cell.temperature_interval_k;
        let mass = thermal.cells().get(&cell.cell)
            .ok_or_else(|| invalid("expansion caloric receiver disappeared"))?.mass_kg;
        for segment in cell.caloric.curve().knots().windows(2) {
            let a = segment[0];
            let b = segment[1];
            if b.temperature_k < low || a.temperature_k > high {
                continue;
            }
            let capacity = mass * (b.specific_enthalpy_j_kg - a.specific_enthalpy_j_kg)
                / (b.temperature_k - a.temperature_k)
                - high.min(b.temperature_k) * coefficient;
            if !capacity.is_finite() || capacity <= 0.0 {
                return Err(invalid(
                    "thermal expansion makes the total fixed-strain heat capacity nonpositive",
                ));
            }
        }
    }
    check(gate)?;
    Ok(prepared)
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
