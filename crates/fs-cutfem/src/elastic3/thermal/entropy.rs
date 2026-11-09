//! Physical reciprocal thermal/mechanical coupling on the retained bulk rules.
//!
//! The original free energy supplies S = -Psi_T and beta = sigma_T.
//! Its displacement derivative is S_u = -integral(beta : B), with the native
//! B-bar and hanging-node maps. Numerical ghost and support energies do not
//! enter physical entropy, caloric storage, or reversible heat exchange.

use std::ops::ControlFlow;
use std::sync::Arc;

use fs_material::tensor::from_mandel;

use super::jacobian::{ThermalCoordinateMap3, strain_basis};
use crate::elastic3::eigenstrain::physical_response;
use crate::elastic3::{CutElasticity3, ElasticityError3, q1};

/// Explicit bounds for a physical bulk coupling preparation and its reuse.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ThermoelasticEntropyLimits3 {
    /// Retained bulk points; numerical ghost/support points are not integrated.
    pub max_points: usize,
    /// Peak cached scalar slots, each charged as eight bytes. Fixed stack
    /// scratch, Vec/Arc metadata and allocator overhead are excluded.
    pub max_workspace_components: usize,
    /// Scalar stencil work per action, including original hanging rows and
    /// reversible-heat multiplication. Linear field checks are additional.
    pub max_action_terms: usize,
}

/// Extensive elastic-field potential of one original reference cell.
///
/// Every value includes the current stiffness scale and retained reference
/// volume. These are elastic contributions only. The caloric owner must combine
/// the signed fixed-total-strain heat capacity with its original caloric chart;
/// positivity of this contribution alone is neither required nor sufficient.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermoelasticCellPotential3 {
    pub temperature_k: f64,
    pub free_energy_j: f64,
    pub entropy_j_k: f64,
    pub internal_energy_j: f64,
    pub heat_capacity_j_k: f64,
}

struct CellEntropy3 {
    displacement_derivative: [f64; 24],
    scale_bits: u64,
    potential: ThermoelasticCellPotential3,
}

/// Frozen physical entropy derivative and cell potentials at the SAME native u.
///
/// S_u maps independent displacement increments to integrated cell entropy
/// increments [J/K]; its exact discrete transpose uses the same cached columns.
/// This is the bulk source-potential mixed derivative, not the transpose of the
/// COMPLETE stabilized mechanical residual. Nitsche penalties, ghost terms,
/// prescribed support work, caloric/latent storage and transport are excluded.
///
/// Reference geometry, original mass and cell temperatures are unchanged by
/// preparation and actions. Reprepare after changing the physical displacement.
/// Material/geometry owner identity and exact scales reject stale preparations;
/// no stored displacement reference or separately rebuilt constraint map is used.
pub struct PreparedThermoelasticEntropy3 {
    binding: Arc<()>,
    adaptive: bool,
    n: usize,
    cells: Vec<CellEntropy3>,
    points: usize,
    workspace_components: usize,
    action_terms: usize,
}

impl core::fmt::Debug for PreparedThermoelasticEntropy3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("PreparedThermoelasticEntropy3")
            .field("displacement_dofs", &self.n)
            .field("cells", &self.cells.len())
            .field("points", &self.points)
            .field("workspace_components", &self.workspace_components)
            .field("action_terms", &self.action_terms)
            .finish_non_exhaustive()
    }
}

fn invalid(message: &'static str) -> ElasticityError3 {
    ElasticityError3::Invalid(message)
}

fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ElasticityError3> {
    if checkpoint().is_break() {
        Err(ElasticityError3::Cancelled)
    } else {
        Ok(())
    }
}

fn finite(
    values: &[f64],
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
) -> Result<(), ElasticityError3> {
    for chunk in values.chunks(192) {
        poll(checkpoint)?;
        if chunk.iter().any(|value| !value.is_finite()) {
            return Err(invalid("nonfinite thermoelastic entropy field or action"));
        }
    }
    Ok(())
}

impl PreparedThermoelasticEntropy3 {
    /// Original independent displacement dimension, including clamp rows.
    #[must_use]
    pub const fn n(&self) -> usize {
        self.n
    }

    /// Number of physical cell entropy rows, in original active-cell order.
    #[must_use]
    pub fn cells(&self) -> usize {
        self.cells.len()
    }

    /// Original source potential integrated over this reference cell.
    #[must_use]
    pub fn cell_potential(&self, cell: usize) -> Option<&ThermoelasticCellPotential3> {
        self.cells.get(cell).map(|entry| &entry.potential)
    }

    #[must_use]
    pub const fn points(&self) -> usize {
        self.points
    }

    #[must_use]
    pub const fn workspace_components(&self) -> usize {
        self.workspace_components
    }

    #[must_use]
    pub const fn action_terms(&self) -> usize {
        self.action_terms
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(in crate::elastic3) fn prepare(
        source: &CutElasticity3,
        map: ThermalCoordinateMap3<'_>,
        u: &[f64],
        maximum_kinematic_strain: f64,
        limits: ThermoelasticEntropyLimits3,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Self, ElasticityError3> {
        poll(&mut checkpoint)?;
        if u.len() != map.n()
            || !maximum_kinematic_strain.is_finite()
            || maximum_kinematic_strain <= 0.0
            || maximum_kinematic_strain >= 1.0
        {
            return Err(invalid("invalid thermoelastic entropy displacement or strain budget"));
        }
        finite(u, &mut checkpoint)?;
        // Twenty-four coefficients, one scale bit slot and five potential
        // scalars per cell. No global displacement or point-response scratch.
        let workspace_components = source.cells.len().checked_mul(30)
            .ok_or_else(|| invalid("thermoelastic entropy workspace count overflow"))?;
        if workspace_components > limits.max_workspace_components {
            return Err(invalid("thermoelastic entropy workspace allowance exhausted"));
        }
        let mut points = 0usize;
        let mut action_terms = source.cells.len();
        for cell in &source.cells {
            poll(&mut checkpoint)?;
            if cell.material.thermal_state().is_none() {
                return Err(invalid("entropy coupling requires a native thermal source in every cell"));
            }
            points = points.checked_add(cell.rules.bulk().len())
                .ok_or_else(|| invalid("thermoelastic entropy point count overflow"))?;
            if points > limits.max_points {
                return Err(invalid("thermoelastic entropy point allowance exhausted"));
            }
            for &node in &cell.nodes {
                action_terms = map.terms(node).checked_mul(6)
                    .and_then(|terms| action_terms.checked_add(terms))
                    .ok_or_else(|| invalid("thermoelastic entropy action count overflow"))?;
                if action_terms > limits.max_action_terms {
                    return Err(invalid("thermoelastic entropy action allowance exhausted"));
                }
            }
        }
        if action_terms > limits.max_action_terms {
            return Err(invalid("thermoelastic entropy action allowance exhausted"));
        }

        // All size/work admission precedes allocation and potential evaluation.
        let mut cells = Vec::with_capacity(source.cells.len());
        for (id, cell) in source.cells.iter().enumerate() {
            poll(&mut checkpoint)?;
            let thermal = cell.material.thermal_state().expect("thermal source preflight").source();
            let free_mandel = thermal.free_strain_mandel();
            let free = from_mandel(&free_mandel);
            let free_trace = free_mandel[..3].iter().sum::<f64>();
            let local: [f64; 24] = std::array::from_fn(|i| {
                map.value(u, cell.nodes[i / 3], i % 3)
            });
            finite(&local, &mut checkpoint)?;
            let scale = source.scales[id];
            let mut derivative = [0.0; 24];
            let mut potential = ThermoelasticCellPotential3 {
                temperature_k: thermal.temperature_k(),
                free_energy_j: 0.0,
                entropy_j_k: 0.0,
                internal_energy_j: 0.0,
                heat_capacity_j_k: 0.0,
            };
            for &(position, weight) in cell.rules.bulk() {
                poll(&mut checkpoint)?;
                let physical = physical_response(
                    cell, &local, position, free, free_trace, scale,
                    maximum_kinematic_strain,
                )?;
                let point = thermal.state(physical.constitutive_total_mandel, weight)
                    .map_err(|_| invalid("native entropy potential is outside its admitted domain"))?;
                potential.free_energy_j += scale * point.free_energy_j;
                potential.entropy_j_k += scale * point.entropy_j_k;
                potential.internal_energy_j += scale * point.internal_energy_j;
                potential.heat_capacity_j_k += scale * point.heat_capacity_j_k;
                let (_, gradients) = q1(cell.bounds, position);
                let basis = strain_basis(cell, &gradients);
                // beta is a LOCAL Mandel stress slope, without volume or scale.
                // The source retains its original tensor/native K/G arithmetic.
                for i in 3..24 {
                    let work: f64 = basis[i].iter()
                        .zip(point.stress_temperature_derivative_mandel_pa_k)
                        .map(|(strain, beta)| strain * beta).sum();
                    derivative[i] -= (scale * weight) * work;
                }
            }
            // Differentiate the native relative-displacement reconstruction.
            for component in 0..3 {
                derivative[component] = -(1..8)
                    .map(|node| derivative[3 * node + component]).sum::<f64>();
            }
            finite(&derivative, &mut checkpoint)?;
            finite(&[
                potential.temperature_k, potential.free_energy_j,
                potential.entropy_j_k, potential.internal_energy_j,
                potential.heat_capacity_j_k,
            ], &mut checkpoint)?;
            cells.push(CellEntropy3 {
                displacement_derivative: derivative,
                scale_bits: scale.to_bits(),
                potential,
            });
        }
        poll(&mut checkpoint)?;
        Ok(Self {
            binding: Arc::clone(&source.material_body_binding),
            adaptive: map.rows.is_some(),
            n: map.n(),
            cells,
            points,
            workspace_components,
            action_terms,
        })
    }

    fn validate_for(
        &self,
        source: &CutElasticity3,
        map: ThermalCoordinateMap3<'_>,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        poll(checkpoint)?;
        if !Arc::ptr_eq(&self.binding, &source.material_body_binding)
            || self.adaptive != map.rows.is_some()
            || self.n != map.n()
            || self.cells.len() != source.cells.len()
        {
            return Err(invalid("entropy coupling belongs to another operator or space"));
        }
        for (id, (cell, scale)) in self.cells.iter().zip(&source.scales).enumerate() {
            if id % 64 == 0 {
                poll(checkpoint)?;
            }
            if cell.scale_bits != scale.to_bits() {
                return Err(invalid("entropy coupling stiffness scales changed"));
            }
        }
        poll(checkpoint)
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::elastic3) fn apply_mapped(
        &self,
        source: &CutElasticity3,
        map: ThermalCoordinateMap3<'_>,
        input: &[f64],
        output: &mut [f64],
        transpose: bool,
        reversible_heat: bool,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.validate_for(source, map, &mut checkpoint)?;
        let (rows, columns) = if transpose {
            (self.n, self.cells.len())
        } else {
            (self.cells.len(), self.n)
        };
        if input.len() != columns || output.len() != rows {
            return Err(invalid("thermoelastic entropy action shape mismatch"));
        }
        finite(input, &mut checkpoint)?;
        poll(&mut checkpoint)?;
        for chunk in output.chunks_mut(192) {
            poll(&mut checkpoint)?;
            chunk.fill(0.0);
        }
        for (id, (cell, prepared)) in source.cells.iter().zip(&self.cells).enumerate() {
            poll(&mut checkpoint)?;
            if transpose {
                for i in 0..24 {
                    map.scatter(
                        cell.nodes[i / 3], i % 3,
                        prepared.displacement_derivative[i] * input[id], output,
                    );
                }
            } else {
                for component in 0..3 {
                    let origin = map.value(input, cell.nodes[0], component);
                    for node in 1..8 {
                        let direction = map.value(input, cell.nodes[node], component) - origin;
                        output[id] += prepared.displacement_derivative[3 * node + component] * direction;
                    }
                }
                if reversible_heat {
                    output[id] *= -prepared.potential.temperature_k;
                }
            }
        }
        finite(output, &mut checkpoint)?;
        poll(&mut checkpoint)
    }
}

impl CutElasticity3 {
    /// Prepare physical cell entropy, signed elastic capacity and the S_u block.
    ///
    /// The same native actual/constitutive total and elastic strain domains
    /// govern retained BULK points. This is not a boundary admissibility or
    /// solved-equilibrium certificate. Reuse performs no quadrature, source
    /// evaluation, geometric callbacks or allocation.
    pub fn prepare_thermoelastic_entropy(
        &self,
        u: &[f64],
        maximum_kinematic_strain: f64,
        limits: ThermoelasticEntropyLimits3,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<PreparedThermoelasticEntropy3, ElasticityError3> {
        PreparedThermoelasticEntropy3::prepare(
            self, ThermalCoordinateMap3 { rows: None, fixed: &self.fixed },
            u, maximum_kinematic_strain, limits, checkpoint,
        )
    }

    /// S_u times an algebraic displacement direction, in cell entropy rows.
    ///
    /// Directions have no physical strain-amplitude bound. Owner, exact scales,
    /// dimensions and finite inputs are checked before output writes. Cancellation
    /// or arithmetic overflow during accumulation may leave partial output;
    /// discard it on error. Homogeneous strong-support coefficients are masked.
    pub fn thermoelastic_entropy_displacement_apply_into(
        &self,
        prepared: &PreparedThermoelasticEntropy3,
        direction: &[f64],
        output: &mut [f64],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        prepared.apply_mapped(
            self, ThermalCoordinateMap3 { rows: None, fixed: &self.fixed },
            direction, output, false, false, checkpoint,
        )
    }

    /// S_u^T times a cell entropy cotangent using the identical cached stencil.
    /// The native support mask leaves all strong-support output rows zero.
    /// Validation and partial-output semantics match the forward action.
    pub fn thermoelastic_entropy_displacement_transpose_into(
        &self,
        prepared: &PreparedThermoelasticEntropy3,
        entropy_cotangent: &[f64],
        output: &mut [f64],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        prepared.apply_mapped(
            self, ThermalCoordinateMap3 { rows: None, fixed: &self.fixed },
            entropy_cotangent, output, true, false, checkpoint,
        )
    }

    /// Reversible mechanical heat rate -T S_u v = integral(T beta : strain_rate).
    ///
    /// Output is watts per original cell. This is the mechanical term on the
    /// RIGHT of a fixed-total-strain heat-capacity temperature equation. Positive
    /// means heating; expansion commonly cools. The caloric owner must also add
    /// each cell's signed elastic heat capacity exactly once. This is reversible
    /// exchange, not dissipation, a damping law, or a complete thermal step.
    pub fn thermoelastic_reversible_heat_rate_into(
        &self,
        prepared: &PreparedThermoelasticEntropy3,
        velocity: &[f64],
        output: &mut [f64],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        prepared.apply_mapped(
            self, ThermalCoordinateMap3 { rows: None, fixed: &self.fixed },
            velocity, output, false, true, checkpoint,
        )
    }
}
