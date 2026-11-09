//! Finite physical thermoelastic exchange on original retained bulk quadrature.
//!
//! Source intervals are prepared once per cell. Native strain reconstruction,
//! B-bar projection and hanging-node maps remain authoritative; no geometry or
//! global physical displacement is rebuilt. Stable source increments are summed
//! in double-double precision before applying the fixed cell stiffness scale.
//! Ghost/support numerical work is deliberately outside this physical account.

use std::ops::ControlFlow;

use fs_material::elastic::thermal::{AffineThermoelastic, ResolvedThermoelastic};
use fs_math::dd::Dd;
use fs_solver::op::LinearOp;

use crate::HexCell;
use crate::elastic3::strain_domain::FrobeniusStrainBall;
use crate::elastic3::{Cell3, CutElasticity3, ElasticityError3, q1};

/// Explicit cell, bulk-work and retained-output allowances.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThermoelasticExchangeLimits3 {
    pub maximum_cells: usize,
    pub maximum_points: usize,
    /// Cell output payload in bytes, checked before allocation. Borrowed fields,
    /// fixed stack scratch, allocator/Vec metadata and the original material
    /// source payload resolved for one cell at a time are excluded.
    pub maximum_workspace_bytes: usize,
}

/// Finite physical account for one original active cell.
///
/// Temperature and displacement endpoints belong to one fixed reference space
/// and the same original law. The current operator supplies the AFTER source.
/// Stiffness scales must be unchanged over the interval; changing design scales
/// requires its own work account. Numerical ghost/support energies, caloric
/// storage and heat transport are excluded. No path-exact or entropy-conserving
/// time integration is asserted: endpoint-average stress requires refinement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoelasticCellExchange3 {
    pub cell: usize,
    pub reference_volume_m3: f64,
    pub temperature_k: [f64; 2],
    pub stiffness_scale: f64,
    pub mechanical_work_j: f64,
    pub free_energy_change_j: f64,
    pub entropy_change_j_k: f64,
    pub internal_energy_change_j: f64,
    /// Positive heat delivered to the caloric owner: W_bulk - Delta U_elastic.
    pub reversible_heat_j: f64,
    pub expansion_work_j: f64,
    pub stiffness_work_j: f64,
    retained: [Dd; 7],
}

impl ThermoelasticCellExchange3 {
    pub const fn retained_mechanical_work_j(&self) -> Dd {
        self.retained[0]
    }

    pub const fn retained_free_energy_change_j(&self) -> Dd {
        self.retained[1]
    }

    pub const fn retained_entropy_change_j_k(&self) -> Dd {
        self.retained[2]
    }

    pub const fn retained_internal_energy_change_j(&self) -> Dd {
        self.retained[3]
    }

    pub const fn retained_reversible_heat_j(&self) -> Dd {
        self.retained[4]
    }

    pub const fn retained_expansion_work_j(&self) -> Dd {
        self.retained[5]
    }

    pub const fn retained_stiffness_work_j(&self) -> Dd {
        self.retained[6]
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

// The actual and constitutive total/elastic component boxes are both required.
// Evaluating only the projected potential would miss a forbidden actual strain
// that B-bar can remove from its cell-mean volumetric constitutive state.
fn admit(
    source: &ResolvedThermoelastic,
    actual: [f64; 6],
    assumed: [f64; 6],
    maximum_kinematic_strain: f64,
) -> Result<(), ElasticityError3> {
    let limit = source.maximum_absolute_mandel_strain();
    let free = source.free_strain_mandel();
    for total in [actual, assumed] {
        if FrobeniusStrainBall::new(total, maximum_kinematic_strain).is_none() {
            return Err(invalid(
                "finite thermal exchange exceeds small-deformation budget; finite-deformation model required",
            ));
        }
        for i in 0..6 {
            let elastic = total[i] - free[i];
            if total[i].abs() > limit || !elastic.is_finite() || elastic.abs() > limit {
                return Err(invalid(
                    "finite thermal exchange leaves original total or elastic strain domain",
                ));
            }
        }
    }
    Ok(())
}

impl CutElasticity3 {
    /// Bounds of the original active reference cell; no geometry is resampled.
    #[must_use]
    pub fn reference_cell_bounds(&self, cell: usize) -> Option<HexCell> {
        self.cells.get(cell).map(|cell| cell.bounds)
    }

    /// Physical cut volume from the original retained integration rule.
    #[must_use]
    pub fn reference_cell_volume_m3(&self, cell: usize) -> Option<f64> {
        self.cells.get(cell).map(|cell| cell.rules.volume())
    }

    /// Whether this operator owns weak embedded displacement supports.
    #[must_use]
    pub fn has_embedded_supports(&self) -> bool {
        self.embedded.is_some()
    }

    /// Integrate stable physical exchange between two native displacement fields.
    ///
    /// Laws and old temperatures use original active-cell order. Each law's
    /// resolved AFTER source must equal this operator's complete retained source,
    /// including axes, domain and claims. The caller retains the same geometry,
    /// fixed stiffness scales and original law at the BEFORE endpoint.
    ///
    /// Strong support coefficients are masked exactly as in the native operator.
    /// Both actual and projected total/elastic strains are admitted at both
    /// endpoints. Only bulk quadrature contributes; numerical support and ghost
    /// work must be accounted for separately by the coupled owner. Cancellation
    /// or any refusal returns no partial cell array and changes no owner.
    #[allow(clippy::too_many_arguments)]
    pub fn thermoelastic_exchange(
        &self,
        laws: &[AffineThermoelastic],
        before_temperatures_k: &[f64],
        before_u: &[f64],
        after_u: &[f64],
        maximum_kinematic_strain: f64,
        limits: ThermoelasticExchangeLimits3,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<ThermoelasticCellExchange3>, ElasticityError3> {
        self.thermoelastic_exchange_mapped(
            laws, before_temperatures_k, before_u, after_u, self.n(),
            maximum_kinematic_strain, limits,
            |cell, field| std::array::from_fn(|i| {
                let node = cell.nodes[i / 3];
                if self.fixed[node] { 0.0 } else { field[3 * node + i % 3] }
            }),
            checkpoint,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(in crate::elastic3) fn thermoelastic_exchange_mapped(
        &self,
        laws: &[AffineThermoelastic],
        before_temperatures_k: &[f64],
        before_u: &[f64],
        after_u: &[f64],
        n: usize,
        maximum_kinematic_strain: f64,
        limits: ThermoelasticExchangeLimits3,
        mut local: impl FnMut(&Cell3, &[f64]) -> [f64; 24],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<ThermoelasticCellExchange3>, ElasticityError3> {
        poll(&mut checkpoint)?;
        if laws.len() != self.cells.len()
            || before_temperatures_k.len() != self.cells.len()
            || before_u.len() != n || after_u.len() != n
            || !maximum_kinematic_strain.is_finite()
            || maximum_kinematic_strain <= 0.0 || maximum_kinematic_strain >= 1.0
        {
            return Err(invalid("invalid finite thermoelastic exchange fields or strain budget"));
        }
        let bytes = self.cells.len()
            .checked_mul(core::mem::size_of::<ThermoelasticCellExchange3>())
            .ok_or_else(|| invalid("finite thermoelastic exchange output count overflow"))?;
        if self.cells.len() > limits.maximum_cells || bytes > limits.maximum_workspace_bytes {
            return Err(invalid("finite thermoelastic exchange cell/output allowance exhausted"));
        }
        for field in [before_u, after_u, before_temperatures_k] {
            for chunk in field.chunks(192) {
                poll(&mut checkpoint)?;
                if chunk.iter().any(|v| !v.is_finite()) {
                    return Err(invalid("nonfinite finite thermoelastic exchange field"));
                }
            }
        }
        let mut points = 0usize;
        for cell in &self.cells {
            poll(&mut checkpoint)?;
            if cell.material.thermal_state().is_none() {
                return Err(invalid("finite thermoelastic exchange requires every original thermal source"));
            }
            points = points.checked_add(cell.rules.bulk().len())
                .ok_or_else(|| invalid("finite thermoelastic exchange point count overflow"))?;
            if points > limits.maximum_points {
                return Err(invalid("finite thermoelastic exchange point allowance exhausted"));
            }
        }
        let mut output = Vec::new();
        output.try_reserve_exact(self.cells.len())
            .map_err(|_| invalid("finite thermoelastic exchange output allocation failed"))?;
        for (id, cell) in self.cells.iter().enumerate() {
            poll(&mut checkpoint)?;
            let source = cell.material.thermal_state()
                .expect("thermal sources checked before output allocation").source();
            let interval = laws[id].prepare_exchange(
                before_temperatures_k[id], source.temperature_k(),
            ).map_err(|_| invalid("finite thermoelastic exchange source interval is invalid"))?;
            if interval.after_source() != source {
                return Err(invalid("finite thermoelastic exchange law differs from retained native source"));
            }
            let before = local(cell, before_u);
            let after = local(cell, after_u);
            let mut retained = [Dd::ZERO; 7];
            for &(position, weight) in cell.rules.bulk() {
                poll(&mut checkpoint)?;
                let (_, gradients) = q1(cell.bounds, position);
                let actual0 = crate::elastic3::stress::mandel_strain(&gradients, &before);
                let actual1 = crate::elastic3::stress::mandel_strain(&gradients, &after);
                let total0 = cell.project_mandel_strain(actual0, &before);
                let total1 = cell.project_mandel_strain(actual1, &after);
                admit(interval.before_source(), actual0, total0, maximum_kinematic_strain)?;
                admit(interval.after_source(), actual1, total1, maximum_kinematic_strain)?;
                let point = interval.evaluate(total0, total1, weight)
                    .map_err(|_| invalid("finite thermoelastic exchange point is unrepresentable"))?;
                let increments = [
                    point.retained_mechanical_work_j(),
                    point.retained_free_energy_change_j(),
                    point.retained_entropy_change_j_k(),
                    point.retained_internal_energy_change_j(),
                    point.retained_reversible_heat_j(),
                    point.retained_expansion_work_j(),
                    point.retained_stiffness_work_j(),
                ];
                for (sum, increment) in retained.iter_mut().zip(increments) {
                    *sum = *sum + increment;
                }
            }
            let scale = self.scales[id];
            for value in &mut retained {
                *value = *value * Dd::from_f64(scale);
            }
            let values = retained.map(Dd::to_f64);
            if values.iter().any(|value| !value.is_finite()) {
                return Err(invalid("scaled finite thermoelastic exchange overflow"));
            }
            output.push(ThermoelasticCellExchange3 {
                cell: id,
                reference_volume_m3: cell.rules.volume(),
                temperature_k: [before_temperatures_k[id], source.temperature_k()],
                stiffness_scale: scale,
                mechanical_work_j: values[0],
                free_energy_change_j: values[1],
                entropy_change_j_k: values[2],
                internal_energy_change_j: values[3],
                reversible_heat_j: values[4],
                expansion_work_j: values[5],
                stiffness_work_j: values[6],
                retained,
            });
        }
        poll(&mut checkpoint)?;
        Ok(output)
    }
}
