//! Passive intrinsic material memory on the original tetrahedral midpoint maps.
//!
//! Branch coefficients are eliminated locally before the one global solve.
//! The accepted Gram/polarization states remain private and source-bound. Field
//! recovery independently integrates their stress and free-charge increments;
//! neither the assembled matrix nor a fitted modal damping is a heat oracle.

use super::super::strain_matrix;
use super::*;
use fs_material::piezoelectric::tensor::relaxation::transient::{
    MemoryResponse, PiezoelectricMemory,
};

/// Intrinsic account integrated over one original reference tetrahedron.
#[derive(Clone, Copy, Debug, Default)]
pub struct CellRelaxationStep {
    /// Nonnegative heat from the source-bound Maxwell branches [J].
    pub mechanical_heat_j: f64,
    /// Nonnegative heat from the source-bound Debye poles [J].
    pub dielectric_heat_j: f64,
    /// Material memory storage before/after the window [J], already included
    /// in the body's cell and total stored energy.
    pub memory_energy_j: [f64; 2],
    /// Native endpoint increment, retained separately from rounded inventories.
    pub memory_energy_change_j: f64,
}

#[derive(Clone)]
pub(super) struct Relaxation {
    materials: Arc<[RelaxingTensorPiezoelectric]>,
    memories: Vec<PiezoelectricMemory>,
    steps: Vec<CellRelaxationStep>,
}

impl Relaxation {
    pub(super) fn new(
        model: &TetPiezoelectricModel<'_>,
        materials: Vec<RelaxingTensorPiezoelectric>,
        maximum_internal_variables: usize,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        if materials.len() != 1 && materials.len() != model.tetrahedra.len() {
            return Err(invalid(
                "intrinsic relaxation requires one uniform law or one per original cell",
            ));
        }
        let mut variables = 0usize;
        for cell in 0..model.tetrahedra.len() {
            poll(cx)?;
            let original = &model.materials[if model.materials.len() == 1 { 0 } else { cell }];
            let material = &materials[if materials.len() == 1 { 0 } else { cell }];
            if material.base().identity() != original.identity()
                || material.base().frame_identity() != model.mesh_frame
            {
                return Err(invalid(
                    "relaxation must retain the original cell tensor source and mesh frame",
                ));
            }
            if material.frequency_band_rad_s() != [0.0, f64::MAX] {
                return Err(invalid(
                    "frequency-limited material claims cannot authorize unrestricted transients",
                ));
            }
            variables = variables
                .checked_add(material.internal_variables())
                .ok_or_else(|| invalid("intrinsic cell memory budget overflows"))?;
            if variables > maximum_internal_variables {
                return Err(invalid(
                    "intrinsic cell memory exceeds the explicit scalar budget",
                ));
            }
        }
        let mut memories = Vec::with_capacity(model.tetrahedra.len());
        for cell in 0..model.tetrahedra.len() {
            poll(cx)?;
            memories.push(materials[if materials.len() == 1 { 0 } else { cell }].zero_memory()?);
        }
        Ok(Self {
            materials: materials.into(),
            memories,
            steps: vec![CellRelaxationStep::default(); model.tetrahedra.len()],
        })
    }

    pub(super) fn steps(&self) -> &[CellRelaxationStep] {
        &self.steps
    }

    /// Original source index and law for one original tetrahedron.
    pub(super) fn material(&self, cell: usize) -> (usize, &RelaxingTensorPiezoelectric) {
        let index = if self.materials.len() == 1 { 0 } else { cell };
        (index, &self.materials[index])
    }
}

/// Add local tangent and history terms to the existing coupled increment.
pub(super) fn stamp(
    body: &TetPiezoelectricTransient<'_, '_>,
    dt_s: f64,
    boundary: &mut BoundaryMap,
    entries: &mut Vec<(usize, usize, f64)>,
    cx: &Cx<'_>,
) -> Result<()> {
    let Some(intrinsic) = &body.intrinsic else {
        return Ok(());
    };
    let model = body.model;
    let nu = 3 * model.nodes_m.len();
    let coefficients = intrinsic
        .materials
        .iter()
        .map(|material| {
            poll(cx)?;
            Ok(material.midpoint_coefficients(dt_s)?)
        })
        .collect::<Result<Vec<_>>>()?;
    for (cell, tet) in model.tetrahedra.iter().enumerate() {
        poll(cx)?;
        let material_index = if intrinsic.materials.len() == 1 {
            0
        } else {
            cell
        };
        let coefficient = &coefficients[material_index];
        let history = intrinsic.materials[material_index].history_correction(
            &intrinsic.memories[cell],
            dt_s,
            body.state.field.cells[cell].response.electric_field_v_m,
        )?;
        let gradients = &model.dielectric.gradients()[cell];
        let b = strain_matrix(gradients);
        let volume = model.dielectric.volumes_m3()[cell];
        let cb: [[f64; 12]; 6] = core::array::from_fn(|row| {
            core::array::from_fn(|column| {
                (0..6)
                    .map(|k| coefficient.stiffness_increment_mandel_pa[row][k] * b[k][column])
                    .sum()
            })
        });
        for row in 0..12 {
            let global_row = 3 * tet[row / 3] + row % 3;
            if let Some(index) = boundary.index[global_row] {
                boundary.rhs[index] += volume
                    * (0..6)
                        .map(|k| b[k][row] * history.mechanical_rhs_stress_mandel_pa[k])
                        .sum::<f64>();
            }
            for column in row..12 {
                let global_column = 3 * tet[column / 3] + column % 3;
                let value = volume * (0..6).map(|k| b[k][row] * cb[k][column]).sum::<f64>();
                push_symmetric(entries, global_row, global_column, value)?;
            }
        }
        let eg: [[f64; 3]; 4] = core::array::from_fn(|node| {
            core::array::from_fn(|axis| {
                (0..3)
                    .map(|k| coefficient.permittivity_increment_f_m[axis][k] * gradients[node][k])
                    .sum()
            })
        });
        for row in 0..4 {
            if let Some(index) = boundary.index[nu + tet[row]] {
                boundary.rhs[index] -= volume
                    * (0..3)
                        .map(|axis| {
                            gradients[row][axis] * history.polarization_history_increment_c_m2[axis]
                        })
                        .sum::<f64>();
            }
            for column in row..4 {
                let value = -volume
                    * (0..3)
                        .map(|axis| gradients[row][axis] * eg[column][axis])
                        .sum::<f64>();
                push_symmetric(entries, nu + tet[row], nu + tet[column], value)?;
            }
        }
    }
    poll(cx)
}

fn push_symmetric(
    entries: &mut Vec<(usize, usize, f64)>,
    row: usize,
    column: usize,
    value: f64,
) -> Result<()> {
    if !value.is_finite() {
        return Err(invalid(
            "intrinsic midpoint coefficient is not representable",
        ));
    }
    if value != 0.0 {
        entries.push((row, column, value));
        if row != column {
            entries.push((column, row, value));
        }
    }
    Ok(())
}

pub(super) struct Recovery {
    pub(super) force_n: Vec<[f64; 3]>,
    pub(super) charge_c: Vec<f64>,
    cells: Vec<MemoryResponse>,
}

impl Recovery {
    fn empty(nodes: usize, cells: usize) -> Self {
        Self {
            force_n: vec![[0.0; 3]; nodes],
            charge_c: vec![0.0; nodes],
            cells: Vec::with_capacity(cells),
        }
    }

    fn push(&mut self, model: &TetPiezoelectricModel<'_>, cell: usize, response: MemoryResponse) {
        let tet = model.tetrahedra[cell];
        let gradients = &model.dielectric.gradients()[cell];
        let b = strain_matrix(gradients);
        let volume = model.dielectric.volumes_m3()[cell];
        for local in 0..12 {
            self.force_n[tet[local / 3]][local % 3] += volume
                * (0..6)
                    .map(|axis| b[axis][local] * response.additional_stress_mandel_pa[axis])
                    .sum::<f64>();
        }
        for local in 0..4 {
            self.charge_c[tet[local]] -= volume
                * (0..3)
                    .map(|axis| gradients[local][axis] * response.polarization_c_m2[axis])
                    .sum::<f64>();
        }
        self.cells.push(response);
    }

    fn validate(&self) -> Result<()> {
        if self
            .force_n
            .iter()
            .flatten()
            .chain(&self.charge_c)
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "intrinsic stress force or polarization charge is not representable",
            ));
        }
        Ok(())
    }

    pub(super) fn enrich(
        &self,
        model: &TetPiezoelectricModel<'_>,
        field: &mut TetPiezoelectricSolution,
        cx: &Cx<'_>,
    ) -> Result<()> {
        for (id, (cell, memory)) in field.cells.iter_mut().zip(&self.cells).enumerate() {
            poll(cx)?;
            for axis in 0..6 {
                cell.response.stress_mandel_pa[axis] += memory.additional_stress_mandel_pa[axis];
            }
            for axis in 0..3 {
                cell.response.displacement_c_m2[axis] += memory.polarization_c_m2[axis];
            }
            cell.response.stored_energy_j_m3 += memory.memory_energy_j_m3;
            let energy = model.dielectric.volumes_m3()[id] * memory.memory_energy_j_m3;
            cell.stored_energy_j += energy;
            field.stored_energy_j += energy;
            if cell
                .response
                .stress_mandel_pa
                .iter()
                .chain(&cell.response.displacement_c_m2)
                .chain([
                    &cell.response.stored_energy_j_m3,
                    &cell.stored_energy_j,
                    &field.stored_energy_j,
                ])
                .any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "complete intrinsic material response is not representable",
                ));
            }
        }
        Ok(())
    }
}

pub(super) struct Update {
    pub(super) next: Relaxation,
    pub(super) midpoint: Recovery,
    pub(super) endpoint: Recovery,
    pub(super) mechanical_heat_j: f64,
    pub(super) dielectric_heat_j: f64,
    pub(super) memory_energy_j: [f64; 2],
    pub(super) memory_energy_change: Dd,
}

pub(super) fn advance(
    body: &TetPiezoelectricTransient<'_, '_>,
    dt_s: f64,
    full_increment: &[f64],
    full_midpoint: &[f64],
    controls: TetPiezoelectricControls,
    cx: &Cx<'_>,
) -> Result<Option<Update>> {
    let Some(intrinsic) = &body.intrinsic else {
        return Ok(None);
    };
    let model = body.model;
    let n = model.nodes_m.len();
    let nu = 3 * n;
    let cells = model.tetrahedra.len();
    let mut midpoint = Recovery::empty(n, cells);
    let mut endpoint = Recovery::empty(n, cells);
    let mut memories = Vec::with_capacity(cells);
    let mut steps = Vec::with_capacity(cells);
    let mut mechanical_heat_j = 0.0;
    let mut dielectric_heat_j = 0.0;
    let mut memory_energy_j = [0.0; 2];
    let mut memory_energy_change = Products::default();
    let mut absolute_memory_residual_j = 0.0;
    for (cell, tet) in model.tetrahedra.iter().enumerate() {
        poll(cx)?;
        let material = &intrinsic.materials[if intrinsic.materials.len() == 1 {
            0
        } else {
            cell
        }];
        let gradients = &model.dielectric.gradients()[cell];
        let b = strain_matrix(gradients);
        let volume = model.dielectric.volumes_m3()[cell];
        let mut delta_strain = [0.0; 6];
        let mut electric_field = [0.0; 3];
        for axis in 0..6 {
            for local in 0..12 {
                delta_strain[axis] += b[axis][local]
                    * (full_increment[3 * tet[local / 3] + local % 3]
                        - full_increment[3 * tet[0] + local % 3]);
            }
        }
        for local in 1..4 {
            let potential = full_midpoint[nu + tet[local]] - full_midpoint[nu + tet[0]];
            for axis in 0..3 {
                electric_field[axis] -= potential * gradients[local][axis];
            }
        }
        let (memory, step) = material
            .midpoint_candidate(
                &intrinsic.memories[cell],
                dt_s,
                delta_strain,
                electric_field,
            )?
            .into_parts();
        midpoint.push(model, cell, step.midpoint);
        endpoint.push(model, cell, step.endpoint);
        let report = CellRelaxationStep {
            mechanical_heat_j: volume * step.mechanical_heat_j_m3,
            dielectric_heat_j: volume * step.dielectric_heat_j_m3,
            memory_energy_j: [
                volume * step.memory_energy_before_j_m3,
                volume * step.memory_energy_after_j_m3,
            ],
            memory_energy_change_j: volume * step.memory_energy_change_j_m3,
        };
        memory_energy_change.add(
            Dd::from_f64(volume), Dd::from_f64(step.memory_energy_change_j_m3),
        );
        mechanical_heat_j += report.mechanical_heat_j;
        dielectric_heat_j += report.dielectric_heat_j;
        for (sum, value) in memory_energy_j.iter_mut().zip(report.memory_energy_j) {
            *sum += value;
        }
        absolute_memory_residual_j += volume * step.energy_residual_j_m3.abs();
        memories.push(memory);
        steps.push(report);
    }
    midpoint.validate()?;
    endpoint.validate()?;
    if [
        mechanical_heat_j,
        dielectric_heat_j,
        memory_energy_j[0],
        memory_energy_j[1],
        absolute_memory_residual_j,
    ]
    .iter()
    .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(invalid(
            "intrinsic material energy or heat is not representable",
        ));
    }
    if absolute_memory_residual_j > controls.energy_tolerance_j
        || (memory_energy_j[0] - body.state.material_memory_energy_j).abs()
            > controls.energy_tolerance_j
    {
        return Err(TetPiezoelectricError::Balance(
            "independent material memory work and heat do not close",
        ));
    }
    poll(cx)?;
    Ok(Some(Update {
        next: Relaxation {
            materials: Arc::clone(&intrinsic.materials),
            memories,
            steps,
        },
        midpoint,
        endpoint,
        mechanical_heat_j,
        dielectric_heat_j,
        memory_energy_j,
        memory_energy_change: memory_energy_change.value(),
    }))
}
