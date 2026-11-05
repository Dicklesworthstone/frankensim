//! Spatial, source-backed multipole Debye polarization on the native P1 mesh.
//!
//! Each cell retains `tau_k Pdot_k + P_k = delta_epsilon_k E` and
//! `D = epsilon_infinity E + sum(P_k)`. Eliminating the LOCAL midpoint pole
//! updates gives one positive scalar dielectric operator for the endpoint
//! potential. There are no extra global polarization unknowns. A prepared
//! time step reuses that operator and the immutable original mesh/materials.
//!
//! Free nodal charge is retained on unprescribed nodes; voltage traces exchange
//! charge with their sources. The independently evaluated balance is
//! `Delta U + heat = V_mid . Delta Q_boundary`, where
//! `u = epsilon_infinity |E|^2/2 + sum(|P_k|^2/(2 delta_epsilon_k))`.
//! Polarization heat is distinct from Ohmic leakage. Geometry, coefficients
//! and external material conditions remain fixed. This electroquasistatic
//! lane does not include induction, carrier transport or dielectric breakdown.
//! Midpoint is second order, but large steps can oscillate: resolve the
//! relaxation times and refine time separately from algebraic tolerances.
//!
//! Model reference: COMSOL AC/DC, Dispersion Theory,
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.acdc/acdc_ug_theory.05.36.html>.

use std::sync::Arc;

use fs_feec::{
    assembly::{incidence_to_csr, stiffness},
    weighted::{WeightedAssemblyLimits, WeightedError},
    weighted2::{CellWeight2, weighted_mass_matrix_2d},
};
use fs_material::field_storage::{
    ELECTRIC_PERMITTIVITY_DIMS, dielectric::bulk::ResolvedDebyeMaterial,
};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::SolveReport;
use fs_sparse::Csr;

use super::{
    ElectrostaticControls, ElectrostaticError, accumulate, finite, geometry, poll, positive,
};
use crate::{planar::gradients, scalar_potential};

type Result<T> = core::result::Result<T, ElectrostaticError>;

struct Problem {
    mesh: TriComplex2,
    materials: Vec<ResolvedDebyeMaterial>,
    gradients: Vec<[[f64; 2]; 3]>,
    volumes_m3: Vec<f64>,
    instantaneous_stiffness: Csr,
    limits: WeightedAssemblyLimits,
}

/// Accepted physical field and cell-local polarization history. Clones share
/// immutable geometry/materials; trial steps never modify accepted history.
#[derive(Clone)]
pub struct PlanarDebyeDielectric {
    problem: Arc<Problem>,
    voltage_v: Vec<f64>,
    relative_voltage_v: Vec<f64>,
    electric_field_v_m: Vec<[f64; 2]>,
    displacement_c_m2: Vec<[f64; 2]>,
    polarization_c_m2: Vec<Vec<[f64; 2]>>,
    nodal_charge_c: Vec<f64>,
    cell_energy_j: Vec<f64>,
    stored_energy_j: f64,
    voltage_vertices: Vec<usize>,
    time_s: f64,
}

/// Reusable condensed endpoint operator for ONE fixed material/geometry owner
/// and positive time step. It cannot be used with an unrelated field history.
pub struct PreparedDebyeStep {
    problem: Arc<Problem>,
    dt_s: f64,
    stiffness: Csr,
}
impl PreparedDebyeStep {
    pub fn time_step_s(&self) -> f64 {
        self.dt_s
    }
}

/// Physical exchanges on one accepted-candidate interval. Boundary charge is
/// positive into the dielectric from the voltage source. Local heat can be
/// transferred to a thermal owner without treating stored field energy as heat.
#[derive(Clone, Debug)]
pub struct DebyeFieldStep {
    pub interval_s: [f64; 2],
    pub nodal_boundary_charge_increment_c: Vec<Option<f64>>,
    pub polarization_current_density_a_m2: Vec<[f64; 2]>,
    pub cell_heat_j: Vec<f64>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Sum of local endpoint storage increments, without subtracting two large
    /// total energies. This is the change used in the physical balance.
    pub stored_change_j: f64,
    pub heat_j: f64,
    pub supplied_work_j: f64,
    pub maximum_free_charge_residual_c: f64,
    pub boundary_charge_residual_c: f64,
    pub energy_residual_j: f64,
    pub linear: Option<SolveReport>,
}

/// Transactional result: publish its field only after all physical budgets pass.
pub struct DebyeFieldCandidate {
    field: PlanarDebyeDielectric,
    step: DebyeFieldStep,
}
impl DebyeFieldCandidate {
    pub fn field(&self) -> &PlanarDebyeDielectric {
        &self.field
    }
    pub fn step(&self) -> &DebyeFieldStep {
        &self.step
    }
    pub fn into_parts(self) -> (PlanarDebyeDielectric, DebyeFieldStep) {
        (self.field, self.step)
    }
}

fn invalid(reason: &'static str) -> ElectrostaticError {
    ElectrostaticError::Invalid(reason)
}
fn admit_controls(controls: ElectrostaticControls) -> Result<()> {
    if !positive(controls.linear.tolerance)
        || controls.linear.tolerance >= 1.0
        || !positive(controls.charge_tolerance_c)
        || !positive(controls.energy_tolerance_j)
    {
        return Err(invalid(
            "positive finite Debye field solve and balance budgets required",
        ));
    }
    Ok(())
}

fn assemble(
    mesh: &TriComplex2,
    epsilon: &[f64],
    limits: WeightedAssemblyLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Csr> {
    let weights: Vec<_> = epsilon.iter().copied().map(CellWeight2::Scalar).collect();
    let mass = weighted_mass_matrix_2d(
        mesh,
        1,
        &weights,
        ELECTRIC_PERMITTIVITY_DIMS,
        limits,
        cancelled,
    )
    .map_err(|error| match error {
        WeightedError::Cancelled => ElectrostaticError::Cancelled,
        other => ElectrostaticError::Assembly(other),
    })?;
    poll(cancelled)?;
    let matrix = stiffness(&incidence_to_csr(&mesh.d0()), mass.matrix());
    for row in 0..matrix.nrows() {
        poll(cancelled)?;
        if matrix.row(row).1.iter().any(|value| !value.is_finite()) {
            return Err(invalid("Debye field operator is not representable"));
        }
    }
    Ok(matrix)
}

// Constant-potential nullspace and compensated sums are retained before any
// multiplication, including independently biased disconnected components.
fn apply(matrix: &Csr, voltage: &[f64], cancelled: &mut impl FnMut() -> bool) -> Result<Vec<f64>> {
    let mut result = Vec::with_capacity(matrix.nrows());
    for row in 0..matrix.nrows() {
        poll(cancelled)?;
        let (columns, values) = matrix.row(row);
        let (mut sum, mut correction) = (0.0, 0.0);
        for (&column, &value) in columns.iter().zip(values) {
            accumulate(
                &mut sum,
                &mut correction,
                finite(value * finite(voltage[column] - voltage[row])?)?,
            )?;
        }
        result.push(finite(sum + correction)?);
    }
    Ok(result)
}

impl Problem {
    fn fields(
        &self,
        voltage: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 2]>> {
        let mut fields = Vec::with_capacity(self.materials.len());
        for (face, g) in self.mesh.faces().iter().zip(&self.gradients) {
            poll(cancelled)?;
            let mut e = [0.0; 2];
            for local in 1..3 {
                let difference = finite(voltage[face[local] as usize] - voltage[face[0] as usize])?;
                for axis in 0..2 {
                    e[axis] = finite((-difference).mul_add(g[local][axis], e[axis]))?;
                }
            }
            fields.push(e);
        }
        Ok(fields)
    }

    // Same constant-cell P1 moment as the impressed-polarization weak load.
    // The first local entry closes partition of unity exactly before reduction.
    fn polarization_load(
        &self,
        polarization: &[[f64; 2]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>> {
        let mut rhs = vec![0.0; self.mesh.vertices().len()];
        let mut corrections = vec![0.0; rhs.len()];
        for (cell, face) in self.mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            let g = self.gradients[cell];
            let p = polarization[cell];
            let volume = self.volumes_m3[cell];
            let b1 = finite(volume * finite(g[1][0].mul_add(p[0], g[1][1] * p[1]))?)?;
            let b2 = finite(volume * finite(g[2][0].mul_add(p[0], g[2][1] * p[1]))?)?;
            for (&vertex, load) in face.iter().zip([-finite(b1 + b2)?, b1, b2]) {
                accumulate(
                    &mut rhs[vertex as usize],
                    &mut corrections[vertex as usize],
                    load,
                )?;
            }
        }
        for (value, correction) in rhs.iter_mut().zip(corrections) {
            *value = finite(*value + correction)?;
        }
        Ok(rhs)
    }

    fn charges(
        &self,
        relative_voltage: &[f64],
        polarization: &[Vec<[f64; 2]>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>> {
        let mut total_p = Vec::with_capacity(polarization.len());
        for poles in polarization {
            poll(cancelled)?;
            let mut sum = [0.0; 2];
            let mut corrections = [0.0; 2];
            for p in poles {
                for axis in 0..2 {
                    accumulate(&mut sum[axis], &mut corrections[axis], p[axis])?;
                }
            }
            for axis in 0..2 {
                sum[axis] = finite(sum[axis] + corrections[axis])?;
            }
            total_p.push(sum);
        }
        let load = self.polarization_load(&total_p, cancelled)?;
        apply(&self.instantaneous_stiffness, relative_voltage, cancelled)?
            .iter()
            .zip(load)
            .map(|(&applied, p)| finite(applied - p))
            .collect()
    }
}

impl PlanarDebyeDielectric {
    /// Initialize explicit cell/pole polarization and native nodal potentials.
    /// The resulting weak free charges, including any precharge, are retained;
    /// they are not added again as an impressed volume source. Each component
    /// requires a complete nonzero-length voltage trace. Empty pole vectors
    /// are valid for cells with a purely instantaneous dielectric law.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mesh: TriComplex2,
        materials: Vec<ResolvedDebyeMaterial>,
        voltage_v: Vec<f64>,
        polarization_c_m2: Vec<Vec<[f64; 2]>>,
        voltage_vertices: &[usize],
        time_s: f64,
        limits: WeightedAssemblyLimits,
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        admit_controls(controls)?;
        if !matches!(mesh.metric(), Metric2::Planar { .. })
            || materials.len() != mesh.faces().len()
            || polarization_c_m2.len() != materials.len()
            || materials.len() > limits.max_cells
            || mesh.vertices().len() > limits.max_dofs
            || mesh.edges().len() > limits.max_dofs
            || materials
                .len()
                .checked_mul(9)
                .is_none_or(|n| n > limits.max_triplets)
            || voltage_v.len() != mesh.vertices().len()
            || voltage_v.iter().any(|value| !value.is_finite())
            || voltage_vertices.len() > voltage_v.len()
            || voltage_vertices.iter().any(|v| *v >= voltage_v.len())
            || !time_s.is_finite()
            || time_s < 0.0
        {
            return Err(invalid(
                "bounded planar mesh, complete Debye states, voltages and finite time required",
            ));
        }
        let boundary: Vec<_> = voltage_vertices
            .iter()
            .map(|&v| (v, voltage_v[v]))
            .collect();
        let prescribed = scalar_potential::prescribed_boundary(&mesh, &boundary, cancelled)?;
        let epsilon: Vec<_> = materials
            .iter()
            .map(ResolvedDebyeMaterial::high_frequency_permittivity_f_m)
            .collect();
        let instantaneous_stiffness = assemble(&mesh, &epsilon, limits, cancelled)?;
        let references = scalar_potential::component_references(
            &instantaneous_stiffness,
            &prescribed,
            cancelled,
        )?;
        let relative_voltage_v: Vec<_> = voltage_v
            .iter()
            .zip(&references)
            .map(|(&v, &r)| finite(v - r))
            .collect::<Result<_>>()?;
        let mut cell_gradients = Vec::with_capacity(materials.len());
        let mut volumes_m3 = Vec::with_capacity(materials.len());
        for (cell, face) in mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            cell_gradients
                .push(gradients(face.map(|v| mesh.vertices()[v as usize])).map_err(geometry)?);
            volumes_m3.push(
                mesh.face_measure(cell)
                    .ok_or_else(|| invalid("Debye cell lacks its physical volume"))?,
            );
        }
        let problem = Arc::new(Problem {
            mesh,
            materials,
            gradients: cell_gradients,
            volumes_m3,
            instantaneous_stiffness,
            limits,
        });
        let electric_field_v_m = problem.fields(&relative_voltage_v, cancelled)?;
        let mut displacement_c_m2 = Vec::with_capacity(problem.materials.len());
        let mut cell_energy_j = Vec::with_capacity(problem.materials.len());
        let (mut stored_energy_j, mut correction) = (0.0, 0.0);
        for (cell, law) in problem.materials.iter().enumerate() {
            poll(cancelled)?;
            let e = electric_field_v_m[cell];
            let p = &polarization_c_m2[cell];
            displacement_c_m2.push(
                law.displacement_c_m2(e, p)
                    .map_err(ElectrostaticError::Material)?,
            );
            let energy = finite(
                problem.volumes_m3[cell]
                    * law
                        .stored_energy_density_j_m3(e, p)
                        .map_err(ElectrostaticError::Material)?,
            )?;
            cell_energy_j.push(energy);
            accumulate(&mut stored_energy_j, &mut correction, energy)?;
        }
        let stored_energy_j = finite(stored_energy_j + correction)?;
        let nodal_charge_c = problem.charges(&relative_voltage_v, &polarization_c_m2, cancelled)?;
        let charge_sum = compensated_sum(&nodal_charge_c)?;
        if charge_sum.abs() > controls.charge_tolerance_c {
            return Err(invalid("initial Debye weak charges do not close"));
        }
        let mut voltage_vertices = voltage_vertices.to_vec();
        voltage_vertices.sort_unstable();
        poll(cancelled)?;
        Ok(Self {
            problem,
            voltage_v,
            relative_voltage_v,
            electric_field_v_m,
            displacement_c_m2,
            polarization_c_m2,
            nodal_charge_c,
            cell_energy_j,
            stored_energy_j,
            voltage_vertices,
            time_s,
        })
    }

    pub fn mesh(&self) -> &TriComplex2 {
        &self.problem.mesh
    }
    pub fn materials(&self) -> &[ResolvedDebyeMaterial] {
        &self.problem.materials
    }
    pub fn voltage_v(&self) -> &[f64] {
        &self.voltage_v
    }
    /// Gauge-separated physical voltages used in every field/energy calculation.
    pub fn relative_voltage_v(&self) -> &[f64] {
        &self.relative_voltage_v
    }
    pub fn electric_field_v_m(&self) -> &[[f64; 2]] {
        &self.electric_field_v_m
    }
    pub fn displacement_c_m2(&self) -> &[[f64; 2]] {
        &self.displacement_c_m2
    }
    /// Original cell order, then the original material's pole order.
    pub fn polarization_c_m2(&self) -> &[Vec<[f64; 2]>] {
        &self.polarization_c_m2
    }
    pub fn nodal_charge_c(&self) -> &[f64] {
        &self.nodal_charge_c
    }
    pub fn cell_energy_j(&self) -> &[f64] {
        &self.cell_energy_j
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Condense all local pole updates into a positive endpoint permittivity.
    /// Reuse this object for successive fixed-size steps and different voltage
    /// drives of this same mesh/material owner; no Hodge reassembly is then needed.
    pub fn prepare_step(
        &self,
        dt_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedDebyeStep> {
        poll(cancelled)?;
        if !positive(dt_s) {
            return Err(invalid("positive finite Debye time step required"));
        }
        let mut epsilon = Vec::with_capacity(self.problem.materials.len());
        for law in &self.problem.materials {
            poll(cancelled)?;
            epsilon.push(
                law.midpoint_permittivity_f_m(dt_s)
                    .map_err(ElectrostaticError::Material)?,
            );
        }
        let stiffness = assemble(&self.problem.mesh, &epsilon, self.problem.limits, cancelled)?;
        poll(cancelled)?;
        Ok(PreparedDebyeStep {
            problem: Arc::clone(&self.problem),
            dt_s,
            stiffness,
        })
    }

    /// Advance the original accepted history under linearly ramped prescribed
    /// endpoint voltages. Unprescribed nodes retain their free charge exactly
    /// up to the declared algebraic budget. Failure/cancellation returns no
    /// partial field, polarization update, heat or advanced clock.
    pub fn step_candidate(
        &self,
        prepared: &PreparedDebyeStep,
        endpoint_voltage_v: &[(usize, f64)],
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DebyeFieldCandidate> {
        poll(cancelled)?;
        admit_controls(controls)?;
        if !Arc::ptr_eq(&self.problem, &prepared.problem) {
            return Err(invalid(
                "prepared Debye operator belongs to a different field owner",
            ));
        }
        let dt = prepared.dt_s;
        let time_s = finite(self.time_s + dt)?;
        if time_s <= self.time_s {
            return Err(invalid("Debye time step cannot advance the clock"));
        }
        let prescribed = scalar_potential::prescribed_boundary(
            &self.problem.mesh,
            endpoint_voltage_v,
            cancelled,
        )?;
        if prescribed
            .iter()
            .enumerate()
            .filter_map(|(v, x)| x.map(|_| v))
            .ne(self.voltage_vertices.iter().copied())
        {
            return Err(invalid("Debye voltage-trace topology must remain fixed"));
        }
        let mut history = Vec::with_capacity(self.problem.materials.len());
        for (cell, law) in self.problem.materials.iter().enumerate() {
            poll(cancelled)?;
            history.push(
                law.midpoint_history_c_m2(
                    dt,
                    self.electric_field_v_m[cell],
                    &self.polarization_c_m2[cell],
                )
                .map_err(ElectrostaticError::Material)?,
            );
        }
        let history_load = self.problem.polarization_load(&history, cancelled)?;
        let rhs: Vec<_> = self
            .nodal_charge_c
            .iter()
            .zip(history_load)
            .map(|(&q, p)| finite(q + p))
            .collect::<Result<_>>()?;
        let solved = scalar_potential::solve_relative(
            &prepared.stiffness,
            &rhs,
            &prescribed,
            controls.linear,
            cancelled,
        )?;
        let electric_field_v_m = self.problem.fields(&solved.relative, cancelled)?;
        let cells = self.problem.materials.len();
        let mut polarization_c_m2 = Vec::with_capacity(cells);
        let mut displacement_c_m2 = Vec::with_capacity(cells);
        let mut cell_energy_j = Vec::with_capacity(cells);
        let mut cell_heat_j = Vec::with_capacity(cells);
        let mut cell_storage_change_j = Vec::with_capacity(cells);
        let mut polarization_current_density_a_m2 = Vec::with_capacity(cells);
        let mut polarization_change_c_m2 = Vec::with_capacity(cells);
        for (cell, law) in self.problem.materials.iter().enumerate() {
            poll(cancelled)?;
            let response = law
                .step_midpoint(
                    dt,
                    self.electric_field_v_m[cell],
                    electric_field_v_m[cell],
                    &self.polarization_c_m2[cell],
                )
                .map_err(ElectrostaticError::Material)?;
            let volume = self.problem.volumes_m3[cell];
            cell_energy_j.push(finite(volume * response.stored_after_j_m3)?);
            cell_heat_j.push(finite(volume * response.heat_j_m3)?);
            cell_storage_change_j.push(finite(volume * response.stored_change_j_m3)?);
            let mut change = [0.0; 2];
            let mut corrections = [0.0; 2];
            for (p0, p1) in self.polarization_c_m2[cell]
                .iter()
                .zip(&response.polarization_after_c_m2)
            {
                for axis in 0..2 {
                    accumulate(
                        &mut change[axis],
                        &mut corrections[axis],
                        finite(p1[axis] - p0[axis])?,
                    )?;
                }
            }
            for axis in 0..2 {
                change[axis] = finite(change[axis] + corrections[axis])?;
            }
            let current = [finite(change[0] / dt)?, finite(change[1] / dt)?];
            polarization_change_c_m2.push(change);
            polarization_current_density_a_m2.push(current);
            displacement_c_m2.push(response.displacement_after_c_m2);
            polarization_c_m2.push(response.polarization_after_c_m2);
        }
        let nodal_charge_c =
            self.problem
                .charges(&solved.relative, &polarization_c_m2, cancelled)?;
        // Apply the same native charge functional to endpoint increments.
        // Subtracting two large nodal charges would discard resolved source
        // work even when the local material update retains it.
        let voltage_change: Vec<_> = solved
            .relative
            .iter()
            .zip(&self.relative_voltage_v)
            .map(|(&after, &before)| finite(after - before))
            .collect::<Result<_>>()?;
        let charge_change = apply(
            &self.problem.instantaneous_stiffness,
            &voltage_change,
            cancelled,
        )?;
        let polarization_load_change = self
            .problem
            .polarization_load(&polarization_change_c_m2, cancelled)?;
        let mut nodal_boundary_charge_increment_c = vec![None; nodal_charge_c.len()];
        let mut maximum_free_charge_residual_c = 0.0_f64;
        let (mut boundary_charge_residual_c, mut charge_correction) = (0.0, 0.0);
        let (mut supplied_work_j, mut work_correction) = (0.0, 0.0);
        for vertex in 0..nodal_charge_c.len() {
            poll(cancelled)?;
            let change = finite(charge_change[vertex] - polarization_load_change[vertex])?;
            if prescribed[vertex].is_some() {
                nodal_boundary_charge_increment_c[vertex] = Some(change);
                accumulate(
                    &mut boundary_charge_residual_c,
                    &mut charge_correction,
                    change,
                )?;
                // Each connected component has its own reference; this pairing
                // is invariant under independent endpoint gauge changes.
                let v_mid =
                    finite(0.5 * self.relative_voltage_v[vertex] + 0.5 * solved.relative[vertex])?;
                accumulate(
                    &mut supplied_work_j,
                    &mut work_correction,
                    finite(v_mid * change)?,
                )?;
            } else {
                maximum_free_charge_residual_c = maximum_free_charge_residual_c.max(change.abs());
            }
        }
        let boundary_charge_residual_c = finite(boundary_charge_residual_c + charge_correction)?;
        let supplied_work_j = finite(supplied_work_j + work_correction)?;
        let stored_energy_j = compensated_sum(&cell_energy_j)?;
        let stored_change_j = compensated_sum(&cell_storage_change_j)?;
        let heat_j = compensated_sum(&cell_heat_j)?;
        let energy_residual_j = compensated_sum(&[stored_change_j, heat_j, -supplied_work_j])?;
        if maximum_free_charge_residual_c > controls.charge_tolerance_c
            || boundary_charge_residual_c.abs() > controls.charge_tolerance_c
            || energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(ElectrostaticError::Balance {
                charge_residual_c: boundary_charge_residual_c,
                maximum_free_charge_residual_c,
                energy_residual_j,
            });
        }
        let mut voltage_v = Vec::with_capacity(nodal_charge_c.len());
        for vertex in 0..nodal_charge_c.len() {
            poll(cancelled)?;
            voltage_v.push(match prescribed[vertex] {
                Some(value) => value,
                None => finite(solved.references[vertex] + solved.relative[vertex])?,
            });
        }
        let step = DebyeFieldStep {
            interval_s: [self.time_s, time_s],
            nodal_boundary_charge_increment_c,
            polarization_current_density_a_m2,
            cell_heat_j,
            stored_before_j: self.stored_energy_j,
            stored_after_j: stored_energy_j,
            stored_change_j,
            heat_j,
            supplied_work_j,
            maximum_free_charge_residual_c,
            boundary_charge_residual_c,
            energy_residual_j,
            linear: solved.report,
        };
        poll(cancelled)?;
        Ok(DebyeFieldCandidate {
            field: Self {
                problem: Arc::clone(&self.problem),
                voltage_v,
                relative_voltage_v: solved.relative,
                electric_field_v_m,
                displacement_c_m2,
                polarization_c_m2,
                nodal_charge_c,
                cell_energy_j,
                stored_energy_j,
                voltage_vertices: self.voltage_vertices.clone(),
                time_s,
            },
            step,
        })
    }
}

fn compensated_sum(values: &[f64]) -> Result<f64> {
    let (mut sum, mut correction) = (0.0, 0.0);
    for &value in values {
        accumulate(&mut sum, &mut correction, value)?;
    }
    finite(sum + correction)
}
