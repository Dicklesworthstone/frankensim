//! One spatial evolution kernel for triangle and tetrahedron dielectric fields.
//! Geometry adapters own traces and native Hodge assembly; this owner retains
//! constitutive state, gauge-separated charge, physical loss and source work.

use super::{
    DebyeFieldStep, ElectrostaticControls, ElectrostaticError, accumulate, finite, positive,
};
use crate::scalar_potential;
use fs_feec::weighted::WeightedAssemblyLimits;
#[cfg(feature = "resolved-leaky-dielectric")]
use fs_material::field_storage::dielectric::leaky::ResolvedLeakyDielectric;
use fs_material::field_storage::{
    ELECTRIC_PERMITTIVITY_DIMS, dielectric::bulk::ResolvedDebyeMaterial,
};
use fs_qty::Dims;
use fs_sparse::Csr;
use std::sync::Arc;

type Result<T> = core::result::Result<T, ElectrostaticError>;

/// Internal geometry seam: no alternate field solver or numerical material law.
pub(crate) trait DebyeGeometry<const D: usize>: Send + Sync {
    fn vertex_count(&self) -> usize;
    fn cell_count(&self) -> usize;
    fn volumes_m3(&self) -> &[f64];
    fn assemble(
        &self,
        coefficients: &[f64],
        dimensions: Dims,
        limits: WeightedAssemblyLimits,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Csr>;
    fn fields(
        &self,
        relative_voltage_v: &[f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Vec<[f64; D]>>;
    fn polarization_load(
        &self,
        polarization_c_m2: &[[f64; D]],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Vec<f64>>;
}

struct Problem<const D: usize> {
    geometry: Arc<dyn DebyeGeometry<D>>,
    materials: Vec<ResolvedDebyeMaterial>,
    conductivity_s_m: Vec<f64>,
    conductivity_properties: Vec<Option<String>>,
    has_conduction: bool,
    instantaneous_stiffness: Csr,
    limits: WeightedAssemblyLimits,
}

#[derive(Clone)]
pub(crate) struct Field<const D: usize> {
    problem: Arc<Problem<D>>,
    voltage_v: Vec<f64>,
    relative_voltage_v: Vec<f64>,
    electric_field_v_m: Vec<[f64; D]>,
    displacement_c_m2: Vec<[f64; D]>,
    polarization_c_m2: Vec<Vec<[f64; D]>>,
    nodal_charge_c: Vec<f64>,
    cell_energy_j: Vec<f64>,
    stored_energy_j: f64,
    voltage_vertices: Vec<usize>,
    time_s: f64,
}

pub(crate) struct Prepared<const D: usize> {
    problem: Arc<Problem<D>>,
    dt_s: f64,
    stiffness: Csr,
    half_conduction_permittivity_f_m: Vec<f64>,
}
impl<const D: usize> Prepared<D> {
    pub(crate) fn time_step_s(&self) -> f64 {
        self.dt_s
    }
}
fn invalid(reason: &'static str) -> ElectrostaticError {
    ElectrostaticError::Invalid(reason)
}
fn poll(cancelled: &mut dyn FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(ElectrostaticError::Cancelled)
    } else {
        Ok(())
    }
}
fn admit_controls(controls: ElectrostaticControls) -> Result<()> {
    if !positive(controls.linear.tolerance)
        || controls.linear.tolerance >= 1.0
        || !positive(controls.charge_tolerance_c)
        || !positive(controls.energy_tolerance_j)
    {
        return Err(invalid(
            "positive finite dielectric solve and balance budgets required",
        ));
    }
    Ok(())
}

fn apply(matrix: &Csr, voltage: &[f64], cancelled: &mut dyn FnMut() -> bool) -> Result<Vec<f64>> {
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

impl<const D: usize> Problem<D> {
    fn charges(
        &self,
        relative_voltage: &[f64],
        polarization: &[Vec<[f64; D]>],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Vec<f64>> {
        let mut total_p = Vec::with_capacity(polarization.len());
        for poles in polarization {
            poll(cancelled)?;
            let mut sum = [0.0; D];
            let mut corrections = [0.0; D];
            for p in poles {
                for axis in 0..D {
                    accumulate(&mut sum[axis], &mut corrections[axis], p[axis])?;
                }
            }
            for axis in 0..D {
                sum[axis] = finite(sum[axis] + corrections[axis])?;
            }
            total_p.push(sum);
        }
        let load = self.geometry.polarization_load(&total_p, cancelled)?;
        apply(&self.instantaneous_stiffness, relative_voltage, cancelled)?
            .iter()
            .zip(load)
            .map(|(&applied, p)| finite(applied - p))
            .collect()
    }
}

impl<const D: usize> Field<D> {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        geometry: Arc<dyn DebyeGeometry<D>>,
        materials: Vec<ResolvedDebyeMaterial>,
        voltage_v: Vec<f64>,
        polarization_c_m2: Vec<Vec<[f64; D]>>,
        prescribed: Vec<Option<f64>>,
        time_s: f64,
        limits: WeightedAssemblyLimits,
        controls: ElectrostaticControls,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Self> {
        if materials.len() > limits.max_cells {
            return Err(invalid("dielectric material count exceeds cell budget"));
        }
        let count = materials.len();
        Self::new_resolved(
            geometry,
            materials,
            vec![0.0; count],
            vec![None; count],
            voltage_v,
            polarization_c_m2,
            prescribed,
            time_s,
            limits,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    #[cfg(feature = "resolved-leaky-dielectric")]
    pub(crate) fn new_leaky(
        geometry: Arc<dyn DebyeGeometry<D>>,
        materials: Vec<ResolvedLeakyDielectric>,
        voltage_v: Vec<f64>,
        polarization_c_m2: Vec<Vec<[f64; D]>>,
        prescribed: Vec<Option<f64>>,
        time_s: f64,
        limits: WeightedAssemblyLimits,
        controls: ElectrostaticControls,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        if materials.len() > limits.max_cells {
            return Err(invalid("leaky material count exceeds cell budget"));
        }
        let mut dielectric = Vec::with_capacity(materials.len());
        let mut conductivity_s_m = Vec::with_capacity(materials.len());
        let mut conductivity_properties = Vec::with_capacity(materials.len());
        for material in materials {
            poll(cancelled)?;
            let (bulk, name, sigma) = material.into_parts();
            dielectric.push(bulk);
            conductivity_s_m.push(sigma);
            conductivity_properties.push(Some(name));
        }
        Self::new_resolved(
            geometry,
            dielectric,
            conductivity_s_m,
            conductivity_properties,
            voltage_v,
            polarization_c_m2,
            prescribed,
            time_s,
            limits,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn new_resolved(
        geometry: Arc<dyn DebyeGeometry<D>>,
        materials: Vec<ResolvedDebyeMaterial>,
        conductivity_s_m: Vec<f64>,
        conductivity_properties: Vec<Option<String>>,
        voltage_v: Vec<f64>,
        polarization_c_m2: Vec<Vec<[f64; D]>>,
        prescribed: Vec<Option<f64>>,
        time_s: f64,
        limits: WeightedAssemblyLimits,
        controls: ElectrostaticControls,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        admit_controls(controls)?;
        let n = geometry.vertex_count();
        if !(2..=3).contains(&D)
            || materials.len() != geometry.cell_count()
            || geometry.volumes_m3().len() != materials.len()
            || geometry.volumes_m3().iter().any(|v| !positive(*v))
            || polarization_c_m2.len() != materials.len()
            || materials.len() > limits.max_cells
            || n > limits.max_dofs
            || voltage_v.len() != n
            || voltage_v.iter().any(|v| !v.is_finite())
            || prescribed.len() != n
            || prescribed
                .iter()
                .zip(&voltage_v)
                .any(|(p, v)| p.is_some_and(|p| !p.is_finite() || p != *v))
            || !time_s.is_finite()
            || time_s < 0.0
        {
            return Err(invalid(
                "bounded geometry, complete dielectric states and finite initial voltages/time required",
            ));
        }
        let epsilon: Vec<_> = materials
            .iter()
            .map(ResolvedDebyeMaterial::high_frequency_permittivity_f_m)
            .collect();
        let instantaneous_stiffness =
            geometry.assemble(&epsilon, ELECTRIC_PERMITTIVITY_DIMS, limits, cancelled)?;
        let references = scalar_potential::component_references(
            &instantaneous_stiffness,
            &prescribed,
            &mut || cancelled(),
        )?;
        let relative_voltage_v: Vec<_> = voltage_v
            .iter()
            .zip(&references)
            .map(|(&v, &r)| finite(v - r))
            .collect::<Result<_>>()?;
        let has_conduction = conductivity_s_m.iter().any(|sigma| *sigma > 0.0);
        let problem = Arc::new(Problem {
            geometry,
            materials,
            conductivity_s_m,
            conductivity_properties,
            has_conduction,
            instantaneous_stiffness,
            limits,
        });
        let electric_field_v_m = problem.geometry.fields(&relative_voltage_v, cancelled)?;
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
                problem.geometry.volumes_m3()[cell]
                    * law
                        .stored_energy_density_j_m3(e, p)
                        .map_err(ElectrostaticError::Material)?,
            )?;
            cell_energy_j.push(energy);
            accumulate(&mut stored_energy_j, &mut correction, energy)?;
        }
        let stored_energy_j = finite(stored_energy_j + correction)?;
        let nodal_charge_c = problem.charges(&relative_voltage_v, &polarization_c_m2, cancelled)?;
        if compensated_sum(&nodal_charge_c)?.abs() > controls.charge_tolerance_c {
            return Err(invalid("initial dielectric weak charges do not close"));
        }
        let voltage_vertices = prescribed
            .iter()
            .enumerate()
            .filter_map(|(v, p)| p.map(|_| v))
            .collect();
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
    pub(crate) fn cell_volumes_m3(&self) -> &[f64] {
        self.problem.geometry.volumes_m3()
    }
    pub(crate) fn conductivity_s_m(&self) -> &[f64] {
        &self.problem.conductivity_s_m
    }
    pub(crate) fn conductivity_properties(&self) -> &[Option<String>] {
        &self.problem.conductivity_properties
    }
    pub(crate) fn same_problem(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.problem, &other.problem)
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
    pub fn electric_field_v_m(&self) -> &[[f64; D]] {
        &self.electric_field_v_m
    }
    pub fn displacement_c_m2(&self) -> &[[f64; D]] {
        &self.displacement_c_m2
    }
    /// Original cell order, then the original material's pole order.
    pub fn polarization_c_m2(&self) -> &[Vec<[f64; D]>] {
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
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Prepared<D>> {
        poll(cancelled)?;
        if !positive(dt_s) {
            return Err(invalid("positive finite Debye time step required"));
        }
        let mut epsilon = Vec::with_capacity(self.problem.materials.len());
        let mut half_conduction_permittivity_f_m = Vec::with_capacity(self.problem.materials.len());
        for (cell, law) in self.problem.materials.iter().enumerate() {
            poll(cancelled)?;
            let half = half_conductivity(dt_s, self.problem.conductivity_s_m[cell])?;
            epsilon.push(finite(
                law.midpoint_permittivity_f_m(dt_s)
                    .map_err(ElectrostaticError::Material)?
                    + half,
            )?);
            half_conduction_permittivity_f_m.push(half);
        }
        let stiffness = self.problem.geometry.assemble(
            &epsilon,
            ELECTRIC_PERMITTIVITY_DIMS,
            self.problem.limits,
            cancelled,
        )?;
        poll(cancelled)?;
        Ok(Prepared {
            problem: Arc::clone(&self.problem),
            dt_s,
            stiffness,
            half_conduction_permittivity_f_m,
        })
    }

    /// Advance accepted history under linearly ramped voltage traces. Free-node
    /// charge change balances its native integrated conduction outflow; the
    /// insulating limit retains free charge. Failure/cancellation returns no
    /// partial field, polarization update, heat or advanced clock.
    pub fn step_candidate(
        &self,
        prepared: &Prepared<D>,
        prescribed: Vec<Option<f64>>,
        controls: ElectrostaticControls,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(Self, DebyeFieldStep<D>)> {
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
        if prescribed.len() != self.voltage_v.len()
            || prescribed.iter().flatten().any(|value| !value.is_finite())
        {
            return Err(invalid(
                "complete finite endpoint voltage partition required",
            ));
        }
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
            let mut local_history = law
                .midpoint_history_c_m2(
                    dt,
                    self.electric_field_v_m[cell],
                    &self.polarization_c_m2[cell],
                )
                .map_err(ElectrostaticError::Material)?;
            let half = prepared.half_conduction_permittivity_f_m[cell];
            if half != 0.0 {
                for (axis, value) in local_history.iter_mut().enumerate() {
                    *value = finite(half.mul_add(self.electric_field_v_m[cell][axis], *value))?;
                }
            }
            history.push(local_history);
        }
        let history_load = self
            .problem
            .geometry
            .polarization_load(&history, cancelled)?;
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
            &mut || cancelled(),
        )?;
        let electric_field_v_m = self.problem.geometry.fields(&solved.relative, cancelled)?;
        let cells = self.problem.materials.len();
        let mut polarization_c_m2 = Vec::with_capacity(cells);
        let mut displacement_c_m2 = Vec::with_capacity(cells);
        let mut cell_energy_j = Vec::with_capacity(cells);
        let mut cell_heat_j = Vec::with_capacity(cells);
        let mut cell_polarization_heat_j = Vec::with_capacity(cells);
        let mut cell_joule_heat_j = Vec::with_capacity(cells);
        let mut cell_storage_change_j = Vec::with_capacity(cells);
        let mut polarization_current_density_a_m2 = Vec::with_capacity(cells);
        let mut polarization_change_c_m2 = Vec::with_capacity(cells);
        let mut conduction_current_density_a_m2 = Vec::with_capacity(cells);
        let mut integrated_conduction_c_m2 = Vec::with_capacity(cells);
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
            let volume = self.problem.geometry.volumes_m3()[cell];
            cell_energy_j.push(finite(volume * response.stored_after_j_m3)?);
            let polarization_heat = finite(volume * response.heat_j_m3)?;
            let sigma = self.problem.conductivity_s_m[cell];
            let mut conduction_current = [0.0; D];
            let mut integrated_conduction = [0.0; D];
            let mut joule_terms = [0.0; D];
            if sigma > 0.0 {
                for axis in 0..D {
                    let midpoint = f64::midpoint(
                        self.electric_field_v_m[cell][axis],
                        electric_field_v_m[cell][axis],
                    );
                    conduction_current[axis] = finite(sigma * midpoint)?;
                    integrated_conduction[axis] = product3(dt, sigma, midpoint)?;
                    joule_terms[axis] = finite(midpoint * integrated_conduction[axis])?;
                    if midpoint != 0.0
                        && (conduction_current[axis] == 0.0
                            || integrated_conduction[axis] == 0.0
                            || joule_terms[axis] == 0.0)
                    {
                        return Err(invalid("nonzero Ohmic current or Joule heat underflows"));
                    }
                }
            }
            let joule_heat = finite(volume * compensated_sum(&joule_terms)?)?;
            cell_polarization_heat_j.push(polarization_heat);
            cell_joule_heat_j.push(joule_heat);
            cell_heat_j.push(finite(polarization_heat + joule_heat)?);
            conduction_current_density_a_m2.push(conduction_current);
            integrated_conduction_c_m2.push(integrated_conduction);
            cell_storage_change_j.push(finite(volume * response.stored_change_j_m3)?);
            let mut change = [0.0; D];
            let mut corrections = [0.0; D];
            for (p0, p1) in self.polarization_c_m2[cell]
                .iter()
                .zip(&response.polarization_after_c_m2)
            {
                for axis in 0..D {
                    accumulate(
                        &mut change[axis],
                        &mut corrections[axis],
                        finite(p1[axis] - p0[axis])?,
                    )?;
                }
            }
            for axis in 0..D {
                change[axis] = finite(change[axis] + corrections[axis])?;
            }
            let mut current = [0.0; D];
            for axis in 0..D {
                current[axis] = finite(change[axis] / dt)?;
            }
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
            .geometry
            .polarization_load(&polarization_change_c_m2, cancelled)?;
        // -integral grad(N).(dt J_mid) is the weak conduction outflow.
        // It permits zero conductivity without constructing a singular or
        // fictitiously positive conductivity Hodge.
        let nodal_conduction_outflow_c = if self.problem.has_conduction {
            self.problem
                .geometry
                .polarization_load(&integrated_conduction_c_m2, cancelled)?
                .into_iter()
                .map(|value| finite(-value))
                .collect::<Result<Vec<_>>>()?
        } else {
            vec![0.0; nodal_charge_c.len()]
        };
        let mut nodal_boundary_charge_increment_c = vec![None; nodal_charge_c.len()];
        let mut nodal_boundary_supplied_charge_c = vec![None; nodal_charge_c.len()];
        let mut maximum_free_charge_residual_c = 0.0_f64;
        let (mut boundary_charge_residual_c, mut charge_correction) = (0.0, 0.0);
        let (mut supplied_displacement_work_j, mut displacement_work_correction) = (0.0, 0.0);
        let (mut supplied_conduction_work_j, mut conduction_work_correction) = (0.0, 0.0);
        for vertex in 0..nodal_charge_c.len() {
            poll(cancelled)?;
            let change = finite(charge_change[vertex] - polarization_load_change[vertex])?;
            let outflow = nodal_conduction_outflow_c[vertex];
            let supplied_charge = finite(change + outflow)?;
            if prescribed[vertex].is_some() {
                nodal_boundary_charge_increment_c[vertex] = Some(change);
                nodal_boundary_supplied_charge_c[vertex] = Some(supplied_charge);
                accumulate(
                    &mut boundary_charge_residual_c,
                    &mut charge_correction,
                    supplied_charge,
                )?;
                // Each connected component has its own reference; this pairing
                // is invariant under independent endpoint gauge changes.
                let v_mid =
                    finite(0.5 * self.relative_voltage_v[vertex] + 0.5 * solved.relative[vertex])?;
                accumulate(
                    &mut supplied_displacement_work_j,
                    &mut displacement_work_correction,
                    finite(v_mid * change)?,
                )?;
                accumulate(
                    &mut supplied_conduction_work_j,
                    &mut conduction_work_correction,
                    finite(v_mid * outflow)?,
                )?;
            } else {
                maximum_free_charge_residual_c =
                    maximum_free_charge_residual_c.max(supplied_charge.abs());
            }
        }
        let boundary_charge_residual_c = finite(boundary_charge_residual_c + charge_correction)?;
        let supplied_displacement_work_j =
            finite(supplied_displacement_work_j + displacement_work_correction)?;
        let supplied_conduction_work_j =
            finite(supplied_conduction_work_j + conduction_work_correction)?;
        let supplied_work_j =
            compensated_sum(&[supplied_displacement_work_j, supplied_conduction_work_j])?;
        let stored_energy_j = compensated_sum(&cell_energy_j)?;
        let stored_change_j = compensated_sum(&cell_storage_change_j)?;
        let heat_j = compensated_sum(&cell_heat_j)?;
        let polarization_heat_j = compensated_sum(&cell_polarization_heat_j)?;
        let joule_heat_j = compensated_sum(&cell_joule_heat_j)?;
        // Retain small dielectric changes under large DC through-power by
        // reducing the original independent contributions before rounding a
        // total heat/work scalar. No balance defect is converted into heat.
        let energy_residual_j = compensated_sum(&[
            stored_change_j,
            polarization_heat_j,
            joule_heat_j,
            -supplied_displacement_work_j,
            -supplied_conduction_work_j,
        ])?;
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
            nodal_boundary_supplied_charge_c,
            nodal_conduction_outflow_c,
            polarization_current_density_a_m2,
            conduction_current_density_a_m2,
            cell_polarization_heat_j,
            cell_joule_heat_j,
            cell_heat_j,
            stored_before_j: self.stored_energy_j,
            stored_after_j: stored_energy_j,
            stored_change_j,
            polarization_heat_j,
            joule_heat_j,
            heat_j,
            supplied_work_j,
            supplied_displacement_work_j,
            supplied_conduction_work_j,
            maximum_free_charge_residual_c,
            boundary_charge_residual_c,
            energy_residual_j,
            linear: solved.report,
        };
        poll(cancelled)?;
        Ok((
            Self {
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
        ))
    }
}

fn compensated_sum(values: &[f64]) -> Result<f64> {
    let (mut sum, mut correction) = (0.0, 0.0);
    for &value in values {
        accumulate(&mut sum, &mut correction, value)?;
    }
    finite(sum + correction)
}

fn half_conductivity(dt: f64, sigma: f64) -> Result<f64> {
    if sigma == 0.0 {
        return Ok(0.0);
    }
    // Halve the larger factor first: avoid both unnecessary overflow and
    // double rounding of representable subnormal midpoint coefficients.
    let half = if dt >= sigma {
        (0.5 * dt) * sigma
    } else {
        dt * (0.5 * sigma)
    };
    if positive(half) {
        Ok(half)
    } else {
        Err(invalid("midpoint Ohmic coefficient is unrepresentable"))
    }
}

fn product3(a: f64, b: f64, c: f64) -> Result<f64> {
    let mut factors = [a, b, c];
    factors.sort_by(|x, y| x.abs().total_cmp(&y.abs()));
    finite((factors[0] * factors[2]) * factors[1])
}
