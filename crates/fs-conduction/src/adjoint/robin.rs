//! Linear FEM response to Robin references, log(h), and assembled nodal loads.
//!
//! A tangent solves A dT = db - dA T; a pullback solves A^T lambda = dJ/dT.
//! The same production P1 operator and preconditioned CG serve both. Uniform
//! Robin faces use their consistent mass matrix, NOT a lumped approximation.
//! Prescribed temperatures enter dA T but their tangents/adjoints are zero.
//! This owns the solid half of a conjugate derivative, not the air fixed point.
//! Geometry, conductivity and prescribed temperatures are held fixed; k(T),
//! radiation coupling and continuum-error certification are not covered.
//! Explicit matching-P1 contact is supported with its resistance held fixed.

use std::collections::BTreeSet;

use fs_exec::Cx;
use fs_solver::{CgState, CsrOp, norm2};
use fs_sparse::Csr;

use crate::assemble::{DofMap, assemble_operator_scaled_with_interfaces, reduce};
use crate::{ConductionError, ConductionProblem, ConductionSolution, LinearConfig, ScalarField,
    SolveConfig, ThermalBc, ThermalInterfaces};

/// One selected uniform Robin region. Ordering is chosen by the caller.
#[derive(Debug)]
pub struct RobinPort {
    /// Exact boundary-region name.
    pub name: String,
    /// Integrated boundary area, m^2.
    pub area_m2: f64,
    /// Nominal transfer coefficient, W/(m^2 K).
    pub htc_w_m2_k: f64,
    /// Nominal Robin reference, K.
    pub reference_k: f64,
    faces: Vec<([usize; 3], f64)>,
}

/// Tangent inputs in port order, plus a full assembled-load perturbation.
#[derive(Debug, Clone, PartialEq)]
pub struct RobinDirection {
    /// Changes in uniform references, K.
    pub references_k: Vec<f64>,
    /// Changes in ln(h), equivalent to ln(hA) only at fixed area.
    pub log_htc: Vec<f64>,
    /// Changes to assembled nodal loads, W; prescribed-node entries do not
    /// affect temperature. This is not a volumetric-source density vector.
    pub nodal_load_w: Vec<f64>,
}

/// One tangent, with the explicitly recomputed algebraic solve residual.
#[derive(Debug, Clone, PartialEq)]
pub struct RobinDifferential {
    /// Full nodal temperature changes, K; zero on prescribed vertices.
    pub temperature_k: Vec<f64>,
    /// Area-mean wall-temperature changes, K, in port order.
    pub mean_wall_temperatures_k: Vec<f64>,
    /// Outward heat-rate changes, W, in port order.
    pub heat_rates_w: Vec<f64>,
    /// Recomputed Euclidean relative residual of the normalized linear solve.
    pub relative_residual: f64,
    /// Krylov iterations, excluding the primal solve.
    pub iterations: usize,
}

/// Gradient of a linear functional of nodal T, mean wall T and outward heat.
#[derive(Debug, Clone, PartialEq)]
pub struct RobinGradient {
    /// Derivatives with respect to reference temperatures, in port order.
    pub references: Vec<f64>,
    /// Derivatives with respect to ln(h), in port order.
    pub log_htc: Vec<f64>,
    /// Derivatives with respect to assembled nodal loads; zero on fixed nodes.
    pub nodal_load: Vec<f64>,
    /// Recomputed Euclidean relative residual of the normalized adjoint solve.
    pub relative_residual: f64,
    /// Krylov iterations, excluding the primal solve.
    pub iterations: usize,
}

/// An owned operator bound to a newly solved and residual-checked primal.
/// Publicly assembled ConductionSolution values cannot be substituted for it.
pub struct RobinLinearization {
    primal: ConductionSolution,
    matrix: Csr,
    dofs: DofMap,
    ports: Vec<RobinPort>,
    linear: LinearConfig,
}

impl RobinLinearization {
    /// Solve a temperature-independent conduction problem and bind selected
    /// uniform Robin regions. Non-selected boundaries remain in the operator.
    /// Refuses unknown/duplicate/nonuniform ports, nonlinear materials, invalid
    /// solve budgets, a failed true-residual gate, and cancellation.
    pub fn new(
        cx: &Cx<'_>, problem: ConductionProblem<'_>, config: SolveConfig, regions: &[&str],
    ) -> Result<Self, ConductionError> {
        Self::new_inner(cx, problem, None, config, regions)
    }

    /// Bind the same Robin derivatives through explicitly declared matching-P1
    /// contacts. The primal, tangent and adjoint all contain the SAME contact
    /// operator; temperatures on the two traces remain distinct unknowns.
    /// Contact resistance is fixed, not a Robin log(h) control. This does not
    /// differentiate contact geometry or resistance and does not infer perfect
    /// contact. Interface ownership, geometry and missing-pair checks remain
    /// the production ThermalInterfaces/solve_with_interfaces admission rules.
    ///
    /// # Errors
    /// All errors from `new`, plus the production contact-binding refusals.
    pub fn new_with_interfaces(
        cx: &Cx<'_>, problem: ConductionProblem<'_>, interfaces: &ThermalInterfaces,
        config: SolveConfig, regions: &[&str],
    ) -> Result<Self, ConductionError> {
        Self::new_inner(cx, problem, Some(interfaces), config, regions)
    }

    fn new_inner(
        cx: &Cx<'_>, problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
        config: SolveConfig, regions: &[&str],
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        if !(config.linear.tolerance.is_finite() && config.linear.tolerance > 0.0
            && config.linear.tolerance < 1.0 && config.linear.max_iterations > 0)
        {
            return Err(invalid("positive Krylov budget and a relative tolerance in (0, 1) required"));
        }
        if let Some(materials) = problem.element_materials {
            materials.validate_for(problem.mesh)?;
            for element in 0..problem.mesh.element_count() {
                poll(cx, element)?;
                if materials.model_for(element)?.is_temperature_dependent() {
                    return Err(invalid("Robin sensitivities require temperature-independent materials"));
                }
            }
        } else if problem.material.is_temperature_dependent() {
            return Err(invalid("Robin sensitivities require temperature-independent materials"));
        }
        let mut seen = BTreeSet::new();
        let mut ports = Vec::with_capacity(regions.len());
        for &name in regions {
            poll(cx, ports.len())?;
            if !seen.insert(name) { return Err(invalid("duplicate Robin sensitivity region")); }
            let region = problem.boundary.region_names().iter().position(|n| n == name)
                .ok_or_else(|| invalid("unknown Robin sensitivity region"))?;
            let ThermalBc::Robin { htc: ScalarField::Uniform(h), t_ref: ScalarField::Uniform(r) }
                = &problem.boundary.conditions()[region]
            else { return Err(invalid("selected sensitivity regions require uniform Robin h and reference")); };
            let mut port = RobinPort { name: name.to_string(), area_m2: 0.0,
                htc_w_m2_k: *h, reference_k: *r, faces: Vec::new() };
            for (slot, face) in problem.mesh.boundary().iter().enumerate() {
                poll(cx, slot)?;
                if problem.boundary.region_for(slot) == Some(region) {
                    add(&mut port.area_m2, face.area)?;
                    port.faces.push((face.vertices.map(|v| v as usize), face.area));
                }
            }
            if port.area_m2 <= 0.0 { return Err(invalid("a Robin sensitivity region has no area")); }
            ports.push(port);
        }
        let linear = config.linear;
        let primal = match interfaces {
            Some(interfaces) => crate::solve::solve_with_interfaces(cx, problem, interfaces, config)?,
            None => crate::solve::solve(cx, problem, config)?,
        };
        let system = assemble_operator_scaled_with_interfaces(cx, problem.mesh, problem.boundary,
            problem.material, problem.source, &primal.temperature, None, interfaces, problem.element_materials)?;
        let dofs = DofMap::new(problem.boundary, problem.mesh.vertex_count())?;
        let (matrix, rhs) = reduce(&system, &dofs);
        let relative = true_residual(&matrix, &dofs.gather(&primal.temperature), &rhs)?;
        if relative >= linear.tolerance {
            return Err(failed(0, relative, linear));
        }
        poll(cx, 0)?;
        Ok(Self { primal, matrix, dofs, ports, linear })
    }

    /// Freshly executed primal; its references and coefficients define this map.
    #[must_use]
    pub const fn primal(&self) -> &ConductionSolution { &self.primal }

    /// Selected regions in the exact order used by every input/output vector.
    #[must_use]
    pub fn ports(&self) -> &[RobinPort] { &self.ports }

    /// Correctly sized zero perturbation.
    #[must_use]
    pub fn zero_direction(&self) -> RobinDirection {
        RobinDirection { references_k: vec![0.0; self.ports.len()], log_htc: vec![0.0; self.ports.len()],
            nodal_load_w: vec![0.0; self.primal.temperature.len()] }
    }

    /// Area-average a full nodal field on the selected Robin traces.
    /// This accepts differences or absolute temperatures; it never adds a lift.
    pub fn wall_means(&self, cx: &Cx<'_>, field: &[f64]) -> Result<Vec<f64>, ConductionError> {
        vector(cx, field, self.primal.temperature.len())?;
        self.ports.iter().map(|port| {
            let mut mean = 0.0;
            for (vertices, area) in &port.faces {
                poll(cx, 0)?;
                for &vertex in vertices { add(&mut mean, (area / port.area_m2 / 3.0) * field[vertex])?; }
            }
            Ok(mean)
        }).collect()
    }

    /// Apply db-dA*T and solve once; no perturbed primal or solver-iteration AD.
    pub fn apply(&self, cx: &Cx<'_>, direction: &RobinDirection) -> Result<RobinDifferential, ConductionError> {
        vector(cx, &direction.references_k, self.ports.len())?;
        vector(cx, &direction.log_htc, self.ports.len())?;
        vector(cx, &direction.nodal_load_w, self.primal.temperature.len())?;
        let mut rhs = direction.nodal_load_w.clone();
        for (index, port) in self.ports.iter().enumerate() {
            for (vertices, area) in &port.faces {
                poll(cx, index)?;
                for (a, &vertex) in vertices.iter().enumerate() {
                    let mut load = port.htc_w_m2_k * (area / 3.0) * direction.references_k[index];
                    for (b, &other) in vertices.iter().enumerate() {
                        let mass = port.htc_w_m2_k * (area / 12.0) * if a == b { 2.0 } else { 1.0 };
                        add(&mut load, mass * (port.reference_k - self.primal.temperature[other]) * direction.log_htc[index])?;
                    }
                    add(&mut rhs[vertex], load)?;
                }
            }
        }
        let (temperature_k, relative_residual, iterations) = self.solve_rhs(cx, &rhs)?;
        let mean_wall_temperatures_k = self.wall_means(cx, &temperature_k)?;
        let means = self.wall_means(cx, &self.primal.temperature)?;
        let heat_rates_w = self.ports.iter().enumerate().map(|(i, port)| {
            checked(port.htc_w_m2_k * port.area_m2 * (mean_wall_temperatures_k[i]
                - direction.references_k[i] + (means[i] - port.reference_k) * direction.log_htc[i]))
        }).collect::<Result<Vec<_>, _>>()?;
        poll(cx, iterations)?;
        Ok(RobinDifferential { temperature_k, mean_wall_temperatures_k, heat_rates_w, relative_residual, iterations })
    }

    /// One transposed solve for weights on full nodal temperature, selected
    /// area-mean wall temperatures, and selected outward heat rates.
    pub fn pullback(
        &self, cx: &Cx<'_>, nodal_weights: &[f64], wall_weights: &[f64], heat_weights: &[f64],
    ) -> Result<RobinGradient, ConductionError> {
        vector(cx, nodal_weights, self.primal.temperature.len())?;
        vector(cx, wall_weights, self.ports.len())?;
        vector(cx, heat_weights, self.ports.len())?;
        let mut rhs = nodal_weights.to_vec();
        for (i, port) in self.ports.iter().enumerate() {
            for (vertices, area) in &port.faces {
                poll(cx, i)?;
                for &vertex in vertices {
                    add(&mut rhs[vertex], (area / 3.0) * (wall_weights[i] / port.area_m2
                        + heat_weights[i] * port.htc_w_m2_k))?;
                }
            }
        }
        let (lambda, relative_residual, iterations) = self.solve_rhs(cx, &rhs)?;
        let means = self.wall_means(cx, &self.primal.temperature)?;
        let mut references = vec![0.0; self.ports.len()];
        let mut log_htc = vec![0.0; self.ports.len()];
        for (i, port) in self.ports.iter().enumerate() {
            let conductance = checked(port.htc_w_m2_k * port.area_m2)?;
            references[i] = checked(-heat_weights[i] * conductance)?;
            log_htc[i] = checked(heat_weights[i] * conductance * (means[i] - port.reference_k))?;
            for (vertices, area) in &port.faces {
                poll(cx, i)?;
                for (a, &vertex) in vertices.iter().enumerate() {
                    add(&mut references[i], lambda[vertex] * port.htc_w_m2_k * (area / 3.0))?;
                    for (b, &other) in vertices.iter().enumerate() {
                        let mass = port.htc_w_m2_k * (area / 12.0) * if a == b { 2.0 } else { 1.0 };
                        add(&mut log_htc[i], lambda[vertex] * mass * (port.reference_k - self.primal.temperature[other]))?;
                    }
                }
            }
        }
        poll(cx, iterations)?;
        Ok(RobinGradient { references, log_htc, nodal_load: lambda, relative_residual, iterations })
    }

    fn solve_rhs(&self, cx: &Cx<'_>, rhs: &[f64]) -> Result<(Vec<f64>, f64, usize), ConductionError> {
        poll(cx, 0)?;
        let rhs = self.dofs.gather(rhs);
        let scale = rhs.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        if scale == 0.0 { return Ok((vec![0.0; self.primal.temperature.len()], 0.0, 0)); }
        let normalized: Vec<f64> = rhs.iter().map(|v| v / scale).collect();
        let op = CsrOp::symmetric(self.matrix.clone());
        let pre = crate::solve::spd_preconditioner(&self.matrix);
        let mut state = CgState::new(&op, &pre, &normalized);
        while state.rel_residual() >= self.linear.tolerance && state.iters < self.linear.max_iterations {
            poll(cx, state.iters)?;
            let before = state.iters;
            let batch = (self.linear.max_iterations - before).min(16);
            state.run(&op, &pre, self.linear.tolerance, batch);
            if state.iters == before { break; }
        }
        poll(cx, state.iters)?;
        let residual = true_residual(&self.matrix, &state.x, &normalized)?;
        if residual >= self.linear.tolerance { return Err(failed(state.iters, residual, self.linear)); }
        let mut full = vec![0.0; self.primal.temperature.len()];
        for (i, &vertex) in self.dofs.free().iter().enumerate() {
            poll(cx, i)?;
            full[vertex] = checked(state.x[i] * scale)?;
        }
        Ok((full, residual, state.iters))
    }
}

fn true_residual(matrix: &Csr, x: &[f64], rhs: &[f64]) -> Result<f64, ConductionError> {
    let mut ax = vec![0.0; rhs.len()];
    matrix.spmv(x, &mut ax);
    let scale = rhs.iter().chain(&ax).map(|v| v.abs()).fold(0.0_f64, f64::max);
    for &v in rhs.iter().chain(&ax) { checked(v)?; }
    if scale == 0.0 { return Ok(0.0); }
    let residual: Vec<f64> = rhs.iter().zip(&ax).map(|(b, a)| b / scale - a / scale).collect();
    let normalized: Vec<f64> = rhs.iter().map(|b| b / scale).collect();
    checked(norm2(&residual) / norm2(&normalized).max(f64::MIN_POSITIVE))
}
fn vector(cx: &Cx<'_>, values: &[f64], expected: usize) -> Result<(), ConductionError> {
    poll(cx, 0)?;
    if values.len() != expected { return Err(ConductionError::FieldLength { field: "Robin sensitivity vector", expected, found: values.len() }); }
    for (i, &value) in values.iter().enumerate() { if i % 512 == 0 { poll(cx, i)?; } checked(value)?; }
    Ok(())
}
fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled { stage: "robin-sensitivity", at })
}
fn checked(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() { Ok(value) }
    else { Err(ConductionError::NonFinite { field: "Robin sensitivity", bits: value.to_bits() }) }
}
fn add(sum: &mut f64, value: f64) -> Result<(), ConductionError> { *sum = checked(*sum + value)?; Ok(()) }
fn invalid(what: &str) -> ConductionError { ConductionError::Config { parameter: "Robin sensitivity", what: what.to_string() } }
fn failed(iterations: usize, residual: f64, config: LinearConfig) -> ConductionError {
    ConductionError::LinearSolveFailed { iteration: 0, krylov_iterations: iterations,
        true_relative_residual: residual, tolerance: config.tolerance }
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod contact_tests;
