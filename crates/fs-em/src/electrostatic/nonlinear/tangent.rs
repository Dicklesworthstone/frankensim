//! Matrix-free analytic field tangent and local differential terminal response.
use super::electrodes::Dofs;
use super::*;
use fs_solver::{CheckedCgError, LinearOp, checked_cg};
use fs_sparse::precond::IdentityPrecond;
use std::cell::{Cell, RefCell};

fn action(element: &Element, epsilon: [[f64; 2]; 2], values: [f64; 3]) -> Result<[f64; 3], Error> {
    let g = element.gradient_of(values);
    let d = [
        -epsilon[0][0].mul_add(g[0], epsilon[0][1] * g[1]),
        -epsilon[1][0].mul_add(g[0], epsilon[1][1] * g[1]),
    ];
    element.reactions(d)
}

struct PreparedTangent<'a> {
    problem: &'a PlanarNonlinearElectrostatics,
    weights: &'a [[[f64; 2]; 2]],
    dofs: &'a Dofs,
    scales: Vec<f64>,
}
impl<'a> PreparedTangent<'a> {
    fn new(
        problem: &'a PlanarNonlinearElectrostatics,
        weights: &'a [[[f64; 2]; 2]],
        dofs: &'a Dofs,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let mut diagonal = vec![0.0; dofs.len()];
        let mut corrections = vec![0.0; dofs.len()];
        for (cell, (element, &epsilon)) in problem.elements.iter().zip(weights).enumerate() {
            poll(cancelled)?;
            CellWeight2::Tensor(epsilon)
                .validate(cell, 1)
                .map_err(weighted)?;
            let indices = element.vertices.map(|v| dofs.index[v]);
            for local in 0..3 {
                let Some(i) = indices[local] else { continue };
                if indices[..local].contains(&Some(i)) {
                    continue;
                }
                // Collapsed terminal coordinates need their whole cell trace,
                // including off-diagonal terms, in the actual Jacobi diagonal.
                let g = element
                    .gradient_of(indices.map(|index| if index == Some(i) { 1.0 } else { 0.0 }));
                let d = [
                    epsilon[0][0].mul_add(g[0], epsilon[0][1] * g[1]),
                    epsilon[1][0].mul_add(g[0], epsilon[1][1] * g[1]),
                ];
                let value = finite(element.volume_m3 * g[0].mul_add(d[0], g[1] * d[1]))?;
                if value < 0.0 || (g != [0.0; 2] && value == 0.0) {
                    return Err(invalid(
                        "electric tangent diagonal is not representably positive",
                    ));
                }
                add(&mut diagonal[i], &mut corrections[i], value)?;
            }
        }
        for (d, c) in diagonal.iter_mut().zip(corrections) {
            *d = finite(*d + c)?.sqrt();
            if !positive(*d) {
                return Err(invalid(
                    "electric tangent has an unanchored or unresolved coordinate",
                ));
            }
        }
        Ok(Self {
            problem,
            weights,
            dofs,
            scales: diagonal,
        })
    }
    fn solve(
        &self,
        rhs: &[f64],
        controls: CheckedCgConfig,
        maximum_applications: usize,
        applications: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, SolveReport), Error> {
        let cancel = RefCell::new(cancelled);
        let interrupted = Cell::new(false);
        let arithmetic = Cell::new(false);
        let budget = Cell::new(false);
        let count = Cell::new(*applications);
        let stop = || {
            if !interrupted.get() && (**cancel.borrow_mut())() {
                interrupted.set(true);
            }
            interrupted.get()
        };
        let operator = Operator {
            prepared: self,
            stop: &stop,
            arithmetic: &arithmetic,
            budget: &budget,
            count: &count,
            maximum: maximum_applications,
        };
        let scaled: Vec<_> = rhs
            .iter()
            .zip(&self.scales)
            .map(|(&v, &s)| finite(v / s))
            .collect::<Result<_, _>>()?;
        let result = checked_cg(&operator, &IdentityPrecond, &scaled, controls, |_| {
            if stop() {
                Err(ElectrostaticError::Cancelled.into())
            } else if budget.get() {
                Err(Error::OperatorBudget)
            } else if arithmetic.get() {
                Err(invalid("electric tangent application is unrepresentable"))
            } else {
                Ok(())
            }
        });
        *applications = count.get();
        let result = result.map_err(|error| match error {
            CheckedCgError::InvalidInput(reason) => invalid(reason),
            CheckedCgError::Interrupted(error) => error,
        })?;
        if !result.report.converged_euclidean() {
            return Err(ElectrostaticError::Nonconvergence(result.report).into());
        }
        let x = result
            .x
            .iter()
            .zip(&self.scales)
            .map(|(&v, &s)| finite(v / s))
            .collect::<Result<Vec<_>, _>>()?;
        Ok((x, result.report))
    }
}

struct Operator<'a> {
    prepared: &'a PreparedTangent<'a>,
    stop: &'a dyn Fn() -> bool,
    arithmetic: &'a Cell<bool>,
    budget: &'a Cell<bool>,
    count: &'a Cell<usize>,
    maximum: usize,
}
impl LinearOp for Operator<'_> {
    fn n(&self) -> usize {
        self.prepared.dofs.len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        if self.count.get() >= self.maximum {
            self.budget.set(true);
            y.fill(f64::NAN);
            return;
        }
        self.count.set(self.count.get() + 1);
        y.fill(0.0);
        let p = self.prepared;
        for (cell, (element, &epsilon)) in p.problem.elements.iter().zip(p.weights).enumerate() {
            if cell % 64 == 0 && (self.stop)() {
                y.fill(f64::NAN);
                return;
            }
            let values = element
                .vertices
                .map(|v| p.dofs.index[v].map_or(0.0, |i| x[i] / p.scales[i]));
            let Ok(reactions) = action(element, epsilon, values) else {
                self.arithmetic.set(true);
                y.fill(f64::NAN);
                return;
            };
            for (v, value) in element.vertices.into_iter().zip(reactions) {
                if let Some(i) = p.dofs.index[v] {
                    y[i] += value / p.scales[i];
                    if !y[i].is_finite() {
                        self.arithmetic.set(true);
                        y.fill(f64::NAN);
                        return;
                    }
                }
            }
        }
    }
}

pub(super) fn direction(
    problem: &PlanarNonlinearElectrostatics,
    weights: &[[[f64; 2]; 2]],
    dofs: &Dofs,
    rhs: &[f64],
    controls: CheckedCgConfig,
    maximum: usize,
    applications: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, SolveReport), Error> {
    PreparedTangent::new(problem, weights, dofs, cancelled)?.solve(
        rhs,
        controls,
        maximum,
        applications,
        cancelled,
    )
}

fn full_action(
    problem: &PlanarNonlinearElectrostatics,
    weights: &[[[f64; 2]; 2]],
    values: &[f64],
    out: &mut [f64],
    correction: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    out.fill(0.0);
    correction.fill(0.0);
    for (element, &epsilon) in problem.elements.iter().zip(weights) {
        poll(cancelled)?;
        let reactions = action(element, epsilon, element.vertices.map(|v| values[v]))?;
        for (v, value) in element.vertices.into_iter().zip(reactions) {
            add(&mut out[v], &mut correction[v], value)?;
        }
    }
    for (value, c) in out.iter_mut().zip(correction) {
        *value = finite(*value + *c)?;
    }
    Ok(())
}

/// Finite extraction work and an absolute capacitance residual [F].
#[derive(Clone, Copy, Debug)]
pub struct DifferentialCapacitanceControls {
    pub maximum_terminals: usize,
    pub maximum_matrix_entries: usize,
    pub max_operator_applications: usize,
    pub coefficient_tolerance_f: f64,
    pub linear: CheckedCgConfig,
}

/// Local dQ/dV on voltage-controlled terminals, with every floating charge
/// held fixed. It remains attached to the actual nonlinear operating point;
/// this is not a global capacitor network or a replacement for its trapped charge.
pub struct DifferentialCapacitance<'a, 'm> {
    field: &'a NonlinearElectricField<'m>,
    terminals: Vec<usize>,
    matrix_f: Vec<f64>,
    maximum_residual_f: f64,
    operator_applications: usize,
}
impl<'a, 'm> DifferentialCapacitance<'a, 'm> {
    pub fn field(&self) -> &'a NonlinearElectricField<'m> {
        self.field
    }
    pub fn terminals(&self) -> &[usize] {
        &self.terminals
    }
    pub fn matrix_f(&self) -> &[f64] {
        &self.matrix_f
    }
    pub fn maximum_residual_f(&self) -> f64 {
        self.maximum_residual_f
    }
    pub fn operator_applications(&self) -> usize {
        self.operator_applications
    }
}

impl<'m> NonlinearElectricField<'m> {
    /// Differentiate the admitted original residual, including equipotential
    /// quotient and zero floating-charge increments. The analytic tangent and
    /// diagonal preparation are reused for all columns; source charge is fixed.
    /// At a genuine material kink a unique differential capacitance does not
    /// exist, so extraction refuses rather than selecting an unstated side.
    pub fn differential_capacitance(
        &self,
        controls: DifferentialCapacitanceControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DifferentialCapacitance<'_, 'm>, Error> {
        poll(cancelled)?;
        let terminals: Vec<_> = self
            .controls
            .iter()
            .enumerate()
            .filter_map(|(t, c)| matches!(c, ElectrodeControl::Voltage(_)).then_some(t))
            .collect();
        let count = terminals.len();
        if count > controls.maximum_terminals
            || count
                .checked_mul(count)
                .is_none_or(|n| n > controls.maximum_matrix_entries)
            || !positive(controls.coefficient_tolerance_f)
            || !positive(controls.linear.tolerance)
            || controls.linear.tolerance >= 1.0
        {
            return Err(invalid(
                "bounded differential terminal matrix and positive tolerances required",
            ));
        }
        for (cell, material) in self.problem.materials.iter().enumerate() {
            poll(cancelled)?;
            let DielectricMaterial::Reversible(law) = material else {
                continue;
            };
            let e = self.evaluation.electric[cell][0].hypot(self.evaluation.electric[cell][1]);
            if e >= law.maximum_field_strength_v_m() {
                return Err(invalid(
                    "differential capacitance needs interior material field support",
                ));
            }
            let knots = &law.source_path().knots;
            if let Ok(i) = knots.binary_search_by(|(v, _)| v.total_cmp(&e)) {
                if i > 0 && i + 1 < knots.len() {
                    let value = |j: usize| knots[j].1.evidence.value.value;
                    let left = (value(i) - value(i - 1)) / (knots[i].0 - knots[i - 1].0);
                    let right = (value(i + 1) - value(i)) / (knots[i + 1].0 - knots[i].0);
                    if left != right {
                        return Err(invalid(
                            "differential capacitance is not unique at a material kink",
                        ));
                    }
                }
            }
        }
        let prepared = PreparedTangent::new(
            self.problem,
            &self.evaluation.tangent,
            &self.dofs,
            cancelled,
        )?;
        let n = self.relative_voltage_v.len();
        let mut values = vec![0.0; n];
        let mut applied = vec![0.0; n];
        let mut corrections = vec![0.0; n];
        let mut rhs = vec![0.0; self.dofs.len()];
        let mut rhs_c = vec![0.0; self.dofs.len()];
        let mut matrix_f = vec![0.0; count * count];
        let mut maximum_residual_f = 0.0_f64;
        let mut operator_applications = 0;
        for (column, &terminal) in terminals.iter().enumerate() {
            poll(cancelled)?;
            values.fill(0.0);
            for &v in &self.electrodes.terminals[terminal] {
                values[v] = 1.0;
            }
            full_action(
                self.problem,
                &self.evaluation.tangent,
                &values,
                &mut applied,
                &mut corrections,
                cancelled,
            )?;
            rhs.fill(0.0);
            rhs_c.fill(0.0);
            for (v, index) in self.dofs.index.iter().enumerate() {
                if let Some(i) = *index {
                    add(&mut rhs[i], &mut rhs_c[i], -applied[v])?;
                }
            }
            for (v, c) in rhs.iter_mut().zip(&rhs_c) {
                *v = finite(*v + c)?;
            }
            if !rhs.is_empty() {
                let (unknown, _) = prepared.solve(
                    &rhs,
                    controls.linear,
                    controls.max_operator_applications,
                    &mut operator_applications,
                    cancelled,
                )?;
                for (v, index) in self.dofs.index.iter().enumerate() {
                    if let Some(i) = *index {
                        values[v] = unknown[i];
                    }
                }
            }
            full_action(
                self.problem,
                &self.evaluation.tangent,
                &values,
                &mut applied,
                &mut corrections,
                cancelled,
            )?;
            for v in 0..n {
                if self.electrodes.node_terminal[v].is_none() {
                    maximum_residual_f = maximum_residual_f.max(applied[v].abs());
                }
            }
            for (t, trace) in self.electrodes.terminals.iter().enumerate() {
                let (mut charge, mut correction) = (0.0, 0.0);
                for &v in trace {
                    add(&mut charge, &mut correction, applied[v])?;
                }
                let charge = finite(charge + correction)?;
                if let Some(row) = terminals.iter().position(|&retained| retained == t) {
                    matrix_f[row * count + column] = charge;
                } else {
                    maximum_residual_f = maximum_residual_f.max(charge.abs());
                }
            }
        }
        for row in 0..count {
            poll(cancelled)?;
            let (mut gauge, mut correction) = (0.0, 0.0);
            for col in 0..count {
                add(&mut gauge, &mut correction, matrix_f[row * count + col])?;
                maximum_residual_f = maximum_residual_f
                    .max((matrix_f[row * count + col] - matrix_f[col * count + row]).abs());
            }
            maximum_residual_f = maximum_residual_f.max(finite(gauge + correction)?.abs());
            if matrix_f[row * count + row] < -controls.coefficient_tolerance_f {
                return Err(invalid(
                    "differential capacitance lost passive diagonal storage",
                ));
            }
        }
        if maximum_residual_f > controls.coefficient_tolerance_f {
            return Err(invalid(
                "differential terminal charge, reciprocity or gauge residual failed",
            ));
        }
        poll(cancelled)?;
        Ok(DifferentialCapacitance {
            field: self,
            terminals,
            matrix_f,
            maximum_residual_f,
            operator_applications,
        })
    }
}
