//! Matrix-free tetrahedral material tangent and local terminal capacitance.
use super::*;
use fs_solver::{CheckedCgError, LinearOp, checked_cg};
use fs_sparse::precond::IdentityPrecond;
use std::cell::{Cell, RefCell};

pub use super::super::tangent::DifferentialCapacitanceControls;

fn action(
    element: Element<'_>,
    epsilon: [[f64; 3]; 3],
    values: [f64; 4],
) -> Result<[f64; 4], Error> {
    let gradient = element.gradient_of(values);
    element.reactions(epsilon.map(|row| -dot(row, gradient)))
}

struct PreparedTangent<'a> {
    problem: &'a SpatialNonlinearElectrostatics,
    weights: &'a [[[f64; 3]; 3]],
    dofs: &'a Dofs,
    scales: Vec<f64>,
}
impl<'a> PreparedTangent<'a> {
    fn new(
        problem: &'a SpatialNonlinearElectrostatics,
        weights: &'a [[[f64; 3]; 3]],
        dofs: &'a Dofs,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let mut diagonal = vec![0.0; dofs.len()];
        let mut correction = vec![0.0; dofs.len()];
        for (cell, &epsilon) in weights.iter().enumerate() {
            poll(cancelled)?;
            let element = problem.element(cell);
            let indices = element.vertices.map(|v| dofs.index[v as usize]);
            for local in 0..4 {
                let Some(i) = indices[local] else { continue };
                if indices[..local].contains(&Some(i)) {
                    continue;
                }
                // A floating metal is one coordinate, including the cell's
                // off-diagonal couplings between all of its trace vertices.
                let gradient = element
                    .gradient_of(indices.map(|index| if index == Some(i) { 1.0 } else { 0.0 }));
                let displacement = epsilon.map(|row| dot(row, gradient));
                let value = finite(element.volume_m3 * dot(gradient, displacement))?;
                if value < 0.0 || (gradient != [0.0; 3] && value == 0.0) {
                    return Err(invalid(
                        "tetrahedral tangent diagonal is not representably positive",
                    ));
                }
                add(&mut diagonal[i], &mut correction[i], value)?;
            }
        }
        for (value, correction) in diagonal.iter_mut().zip(correction) {
            *value = finite(*value + correction)?.sqrt();
            if !positive(*value) {
                return Err(invalid(
                    "tetrahedral tangent has an unresolved or unanchored coordinate",
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
        settings: CheckedCgConfig,
        maximum: usize,
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
            maximum,
            correction: RefCell::new(vec![0.0; self.dofs.len()]),
        };
        let scaled: Vec<_> = rhs
            .iter()
            .zip(&self.scales)
            .map(|(&value, &scale)| finite(value / scale))
            .collect::<Result<_, _>>()?;
        let result = checked_cg(&operator, &IdentityPrecond, &scaled, settings, |_| {
            if stop() {
                Err(ElectrostaticError::Cancelled.into())
            } else if budget.get() {
                Err(Error::OperatorBudget)
            } else if arithmetic.get() {
                Err(invalid(
                    "tetrahedral tangent application is unrepresentable",
                ))
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
        let values = result
            .x
            .iter()
            .zip(&self.scales)
            .map(|(&value, &scale)| finite(value / scale))
            .collect::<Result<Vec<_>, _>>()?;
        Ok((values, result.report))
    }
}

struct Operator<'a> {
    prepared: &'a PreparedTangent<'a>,
    stop: &'a dyn Fn() -> bool,
    arithmetic: &'a Cell<bool>,
    budget: &'a Cell<bool>,
    count: &'a Cell<usize>,
    maximum: usize,
    correction: RefCell<Vec<f64>>,
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
        let mut correction = self.correction.borrow_mut();
        correction.fill(0.0);
        let p = self.prepared;
        for (cell, &epsilon) in p.weights.iter().enumerate() {
            if cell % 64 == 0 && (self.stop)() {
                y.fill(f64::NAN);
                return;
            }
            let element = p.problem.element(cell);
            let values = element
                .vertices
                .map(|v| p.dofs.index[v as usize].map_or(0.0, |i| x[i] / p.scales[i]));
            let Ok(reactions) = action(element, epsilon, values) else {
                self.arithmetic.set(true);
                y.fill(f64::NAN);
                return;
            };
            for (&vertex, value) in element.vertices.iter().zip(reactions) {
                if let Some(i) = p.dofs.index[vertex as usize] {
                    if add(&mut y[i], &mut correction[i], value / p.scales[i]).is_err() {
                        self.arithmetic.set(true);
                        y.fill(f64::NAN);
                        return;
                    }
                }
            }
        }
        for i in 0..y.len() {
            y[i] += correction[i];
            if !y[i].is_finite() {
                self.arithmetic.set(true);
                y.fill(f64::NAN);
                return;
            }
        }
    }

    // The scaled tangent of an energy (cell Hessians assembled and scaled by
    // the same diagonal on both sides): symmetric; checked_cg solves it.
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.apply(x, y);
    }
}

pub(super) fn direction(
    problem: &SpatialNonlinearElectrostatics,
    weights: &[[[f64; 3]; 3]],
    dofs: &Dofs,
    rhs: &[f64],
    settings: CheckedCgConfig,
    maximum: usize,
    applications: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, SolveReport), Error> {
    PreparedTangent::new(problem, weights, dofs, cancelled)?.solve(
        rhs,
        settings,
        maximum,
        applications,
        cancelled,
    )
}

fn full_action(
    problem: &SpatialNonlinearElectrostatics,
    weights: &[[[f64; 3]; 3]],
    values: &[f64],
    out: &mut [f64],
    correction: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    out.fill(0.0);
    correction.fill(0.0);
    for (cell, &epsilon) in weights.iter().enumerate() {
        poll(cancelled)?;
        let element = problem.element(cell);
        let reactions = action(
            element,
            epsilon,
            element.vertices.map(|v| values[v as usize]),
        )?;
        for (&vertex, charge) in element.vertices.iter().zip(reactions) {
            let v = vertex as usize;
            add(&mut out[v], &mut correction[v], charge)?;
        }
    }
    for (value, correction) in out.iter_mut().zip(correction) {
        *value = finite(*value + *correction)?;
    }
    Ok(())
}

/// Local dQ/dV on voltage-controlled terminals while all floating charges and
/// volume sources stay fixed. This response belongs to the accepted 3D bias;
/// integrating a constant tangent would not reproduce its nonlinear storage.
pub struct SpatialDifferentialCapacitance<'a, 'm> {
    field: &'a SpatialNonlinearElectricField<'m>,
    terminals: Vec<usize>,
    matrix_f: Vec<f64>,
    maximum_residual_f: f64,
    operator_applications: usize,
}
impl<'a, 'm> SpatialDifferentialCapacitance<'a, 'm> {
    pub fn field(&self) -> &'a SpatialNonlinearElectricField<'m> {
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

impl<'m> SpatialNonlinearElectricField<'m> {
    /// Reuse one prepared analytic tangent for all terminal columns. Source
    /// endpoints and genuine constitutive kinks refuse an unspecified two-sided
    /// derivative; the field itself still retains its admitted one-sided law.
    pub fn differential_capacitance(
        &self,
        settings: DifferentialCapacitanceControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialDifferentialCapacitance<'_, 'm>, Error> {
        poll(cancelled)?;
        let terminals: Vec<_> = self
            .controls
            .iter()
            .enumerate()
            .filter_map(|(terminal, control)| {
                matches!(control, ElectrodeControl::Voltage(_)).then_some(terminal)
            })
            .collect();
        let count = terminals.len();
        if count > settings.maximum_terminals
            || count
                .checked_mul(count)
                .is_none_or(|n| n > settings.maximum_matrix_entries)
            || !positive(settings.coefficient_tolerance_f)
            || !positive(settings.linear.tolerance)
            || settings.linear.tolerance >= 1.0
        {
            return Err(invalid(
                "bounded differential terminal matrix and positive tolerances required",
            ));
        }
        for (cell, material) in self.problem.materials.iter().enumerate() {
            poll(cancelled)?;
            let SpatialDielectricMaterial::Reversible(law) = material else {
                continue;
            };
            let e = self.evaluation.electric[cell];
            let magnitude = e[0].hypot(e[1]).hypot(e[2]);
            if magnitude >= law.maximum_field_strength_v_m() {
                return Err(invalid(
                    "differential capacitance needs interior material field support",
                ));
            }
            let knots = &law.source_path().knots;
            if let Ok(i) = knots.binary_search_by(|(value, _)| value.total_cmp(&magnitude)) {
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
        let mut correction = vec![0.0; n];
        let mut rhs = vec![0.0; self.dofs.len()];
        let mut rhs_c = vec![0.0; self.dofs.len()];
        let mut matrix_f = vec![0.0; count * count];
        let mut row_for_terminal = vec![None; self.controls.len()];
        for (row, &terminal) in terminals.iter().enumerate() {
            row_for_terminal[terminal] = Some(row);
        }
        let mut component_charge = vec![0.0; self.controls.len()];
        let mut component_correction = vec![0.0; self.controls.len()];
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
                &mut correction,
                cancelled,
            )?;
            rhs.fill(0.0);
            rhs_c.fill(0.0);
            for (v, index) in self.dofs.index.iter().enumerate() {
                if let Some(i) = *index {
                    add(&mut rhs[i], &mut rhs_c[i], -applied[v])?;
                }
            }
            for (value, correction) in rhs.iter_mut().zip(&rhs_c) {
                *value = finite(*value + *correction)?;
            }
            if !rhs.is_empty() {
                let (unknown, _) = prepared.solve(
                    &rhs,
                    settings.linear,
                    settings.max_operator_applications,
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
                &mut correction,
                cancelled,
            )?;
            for v in 0..n {
                if self.electrodes.node_terminal[v].is_none() {
                    maximum_residual_f = maximum_residual_f.max(applied[v].abs());
                }
            }
            component_charge.fill(0.0);
            component_correction.fill(0.0);
            for (terminal, trace) in self.electrodes.terminals.iter().enumerate() {
                let (mut charge, mut cq) = (0.0, 0.0);
                for &v in trace {
                    add(&mut charge, &mut cq, applied[v])?;
                }
                let charge = finite(charge + cq)?;
                if let Some(row) = row_for_terminal[terminal] {
                    matrix_f[row * count + column] = charge;
                } else {
                    maximum_residual_f = maximum_residual_f.max(charge.abs());
                }
                let reference = self.references.terminal_by_vertex[trace[0]];
                add(
                    &mut component_charge[reference],
                    &mut component_correction[reference],
                    charge,
                )?;
            }
            for (&charge, &correction) in component_charge.iter().zip(&component_correction) {
                maximum_residual_f = maximum_residual_f.max(finite(charge + correction)?.abs());
            }
        }
        for row in 0..count {
            poll(cancelled)?;
            let (mut gauge, mut correction) = (0.0, 0.0);
            for column in 0..count {
                add(&mut gauge, &mut correction, matrix_f[row * count + column])?;
                maximum_residual_f = maximum_residual_f
                    .max((matrix_f[row * count + column] - matrix_f[column * count + row]).abs());
            }
            maximum_residual_f = maximum_residual_f.max(finite(gauge + correction)?.abs());
            if matrix_f[row * count + row] < -settings.coefficient_tolerance_f {
                return Err(invalid(
                    "differential capacitance lost passive diagonal storage",
                ));
            }
        }
        if maximum_residual_f > settings.coefficient_tolerance_f {
            return Err(invalid(
                "differential terminal charge, reciprocity or gauge residual failed",
            ));
        }
        poll(cancelled)?;
        Ok(SpatialDifferentialCapacitance {
            field: self,
            terminals,
            matrix_f,
            maximum_residual_f,
            operator_applications,
        })
    }
}
