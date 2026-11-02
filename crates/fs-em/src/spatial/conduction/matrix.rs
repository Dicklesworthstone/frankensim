//! Native multi-electrode Dirichlet-to-current response on the same tetrahedra.
//! Only potentials are superposed. Actual E/J and bulk/contact Joule work are
//! recovered together, retaining every cross term between simultaneous drives.

use super::terminals::{positive, terminal_current, voltage_boundary};
use super::*;
use fs_la::factor::{Cholesky, cholesky};

/// Additional bounds for the reusable terminal response, beyond native field
/// tolerances and the already admitted volume/contact assembly.
#[derive(Clone, Copy, Debug)]
pub struct SpatialConductanceControls {
    pub maximum_terminals: usize,
    /// Exactly N-1 native solves; the first terminal fixes the voltage gauge.
    pub maximum_field_solves: usize,
    /// Conservative scalar storage: (N+1)*vertices + 4*N*N, including retained
    /// relative basis potentials, terminal matrices and factorization buffers.
    pub maximum_response_values: usize,
    /// Conservative dense products: (N-1)^3 + (N-1)^2*(vertices+2).
    /// The native Cholesky owner is polled immediately before and after.
    pub maximum_dense_work: usize,
    /// Positive numerical probe amplitude; never an additional heat source.
    pub reference_voltage_v: f64,
    /// Absolute [S] allowance for energy/reaction coefficients and gauge closure.
    pub coefficient_tolerance_s: f64,
}

/// Frozen linear source and geometry, with a reusable grounded terminal factor.
/// Terminal ordering is caller-owned; terminal zero is the voltage reference.
/// Reprepare after any material, contact or geometry change.
pub struct SpatialConductanceMatrix<'a, C = ResolvedConductorState> {
    problem: &'a SpatialConduction<C>,
    faces: Vec<Vec<usize>>,
    vertices: Vec<Vec<usize>>,
    basis: Vec<Vec<f64>>,
    matrix_s: Vec<f64>,
    factor: Cholesky,
    scales: Vec<f64>,
    maximum_coefficient_residual_s: f64,
}

fn invalid(reason: &'static str) -> SpatialConductionError {
    SpatialConductionError::Invalid(reason)
}

fn add(sum: &mut f64, correction: &mut f64, value: f64) -> Result<()> {
    let adjusted = finite(value - *correction)?;
    let next = finite(*sum + adjusted)?;
    *correction = finite((next - *sum) - adjusted)?;
    *sum = next;
    Ok(())
}

impl<C: SpatialConductorMaterial> SpatialConduction<C> {
    /// Extract the full reciprocal terminal conductance matrix. Each electrode
    /// is a disjoint nonempty union of complete exterior faces. Other faces
    /// insulate; original material contacts remain in the native operator.
    ///
    /// All declared electrodes must connect through bulk, finite contacts or a
    /// shared equipotential electrode. Untouched disconnected bodies retain
    /// their native zero-current gauges. Exactly N-1 basis solves are performed.
    /// The matrix is the original stiffness energy Gram matrix, independently
    /// checked against every basis terminal reaction. No fitted pair-resistor
    /// graph or assumed uniform current injection is introduced.
    pub fn conductance_matrix(
        &self,
        terminals: &[&[usize]],
        controls: SpatialConductanceControls,
        field_controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialConductanceMatrix<'_, C>> {
        check(cancelled)?;
        let n = terminals.len();
        let m = n.saturating_sub(1);
        let vertices = self.mesh().vertex_count;
        let entries = n.checked_mul(n);
        let values = n
            .checked_add(1)
            .and_then(|count| count.checked_mul(vertices))
            .and_then(|basis| entries?.checked_mul(4)?.checked_add(basis));
        let work = m.checked_mul(m).and_then(|square| {
            square
                .checked_mul(m)?
                .checked_add(square.checked_mul(vertices.checked_add(2)?)?)
        });
        if n < 2
            || n > controls.maximum_terminals
            || m > controls.maximum_field_solves
            || values.is_none_or(|v| v > controls.maximum_response_values)
            || work.is_none_or(|w| w > controls.maximum_dense_work)
            || !positive(controls.reference_voltage_v)
            || !positive(controls.coefficient_tolerance_s)
        {
            return Err(invalid(
                "bounded terminal basis/storage/work and positive reference/tolerance required",
            ));
        }
        let (faces, terminal_vertices) = self.electrode_faces(terminals, cancelled)?;
        let mut basis = Vec::with_capacity(m);
        let mut reactions = vec![0.0; n * m];
        for column in 0..m {
            check(cancelled)?;
            let mut voltages = vec![0.0; n];
            voltages[column + 1] = controls.reference_voltage_v;
            let boundary = voltage_boundary(&faces, &voltages);
            let traces = self.prepare_boundary(&boundary, field_controls, cancelled)?;
            let (relative, report) = solve_scalar_system(
                &self.stiffness,
                &traces.flux_load_wb,
                &traces.prescribed,
                field_controls.linear,
                cancelled,
            )?;
            let normalized = relative
                .iter()
                .map(|v| finite(v / controls.reference_voltage_v))
                .collect::<Result<Vec<_>>>()?;
            let field = self.recover_field(
                relative,
                traces,
                field_controls.current_tolerance_a,
                field_controls.power_tolerance_w,
                report,
                cancelled,
            )?;
            for (row, vertices) in terminal_vertices.iter().enumerate() {
                reactions[row * m + column] = finite(
                    terminal_current(&field, vertices, cancelled)? / controls.reference_voltage_v,
                )?;
            }
            basis.push(normalized);
        }

        // K contains both the original bulk Hodge and exact triangular contact
        // operator. This Gram form preserves their cross terms and passivity.
        let mut matrix_s = vec![0.0; n * n];
        for column in 0..m {
            let mut applied = vec![0.0; vertices];
            for (v, out) in applied.iter_mut().enumerate() {
                check(cancelled)?;
                let (indices, values) = self.stiffness.row(v);
                let mut correction = 0.0;
                for (&j, &value) in indices.iter().zip(values) {
                    add(out, &mut correction, finite(value * basis[column][j])?)?;
                }
            }
            for row in 0..=column {
                let (mut value, mut correction) = (0.0, 0.0);
                for (&potential, &current) in basis[row].iter().zip(&applied) {
                    check(cancelled)?;
                    add(&mut value, &mut correction, finite(potential * current)?)?;
                }
                matrix_s[(row + 1) * n + column + 1] = value;
                matrix_s[(column + 1) * n + row + 1] = value;
            }
        }
        for row in 1..n {
            let (mut sum, mut correction) = (0.0, 0.0);
            for column in 1..n {
                add(&mut sum, &mut correction, matrix_s[row * n + column])?;
            }
            matrix_s[row * n] = -sum;
            matrix_s[row] = -sum;
        }
        let (mut reference_diagonal, mut correction) = (0.0, 0.0);
        for &value in &matrix_s[1..n] {
            add(&mut reference_diagonal, &mut correction, -value)?;
        }
        matrix_s[0] = reference_diagonal;
        let mut maximum_coefficient_residual_s = 0.0_f64;
        for row in 0..n {
            check(cancelled)?;
            if !positive(matrix_s[row * n + row]) {
                return Err(invalid(
                    "terminal has no representably positive dissipative connection",
                ));
            }
            let (mut gauge, mut correction) = (0.0, 0.0);
            for column in 0..n {
                add(&mut gauge, &mut correction, matrix_s[row * n + column])?;
                if column != 0 {
                    maximum_coefficient_residual_s = maximum_coefficient_residual_s.max(
                        finite(matrix_s[row * n + column] - reactions[row * m + column - 1])?.abs(),
                    );
                }
            }
            maximum_coefficient_residual_s = maximum_coefficient_residual_s.max(gauge.abs());
        }
        if maximum_coefficient_residual_s > controls.coefficient_tolerance_s {
            return Err(invalid(
                "native terminal reactions, energy Gram matrix or gauge do not agree",
            ));
        }
        let scales: Vec<_> = (1..n).map(|i| matrix_s[i * n + i].sqrt()).collect();
        let mut block = vec![0.0; m * m];
        for row in 0..m {
            check(cancelled)?;
            for column in 0..m {
                block[row * m + column] =
                    finite((matrix_s[(row + 1) * n + column + 1] / scales[row]) / scales[column])?;
            }
        }
        check(cancelled)?;
        let factor = cholesky(&block, m)
            .map_err(|_| invalid("grounded terminal conductance is not positive definite"))?;
        check(cancelled)?;
        Ok(SpatialConductanceMatrix {
            problem: self,
            faces,
            vertices: terminal_vertices,
            basis,
            matrix_s,
            factor,
            scales,
            maximum_coefficient_residual_s,
        })
    }
}

impl<C: SpatialConductorMaterial> SpatialConductanceMatrix<'_, C> {
    /// Canonical face sets in the original terminal ordering.
    pub fn terminal_faces(&self) -> &[Vec<usize>] {
        &self.faces
    }
    /// Full row-major [S] matrix I=G V, including the reference terminal.
    pub fn conductance_s(&self) -> &[f64] {
        &self.matrix_s
    }
    pub fn field_solves(&self) -> usize {
        self.basis.len()
    }
    pub fn maximum_coefficient_residual_s(&self) -> f64 {
        self.maximum_coefficient_residual_s
    }

    /// Recover the original field at arbitrary simultaneous terminal voltages,
    /// without another FE solve. Nodal voltages are relative to terminal zero;
    /// a common potential is removed before all field and power arithmetic.
    /// Actual-amplitude current and power gates remain independent of extraction.
    pub fn field_at_voltages(
        &self,
        voltages_v: &[f64],
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialConductionField> {
        check(cancelled)?;
        let n = self.faces.len();
        if voltages_v.len() != n
            || voltages_v.iter().any(|v| !v.is_finite())
            || !positive(current_tolerance_a)
            || !positive(power_tolerance_w)
        {
            return Err(invalid(
                "one finite voltage per terminal and positive actual-field tolerances required",
            ));
        }
        let voltages = voltages_v
            .iter()
            .map(|v| finite(v - voltages_v[0]))
            .collect::<Result<Vec<_>>>()?;
        let mut relative = vec![0.0; self.problem.mesh().vertex_count];
        for (vertex, out) in relative.iter_mut().enumerate() {
            check(cancelled)?;
            let mut correction = 0.0;
            for (basis, &voltage) in self.basis.iter().zip(&voltages[1..]) {
                add(out, &mut correction, finite(basis[vertex] * voltage)?)?;
            }
        }
        let boundary = voltage_boundary(&self.faces, &voltages);
        let traces = prepare_scalar_boundary_traces(
            &self.problem.geometry.mesh,
            &self.problem.geometry.faces,
            &self.problem.components,
            boundary.voltage.iter().map(|t| (t.face, t.voltage_v)),
            std::iter::empty(),
            current_tolerance_a,
            cancelled,
        )?;
        // Preserve exact authored Dirichlet data even when cancellation among
        // basis potentials loses a few ulps. Native free-node residuals still
        // check the resulting actual field on the full original operator.
        for (vertex, prescribed) in traces.prescribed.iter().enumerate() {
            if let Some(value) = prescribed {
                relative[vertex] = *value;
            }
        }
        let field = self.problem.recover_field(
            relative,
            traces,
            current_tolerance_a,
            power_tolerance_w,
            None,
            cancelled,
        )?;
        let mut currents = vec![0.0; n];
        for (row, current) in currents.iter_mut().enumerate() {
            check(cancelled)?;
            let mut correction = 0.0;
            for (column, &voltage) in voltages.iter().enumerate() {
                add(
                    current,
                    &mut correction,
                    finite(self.matrix_s[row * n + column] * voltage)?,
                )?;
            }
        }
        self.admit_response(
            &field,
            &voltages,
            &currents,
            current_tolerance_a,
            power_tolerance_w,
            cancelled,
        )?;
        Ok(field)
    }

    /// Prescribe signed TOTAL currents into all original equipotential terminals.
    /// Their sum must vanish within the declared current tolerance. A zero-net
    /// electrode floats to its solved potential and may carry nonzero local
    /// current between its faces; it is neither grounded nor locally insulated.
    /// The retained equilibrated native factor is reused for every right-hand
    /// side, followed by original-field recovery and independent I.V checks.
    pub fn field_at_currents(
        &self,
        currents_a: &[f64],
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialConductionField> {
        check(cancelled)?;
        if currents_a.len() != self.faces.len()
            || currents_a.iter().any(|i| !i.is_finite())
            || !positive(current_tolerance_a)
            || !positive(power_tolerance_w)
        {
            return Err(invalid(
                "one finite current per terminal and positive actual-field tolerances required",
            ));
        }
        let (mut sum, mut correction) = (0.0, 0.0);
        for &current in currents_a {
            add(&mut sum, &mut correction, current)?;
        }
        if sum.abs() > current_tolerance_a {
            return Err(invalid("prescribed terminal currents must sum to zero"));
        }
        let mut rhs = currents_a[1..]
            .iter()
            .zip(&self.scales)
            .map(|(&current, &scale)| finite(current / scale))
            .collect::<Result<Vec<_>>>()?;
        check(cancelled)?;
        self.factor.solve(&mut rhs);
        check(cancelled)?;
        let mut voltages = vec![0.0; currents_a.len()];
        for ((voltage, value), &scale) in voltages[1..].iter_mut().zip(rhs).zip(&self.scales) {
            *voltage = finite(value / scale)?;
        }
        let field =
            self.field_at_voltages(&voltages, current_tolerance_a, power_tolerance_w, cancelled)?;
        self.admit_response(
            &field,
            &voltages,
            currents_a,
            current_tolerance_a,
            power_tolerance_w,
            cancelled,
        )?;
        if currents_a.iter().any(|&i| i != 0.0) && !positive(field.power_w) {
            return Err(invalid(
                "nonzero terminal current lost representable Joule work",
            ));
        }
        Ok(field)
    }

    fn admit_response(
        &self,
        field: &SpatialConductionField,
        voltages: &[f64],
        currents: &[f64],
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        let mut current_residual_a = 0.0_f64;
        let (mut work, mut correction) = (0.0, 0.0);
        for ((vertices, &current), &voltage) in self.vertices.iter().zip(currents).zip(voltages) {
            check(cancelled)?;
            let residual = finite(terminal_current(field, vertices, cancelled)? - current)?;
            if residual.abs() > current_residual_a.abs() {
                current_residual_a = residual;
            }
            add(&mut work, &mut correction, finite(voltage * current)?)?;
        }
        let power_residual_w = finite(field.power_w - work)?;
        if current_residual_a.abs() > current_tolerance_a
            || power_residual_w.abs() > power_tolerance_w
        {
            return Err(SpatialConductionError::Balance {
                current_residual_a,
                maximum_free_current_residual_a: field.maximum_free_current_residual_a,
                power_residual_w,
                hodge_power_residual_w: field.hodge_power_residual_w,
            });
        }
        check(cancelled)
    }
}
