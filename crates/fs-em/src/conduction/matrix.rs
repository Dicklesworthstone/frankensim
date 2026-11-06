//! Multi-electrode resistance on the original planar conductivity operator.
//!
//! Exactly N-1 native basis solves prepare a reciprocal energy Gram matrix.
//! Subsequent voltage, total-current and mixed drives superpose POTENTIALS and
//! recover the actual native cell/contact heat, including all cross terms.

use super::*;
use fs_la::factor::{Cholesky, FactorRunError, cholesky_cancellable};

/// Bounded preparation of a frozen planar terminal response.
#[derive(Clone, Copy, Debug)]
pub struct PlanarConductanceControls {
    pub maximum_terminals: usize,
    /// Exactly N-1 native solves; terminal zero fixes the reporting gauge.
    pub maximum_field_solves: usize,
    /// Scalar response storage bound: (N+4)*vertices + 4*N*N. Native field
    /// assembly/recovery retains its separately admitted geometry bounds.
    pub maximum_response_values: usize,
    /// Conservative dense work: (N-1)^3 + (N-1)^2*(vertices+2).
    pub maximum_dense_work: usize,
    /// Positive physical probe amplitude, inside the original source domain.
    pub reference_voltage_v: f64,
    /// Absolute [S] reaction/energy coefficient and gauge allowance.
    pub coefficient_tolerance_s: f64,
}

/// Frozen source/geometry with an equilibrated grounded terminal factor.
/// Reprepare after a material, contact, geometry or nonlinear state change.
pub struct PlanarConductanceMatrix<'m, C = ResolvedConductorState> {
    problem: &'m PlanarConduction<C>,
    terminals: Vec<Vec<usize>>,
    basis: Vec<Vec<f64>>,
    matrix_s: Vec<f64>,
    factor: Cholesky,
    scales: Vec<f64>,
    maximum_coefficient_residual_s: f64,
}

/// Bounds for preparing a reusable mixed voltage/current terminal partition.
#[derive(Clone, Copy, Debug)]
pub struct PlanarMixedConductanceControls {
    pub maximum_current_terminals: usize,
    /// Scalar storage for current block/factor and vectors: 2*F*F + 3*N.
    pub maximum_response_values: usize,
    /// Conservative dense factorization work: F*F*(F+N).
    pub maximum_dense_work: usize,
}

/// Current-controlled terminal block factored once on the original response.
/// A zero-total-current electrode floats and may carry nonzero local current.
pub struct PlanarMixedConductance<'s, 'm, C = ResolvedConductorState> {
    source: &'s PlanarConductanceMatrix<'m, C>,
    voltage_terminals: Vec<usize>,
    current_terminals: Vec<usize>,
    factor: Cholesky,
    scales: Vec<f64>,
}

type Result<T> = std::result::Result<T, ConductionError>;

fn invalid(reason: &'static str) -> ConductionError {
    ConductionError::Invalid(reason)
}

fn factor_error(error: FactorRunError) -> ConductionError {
    match error {
        FactorRunError::Cancelled => ConductionError::Cancelled,
        FactorRunError::Factor(_) => invalid("terminal conductance block is not positive definite"),
    }
}

fn add(sum: &mut f64, correction: &mut f64, value: f64) -> Result<()> {
    let adjusted = finite(value - *correction)?;
    let next = finite(*sum + adjusted)?;
    *correction = finite((next - *sum) - adjusted)?;
    *sum = next;
    Ok(())
}

fn terminal_current<C>(
    field: &ConductionField<'_, C>,
    vertices: &[usize],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let (mut value, mut correction) = (0.0, 0.0);
    for &vertex in vertices {
        poll(cancelled)?;
        add(
            &mut value,
            &mut correction,
            field.nodal_boundary_current_a[vertex]
                .ok_or_else(|| invalid("lost original terminal current reaction"))?,
        )?;
    }
    Ok(value)
}

impl<C: PlanarOhmicResponse> PlanarConduction<C> {
    /// Prepare I=G V for disjoint complete exterior-edge electrodes. Other
    /// edges insulate. Every original conducting component needs a terminal,
    /// and all terminals must connect through bulk, finite contacts or a shared
    /// ideal electrode. A terminal can join disconnected mesh pieces.
    ///
    /// The same scalar/tensor conductivity and native contact operator supply
    /// the basis, energy Gram matrix, terminal reactions and local heat. No
    /// pair-resistor fit or uniform electrode current density is assumed.
    pub fn conductance_matrix(
        &self,
        terminals: &[&[usize]],
        controls: PlanarConductanceControls,
        field_controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarConductanceMatrix<'_, C>> {
        poll(cancelled)?;
        let n = terminals.len();
        let m = n.saturating_sub(1);
        let vertices = self.mesh.vertices().len();
        let entries = n.checked_mul(n);
        let values = n
            .checked_add(4)
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
                "bounded terminal response and positive probe/tolerance required",
            ));
        }
        // Original trace admission rejects overlapping/point electrodes and
        // internal contact edges, and checks an anchor on every physical body.
        let zero =
            scalar_potential::terminal_boundary(&self.mesh, terminals, &vec![0.0; n], cancelled)?;
        self.validate_controls(zero.len(), field_controls)?;
        self.prescribed_boundary(&zero, cancelled)?;
        let contact_links: Vec<_> = self
            .contacts
            .iter()
            .map(|contact| [contact.first_vertices[0], contact.second_vertices[0]])
            .collect();
        let mut links = terminals.to_vec();
        links.extend(contact_links.iter().map(|pair| pair.as_slice()));
        let components = scalar_potential::terminal_components(&self.mesh, &links, cancelled)?;
        if terminals
            .iter()
            .any(|terminal| components[terminal[0]] != components[terminals[0][0]])
        {
            return Err(invalid(
                "terminals have no common bulk/contact/electrode connection",
            ));
        }
        let terminals: Vec<Vec<usize>> = terminals
            .iter()
            .map(|vertices| {
                let mut canonical = vertices.to_vec();
                canonical.sort_unstable();
                canonical
            })
            .collect();
        let views: Vec<_> = terminals.iter().map(Vec::as_slice).collect();
        let mut basis = Vec::with_capacity(m);
        let mut reactions = vec![0.0; n * m];
        for column in 0..m {
            poll(cancelled)?;
            let mut voltages = vec![0.0; n];
            voltages[column + 1] = controls.reference_voltage_v;
            let boundary =
                scalar_potential::terminal_boundary(&self.mesh, &views, &voltages, cancelled)?;
            let field = self.solve(&boundary, field_controls, cancelled)?;
            for (row, vertices) in terminals.iter().enumerate() {
                reactions[row * m + column] = finite(
                    terminal_current(&field, vertices, cancelled)? / controls.reference_voltage_v,
                )?;
            }
            // These probe potentials have no large common bias. Retain their
            // common terminal gauge, including bodies joined by an electrode;
            // independent native per-component reference offsets are not lost.
            basis.push(
                field
                    .voltage_v
                    .iter()
                    .map(|value| finite(value / controls.reference_voltage_v))
                    .collect::<Result<Vec<_>>>()?,
            );
        }
        let mut matrix_s = vec![0.0; n * n];
        for column in 0..m {
            let mut applied = vec![0.0; vertices];
            for (vertex, out) in applied.iter_mut().enumerate() {
                poll(cancelled)?;
                let (indices, values) = self.stiffness.row(vertex);
                let mut correction = 0.0;
                for (&index, &value) in indices.iter().zip(values) {
                    add(out, &mut correction, finite(value * basis[column][index])?)?;
                }
            }
            for row in 0..=column {
                let (mut value, mut correction) = (0.0, 0.0);
                for (&potential, &current) in basis[row].iter().zip(&applied) {
                    poll(cancelled)?;
                    add(&mut value, &mut correction, finite(potential * current)?)?;
                }
                matrix_s[(row + 1) * n + column + 1] = value;
                matrix_s[(column + 1) * n + row + 1] = value;
            }
        }
        // Complete the common-voltage nullspace exactly in the same symmetric
        // energy matrix, then audit it against independently recovered currents.
        for row in 1..n {
            let (mut sum, mut correction) = (0.0, 0.0);
            for column in 1..n {
                add(&mut sum, &mut correction, matrix_s[row * n + column])?;
            }
            matrix_s[row * n] = -sum;
            matrix_s[row] = -sum;
        }
        let (mut diagonal, mut correction) = (0.0, 0.0);
        for &value in &matrix_s[1..n] {
            add(&mut diagonal, &mut correction, -value)?;
        }
        matrix_s[0] = diagonal;
        let mut maximum_coefficient_residual_s = 0.0_f64;
        for row in 0..n {
            poll(cancelled)?;
            if !positive(matrix_s[row * n + row]) {
                return Err(invalid(
                    "terminal lacks representably positive dissipative connection",
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
                "terminal reactions, energy matrix or gauge do not agree",
            ));
        }
        let scales: Vec<_> = (1..n).map(|i| matrix_s[i * n + i].sqrt()).collect();
        let mut block = vec![0.0; m * m];
        for row in 0..m {
            poll(cancelled)?;
            for column in 0..m {
                block[row * m + column] =
                    finite((matrix_s[(row + 1) * n + column + 1] / scales[row]) / scales[column])?;
            }
        }
        poll(cancelled)?;
        let factor = cholesky_cancellable(&block, m, &mut *cancelled).map_err(factor_error)?;
        poll(cancelled)?;
        Ok(PlanarConductanceMatrix {
            problem: self,
            terminals,
            basis,
            matrix_s,
            factor,
            scales,
            maximum_coefficient_residual_s,
        })
    }
}

impl<'m, C: PlanarOhmicResponse> PlanarConductanceMatrix<'m, C> {
    /// Canonical vertex sets in the original terminal order.
    pub fn terminal_vertices(&self) -> &[Vec<usize>] {
        &self.terminals
    }
    /// Full row-major matrix [S], including the reference terminal: I=G V.
    pub fn conductance_s(&self) -> &[f64] {
        &self.matrix_s
    }
    pub fn field_solves(&self) -> usize {
        self.basis.len()
    }
    pub fn maximum_coefficient_residual_s(&self) -> f64 {
        self.maximum_coefficient_residual_s
    }

    /// Recover simultaneous terminal voltages without another field solve.
    /// Common bias is removed before all field/current/power arithmetic; the
    /// returned voltages are relative to original terminal zero. Original
    /// material field domains and native actual-amplitude residuals still apply.
    pub fn field_at_voltages(
        &self,
        voltages_v: &[f64],
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionField<'m, C>> {
        poll(cancelled)?;
        let n = self.terminals.len();
        if voltages_v.len() != n
            || voltages_v.iter().any(|v| !v.is_finite())
            || !positive(current_tolerance_a)
            || !positive(power_tolerance_w)
        {
            return Err(invalid(
                "complete finite terminal voltages and positive tolerances required",
            ));
        }
        let voltages = voltages_v
            .iter()
            .map(|v| finite(v - voltages_v[0]))
            .collect::<Result<Vec<_>>>()?;
        let mut potentials = vec![0.0; self.problem.mesh.vertices().len()];
        for (vertex, out) in potentials.iter_mut().enumerate() {
            poll(cancelled)?;
            let mut correction = 0.0;
            for (basis, &voltage) in self.basis.iter().zip(&voltages[1..]) {
                add(out, &mut correction, finite(basis[vertex] * voltage)?)?;
            }
        }
        let views: Vec<_> = self.terminals.iter().map(Vec::as_slice).collect();
        let boundary =
            scalar_potential::terminal_boundary(&self.problem.mesh, &views, &voltages, cancelled)?;
        let prescribed = self.problem.prescribed_boundary(&boundary, cancelled)?;
        let references = scalar_potential::component_references(
            &self.problem.stiffness,
            &prescribed,
            cancelled,
        )?;
        for (vertex, value) in potentials.iter_mut().enumerate() {
            poll(cancelled)?;
            // Restore exact authored Dirichlet values after superposition;
            // native residuals independently audit the resulting free vertices.
            *value = finite(prescribed[vertex].unwrap_or(*value) - references[vertex])?;
        }
        let field = self.problem.recover_field(
            potentials,
            references,
            &prescribed,
            vec![None; self.problem.mesh.vertices().len()],
            current_tolerance_a,
            power_tolerance_w,
            None,
            cancelled,
        )?;
        let mut currents = vec![0.0; n];
        for (row, out) in currents.iter_mut().enumerate() {
            poll(cancelled)?;
            let mut correction = 0.0;
            for (column, &voltage) in voltages.iter().enumerate() {
                add(
                    out,
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
        if voltages.iter().any(|&v| v != 0.0) && !positive(field.power_w) {
            return Err(invalid(
                "nonzero terminal voltage lost representable Joule work",
            ));
        }
        Ok(field)
    }

    /// Prescribe total currents INTO every original equipotential electrode.
    /// Their sum must vanish within the authored tolerance. Zero total current
    /// floats an electrode; it does not force every local surface current zero.
    /// The equilibrated grounded factor is reused, then native fields and I.V
    /// are independently recovered and checked at the actual drive amplitude.
    pub fn field_at_currents(
        &self,
        currents_a: &[f64],
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionField<'m, C>> {
        poll(cancelled)?;
        if currents_a.len() != self.terminals.len()
            || currents_a.iter().any(|v| !v.is_finite())
            || !positive(current_tolerance_a)
            || !positive(power_tolerance_w)
        {
            return Err(invalid(
                "complete finite terminal currents and positive tolerances required",
            ));
        }
        let (mut sum, mut correction) = (0.0, 0.0);
        for &current in currents_a {
            add(&mut sum, &mut correction, current)?;
        }
        if sum.abs() > current_tolerance_a {
            return Err(invalid("terminal currents must sum to zero"));
        }
        let mut rhs = currents_a[1..]
            .iter()
            .zip(&self.scales)
            .map(|(&current, &scale)| finite(current / scale))
            .collect::<Result<Vec<_>>>()?;
        poll(cancelled)?;
        self.factor
            .solve_cancellable(&mut rhs, &mut *cancelled)
            .map_err(factor_error)?;
        poll(cancelled)?;
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
        if currents_a.iter().any(|&v| v != 0.0) && !positive(field.power_w) {
            return Err(invalid(
                "nonzero terminal current lost representable Joule work",
            ));
        }
        Ok(field)
    }

    /// Factor an arbitrary nonempty proper subset of current electrodes once.
    /// The given order fixes future current input order; the complement's
    /// original order fixes voltage input order. At least one voltage electrode
    /// removes the sole connected terminal gauge.
    pub fn prepare_current_control(
        &self,
        current_terminals: &[usize],
        controls: PlanarMixedConductanceControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarMixedConductance<'_, 'm, C>> {
        poll(cancelled)?;
        let n = self.terminals.len();
        let f = current_terminals.len();
        let square = f.checked_mul(f);
        let work = square.and_then(|s| s.checked_mul(f.checked_add(n)?));
        let values = square.and_then(|s| s.checked_mul(2)?.checked_add(n.checked_mul(3)?));
        if f == 0
            || f >= n
            || f > controls.maximum_current_terminals
            || work.is_none_or(|w| w > controls.maximum_dense_work)
            || values.is_none_or(|v| v > controls.maximum_response_values)
            || current_terminals
                .iter()
                .enumerate()
                .any(|(i, &t)| t >= n || current_terminals[..i].contains(&t))
        {
            return Err(invalid(
                "bounded unique current terminals and a voltage terminal required",
            ));
        }
        let voltage_terminals = (0..n).filter(|t| !current_terminals.contains(t)).collect();
        let scales: Vec<_> = current_terminals
            .iter()
            .map(|&t| self.matrix_s[t * n + t].sqrt())
            .collect();
        let mut block = vec![0.0; f * f];
        for (row, &a) in current_terminals.iter().enumerate() {
            poll(cancelled)?;
            for (column, &b) in current_terminals.iter().enumerate() {
                block[row * f + column] =
                    finite((self.matrix_s[a * n + b] / scales[row]) / scales[column])?;
            }
        }
        poll(cancelled)?;
        let factor = cholesky_cancellable(&block, f, &mut *cancelled).map_err(factor_error)?;
        poll(cancelled)?;
        Ok(PlanarMixedConductance {
            source: self,
            voltage_terminals,
            current_terminals: current_terminals.to_vec(),
            factor,
            scales,
        })
    }

    fn admit_response(
        &self,
        field: &ConductionField<'_, C>,
        voltages: &[f64],
        currents: &[f64],
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        let mut current_residual_a = 0.0_f64;
        let (mut work, mut correction) = (0.0, 0.0);
        for ((vertices, &current), &voltage) in self.terminals.iter().zip(currents).zip(voltages) {
            poll(cancelled)?;
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
            return Err(ConductionError::Balance {
                current_residual_a,
                maximum_free_current_residual_a: field.maximum_free_current_residual_a,
                power_residual_w,
            });
        }
        poll(cancelled)
    }
}

impl<'m, C: PlanarOhmicResponse> PlanarMixedConductance<'_, 'm, C> {
    pub fn voltage_terminals(&self) -> &[usize] {
        &self.voltage_terminals
    }
    pub fn current_terminals(&self) -> &[usize] {
        &self.current_terminals
    }

    /// Solve G_FF V_F = I_F - G_FD V_D with the retained factor. Voltage-source
    /// reactions supply the remaining current/work. Common voltage bias is
    /// removed before arithmetic; returned potentials use original terminal zero.
    pub fn field_at_controls(
        &self,
        voltages_v: &[f64],
        currents_a: &[f64],
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionField<'m, C>> {
        poll(cancelled)?;
        if voltages_v.len() != self.voltage_terminals.len()
            || currents_a.len() != self.current_terminals.len()
            || voltages_v.iter().chain(currents_a).any(|v| !v.is_finite())
            || !positive(current_tolerance_a)
            || !positive(power_tolerance_w)
        {
            return Err(invalid(
                "complete finite mixed controls and positive tolerances required",
            ));
        }
        let n = self.source.terminals.len();
        let mut voltages = vec![0.0; n];
        for (&terminal, &value) in self.voltage_terminals.iter().zip(voltages_v) {
            voltages[terminal] = finite(value - voltages_v[0])?;
        }
        let mut rhs = Vec::with_capacity(currents_a.len());
        for ((&terminal, &current), &scale) in self
            .current_terminals
            .iter()
            .zip(currents_a)
            .zip(&self.scales)
        {
            poll(cancelled)?;
            let (mut load, mut correction) = (current, 0.0);
            for &fixed in &self.voltage_terminals {
                add(
                    &mut load,
                    &mut correction,
                    finite(-self.source.matrix_s[terminal * n + fixed] * voltages[fixed])?,
                )?;
            }
            rhs.push(finite(load / scale)?);
        }
        poll(cancelled)?;
        self.factor
            .solve_cancellable(&mut rhs, &mut *cancelled)
            .map_err(factor_error)?;
        poll(cancelled)?;
        for ((&terminal, value), &scale) in self.current_terminals.iter().zip(rhs).zip(&self.scales)
        {
            voltages[terminal] = finite(value / scale)?;
        }
        let field = self.source.field_at_voltages(
            &voltages,
            current_tolerance_a,
            power_tolerance_w,
            cancelled,
        )?;
        let mut currents = self
            .source
            .terminals
            .iter()
            .map(|vertices| terminal_current(&field, vertices, cancelled))
            .collect::<Result<Vec<_>>>()?;
        for (&terminal, &current) in self.current_terminals.iter().zip(currents_a) {
            currents[terminal] = current;
        }
        let gauge = voltages[0];
        for voltage in &mut voltages {
            *voltage = finite(*voltage - gauge)?;
        }
        self.source.admit_response(
            &field,
            &voltages,
            &currents,
            current_tolerance_a,
            power_tolerance_w,
            cancelled,
        )?;
        if currents_a.iter().any(|&v| v != 0.0) && !positive(field.power_w) {
            return Err(invalid(
                "nonzero mixed terminal current lost representable Joule work",
            ));
        }
        Ok(field)
    }
}
