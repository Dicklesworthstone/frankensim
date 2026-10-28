//! Pressure-continuous star junctions in the ORIGINAL DG energy coordinates.
//!
//! With H=[-B_bore,+B_branch], the source equations are
//! x'=Lx+Hq+Bu and H^T x=0. Eliminating only q gives the orthogonal
//! projector P=I-H(H^T H)^-1 H^T. PLP and PB extend the same constrained
//! dynamics to the full original coordinate space without discarding states.
//! Every physical loss and pressure observation still uses its original owner.
use super::*;

pub(crate) struct DgContinuousOwner {
    pub offset: usize,
    pub dimension: usize,
    /// Original local -sym(L), before junction projection.
    pub loss: Vec<f64>,
}

pub(crate) struct DgContinuousJunction {
    pub port: PortId,
    pub flow_column: Vec<f64>,
    pub inlet_flow_factor: f64,
    /// Coefficients of OUTWARD terminal flows in explicit source channel order.
    pub terminal_flow_factors: Vec<f64>,
    pub bore_pressure_column: Vec<f64>,
    pub branch_pressure_column: Vec<f64>,
}

pub(crate) struct DgContinuousNetwork {
    pub dimension: usize,
    pub generator: Vec<f64>,
    pub inlet_column: Vec<f64>,
    pub inlet_pressure_column: Vec<f64>,
    pub terminal_pressure_columns: Vec<Vec<f64>>,
    pub terminal_input_columns: Vec<Vec<f64>>,
    pub owners: Vec<DgContinuousOwner>,
    pub junctions: Vec<DgContinuousJunction>,
}

fn poll(gate: &fs_exec::CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() {
        Err(AcousticRealizeError::Cancelled)
    } else {
        Ok(())
    }
}

fn trace(
    view: &DgContinuous<'_>,
    id: PortId,
    offset: usize,
    dimension: usize,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let (_, original) = view
        .ports()
        .iter()
        .find(|(port, _)| *port == id)
        .ok_or_else(|| invalid("continuous network lost an original pressure trace"))?;
    let mut column = vec![0.0; dimension];
    column[offset..offset + original.len()].copy_from_slice(original);
    Ok(column)
}

fn apply(matrix: &[f64], column: &[f64]) -> Vec<f64> {
    matrix
        .chunks_exact(column.len())
        .map(|row| row.iter().zip(column).map(|(a, x)| a * x).sum())
        .collect()
}

fn sandwich(
    projection: &[f64],
    matrix: &[f64],
    n: usize,
    gate: &fs_exec::CancelGate,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let mut right = vec![0.0; n * n];
    let mut out = vec![0.0; n * n];
    for i in 0..n {
        poll(gate)?;
        for k in 0..n {
            for j in 0..n {
                right[i * n + j] += matrix[i * n + k] * projection[k * n + j];
            }
        }
    }
    for i in 0..n {
        poll(gate)?;
        for k in 0..n {
            for j in 0..n {
                out[i * n + j] += projection[i * n + k] * right[k * n + j];
            }
        }
    }
    Ok(out)
}

impl DgPortLoad {
    pub(crate) fn continuous_energy_coordinates(&self) -> Result<Vec<f64>, AcousticRealizeError> {
        self.admit_continuous_network()?;
        let mut state = self
            .duct
            .continuous_operator()
            .energy_coordinates(&self.state)
            .map_err(AcousticRealizeError::Duct)?
            .to_vec();
        for branch in &self.branches {
            state.extend_from_slice(
                branch
                    .duct()
                    .continuous_operator()
                    .energy_coordinates(branch.state())
                    .map_err(AcousticRealizeError::Duct)?,
            );
        }
        Ok(state)
    }

    pub(super) fn admit_continuous_network(&self) -> Result<usize, AcousticRealizeError> {
        if !self.terminal_channels.is_empty() || self.has_driven_mouth() {
            return Err(invalid(
                "continuous DG network must retain finite gap and driven-mouth states before admitting them",
            ));
        }
        let original = self.duct.continuous_operator();
        if original.ports().len() != 1 + self.terminal_ports().len()
            || self
                .terminal
                .as_ref()
                .is_some_and(|t| t.ports() != self.terminal_ports().len())
        {
            return Err(invalid(
                "continuous network must retain every original terminal port",
            ));
        }
        let mut dimension = original.dimension();
        for branch in &self.branches {
            if branch.propagator().dt().to_bits() != self.propagator.dt().to_bits()
                || branch.state().step_count() != self.state.step_count()
                || branch.duct().continuous_operator().ports().len() != 2
                || !self.terminal_ports().contains(&branch.port())
            {
                return Err(invalid(
                    "continuous branch lost its original ports or shared accepted clock",
                ));
            }
            dimension = dimension
                .checked_add(branch.duct().total_dofs())
                .ok_or_else(|| invalid("continuous DG network dimension overflow"))?;
        }
        Ok(dimension)
    }

    pub(crate) fn prepare_continuous_network(
        &self,
        maximum_dense_bytes: usize,
        gate: &fs_exec::CancelGate,
    ) -> Result<DgContinuousNetwork, AcousticRealizeError> {
        poll(gate)?;
        let n = self.admit_continuous_network()?;
        if n == 0
            || n > 256
            || n.checked_mul(n)
                .and_then(|v| v.checked_mul(128))
                .is_none_or(|v| v > maximum_dense_bytes)
        {
            return Err(invalid(
                "complete original continuous network exceeds dense allowance",
            ));
        }
        let main = self.duct.continuous_operator();
        let mut source = vec![0.0; n * n];
        let mut owners = Vec::with_capacity(1 + self.branches.len());
        let mut offset = 0;
        for duct in std::iter::once(&self.duct).chain(self.branches.iter().map(|b| b.duct())) {
            poll(gate)?;
            let view = duct.continuous_operator();
            let d = view.dimension();
            let a = view.generator_per_s();
            let mut loss = vec![0.0; d * d];
            for i in 0..d {
                source[(offset + i) * n + offset..(offset + i) * n + offset + d]
                    .copy_from_slice(&a[i * d..(i + 1) * d]);
                for j in 0..d {
                    loss[i * d + j] = -0.5 * (a[i * d + j] + a[j * d + i]);
                }
            }
            owners.push(DgContinuousOwner {
                offset,
                dimension: d,
                loss,
            });
            offset += d;
        }
        let inlet_pressure_column = trace(&main, self.inlet, 0, n)?;
        let terminal_pressure_columns: Vec<_> = self
            .terminal_ports()
            .iter()
            .map(|&id| {
                if let Some((index, branch)) = self
                    .branches
                    .iter()
                    .enumerate()
                    .find(|(_, b)| b.port() == id)
                {
                    let DgBoundary::Port(outlet) = branch.duct().outlet() else {
                        return Err(invalid("continuous branch lost its original outlet"));
                    };
                    trace(
                        &branch.duct().continuous_operator(),
                        outlet,
                        owners[index + 1].offset,
                        n,
                    )
                } else {
                    trace(&main, id, 0, n)
                }
            })
            .collect::<Result<_, _>>()?;
        if self.branches.is_empty() {
            return Ok(DgContinuousNetwork {
                dimension: n,
                generator: source,
                inlet_column: inlet_pressure_column.clone(),
                inlet_pressure_column,
                terminal_input_columns: terminal_pressure_columns.clone(),
                terminal_pressure_columns,
                owners,
                junctions: Vec::new(),
            });
        }

        // Normalize each constraint before using the existing conditioned LU.
        // Distinct branch inlet blocks make this complete Gram matrix SPD even
        // when bore-side pressure traces overlap. Scalar independent solves
        // would lose precisely that cross-junction coupling.
        let b = self.branches.len();
        let mut constraints = Vec::with_capacity(b);
        let mut scales = Vec::with_capacity(b);
        let mut junctions = Vec::with_capacity(b);
        for (index, branch) in self.branches.iter().enumerate() {
            poll(gate)?;
            let DgBoundary::Port(inlet) = branch.duct().inlet() else {
                return Err(invalid("continuous branch lost its original inlet"));
            };
            let bore = trace(&main, branch.port(), 0, n)?;
            let branch_inlet = trace(
                &branch.duct().continuous_operator(),
                inlet,
                owners[index + 1].offset,
                n,
            )?;
            let mut h: Vec<_> = branch_inlet.iter().zip(&bore).map(|(a, b)| a - b).collect();
            let maximum = h.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
            let norm = maximum
                * fs_math::det::sqrt(
                    h.iter()
                        .map(|v| {
                            let scaled = v / maximum;
                            scaled * scaled
                        })
                        .sum(),
                );
            if !norm.is_finite() || norm <= 0.0 {
                return Err(invalid(
                    "continuous junction needs a finite nonzero original pressure constraint",
                ));
            }
            for value in &mut h {
                *value /= norm;
            }
            constraints.push(h);
            scales.push(norm);
            junctions.push(DgContinuousJunction {
                port: branch.port(),
                flow_column: vec![0.0; n],
                inlet_flow_factor: 0.0,
                terminal_flow_factors: Vec::new(),
                bore_pressure_column: bore,
                branch_pressure_column: branch_inlet,
            });
        }
        let mut gram = vec![0.0; b * b];
        for i in 0..b {
            poll(gate)?;
            for j in 0..b {
                gram[i * b + j] = constraints[i]
                    .iter()
                    .zip(&constraints[j])
                    .map(|(a, b)| a * b)
                    .sum();
            }
        }
        let solve = super::super::terminal::TerminalSolve::new(&gram, b)?;
        let mut dual = vec![vec![0.0; n]; b];
        for col in 0..n {
            poll(gate)?;
            let rhs: Vec<_> = constraints.iter().map(|h| h[col]).collect();
            for (row, value) in solve.apply(&rhs)?.into_iter().enumerate() {
                dual[row][col] = value;
            }
        }
        let mut projection = vec![0.0; n * n];
        for i in 0..n {
            poll(gate)?;
            for j in 0..=i {
                let ij: f64 = (0..b).map(|k| constraints[k][i] * dual[k][j]).sum();
                let ji: f64 = (0..b).map(|k| constraints[k][j] * dual[k][i]).sum();
                let value = if i == j { 1.0 } else { 0.0 } - 0.5 * (ij + ji);
                projection[i * n + j] = value;
                projection[j * n + i] = value;
            }
        }
        let mut constraint_residual = 0.0_f64;
        for h in &constraints {
            poll(gate)?;
            for j in 0..n {
                let value: f64 = (0..n).map(|i| h[i] * projection[i * n + j]).sum();
                constraint_residual = constraint_residual.max(value.abs());
            }
        }
        if !constraint_residual.is_finite() || constraint_residual > 256.0 * f64::EPSILON * n as f64
        {
            return Err(invalid(
                "continuous pressure-junction projection fails its original constraints",
            ));
        }
        for (index, junction) in junctions.iter_mut().enumerate() {
            poll(gate)?;
            for col in 0..n {
                junction.flow_column[col] = -(0..n)
                    .map(|row| dual[index][row] * source[row * n + col])
                    .sum::<f64>()
                    / scales[index];
            }
            junction.inlet_flow_factor = -dual[index]
                .iter()
                .zip(&inlet_pressure_column)
                .map(|(a, b)| a * b)
                .sum::<f64>()
                / scales[index];
            junction.terminal_flow_factors = terminal_pressure_columns
                .iter()
                .map(|column| {
                    dual[index]
                        .iter()
                        .zip(column)
                        .map(|(a, b)| a * b)
                        .sum::<f64>()
                        / scales[index]
                })
                .collect();
        }

        // Project the physical skew and dissipative parts separately. Explicit
        // symmetry preserves an exactly skew lossless extension in floating
        // point; no damping is introduced on a pressure constraint.
        let mut skew = vec![0.0; n * n];
        let mut loss = vec![0.0; n * n];
        for i in 0..n {
            poll(gate)?;
            for j in 0..n {
                skew[i * n + j] = 0.5 * (source[i * n + j] - source[j * n + i]);
                loss[i * n + j] = -0.5 * (source[i * n + j] + source[j * n + i]);
            }
        }
        let skew = sandwich(&projection, &skew, n, gate)?;
        let loss = sandwich(&projection, &loss, n, gate)?;
        let mut generator = vec![0.0; n * n];
        for i in 0..n {
            poll(gate)?;
            for j in 0..n {
                generator[i * n + j] = 0.5 * (skew[i * n + j] - skew[j * n + i])
                    - 0.5 * (loss[i * n + j] + loss[j * n + i]);
            }
        }
        let inlet_column = apply(&projection, &inlet_pressure_column);
        let terminal_input_columns: Vec<_> = terminal_pressure_columns
            .iter()
            .map(|p| apply(&projection, p))
            .collect();
        if generator
            .iter()
            .chain(inlet_column.iter())
            .chain(terminal_input_columns.iter().flatten())
            .chain(junctions.iter().flat_map(|j| {
                j.flow_column
                    .iter()
                    .chain(std::iter::once(&j.inlet_flow_factor))
                    .chain(&j.terminal_flow_factors)
            }))
            .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "continuous junction composition left its finite source domain",
            ));
        }
        poll(gate)?;
        Ok(DgContinuousNetwork {
            dimension: n,
            generator,
            inlet_column,
            inlet_pressure_column,
            terminal_pressure_columns,
            terminal_input_columns,
            owners,
            junctions,
        })
    }
}
