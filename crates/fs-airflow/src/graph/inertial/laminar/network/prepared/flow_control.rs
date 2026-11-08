//! Required held pressures for actual mean flows through original network ports.
use super::*;
use fs_la::factor::{Cholesky, FactorRunError, cholesky_cancellable};
use fs_math::dd::Dd;
use fs_qty::VolumetricFlowRate;

/// One original power-conjugate port; indices follow the authored topology.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LaminarNetworkFlowPort {
    /// Index in `topology().pressure_nodes`; positive flow enters the network.
    Reservoir(usize),
    /// Original branch index; positive flow follows its authored direction.
    BranchHead(usize),
}

/// Explicit absolute and relative tolerance for each actual interval-mean flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminarNetworkFlowTolerance {
    pub absolute_m3_s: f64,
    pub relative: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct LaminarNetworkFlowControlStep {
    /// Original trajectory, independently integrated heat and physical work.
    pub network: LaminarNetworkStep,
    pub maximum_flow_residual_m3_s: f64,
}

/// Original exact dynamics plus one prepared, scaled pressure-response factor.
/// Warm calls reuse every buffer and perform no new matrix factorization.
pub struct PreparedLaminarNetworkFlowControl<'a, M> {
    inner: PreparedLaminarNetwork<'a, M>,
    ports: Vec<LaminarNetworkFlowPort>,
    slots: Vec<usize>,
    // Row-major complete original input response to selected physical pressure
    // increments. A controlled reference reservoir moves the OTHER gauges.
    pressure_columns: Vec<f64>,
    mobility: Vec<f64>,
    scales: Vec<f64>,
    factor: Cholesky,
    base_input: Vec<f64>,
    candidate_input: Vec<f64>,
    accepted_input: Vec<f64>,
    pressures: Vec<Pressure>,
    free: Vec<f64>,
    pending_offsets: Vec<f64>,
    offsets: Vec<f64>,
    correction: Vec<f64>,
    pending_flows: Vec<f64>,
    flows: Vec<f64>,
    last: Option<LaminarNetworkFlowControlStep>,
}

impl<M: FluidTransport> LaminarNetwork<M> {
    /// Prepare pressures that deliver prescribed mean flows over one interval.
    ///
    /// Each controlled reservoir pressure or directed branch head receives an
    /// unknown held increment. Other supplied pressures/heads stay fixed. The
    /// true interval average includes the original initial momentum and every
    /// radial-memory state; delivered volume is that average times duration.
    /// Neither endpoint flow nor steady Poiseuille resistance replaces it.
    ///
    /// Ports must be distinct and their physical flow observations independent.
    /// Complete constrained flow rank, symmetric positive response and scaled
    /// Cholesky pivots are checked. Redundant or numerically unresolved port
    /// combinations refuse, including attempts to control both ends of one
    /// incompressible branch independently. The original bounded network and
    /// fixed-clock preparation limits still apply.
    pub fn prepare_flow_control(
        &self,
        cx: &Cx<'_>,
        duration_s: f64,
        ports: &[LaminarNetworkFlowPort],
    ) -> Result<PreparedLaminarNetworkFlowControl<'_, M>, LaminarPipeError> {
        checkpoint(cx)?;
        let (n, m, k) = (self.system.state_dim(), self.port_count(), ports.len());
        let boundary = self.topology.pressure_nodes.len();
        if k == 0 || k > self.basis.columns {
            return Err(invalid(
                "controlled flows must fit the original independent mean-flow rank",
            ));
        }
        let mut slots = Vec::with_capacity(k);
        for &port in ports {
            let slot = match port {
                LaminarNetworkFlowPort::Reservoir(p) if p < boundary => p,
                LaminarNetworkFlowPort::BranchHead(e) if e < self.branches.len() => boundary + e,
                _ => {
                    return Err(invalid(
                        "controlled flow port is outside the original network",
                    ));
                }
            };
            if slots.contains(&slot) {
                return Err(invalid("controlled network ports must be distinct"));
            }
            slots.push(slot);
        }
        let mut columns = vec![0.0; m * k];
        for (column, &slot) in slots.iter().enumerate() {
            if slot < boundary && self.basis.pressure_reference[slot] == slot {
                for row in 0..boundary {
                    if row != slot && self.basis.pressure_reference[row] == slot {
                        columns[row * k + column] = -1.0;
                    }
                }
            } else {
                columns[slot * k + column] = 1.0;
            }
        }
        let inner = self.prepare(cx, duration_s)?;
        let mut mobility = vec![0.0; k * k];
        let mut scales = vec![0.0; k];
        for row in 0..k {
            checkpoint(cx)?;
            let (_, input) = inner.native.average_output_map(slots[row])?;
            for column in 0..k {
                mobility[row * k + column] =
                    sum((0..m)
                        .map(|p| Dd::from_f64(input[p]) * Dd::from_f64(columns[p * k + column])));
            }
            let diagonal = mobility[row * k + row];
            if !positive(diagonal) {
                return Err(invalid(
                    "controlled flow has no positive finite pressure mobility",
                ));
            }
            scales[row] = det::sqrt(diagonal);
        }
        let mut scaled = vec![0.0; k * k];
        for row in 0..k {
            for column in 0..k {
                scaled[row * k + column] =
                    (mobility[row * k + column] / scales[row]) / scales[column];
            }
        }
        let tolerance = ROUNDOFF * (n + m + k) as f64;
        for row in 0..k {
            checkpoint(cx)?;
            for column in 0..k {
                let (a, b) = (scaled[row * k + column], scaled[column * k + row]);
                if !a.is_finite()
                    || !b.is_finite()
                    || (a - b).abs() > tolerance * a.abs().max(b.abs()).max(1.0)
                {
                    return Err(invalid(
                        "original coupled pressure mobility is numerically unresolved",
                    ));
                }
            }
        }
        // Retain the actual lower triangle, never an averaged/fitted mobility.
        // Every warm solve is checked against the COMPLETE original response.
        let factor =
            cholesky_cancellable(&scaled, k, || cx.checkpoint().is_err()).map_err(|error| {
                match error {
                    FactorRunError::Cancelled => {
                        LaminarPipeError::Execution(super::super::super::GraphError::Interrupted)
                    }
                    FactorRunError::Factor(_) => invalid(
                        "controlled pressure responses are not independent positive directions",
                    ),
                }
            })?;
        if (0..k).any(|i| !positive(factor.l(i, i)) || factor.l(i, i).powi(2) <= tolerance) {
            return Err(invalid(
                "controlled pressure directions cannot be resolved independently",
            ));
        }
        checkpoint(cx)?;
        Ok(PreparedLaminarNetworkFlowControl {
            inner,
            ports: ports.to_vec(),
            slots,
            pressure_columns: columns,
            mobility,
            scales,
            factor,
            base_input: vec![0.0; m],
            candidate_input: vec![0.0; m],
            accepted_input: vec![0.0; m],
            pressures: vec![Pressure::new(0.0); m],
            free: vec![0.0; k],
            pending_offsets: vec![0.0; k],
            offsets: vec![0.0; k],
            correction: vec![0.0; k],
            pending_flows: vec![0.0; k],
            flows: vec![0.0; k],
            last: None,
        })
    }
}

impl<M: FluidTransport> PreparedLaminarNetworkFlowControl<'_, M> {
    pub fn duration_s(&self) -> f64 {
        self.inner.duration_s()
    }
    pub fn ports(&self) -> &[LaminarNetworkFlowPort] {
        &self.ports
    }
    pub fn branch_steps(&self) -> &[LaminarNetworkBranchStep] {
        self.inner.branch_steps()
    }
    pub fn boundary_volume_flows_m3_s(&self) -> &[f64] {
        self.inner.boundary_volume_flows_m3_s()
    }
    /// Actual means in controlled-port order; empty until first acceptance.
    pub fn controlled_volume_flows_m3_s(&self) -> &[f64] {
        if self.last.is_some() {
            &self.flows
        } else {
            &[]
        }
    }
    /// Required increments relative to the supplied original pressures/heads.
    pub fn pressure_offsets_pa(&self) -> &[f64] {
        if self.last.is_some() {
            &self.offsets
        } else {
            &[]
        }
    }
    /// Accepted complete native pressure inputs: reservoirs relative to each
    /// component's original reference, then directed branch heads. Keeping this
    /// gauge avoids losing a small solved pressure beside huge absolute values.
    pub fn gauged_port_pressures_pa(&self) -> &[f64] {
        if self.last.is_some() {
            &self.accepted_input
        } else {
            &[]
        }
    }
    pub const fn last_step(&self) -> Option<LaminarNetworkFlowControlStep> {
        self.last
    }

    /// Deliver requested signed interval-mean flows with solved held pressures.
    ///
    /// `reservoir_pressures` and `branch_heads` provide original baseline loads;
    /// only selected physical sources receive a solved increment. Flow is NOT
    /// prescribed pointwise in time: its integral is the volume delivered during
    /// this window, while the original fluid/profile inertia determines the
    /// complete trajectory. Original whole-interval Reynolds and energy/heat/
    /// continuity audits still gate acceptance. Unrepresentable or unresolved
    /// pressures/targets refuse without changing state, clock or public reports.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &mut self,
        cx: &Cx<'_>,
        state: &mut LaminarNetworkState,
        reservoir_pressures: &[Pressure],
        branch_heads: &[Pressure],
        targets: &[VolumetricFlowRate],
        tolerance: LaminarNetworkFlowTolerance,
    ) -> Result<LaminarNetworkFlowControlStep, LaminarPipeError> {
        checkpoint(cx)?;
        self.inner.network.admit_owner(state)?;
        let k = self.slots.len();
        if targets.len() != k
            || targets.iter().any(|q| !q.value().is_finite())
            || !tolerance.absolute_m3_s.is_finite()
            || tolerance.absolute_m3_s < 0.0
            || !tolerance.relative.is_finite()
            || !(0.0..1.0).contains(&tolerance.relative)
        {
            return Err(invalid(
                "complete finite mean-flow targets and nonnegative bounded tolerances required",
            ));
        }
        self.inner.network.fill_inputs(
            reservoir_pressures,
            Some(branch_heads),
            &mut self.base_input,
        )?;
        for row in 0..k {
            checkpoint(cx)?;
            let (initial, input) = self.inner.native.average_output_map(self.slots[row])?;
            self.free[row] = sum(initial
                .iter()
                .zip(&state.coordinates)
                .map(|(&a, &x)| Dd::from_f64(a) * Dd::from_f64(x))
                .chain(
                    input
                        .iter()
                        .zip(&self.base_input)
                        .map(|(&a, &u)| Dd::from_f64(a) * Dd::from_f64(u)),
                ));
            self.pending_offsets[row] =
                difference(targets[row].value(), self.free[row]) / self.scales[row];
        }
        if self
            .free
            .iter()
            .chain(&self.pending_offsets)
            .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "controlled mean-flow history or pressure load is unrepresentable",
            ));
        }
        self.factor.solve(&mut self.pending_offsets);
        for (pressure, scale) in self.pending_offsets.iter_mut().zip(&self.scales) {
            *pressure /= scale;
        }
        // Two bounded refinements retain the full unsymmetrized original map.
        // Cholesky triangular solves allocate no warm scratch.
        for _ in 0..2 {
            checkpoint(cx)?;
            for row in 0..k {
                let mut residual =
                    Dd::from_f64(targets[row].value()) - Dd::from_f64(self.free[row]);
                for column in 0..k {
                    residual = residual
                        - Dd::from_f64(self.mobility[row * k + column])
                            * Dd::from_f64(self.pending_offsets[column]);
                }
                self.correction[row] = residual.to_f64() / self.scales[row];
            }
            self.factor.solve(&mut self.correction);
            for row in 0..k {
                self.pending_offsets[row] = (Dd::from_f64(self.pending_offsets[row])
                    + Dd::from_f64(self.correction[row] / self.scales[row]))
                .to_f64();
            }
        }
        if self.pending_offsets.iter().any(|x| !x.is_finite()) {
            return Err(invalid("required network pressures are unrepresentable"));
        }
        for port in 0..self.base_input.len() {
            let mut pressure = Dd::from_f64(self.base_input[port]);
            for column in 0..k {
                pressure = pressure
                    + Dd::from_f64(self.pressure_columns[port * k + column])
                        * Dd::from_f64(self.pending_offsets[column]);
            }
            self.candidate_input[port] = pressure.to_f64();
            if !self.candidate_input[port].is_finite() {
                return Err(invalid("controlled gauged pressure is unrepresentable"));
            }
            self.pressures[port] = Pressure::new(self.candidate_input[port]);
        }
        let boundary = self.inner.network.topology.pressure_nodes.len();
        let network = self.inner.preview(
            cx,
            state,
            &self.pressures[..boundary],
            Some(&self.pressures[boundary..]),
        )?;
        let mut maximum_flow_residual_m3_s = 0.0_f64;
        for row in 0..k {
            checkpoint(cx)?;
            let actual = self.inner.output[self.slots[row]];
            let residual = difference(actual, targets[row].value()).abs();
            let allowance = tolerance
                .relative
                .mul_add(targets[row].value().abs(), tolerance.absolute_m3_s);
            if !actual.is_finite() || !allowance.is_finite() || residual > allowance {
                return Err(invalid(
                    "original network trajectory misses a prescribed mean flow",
                ));
            }
            self.pending_flows[row] = actual;
            maximum_flow_residual_m3_s = maximum_flow_residual_m3_s.max(residual);
        }
        let report = LaminarNetworkFlowControlStep {
            network,
            maximum_flow_residual_m3_s,
        };
        checkpoint(cx)?;
        self.inner.publish(state, network);
        core::mem::swap(&mut self.offsets, &mut self.pending_offsets);
        core::mem::swap(&mut self.flows, &mut self.pending_flows);
        core::mem::swap(&mut self.accepted_input, &mut self.candidate_input);
        self.last = Some(report);
        Ok(report)
    }
}

fn invalid(message: &'static str) -> LaminarPipeError {
    LaminarPipeError::Invalid(message)
}
fn difference(a: f64, b: f64) -> f64 {
    (Dd::from_f64(a) - Dd::from_f64(b)).to_f64()
}
fn sum(values: impl Iterator<Item = Dd>) -> f64 {
    values
        .fold(Dd::from_f64(0.0), |sum, value| sum + value)
        .to_f64()
}

#[cfg(test)]
mod tests;
