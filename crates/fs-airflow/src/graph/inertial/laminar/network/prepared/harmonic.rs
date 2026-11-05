//! Continuous pressure/head response of the already projected laminar network.
//!
//! The original port-Hamiltonian coordinates retain every independent mean-flow
//! mode and every pipe's viscous-memory arm. Harmonic queries reuse that system
//! directly under `exp(-i omega t)`; no steady resistance network or isolated
//! impedance assembly replaces it. Pressure nodes are gauge-normalized exactly
//! as in time-domain stepping. Outputs are pressure-node flows first, followed by
//! one original branch flow per explicit head port.

use super::*;
use fs_math::{c64::C64, det};
use fs_phs::harmonic::{HarmonicBudget, HarmonicError, HarmonicModel, HarmonicResponse};

impl<M: FluidTransport> LaminarNetwork<M> {
    /// Prepare the complete native continuous-frequency model once. External
    /// composition may reuse it, but each queried frequency remains subject to
    /// the original pipes' wall-memory bands; the checked one-shot methods below
    /// perform that admission and the network continuity audit automatically.
    pub fn harmonic_model(
        &self,
        cx: &Cx<'_>,
        budget: HarmonicBudget,
    ) -> Result<HarmonicModel, HarmonicError> {
        poll(cx)?;
        let model = HarmonicModel::from_ode(&self.system, budget, || cx.checkpoint().is_err())?;
        poll(cx)?;
        Ok(model)
    }

    /// Gauge-normalize reservoir pressure phasors and append one signed pressure
    /// rise per branch. Positive heads drive each authored branch direction.
    /// The returned port order exactly matches [`Self::system`].
    pub fn harmonic_input_phasors(
        &self,
        reservoir_pressures_pa: &[C64],
        branch_heads_pa: Option<&[C64]>,
    ) -> Result<Vec<C64>, HarmonicError> {
        let boundary = self.topology.pressure_nodes.len();
        let edges = self.branches.len();
        if reservoir_pressures_pa.len() != boundary
            || branch_heads_pa.is_some_and(|heads| heads.len() != edges)
            || reservoir_pressures_pa.iter().any(|value| !finite(*value))
            || branch_heads_pa.is_some_and(|heads| heads.iter().any(|value| !finite(*value)))
        {
            return Err(invalid(
                "complete finite reservoir-pressure and branch-head phasors required",
            ));
        }
        let mut input = vec![C64::ZERO; boundary + edges];
        for port in 0..boundary {
            input[port] = reservoir_pressures_pa[port]
                - reservoir_pressures_pa[self.basis.pressure_reference[port]];
        }
        for edge in 0..edges {
            input[boundary + edge] = branch_heads_pa.map_or(C64::ZERO, |heads| heads[edge]);
        }
        if input.iter().any(|value| !finite(*value)) {
            return Err(invalid("gauge-normalized pressure phasors are not representable"));
        }
        Ok(input)
    }

    /// Solve one reservoir-pressure excitation with zero explicit branch heads.
    pub fn harmonic_response(
        &self,
        cx: &Cx<'_>,
        omega_rad_s: f64,
        reservoir_pressures_pa: &[C64],
        budget: HarmonicBudget,
    ) -> Result<HarmonicResponse, HarmonicError> {
        self.harmonic_response_impl(cx, omega_rad_s, reservoir_pressures_pa, None, budget)
    }

    /// Solve one reservoir-pressure and explicit branch-head excitation.
    pub fn harmonic_response_with_heads(
        &self,
        cx: &Cx<'_>,
        omega_rad_s: f64,
        reservoir_pressures_pa: &[C64],
        branch_heads_pa: &[C64],
        budget: HarmonicBudget,
    ) -> Result<HarmonicResponse, HarmonicError> {
        self.harmonic_response_impl(
            cx,
            omega_rad_s,
            reservoir_pressures_pa,
            Some(branch_heads_pa),
            budget,
        )
    }

    fn harmonic_response_impl(
        &self,
        cx: &Cx<'_>,
        omega_rad_s: f64,
        reservoir_pressures_pa: &[C64],
        branch_heads_pa: Option<&[C64]>,
        budget: HarmonicBudget,
    ) -> Result<HarmonicResponse, HarmonicError> {
        poll(cx)?;
        // Keep every original pipe's sampled wall-memory band authoritative.
        for pipe in &self.branches {
            pipe.impedance_at(omega_rad_s).map_err(|_| {
                invalid("frequency lies outside an original pipe wall-memory admission")
            })?;
            poll(cx)?;
        }
        let input = self.harmonic_input_phasors(reservoir_pressures_pa, branch_heads_pa)?;
        let tolerance = budget.relative_residual_tolerance
            + ROUNDOFF * self.system.state_dim().max(input.len()) as f64;
        let model = self.harmonic_model(cx, budget)?;
        let factor = model.factor(omega_rad_s, || cx.checkpoint().is_err())?;
        let response = factor.response(&input, || cx.checkpoint().is_err())?;
        self.audit_harmonic_response(&response, tolerance, cx)?;
        poll(cx)?;
        Ok(response)
    }

    fn audit_harmonic_response(
        &self,
        response: &HarmonicResponse,
        tolerance: f64,
        cx: &Cx<'_>,
    ) -> Result<(), HarmonicError> {
        let boundary = self.topology.pressure_nodes.len();
        if response.state.len() != self.system.state_dim()
            || response.output.len() != boundary + self.branches.len()
        {
            return Err(invalid("native harmonic response shape disagrees with the network"));
        }
        let state_norm = response
            .state
            .iter()
            .fold(0.0, |total, value| det::hypot(total, value.abs()));
        if !state_norm.is_finite() {
            return Err(invalid("harmonic network state norm is not representable"));
        }
        let branch_flows = &response.output[boundary..];
        let mut node_flow = vec![C64::ZERO; self.topology.node_count];
        let mut node_scale = vec![0.0; self.topology.node_count];
        for (edge, &[from, to]) in self.topology.branch_nodes.iter().enumerate() {
            poll(cx)?;
            let mapped = complex_dot(&self.flow_rows[edge], &response.state)?;
            let uncancelled = self.flow_norms[edge] * state_norm;
            if !uncancelled.is_finite() {
                return Err(invalid("harmonic branch-flow scale is not representable"));
            }
            close(
                branch_flows[edge] - mapped,
                branch_flows[edge].abs() + mapped.abs() + uncancelled,
                tolerance,
                "head-port output disagrees with the original branch flow map",
            )?;
            node_flow[from] = node_flow[from] + mapped;
            node_flow[to] = node_flow[to] - mapped;
            let scale = mapped.abs() + uncancelled;
            node_scale[from] += scale;
            node_scale[to] += scale;
        }
        let mut component_flow = vec![C64::ZERO; boundary];
        let mut component_scale = vec![0.0; boundary];
        for node in 0..node_flow.len() {
            poll(cx)?;
            if let Some(port) = self.basis.boundary_slot[node] {
                close(
                    response.output[port] - node_flow[node],
                    response.output[port].abs() + node_scale[node],
                    tolerance,
                    "reservoir flow disagrees with original branch continuity",
                )?;
                let root = self.basis.pressure_reference[port];
                component_flow[root] = component_flow[root] + response.output[port];
                component_scale[root] += response.output[port].abs() + node_scale[node];
            } else {
                close(
                    node_flow[node],
                    node_scale[node],
                    tolerance,
                    "harmonic trajectory violates internal volume continuity",
                )?;
            }
        }
        for (&flow, &scale) in component_flow.iter().zip(&component_scale) {
            close(
                flow,
                scale,
                tolerance,
                "harmonic component reservoir flows do not conserve volume",
            )?;
        }
        Ok(())
    }
}

fn complex_dot(row: &[f64], values: &[C64]) -> Result<C64, HarmonicError> {
    if row.len() != values.len() {
        return Err(invalid("branch flow observation row has invalid shape"));
    }
    let mut real = 0.0;
    let mut imaginary = 0.0;
    for (&coefficient, &value) in row.iter().zip(values) {
        real += coefficient * value.re;
        imaginary += coefficient * value.im;
    }
    let result = C64::new(real, imaginary);
    if !finite(result) {
        return Err(invalid("branch flow observation is not representable"));
    }
    Ok(result)
}

fn close(
    residual: C64,
    scale: f64,
    tolerance: f64,
    what: &'static str,
) -> Result<(), HarmonicError> {
    if !finite(residual)
        || !scale.is_finite()
        || residual.abs() > tolerance * scale.max(1.0e-300)
    {
        return Err(invalid(what));
    }
    Ok(())
}

fn poll(cx: &Cx<'_>) -> Result<(), HarmonicError> {
    cx.checkpoint().map_err(|_| HarmonicError::Cancelled)
}
fn invalid(what: &'static str) -> HarmonicError {
    HarmonicError::Invalid { what }
}
fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::inertial::laminar::LaminarPipeSpec;
    use crate::graph::tests::with_cx;
    use fs_material::fluid_transport::FluidTransportState;
    use fs_phs::ZkDiffusiveBasis;
    use fs_qty::Length;

    fn pipe(cx: &Cx<'_>, length: f64, viscosity: f64) -> LaminarPipe<FluidTransportState> {
        LaminarPipe::try_new(
            cx,
            FluidTransportState::try_new(1000.0, 4200.0, 0.6, viscosity).unwrap(),
            LaminarPipeSpec {
                length: Length::new(length),
                radius: Length::new(0.001),
                maximum_reynolds: 1000.0,
                band_rad_s: [0.01, 900.0],
                maximum_relative_wall_error: 1.0e-4,
            },
            &ZkDiffusiveBasis::new(16).unwrap(),
        )
        .unwrap()
    }

    fn topology(nodes: usize, edges: &[[usize; 2]], ports: &[usize]) -> LaminarNetworkTopology {
        LaminarNetworkTopology {
            node_count: nodes,
            branch_nodes: edges.to_vec(),
            pressure_nodes: ports.to_vec(),
        }
    }

    fn budget() -> HarmonicBudget {
        HarmonicBudget {
            max_dimension: 128,
            max_dense_bytes: 64 << 20,
            max_dense_work: 200_000_000,
            relative_residual_tolerance: 1.0e-10,
        }
    }

    fn near(actual: C64, expected: C64) {
        let scale = actual.abs() + expected.abs();
        assert!((actual - expected).abs() <= 3.0e-9 * scale.max(1.0e-300));
    }

    #[test]
    fn single_and_series_networks_retain_original_complex_impedance() {
        with_cx(|cx| {
            let omega = 17.0;
            let first = pipe(cx, 0.15, 0.001);
            let z0 = first.impedance_at(omega).unwrap().conj();
            let single = LaminarNetwork::try_new(
                cx,
                topology(2, &[[0, 1]], &[0, 1]),
                vec![first],
            )
            .unwrap();
            let pressure = C64::new(80.0, -25.0);
            let response = single
                .harmonic_response(cx, omega, &[pressure, C64::ZERO], budget())
                .unwrap();
            near(response.output[2], pressure / z0);
            near(response.output[0], response.output[2]);
            near(response.output[1], response.output[2].scale(-1.0));

            let a = pipe(cx, 0.15, 0.001);
            let b = pipe(cx, 0.25, 0.003);
            let impedance =
                a.impedance_at(omega).unwrap().conj() + b.impedance_at(omega).unwrap().conj();
            let series = LaminarNetwork::try_new(
                cx,
                topology(3, &[[0, 1], [1, 2]], &[0, 2]),
                vec![a, b],
            )
            .unwrap();
            let response = series
                .harmonic_response(cx, omega, &[pressure, C64::ZERO], budget())
                .unwrap();
            near(response.output[2], pressure / impedance);
            near(response.output[3], response.output[2]);
            assert!(
                response.power_residual_w.abs()
                    < 2.0e-9 * response.supplied_power_w.abs().max(1.0e-300)
            );
        });
    }

    #[test]
    fn explicit_head_drives_a_closed_loop_without_reservoir_flow() {
        with_cx(|cx| {
            let omega = 9.0;
            let a = pipe(cx, 0.10, 0.001);
            let b = pipe(cx, 0.20, 0.002);
            let c = pipe(cx, 0.12, 0.0015);
            let impedance = a.impedance_at(omega).unwrap().conj()
                + b.impedance_at(omega).unwrap().conj()
                + c.impedance_at(omega).unwrap().conj();
            let network = LaminarNetwork::try_new(
                cx,
                topology(3, &[[0, 1], [1, 2], [2, 0]], &[0]),
                vec![a, b, c],
            )
            .unwrap();
            let head = C64::new(120.0, 30.0);
            let response = network
                .harmonic_response_with_heads(
                    cx,
                    omega,
                    &[C64::new(1.0e9, -2.0e8)],
                    &[head, C64::ZERO, C64::ZERO],
                    budget(),
                )
                .unwrap();
            assert!(response.output[0].abs() < 1.0e-18);
            for flow in &response.output[1..] {
                near(*flow, head / impedance);
            }
            assert!(response.supplied_power_w > 0.0);
        });
    }

    #[test]
    fn malformed_ports_and_out_of_band_frequency_refuse() {
        with_cx(|cx| {
            let network = LaminarNetwork::try_new(
                cx,
                topology(2, &[[0, 1]], &[0, 1]),
                vec![pipe(cx, 0.2, 0.001)],
            )
            .unwrap();
            assert!(
                network
                    .harmonic_response(cx, 10.0, &[C64::ZERO], budget())
                    .is_err()
            );
            assert!(
                network
                    .harmonic_response_with_heads(
                        cx,
                        10.0,
                        &[C64::ZERO; 2],
                        &[C64::ZERO; 2],
                        budget(),
                    )
                    .is_err()
            );
            assert!(
                network
                    .harmonic_response(cx, 1000.0, &[C64::ONE, C64::ZERO], budget())
                    .is_err()
            );
        });
    }
}
