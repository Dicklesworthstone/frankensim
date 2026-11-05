//! Original per-pipe applicability, storage, heat and continuity audits.

use super::*;
use fs_math::det;
use fs_phs::harmonic::HarmonicFactor;

impl<M: FluidTransport> LaminarNetwork<M> {
    pub(super) fn response_with_factor(
        &self,
        cx: &Cx<'_>,
        factor: &HarmonicFactor<'_>,
        input: &[C64],
        tolerance: f64,
    ) -> Result<HarmonicResponse, HarmonicError> {
        let response = factor.response(input, || cx.checkpoint().is_err())?;
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
        let n = self.system.state_dim();
        if response.state.len() != n
            || response.effort.len() != n
            || response.output.len() != boundary + self.branches.len()
        {
            return Err(invalid(
                "native harmonic response shape disagrees with the network",
            ));
        }
        let state_norm = response
            .state
            .iter()
            .fold(0.0, |total, value| det::hypot(total, value.abs()));
        if !state_norm.is_finite() {
            return Err(invalid(
                "harmonic network state norm is not representable",
            ));
        }
        let branch_flows = &response.output[boundary..];
        let mut node_flow = vec![C64::ZERO; self.topology.node_count];
        let mut node_scale = vec![0.0; self.topology.node_count];
        let mut local_heat_w = 0.0;
        let mut local_storage_j = 0.0;
        let mean_count = self.basis.columns;
        for (edge, &[from, to]) in self.topology.branch_nodes.iter().enumerate() {
            poll(cx)?;
            let mapped = complex_dot(&self.flow_rows[edge], &response.state)?;
            let uncancelled = self.flow_norms[edge] * state_norm;
            if !uncancelled.is_finite() {
                return Err(invalid(
                    "harmonic branch-flow scale is not representable",
                ));
            }
            close(
                branch_flows[edge] - mapped,
                branch_flows[edge].abs() + mapped.abs() + uncancelled,
                tolerance,
                "head-port output disagrees with the original branch flow map",
            )?;
            let reynolds = self.branches[edge]
                .admit_flow(branch_flows[edge].abs())
                .map_err(|_| {
                    invalid("harmonic branch peak flow exceeds original laminar admission")
                })?;
            if !reynolds.is_finite() {
                return Err(invalid(
                    "harmonic branch Reynolds amplitude is not representable",
                ));
            }

            let mean = complex_dot(
                &self.basis.energy_rows
                    [edge * mean_count..(edge + 1) * mean_count],
                &response.state[..mean_count],
            )?;
            let memory_start = self.offsets[edge];
            let memory_end = self.offsets[edge + 1];
            let memory_energy: f64 = response.state[memory_start..memory_end]
                .iter()
                .map(|value| 0.25 * value.norm_sq())
                .sum();
            let mean_energy = 0.25 * mean.norm_sq();
            if !mean_energy.is_finite() || !memory_energy.is_finite() {
                return Err(invalid(
                    "harmonic branch kinetic storage is not representable",
                ));
            }
            local_storage_j += mean_energy + memory_energy;

            let source = &self.source_rows[edge];
            if source.len() % n != 0 {
                return Err(invalid(
                    "harmonic branch source-square layout is inconsistent",
                ));
            }
            for row in source.chunks_exact(n) {
                let amplitude = complex_dot(row, &response.effort)?;
                local_heat_w += 0.5 * amplitude.norm_sq();
            }
            if !local_heat_w.is_finite() || !local_storage_j.is_finite() {
                return Err(invalid(
                    "harmonic local heat or storage is not representable",
                ));
            }

            node_flow[from] = node_flow[from] + mapped;
            node_flow[to] = node_flow[to] - mapped;
            let scale = mapped.abs() + uncancelled;
            node_scale[from] += scale;
            node_scale[to] += scale;
        }

        scalar_close(
            local_heat_w - response.dissipated_power_w,
            local_heat_w.abs() + response.dissipated_power_w.abs(),
            tolerance + ROUNDOFF * n as f64,
            "original per-pipe viscous heat disagrees with native network loss",
        )?;
        scalar_close(
            local_storage_j - response.mean_stored_energy_j,
            local_storage_j.abs() + response.mean_stored_energy_j.abs(),
            tolerance + ROUNDOFF * n as f64,
            "original per-pipe kinetic storage disagrees with native network storage",
        )?;

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

fn scalar_close(
    residual: f64,
    scale: f64,
    tolerance: f64,
    what: &'static str,
) -> Result<(), HarmonicError> {
    if !residual.is_finite()
        || !scale.is_finite()
        || residual.abs() > tolerance * scale.max(1.0e-300)
    {
        return Err(invalid(what));
    }
    Ok(())
}
