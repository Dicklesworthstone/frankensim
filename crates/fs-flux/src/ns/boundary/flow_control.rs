//! Integral port flow with an unknown, power-conjugate pressure.

use fs_la::factor::{Lu, lu};
use fs_math::dd::Dd;

use super::{StokesAssembly, StokesBoundaryError, StokesSolution};
use crate::ns::FluxSolution;

#[path = "flow_control/resistance.rs"]
mod resistance;
pub use resistance::{StokesResistancePort, StokesResistanceSolution};
#[path = "flow_control/windkessel.rs"]
mod windkessel;
pub(crate) use windkessel::solve_windkessel_assembly;
pub use windkessel::{
    StokesWindkesselBalance, StokesWindkesselHistory, StokesWindkesselPort,
    StokesWindkesselSolution,
};

/// Resolved flow and the uniform pressure increments required at controlled ports.
pub struct StokesFlowControlSolution<S = StokesSolution> {
    /// Accepted linear or rheological flow, including its boundary reactions.
    pub flow: S,
    /// Pressure added to each port's original traction as -p*n [m²/s²].
    /// Multiply by the SAME reference density as the assembled equations for Pa.
    pub pressure_offsets_m2_s2: Vec<f64>,
    /// Actual signed outward integral of velocity at each port [m²/s].
    /// This is flow per unit out-of-plane depth; multiply by depth for m³/s.
    pub outward_volume_rates_m2_s: Vec<f64>,
    /// Largest absolute difference between the solved and requested port flows.
    pub maximum_flow_residual_m2_s: f64,
}

/// Cached native responses for repeated flow targets on one frozen linear system.
///
/// The original operator, body force, wall velocities, traction shapes, density
/// metric and (for transient assemblies) history stay fixed. Reprepare when any
/// of them changes. Preparation performs one original solve plus one pressure
/// response solve per port and stores one complete response vector per port.
/// Each subsequent solve uses the small cached pressure factorization and
/// superposes those fields; it does not reassemble, query a callback or solve
/// another full flow system. Both original equation residuals and each integral
/// flow constraint are checked before a result is returned.
pub struct PreparedStokesFlowControl<'a, 'm> {
    assembly: StokesAssembly<'a, 'm>,
    ports: Vec<Vec<usize>>,
    baseline: Vec<f64>,
    baseline_flows: Vec<f64>,
    responses: Vec<Vec<f64>>,
    mobility: Vec<f64>,
    row_scales: Vec<f64>,
    pressure_factor: Lu,
}

impl<'a, 'm> StokesAssembly<'a, 'm> {
    /// Control the TOTAL outward flow of each supplied group of traction edges.
    ///
    /// Each nonempty port must lie in one connected fluid component. Ports are
    /// disjoint and may contain only original natural-traction edges. Every
    /// controlled component must retain at least one uncontrolled traction edge:
    /// its original load supplies the pressure datum. All-velocity components
    /// retain their independent gauges and are otherwise untouched.
    ///
    /// The added load is one spatially uniform pressure increment per port,
    /// leaving its original tangential traction and pressure shape intact.
    /// No inlet velocity profile is prescribed. The exact BDM mean-normal
    /// moments supply BOTH the integral constraint and its conjugate pressure
    /// load, including curved/piecewise-straight groups of boundary edges.
    ///
    /// Works with the existing variable-viscosity and fixed-history linear
    /// Stokes assemblies. A frozen Oseen/rheological operator remains frozen;
    /// this entry does not add a nonlinear iteration.
    ///
    /// # Errors
    /// Refuses invalid/repeated/non-traction edges, missing pressure datums,
    /// singular port response or failure of any original checked linear solve.
    pub fn prepare_flow_control(
        self,
        ports: &[&[usize]],
    ) -> Result<PreparedStokesFlowControl<'a, 'm>, StokesBoundaryError> {
        if ports.is_empty() || ports.iter().any(|port| port.is_empty()) {
            return Err(invalid("at least one nonempty controlled port required"));
        }
        let mesh = self.system.mesh;
        let mut selected = vec![false; mesh.edges.len()];
        let mut controlled_components = vec![false; self.has_traction.len()];
        for port in ports {
            let mut component = None;
            for &edge in *port {
                if edge >= mesh.edges.len() {
                    return Err(StokesBoundaryError::EdgeOutOfRange { edge });
                }
                if !mesh.is_boundary(edge) {
                    return Err(StokesBoundaryError::InteriorEdge { edge });
                }
                if selected[edge] {
                    return Err(StokesBoundaryError::DuplicateEdge { edge });
                }
                if !self.reactions.is_traction(edge) {
                    return Err(invalid(
                        "flow-controlled edges must retain natural traction",
                    ));
                }
                let current = self.system.pressure_component[mesh.edges[edge].tris.0];
                if component.is_some_and(|previous| previous != current) {
                    return Err(invalid(
                        "each controlled port must belong to one fluid component",
                    ));
                }
                component = Some(current);
                controlled_components[current] = true;
                selected[edge] = true;
            }
        }
        let mut pressure_datum = vec![false; self.has_traction.len()];
        for (edge, geometry) in mesh.edges.iter().enumerate() {
            if !selected[edge] && self.reactions.is_traction(edge) {
                pressure_datum[self.system.pressure_component[geometry.tris.0]] = true;
            }
        }
        if controlled_components
            .iter()
            .zip(&pressure_datum)
            .any(|(&controlled, &datum)| controlled && !datum)
        {
            return Err(invalid(
                "each controlled component needs an uncontrolled traction pressure datum",
            ));
        }
        let count = ports.len();
        let entries = count
            .checked_mul(count)
            .ok_or_else(|| invalid("port-response dimensions overflow"))?;
        let baseline = checked_solve(&self, &self.rhs)?.x;
        let baseline_flows = ports
            .iter()
            .map(|port| port_flow(&self, port, &baseline))
            .collect();
        let mut responses = Vec::with_capacity(count);
        for port in ports {
            let mut load = vec![0.0; self.system.n];
            for &edge in *port {
                // Mean-normal BDM moment: integral(v.n)=length*v[2*edge].
                // The same functional makes pressure and volume flow conjugate.
                load[2 * edge] = -mesh.edges[edge].len;
            }
            responses.push(checked_solve(&self, &load)?.x);
        }
        let mut mobility = vec![0.0; entries];
        let mut scaled = vec![0.0; entries];
        let mut row_scales = Vec::with_capacity(count);
        for (i, port) in ports.iter().enumerate() {
            for (j, response) in responses.iter().enumerate() {
                mobility[i * count + j] = port_flow(&self, port, response);
            }
            let scale = mobility[i * count..(i + 1) * count]
                .iter()
                .fold(0.0_f64, |largest, value| largest.max(value.abs()));
            if !scale.is_finite() || scale == 0.0 {
                return Err(StokesBoundaryError::SingularSystem);
            }
            row_scales.push(scale);
            for j in 0..count {
                scaled[i * count + j] = mobility[i * count + j] / scale;
            }
        }
        if mobility.iter().any(|value| !value.is_finite()) {
            return Err(invalid("unrepresentable hydraulic port response"));
        }
        let pressure_factor =
            lu(&scaled, count).map_err(|_| StokesBoundaryError::SingularSystem)?;
        Ok(PreparedStokesFlowControl {
            assembly: self,
            ports: ports.iter().map(|port| port.to_vec()).collect(),
            baseline,
            baseline_flows,
            responses,
            mobility,
            row_scales,
            pressure_factor,
        })
    }
}

impl PreparedStokesFlowControl<'_, '_> {
    /// Solve for pressures producing signed OUTWARD flows per unit depth [m²/s].
    ///
    /// Positive target means outflow; inflow is negative. The flow allowance is
    /// `absolute_tolerance_m2_s + native_relative_tolerance*abs(target)`.
    /// The caller supplies a finite nonnegative absolute allowance, including
    /// for zero-flow targets. No pressure sign or flow direction is clipped.
    ///
    /// # Errors
    /// Refuses invalid targets/tolerance, unrepresentable pressure or field,
    /// failed original equation residuals, or a missed integral flow target.
    pub fn solve(
        &self,
        outward_volume_rates_m2_s: &[f64],
        absolute_tolerance_m2_s: f64,
    ) -> Result<StokesFlowControlSolution, StokesBoundaryError> {
        let count = self.ports.len();
        admit_flow_targets(count, outward_volume_rates_m2_s, absolute_tolerance_m2_s)?;
        let mut pressures: Vec<_> = outward_volume_rates_m2_s
            .iter()
            .zip(&self.baseline_flows)
            .zip(&self.row_scales)
            .map(|((&target, &base), &scale)| (target - base) / scale)
            .collect();
        if pressures.iter().any(|value| !value.is_finite()) {
            return Err(invalid("unrepresentable controlled-flow load"));
        }
        self.pressure_factor.solve(&mut pressures);
        // Refine only the small pressure system; native responses stay cached.
        for _ in 0..2 {
            let mut correction: Vec<_> = outward_volume_rates_m2_s
                .iter()
                .enumerate()
                .map(|(i, &target)| {
                    let mut residual = Dd::from_f64(target) - Dd::from_f64(self.baseline_flows[i]);
                    for (j, &pressure) in pressures.iter().enumerate() {
                        residual = residual
                            - Dd::from_f64(self.mobility[i * count + j]) * Dd::from_f64(pressure);
                    }
                    residual.to_f64() / self.row_scales[i]
                })
                .collect();
            self.pressure_factor.solve(&mut correction);
            for (pressure, delta) in pressures.iter_mut().zip(correction) {
                *pressure += delta;
            }
        }
        if pressures.iter().any(|value| !value.is_finite()) {
            return Err(invalid("unrepresentable controlled port pressure"));
        }
        let flow = self.report_pressures(&pressures)?;
        let (flows, maximum_flow_residual_m2_s) = checked_port_flows(
            &self.assembly,
            &self.ports,
            outward_volume_rates_m2_s,
            absolute_tolerance_m2_s,
            &flow.field.x,
        )?;
        Ok(StokesFlowControlSolution {
            flow,
            pressure_offsets_m2_s2: pressures,
            outward_volume_rates_m2_s: flows,
            maximum_flow_residual_m2_s,
        })
    }

    // Every cached port solve uses the same compensated superposition,
    // original-equation check and pressure-aware boundary reactions.
    fn report_pressures(&self, pressures: &[f64]) -> Result<StokesSolution, StokesBoundaryError> {
        let mut x = self.baseline.clone();
        for (i, value) in x.iter_mut().enumerate() {
            let mut sum = Dd::from_f64(*value);
            for (&pressure, response) in pressures.iter().zip(&self.responses) {
                sum = sum + Dd::from_f64(pressure) * Dd::from_f64(response[i]);
            }
            *value = sum.to_f64();
        }
        let rhs = pressure_rhs(&self.assembly, &self.ports, pressures)?;
        if x.iter().any(|value| !value.is_finite()) {
            return Err(invalid("unrepresentable controlled flow field or load"));
        }
        let mut work = vec![0.0; x.len()];
        let (rel_residual, scaled_rel_residual) =
            crate::ns::sparse::residuals(&self.assembly.matrix, &rhs, &x, &mut work);
        let field = FluxSolution {
            x,
            n_u: self.assembly.system.n_u,
            iters: 0,
            rel_residual,
            scaled_rel_residual,
        };
        check_convergence(&self.assembly, &field)?;
        let mut flow = self.assembly.report_solution(field)?;
        flow.reactions.add_port_pressures(&self.ports, pressures)?;
        Ok(flow)
    }

    /// Kinetic energy per depth divided by the same reference density.
    /// Uses the original retained density metric without source callbacks.
    ///
    /// # Errors
    /// Preserves the original assembly's velocity and energy admission.
    pub fn kinetic_energy(&self, velocity: &[f64]) -> Result<f64, StokesBoundaryError> {
        self.assembly.kinetic_energy(velocity)
    }
}

// Borrowed nonlinear constraints; geometry is admitted by each native
// response preparation, while the same load/flow functionals check endpoints.
#[derive(Clone, Copy)]
pub(super) struct FlowControl<'a> {
    pub(super) ports: &'a [&'a [usize]],
    pub(super) targets: &'a [f64],
    pub(super) absolute_tolerance_m2_s: f64,
}

impl FlowControl<'_> {
    pub(super) fn admit(&self) -> Result<(), StokesBoundaryError> {
        if self.ports.is_empty() || self.ports.iter().any(|port| port.is_empty()) {
            return Err(invalid("at least one nonempty controlled port required"));
        }
        admit_flow_targets(self.ports.len(), self.targets, self.absolute_tolerance_m2_s)
    }
}

fn admit_flow_targets(
    count: usize,
    targets: &[f64],
    absolute_tolerance_m2_s: f64,
) -> Result<(), StokesBoundaryError> {
    if targets.len() != count
        || targets.iter().any(|value| !value.is_finite())
        || !absolute_tolerance_m2_s.is_finite()
        || absolute_tolerance_m2_s < 0.0
    {
        return Err(invalid(
            "finite complete port flows and nonnegative absolute flow tolerance required",
        ));
    }
    Ok(())
}

// The frozen linear solve and the reassembled nonlinear residual use exactly
// the same conjugate pressure load, including its compensated arithmetic.
pub(super) fn pressure_rhs(
    assembly: &StokesAssembly<'_, '_>,
    ports: &[impl AsRef<[usize]>],
    pressures: &[f64],
) -> Result<Vec<f64>, StokesBoundaryError> {
    if pressures.len() != ports.len() || pressures.iter().any(|value| !value.is_finite()) {
        return Err(invalid(
            "finite complete controlled port pressures required",
        ));
    }
    let mut rhs = assembly.rhs.clone();
    for (port, &pressure) in ports.iter().zip(pressures) {
        for &edge in port.as_ref() {
            rhs[2 * edge] = (Dd::from_f64(rhs[2 * edge])
                - Dd::from_f64(pressure) * Dd::from_f64(assembly.system.mesh.edges[edge].len))
            .to_f64();
        }
    }
    if rhs.iter().any(|value| !value.is_finite()) {
        return Err(invalid("unrepresentable controlled flow field or load"));
    }
    Ok(rhs)
}

pub(super) fn checked_port_flows(
    assembly: &StokesAssembly<'_, '_>,
    ports: &[impl AsRef<[usize]>],
    targets: &[f64],
    absolute_tolerance_m2_s: f64,
    velocity: &[f64],
) -> Result<(Vec<f64>, f64), StokesBoundaryError> {
    admit_flow_targets(ports.len(), targets, absolute_tolerance_m2_s)?;
    let flows: Vec<_> = ports
        .iter()
        .map(|port| port_flow(assembly, port.as_ref(), velocity))
        .collect();
    let mut maximum_flow_residual_m2_s = 0.0_f64;
    for (port, (&requested, &actual)) in targets.iter().zip(&flows).enumerate() {
        let tolerance = absolute_tolerance_m2_s + assembly.params.tol * requested.abs();
        let residual = (actual - requested).abs();
        if !actual.is_finite() || !tolerance.is_finite() || residual > tolerance {
            return Err(StokesBoundaryError::FlowRateNotConverged {
                port,
                requested_m2_s: requested,
                actual_m2_s: actual,
                tolerance_m2_s: tolerance,
            });
        }
        maximum_flow_residual_m2_s = maximum_flow_residual_m2_s.max(residual);
    }
    Ok((flows, maximum_flow_residual_m2_s))
}

fn port_flow(assembly: &StokesAssembly<'_, '_>, edges: &[usize], x: &[f64]) -> f64 {
    edges
        .iter()
        .fold(Dd::from_f64(0.0), |sum, &edge| {
            sum + Dd::from_f64(assembly.system.mesh.edges[edge].len) * Dd::from_f64(x[2 * edge])
        })
        .to_f64()
}

fn checked_solve(
    assembly: &StokesAssembly<'_, '_>,
    rhs: &[f64],
) -> Result<FluxSolution, StokesBoundaryError> {
    let field = assembly
        .system
        .solve_linear_result(&assembly.matrix, rhs, assembly.params)?;
    check_convergence(assembly, &field)?;
    Ok(field)
}

fn check_convergence(
    assembly: &StokesAssembly<'_, '_>,
    field: &FluxSolution,
) -> Result<(), StokesBoundaryError> {
    if field.x.iter().any(|value| !value.is_finite()) || !field.converged(assembly.params.tol) {
        return Err(StokesBoundaryError::NotConverged {
            iterations: field.iters,
            relative_residual: field.rel_residual,
            scaled_relative_residual: field.scaled_rel_residual,
        });
    }
    Ok(())
}

fn invalid(reason: &'static str) -> StokesBoundaryError {
    StokesBoundaryError::InvalidInput(reason)
}
