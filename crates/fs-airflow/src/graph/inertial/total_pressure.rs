//! Opt-in incompressible networks with ideal reversible total-pressure nodes.
//!
//! Each uniform-area branch retains the original native impulse L Q and its
//! kinetic storage. At a node, P=p+rho*(Q/A)^2/2 is common, while static pressure
//! depends on the incident branch. Equal dynamic heads at a pipe's two ends
//! cancel in its impulse equation, so the existing descriptor, original drag,
//! gravity, nonlinear residual and acceptance gates remain the numerical owner.
//!
//! This is an authored reversible junction/reducer approximation. It supplies
//! neither mixing loss nor compressible/acoustic storage, cavitation, momentum
//! profiles, swirl, or a moving-pipe law. Fixed total pressures are reservoir
//! stagnation pressures or explicitly authored total-head stream cuts; a
//! finite-speed STATIC-pressure boundary cannot be silently substituted.

use super::loads::{HydraulicAssemblyLoads, JunctionLoad, PipeLoad, PipeLoadGeometry};
use super::{
    FixedPressure, InertialBranch, InertialGraph, InertialGraphControls, InertialGraphError,
    InertialGraphReport, InertialGraphTrial, LossGraph, MAX_INERTIAL_GRAPH_STATES, PressureChart,
    Sum, checkpoint, finite, nonnegative,
};
use fs_exec::Cx;
use fs_qty::{Density, Pressure};

mod transport;
pub use transport::{BranchStreamWork, NodeStreamWork};

/// Prescribed P=p+rho*v^2/2, in the same pressure datum as the source graph.
/// This type is intentionally distinct from the static FixedPressure boundary.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FixedTotalPressure {
    pub node: usize,
    pub total_pressure: Pressure,
}

/// The sole retained physical state is the original branch-momentum owner.
/// Its private pressure chart cannot escape into the static thermal adapters.
#[derive(Clone, Debug)]
pub struct TotalPressureGraph {
    native: InertialGraph,
}

#[derive(Clone)]
pub struct TotalPressureGraphTrial {
    native: InertialGraphTrial,
    report: TotalPressureGraphReport,
}

/// Exact affine stream work accompanying the same native midpoint trajectory.
/// The native report is private so total multipliers cannot be consumed as
/// static nodal pressures by an unrelated source/thermal adapter.
#[derive(Clone, Debug)]
pub struct TotalPressureGraphReport {
    native: InertialGraphReport,
    /// Original branch order, with from/to port orientation.
    pub branch_streams: Vec<BranchStreamWork>,
    /// Original node order, retaining every local pressure/kinetic balance.
    pub node_streams: Vec<NodeStreamWork>,
    /// Positive work/kinetic transport ENTERS the graph at fixed pressure cuts.
    pub boundary_static_pressure_work_j: f64,
    pub boundary_kinetic_energy_inflow_j: f64,
    /// Full physical gravity, rather than only its anomalous reduced component.
    pub gravity_work_j: f64,
    /// Delta K + loss - drive - static pressure work - kinetic inflow - gravity.
    pub mechanical_balance_residual_j: f64,
    pub maximum_internal_node_work_residual_j: f64,
    pub maximum_stream_identity_residual_j: f64,
}

impl TotalPressureGraphReport {
    pub fn duration_s(&self) -> f64 {
        self.native.duration_s
    }
    pub fn average_total_pressures_pa(&self) -> &[f64] {
        &self.native.average_pressures_pa
    }
    pub fn average_volume_flows_m3_s(&self) -> &[f64] {
        &self.native.average_volume_flows_m3_s
    }
    pub fn initial_flows_m3_s(&self) -> &[f64] {
        &self.native.initial_flows_m3_s
    }
    pub fn final_flows_m3_s(&self) -> &[f64] {
        &self.native.final_flows_m3_s
    }
    pub fn average_flow_residual_m3_s(&self) -> &[f64] {
        &self.native.average_flow_residual_m3_s
    }
    pub fn original_impulse_residuals_pa_s(&self) -> &[f64] {
        &self.native.original_impulse_residuals_pa_s
    }
    pub fn maximum_impulse_residual_pa_s(&self) -> f64 {
        self.native.maximum_impulse_residual_pa_s
    }
    pub fn hydraulic_loss_j(&self) -> f64 {
        self.native.hydraulic_loss_j
    }
    pub fn driving_work_j(&self) -> f64 {
        self.native.driving_work_j
    }
    pub fn neglected_energy_bound_j(&self) -> f64 {
        self.native.neglected_energy_bound_j
    }
    /// Unmodified native integrator result; algebraic tail entries encode total
    /// pressure in this owner. They are not additional stored fluid history.
    pub fn momentum(&self) -> &fs_phs::StepRecord {
        &self.native.momentum
    }
}

fn boundaries(values: &[FixedTotalPressure]) -> Result<Vec<FixedPressure>, InertialGraphError> {
    if values.len() > MAX_INERTIAL_GRAPH_STATES {
        return Err(InertialGraphError::Invalid(
            "total-pressure boundary count exceeds native state cap",
        ));
    }
    Ok(values
        .iter()
        .map(|value| FixedPressure {
            node: value.node,
            pressure: value.total_pressure,
        })
        .collect())
}

impl TotalPressureGraph {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        cx: &Cx<'_>,
        graph: &LossGraph,
        reference_density: Density,
        branches: &[InertialBranch],
        pressure_boundaries: &[FixedTotalPressure],
        initial_flows_m3_s: &[f64],
        accepted_time_s: f64,
        maximum_states: usize,
    ) -> Result<Self, InertialGraphError> {
        Ok(Self {
            native: InertialGraph::try_new_in_chart(
                cx,
                graph,
                reference_density,
                branches,
                &boundaries(pressure_boundaries)?,
                initial_flows_m3_s,
                accepted_time_s,
                maximum_states,
                PressureChart::Total,
            )?,
        })
    }

    pub fn node_count(&self) -> usize {
        self.native.node_count()
    }
    pub fn branch_count(&self) -> usize {
        self.native.branch_count()
    }
    pub fn accepted_time_s(&self) -> f64 {
        self.native.accepted_time_s()
    }
    pub fn reference_density(&self) -> Density {
        self.native.reference_density()
    }
    pub fn inertances_pa_s2_m3(&self) -> &[f64] {
        self.native.inertances_pa_s2_m3()
    }
    pub fn volume_flows_m3_s(&self) -> Vec<f64> {
        self.native.volume_flows_m3_s()
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        self.native.kinetic_energy_j()
    }

    /// The original held-coefficient proposal, with explicit total-pressure
    /// ports. Source pressure admission uses the full STATIC-pressure extrema
    /// of the affine flow, including its interior maximum during reversal.
    /// Nonlinear drag still requires original candidate/actual flow closure.
    #[allow(clippy::too_many_arguments)]
    pub fn propose(
        &self,
        cx: &Cx<'_>,
        graph: &LossGraph,
        pressure_boundaries: &[FixedTotalPressure],
        driving_pressure_pa: &[f64],
        candidate_average_flows_m3_s: &[f64],
        target_s: f64,
        controls: InertialGraphControls,
    ) -> Result<TotalPressureGraphTrial, InertialGraphError> {
        TotalPressureGraphTrial::bind(
            cx,
            self.native.propose(
                cx,
                graph,
                &boundaries(pressure_boundaries)?,
                driving_pressure_pa,
                candidate_average_flows_m3_s,
                target_s,
                controls,
            )?,
        )
    }
}

impl TotalPressureGraphTrial {
    fn bind(cx: &Cx<'_>, native: InertialGraphTrial) -> Result<Self, InertialGraphError> {
        let report = stream_report(cx, &native)?;
        Ok(Self { native, report })
    }

    pub fn report(&self) -> &TotalPressureGraphReport {
        &self.report
    }

    /// Re-admit original source coefficients on the SAME native trajectory.
    pub fn recheck(
        self,
        cx: &Cx<'_>,
        graph: &LossGraph,
        pressure_boundaries: &[FixedTotalPressure],
        driving_pressure_pa: &[f64],
        controls: InertialGraphControls,
    ) -> Result<Self, InertialGraphError> {
        Self::bind(
            cx,
            self.native.recheck(
                cx,
                graph,
                &boundaries(pressure_boundaries)?,
                driving_pressure_pa,
                controls,
            )?,
        )
    }

    /// Static pressure traction on the actual fixed path. Mean Q^2 sets this
    /// load; mean Q^3 independently sets stream work in report(). No second
    /// pressure-loss force or thermal credit is added.
    pub fn pipe_load(
        &self,
        cx: &Cx<'_>,
        branch: usize,
        geometry: &PipeLoadGeometry,
        exterior_pressure: Pressure,
        about_m: [f64; 3],
    ) -> Result<PipeLoad, InertialGraphError> {
        self.native
            .pipe_load(cx, branch, geometry, exterior_pressure, about_m)
    }

    /// All incident port tractions use their own physical static pressure.
    /// report().node_streams[node] supplies the corresponding pressure and
    /// kinetic work; no common-static-pressure assumption is made here.
    pub fn junction_load(
        &self,
        cx: &Cx<'_>,
        node: usize,
        geometries: &[PipeLoadGeometry],
        exterior_pressure: Pressure,
        about_m: [f64; 3],
    ) -> Result<JunctionLoad, InertialGraphError> {
        self.native
            .junction_load(cx, node, geometries, exterior_pressure, about_m)
    }

    /// Original component and whole-assembly momentum balances, evaluated with
    /// static tractions from this total-pressure chart. Internal interfaces
    /// still cancel analytically, retaining actual fluid storage and weight.
    pub fn assembly_loads(
        &self,
        cx: &Cx<'_>,
        geometries: &[PipeLoadGeometry],
        enclosed_pressure_nodes: &[usize],
        exterior_pressure: Pressure,
        about_m: [f64; 3],
    ) -> Result<HydraulicAssemblyLoads, InertialGraphError> {
        self.native.assembly_loads(
            cx,
            geometries,
            enclosed_pressure_nodes,
            exterior_pressure,
            about_m,
        )
    }

    pub fn accept(
        mut self,
        cx: &Cx<'_>,
        impulse_tolerance_pa_s: f64,
        mechanical_tolerance_j: f64,
    ) -> Result<(TotalPressureGraph, TotalPressureGraphReport), InertialGraphError> {
        checkpoint(cx)?;
        self.admit_work(mechanical_tolerance_j)?;
        let (native, report) =
            self.native
                .accept(cx, impulse_tolerance_pa_s, mechanical_tolerance_j)?;
        self.report.native = report;
        Ok((TotalPressureGraph { native }, self.report))
    }

    #[allow(clippy::too_many_arguments)]
    pub fn accept_with_assembly_loads(
        self,
        cx: &Cx<'_>,
        geometries: &[PipeLoadGeometry],
        enclosed_pressure_nodes: &[usize],
        exterior_pressure: Pressure,
        about_m: [f64; 3],
        impulse_tolerance_pa_s: f64,
        mechanical_tolerance_j: f64,
    ) -> Result<
        (
            TotalPressureGraph,
            TotalPressureGraphReport,
            HydraulicAssemblyLoads,
        ),
        InertialGraphError,
    > {
        let loads = self.assembly_loads(
            cx,
            geometries,
            enclosed_pressure_nodes,
            exterior_pressure,
            about_m,
        )?;
        self.admit_work(mechanical_tolerance_j)?;
        for junction in &loads.junctions {
            checkpoint(cx)?;
            let energy = self.report.node_streams[junction.node].zero_storage_balance_residual_j;
            if energy.abs() > mechanical_tolerance_j {
                return Err(InertialGraphError::MomentumResidual {
                    impulse_pa_s: self.report.maximum_impulse_residual_pa_s(),
                    energy_j: energy,
                });
            }
        }
        let (next, report) = self.accept(cx, impulse_tolerance_pa_s, mechanical_tolerance_j)?;
        Ok((next, report, loads))
    }

    fn admit_work(&self, tolerance: f64) -> Result<(), InertialGraphError> {
        if !nonnegative(tolerance) {
            return Err(InertialGraphError::Invalid(
                "finite nonnegative work tolerance required",
            ));
        }
        let residual = self
            .report
            .mechanical_balance_residual_j
            .abs()
            .max(self.report.maximum_internal_node_work_residual_j)
            .max(self.report.maximum_stream_identity_residual_j);
        if residual > tolerance {
            return Err(InertialGraphError::MomentumResidual {
                impulse_pa_s: self.report.maximum_impulse_residual_pa_s(),
                energy_j: residual,
            });
        }
        Ok(())
    }
}

impl PressureChart {
    pub(super) fn static_mean_pressures(
        self,
        pressures: [f64; 2],
        flows: [f64; 2],
        density: f64,
        area: f64,
    ) -> Result<[f64; 2], InertialGraphError> {
        if self == Self::Static {
            return Ok(pressures);
        }
        let head = transport::mean_dynamic_pressure(flows, density, area)?;
        Ok([
            finite(pressures[0] - head, "static inlet pressure")?,
            finite(pressures[1] - head, "static outlet pressure")?,
        ])
    }

    pub(super) fn source_pressure_bounds(
        self,
        pressures: [f64; 2],
        flows: [f64; 2],
        density: f64,
        area: f64,
    ) -> Result<[[f64; 2]; 2], InertialGraphError> {
        if self == Self::Static {
            return Ok([pressures; 2]);
        }
        Ok(transport::static_pressure_bounds(
            pressures, flows, density, area,
        )?)
    }
}

fn stream_report(
    cx: &Cx<'_>,
    trial: &InertialGraphTrial,
) -> Result<TotalPressureGraphReport, InertialGraphError> {
    checkpoint(cx)?;
    let c = &trial.initial.compiled;
    let r = &trial.report;
    let mut branch_streams = Vec::with_capacity(c.branches.len());
    let mut pressure = vec![Sum::default(); c.free.len()];
    let mut kinetic = vec![Sum::default(); c.free.len()];
    let mut total = vec![Sum::default(); c.free.len()];
    let mut gravity = Sum::default();
    for (branch, edge) in c.topology.branches().iter().enumerate() {
        checkpoint(cx)?;
        let work = transport::branch_work(
            r.duration_s,
            c.reference_density.value(),
            c.branches[branch].flow_area_m2,
            [
                r.average_pressures_pa[edge.from],
                r.average_pressures_pa[edge.to],
            ],
            [r.initial_flows_m3_s[branch], r.final_flows_m3_s[branch]],
        )?;
        for (port, node) in [edge.from, edge.to].into_iter().enumerate() {
            pressure[node].add(work.static_pressure_work_into_branch_j[port])?;
            kinetic[node].add(work.kinetic_energy_into_branch_j[port])?;
            total[node].add(work.total_pressure_work_into_branch_j[port])?;
        }
        gravity.add(finite(
            r.duration_s
                * r.gravity_pressure_heads_pa[branch]
                * r.average_volume_flows_m3_s[branch],
            "full gravity work",
        )?)?;
        branch_streams.push(work);
    }
    let mut boundary_pressure = Sum::default();
    let mut boundary_kinetic = Sum::default();
    let mut node_streams = Vec::with_capacity(c.free.len());
    let mut internal_residual = 0.0_f64;
    let mut stream_residual = 0.0_f64;
    for node in 0..c.free.len() {
        checkpoint(cx)?;
        let p = pressure[node].value()?;
        let k = kinetic[node].value()?;
        let t = total[node].value()?;
        let balance = finite(p + k, "node static-pressure and kinetic work")?;
        let identity = finite(balance - t, "node total-pressure identity")?;
        stream_residual = stream_residual.max(identity.abs());
        if c.free[node].is_some() {
            internal_residual = internal_residual.max(balance.abs());
        } else {
            boundary_pressure.add(p)?;
            boundary_kinetic.add(k)?;
        }
        node_streams.push(NodeStreamWork {
            node,
            static_pressure_work_out_j: p,
            kinetic_energy_out_j: k,
            total_pressure_work_out_j: t,
            pressure_kinetic_identity_residual_j: identity,
            zero_storage_balance_residual_j: balance,
        });
    }
    let boundary_static_pressure_work_j = boundary_pressure.value()?;
    let boundary_kinetic_energy_inflow_j = boundary_kinetic.value()?;
    let gravity_work_j = gravity.value()?;
    let mechanical_balance_residual_j = finite(
        r.momentum.delta_h + r.hydraulic_loss_j
            - r.driving_work_j
            - boundary_static_pressure_work_j
            - boundary_kinetic_energy_inflow_j
            - gravity_work_j,
        "full pressure and kinetic stream balance",
    )?;
    checkpoint(cx)?;
    Ok(TotalPressureGraphReport {
        native: r.clone(),
        branch_streams,
        node_streams,
        boundary_static_pressure_work_j,
        boundary_kinetic_energy_inflow_j,
        gravity_work_j,
        mechanical_balance_residual_j,
        maximum_internal_node_work_residual_j: internal_residual,
        maximum_stream_identity_residual_j: stream_residual,
    })
}

#[cfg(test)]
mod tests;
