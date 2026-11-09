//! Resolved Carreau momentum and sensible heat in one shared inventory.
//!
//! All pipes, hydraulic storage and returning reservoir fluid are declared to
//! share one uniform temperature and one original caloric storage law. Every
//! trial advances the original radial histories at the trial midpoint
//! temperature. Only their actual integrated viscous heat enters sensible
//! storage; pressure work, kinetic energy and compliance energy remain distinct.
//!
//! The scalar energy equation uses the same bounded midpoint temperature solve
//! as the quasisteady thermal owner. This closes viscosity feedback during
//! startup, coastdown and pressure reversal without fitting a mean-flow inertia.
//! Fixed density/compliance and rapid thermal equilibration are explicit model
//! assumptions. This is not spatial thermal advection, an open stream, thermal
//! expansion, an inferred heat capacity, or a radial temperature-field solve.

use super::{
    Cx, FixedPressure, FlowInjection, GraphError, MAX_THERMAL_HYDRAULIC_SOLVES, Point,
    TabulatedHeatStorage, ThermalCarreauBoundary, ThermalCarreauError, Trial, checkpoint, finite,
    nonnegative, positive, shared_temperature_span, solve_temperature, sum,
};
use crate::graph::rheology::transient::{
    TransientCarreauControls, TransientCarreauError, TransientCarreauGraph, TransientCarreauState,
    TransientCarreauStep, TransientCarreauWorkspace,
};
use fs_qty::Pressure;

/// One shared thermal inventory around the original retained momentum owner.
/// The accepted hydraulic history is supplied immutably to each proposal.
#[derive(Clone, Debug)]
pub struct TransientThermalCarreauGraph<'a> {
    graph: TransientCarreauGraph,
    storage: &'a TabulatedHeatStorage,
    temperature_span_k: [f64; 2],
}

/// Thermal accuracy and aggregate work bounds, alongside original hydraulics.
#[derive(Clone, Copy, Debug)]
pub struct TransientThermalCarreauControls {
    pub hydraulics: TransientCarreauControls,
    /// All complete native proposals, including initial-temperature evaluation.
    pub maximum_hydraulic_solves: usize,
    /// Shared across every thermal trial and every native pressure iteration.
    pub maximum_total_branch_solves: u64,
    pub maximum_temperature_change_k: f64,
    pub temperature_tolerance_k: f64,
    pub absolute_energy_tolerance_j: f64,
    pub relative_energy_tolerance: f64,
}

impl Default for TransientThermalCarreauControls {
    fn default() -> Self {
        Self {
            hydraulics: TransientCarreauControls::default(),
            maximum_hydraulic_solves: 96,
            maximum_total_branch_solves: 65_536,
            maximum_temperature_change_k: 5.0,
            temperature_tolerance_k: 1e-8,
            absolute_energy_tolerance_j: 1e-10,
            relative_energy_tolerance: 1e-9,
        }
    }
}

impl TransientThermalCarreauControls {
    fn validate(self) -> Result<(), ThermalCarreauError> {
        if !(1..=MAX_THERMAL_HYDRAULIC_SOLVES).contains(&self.maximum_hydraulic_solves)
            || self.maximum_total_branch_solves == 0
            || !positive(self.maximum_temperature_change_k)
            || !positive(self.temperature_tolerance_k)
            || !nonnegative(self.absolute_energy_tolerance_j)
            || !nonnegative(self.relative_energy_tolerance)
            || self.relative_energy_tolerance >= 1.0
            || self.absolute_energy_tolerance_j == 0.0 && self.relative_energy_tolerance == 0.0
        {
            return Err(
                GraphError::InvalidInput("invalid transient thermal Carreau controls").into(),
            );
        }
        // The original native proposal admits all hydraulic controls before
        // performing its first pressure or radial solve.
        Ok(())
    }
}

/// The same accepted trajectory owns all momentum, work and thermal quantities.
#[derive(Clone, Debug)]
pub struct TransientThermalCarreauReport {
    pub initial_temperature_k: f64,
    pub midpoint_temperature_k: f64,
    pub final_temperature_k: f64,
    pub stored_energy_change_j: f64,
    pub viscous_heat_j: f64,
    pub prescribed_heat_j: f64,
    pub bath_heat_j: f64,
    /// Sensible storage minus actual viscous, imposed and bath heat.
    pub thermal_residual_j: f64,
    pub thermal_residual_allowance_j: f64,
    /// Sensible + kinetic + hydraulic storage changes minus pressure/head,
    /// imposed and bath work. Original hydraulic defects remain in this value.
    pub combined_energy_residual_j: f64,
    pub combined_residual_allowance_j: f64,
    pub temperature_bracket_width_k: f64,
    pub hydraulic_solves: usize,
    pub total_branch_solves: u64,
    /// Original native report at the accepted midpoint viscosity, unchanged.
    pub hydraulics: TransientCarreauStep,
}

/// Publish temperature and retained radial/storage histories together.
#[derive(Clone, Debug)]
pub struct TransientThermalCarreauStep<'a> {
    pub next: TransientThermalCarreauGraph<'a>,
    pub state: TransientCarreauState,
    pub report: TransientThermalCarreauReport,
}

impl<'a> TransientThermalCarreauGraph<'a> {
    /// Bind one explicit well-mixed inventory, including every returning stream.
    /// Unequal regional temperatures require a spatial thermal owner and refuse.
    pub fn new(
        graph: TransientCarreauGraph,
        storage: &'a TabulatedHeatStorage,
    ) -> Result<Self, ThermalCarreauError> {
        let temperature_span_k = shared_temperature_span(graph.graph(), storage)?;
        Ok(Self {
            graph,
            storage,
            temperature_span_k,
        })
    }

    pub const fn graph(&self) -> &TransientCarreauGraph {
        &self.graph
    }
    pub const fn storage(&self) -> &TabulatedHeatStorage {
        self.storage
    }
    pub const fn temperature_span_k(&self) -> [f64; 2] {
        self.temperature_span_k
    }
    pub fn temperature_k(&self) -> f64 {
        self.graph.graph().branches()[0].pipe.temperature_k()
    }

    /// Propose one interval without changing the accepted temperature or history.
    ///
    /// Existing states/workspaces retain their native compatibility identities.
    /// Workspaces are scratch and may change on refusal. Every thermal trial
    /// restarts from the SAME accepted radial/storage history; no trial advances
    /// the clock of a later trial. Final source support and both energy gates
    /// must pass before the complete next owner and history are returned.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        state: &TransientCarreauState,
        workspace: &mut TransientCarreauWorkspace,
        duration_s: f64,
        boundaries: &[FixedPressure],
        injections: &[FlowInjection],
        boundary: ThermalCarreauBoundary,
        controls: TransientThermalCarreauControls,
        cx: &Cx<'_>,
    ) -> Result<TransientThermalCarreauStep<'a>, ThermalCarreauError> {
        self.advance_impl(
            state, workspace, duration_s, boundaries, injections, None, boundary, controls, cx,
        )
    }

    /// The same coupled proposal with one held signed pressure rise per branch.
    /// Head work can be negative during inertial reversal; it is not heat.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_heads(
        &self,
        state: &TransientCarreauState,
        workspace: &mut TransientCarreauWorkspace,
        duration_s: f64,
        boundaries: &[FixedPressure],
        injections: &[FlowInjection],
        branch_heads: &[Pressure],
        boundary: ThermalCarreauBoundary,
        controls: TransientThermalCarreauControls,
        cx: &Cx<'_>,
    ) -> Result<TransientThermalCarreauStep<'a>, ThermalCarreauError> {
        self.advance_impl(
            state,
            workspace,
            duration_s,
            boundaries,
            injections,
            Some(branch_heads),
            boundary,
            controls,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_impl(
        &self,
        state: &TransientCarreauState,
        workspace: &mut TransientCarreauWorkspace,
        duration_s: f64,
        boundaries: &[FixedPressure],
        injections: &[FlowInjection],
        branch_heads: Option<&[Pressure]>,
        boundary: ThermalCarreauBoundary,
        controls: TransientThermalCarreauControls,
        cx: &Cx<'_>,
    ) -> Result<TransientThermalCarreauStep<'a>, ThermalCarreauError> {
        checkpoint(cx)?;
        controls.validate()?;
        boundary.validate()?;
        if !positive(duration_s) {
            return Err(
                GraphError::InvalidInput("positive transient thermal duration required").into(),
            );
        }
        let initial_temperature = self.temperature_k();
        let initial_slope = finite(
            self.storage.heat_capacity_j_k(initial_temperature)?
                + (0.5 * duration_s) * boundary.bath_conductance_w_k,
            "transient Carreau thermal initial slope",
        )?;
        let mut evaluations = 0;
        let mut total_branch_solves = 0_u64;
        let mut solve = |temperature: f64| -> Result<
            Trial<(TransientCarreauState, TransientCarreauStep)>,
            ThermalCarreauError,
        > {
            checkpoint(cx)?;
            if evaluations == controls.maximum_hydraulic_solves {
                return Err(ThermalCarreauError::HydraulicSolveBudget);
            }
            let remaining = controls.maximum_total_branch_solves - total_branch_solves;
            if remaining == 0 {
                return Err(TransientCarreauError::BranchSolveBudget.into());
            }
            if !temperature.is_finite()
                || temperature < self.temperature_span_k[0]
                || temperature > self.temperature_span_k[1]
            {
                return Err(ThermalCarreauError::TemperatureDomain {
                    temperature_k: temperature,
                });
            }
            let midpoint = 0.5 * initial_temperature + 0.5 * temperature;
            let model = self.graph.at_temperature(midpoint, cx)?;
            let mut candidate = state.clone();
            let hydraulics = TransientCarreauControls {
                maximum_branch_solves: controls.hydraulics.maximum_branch_solves.min(remaining),
                ..controls.hydraulics
            };
            evaluations += 1;
            let report = match branch_heads {
                Some(heads) => model.step_with_heads(
                    &mut candidate,
                    workspace,
                    boundaries,
                    injections,
                    heads,
                    duration_s,
                    hydraulics,
                    cx,
                )?,
                None => model.step(
                    &mut candidate,
                    workspace,
                    boundaries,
                    injections,
                    duration_s,
                    hydraulics,
                    cx,
                )?,
            };
            total_branch_solves += report.branch_solves;
            let heat = report.viscous_heat_j;
            if !nonnegative(heat) {
                return Err(ThermalCarreauError::EnergyBalance);
            }
            let stored = self
                .storage
                .energy_change_j(initial_temperature, temperature)?;
            let imposed = finite(
                duration_s * boundary.prescribed_power_w,
                "transient Carreau imposed heat",
            )?;
            let bath = finite(
                duration_s * boundary.bath_power(midpoint)?,
                "transient Carreau bath heat",
            )?;
            let residual = sum(&[stored, -heat, -imposed, -bath])?;
            let scale = stored.abs().max(heat).max(imposed.abs()).max(bath.abs());
            let allowance = finite(
                controls.absolute_energy_tolerance_j + controls.relative_energy_tolerance * scale,
                "transient Carreau thermal tolerance",
            )?;
            Ok(Trial {
                point: Point {
                    temperature,
                    residual,
                    allowance,
                },
                stored,
                heat,
                imposed,
                bath,
                hydraulics: (candidate, report),
            })
        };
        let (accepted, width, _) = solve_temperature(
            initial_temperature,
            self.temperature_span_k,
            controls.maximum_temperature_change_k,
            controls.temperature_tolerance_k,
            initial_slope,
            &mut solve,
        )?;
        drop(solve);
        if accepted.point.residual.abs() > accepted.point.allowance {
            return Err(ThermalCarreauError::EnergyBalance);
        }
        // Admit the endpoint's complete original physical rate envelope without
        // taking another momentum step or replacing the retained radial field.
        let next_graph = self.graph.at_temperature(accepted.point.temperature, cx)?;
        let (next_state, hydraulics) = accepted.hydraulics;
        let terms = [
            accepted.stored,
            hydraulics.kinetic_energy_change_j,
            hydraulics.storage_energy_change_j,
            -hydraulics.reservoir_work_j,
            -hydraulics.injection_work_j,
            -hydraulics.branch_head_work_j,
            -accepted.imposed,
            -accepted.bath,
        ];
        let combined = sum(&terms)?;
        let scale = sum(&terms.map(f64::abs))?;
        let combined_allowance = finite(
            controls.absolute_energy_tolerance_j
                + controls.relative_energy_tolerance * scale
                + 64.0 * f64::EPSILON * scale,
            "transient Carreau combined energy tolerance",
        )?;
        if combined.abs() > combined_allowance {
            return Err(ThermalCarreauError::EnergyBalance);
        }
        checkpoint(cx)?;
        Ok(TransientThermalCarreauStep {
            next: Self {
                graph: next_graph,
                storage: self.storage,
                temperature_span_k: self.temperature_span_k,
            },
            state: next_state,
            report: TransientThermalCarreauReport {
                initial_temperature_k: initial_temperature,
                midpoint_temperature_k: 0.5 * initial_temperature
                    + 0.5 * accepted.point.temperature,
                final_temperature_k: accepted.point.temperature,
                stored_energy_change_j: accepted.stored,
                viscous_heat_j: accepted.heat,
                prescribed_heat_j: accepted.imposed,
                bath_heat_j: accepted.bath,
                thermal_residual_j: accepted.point.residual,
                thermal_residual_allowance_j: accepted.point.allowance,
                combined_energy_residual_j: combined,
                combined_residual_allowance_j: combined_allowance,
                temperature_bracket_width_k: width,
                hydraulic_solves: evaluations,
                total_branch_solves,
                hydraulics,
            },
        })
    }
}

#[cfg(test)]
mod tests;
