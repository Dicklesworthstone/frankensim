//! Thermal-only nonlinear endpoints through the existing complete heat owner.
//!
//! All native enthalpy residuals, Newton work and acceptance remain in
//! `cell_driver_candidate`. No circuit or independent thermal state is created.

use super::{HeatAssembly, HeatNetworkError, as_circuit_error, closed_balance};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, DrivenCoefficients, DriverAccount, ThermalDriver,
    cell_driver_candidate,
};
use crate::material_circuit::MaterialCircuitSolve;
use fs_exec::Cx;
use fs_solver::{NewtonKrylovConfig, NewtonReport};
use fs_thermal::enthalpy::{DepositionSimulation, EnergyAudit};
use std::collections::BTreeMap;

/// Bounded nonlinear work and independent physical gates for one thermal step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatAssemblySolve {
    pub maximum_cells: usize,
    /// Every assembly query, including Jacobian probes and fresh final admission.
    pub maximum_evaluations: usize,
    pub maximum_newton_iterations: usize,
    /// Original numerical Jacobian perturbation in cell energy coordinates [J].
    pub difference_step_j: f64,
    /// Largest local cell residual; opposite errors cannot cancel this gate [J].
    pub cell_energy_tolerance_j: f64,
    /// Native total thermal residual, including numerical enclosure net heat [J].
    pub thermal_energy_tolerance_j: f64,
    /// Thermal storage minus actual external supply, excluding enclosure net [J].
    pub closed_energy_tolerance_j: f64,
    /// The existing Newton--Krylov policy; residual coordinates are joules.
    pub newton: NewtonKrylovConfig,
}

impl HeatAssemblySolve {
    fn native(self) -> CircuitEnthalpySolve {
        CircuitEnthalpySolve {
            maximum_cells: self.maximum_cells,
            maximum_circuit_trials: self.maximum_evaluations,
            maximum_newton_iterations: self.maximum_newton_iterations,
            difference_step_j: self.difference_step_j,
            cell_energy_tolerance_j: self.cell_energy_tolerance_j,
            thermal_energy_tolerance_j: self.thermal_energy_tolerance_j,
            total_energy_tolerance_j: self.closed_energy_tolerance_j,
            // The shared cell driver does not read circuit controls or invoke
            // a circuit solve. These inactive fields never enter the public API.
            circuit: MaterialCircuitSolve {
                energy_tolerance_j: 0.0,
                power_tolerance_w: 0.0,
            },
            newton: self.newton,
        }
    }
}

/// Actual endpoint material sources and native thermal/nonlinear accounts.
pub struct HeatAssemblyStepReport<S> {
    /// Original finite-cell audit; external heat includes numerical enclosure net.
    pub thermal: EnergyAudit,
    /// Complete fresh assembly at the returned physical endpoint.
    pub source: S,
    pub nonlinear: NewtonReport,
    pub maximum_cell_residual_j: f64,
    /// All source queries, including discarded probes and final admission.
    pub evaluations: usize,
    /// Authored source/boundary heat with closed-enclosure net removed once [J].
    pub external_heat_j: f64,
    /// Thermal storage change minus actual external heat [J].
    pub energy_residual_j: f64,
}

// A pure thermal consumer has no external driver store, supplied work or heat.
// This private empty account selects the unchanged native thermal residual and
// does not construct, advance or impersonate an electrical/mechanical model.
struct ThermalOnly;
impl ThermalDriver for ThermalOnly {
    type Model = ();
    type Step = ();

    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: 0.0,
            numerical_dissipation_j: 0.0,
            stored_change_j: 0.0,
            supplied_j: 0.0,
        }
    }

    fn into_parts(self) -> ((), ()) {
        ((), ())
    }
}

/// Solve a complete material heat assembly on its original finite-cell stores.
///
/// Every nonlinear trial, numerical Jacobian probe and final acceptance queries
/// the SAME assembly at that trial's actual enthalpies. Source-dependent wall
/// conduction, resolved fluid-facing surfaces, convection, contacts and radiation
/// therefore respond to both finite temperatures inside the native endpoint
/// solve. Reciprocal link heat enters the existing conductance map exactly once;
/// prescribed reservoirs remain external source heat.
///
/// The prepared owner must be thermal-only. Electrical bindings and unexpected
/// electrical coefficients refuse, so no authored work is silently discarded.
/// Other coupled heat owners retain their own circuit/field/mechanical drivers.
///
/// Guesses are cell energy increments [J] in sorted CellId order. Native phase,
/// source, pending-event and time admission remain authoritative. Process start
/// events first; interior events refuse, and target events remain pending for
/// the existing timeline after acceptance. Geometry, mass and material history
/// stay fixed through this one interval. Failed work publishes no endpoint.
///
/// This is the existing first-order backward-Euler/Newton--Krylov path, with its
/// bounded numerical Jacobian. It supplies no new temporal-error estimator or
/// unconditional convergence claim; refine steps and source tolerances for the
/// physical problem. Local, native global and corrected closed-energy gates are
/// independent of the numerical solver's convergence criterion.
pub fn step_heat_assembly_candidate<'m, A: HeatAssembly + ?Sized>(
    assembly: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    controls: HeatAssemblySolve,
    guess_j: &[f64],
) -> Result<
    (DepositionSimulation<'m>, HeatAssemblyStepReport<A::Sample>),
    CircuitEnthalpyError<A::Error>,
> {
    cx.checkpoint()
        .map_err(|_| CircuitEnthalpyError::Cancelled)?;
    assembly
        .prepared_heat_network()
        .admit_thermal_only()
        .map_err(|error| as_circuit_error(error).map_source(A::Error::from))?;
    let (next, (), report) = cell_driver_candidate(
        thermal,
        target_s,
        controls.native(),
        guess_j,
        |state, _| {
            let coefficients = assembly.resolve_heat(cx, state)?;
            if !coefficients.bulk.resistors.is_empty() || !coefficients.contacts.is_empty() {
                return Err(CircuitEnthalpyError::Source(
                    HeatNetworkError::Invalid(
                        "thermal-only assembly cannot discard electrical coefficients",
                    )
                    .into(),
                ));
            }
            Ok(DrivenCoefficients {
                circuit: ThermalOnly,
                sources_w: coefficients.bulk.sources_w,
                conductances_w_k: coefficients.bulk.conductances_w_k,
                source: coefficients.bulk.source,
                circuit_trials: 1,
            })
        },
        |_, _, checkpoint| {
            checkpoint()?;
            Ok(BTreeMap::new())
        },
        || cx.checkpoint().is_err(),
    )?;
    let duration_s = target_s - thermal.time_s();
    closed_balance(
        duration_s,
        report.total_energy_residual_j,
        report.source.as_ref(),
        controls.closed_energy_tolerance_j,
    )
    .map_err(|error| as_circuit_error(error).map_source(A::Error::from))?;
    let internal_radiation_j = duration_s * report.source.as_ref().radiation_net_power_w();
    let external_heat_j = report.thermal.external_heat_j - internal_radiation_j;
    let energy_residual_j = report.total_energy_residual_j + internal_radiation_j;
    if !external_heat_j.is_finite() || !energy_residual_j.is_finite() {
        return Err(CircuitEnthalpyError::Invalid(
            "thermal-only source energy account is unrepresentable",
        ));
    }
    cx.checkpoint()
        .map_err(|_| CircuitEnthalpyError::Cancelled)?;
    Ok((
        next,
        HeatAssemblyStepReport {
            thermal: report.thermal,
            source: report.source,
            nonlinear: report.nonlinear,
            maximum_cell_residual_j: report.maximum_cell_residual_j,
            evaluations: report.circuit_trials,
            external_heat_j,
            energy_residual_j,
        },
    ))
}
