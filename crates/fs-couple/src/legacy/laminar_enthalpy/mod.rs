//! Source-owned BDM laminar momentum, enthalpy advection and finite thermal feedback.
//!
//! Every implicit thermal trial resolves the original convective/rheological
//! BDF1 momentum equations from the SAME accepted velocity. Its exact sign-split
//! BDM face rates drive the native upwind enthalpy owner. Original bulk strain
//! heat and explicitly routed Navier-slip heat are deposited once. SIP, upwind,
//! backflow and BDF numerical losses are not material heat.
//!
//! This is a fixed-mesh, fixed-density incompressible composition. It preserves
//! native caloric laws and authored conduction/contact/radiation links; it does
//! not manufacture enthalpy from a heat-capacity sample, infer thermal expansion,
//! or qualify turbulent flow. Coupled accuracy is first order in time. A thermal
//! balance and momentum residual are not a whole-system passivity certificate.

use std::cell::Cell;

use fs_exec::Cx;
use fs_flux::{
    FluxSystem, MaterialFluxError, MaterialFluxNumerics, MaterialFluxRheologyPoint,
    MaterialFluxSystem, MaterialNavierStokesBoundary, MaterialNavierStokesRheologyStep,
    NavierStokesControls,
};
use fs_material::MaterialError;
use fs_material::fluid_transport::{
    FluidTransport, FluidTransportState, ThermalFluidTransport, ThermalFluidTransportState,
};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport};
use fs_thermal::enthalpy::backward_euler::{
    BackwardEulerStep, BoundaryStepCoefficients, StepCoefficients, TransportStepCoefficients,
};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionSimulation, EnergyAudit, MaterialId,
};

use crate::stokes_motion::heating::{StokesHeatFailure, StokesSlipHeatRoute, StokesThermalCells};

mod binding;
use binding::ThermalBinding;

/// Joint native momentum/enthalpy temporal accuracy and whole-window acceptance.
pub mod adaptive;

/// Source-owned spatial mechanical drives resolved at the actual thermal endpoint.
pub mod drive;

/// Original RCR/inertance ports resolved at the same thermal endpoint.
pub mod windkessel;

/// Source-density buoyancy on the same coupled native thermal/momentum owner.
pub mod boussinesq;
use boussinesq::LaminarBuoyancyReport;

/// The accepted BDM velocity and its physical clock; neither is mutated.
#[derive(Clone, Copy)]
pub struct LaminarEnthalpyHistory<'a> {
    pub velocity_m_s: &'a [f64],
    pub time_s: f64,
}

/// Original physical momentum conditions held over this thermal interval.
#[derive(Clone, Copy)]
pub struct LaminarEnthalpyInput<'a> {
    pub history: LaminarEnthalpyHistory<'a>,
    /// N/m³; this additional force is not inferred from temperature.
    pub force_density: &'a dyn Fn([f64; 2]) -> [f64; 2],
    pub boundaries: &'a [MaterialNavierStokesBoundary<'a>],
    /// One explicit route for every Navier-slip edge, including perfect slip.
    /// The original receiver/fraction admission is shared with Stokes heating.
    pub slip_heat_routes: &'a [StokesSlipHeatRoute],
}

/// Explicit incoming caloric state on an original boundary edge.
/// Outlet enthalpy is always the actual local donor's original value.
#[derive(Clone, Copy, Debug)]
pub struct LaminarInletEnthalpy {
    pub material_id: MaterialId,
    pub specific_h_j_kg: f64,
}

/// Nonadvective heat and source-owned external streams at the SAME endpoint.
pub struct LaminarHeatExchange<T> {
    /// Original edge order. Interior edges require None. Every edge carrying
    /// any inflow requires Some, even if its NET flow is outward or zero.
    /// A dormant boundary entry may remain for flow reversal during iteration.
    pub inlet_enthalpies: Vec<Option<LaminarInletEnthalpy>>,
    /// Additional sources and replacements of original physical heat links.
    /// Exclude advection and fluid/slip dissipation, which this owner supplies.
    pub heat: BoundaryStepCoefficients<T>,
}

/// Original source object and its checked, coherent endpoint transport view.
/// Native momentum reads this immutable snapshot without re-querying sources.
pub struct LaminarThermalMaterial<G> {
    source: G,
    state: ThermalFluidTransportState,
}
impl<G> LaminarThermalMaterial<G> {
    pub fn source(&self) -> &G {
        &self.source
    }
    pub const fn thermal_state(&self) -> &ThermalFluidTransportState {
        &self.state
    }
}
impl<G> FluidTransport for LaminarThermalMaterial<G> {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.state.transport())
    }
}

/// Bounded native momentum and nonlinear thermal work, with independent
/// physical thermal acceptance. The finite-difference coordinate is joules.
#[derive(Clone, Copy, Debug)]
pub struct LaminarEnthalpyControls {
    pub momentum: NavierStokesControls,
    pub maximum_active_cells: usize,
    pub maximum_endpoint_evaluations: usize,
    pub maximum_thermal_iterations: usize,
    pub difference_step_j: f64,
    pub cell_tolerance_j: f64,
    pub total_tolerance_j: f64,
    /// Original nonlinear policy. Use an inexact linear forcing floor suited
    /// to numerical directional derivatives (for example `forcing_minimum =
    /// 1e-4`), rather than requesting machine-precision Krylov convergence of
    /// a finite-difference action. Final physical heat tolerances remain exact.
    pub newton: NewtonKrylovConfig,
}

/// Original refusal before any replacement physical history is returned.
#[derive(Debug)]
pub enum LaminarEnthalpyError<E = DepositionError> {
    Invalid(&'static str),
    Momentum(MaterialFluxError),
    Thermal(DepositionError),
    Exchange(E),
    Drive(E),
    FaceTransport(fs_flux::ns::transport::FluxTransportError),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    EndpointBudget,
    Cancelled,
}
impl<E> From<DepositionError> for LaminarEnthalpyError<E> {
    fn from(error: DepositionError) -> Self {
        match error {
            DepositionError::Cancelled => Self::Cancelled,
            other => Self::Thermal(other),
        }
    }
}
impl<E> From<StokesHeatFailure> for LaminarEnthalpyError<E> {
    fn from(error: StokesHeatFailure) -> Self {
        match error {
            StokesHeatFailure::Binding(reason) => Self::Invalid(reason),
            StokesHeatFailure::Flow(error) => Self::Momentum(error),
            StokesHeatFailure::Thermal(error) => error.into(),
            StokesHeatFailure::Cancelled => Self::Cancelled,
        }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for LaminarEnthalpyError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "laminar enthalpy coupling: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for LaminarEnthalpyError<E> {}

/// Independent thermal storage/source observations at the accepted endpoint.
/// Rates use backward-Euler endpoint quadrature over the actual interval.
pub struct LaminarEnthalpyReport<T> {
    pub thermal: EnergyAudit,
    pub exchange: T,
    pub nonlinear: NewtonReport,
    pub endpoint_evaluations: usize,
    pub maximum_cell_residual_j: f64,
    pub bulk_viscous_heat_j: f64,
    pub slip_heat_j: f64,
    pub external_stream_enthalpy_j: f64,
    pub other_external_heat_j: f64,
    /// Delta H - bulk heat - slip heat - net stream enthalpy - other heat.
    /// Internal advective/conductive transfers cancel; nothing is repaired.
    pub thermal_balance_residual_j: f64,
    /// Present only for explicitly requested source-density buoyancy. Its
    /// mechanical work is separate from every heat term above.
    pub buoyancy: Option<LaminarBuoyancyReport>,
    /// Selected actual mechanical exchanges for explicitly resolved drives.
    /// These are separate from thermal storage and never deposited as heat.
    pub mechanical: Option<drive::LaminarMechanicalWork>,
}

/// Joint accepted flow/material and thermal histories. A failed attempt never
/// changes either input, including after successful provisional flow solves.
pub struct LaminarEnthalpyCandidate<'m, G, T> {
    thermal: DepositionSimulation<'m>,
    flow: MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>>,
    report: LaminarEnthalpyReport<T>,
}
impl<'m, G, T> LaminarEnthalpyCandidate<'m, G, T> {
    pub fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub fn flow(&self) -> &MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>> {
        &self.flow
    }
    pub fn report(&self) -> &LaminarEnthalpyReport<T> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'m>,
        MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>>,
        LaminarEnthalpyReport<T>,
    ) {
        (self.thermal, self.flow, self.report)
    }
}

/// One original continuum source and its fixed finite caloric inventories.
pub struct LaminarEnthalpy<'s, 'mesh, 'source, F: FluidTransport + ?Sized> {
    system: &'s FluxSystem<'mesh>,
    flow: MaterialFluxSystem<'s, 'mesh, 'source, F>,
    binding: ThermalBinding,
}

struct Endpoint<G, T> {
    flow: MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>>,
    exchange: T,
    powers_w: [f64; 4],
    buoyancy: Option<LaminarBuoyancyReport>,
}

impl<'s, 'mesh, 'source, F: FluidTransport + ?Sized> LaminarEnthalpy<'s, 'mesh, 'source, F> {
    /// Bind one active triangular prism per original fluid cell. Additional
    /// native thermal cells may carry walls, solids and external heat exchange.
    /// All fluid cells must retain one original identity-bearing ENTHALPY law.
    /// The caller associates that caloric law with the original fluid source.
    /// No caloric model is inferred from the momentum transport coefficients.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        cx: &Cx<'_>,
        system: &'s FluxSystem<'mesh>,
        material: &'source F,
        numerics: MaterialFluxNumerics,
        thermal: &DepositionSimulation<'_>,
        cells: StokesThermalCells<'_>,
        depth_m: f64,
        maximum_connections: usize,
    ) -> Result<Self, LaminarEnthalpyError> {
        checkpoint(cx)?;
        let flow = MaterialFluxSystem::new(system, material, numerics)
            .map_err(LaminarEnthalpyError::Momentum)?;
        let binding = ThermalBinding::new(
            cx,
            system,
            thermal,
            cells,
            depth_m,
            flow.transport().density_kg_m3(),
            maximum_connections,
        )?;
        Ok(Self {
            system,
            flow,
            binding,
        })
    }

    pub fn triangle_cells(&self) -> &[CellId] {
        self.binding.cells()
    }

    /// Resolve thermal advection, source viscosity, strain and convection at
    /// one fully coupled BDF1 endpoint. Every flow solve uses the SAME accepted
    /// velocity and held physical boundary data; no provisional field becomes
    /// history. Source queries receive actual native caloric temperature and
    /// current objective shear, and retain their exact original material owner.
    ///
    /// The exchange callback sees that converged physical flow and supplies
    /// explicit incoming enthalpies plus nonadvective thermal laws. Conductivity
    /// may be resolved through the returned cell sources and original thermal
    /// links; no centre-distance Fourier approximation is manufactured here.
    ///
    /// Fixed density, geometry and caloric identities remain exact invariants.
    /// Native fixed-mass advection can refuse a continuity error admitted by a
    /// looser momentum tolerance: no mass-flux correction is applied. Scheduled
    /// thermal events retain their original interval admission. Callbacks must
    /// be pure and valid on nonlinear derivative probes. Cancellation is checked
    /// at composition/source boundaries; native sparse solves keep their own
    /// bounded work contract. No within-interval flow-reversal location or
    /// timestep-accuracy certificate is implied by nonlinear convergence.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'m, G: ThermalFluidTransport, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        input: LaminarEnthalpyInput<'_>,
        controls: LaminarEnthalpyControls,
        energy_guess_j: &[f64],
        resolve_material: impl Fn(MaterialFluxRheologyPoint, f64, &F) -> Result<G, MaterialError>,
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>>,
        ) -> Result<LaminarHeatExchange<T>, E>,
    ) -> Result<LaminarEnthalpyCandidate<'m, G, T>, LaminarEnthalpyError<E>> {
        self.advance_with_flow(
            cx,
            thermal,
            target_s,
            input,
            controls,
            energy_guess_j,
            |temperatures, duration_s| {
                self.resolve_laminar_flow(
                    cx,
                    input,
                    controls,
                    temperatures,
                    duration_s,
                    &resolve_material,
                )
                .map(|flow| (flow, None))
            },
            exchange,
        )
    }

    // Both held and resolved drives use this exact original momentum source path.
    #[allow(clippy::too_many_arguments)]
    fn resolve_laminar_flow<G: ThermalFluidTransport>(
        &self,
        cx: &Cx<'_>,
        input: LaminarEnthalpyInput<'_>,
        controls: LaminarEnthalpyControls,
        temperatures: &[f64],
        duration_s: f64,
        resolve_material: &impl Fn(MaterialFluxRheologyPoint, f64, &F) -> Result<G, MaterialError>,
    ) -> Result<MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>>, MaterialFluxError>
    {
        self.flow.bdf1_navier_stokes_mixed_rheology_step(
            input.force_density,
            input.boundaries,
            input.history.velocity_m_s,
            duration_s,
            controls.momentum,
            |point, reference| {
                self.resolve_laminar_material(
                    cx,
                    point,
                    temperatures[point.cell],
                    reference,
                    resolve_material,
                )
            },
        )
    }

    // The original thermal source and density admission are shared by every
    // fixed-density endpoint, including dynamically loaded hydraulic ports.
    fn resolve_laminar_material<G: ThermalFluidTransport>(
        &self,
        cx: &Cx<'_>,
        point: MaterialFluxRheologyPoint,
        temperature: f64,
        reference: &FluidTransportState,
        resolve_material: &impl Fn(MaterialFluxRheologyPoint, f64, &F) -> Result<G, MaterialError>,
    ) -> Result<LaminarThermalMaterial<G>, MaterialError> {
        if cx.checkpoint().is_err() {
            return Err(material_error("coupled flow/enthalpy operation cancelled"));
        }
        let source = resolve_material(point, temperature, self.flow.material())?;
        let material = thermal_source(source, temperature)?;
        if material.state.transport().density_kg_m3() != reference.density_kg_m3() {
            return Err(material_error(
                "original fixed-density thermal flow must retain reference density",
            ));
        }
        Ok(material)
    }

    // Held inputs retain their existing entry point and original routes.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn advance_with_flow<'m, G, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        input: LaminarEnthalpyInput<'_>,
        controls: LaminarEnthalpyControls,
        energy_guess_j: &[f64],
        resolve_flow: impl Fn(
            &[f64],
            f64,
        ) -> Result<
            (
                MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>>,
                Option<LaminarBuoyancyReport>,
            ),
            MaterialFluxError,
        >,
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>>,
        ) -> Result<LaminarHeatExchange<T>, E>,
    ) -> Result<LaminarEnthalpyCandidate<'m, G, T>, LaminarEnthalpyError<E>> {
        checkpoint(cx)?;
        let routes = self.binding.slip_routes(
            cx,
            self.system,
            thermal,
            input.boundaries,
            input.slip_heat_routes,
        )?;
        self.advance_with_endpoint(
            cx,
            thermal,
            target_s,
            input.history,
            controls,
            energy_guess_j,
            |candidate, temperatures, duration_s| {
                let flow = resolve_flow(temperatures, duration_s);
                checkpoint(cx)?;
                let (flow, buoyancy) = flow.map_err(LaminarEnthalpyError::Momentum)?;
                let response =
                    exchange(cx, candidate, &flow).map_err(LaminarEnthalpyError::Exchange)?;
                checkpoint(cx)?;
                self.thermal_endpoint(cx, candidate, flow, buoyancy, &routes, response)
            },
        )
    }

    // Same original heat/advection construction for static and resolved drives.
    #[allow(clippy::too_many_arguments)]
    fn thermal_endpoint<G, T, E>(
        &self,
        cx: &Cx<'_>,
        candidate: &DepositionSimulation<'_>,
        flow: MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>>,
        buoyancy: Option<LaminarBuoyancyReport>,
        routes: &std::collections::BTreeMap<usize, StokesSlipHeatRoute>,
        response: LaminarHeatExchange<T>,
    ) -> Result<TransportStepCoefficients<Endpoint<G, T>>, LaminarEnthalpyError<E>> {
        let (coefficients, powers_w) = self.binding.coefficients(
            cx,
            self.system,
            candidate,
            &flow.step.solution,
            routes,
            response,
        )?;
        let TransportStepCoefficients { heat, transfers } = coefficients;
        let BoundaryStepCoefficients { bulk, boundaries } = heat;
        let StepCoefficients {
            sources_w,
            conductances_w_k,
            payload,
        } = bulk;
        Ok(TransportStepCoefficients {
            heat: BoundaryStepCoefficients {
                bulk: StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload: Endpoint {
                        flow,
                        exchange: payload,
                        powers_w,
                        buoyancy,
                    },
                },
                boundaries,
            },
            transfers,
        })
    }

    // The thermal residual, original storage and final transactional gates
    // are shared by every endpoint mechanic. No additional thermal solve.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn advance_with_endpoint<'m, G, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        history: LaminarEnthalpyHistory<'_>,
        controls: LaminarEnthalpyControls,
        energy_guess_j: &[f64],
        resolve_endpoint: impl Fn(
            &DepositionSimulation<'m>,
            &[f64],
            f64,
        ) -> Result<
            TransportStepCoefficients<Endpoint<G, T>>,
            LaminarEnthalpyError<E>,
        >,
    ) -> Result<LaminarEnthalpyCandidate<'m, G, T>, LaminarEnthalpyError<E>> {
        checkpoint(cx)?;
        let duration_s = target_s - thermal.time_s();
        if history.time_s != thermal.time_s()
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || controls.maximum_active_cells == 0
            || thermal.cells().len() > controls.maximum_active_cells
            || controls.maximum_endpoint_evaluations == 0
            || controls.maximum_thermal_iterations == 0
            || !controls.difference_step_j.is_finite()
            || controls.difference_step_j <= 0.0
            || !nonnegative(controls.cell_tolerance_j)
            || !nonnegative(controls.total_tolerance_j)
        {
            return Err(LaminarEnthalpyError::Invalid(
                "matching clocks, positive interval and finite bounded coupled controls required",
            ));
        }
        controls
            .newton
            .validate()
            .map_err(LaminarEnthalpyError::Newton)?;
        self.binding.check(cx, self.system, thermal)?;
        let step = BackwardEulerStep::new(thermal, target_s, || thermal_checkpoint(cx))?;
        if energy_guess_j.len() != step.cell_ids().len()
            || energy_guess_j.iter().any(|x| !x.is_finite())
        {
            return Err(LaminarEnthalpyError::Invalid(
                "complete finite original cell-energy guess required",
            ));
        }
        let evaluations = Cell::new(0usize);
        let evaluate = |energies: &[f64]| {
            checkpoint(cx)?;
            if evaluations.get() == controls.maximum_endpoint_evaluations {
                return Err(LaminarEnthalpyError::EndpointBudget);
            }
            evaluations.set(evaluations.get() + 1);
            let energies = step.resolvable_increments_j(energies);
            step.evaluate_with_transport(
                &energies,
                |candidate| {
                    self.binding.check(cx, self.system, candidate)?;
                    let temperatures = self.binding.temperatures(cx, candidate)?;
                    resolve_endpoint(candidate, &temperatures, duration_s)
                },
                || thermal_checkpoint(cx),
            )
        };
        let difference = FallibleFiniteDifference::new(
            step.cell_ids().len(),
            controls.difference_step_j,
            |energies: &[f64]| evaluate(energies).map(|trial| trial.residual_j().to_vec()),
        )
        .map_err(difference_error)?;
        let initialized =
            NewtonKrylovState::new(&difference, energy_guess_j.to_vec(), controls.newton);
        if let Some(error) = difference.take_failure() {
            return Err(difference_error(error));
        }
        let mut solver = initialized.map_err(LaminarEnthalpyError::Newton)?;
        let run = solver.run_controlled(
            &difference,
            controls.maximum_thermal_iterations,
            &mut || difference.has_failed() || cx.checkpoint().is_err(),
        );
        if let Some(error) = difference.take_failure() {
            return Err(difference_error(error));
        }
        checkpoint(cx)?;
        if !run.report.converged {
            return Err(LaminarEnthalpyError::NotConverged(run.report));
        }
        let final_trial = evaluate(&solver.x)?;
        let maximum_cell_residual_j = final_trial.maximum_cell_residual_j();
        let (next, endpoint, audit) = final_trial.accept(
            controls.cell_tolerance_j,
            controls.total_tolerance_j,
            || thermal_checkpoint(cx),
        )?;
        let energies = endpoint.powers_w.map(|power| power * duration_s);
        let residual = binding::checked_sum(&[
            audit.stored_energy_change_j,
            -energies[0],
            -energies[1],
            -energies[2],
            -energies[3],
        ])?;
        let scale = binding::checked_sum(&[
            audit.stored_energy_change_j.abs(),
            energies[0].abs(),
            energies[1].abs(),
            energies[2].abs(),
            energies[3].abs(),
        ])?;
        let tolerance = controls.total_tolerance_j + 64.0 * f64::EPSILON * scale;
        if !tolerance.is_finite() || residual.abs() > tolerance {
            return Err(DepositionError::EnergyBalance(residual).into());
        }
        checkpoint(cx)?;
        Ok(LaminarEnthalpyCandidate {
            thermal: next,
            flow: endpoint.flow,
            report: LaminarEnthalpyReport {
                thermal: audit,
                exchange: endpoint.exchange,
                nonlinear: run.report,
                endpoint_evaluations: evaluations.get(),
                maximum_cell_residual_j,
                bulk_viscous_heat_j: energies[0],
                slip_heat_j: energies[1],
                external_stream_enthalpy_j: energies[2],
                other_external_heat_j: energies[3],
                thermal_balance_residual_j: residual,
                buoyancy: endpoint.buoyancy,
                mechanical: None,
            },
        })
    }
}

// Freeze coherent original views once. Boussinesq projects THIS immutable
// snapshot, not a second source query that could change density or viscosity.
fn thermal_source<G: ThermalFluidTransport>(
    source: G,
    temperature_k: f64,
) -> Result<LaminarThermalMaterial<G>, MaterialError> {
    let state = source.thermal_transport_state()?;
    if state.transport() != source.transport_state()?
        || (state.temperature_k() - temperature_k).abs()
            > 64.0 * f64::EPSILON * temperature_k.abs().max(state.temperature_k().abs())
    {
        return Err(material_error(
            "source temperature and thermal/momentum transport must match the actual endpoint",
        ));
    }
    Ok(LaminarThermalMaterial { source, state })
}

fn material_error(what: &str) -> MaterialError {
    MaterialError::State { what: what.into() }
}
fn thermal_checkpoint(cx: &Cx<'_>) -> Result<(), DepositionError> {
    cx.checkpoint().map_err(|_| DepositionError::Cancelled)
}
fn checkpoint<E>(cx: &Cx<'_>) -> Result<(), LaminarEnthalpyError<E>> {
    thermal_checkpoint(cx).map_err(Into::into)
}
fn difference_error<E>(
    error: FiniteDifferenceError<LaminarEnthalpyError<E>>,
) -> LaminarEnthalpyError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => LaminarEnthalpyError::Invalid(message),
    }
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}

#[cfg(test)]
mod tests;
