//! Passive external hydraulic storage on the original implicit thermal endpoint.
//!
//! Every thermal residual resolves native convection, source rheology and the
//! RCR/inertance equations from the same accepted fluid/compliance history.
//! Resistor loss and numerical BDF defect stay outside the caloric source map.

use super::*;
use fs_flux::{
    MaterialStokesWindkesselBalance, MaterialStokesWindkesselHistory, MaterialStokesWindkesselPort,
    MaterialStokesWindkesselSolution,
};

/// Joint adaptive momentum, caloric and external hydraulic continuation.
pub mod adaptive;

/// Held fluid data and external circuit data for one physical time interval.
///
/// Compliance pressures share the clock of `fluid.history` and are measured
/// relative to the SAME reservoir references as the accepted prior state.
/// Port order, geometry, reference density and C/L/R remain fixed. Reservoir
/// increments are held over this interval; changing a compliance reference
/// requires an explicitly transformed physical history, never a silent reset.
#[derive(Clone, Copy)]
pub struct LaminarWindkesselInput<'a> {
    pub fluid: LaminarEnthalpyInput<'a>,
    /// Original disjoint open-boundary edge sets, in native port order.
    pub ports: &'a [&'a [usize]],
    /// Original SI RCR loads, with optional nonnegative series inertance.
    pub loads: &'a [MaterialStokesWindkesselPort],
    /// Accepted previous compliance pressures relative to their reservoirs [Pa].
    /// Previous Q is integrated canonically from `fluid.history.velocity_m_s`.
    pub previous_compliance_pressures_pa: &'a [f64],
}

#[derive(Clone, Copy, Debug)]
pub struct LaminarWindkesselControls {
    pub coupled: LaminarEnthalpyControls,
    pub absolute_pressure_tolerance_pa: f64,
    pub absolute_flow_tolerance_m2_s: f64,
}

/// Exact native accepted circuit state and per-depth physical receipts.
///
/// Positive Q leaves the resolved fluid and enters its external load. Native
/// port power uses the added pressure increment times that signed Q; it does
/// not include the caller's original traction offset. Multiply native energies
/// by `depth_m` and native powers by `depth_m * duration_s` for J. These are
/// external-circuit accounts, not a certificate for the entire fluid/heat system.
#[derive(Debug)]
pub struct LaminarWindkesselReport {
    pub duration_s: f64,
    pub depth_m: f64,
    pub pressure_offsets_pa: Vec<f64>,
    pub outward_volume_rates_m2_s: Vec<f64>,
    /// Retains accepted compliance, storage, resistor losses, reservoir/port
    /// work, original dynamic-law residuals and the separate BDF energy defect.
    pub balances: Vec<MaterialStokesWindkesselBalance>,
    pub maximum_pressure_residual_pa: f64,
    pub maximum_flow_residual_m2_s: f64,
}
impl LaminarWindkesselReport {
    /// Accepted compliance coordinates for the next paired physical history.
    /// Pair these with the SAME candidate's velocity and thermal clock.
    pub fn compliance_pressures_pa(&self) -> Vec<f64> {
        self.balances
            .iter()
            .map(|balance| balance.compliance_pressure_pa)
            .collect()
    }
}

/// Both source-owned thermal payload and the actual accepted hydraulic state.
pub struct LaminarWindkesselExchange<T> {
    pub thermal: T,
    pub circuit: LaminarWindkesselReport,
}

pub type LaminarWindkesselCandidate<'m, G, T> =
    LaminarEnthalpyCandidate<'m, G, LaminarWindkesselExchange<T>>;

impl<'s, 'mesh, 'source, F: FluidTransport + ?Sized> LaminarEnthalpy<'s, 'mesh, 'source, F> {
    /// Resolve original nonlinear laminar momentum, RCR/inertance storage and
    /// finite caloric transport together at one backward-Euler endpoint.
    ///
    /// Every thermal trial uses its actual original temperature/shear-dependent
    /// material and the SAME accepted velocity/compliance history. The native
    /// load owner enforces both original dynamic laws and momentum residual;
    /// its pressure response admits an entirely circuit-connected component.
    /// The original sign-split BDM rates drive enthalpy transfer, including
    /// reversal. Supply explicit inlet enthalpy for every incoming branch.
    ///
    /// `exchange` borrows the full native candidate, including its original
    /// material sources and circuit state. It may author a physical heat route;
    /// this owner never guesses a recipient for external resistor heat. Native
    /// bulk/slip heat remains deposited once; external storage and BDF defects
    /// are never added as heat. Callbacks must be pure and cannot accept or
    /// advance external physical history during thermal/Newton probes.
    ///
    /// Native momentum, endpoint budgets, cancellation and independent thermal
    /// acceptance are transactional: no replacement flow, compliance state or
    /// thermal inventory is returned until every original gate succeeds. This
    /// entry is BDF1 consistently for all three owners; it does not advance
    /// Boussinesq forcing, variable density, a moving mesh or an adaptive window.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_windkessel<'m, G: ThermalFluidTransport, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        input: LaminarWindkesselInput<'_>,
        controls: LaminarWindkesselControls,
        energy_guess_j: &[f64],
        resolve_material: impl Fn(MaterialFluxRheologyPoint, f64, &F) -> Result<G, MaterialError>,
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &MaterialStokesWindkesselSolution<
                MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>>,
            >,
        ) -> Result<LaminarHeatExchange<T>, E>,
    ) -> Result<LaminarWindkesselCandidate<'m, G, T>, LaminarEnthalpyError<E>> {
        checkpoint(cx)?;
        let routes = self.binding.slip_routes(
            cx,
            self.system,
            thermal,
            input.fluid.boundaries,
            input.fluid.slip_heat_routes,
        )?;
        let previous = MaterialStokesWindkesselHistory {
            velocity: input.fluid.history.velocity_m_s,
            compliance_pressures_pa: input.previous_compliance_pressures_pa,
        };
        self.advance_with_endpoint(
            cx,
            thermal,
            target_s,
            input.fluid.history,
            controls.coupled,
            energy_guess_j,
            |candidate, temperatures, duration_s| {
                let solved = self
                    .flow
                    .bdf2_navier_stokes_mixed_rheology_variable_step_windkessel_ports(
                        input.fluid.force_density,
                        input.fluid.boundaries,
                        input.ports,
                        input.loads,
                        previous,
                        None,
                        duration_s,
                        controls.coupled.momentum,
                        controls.absolute_pressure_tolerance_pa,
                        controls.absolute_flow_tolerance_m2_s,
                        |point, reference| {
                            self.resolve_laminar_material(
                                cx,
                                point,
                                temperatures[point.cell],
                                reference,
                                &resolve_material,
                            )
                        },
                    );
                checkpoint(cx)?;
                let solved = solved.map_err(LaminarEnthalpyError::Momentum)?;
                let response =
                    exchange(cx, candidate, &solved).map_err(LaminarEnthalpyError::Exchange)?;
                checkpoint(cx)?;
                let MaterialStokesWindkesselSolution {
                    flow,
                    pressure_offsets_pa,
                    outward_volume_rates_m2_s,
                    balances,
                    maximum_pressure_residual_pa,
                    maximum_flow_residual_m2_s,
                } = solved;
                let LaminarHeatExchange {
                    inlet_enthalpies,
                    heat,
                } = response;
                let BoundaryStepCoefficients { bulk, boundaries } = heat;
                let StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload,
                } = bulk;
                self.thermal_endpoint(
                    cx,
                    candidate,
                    flow,
                    None,
                    &routes,
                    LaminarHeatExchange {
                        inlet_enthalpies,
                        heat: BoundaryStepCoefficients {
                            bulk: StepCoefficients {
                                sources_w,
                                conductances_w_k,
                                payload: LaminarWindkesselExchange {
                                    thermal: payload,
                                    circuit: LaminarWindkesselReport {
                                        duration_s,
                                        depth_m: self.binding.depth_m(),
                                        pressure_offsets_pa,
                                        outward_volume_rates_m2_s,
                                        balances,
                                        maximum_pressure_residual_pa,
                                        maximum_flow_residual_m2_s,
                                    },
                                },
                            },
                            boundaries,
                        },
                    },
                )
            },
        )
    }
}
