//! Spatial mechanical snapshots on the SAME original thermal endpoint.
//!
//! The resolver runs inside every native thermal residual, including numerical
//! derivative probes. Its returned source object supplies all force/boundary
//! values for that momentum solve and is retained at acceptance. This supports
//! prescribed time laws and finite thermal actuator feedback without lagging
//! the drive or introducing a second momentum or thermal solver.

use super::*;
use fs_flux::OpenBoundaryPolicy;

mod mechanics;
pub use mechanics::LaminarMechanicalWork;

mod scheduled;
pub use scheduled::LaminarScheduledCandidate;

/// Absolute endpoint time and authored continuous segment. A closing event
/// belongs to its preceding segment; no floating-time perturbation is needed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminarDriveQuery {
    pub time_s: f64,
    pub segment: usize,
}

/// Original mechanical boundary kinds, with their original SI parameters.
/// Values are evaluated through the SAME retained mechanical snapshot.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LaminarDrivenBoundary {
    Velocity {
        edge: usize,
    },
    NavierSlip {
        edge: usize,
        tangential_resistance_pa_s_per_m: f64,
    },
    Open {
        edge: usize,
        backflow: OpenBoundaryPolicy,
    },
    OpenReservoir {
        edge: usize,
    },
}
impl LaminarDrivenBoundary {
    fn edge(self) -> usize {
        match self {
            Self::Velocity { edge }
            | Self::NavierSlip { edge, .. }
            | Self::Open { edge, .. }
            | Self::OpenReservoir { edge } => edge,
        }
    }
}

/// One coherent, immutable mechanical source. Implementations must be pure:
/// repeated spatial queries on a snapshot return the same finite SI data.
/// Native quadrature, boundary admission and nonlinear residuals remain
/// authoritative. No interpolation, nearest-cell lookup or inferred traction.
///
/// A resolver may own temperature-dependent parameters and borrow original
/// authored source laws. It must not retain a provisional thermal history as
/// externally accepted state. Boundary kinds may change only at authored
/// events; continuous within-segment feedback may change their values.
pub trait LaminarMechanicalDrive {
    fn boundaries(&self) -> &[LaminarDrivenBoundary];
    /// Additional force [N/m³]; excludes the independently requested buoyancy.
    fn force_density(&self, point_m: [f64; 2]) -> [f64; 2];
    /// Prescribed fluid/wall velocity, or carried incoming reservoir velocity.
    fn velocity(&self, edge: usize, point_m: [f64; 2]) -> [f64; 2];
    /// Literal supplied stress [Pa], including its original pressure datum.
    /// Boussinesq callers provide the existing hydrostatically reduced stress.
    fn traction(&self, edge: usize, point_m: [f64; 2], normal: [f64; 2]) -> [f64; 2];
    fn slip_heat_routes(&self) -> &[StokesSlipHeatRoute] {
        &[]
    }
}

/// Original accepted mechanical source, time/segment and thermal source payload.
pub struct LaminarDrivenExchange<D, T> {
    pub query: LaminarDriveQuery,
    pub drive: D,
    pub thermal: T,
}

pub type LaminarDrivenCandidate<'m, D, G, T> =
    LaminarEnthalpyCandidate<'m, G, LaminarDrivenExchange<D, T>>;

impl<'s, 'mesh, 'source, F: FluidTransport + ?Sized> LaminarEnthalpy<'s, 'mesh, 'source, F> {
    /// Resolve original momentum, thermal rheology and mechanical actuator
    /// feedback at one coupled endpoint. The mechanical resolver observes the
    /// actual original caloric candidate, not the start of the time interval.
    /// Both source callbacks share that accepted mechanical snapshot; neither
    /// may advance external physical history. Native fixed-density admission,
    /// sign-split streams, slip routes and final conservation gates are retained.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_driven<'m, D: LaminarMechanicalDrive, G: ThermalFluidTransport, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        query: LaminarDriveQuery,
        history: LaminarEnthalpyHistory<'_>,
        controls: LaminarEnthalpyControls,
        energy_guess_j: &[f64],
        resolve_drive: impl Fn(&Cx<'_>, LaminarDriveQuery, &DepositionSimulation<'m>) -> Result<D, E>,
        resolve_material: impl Fn(MaterialFluxRheologyPoint, f64, &F, &D) -> Result<G, MaterialError>,
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>>,
            &D,
        ) -> Result<LaminarHeatExchange<T>, E>,
    ) -> Result<LaminarDrivenCandidate<'m, D, G, T>, LaminarEnthalpyError<E>> {
        self.advance_driven_with_flow(
            cx,
            thermal,
            query,
            history,
            controls,
            energy_guess_j,
            resolve_drive,
            |input, temperatures, duration, drive| {
                self.resolve_laminar_flow(
                    cx,
                    input,
                    controls,
                    temperatures,
                    duration,
                    &|point, temperature, original| {
                        resolve_material(point, temperature, original, drive)
                    },
                )
                .map(|flow| (flow, None))
            },
            exchange,
        )
    }

    /// The same mechanical feedback with the original SAME-sample bounded
    /// density buoyancy. The source projection, caloric association and reduced
    /// pressure convention are unchanged from advance_boussinesq.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_boussinesq_driven<
        'm,
        D: LaminarMechanicalDrive,
        G: ThermalFluidTransport,
        T,
        E,
    >(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        query: LaminarDriveQuery,
        history: LaminarEnthalpyHistory<'_>,
        config: boussinesq::BoussinesqConfig<2>,
        controls: LaminarEnthalpyControls,
        energy_guess_j: &[f64],
        resolve_drive: impl Fn(&Cx<'_>, LaminarDriveQuery, &DepositionSimulation<'m>) -> Result<D, E>,
        resolve_material: impl Fn(MaterialFluxRheologyPoint, f64, &F, &D) -> Result<G, MaterialError>,
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &MaterialNavierStokesRheologyStep<
                LaminarThermalMaterial<boussinesq::LaminarBoussinesqProjection<G>>,
            >,
            &D,
        ) -> Result<LaminarHeatExchange<T>, E>,
    ) -> Result<
        LaminarDrivenCandidate<'m, D, boussinesq::LaminarBoussinesqProjection<G>, T>,
        LaminarEnthalpyError<E>,
    > {
        checkpoint(cx)?;
        config.validate().map_err(|error| {
            LaminarEnthalpyError::Momentum(MaterialFluxError::Material(MaterialError::State {
                what: error.to_string(),
            }))
        })?;
        self.advance_driven_with_flow(
            cx,
            thermal,
            query,
            history,
            controls,
            energy_guess_j,
            resolve_drive,
            |input, temperatures, duration, drive| {
                self.resolve_boussinesq_flow(
                    cx,
                    input,
                    config,
                    controls,
                    temperatures,
                    duration,
                    &|point, temperature, original| {
                        resolve_material(point, temperature, original, drive)
                    },
                )
            },
            exchange,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn advance_driven_with_flow<'m, D: LaminarMechanicalDrive, G, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        query: LaminarDriveQuery,
        history: LaminarEnthalpyHistory<'_>,
        controls: LaminarEnthalpyControls,
        energy_guess_j: &[f64],
        resolve_drive: impl Fn(&Cx<'_>, LaminarDriveQuery, &DepositionSimulation<'m>) -> Result<D, E>,
        resolve_flow: impl Fn(
            LaminarEnthalpyInput<'_>,
            &[f64],
            f64,
            &D,
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
            &D,
        ) -> Result<LaminarHeatExchange<T>, E>,
    ) -> Result<LaminarDrivenCandidate<'m, D, G, T>, LaminarEnthalpyError<E>> {
        let mut solved = self.advance_with_endpoint(
            cx,
            thermal,
            query.time_s,
            history,
            controls,
            energy_guess_j,
            |candidate, temperatures, duration| {
                let drive =
                    resolve_drive(cx, query, candidate).map_err(LaminarEnthalpyError::Drive)?;
                checkpoint(cx)?;
                if drive.boundaries().len() > self.system.mesh.edges.len() {
                    return Err(LaminarEnthalpyError::Invalid(
                        "mechanical boundary snapshot exceeds original edge capacity",
                    ));
                }
                let (flow, buoyancy, routes) = with_native(&drive, history, |input| {
                    let routes = self.binding.slip_routes(
                        cx,
                        self.system,
                        candidate,
                        input.boundaries,
                        input.slip_heat_routes,
                    )?;
                    let result = resolve_flow(input, temperatures, duration, &drive);
                    checkpoint(cx)?;
                    let (flow, buoyancy) = result.map_err(LaminarEnthalpyError::Momentum)?;
                    Ok::<_, LaminarEnthalpyError<E>>((flow, buoyancy, routes))
                })?;
                let response = exchange(cx, candidate, &flow, &drive)
                    .map_err(LaminarEnthalpyError::Exchange)?;
                checkpoint(cx)?;
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
                    buoyancy,
                    &routes,
                    LaminarHeatExchange {
                        inlet_enthalpies,
                        heat: BoundaryStepCoefficients {
                            bulk: StepCoefficients {
                                sources_w,
                                conductances_w_k,
                                payload: LaminarDrivenExchange {
                                    query,
                                    drive,
                                    thermal: payload,
                                },
                            },
                            boundaries,
                        },
                    },
                )
            },
        )?;
        solved.report.mechanical = Some(mechanics::measure(
            cx,
            self.system,
            &solved.flow.step.solution,
            &solved.report.exchange.drive,
            self.binding.depth_m(),
            query.time_s - thermal.time_s(),
        )?);
        checkpoint(cx)?;
        Ok(solved)
    }
}

// Original callback views borrow the snapshot only during the native call.
// Homogeneous closure vectors avoid one heap-allocated trait object per edge.
fn with_native<D: LaminarMechanicalDrive, R>(
    drive: &D,
    history: LaminarEnthalpyHistory<'_>,
    run: impl FnOnce(LaminarEnthalpyInput<'_>) -> R,
) -> R {
    let velocities: Vec<_> = drive
        .boundaries()
        .iter()
        .map(|boundary| {
            let edge = boundary.edge();
            move |point| drive.velocity(edge, point)
        })
        .collect();
    let tractions: Vec<_> = drive
        .boundaries()
        .iter()
        .map(|boundary| {
            let edge = boundary.edge();
            move |point, normal| drive.traction(edge, point, normal)
        })
        .collect();
    let boundaries: Vec<_> = drive
        .boundaries()
        .iter()
        .enumerate()
        .map(|(i, boundary)| match *boundary {
            LaminarDrivenBoundary::Velocity { edge } => MaterialNavierStokesBoundary::Velocity {
                edge,
                value: &velocities[i],
            },
            LaminarDrivenBoundary::NavierSlip {
                edge,
                tangential_resistance_pa_s_per_m,
            } => MaterialNavierStokesBoundary::NavierSlip {
                edge,
                value: &velocities[i],
                tangential_resistance_pa_s_per_m,
            },
            LaminarDrivenBoundary::Open { edge, backflow } => MaterialNavierStokesBoundary::Open {
                edge,
                traction: &tractions[i],
                backflow,
            },
            LaminarDrivenBoundary::OpenReservoir { edge } => {
                MaterialNavierStokesBoundary::OpenReservoir {
                    edge,
                    traction: &tractions[i],
                    velocity: &velocities[i],
                }
            }
        })
        .collect();
    let force = |point| drive.force_density(point);
    run(LaminarEnthalpyInput {
        history,
        force_density: &force,
        boundaries: &boundaries,
        slip_heat_routes: drive.slip_heat_routes(),
    })
}
