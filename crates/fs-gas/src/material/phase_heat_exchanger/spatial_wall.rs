//! Original phase streams coupled to the existing spatial wall-energy owner.
//!
//! Attach this local source with the native implicit chain/graph's
//! `with_heat_sources`. Reference mass, original h/u storage, nonlinear
//! conductivity, contacts, geometry, time integration and thermal acceptance
//! remain there. This adapter owns no second wall state or time integrator.
//!
//! Each port is an independently supplied, perfectly mixed, negligible-storage
//! passage at fixed feed pressure. Repeating an inlet across cells does not
//! create a serial or shared finite stream. Conductance is the complete local
//! wall-cell-to-mixed-outlet coupling; do not also apply that physical heat as
//! an ambient boundary or authored source. Wall mechanical/storage constraints
//! remain those of the original thermal law; no pressure work is added here.
//!
//! The native thermal audit correctly leaves this source's entropy unclassified.
//! Evaluations retain the actual stream mdot*Delta s, not a fictitious bath at
//! the inlet or outlet temperature. Integrate receipts with the accepted native
//! stage quadrature; one endpoint times a multistep interval is not its receipt.

use super::{
    MixtureHeatExchangerError as StreamError,
    wall::MixtureWallPort,
    wall_stream::{MixtureWallStreamStep, PreparedMixtureWallStream},
};
use fs_material::thermochemical_gas::nonideal::flash::{
    ThermochemicalFlashError,
    equilibrium::acoustics::{PhaseEquilibriumAcousticControls, PhaseEquilibriumAcoustics},
};
use fs_qty::Temperature;
use fs_thermal::enthalpy::{
    DepositionError, ThermalLaw,
    conduction::implicit::{
        ImplicitConductionLaw,
        source::{ImplicitHeatSourceLaw, ImplicitHeatSourceState},
    },
};
use std::cell::{Cell, RefCell};

/// One immutable feed and the explicitly chosen original response neighborhood.
#[derive(Clone, Copy, Debug)]
pub struct MixtureStreamWallPort<'m> {
    pub stream: MixtureWallPort<'m>,
    pub response: PhaseEquilibriumAcousticControls,
}

/// Aggregate work across preparation, attachment, all trials and all advances.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MixtureStreamWallBudget {
    pub maximum_duty_evaluations: usize,
    pub maximum_pt_evaluations: usize,
}

/// Charged work ceilings, not physical heat or committed timestep history.
///
/// Successful queries charge their native counters. A refused stream query
/// retains its whole granted reservation because some native errors do not
/// expose partially consumed PT work. No retry or failed advance refunds it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MixtureStreamWallWork {
    pub charged_duty_evaluations: usize,
    pub charged_pt_evaluations: usize,
    pub refused_stream_queries: usize,
}

/// A typed source refusal; the native trait also retains this in last_failure().
#[derive(Clone, Debug, PartialEq)]
pub enum MixtureStreamWallError {
    Thermal(DepositionError),
    Stream { port: usize, source: StreamError },
}
impl From<DepositionError> for MixtureStreamWallError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl core::fmt::Display for MixtureStreamWallError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial wall-stream source: {self:?}")
    }
}
impl std::error::Error for MixtureStreamWallError {}

/// Instantaneous original source receipts, with no physical history mutation.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureStreamWallEvaluation {
    /// Complete power INTO the wall and its derivative on the original h chart.
    pub source: ImplicitHeatSourceState,
    /// Actual source-integrated enthalpy/entropy for each independent passage.
    pub streams: Vec<MixtureWallStreamStep>,
    /// Native outlet response evidence. None when zero G or a wall latent
    /// plateau makes a fluid temperature derivative unnecessary.
    pub responses: Vec<Option<PhaseEquilibriumAcoustics>>,
}

struct PreparedPort<'m> {
    stream: PreparedMixtureWallStream<'m>,
    declaration: MixtureStreamWallPort<'m>,
}

/// One or two mixed-flow ports on the SAME original wall caloric/conduction law.
///
/// The local native source trait is synchronous and has no Send/Sync bounds.
/// Interior state contains only charged work and refusal diagnostics, never
/// accepted heat, phase history or wall inventory. A new source declares a new
/// work allowance; sharing this instance across calls shares its finite budget.
pub struct MixtureStreamWallSource<'w, 'm> {
    wall: &'w dyn ImplicitConductionLaw,
    ports: Vec<PreparedPort<'m>>,
    budget: MixtureStreamWallBudget,
    work: Cell<MixtureStreamWallWork>,
    failure: RefCell<Option<MixtureStreamWallError>>,
}

impl<'w, 'm> MixtureStreamWallSource<'w, 'm> {
    /// Prepare each original inlet once. Stream policies must already match
    /// the response policy; tighter PH enthalpy accuracy does not change it.
    pub fn new(
        wall: &'w dyn ImplicitConductionLaw,
        ports: &[MixtureStreamWallPort<'m>],
        budget: MixtureStreamWallBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, MixtureStreamWallError> {
        checkpoint()?;
        if !(1..=2).contains(&ports.len())
            || budget.maximum_duty_evaluations == 0
            || budget.maximum_pt_evaluations == 0
        {
            return Err(DepositionError::InvalidInput(
                "one or two wall streams and positive work limits required",
            )
            .into());
        }
        let mut prepared = Vec::with_capacity(ports.len());
        let mut work = MixtureStreamWallWork::default();
        for (port, original) in ports.iter().enumerate() {
            checkpoint()?;
            original
                .response
                .validate()
                .map_err(|source| MixtureStreamWallError::Stream {
                    port,
                    source: super::phase_error(port, work.charged_pt_evaluations, source),
                })?;
            if original.response.phase != original.stream.controls.equilibrium.phase_options()
                || !original.stream.conductance.value().is_finite()
                || original.stream.conductance.value() < 0.0
            {
                return Err(DepositionError::InvalidInput(
                    "wall-stream conductance or original phase policy mismatch",
                )
                .into());
            }
            let remaining_pt = budget.maximum_pt_evaluations - work.charged_pt_evaluations;
            if remaining_pt == 0 {
                return Err(DepositionError::Budget("spatial wall-stream preparation").into());
            }
            let mut declaration = *original;
            declaration.stream.controls.maximum_pt_evaluations = declaration
                .stream
                .controls
                .maximum_pt_evaluations
                .min(remaining_pt);
            let mut cancellation = None;
            let result = PreparedMixtureWallStream::new(
                declaration.stream.inlet,
                declaration.stream.controls,
                || checkpoint_cancelled(&mut checkpoint, &mut cancellation),
            );
            if let Some(error) = cancellation {
                return Err(error.into());
            }
            let stream =
                result.map_err(|source| MixtureStreamWallError::Stream { port, source })?;
            work.charged_pt_evaluations += stream.preparation_pt_evaluations();
            prepared.push(PreparedPort {
                stream,
                declaration,
            });
        }
        checkpoint()?;
        Ok(Self {
            wall,
            ports: prepared,
            budget,
            work: Cell::new(work),
            failure: RefCell::new(None),
        })
    }

    #[must_use]
    pub fn work(&self) -> MixtureStreamWallWork {
        self.work.get()
    }

    /// Last refusal, retained even if a later query succeeds.
    #[must_use]
    pub fn last_failure(&self) -> Option<MixtureStreamWallError> {
        self.failure.borrow().clone()
    }

    /// Query original stream receipts and the complete wall-energy Jacobian.
    ///
    /// On a latent wall plateau, dT_wall/dh=0 exactly: evaluate actual flow
    /// heat but skip every otherwise unnecessary fluid response probe.
    pub fn evaluate(
        &self,
        specific_h_j_kg: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<MixtureStreamWallEvaluation, MixtureStreamWallError> {
        let result = self.evaluate_inner(specific_h_j_kg, &mut checkpoint);
        if let Err(error) = &result {
            self.failure.replace(Some(error.clone()));
        }
        result
    }

    fn evaluate_inner(
        &self,
        h: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<MixtureStreamWallEvaluation, MixtureStreamWallError> {
        checkpoint()?;
        let wall = self.wall.state_at_specific_enthalpy(h)?;
        if !super::positive(wall.temperature_k)
            || !wall.temperature_slope_kg_k_j.is_finite()
            || wall.temperature_slope_kg_k_j < 0.0
        {
            return Err(DepositionError::MaterialDomain.into());
        }
        let mut streams = Vec::with_capacity(self.ports.len());
        let mut responses = Vec::with_capacity(self.ports.len());
        let mut source = ImplicitHeatSourceState {
            power_w: 0.0,
            power_enthalpy_derivative_w_kg_j: 0.0,
        };
        for (port, prepared) in self.ports.iter().enumerate() {
            checkpoint()?;
            let mut work = self.work.get();
            let limits = prepared.declaration.stream.controls;
            let duty = limits
                .maximum_duty_evaluations
                .min(self.budget.maximum_duty_evaluations - work.charged_duty_evaluations);
            let pt = limits
                .maximum_pt_evaluations
                .min(self.budget.maximum_pt_evaluations - work.charged_pt_evaluations);
            if duty == 0 {
                return Err(DepositionError::Budget("spatial wall-stream duty evaluations").into());
            }
            // Reserve before any user callback. Success returns unused allowance;
            // opaque failure keeps this finite ceiling across every native retry.
            work.charged_duty_evaluations += duty;
            work.charged_pt_evaluations += pt;
            self.work.set(work);
            let mut cancellation = None;
            let temperature = Temperature::new(wall.temperature_k);
            let conductance = prepared.declaration.stream.conductance;
            let result = if wall.temperature_slope_kg_k_j == 0.0 {
                prepared
                    .stream
                    .evaluate(temperature, conductance, duty, pt, || {
                        checkpoint_cancelled(checkpoint, &mut cancellation)
                    })
                    .map(|stream| {
                        let used_pt = stream.pt_evaluations;
                        (stream, None, 0.0, used_pt)
                    })
            } else {
                prepared
                    .stream
                    .evaluate_response(
                        temperature,
                        conductance,
                        prepared.declaration.response,
                        duty,
                        pt,
                        || checkpoint_cancelled(checkpoint, &mut cancellation),
                    )
                    .map(|response| {
                        let used_pt = response.pt_evaluations();
                        let tangent = -response.heat_power_temperature_derivative_w_k
                            * wall.temperature_slope_kg_k_j;
                        (response.stream, response.outlet_response, tangent, used_pt)
                    })
            };
            let result = match cancellation {
                Some(error) => Err(MixtureStreamWallError::Thermal(error)),
                None => result.map_err(|source| MixtureStreamWallError::Stream { port, source }),
            };
            let mut work = self.work.get();
            let (stream, response, tangent, used_pt) = match result {
                Ok(value) => value,
                Err(error) => {
                    work.refused_stream_queries += 1;
                    self.work.set(work);
                    return Err(error);
                }
            };
            work.charged_duty_evaluations -= duty - stream.duty_evaluations;
            work.charged_pt_evaluations -= pt - used_pt;
            self.work.set(work);
            source.power_w -= stream.enthalpy_flow_change_w;
            source.power_enthalpy_derivative_w_kg_j += tangent;
            streams.push(stream);
            responses.push(response);
        }
        if !source.power_w.is_finite() || !source.power_enthalpy_derivative_w_kg_j.is_finite() {
            return Err(DepositionError::Arithmetic.into());
        }
        checkpoint()?;
        Ok(MixtureStreamWallEvaluation {
            source,
            streams,
            responses,
        })
    }
}

impl ImplicitHeatSourceLaw for MixtureStreamWallSource<'_, '_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.wall.caloric_law()
    }

    fn state_at_specific_enthalpy(
        &self,
        h: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitHeatSourceState, DepositionError> {
        self.evaluate(h, checkpoint)
            .map(|value| value.source)
            .map_err(|error| match error {
                MixtureStreamWallError::Thermal(error) => error,
                MixtureStreamWallError::Stream {
                    source: StreamError::Cancelled,
                    ..
                } => DepositionError::Cancelled,
                MixtureStreamWallError::Stream {
                    source: StreamError::WorkLimit { .. },
                    ..
                }
                | MixtureStreamWallError::Stream {
                    source:
                        StreamError::Equilibrium {
                            source: ThermochemicalFlashError::WorkLimit { .. },
                            ..
                        },
                    ..
                } => DepositionError::Budget("spatial wall-stream phase work"),
                MixtureStreamWallError::Stream { .. } => DepositionError::MaterialDomain,
            })
    }
}

fn checkpoint_cancelled(
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    failure: &mut Option<DepositionError>,
) -> bool {
    if failure.is_none() {
        *failure = checkpoint().err();
    }
    failure.is_some()
}
