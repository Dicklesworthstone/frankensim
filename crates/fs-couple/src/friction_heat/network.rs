//! Simultaneous fixed interfaces on one original finite enthalpy network.
//!
//! Each physical contact has its own massless common-temperature coordinate.
//! All interface powers enter one backward-Euler residual, including when
//! several contacts share a finite cell. No intermediate thermal state is
//! accepted, no capacitance is added and no heat share is prescribed.

use super::{FrictionHeatRoute, FrictionHeatTolerances, finite, invalid};
use fs_thermal::enthalpy::backward_euler::{
    BackwardEulerStep, BackwardEulerTrial, StepCoefficients,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, EnergyAudit};
use std::cell::RefCell;

/// Heat and original material data for one interface in canonical link order.
#[derive(Debug, Clone, Copy)]
pub struct InterfaceHeatInput {
    pub dissipated_j: f64,
    pub conductivity_w_m_k: [f64; 2],
}

/// One joint mechanics/material evaluation, before any owner publishes state.
pub struct FrictionHeatNetworkInput<P> {
    /// Exactly one item per `link_ids()`, in that order.
    pub interfaces: Vec<InterfaceHeatInput>,
    /// Other heat only. Bound interface conductances must not be overridden.
    pub transport: StepCoefficients<P>,
}

#[derive(Debug, Clone)]
pub struct InterfaceHeatReport {
    pub link_id: u64,
    pub cells: [CellId; 2],
    pub interface_temperature_k: f64,
    pub conductances_w_k: [f64; 2],
    /// Signed net heat, including ordinary hot-to-cold exchange.
    pub body_heat_j: [f64; 2],
    pub dissipated_j: f64,
    pub interface_residual_j: f64,
}

#[derive(Debug, Clone)]
pub struct FrictionHeatNetworkReport {
    /// Canonical ascending original link ID, independently audited per contact.
    pub interfaces: Vec<InterfaceHeatReport>,
    pub dissipated_j: f64,
    pub external_heat_j: f64,
    pub thermal: EnergyAudit,
    pub total_energy_residual_j: f64,
}

/// N original cell energy increments and M distinct common temperatures.
pub struct FrictionHeatNetworkStep<'s, 'm> {
    thermal: BackwardEulerStep<'s, 'm>,
    routes: Vec<&'s FrictionHeatRoute>,
    link_ids: Vec<u64>,
}
impl<'s, 'm> FrictionHeatNetworkStep<'s, 'm> {
    pub fn new(
        thermal: &'s DepositionSimulation<'m>,
        target_s: f64,
        routes: &[&'s FrictionHeatRoute],
        maximum_interfaces: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if routes.is_empty() || routes.len() > maximum_interfaces {
            return Err(invalid(
                "bounded nonempty friction heat interfaces required",
            ));
        }
        let mut routes = routes.to_vec();
        routes.sort_by_key(|route| route.patch.link.id);
        for (i, route) in routes.iter().enumerate() {
            checkpoint()?;
            route.validate(thermal)?;
            if i > 0 && routes[i - 1].patch.link.id == route.patch.link.id {
                return Err(invalid("one heat owner per original contact link required"));
            }
        }
        let link_ids = routes.iter().map(|route| route.patch.link.id).collect();
        Ok(Self {
            thermal: BackwardEulerStep::new(thermal, target_s, checkpoint)?,
            routes,
            link_ids,
        })
    }
    pub fn cell_ids(&self) -> &[CellId] {
        self.thermal.cell_ids()
    }
    pub fn link_ids(&self) -> &[u64] {
        &self.link_ids
    }
    pub fn duration_s(&self) -> f64 {
        self.thermal.duration_s()
    }

    /// Resolve all laws and transport at one candidate thermal endpoint.
    /// Temperatures and returned interface inputs follow `link_ids()`.
    pub fn evaluate<P, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        interface_temperatures_k: &[f64],
        resolve: impl FnOnce(
            &DepositionSimulation<'m>,
            &[f64],
        ) -> Result<FrictionHeatNetworkInput<P>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<FrictionHeatNetworkTrial<'m, P>, E> {
        if interface_temperatures_k.len() != self.routes.len()
            || interface_temperatures_k
                .iter()
                .any(|t| !t.is_finite() || *t <= 0.0)
        {
            return Err(invalid("one positive Kelvin coordinate per interface required").into());
        }
        let checkpoint = RefCell::new(checkpoint);
        let check = || (*checkpoint.borrow_mut())();
        let dt = self.duration_s();
        let trial = self.thermal.evaluate(
            energy_increments_j,
            |state| -> Result<StepCoefficients<Payload<P>>, E> {
                let FrictionHeatNetworkInput {
                    interfaces,
                    transport,
                } = resolve(state, interface_temperatures_k)?;
                if interfaces.len() != self.routes.len() {
                    return Err(invalid("incomplete friction heat interface inputs").into());
                }
                let StepCoefficients {
                    mut sources_w,
                    mut conductances_w_k,
                    payload,
                } = transport;
                // Suppress all original series links together. Surface coverage
                // stays in the original network; no exposed face is invented.
                for id in &self.link_ids {
                    check()?;
                    if conductances_w_k.insert(*id, 0.0).is_some() {
                        return Err(invalid("friction interface conduction has two owners").into());
                    }
                }
                let external_heat_j = finite(
                    dt * state.external_step_power_w(&sources_w, &conductances_w_k, &check)?,
                )?;
                let mut reports = Vec::with_capacity(self.routes.len());
                let mut dissipated_j = 0.0;
                for (i, route) in self.routes.iter().enumerate() {
                    check()?;
                    let input = interfaces[i];
                    if !input.dissipated_j.is_finite() || input.dissipated_j < 0.0 {
                        return Err(invalid(
                            "only nonnegative physical friction heat may enter cells",
                        )
                        .into());
                    }
                    let conductances = route.conductances(input.conductivity_w_m_k)?;
                    let mut body_heat_j = [0.0; 2];
                    for (side, id) in route.cells().into_iter().enumerate() {
                        let power = finite(
                            conductances[side]
                                * (interface_temperatures_k[i] - state.temperature_k(id)?),
                        )?;
                        body_heat_j[side] = finite(dt * power)?;
                        if power != 0.0 && body_heat_j[side] == 0.0 {
                            return Err(DepositionError::Arithmetic.into());
                        }
                        let source = sources_w.entry(id).or_insert(0.0);
                        *source = finite(*source + power)?;
                    }
                    reports.push(InterfaceHeatReport {
                        link_id: route.patch.link.id,
                        cells: route.cells(),
                        interface_temperature_k: interface_temperatures_k[i],
                        conductances_w_k: conductances,
                        body_heat_j,
                        dissipated_j: input.dissipated_j,
                        interface_residual_j: finite(
                            finite(body_heat_j[0] + body_heat_j[1])? - input.dissipated_j,
                        )?,
                    });
                    dissipated_j = finite(dissipated_j + input.dissipated_j)?;
                }
                Ok(StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload: Payload {
                        source: payload,
                        reports,
                        dissipated_j,
                        external_heat_j,
                    },
                })
            },
            &check,
        )?;
        let payload = trial.payload();
        let report = FrictionHeatNetworkReport {
            interfaces: payload.reports.clone(),
            dissipated_j: payload.dissipated_j,
            external_heat_j: payload.external_heat_j,
            thermal: *trial.audit(),
            total_energy_residual_j: finite(
                trial.audit().stored_energy_change_j
                    - payload.external_heat_j
                    - payload.dissipated_j,
            )?,
        };
        let mut residual = trial.residual_j().to_vec();
        residual.extend(report.interfaces.iter().map(|r| r.interface_residual_j));
        Ok(FrictionHeatNetworkTrial {
            trial,
            residual,
            report,
        })
    }
}

struct Payload<P> {
    source: P,
    reports: Vec<InterfaceHeatReport>,
    dissipated_j: f64,
    external_heat_j: f64,
}

pub struct FrictionHeatNetworkTrial<'m, P> {
    trial: BackwardEulerTrial<'m, Payload<P>>,
    residual: Vec<f64>,
    report: FrictionHeatNetworkReport,
}
impl<'m, P> FrictionHeatNetworkTrial<'m, P> {
    pub fn residual_j(&self) -> &[f64] {
        &self.residual
    }
    pub fn payload(&self) -> &P {
        &self.trial.payload().source
    }
    pub fn report(&self) -> &FrictionHeatNetworkReport {
        &self.report
    }
    pub fn accept(
        self,
        tolerance: FrictionHeatTolerances,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, P, FrictionHeatNetworkReport), DepositionError> {
        checkpoint()?;
        if [
            tolerance.cell_j,
            tolerance.thermal_j,
            tolerance.interface_j,
            tolerance.total_j,
        ]
        .iter()
        .any(|t| !t.is_finite() || *t < 0.0)
        {
            return Err(invalid(
                "finite nonnegative friction heat tolerances required",
            ));
        }
        for interface in &self.report.interfaces {
            checkpoint()?;
            if interface.interface_residual_j.abs() > tolerance.interface_j {
                return Err(DepositionError::EnergyBalance(
                    interface.interface_residual_j,
                ));
            }
        }
        if self.report.total_energy_residual_j.abs() > tolerance.total_j {
            return Err(DepositionError::EnergyBalance(
                self.report.total_energy_residual_j,
            ));
        }
        let (thermal, payload, _) =
            self.trial
                .accept(tolerance.cell_j, tolerance.thermal_j, &mut checkpoint)?;
        checkpoint()?;
        Ok((thermal, payload.source, self.report))
    }
}
