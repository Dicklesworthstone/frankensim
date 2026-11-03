//! Endpoint material closure on the original unsteady SIMPLEC/heat owner.
//!
//! Only the nonlinear iterate changes inside a coupling pass. The native BE/
//! BDF2 previous and older velocity levels and the accepted old temperature
//! remain fixed until the enclosing march accepts one physical step.

use super::{
    ChtError, EnergyConfig, FluidProperties, Marcher, MaterialCoefficients, MaterialCoupling,
    MaterialHeatStep, RadiationConfig, Radiators, SolidMaterial, ThermalSetup, TransientRecord,
    UnsteadyConfig, finite, poll,
};
use crate::conjugate::turbulence::TURBULENT_PRANDTL;
#[cfg(feature = "enthalpy-coupling")]
use super::MaterialViscousHeat;
use core::cell::Cell;
use fs_exec::CancelGate;

/// Explicit nonlinear and whole-march work bounds for implicit material feedback.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialEndpointConfig {
    /// Maximum complete native momentum/heat pairs at each physical time level.
    pub maximum_passes_per_step: usize,
    /// All native momentum solves, including discarded coupling trials.
    pub maximum_momentum_solves: usize,
    /// Actual linear energy solves across all material, caloric, compact-storage
    /// and radiation trials. Source-caloric controls may impose a tighter bound.
    pub maximum_energy_solves: usize,
    /// All in-domain transport-source epochs, including the opening state,
    /// endpoint checks and separately queried relaxed trial temperatures.
    /// Caloric queries retain their original independent heat-solve budgets.
    pub maximum_source_refreshes: usize,
    /// Maximum unrelaxed output-minus-trial temperature difference, K.
    pub temperature_tolerance_k: f64,
    /// Relative difference of applied and actual endpoint mu and k.
    pub transport_relative_tolerance: f64,
    /// Maximum componentwise endpoint-minus-applied buoyancy force, N/m3.
    pub force_tolerance_n_m3: f64,
    /// Relax only the next temperature iterate, never a convergence defect,
    /// physical velocity, heat record or accepted time level.
    pub temperature_relaxation: f64,
}

/// Closure defects at one accepted original physical time level.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialEndpointStepReport {
    pub time_s: f64,
    pub passes: usize,
    pub temperature_defect_k: f64,
    pub transport_relative_defect: f64,
    pub force_defect_n_m3: f64,
}

/// The original native result and the work spent closing its material feedback.
/// Discarded nonlinear trials never enter its physical histories or averages.
pub struct MaterialEndpointMarch<M> {
    pub solution: M,
    pub config: MaterialEndpointConfig,
    pub records: Vec<MaterialEndpointStepReport>,
    pub momentum_solves: usize,
    pub energy_solves: usize,
    pub source_refreshes: usize,
}

pub(crate) struct MaterialEndpointCoupling {
    config: MaterialEndpointConfig,
    records: Vec<MaterialEndpointStepReport>,
    momentum_solves: usize,
    remaining_energy_solves: Cell<usize>,
    source_refreshes: usize,
    // Reuse both material epochs and the temperature iterate across steps.
    next_fields: Option<MaterialCoefficients>,
    iterate: Vec<f64>,
}

fn invalid(reason: impl Into<String>) -> ChtError {
    ChtError::InvalidInput {
        field: "unsteady.material_endpoint",
        reason: reason.into(),
    }
}

impl MaterialEndpointCoupling {
    pub(crate) fn new(config: MaterialEndpointConfig) -> Result<Self, ChtError> {
        if !(1..=1024).contains(&config.maximum_passes_per_step)
            || !(1..=1_000_000).contains(&config.maximum_momentum_solves)
            || !(1..=1_000_000).contains(&config.maximum_energy_solves)
            || !(1..=2_000_001).contains(&config.maximum_source_refreshes)
            || [
                config.temperature_tolerance_k,
                config.transport_relative_tolerance,
                config.force_tolerance_n_m3,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
            || config.transport_relative_tolerance >= 1.0
            || !config.temperature_relaxation.is_finite()
            || config.temperature_relaxation <= 0.0
            || config.temperature_relaxation > 1.0
        {
            return Err(invalid(
                "bounded positive work, finite nonnegative tolerances, relative tolerance below one and relaxation in (0,1] required",
            ));
        }
        Ok(Self {
            config,
            records: Vec::new(),
            momentum_solves: 0,
            source_refreshes: 0,
            remaining_energy_solves: Cell::new(config.maximum_energy_solves),
            next_fields: None,
            iterate: Vec::new(),
        })
    }

    pub(super) fn admit(&self, unsteady: &UnsteadyConfig) -> Result<(), ChtError> {
        if self.config.maximum_momentum_solves < unsteady.steps
            || self.config.maximum_energy_solves < unsteady.steps
            || self.config.maximum_source_refreshes <= unsteady.steps
        {
            return Err(invalid(
                "whole-march budgets must cover every native time step and its source endpoint, plus the opening source epoch",
            ));
        }
        Ok(())
    }

    fn charge_source(&mut self) -> Result<(), ChtError> {
        if self.source_refreshes == self.config.maximum_source_refreshes {
            return Err(invalid(format!(
                "transport-source refresh budget exhausted after {} of {} epochs",
                self.source_refreshes, self.config.maximum_source_refreshes,
            )));
        }
        self.source_refreshes += 1;
        Ok(())
    }

    pub(super) fn refresh<E: From<ChtError>>(
        &mut self,
        coupling: &mut MaterialCoupling<'_, E>,
        temperature: &[f64],
    ) -> Result<(), E> {
        self.charge_source()?;
        (coupling.refresh)(temperature, coupling.fields)
    }

    pub(crate) fn finish<M>(self, solution: M) -> MaterialEndpointMarch<M> {
        MaterialEndpointMarch {
            solution,
            config: self.config,
            records: self.records,
            momentum_solves: self.momentum_solves,
            source_refreshes: self.source_refreshes,
            energy_solves: self.config.maximum_energy_solves - self.remaining_energy_solves.get(),
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn step<E: From<ChtError>>(
        &mut self,
        marcher: &mut Marcher<'_>,
        coupling: &mut MaterialCoupling<'_, E>,
        material_heat: &mut Option<
            &mut dyn FnMut(
                &MaterialHeatStep<'_>,
                &mut MaterialCoefficients,
            ) -> Result<(Vec<f64>, TransientRecord), E>,
        >,
        #[cfg(feature = "enthalpy-coupling")] mut viscous_heat: Option<&mut MaterialViscousHeat>,
        fluid: &FluidProperties,
        solids: &[SolidMaterial],
        thermal: &mut ThermalSetup,
        energy: &EnergyConfig,
        capacity: &[f64],
        storage: &[f64],
        previous_temperature: &[f64],
        step: usize,
        inflow: f64,
        solids_present: bool,
        radiation: Option<(&Radiators, &RadiationConfig)>,
        gate: &CancelGate,
    ) -> Result<(Vec<f64>, TransientRecord), E> {
        let domain = marcher.domain;
        let dt = marcher.unsteady.time_step_s;
        let time_s = step as f64 * dt;
        let mut iterate = core::mem::take(&mut self.iterate);
        iterate.clear();
        iterate.extend_from_slice(previous_temperature);
        let mut endpoint = self
            .next_fields
            .take()
            .unwrap_or_else(|| coupling.fields.clone());
        endpoint
            .viscosity_pa_s
            .clone_from(&coupling.fields.viscosity_pa_s);
        endpoint
            .conductivity_w_m_k
            .clone_from(&coupling.fields.conductivity_w_m_k);
        endpoint.force_n_m3.clone_from(&coupling.fields.force_n_m3);
        let mut iterations = 0usize;
        #[cfg(feature = "enthalpy-coupling")]
        let authored_power = if viscous_heat.is_some() {
            thermal.power_w.clone()
        } else {
            Vec::new()
        };
        for pass in 1..=self.config.maximum_passes_per_step {
            poll(gate)?;
            if self.momentum_solves == self.config.maximum_momentum_solves {
                return Err(invalid(format!(
                    "native momentum budget exhausted after {} of {} solves",
                    self.momentum_solves, self.config.maximum_momentum_solves,
                ))
                .into());
            }
            marcher.solver.set_material_transport(
                &coupling.fields.viscosity_pa_s,
                Some(&coupling.fields.force_n_m3),
                gate,
            )?;
            self.momentum_solves += 1;
            // The solver's previous/older time levels NEVER change here.
            // Current velocity/pressure are nonlinear guesses only; warming
            // their solve does not integrate another interval of momentum.
            let inner = marcher.solve_step(step, inflow, gate)?;
            iterations = iterations
                .checked_add(inner)
                .ok_or_else(|| invalid("native momentum iteration count overflow"))?;
            let field = marcher.solver.field();
            if marcher.solver.has_sgs() {
                // The just-projected velocity owns this trial's eddy field.
                // Replace it from native nu_t, never compound it into source k
                // or reuse a preceding coupling pass's thermal transport.
                thermal
                    .eddy_conductivity_w_m_k
                    .resize(domain.cell_count(), 0.0);
                for (cell, (conductivity, nu_t)) in thermal
                    .eddy_conductivity_w_m_k
                    .iter_mut()
                    .zip(marcher.solver.heat_eddy_viscosity())
                    .enumerate()
                {
                    if cell.is_multiple_of(256) {
                        poll(gate)?;
                    }
                    *conductivity = fluid.volumetric_heat_capacity() * nu_t / TURBULENT_PRANDTL;
                    finite("thermal.sgs_conductivity_w_m_k", *conductivity)?;
                }
            }
            #[cfg(feature = "enthalpy-coupling")]
            let mechanical_step = viscous_heat.as_deref().map(|heating| {
                marcher.solver.unsteady_mechanical_step(
                    (step - 1) as f64 * dt, time_s,
                    Some(&coupling.fields.force_n_m3), heating.controls, gate,
                )
            }).transpose()?;
            let physical_fixed_heat_j = 0.0;
            #[cfg(feature = "enthalpy-coupling")]
            let physical_fixed_heat_j = if let Some(receipt) = &mechanical_step {
                MaterialViscousHeat::deposit(thermal, &authored_power, receipt, gate)?
            } else {
                physical_fixed_heat_j
            };
            let interval = MaterialHeatStep {
                domain,
                setup: thermal,
                previous_temperature_k: previous_temperature,
                reference_capacity_j_k: capacity,
                time_s,
                time_step_s: dt,
                physical_fixed_heat_j,
                fluid,
                solids,
                flow: &field,
                energy,
                solids_present,
                radiation,
                gate,
                endpoint_energy_budget: Some(&self.remaining_energy_solves),
            };
            let (next, heat) = if let Some(advance) = material_heat.as_deref_mut() {
                advance(&interval, coupling.fields)?
            } else {
                interval.solve(capacity, storage, &coupling.fields.conductivity_w_m_k)?
            };
            // Source refusal at the actual unrelaxed endpoint is terminal.
            // A small temperature update cannot conceal a sharp mu/k/rho change.
            self.charge_source()?;
            // An inner source-solid owner may have closed its own k(T) here.
            // Preserve that accepted solid epoch; the fluid source refresh
            // below overwrites only its own fluid coefficients.
            endpoint
                .conductivity_w_m_k
                .clone_from(&coupling.fields.conductivity_w_m_k);
            (coupling.refresh)(&next, &mut endpoint)?;
            // Physical SGS applicability is independent of loose material
            // convergence tolerances. Observe the already projected nu_t at
            // actual endpoint mu without replacing the applied source epoch.
            marcher
                .solver
                .admit_sgs_material_endpoint(&endpoint.viscosity_pa_s, gate)?;
            let mut report = MaterialEndpointStepReport {
                time_s,
                passes: pass,
                temperature_defect_k: 0.0,
                transport_relative_defect: 0.0,
                force_defect_n_m3: 0.0,
            };
            if next.len() != iterate.len() {
                return Err(invalid("native thermal candidate changed cell count").into());
            }
            for (cell, (&old, &new)) in iterate.iter().zip(&next).enumerate() {
                if cell.is_multiple_of(256) {
                    poll(gate)?;
                }
                let difference = (new - old).abs();
                finite("unsteady.material_temperature_defect_k", difference)?;
                report.temperature_defect_k = report.temperature_defect_k.max(difference);
            }
            let (transport, force) = field_defects(coupling.fields, &endpoint, gate)?;
            report.transport_relative_defect = transport;
            report.force_defect_n_m3 = force;
            let material_closed = report.temperature_defect_k <= self.config.temperature_tolerance_k
                && transport <= self.config.transport_relative_tolerance
                && force <= self.config.force_tolerance_n_m3;
            let physical_ratio = 0.0;
            #[cfg(feature = "enthalpy-coupling")]
            let (mechanical_endpoint, physical_ratio) = if let (Some(heating), Some(applied)) =
                (viscous_heat.as_deref(), mechanical_step.as_ref())
            {
                // Recheck this SAME projected trajectory with actual endpoint
                // mu/force. Histories and velocities are untouched; only the
                // original constitutive rows and physical power are sampled.
                marcher.solver.set_material_transport(
                    &endpoint.viscosity_pa_s, Some(&endpoint.force_n_m3), gate,
                )?;
                let actual = marcher.solver.sample_unsteady_mechanical_step(
                    (step - 1) as f64 * dt, time_s,
                    Some(&endpoint.force_n_m3), heating.controls, gate,
                )?;
                let mut transfer_defect_w = 0.0;
                for (cell, (&credit, &debit)) in applied.power.cell_heat_w.iter()
                    .zip(&actual.power.cell_heat_w).enumerate()
                {
                    if cell.is_multiple_of(256) { poll(gate)?; }
                    // Non-cancelling local credit/debit differences bound both
                    // each cell and the total thermal/mechanical transfer.
                    transfer_defect_w += (credit - debit).abs();
                }
                finite("thermal.viscous_heat.endpoint_defect_w", transfer_defect_w)?;
                let allowance_j = dt * heating.controls.power_tolerance_w;
                finite("thermal.viscous_heat.endpoint_allowance_j", allowance_j)?;
                let mut defect = ratio(actual.residual_j.abs(), allowance_j)
                    .max(ratio(transfer_defect_w, heating.controls.power_tolerance_w));
                if material_closed && defect <= 1.0 {
                    // A global power equality alone can hide local momentum
                    // defects, including a new force orthogonal to velocity.
                    let momentum = marcher.solver.unsteady_material_momentum_residual(gate)?;
                    defect = defect.max(ratio(momentum, marcher.unsteady.inner_tolerance));
                    marcher.residuals.1 = momentum;
                }
                (Some(actual), defect)
            } else {
                (None, physical_ratio)
            };
            if material_closed && physical_ratio <= 1.0 {
                poll(gate)?;
                #[cfg(feature = "enthalpy-coupling")]
                if let (Some(heating), Some(receipt)) =
                    (viscous_heat.as_deref_mut(), mechanical_endpoint)
                {
                    receipt.admit_balance(heating.controls)?;
                    heating.audit.accept(&receipt, heating.controls)?;
                    heating.records.push(receipt);
                }
                marcher.record_step(step, iterations, gate)?;
                self.records.push(report);
                core::mem::swap(coupling.fields, &mut endpoint);
                self.next_fields = Some(endpoint);
                self.iterate = iterate;
                return Ok((next, heat));
            }
            if pass == self.config.maximum_passes_per_step {
                return Err(ChtError::SolverNotConverged {
                    system: "implicit transient material coupling",
                    iterations: pass,
                    relative_residual: ratio(
                        report.temperature_defect_k,
                        self.config.temperature_tolerance_k,
                    )
                    .max(ratio(transport, self.config.transport_relative_tolerance))
                    .max(ratio(force, self.config.force_tolerance_n_m3))
                    .max(physical_ratio),
                    tolerance: 1.0,
                }
                .into());
            }
            if self.config.temperature_relaxation == 1.0 {
                iterate.clone_from(&next);
                core::mem::swap(coupling.fields, &mut endpoint);
            } else {
                for (trial, &value) in iterate.iter_mut().zip(&next) {
                    *trial = self
                        .config
                        .temperature_relaxation
                        .mul_add(value - *trial, *trial);
                }
                self.refresh(coupling, &iterate)?;
            }
        }
        unreachable!("bounded coupling accepts or refuses its final pass")
    }
}

fn field_defects(
    applied: &MaterialCoefficients,
    endpoint: &MaterialCoefficients,
    gate: &CancelGate,
) -> Result<(f64, f64), ChtError> {
    if applied.viscosity_pa_s.len() != endpoint.viscosity_pa_s.len()
        || applied.conductivity_w_m_k.len() != endpoint.conductivity_w_m_k.len()
        || applied.force_n_m3.len() != endpoint.force_n_m3.len()
    {
        return Err(invalid("material refresh changed the original field shape"));
    }
    let mut transport = 0.0_f64;
    let mut force = 0.0_f64;
    for (i, (&a, &b)) in applied
        .viscosity_pa_s
        .iter()
        .zip(&endpoint.viscosity_pa_s)
        .chain(
            applied
                .conductivity_w_m_k
                .iter()
                .flatten()
                .zip(endpoint.conductivity_w_m_k.iter().flatten()),
        )
        .enumerate()
    {
        if i.is_multiple_of(256) {
            poll(gate)?;
        }
        let relative = if a == b {
            0.0
        } else {
            (b - a).abs() / a.abs().max(b.abs())
        };
        finite("unsteady.material_transport_defect", relative)?;
        transport = transport.max(relative);
    }
    for (i, (&a, &b)) in applied
        .force_n_m3
        .iter()
        .flatten()
        .zip(endpoint.force_n_m3.iter().flatten())
        .enumerate()
    {
        if i.is_multiple_of(256) {
            poll(gate)?;
        }
        let difference = (b - a).abs();
        finite("unsteady.material_force_defect_n_m3", difference)?;
        force = force.max(difference);
    }
    Ok((transport, force))
}

fn ratio(difference: f64, tolerance: f64) -> f64 {
    if difference == 0.0 {
        0.0
    } else {
        difference / tolerance
    }
}
