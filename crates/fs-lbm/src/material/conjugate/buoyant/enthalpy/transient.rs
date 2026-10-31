//! Source-caloric storage and transport on the original native CHT clock.

use super::super::{
    MaterialBuoyancyError, MaterialBuoyancyMarch, MaterialChtFluid, admit_open_boundary_materials,
    poll,
};
use super::{
    CaloricCoupling, ConvectionCaloricLaw, MaterialEnthalpyConfig, MaterialEnthalpyReport,
    MaterialEnthalpyState, defect_field,
};
use crate::conjugate::unsteady::{MaterialCoefficients, MaterialHeatStep};
use crate::conjugate::{
    AdvectiveEnthalpyField, ChtError, EnergyConfig, RadiationConfig, SimpleConfig, SolidMaterial,
    ThermalSetup, TransientRecord, UnsteadyConfig, VoxelDomain,
};
use crate::material::transport::boussinesq::BoussinesqConfig;
use fs_exec::CancelGate;
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState},
};

/// Finite nonlinear and work controls for the native source-caloric march.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialEnthalpyTransientConfig {
    /// Common datum and local advective face/cell heat tolerances.
    pub caloric: MaterialEnthalpyConfig,
    /// Maximum caloric linearizations for each accepted thermal time step.
    pub maximum_passes_per_step: usize,
    /// Whole-march ceiling on actual linear energy solves, including every
    /// endpoint-radiation trial inside each caloric pass.
    pub maximum_energy_solves: usize,
    /// Maximum UNRELAXED output-minus-trial temperature change, K.
    pub temperature_tolerance_k: f64,
    /// Maximum per-cell true-minus-linearized storage rate, W.
    pub maximum_storage_defect_w: f64,
    /// Maximum absolute physical heat-ledger closure for each time step, J.
    pub maximum_energy_closure_j: f64,
    /// Next-trial relaxation in (0,1]; never scales an acceptance defect.
    pub temperature_relaxation: f64,
}

/// Actual nonlinear closure and work at one accepted native time level.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialEnthalpyStepReport {
    /// Accepted physical time, s.
    pub time_s: f64,
    /// Caloric linearizations completed during this time step.
    pub passes: usize,
    /// Actual linear heat solves, including endpoint radiation trials.
    pub energy_solves: usize,
    /// Unrelaxed final temperature defect, K.
    pub temperature_defect_k: f64,
    /// Largest per-cell true-minus-linearized storage rate, W.
    pub maximum_storage_defect_w: f64,
    /// Independent local source-enthalpy face/cell heat defects.
    pub advection: MaterialEnthalpyReport,
}

/// Original momentum/material/thermal histories with full source-law storage.
pub struct MaterialEnthalpyMarch<'a, F, C: ConvectionCaloricLaw> {
    /// Native heat records contain actual `rho0 V delta(h)` and corrected
    /// source-enthalpy boundary/reservoir heat, rather than tangent storage.
    pub convection: MaterialBuoyancyMarch<F>,
    /// Same immutable single-phase law that supplied every storage/flux query.
    pub caloric: &'a C,
    /// Applied finite work budgets and nonlinear acceptance tolerances.
    pub config: MaterialEnthalpyTransientConfig,
    /// Final caloric linearization and true endpoint. These are nonlinear
    /// thermal epochs; transport snapshots retain the explicit material lag.
    pub enthalpy: MaterialEnthalpyState,
    /// One nonlinear acceptance/work record per native time step.
    pub records: Vec<MaterialEnthalpyStepReport>,
}

impl MaterialEnthalpyTransientConfig {
    fn validate(self) -> Result<(), MaterialBuoyancyError> {
        if self.maximum_passes_per_step == 0
            || self.maximum_energy_solves == 0
            || [
                self.temperature_tolerance_k,
                self.maximum_storage_defect_w,
                self.maximum_energy_closure_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
            || !self.temperature_relaxation.is_finite()
            || self.temperature_relaxation <= 0.0
            || self.temperature_relaxation > 1.0
        {
            return Err(MaterialBuoyancyError::Invalid(
                "source-caloric transient requires positive work budgets, finite nonnegative tolerances and relaxation in (0,1]",
            ));
        }
        Ok(())
    }
}

/// Exact source increments use the old accepted temperature as a numerical
/// datum. Subtracting two large absolute/formation enthalpies is unnecessary.
#[allow(clippy::too_many_arguments)]
fn stored_change(
    caloric: CaloricCoupling<'_>,
    domain: &VoxelDomain,
    reference: &FluidTransportState,
    bounds: BoussinesqConfig<3>,
    step: &MaterialHeatStep<'_>,
    temperature: &[f64],
    stored: &mut Vec<f64>,
    gate: &CancelGate,
) -> Result<(), MaterialBuoyancyError> {
    stored.clear();
    for (cell, (&old, &new)) in step
        .previous_temperature_k
        .iter()
        .zip(temperature)
        .enumerate()
    {
        poll(gate)?;
        if !new.is_finite() || new <= 0.0 {
            return Err(MaterialBuoyancyError::Invalid(
                "caloric endpoint must be positive finite kelvin",
            ));
        }
        let delta = new - old;
        let value = if domain.is_fluid(cell) {
            let local = CaloricCoupling {
                law: caloric.law,
                config: MaterialEnthalpyConfig {
                    reference_temperature_k: old,
                    ..caloric.config
                },
            };
            let change = local.sample(new, reference, bounds)?.enthalpy_j_kg;
            if (delta == 0.0 && change != 0.0)
                || (delta != 0.0
                    && (change == 0.0 || change.is_sign_positive() != delta.is_sign_positive()))
            {
                return Err(MaterialBuoyancyError::Invalid(
                    "source enthalpy increment must follow positive Cp",
                ));
            }
            reference.density_kg_m3() * domain.volume(cell) * change
        } else {
            step.reference_capacity_j_k[cell] * delta
        };
        if !value.is_finite() {
            return Err(MaterialBuoyancyError::Invalid(
                "source-caloric stored energy is unrepresentable",
            ));
        }
        stored.push(value);
    }
    Ok(())
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn energy_step(
    caloric: CaloricCoupling<'_>,
    controls: MaterialEnthalpyTransientConfig,
    domain: &VoxelDomain,
    reference: &FluidTransportState,
    bounds: BoussinesqConfig<3>,
    step: &MaterialHeatStep<'_>,
    conductivity_w_m_k: &[[f64; 3]],
    remaining_solves: &mut usize,
    gate: &CancelGate,
) -> Result<
    (
        Vec<f64>,
        TransientRecord,
        MaterialEnthalpyState,
        MaterialEnthalpyStepReport,
    ),
    MaterialBuoyancyError,
> {
    let boundary_j_kg = caloric.boundaries(step.setup, reference, bounds, gate)?;
    let mut iterate = step.previous_temperature_k.to_vec();
    let mut applied = caloric.field(domain, &iterate, reference, bounds, gate)?;
    let mut endpoint = Vec::with_capacity(domain.cell_count());
    let mut trial_stored = Vec::with_capacity(domain.cell_count());
    let mut stored = Vec::with_capacity(domain.cell_count());
    let mut storage_w_k = Vec::with_capacity(domain.cell_count());
    let mut cell_defects = vec![0.0; domain.cell_count()];
    stored_change(
        caloric,
        domain,
        reference,
        bounds,
        step,
        &iterate,
        &mut trial_stored,
        gate,
    )?;
    let mut fixed = vec![false; domain.cell_count()];
    for &(cell, _) in &step.setup.fixed_temperature {
        if let Some(value) = fixed.get_mut(cell) {
            *value = true;
        }
    }
    let starting_budget = *remaining_solves;
    let mut total_iterations = 0usize;
    let mut report = MaterialEnthalpyStepReport {
        time_s: step.time_s,
        passes: 0,
        energy_solves: 0,
        temperature_defect_k: f64::INFINITY,
        maximum_storage_defect_w: f64::INFINITY,
        advection: MaterialEnthalpyReport::default(),
    };
    let mut closure_j = f64::INFINITY;
    for pass in 1..=controls.maximum_passes_per_step {
        poll(gate)?;
        storage_w_k.clear();
        for cell in 0..domain.cell_count() {
            poll(gate)?;
            let delta = iterate[cell] - step.previous_temperature_k[cell];
            let capacity = if let Some(tangent) = applied[cell] {
                let mass = reference.density_kg_m3() * domain.volume(cell);
                let cp = if delta == 0.0 {
                    tangent.heat_capacity_j_kg_k
                } else {
                    trial_stored[cell] / mass / delta
                };
                if !cp.is_finite() || cp <= 0.0 {
                    return Err(MaterialBuoyancyError::Invalid(
                        "source storage secant must be finite and positive",
                    ));
                }
                mass * cp
            } else {
                step.reference_capacity_j_k[cell]
            };
            let coefficient = capacity / step.time_step_s;
            if !coefficient.is_finite() || coefficient <= 0.0 {
                return Err(MaterialBuoyancyError::Invalid(
                    "source storage rate is unrepresentable",
                ));
            }
            storage_w_k.push(coefficient);
        }
        let (mut next, radiated_w, radiation_correction_w) = step.solve_with_enthalpy(
            &storage_w_k,
            conductivity_w_m_k,
            AdvectiveEnthalpyField {
                cells: &applied,
                boundary_j_kg,
            },
            remaining_solves,
        )?;
        total_iterations = total_iterations.saturating_add(next.report.iterations);
        caloric.write_field(
            domain,
            &next.temperature,
            reference,
            bounds,
            &mut endpoint,
            gate,
        )?;
        stored_change(
            caloric,
            domain,
            reference,
            bounds,
            step,
            &next.temperature,
            &mut stored,
            gate,
        )?;
        let advection = defect_field(
            domain,
            step.flow,
            reference.density_kg_m3(),
            &next.temperature,
            &applied,
            &endpoint,
            &mut cell_defects,
            gate,
        )?;
        let mut temperature_defect_k = 0.0f64;
        let mut maximum_storage_defect_w = 0.0f64;
        let mut outflow_correction_w = radiation_correction_w;
        for cell in 0..domain.cell_count() {
            poll(gate)?;
            temperature_defect_k =
                temperature_defect_k.max((next.temperature[cell] - iterate[cell]).abs());
            let defect = stored[cell] / step.time_step_s
                - storage_w_k[cell] * (next.temperature[cell] - step.previous_temperature_k[cell]);
            if !defect.is_finite() {
                return Err(MaterialBuoyancyError::Invalid(
                    "source storage defect is unrepresentable",
                ));
            }
            maximum_storage_defect_w = maximum_storage_defect_w.max(defect.abs());
            // Internal faces cancel. At fixed cells their residual is supplied
            // by the native thermostat: exclude that reservoir's matching
            // injection correction, exactly like endpoint radiation does.
            if !fixed[cell] {
                outflow_correction_w += cell_defects[cell];
            }
        }
        next.report.iterations = total_iterations;
        let record = step.record_with_enthalpy(&next, &stored, radiated_w, outflow_correction_w)?;
        closure_j = record.closure_j;
        report = MaterialEnthalpyStepReport {
            time_s: step.time_s,
            passes: pass,
            energy_solves: starting_budget - *remaining_solves,
            temperature_defect_k,
            maximum_storage_defect_w,
            advection,
        };
        if temperature_defect_k <= controls.temperature_tolerance_k
            && maximum_storage_defect_w <= controls.maximum_storage_defect_w
            && advection.maximum_face_flux_defect_w <= caloric.config.maximum_flux_defect_w
            && advection.maximum_cell_flux_defect_w <= caloric.config.maximum_flux_defect_w
            && closure_j.abs() <= controls.maximum_energy_closure_j
        {
            poll(gate)?;
            return Ok((
                next.temperature,
                record,
                MaterialEnthalpyState {
                    applied,
                    endpoint,
                    boundary_j_kg,
                    report: advection,
                },
                report,
            ));
        }
        if pass < controls.maximum_passes_per_step {
            if controls.temperature_relaxation == 1.0 {
                // The just-admitted endpoint is the next trial. Reuse its
                // exact source queries and buffers rather than resolving the
                // same immutable law twice at every nonlinear iteration.
                iterate = next.temperature;
                core::mem::swap(&mut applied, &mut endpoint);
                core::mem::swap(&mut trial_stored, &mut stored);
            } else {
                for (trial, endpoint) in iterate.iter_mut().zip(next.temperature) {
                    *trial = controls
                        .temperature_relaxation
                        .mul_add(endpoint - *trial, *trial);
                }
                caloric.write_field(domain, &iterate, reference, bounds, &mut applied, gate)?;
                stored_change(
                    caloric,
                    domain,
                    reference,
                    bounds,
                    step,
                    &iterate,
                    &mut trial_stored,
                    gate,
                )?;
            }
        }
    }
    Err(MaterialBuoyancyError::EnthalpyStepNotConverged {
        time_s: step.time_s,
        passes: report.passes,
        temperature_defect_k: report.temperature_defect_k,
        storage_defect_w: report.maximum_storage_defect_w,
        enthalpy: report.advection,
        closure_j,
    })
}

impl<F: FluidTransport + ?Sized> MaterialChtFluid<'_, F> {
    /// March natural/mixed convection with actual source h(T) storage and
    /// donor-cell advection on the existing native momentum and thermal clock.
    ///
    /// At every time step the original SIMPLEC BE/BDF2 momentum history uses
    /// the preceding accepted material density, viscosity and conductivity.
    /// On those same projected fluxes, bounded caloric iterations use positive
    /// source secants for `rho0 V delta(h)/dt` and the original upwind h face
    /// operator. Every trial retains the SAME old accepted temperature; no
    /// nonlinear iteration advances time. True local storage/face defects and
    /// the physical heat closure must all pass before accepting an endpoint.
    ///
    /// This remains FIRST ORDER in coupled time accuracy: source transport
    /// and buoyancy retain the explicit native material lag. Reference rho
    /// owns mass/inertia; original density AND raw/chart Cp departure limits
    /// remain enforced. Fixed solid storage, contacts, thermostat reservoirs
    /// and converged endpoint radiation retain their original owners. No new
    /// latent-plateau, variable-mass, pressure-work or viscous-heating model is
    /// implied. The caller binds caloric and transport laws to the same fluid,
    /// pressure, composition and source context, as for the steady entry.
    ///
    /// # Errors
    /// Refuses unsupported source endpoints/boundaries, nonpositive storage
    /// secants, nonconservative or non-upwind transport, exhausted per-step or
    /// whole-march work, physical-closure failure, and native solver errors.
    /// Cancellation or any refusal returns no partial march. Source callbacks
    /// are repeatable queries, never commit physical specimen history.
    #[allow(clippy::too_many_arguments)]
    pub fn march_natural_convection_with_enthalpy<
        'a,
        S: FluidTransport,
        C: ConvectionCaloricLaw,
    >(
        &self,
        domain: &VoxelDomain,
        solids: &[SolidMaterial],
        flow: &SimpleConfig,
        thermal: &ThermalSetup,
        energy: &EnergyConfig,
        initial_temperature_k: &[f64],
        bounds: BoussinesqConfig<3>,
        caloric: &'a C,
        controls: MaterialEnthalpyTransientConfig,
        unsteady: &UnsteadyConfig,
        inlet_schedule: impl Fn(f64) -> f64,
        power_schedule: impl Fn(f64) -> f64,
        radiation: Option<&RadiationConfig>,
        gate: &CancelGate,
        mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<S, MaterialError>,
    ) -> Result<MaterialEnthalpyMarch<'a, S, C>, MaterialBuoyancyError> {
        poll(gate)?;
        controls.validate()?;
        if controls.maximum_energy_solves < unsteady.steps {
            return Err(ChtError::InvalidInput {
                field: "unsteady.material_energy_solves",
                reason: "whole-march budget must cover at least one energy solve per time step"
                    .into(),
            }
            .into());
        }
        let coupling = CaloricCoupling {
            law: caloric,
            config: controls.caloric,
        };
        coupling.validate(energy)?;
        bounds.validate().map_err(|_| {
            MaterialBuoyancyError::Invalid(
                "gravity and reference density/heat-capacity departure bounds must be admitted",
            )
        })?;
        // Admit actual inlet and potential-backflow transport only for this
        // source-caloric entry. The shared source-solid/fixed-Cp APIs retain
        // their existing admission and material-refresh behavior.
        admit_open_boundary_materials(
            domain,
            thermal,
            &self.transport,
            bounds,
            gate,
            &mut resolve,
        )?;
        let mut remaining_solves = controls.maximum_energy_solves;
        let mut records = Vec::with_capacity(unsteady.steps);
        let mut last_state = None;
        let mut thermal_step = |step: &MaterialHeatStep<'_>, fields: &mut MaterialCoefficients| {
            let (next, record, state, defects) = energy_step(
                coupling,
                controls,
                domain,
                &self.transport,
                bounds,
                step,
                &fields.conductivity_w_m_k,
                &mut remaining_solves,
                gate,
            )?;
            last_state = Some(state);
            records.push(defects);
            Ok((next, record))
        };
        let convection = self.march_natural_convection_inner(
            domain,
            solids,
            flow,
            thermal,
            energy,
            initial_temperature_k,
            bounds,
            unsteady,
            inlet_schedule,
            power_schedule,
            radiation,
            Some(&mut thermal_step),
            gate,
            resolve,
        )?;
        poll(gate)?;
        Ok(MaterialEnthalpyMarch {
            convection,
            caloric,
            config: controls,
            enthalpy: last_state.expect("native march admits at least one accepted step"),
            records,
        })
    }
}
