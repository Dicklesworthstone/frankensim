//! Joint source-caloric fluid and solid storage on one native thermal endpoint.
//!
//! This closes the previously separate nonlinear caloric owners in the SAME
//! MaterialHeatStep callback. Fluid h(T), solid integral Cp(T), solid k(T),
//! radiation, thermostats and conservative advection therefore see one trial
//! temperature and one physical energy ledger. No secondary thermal clock or
//! operator splitting is introduced.

use super::{
    MaterialEnthalpyMarch, MaterialEnthalpyStepReport, MaterialEnthalpyTransientConfig,
};
use super::super::{
    CaloricCoupling, ConvectionCaloricLaw, MaterialEnthalpyConfig, MaterialEnthalpyReport,
    MaterialEnthalpyState, defect_field,
};
use super::super::super::{
    MaterialBuoyancyError, MaterialBuoyancyMarch, admit_open_boundary_materials, poll,
};
use crate::conjugate::unsteady::{MaterialCoefficients, MaterialHeatStep};
use crate::conjugate::{
    AdvectiveEnthalpyField, ChtError, EnergyConfig, RadiationConfig, SimpleConfig,
    ThermalSetup, TransientRecord, UnsteadyConfig, Voxel, VoxelDomain,
};
use crate::material::conjugate::MaterialChtFluid;
use crate::material::conjugate::buoyant::{SolidThermalConfig, SolidThermalRecord};
use crate::material::conjugate::solid::{
    SolidConductivitySample, SolidThermalError, SolidThermalLaw,
};
use crate::material::transport::boussinesq::BoussinesqConfig;
use fs_exec::CancelGate;
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState},
};

/// One accepted native march with both original fluid and solid caloric owners.
pub struct MaterialEnthalpySolidMarch<'a, F, C: ConvectionCaloricLaw> {
    /// Original source-enthalpy fluid march and all native flow/heat records.
    pub fluid: MaterialEnthalpyMarch<'a, F, C>,
    /// Actual original solid sources at the returned thermal endpoint.
    pub solid_endpoint: Vec<SolidConductivitySample>,
    /// One joint nonlinear solid record per native thermal interval.
    pub solid_records: Vec<SolidThermalRecord>,
    /// Per voxel sensible energy relative to each solid storage law's lower
    /// declared temperature. Fluid entries are zero.
    pub solid_energy_j: Vec<f64>,
}

fn solid_controls_valid(controls: SolidThermalConfig) -> bool {
    controls.maximum_iterations > 0
        && [
            controls.temperature_tolerance_k,
            controls.conductivity_relative_tolerance,
            controls.storage_tolerance_j,
        ]
        .iter()
        .all(|v| v.is_finite() && *v >= 0.0)
        && controls.conductivity_relative_tolerance < 1.0
        && controls.temperature_relaxation.is_finite()
        && controls.temperature_relaxation > 0.0
        && controls.temperature_relaxation <= 1.0
}

fn solid_failure(
    cell: usize,
    temperature_k: f64,
    source: SolidThermalError,
) -> MaterialBuoyancyError {
    MaterialBuoyancyError::SolidThermal {
        cell,
        temperature_k,
        source,
    }
}

fn solid_samples(
    domain: &VoxelDomain,
    sources: &[SolidThermalLaw<'_>],
    temperature: &[f64],
    gate: &CancelGate,
) -> Result<Vec<SolidConductivitySample>, MaterialBuoyancyError> {
    let mut samples = Vec::new();
    for (cell, &temperature_k) in temperature.iter().enumerate() {
        poll(gate)?;
        let Voxel::Solid(index) = domain.voxel_at(cell) else {
            continue;
        };
        let source = sources
            .get(usize::from(index))
            .ok_or(MaterialBuoyancyError::Invalid(
                "every solid assignment needs its original thermal material law",
            ))?;
        samples.push(
            source
                .at_temperature(cell, index, temperature_k)
                .map_err(|error| solid_failure(cell, temperature_k, error))?,
        );
    }
    Ok(samples)
}

/// Fill exact physical storage increments for both source families. The parent
/// helper already computes exact fluid h(T) and fixed solid storage; solid
/// entries are then replaced by their original source integrals.
#[allow(clippy::too_many_arguments)]
fn stored_change(
    caloric: CaloricCoupling<'_>,
    solid_sources: &[SolidThermalLaw<'_>],
    domain: &VoxelDomain,
    reference: &FluidTransportState,
    bounds: BoussinesqConfig<3>,
    step: &MaterialHeatStep<'_>,
    temperature: &[f64],
    stored: &mut Vec<f64>,
    gate: &CancelGate,
) -> Result<(), MaterialBuoyancyError> {
    super::stored_change(
        caloric,
        domain,
        reference,
        bounds,
        step,
        temperature,
        stored,
        gate,
    )?;
    for cell in 0..domain.cell_count() {
        poll(gate)?;
        let Voxel::Solid(index) = domain.voxel_at(cell) else {
            continue;
        };
        let source = solid_sources
            .get(usize::from(index))
            .ok_or(MaterialBuoyancyError::Invalid(
                "every solid assignment needs its original thermal material law",
            ))?;
        let old = step.previous_temperature_k[cell];
        let new = temperature[cell];
        let actual = domain.volume(cell)
            * source
                .storage()
                .energy_change_j(old, new)
                .map_err(|error| {
                    solid_failure(cell, new, SolidThermalError::Storage(error))
                })?;
        if !actual.is_finite() {
            return Err(MaterialBuoyancyError::Invalid(
                "joint source-caloric stored energy is unrepresentable",
            ));
        }
        stored[cell] = actual;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn joint_energy_step(
    caloric: CaloricCoupling<'_>,
    fluid_controls: MaterialEnthalpyTransientConfig,
    solid_sources: &[SolidThermalLaw<'_>],
    solid_controls: SolidThermalConfig,
    domain: &VoxelDomain,
    reference: &FluidTransportState,
    bounds: BoussinesqConfig<3>,
    step: &MaterialHeatStep<'_>,
    fields: &mut MaterialCoefficients,
    remaining_solves: &mut usize,
    gate: &CancelGate,
) -> Result<
    (
        Vec<f64>,
        TransientRecord,
        MaterialEnthalpyState,
        MaterialEnthalpyStepReport,
        Vec<SolidConductivitySample>,
        SolidThermalRecord,
    ),
    MaterialBuoyancyError,
> {
    let boundary_j_kg = caloric.boundaries(step.setup, reference, bounds, gate)?;
    let mut iterate = step.previous_temperature_k.to_vec();
    let mut applied = caloric.field(domain, &iterate, reference, bounds, gate)?;
    let mut endpoint = Vec::with_capacity(domain.cell_count());
    let mut solid_applied = solid_samples(domain, solid_sources, &iterate, gate)?;
    let mut solid_endpoint = Vec::new();
    let mut trial_stored = Vec::with_capacity(domain.cell_count());
    let mut stored = Vec::with_capacity(domain.cell_count());
    let mut storage_w_k = Vec::with_capacity(domain.cell_count());
    let mut cell_defects = vec![0.0; domain.cell_count()];
    stored_change(
        caloric,
        solid_sources,
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
    let maximum_passes = fluid_controls
        .maximum_passes_per_step
        .min(solid_controls.maximum_iterations);
    let relaxation = fluid_controls
        .temperature_relaxation
        .min(solid_controls.temperature_relaxation);
    let mut total_iterations = 0usize;
    let mut fluid_report = MaterialEnthalpyStepReport {
        time_s: step.time_s,
        passes: 0,
        energy_solves: 0,
        temperature_defect_k: f64::INFINITY,
        maximum_storage_defect_w: f64::INFINITY,
        advection: MaterialEnthalpyReport::default(),
    };
    let mut solid_report = SolidThermalRecord {
        time_s: step.time_s,
        iterations: 0,
        temperature_defect_k: f64::INFINITY,
        conductivity_relative_defect: f64::INFINITY,
        storage_defect_j: f64::INFINITY,
        solid_stored_energy_change_j: 0.0,
    };
    let mut closure_j = f64::INFINITY;

    for pass in 1..=maximum_passes {
        poll(gate)?;
        for sample in &solid_applied {
            fields.conductivity_w_m_k[sample.cell] = sample.conductivity_w_m_k;
        }
        storage_w_k.clear();
        for cell in 0..domain.cell_count() {
            poll(gate)?;
            let old = step.previous_temperature_k[cell];
            let delta = iterate[cell] - old;
            let capacity = if let Some(tangent) = applied[cell] {
                let mass = reference.density_kg_m3() * domain.volume(cell);
                let cp = if delta == 0.0 {
                    tangent.heat_capacity_j_kg_k
                } else {
                    trial_stored[cell] / mass / delta
                };
                mass * cp
            } else {
                let Voxel::Solid(index) = domain.voxel_at(cell) else {
                    return Err(MaterialBuoyancyError::Invalid(
                        "joint caloric field omitted a fluid cell",
                    ));
                };
                let source = &solid_sources[usize::from(index)];
                if delta == 0.0 {
                    domain.volume(cell)
                        * source
                            .storage()
                            .heat_capacity_j_k(old)
                            .map_err(|error| {
                                solid_failure(cell, old, SolidThermalError::Storage(error))
                            })?
                } else {
                    trial_stored[cell] / delta
                }
            };
            let coefficient = capacity / step.time_step_s;
            if !capacity.is_finite()
                || capacity <= 0.0
                || !coefficient.is_finite()
                || coefficient <= 0.0
            {
                return Err(MaterialBuoyancyError::Invalid(
                    "joint source storage secant must be finite and positive",
                ));
            }
            storage_w_k.push(coefficient);
        }

        let (mut next, radiated_w, radiation_correction_w) = step.solve_with_enthalpy(
            &storage_w_k,
            &fields.conductivity_w_m_k,
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
        solid_endpoint = solid_samples(domain, solid_sources, &next.temperature, gate)?;
        stored_change(
            caloric,
            solid_sources,
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

        let mut temperature_defect_k = 0.0_f64;
        let mut largest_temperature = 0.0_f64;
        let mut maximum_fluid_storage_defect_w = 0.0_f64;
        let mut maximum_solid_storage_defect_j = 0.0_f64;
        let mut solid_stored_energy_change_j = 0.0_f64;
        let mut solid_roundoff_j = 0.0_f64;
        let mut outflow_correction_w = radiation_correction_w;
        for cell in 0..domain.cell_count() {
            poll(gate)?;
            let old = step.previous_temperature_k[cell];
            let new = next.temperature[cell];
            temperature_defect_k = temperature_defect_k.max((new - iterate[cell]).abs());
            largest_temperature = largest_temperature.max(iterate[cell].abs()).max(new.abs());
            let linear_j =
                storage_w_k[cell] * step.time_step_s * (new - old);
            let defect_j = stored[cell] - linear_j;
            if !defect_j.is_finite() {
                return Err(MaterialBuoyancyError::Invalid(
                    "joint source storage defect is unrepresentable",
                ));
            }
            if domain.is_fluid(cell) {
                maximum_fluid_storage_defect_w =
                    maximum_fluid_storage_defect_w.max((defect_j / step.time_step_s).abs());
                if !fixed[cell] {
                    outflow_correction_w += cell_defects[cell];
                }
            } else {
                maximum_solid_storage_defect_j =
                    maximum_solid_storage_defect_j.max(defect_j.abs());
                solid_stored_energy_change_j += stored[cell];
                solid_roundoff_j += 64.0
                    * f64::EPSILON
                    * storage_w_k[cell]
                    * step.time_step_s
                    * (old.abs() + new.abs());
            }
        }

        let mut conductivity_relative_defect = 0.0_f64;
        for (before, after) in solid_applied.iter().zip(&solid_endpoint) {
            if before.cell != after.cell || before.material_index != after.material_index {
                return Err(MaterialBuoyancyError::Invalid(
                    "solid source ordering changed inside one thermal interval",
                ));
            }
            for (a, b) in before
                .conductivity_w_m_k
                .iter()
                .copied()
                .zip(after.conductivity_w_m_k)
            {
                conductivity_relative_defect = conductivity_relative_defect
                    .max((b - a).abs() / a.max(b));
            }
        }

        next.report.iterations = total_iterations;
        let record = step.record_with_enthalpy(
            &next,
            &stored,
            radiated_w,
            outflow_correction_w,
        )?;
        closure_j = record.closure_j;
        fluid_report = MaterialEnthalpyStepReport {
            time_s: step.time_s,
            passes: pass,
            energy_solves: starting_budget - *remaining_solves,
            temperature_defect_k,
            maximum_storage_defect_w: maximum_fluid_storage_defect_w,
            advection,
        };
        solid_report = SolidThermalRecord {
            time_s: step.time_s,
            iterations: pass,
            temperature_defect_k,
            conductivity_relative_defect,
            storage_defect_j: maximum_solid_storage_defect_j,
            solid_stored_energy_change_j,
        };

        let temperature_roundoff =
            64.0 * f64::EPSILON * largest_temperature;
        let solid_storage_ok = maximum_solid_storage_defect_j
            <= solid_controls.storage_tolerance_j + solid_roundoff_j;
        let fluid_ok = temperature_defect_k
            <= fluid_controls.temperature_tolerance_k + temperature_roundoff
            && maximum_fluid_storage_defect_w
                <= fluid_controls.maximum_storage_defect_w
            && advection.maximum_face_flux_defect_w
                <= caloric.config.maximum_flux_defect_w
            && advection.maximum_cell_flux_defect_w
                <= caloric.config.maximum_flux_defect_w
            && closure_j.abs() <= fluid_controls.maximum_energy_closure_j;
        let solid_ok = temperature_defect_k
            <= solid_controls.temperature_tolerance_k + temperature_roundoff
            && conductivity_relative_defect
                <= solid_controls.conductivity_relative_tolerance
            && solid_storage_ok
            && closure_j.abs() <= solid_controls.storage_tolerance_j + solid_roundoff_j;
        if fluid_ok && solid_ok {
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
                fluid_report,
                solid_endpoint,
                solid_report,
            ));
        }

        if pass < maximum_passes {
            if relaxation == 1.0 {
                iterate = next.temperature;
                core::mem::swap(&mut applied, &mut endpoint);
                core::mem::swap(&mut solid_applied, &mut solid_endpoint);
                core::mem::swap(&mut trial_stored, &mut stored);
            } else {
                for (trial, target) in iterate.iter_mut().zip(next.temperature) {
                    *trial = relaxation.mul_add(target - *trial, *trial);
                }
                caloric.write_field(
                    domain,
                    &iterate,
                    reference,
                    bounds,
                    &mut applied,
                    gate,
                )?;
                solid_applied = solid_samples(domain, solid_sources, &iterate, gate)?;
                stored_change(
                    caloric,
                    solid_sources,
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

    if solid_report.temperature_defect_k > solid_controls.temperature_tolerance_k
        || solid_report.conductivity_relative_defect
            > solid_controls.conductivity_relative_tolerance
        || solid_report.storage_defect_j > solid_controls.storage_tolerance_j
        || closure_j.abs() > solid_controls.storage_tolerance_j
    {
        Err(MaterialBuoyancyError::SolidNotConverged(solid_report))
    } else {
        Err(MaterialBuoyancyError::EnthalpyStepNotConverged {
            time_s: step.time_s,
            passes: fluid_report.passes,
            temperature_defect_k: fluid_report.temperature_defect_k,
            storage_defect_w: fluid_report.maximum_storage_defect_w,
            enthalpy: fluid_report.advection,
            closure_j,
        })
    }
}

impl<F: FluidTransport + ?Sized> MaterialChtFluid<'_, F> {
    /// Advance source h(T) fluid storage/advection and source-integral solid
    /// storage/conductivity in one bounded native thermal endpoint iteration.
    ///
    /// Momentum/material feedback keeps the existing first-order lag. Inside
    /// each accepted heat interval, however, fluid calorics, solid calorics,
    /// solid conductivity, endpoint radiation and thermostats share the same
    /// trial temperature and energy ledger. The stricter of the two requested
    /// relaxation factors is used; neither acceptance tolerance is weakened.
    ///
    /// # Errors
    /// Refuses unsupported source states, malformed work/tolerance controls,
    /// nonpositive secants, nonconservative enthalpy transport, exhausted
    /// energy work, nonlinear closure failure, native solver errors or
    /// cancellation. No partial march is published.
    #[allow(clippy::too_many_arguments)]
    pub fn march_natural_convection_with_enthalpy_and_material_solids<
        'a,
        S: FluidTransport,
        C: ConvectionCaloricLaw,
    >(
        &self,
        domain: &VoxelDomain,
        solid_sources: &[SolidThermalLaw<'_>],
        flow: &SimpleConfig,
        thermal: &ThermalSetup,
        energy: &EnergyConfig,
        initial_temperature_k: &[f64],
        bounds: BoussinesqConfig<3>,
        caloric: &'a C,
        fluid_controls: MaterialEnthalpyTransientConfig,
        solid_controls: SolidThermalConfig,
        unsteady: &UnsteadyConfig,
        inlet_schedule: impl Fn(f64) -> f64,
        power_schedule: impl Fn(f64) -> f64,
        radiation: Option<&RadiationConfig>,
        gate: &CancelGate,
        mut resolve: impl FnMut(
            usize,
            f64,
            &FluidTransportState,
        ) -> Result<S, MaterialError>,
    ) -> Result<MaterialEnthalpySolidMarch<'a, S, C>, MaterialBuoyancyError> {
        poll(gate)?;
        fluid_controls.validate()?;
        if !solid_controls_valid(solid_controls) {
            return Err(MaterialBuoyancyError::Invalid(
                "joint solid heat iteration needs a positive budget, finite nonnegative tolerances, relative tolerance below one and relaxation in (0,1]",
            ));
        }
        if fluid_controls.maximum_energy_solves < unsteady.steps {
            return Err(ChtError::InvalidInput {
                field: "unsteady.material_energy_solves",
                reason: "whole-march budget must cover at least one energy solve per time step"
                    .into(),
            }
            .into());
        }
        let coupling = CaloricCoupling {
            law: caloric,
            config: fluid_controls.caloric,
        };
        coupling.validate(energy)?;
        bounds.validate().map_err(|_| {
            MaterialBuoyancyError::Invalid(
                "gravity and reference density/heat-capacity departure bounds must be admitted",
            )
        })?;
        admit_open_boundary_materials(
            domain,
            thermal,
            &self.transport,
            bounds,
            gate,
            &mut resolve,
        )?;

        let mut solids = Vec::with_capacity(solid_sources.len());
        for source in solid_sources {
            solids.push(source.reference_material().map_err(|_| {
                MaterialBuoyancyError::Invalid(
                    "solid reference storage is not representable",
                )
            })?);
        }

        let mut remaining_solves = fluid_controls.maximum_energy_solves;
        let mut fluid_records = Vec::with_capacity(unsteady.steps);
        let mut solid_records = Vec::with_capacity(unsteady.steps);
        let mut last_state = None;
        let mut last_solid_endpoint = Vec::new();
        let mut thermal_step =
            |step: &MaterialHeatStep<'_>, fields: &mut MaterialCoefficients| {
                let (next, record, state, fluid_record, solid_endpoint, solid_record) =
                    joint_energy_step(
                        coupling,
                        fluid_controls,
                        solid_sources,
                        solid_controls,
                        domain,
                        &self.transport,
                        bounds,
                        step,
                        fields,
                        &mut remaining_solves,
                        gate,
                    )?;
                last_state = Some(state);
                last_solid_endpoint = solid_endpoint;
                fluid_records.push(fluid_record);
                solid_records.push(solid_record);
                Ok((next, record))
            };

        let convection: MaterialBuoyancyMarch<S> = self.march_natural_convection_inner(
            domain,
            &solids,
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

        let mut solid_energy_j = vec![0.0; domain.cell_count()];
        for sample in &last_solid_endpoint {
            poll(gate)?;
            let source = &solid_sources[usize::from(sample.material_index)];
            let low = source.storage().temperature_span_k()[0];
            let energy_j = domain.volume(sample.cell)
                * source
                    .storage()
                    .energy_change_j(low, sample.temperature_k)
                    .map_err(|error| {
                        solid_failure(
                            sample.cell,
                            sample.temperature_k,
                            SolidThermalError::Storage(error),
                        )
                    })?;
            if !energy_j.is_finite() {
                return Err(MaterialBuoyancyError::Invalid(
                    "joint solid sensible energy is unrepresentable",
                ));
            }
            solid_energy_j[sample.cell] = energy_j;
        }

        Ok(MaterialEnthalpySolidMarch {
            fluid: MaterialEnthalpyMarch {
                convection,
                caloric,
                config: fluid_controls,
                enthalpy: last_state
                    .expect("native march admits at least one accepted step"),
                records: fluid_records,
            },
            solid_endpoint: last_solid_endpoint,
            solid_records,
            solid_energy_j,
        })
    }
}
