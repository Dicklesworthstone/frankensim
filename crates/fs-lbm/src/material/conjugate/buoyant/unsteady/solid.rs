//! Source-integral solid storage on the existing native thermal step.

use super::{MaterialBuoyancyMarch, MaterialChtFluid};
use crate::conjugate::unsteady::{MaterialCoefficients, MaterialHeatStep, compact_owners};
use crate::conjugate::{
    EnergyConfig, RadiationConfig, SimpleConfig, ThermalSetup, TransientRecord, UnsteadyConfig,
    Voxel, VoxelDomain,
};
use crate::material::conjugate::buoyant::{MaterialBuoyancyError, poll};
use crate::material::conjugate::solid::{
    SolidConductivitySample, SolidThermalError, SolidThermalLaw,
};
use crate::material::transport::boussinesq::BoussinesqConfig;
use fs_exec::CancelGate;
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState},
};

/// Bounded implicit solid constitutive iteration inside each native heat step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolidThermalConfig {
    /// Maximum source/conduction/storage iterations within one native step.
    pub maximum_iterations: usize,
    /// Absolute unrelaxed output-minus-trial temperature defect, K.
    pub temperature_tolerance_k: f64,
    /// Relative change in any grid-aligned solid conductivity component.
    pub conductivity_relative_tolerance: f64,
    /// Absolute per-cell storage defect AND corrected whole-step closure limit.
    /// Floating-point roundoff in the existing temperature solve is allowed
    /// separately; it is not evidence of a material approximation error.
    pub storage_tolerance_j: f64,
    /// Relaxation of constitutive trials only, in (0, 1].
    pub temperature_relaxation: f64,
}

impl Default for SolidThermalConfig {
    fn default() -> Self {
        Self {
            maximum_iterations: 100,
            temperature_tolerance_k: 1e-9,
            conductivity_relative_tolerance: 1e-9,
            storage_tolerance_j: 1e-8,
            temperature_relaxation: 1.0,
        }
    }
}

impl SolidThermalConfig {
    fn validate(self) -> Result<(), MaterialBuoyancyError> {
        if self.maximum_iterations == 0
            || [
                self.temperature_tolerance_k,
                self.conductivity_relative_tolerance,
                self.storage_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
            || self.conductivity_relative_tolerance >= 1.0
            || !self.temperature_relaxation.is_finite()
            || self.temperature_relaxation <= 0.0
            || self.temperature_relaxation > 1.0
        {
            return Err(MaterialBuoyancyError::Invalid(
                "solid heat iteration needs a positive budget, finite nonnegative tolerances, relative tolerance below one and relaxation in (0,1]",
            ));
        }
        Ok(())
    }
}

/// Actual endpoint constitutive defects and source-integrated solid storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolidThermalRecord {
    /// End of the original native time interval, s.
    pub time_s: f64,
    /// Constitutive iterations, including any unsupported linear predictors.
    pub iterations: usize,
    /// Maximum unrelaxed output-minus-trial temperature defect, K.
    pub temperature_defect_k: f64,
    /// Maximum relative conductivity difference between applied and endpoint sources.
    pub conductivity_relative_defect: f64,
    /// Maximum absolute per-cell exact-integral minus secant-storage defect, J.
    pub storage_defect_j: f64,
    /// Actual total solid sensible-storage increment over this interval, J.
    pub solid_stored_energy_change_j: f64,
}

/// One native fluid/solid march with original solid source and caloric evidence.
pub struct MaterialSolidBuoyancyMarch<F> {
    /// The original fluid buoyancy march, with the true combined heat ledger.
    pub fluid: MaterialBuoyancyMarch<F>,
    /// Actual source queries at the returned thermal endpoint. Their difference
    /// from the last implicit conductivity iterate passed `solid_records`.
    pub solid_endpoint: Vec<SolidConductivitySample>,
    /// One admitted solid constitutive/energy record per original native step.
    pub solid_records: Vec<SolidThermalRecord>,
    /// Per original voxel: zero for fluid and compact aliases, source-solid
    /// sensible energy relative to that material's declared lower Cp temperature.
    /// This is not formation energy or authored junction storage.
    pub solid_energy_j: Vec<f64>,
}

impl<F: FluidTransport + ?Sized> MaterialChtFluid<'_, F> {
    /// Advance temperature-dependent solid conduction and exact source storage
    /// together with the original material-buoyant native flow march.
    ///
    /// Every ordinary solid assignment names an original `SolidThermalLaw`. Its fixed
    /// reference mass multiplies the exact integral of selected source Cp(T),
    /// while the same source bundle supplies grid-aligned k(T). Inside each
    /// original backward-Euler thermal interval, source-integral secants and
    /// conductivity are iterated to a common endpoint. Radiation, contacts,
    /// imposed-temperature reservoirs and conservative fluxes retain their
    /// existing owner. The returned heat ledger uses actual source integrals,
    /// including the storage supplied by fixed-temperature reservoirs.
    ///
    /// Material fluid forcing/transport remains explicitly first order. Solid
    /// thermodynamics is single-phase at the source's held conditions; no
    /// density-driven mass change, melting plateau, deformation or extra clock
    /// is introduced. Unsupported nonlinear trials may backtrack, but no
    /// clamped or unsupported endpoint is returned. Source, convergence or
    /// cancellation failure publishes none of the candidate march.
    #[allow(clippy::too_many_arguments)]
    pub fn march_natural_convection_with_material_solids<S: FluidTransport>(
        &self,
        domain: &VoxelDomain,
        solid_sources: &[SolidThermalLaw<'_>],
        flow_config: &SimpleConfig,
        thermal: &ThermalSetup,
        energy_config: &EnergyConfig,
        initial_temperature_k: &[f64],
        config: BoussinesqConfig<3>,
        unsteady: &UnsteadyConfig,
        solid_controls: SolidThermalConfig,
        inlet_schedule: impl Fn(f64) -> f64,
        power_schedule: impl Fn(f64) -> f64,
        radiation: Option<&RadiationConfig>,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<S, MaterialError>,
    ) -> Result<MaterialSolidBuoyancyMarch<S>, MaterialBuoyancyError> {
        poll(gate)?;
        solid_controls.validate()?;
        if initial_temperature_k.len() != domain.cell_count() {
            return Err(MaterialBuoyancyError::Invalid(
                "one initial temperature per original voxel required",
            ));
        }
        let mut solids = Vec::with_capacity(solid_sources.len());
        for source in solid_sources {
            solids.push(source.reference_material().map_err(|_| {
                MaterialBuoyancyError::Invalid("solid reference storage is not representable")
            })?);
        }
        let owners = compact_owners(domain, &thermal.compact_components, gate)?;
        let endpoint = samples(
            domain,
            solid_sources,
            initial_temperature_k,
            |cell| owners.get(cell).is_some_and(|&owner| owner != usize::MAX),
            gate,
        )?;
        drop(owners);
        let mut fixed = vec![false; domain.cell_count()];
        for &(cell, _) in &thermal.fixed_temperature {
            if let Some(value) = fixed.get_mut(cell) {
                *value = true;
            }
        }
        let mut solid = SolidHeat {
            sources: solid_sources,
            controls: solid_controls,
            gate,
            endpoint,
            fixed,
            capacity: vec![0.0; domain.cell_count()],
            storage: vec![0.0; domain.cell_count()],
            records: Vec::new(),
        };
        let mut advance = |step: &MaterialHeatStep<'_>, fields: &mut MaterialCoefficients| {
            solid.advance(step, fields)
        };
        let fluid = self.march_natural_convection_inner(
            domain,
            &solids,
            flow_config,
            thermal,
            energy_config,
            initial_temperature_k,
            config,
            unsteady,
            inlet_schedule,
            power_schedule,
            radiation,
            Some(&mut advance),
            gate,
            resolve,
        )?;
        let mut solid_energy_j = vec![0.0; domain.cell_count()];
        for sample in &solid.endpoint {
            poll(gate)?;
            let source = &solid_sources[usize::from(sample.material_index)];
            solid_energy_j[sample.cell] = domain.volume(sample.cell)
                * source
                    .storage()
                    .energy_change_j(
                        source.storage().temperature_span_k()[0],
                        sample.temperature_k,
                    )
                    .map_err(|error| {
                        failure(
                            sample.cell,
                            sample.temperature_k,
                            SolidThermalError::Storage(error),
                        )
                    })?;
            finite(solid_energy_j[sample.cell])?;
        }
        poll(gate)?;
        Ok(MaterialSolidBuoyancyMarch {
            fluid,
            solid_endpoint: solid.endpoint,
            solid_records: solid.records,
            solid_energy_j,
        })
    }
}

fn failure(cell: usize, temperature_k: f64, source: SolidThermalError) -> MaterialBuoyancyError {
    MaterialBuoyancyError::SolidThermal {
        cell,
        temperature_k,
        source,
    }
}

fn finite(value: f64) -> Result<(), MaterialBuoyancyError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(MaterialBuoyancyError::Invalid(
            "solid heat storage or convergence account is unrepresentable",
        ))
    }
}

fn samples(
    domain: &VoxelDomain,
    sources: &[SolidThermalLaw<'_>],
    temperature: &[f64],
    is_compact: impl Fn(usize) -> bool,
    gate: &CancelGate,
) -> Result<Vec<SolidConductivitySample>, MaterialBuoyancyError> {
    let mut output = Vec::new();
    for (cell, &temperature_k) in temperature.iter().enumerate() {
        poll(gate)?;
        if is_compact(cell) {
            continue;
        }
        if let Voxel::Solid(index) = domain.voxel_at(cell) {
            let source = sources
                .get(usize::from(index))
                .ok_or(MaterialBuoyancyError::Invalid(
                    "every solid assignment needs its original thermal material law",
                ))?;
            output.push(
                source
                    .at_temperature(cell, index, temperature_k)
                    .map_err(|error| failure(cell, temperature_k, error))?,
            );
        }
    }
    Ok(output)
}

struct SolidHeat<'a, 'm> {
    sources: &'a [SolidThermalLaw<'m>],
    controls: SolidThermalConfig,
    gate: &'a CancelGate,
    endpoint: Vec<SolidConductivitySample>,
    fixed: Vec<bool>,
    capacity: Vec<f64>,
    storage: Vec<f64>,
    records: Vec<SolidThermalRecord>,
}

impl SolidHeat<'_, '_> {
    fn advance(
        &mut self,
        step: &MaterialHeatStep<'_>,
        fields: &mut MaterialCoefficients,
    ) -> Result<(Vec<f64>, TransientRecord), MaterialBuoyancyError> {
        self.capacity.resize(step.reference_capacity_j_k.len(), 0.0);
        self.capacity.copy_from_slice(step.reference_capacity_j_k);
        self.storage.resize(self.capacity.len(), 0.0);
        let mut trial = step.previous_temperature_k.to_vec();
        let mut applied = core::mem::take(&mut self.endpoint);
        let mut linear_iterations = 0usize;
        for iteration in 1..=self.controls.maximum_iterations {
            poll(self.gate)?;
            for sample in &applied {
                let cell = sample.cell;
                let storage = self.sources[usize::from(sample.material_index)].storage();
                let old = step.previous_temperature_k[cell];
                let delta = trial[cell] - old;
                let capacity = if delta == 0.0 {
                    storage.heat_capacity_j_k(old)
                } else {
                    storage
                        .energy_change_j(old, trial[cell])
                        .map(|energy| energy / delta)
                }
                .map_err(|error| failure(cell, trial[cell], SolidThermalError::Storage(error)))?;
                self.capacity[cell] = step.domain.volume(cell) * capacity;
                fields.conductivity_w_m_k[cell] = sample.conductivity_w_m_k;
            }
            for (slot, (rate, &capacity)) in
                self.storage.iter_mut().zip(&self.capacity).enumerate()
            {
                *rate = capacity / step.time_step_s;
                let alias = slot < step.domain.cell_count() && step.is_compact_cell(slot);
                let invalid_storage = if alias {
                    capacity != 0.0 || *rate != 0.0
                } else {
                    capacity <= 0.0 || *rate <= 0.0
                };
                if !capacity.is_finite() || !rate.is_finite() || invalid_storage {
                    return Err(MaterialBuoyancyError::Invalid(
                        "finite positive storage and rates required outside compact aliases",
                    ));
                }
            }
            let (next, mut heat) =
                step.solve(&self.capacity, &self.storage, &fields.conductivity_w_m_k)?;
            linear_iterations = linear_iterations.saturating_add(heat.iterations);
            let endpoint = match samples(
                step.domain,
                self.sources,
                &next,
                |cell| step.is_compact_cell(cell),
                self.gate,
            ) {
                Ok(endpoint) => endpoint,
                Err(error) => {
                    if iteration == self.controls.maximum_iterations
                        || matches!(&error, MaterialBuoyancyError::SolidThermal { cell, .. } if self.fixed[*cell])
                    {
                        return Err(error);
                    }
                    (trial, applied) = self.backtrack(step, &trial, &next, error)?;
                    continue;
                }
            };
            let mut report = SolidThermalRecord {
                time_s: step.time_s,
                iterations: iteration,
                temperature_defect_k: 0.0,
                conductivity_relative_defect: 0.0,
                storage_defect_j: 0.0,
                solid_stored_energy_change_j: 0.0,
            };
            let mut largest_temperature = 0.0_f64;
            for (&a, &b) in trial.iter().zip(&next) {
                report.temperature_defect_k = report.temperature_defect_k.max((b - a).abs());
                largest_temperature = largest_temperature.max(a.abs()).max(b.abs());
            }
            for (a, b) in applied.iter().zip(&endpoint) {
                for (x, y) in a.conductivity_w_m_k.iter().zip(b.conductivity_w_m_k) {
                    report.conductivity_relative_defect = report
                        .conductivity_relative_defect
                        .max((y - x).abs() / x.max(y));
                }
            }
            let mut true_storage = step.junction_stored_energy_change_j(&next)?;
            let mut fixed_correction = 0.0;
            let mut roundoff = 0.0;
            let mut storage_admitted = true;
            for (cell, (&old, &new)) in step.previous_temperature_k.iter().zip(&next).enumerate() {
                poll(self.gate)?;
                if step.is_compact_cell(cell) {
                    continue;
                }
                let linear = self.capacity[cell] * (new - old);
                let actual = match step.domain.voxel_at(cell) {
                    Voxel::Fluid => linear,
                    Voxel::Solid(index) => {
                        let actual = step.domain.volume(cell)
                            * self.sources[usize::from(index)]
                                .storage()
                                .energy_change_j(old, new)
                                .map_err(|error| {
                                    failure(cell, new, SolidThermalError::Storage(error))
                                })?;
                        report.solid_stored_energy_change_j += actual;
                        actual
                    }
                };
                let allowance = 64.0 * f64::EPSILON * self.capacity[cell] * (old.abs() + new.abs());
                finite(actual)?;
                finite(linear)?;
                finite(allowance)?;
                let defect = actual - linear;
                report.storage_defect_j = report.storage_defect_j.max(defect.abs());
                storage_admitted &= defect.abs() <= self.controls.storage_tolerance_j + allowance;
                if self.fixed[cell] {
                    fixed_correction += defect;
                }
                true_storage += actual;
                roundoff += allowance;
            }
            heat.stored_energy_change_j = true_storage;
            // The fixed reservoir supplies the same TRUE storage change, once.
            heat.boundary_outflow_j -= fixed_correction;
            heat.closure_j = true_storage - (heat.source_j - heat.boundary_outflow_j);
            heat.iterations = linear_iterations;
            roundoff += 64.0
                * f64::EPSILON
                * (heat.source_j.abs() + heat.boundary_outflow_j.abs() + true_storage.abs());
            for value in [
                true_storage,
                fixed_correction,
                roundoff,
                heat.closure_j,
                report.solid_stored_energy_change_j,
            ] {
                finite(value)?;
            }
            if report.temperature_defect_k
                <= self.controls.temperature_tolerance_k + 64.0 * f64::EPSILON * largest_temperature
                && report.conductivity_relative_defect
                    <= self.controls.conductivity_relative_tolerance
                && storage_admitted
                && heat.closure_j.abs() <= self.controls.storage_tolerance_j + roundoff
            {
                self.endpoint = endpoint;
                self.records.push(report);
                return Ok((next, heat));
            }
            if iteration == self.controls.maximum_iterations {
                return Err(MaterialBuoyancyError::SolidNotConverged(report));
            }
            if self.controls.temperature_relaxation == 1.0 {
                trial = next;
                applied = endpoint;
            } else {
                for (value, target) in trial.iter_mut().zip(next) {
                    *value = (1.0 - self.controls.temperature_relaxation) * *value
                        + self.controls.temperature_relaxation * target;
                }
                applied = samples(
                    step.domain,
                    self.sources,
                    &trial,
                    |cell| step.is_compact_cell(cell),
                    self.gate,
                )?;
            }
        }
        unreachable!("positive bounded solid iteration accepts or refuses its last endpoint")
    }

    // A rejected linear predictor is not a physical state. Backtrack only the
    // constitutive iterate; the next native solve must still produce an admitted,
    // unconstrained endpoint and pass the source-integral heat residual.
    fn backtrack(
        &self,
        step: &MaterialHeatStep<'_>,
        trial: &[f64],
        target: &[f64],
        failure: MaterialBuoyancyError,
    ) -> Result<(Vec<f64>, Vec<SolidConductivitySample>), MaterialBuoyancyError> {
        let mut fraction = self.controls.temperature_relaxation.min(0.5);
        for _ in 0..48 {
            poll(self.gate)?;
            let candidate: Vec<_> = trial
                .iter()
                .zip(target)
                .map(|(&a, &b)| a + fraction * (b - a))
                .collect();
            if (0..step.domain.cell_count()).all(|cell| {
                step.domain.is_fluid(cell)
                    || step.is_compact_cell(cell)
                    || candidate[cell] == trial[cell]
            }) {
                return Err(failure);
            }
            match samples(
                step.domain,
                self.sources,
                &candidate,
                |cell| step.is_compact_cell(cell),
                self.gate,
            ) {
                Ok(source) => return Ok((candidate, source)),
                Err(
                    error @ MaterialBuoyancyError::Solver(crate::conjugate::ChtError::Cancelled),
                ) => {
                    return Err(error);
                }
                Err(_) => {}
            }
            fraction *= 0.5;
        }
        Err(failure)
    }
}
