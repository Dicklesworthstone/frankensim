//! Scheduled adaptive electrothermal evolution on the original caloric clock.

use super::{
    CaloricElectrothermalConfig, CaloricElectrothermalReport, admit, advance_with_heat, checkpoint,
    exceeds, number, refusal,
};
use crate::conjugate::tensor_conduction::{
    CaloricAdaptiveForcing, CaloricAdaptiveTrial, advance_caloric_adaptive_with,
    advance_radiating_sdirk2_prepared, caloric_stage_clock, embedded_ratio,
};
use crate::conjugate::{
    ChtError, ElectricSetup, ElectricSolution, ElectricalContact, RadiationConfig, SolidMaterial,
    TENSOR_SDIRK2_GAMMA, TensorCaloricAdaptiveConfig, TensorCaloricAdaptiveError,
    TensorCaloricAdaptiveSolution, TensorCaloricForcing, TensorCaloricForcingSample,
    TensorCaloricRadiationReport, TensorCaloricStorage, TensorConductivityLaw, ThermalSetup, Voxel,
    VoxelDomain, advance_caloric_tensor_conduction_sdirk2, poll,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use std::borrow::Cow;

/// Original physical thermal/electrical data at one native stage sample.
#[derive(Clone)]
pub struct CaloricElectrothermalForcing<'a> {
    pub thermal: TensorCaloricForcing<'a>,
    pub electric: Cow<'a, ElectricSetup>,
    pub contacts: Cow<'a, [ElectricalContact]>,
}

impl<'a> CaloricElectrothermalForcing<'a> {
    pub fn borrowed(
        thermal: &'a ThermalSetup,
        prescribed: &'a [(usize, f64)],
        electric: &'a ElectricSetup,
        contacts: &'a [ElectricalContact],
    ) -> Self {
        Self {
            thermal: TensorCaloricForcing::borrowed(thermal, prescribed),
            electric: Cow::Borrowed(electric),
            contacts: Cow::Borrowed(contacts),
        }
    }
}

impl CaloricAdaptiveForcing for CaloricElectrothermalForcing<'_> {
    fn setup(&self) -> &ThermalSetup {
        self.thermal.setup.as_ref()
    }

    fn prescribed(&self) -> &[(usize, f64)] {
        self.thermal.prescribed_specific_energy_j_kg.as_ref()
    }
}

/// Original per-proposal coupling controls and native adaptive time controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricElectrothermalAdaptiveConfig {
    /// Thermal limits span both stages of each source proposal. Coupling passes
    /// span one attempted interval, including an interval later rejected in time.
    pub step: CaloricElectrothermalConfig,
    pub initial_dt_s: f64,
    pub minimum_dt_s: f64,
    pub maximum_dt_s: f64,
    pub maximum_attempts: usize,
    pub energy_absolute_tolerance_j_kg: f64,
    /// Local embedded discrepancy in APPLIED Joule energy [J], positive.
    /// Includes prescribed cells even when their stored energy never changes.
    pub source_absolute_tolerance_j: f64,
    /// Relative scales use original energy changes and positive deposited heat.
    pub relative_tolerance: f64,
    /// Whole-window allowance for native closure AND the sum of absolute
    /// accepted coupling defects. Per-step allowances do not reset this budget.
    pub energy_tolerance_j: f64,
}

impl Default for CaloricElectrothermalAdaptiveConfig {
    fn default() -> Self {
        let time = TensorCaloricAdaptiveConfig::default();
        Self {
            step: CaloricElectrothermalConfig {
                thermal: time.step,
                ..CaloricElectrothermalConfig::default()
            },
            initial_dt_s: time.initial_dt_s,
            minimum_dt_s: time.minimum_dt_s,
            maximum_dt_s: time.maximum_dt_s,
            maximum_attempts: time.maximum_attempts,
            energy_absolute_tolerance_j_kg: time.energy_absolute_tolerance_j_kg,
            source_absolute_tolerance_j: 1e-3,
            relative_tolerance: time.relative_tolerance,
            energy_tolerance_j: time.energy_tolerance_j,
        }
    }
}

impl CaloricElectrothermalAdaptiveConfig {
    fn time_controls(self) -> TensorCaloricAdaptiveConfig {
        TensorCaloricAdaptiveConfig {
            step: self.step.thermal,
            initial_dt_s: self.initial_dt_s,
            minimum_dt_s: self.minimum_dt_s,
            maximum_dt_s: self.maximum_dt_s,
            maximum_attempts: self.maximum_attempts,
            energy_absolute_tolerance_j_kg: self.energy_absolute_tolerance_j_kg,
            relative_tolerance: self.relative_tolerance,
            energy_tolerance_j: self.energy_tolerance_j,
        }
    }
}

/// Completed numerical work and compensated accepted physical electrical heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricElectrothermalAdaptiveReport {
    /// Complete coupled trials include temporal rejections. Incomplete failures
    /// remain bounded by attempts times the explicit per-attempt work ceilings.
    pub completed_trial_passes: usize,
    pub completed_trial_current_solves: usize,
    pub applied_electrode_energy_j: f64,
    pub applied_joule_energy_j: f64,
    pub thermostat_joule_energy_j: f64,
    pub source_mismatch_j: f64,
    pub absolute_energy_defect_j: f64,
}

/// One accepted original-energy timeline, with small aligned coupling receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct CaloricElectrothermalAdaptiveSolution {
    pub thermal: TensorCaloricAdaptiveSolution,
    /// Exactly one accepted report for each thermal.records interval.
    pub coupling: Vec<CaloricElectrothermalReport>,
    /// Empty without radiation; otherwise aligned with thermal.records.
    pub radiation: Vec<[TensorCaloricRadiationReport; 2]>,
    /// Already included in the native boundary ledger; positive outward [J].
    pub radiated_energy_j: f64,
    pub completed_trial_radiation_solves: usize,
    /// Only the final accepted interval's two stage fields are retained.
    pub applied_electric: [ElectricSolution; 2],
    pub endpoint_electric: [ElectricSolution; 2],
    pub report: CaloricElectrothermalAdaptiveReport,
}

/// Evolve spatial Joule heating, original calorics and optional gray radiation.
///
/// Uses the original adaptive SDIRK2 clock, event left limits, prescribed-mask
/// admission, embedded controller and compensated heat ledger. Both forcing and
/// resistivity receive the same original stage sample, including its left-limit
/// marker. Callbacks may maintain caches or diagnostics but must not publish
/// rejected physical history. Source laws and their query budgets persist
/// across all coupling proposals and accepted/rejected temporal attempts.
///
/// The error gate includes original stored-energy increments AND the embedded
/// discrepancy of each cell's applied Joule energy, including thermostats.
/// Constant-temperature latent storage and changing heat through a prescribed
/// conductor therefore remain visible to time adaptation. These are local
/// asymptotic estimates, not global error certificates or quadrature guarantees
/// for every independent external thermal channel. Declare source jumps as
/// transitions; sharp material corners can reduce temporal order.
///
/// Radiation geometry is built once for the whole window. Fixed mass/geometry,
/// stationary DC, original gray surfaces and fixed radiative surroundings are
/// retained. The same physical radiation must not also be authored as a sink.
/// Any late source/work/energy/cancellation refusal discards the whole private
/// prefix. Only accepted heat and small coupling receipts are published.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_caloric_tensor_electrothermal_adaptive<'f>(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    conductivity_laws: &[&dyn TensorConductivityLaw],
    storage: TensorCaloricStorage<'_>,
    start_time_s: f64,
    target_time_s: f64,
    transition_times_s: &[f64],
    config: &CaloricElectrothermalAdaptiveConfig,
    radiation: Option<&RadiationConfig>,
    forcing: impl FnMut(
        TensorCaloricForcingSample,
    ) -> Result<CaloricElectrothermalForcing<'f>, ChtError>,
    mut resistivity: impl FnMut(
        TensorCaloricForcingSample,
        usize,
        Voxel,
        f64,
        f64,
    ) -> Result<Option<f64>, ChtError>,
    gate: &CancelGate,
) -> Result<CaloricElectrothermalAdaptiveSolution, TensorCaloricAdaptiveError> {
    poll(gate)?;
    admit_window(
        domain.cell_count(),
        start_time_s,
        target_time_s,
        config,
        radiation,
    )?;
    let time_controls = config.time_controls();
    let mut completed_trial_passes = 0usize;
    let mut completed_trial_current_solves = 0usize;
    let mut final_fields = None;
    let run = advance_caloric_adaptive_with(
        domain,
        solids,
        conductivity_laws,
        storage,
        start_time_s,
        target_time_s,
        transition_times_s,
        &time_controls,
        radiation,
        forcing,
        |data, samples, trial_storage, time, dt, radiation| {
            let (a, stage_times) = caloric_stage_clock(time, dt)
                .ok_or_else(|| refusal("adaptive electrothermal stage clock is unrepresentable"))?;
            let mut completed_iterations = 0usize;
            let mut completed_radiation_solves = 0usize;
            let mut accepted_radiation = None;
            let prescribed = [data[0].prescribed(), data[1].prescribed()];
            let coupled = advance_with_heat(
                domain,
                solids,
                [data[0].electric.as_ref(), data[1].electric.as_ref()],
                [data[0].contacts.as_ref(), data[1].contacts.as_ref()],
                [data[0].setup(), data[1].setup()],
                conductivity_laws,
                trial_storage,
                time,
                dt,
                stage_times,
                [1.0 - TENSOR_SDIRK2_GAMMA, TENSOR_SDIRK2_GAMMA],
                &config.step,
                |stage_time, cell, occupancy, energy, temperature| {
                    let stage = usize::from(stage_time == stage_times[1]);
                    resistivity(samples[stage], cell, occupancy, energy, temperature)
                },
                |heated| {
                    let (step, reports) = match radiation {
                        Some((radiation, geometry)) => advance_radiating_sdirk2_prepared(
                            domain,
                            solids,
                            heated,
                            conductivity_laws,
                            trial_storage,
                            prescribed,
                            time,
                            dt,
                            &config.step.thermal,
                            radiation,
                            Some(geometry),
                            gate,
                        )
                        .map(|(step, reports)| (step, Some(reports))),
                        None => advance_caloric_tensor_conduction_sdirk2(
                            domain,
                            solids,
                            heated,
                            conductivity_laws,
                            trial_storage,
                            prescribed,
                            time,
                            dt,
                            &config.step.thermal,
                            gate,
                        )
                        .map(|step| (step, None)),
                    }?;
                    completed_iterations = work_add(completed_iterations, step.record.iterations)?;
                    if let Some(reports) = &reports {
                        for report in reports {
                            completed_radiation_solves =
                                work_add(completed_radiation_solves, report.iterations)?;
                        }
                    }
                    accepted_radiation = reports;
                    Ok(step)
                },
                gate,
            )?;
            completed_trial_passes = work_add(completed_trial_passes, coupled.report.passes)?;
            completed_trial_current_solves = work_add(
                completed_trial_current_solves,
                coupled.report.current_solves,
            )?;
            let mut source_error = 0.0_f64;
            let duration = Dd::from_f64(dt);
            let weights = [1.0 - TENSOR_SDIRK2_GAMMA, TENSOR_SDIRK2_GAMMA];
            for cell in 0..domain.cell_count() {
                checkpoint(gate, cell)?;
                let powers = [
                    Dd::from_f64(coupled.applied_electric[0].joule_w[cell]),
                    Dd::from_f64(coupled.applied_electric[1].joule_w[cell]),
                ];
                let first = number(Dd::from_f64(a) * powers[0])?;
                let second = number(
                    duration
                        * (Dd::from_f64(weights[0]) * powers[0]
                            + Dd::from_f64(weights[1]) * powers[1]),
                )?;
                source_error = source_error.max(
                    embedded_ratio(
                        first,
                        second,
                        a / dt,
                        config.source_absolute_tolerance_j,
                        config.relative_tolerance,
                    )
                    .ok_or_else(|| refusal("unrepresentable embedded Joule-energy comparison"))?,
                );
            }
            // Private diagnostic state may include a rejected completed trial.
            // Successful completion always ends on an accepted trial; every
            // late refusal drops these fields with the entire private window.
            final_fields = Some((coupled.applied_electric, coupled.endpoint_electric));
            Ok(CaloricAdaptiveTrial {
                thermal: coupled.thermal,
                radiation: accepted_radiation,
                receipt: coupled.report,
                completed_iterations,
                completed_radiation_solves,
                additional_normalized_error: source_error,
            })
        },
        gate,
    )?;
    let mut energy = [Dd::ZERO; 5];
    for (index, receipt) in run.receipts.iter().enumerate() {
        checkpoint(gate, index)?;
        for (sum, value) in energy.iter_mut().zip([
            receipt.applied_electrode_energy_j,
            receipt.applied_joule_energy_j,
            receipt.thermostat_joule_energy_j,
            receipt.source_mismatch_j,
            receipt.absolute_energy_defect_j,
        ]) {
            *sum = *sum + Dd::from_f64(value);
        }
    }
    if exceeds(energy[4], config.energy_tolerance_j)? {
        return Err(refusal(&format!(
            "adaptive accepted coupling defect {} J exceeds {} J whole-window allowance",
            number(energy[4])?,
            config.energy_tolerance_j,
        ))
        .into());
    }
    let report = CaloricElectrothermalAdaptiveReport {
        completed_trial_passes,
        completed_trial_current_solves,
        applied_electrode_energy_j: number(energy[0])?,
        applied_joule_energy_j: number(energy[1])?,
        thermostat_joule_energy_j: number(energy[2])?,
        source_mismatch_j: number(energy[3])?,
        absolute_energy_defect_j: number(energy[4])?,
    };
    let (applied_electric, endpoint_electric) = final_fields
        .ok_or_else(|| refusal("accepted adaptive window is missing its final current fields"))?;
    poll(gate)?;
    Ok(CaloricElectrothermalAdaptiveSolution {
        thermal: run.thermal.thermal,
        coupling: run.receipts,
        radiation: run.thermal.radiation,
        radiated_energy_j: run.thermal.radiated_energy_j,
        completed_trial_radiation_solves: run.thermal.completed_trial_radiation_solves,
        applied_electric,
        endpoint_electric,
        report,
    })
}

fn work_add(previous: usize, next: usize) -> Result<usize, ChtError> {
    previous
        .checked_add(next)
        .ok_or_else(|| refusal("adaptive electrothermal completed work counter overflow"))
}

fn admit_window(
    cells: usize,
    start: f64,
    target: f64,
    config: &CaloricElectrothermalAdaptiveConfig,
    radiation: Option<&RadiationConfig>,
) -> Result<(), ChtError> {
    admit(&config.step, start, target - start, cells)?;
    if embedded_ratio(
        0.0,
        0.0,
        TENSOR_SDIRK2_GAMMA,
        config.source_absolute_tolerance_j,
        config.relative_tolerance,
    )
    .is_none()
    {
        return Err(refusal("invalid adaptive Joule-energy accuracy controls"));
    }
    let proposals = config
        .maximum_attempts
        .checked_mul(config.step.maximum_passes)
        .ok_or_else(|| refusal("adaptive coupling proposal bound is unrepresentable"))?;
    for cap in [
        config.step.thermal.max_newton_iterations,
        config.step.thermal.max_backtracks,
        config.step.thermal.linear.max_iterations,
        radiation.map_or(0, |value| value.max_iterations),
    ] {
        proposals
            .checked_mul(cap)
            .ok_or_else(|| refusal("adaptive native work bound is unrepresentable"))?;
    }
    config
        .step
        .maximum_passes
        .checked_add(1)
        .and_then(|count| count.checked_mul(2))
        .and_then(|count| count.checked_mul(config.maximum_attempts))
        .and_then(|count| count.checked_mul(cells))
        .ok_or_else(|| refusal("adaptive current callback bound is unrepresentable"))?;
    Ok(())
}
