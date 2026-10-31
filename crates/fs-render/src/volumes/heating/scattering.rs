//! Multiple-scattering absorption on the native conservative thermal timeline.
//!
//! Only sampled true absorption enters `DepositionSimulation` as a cell source.
//! Ballistic escape, diffuse escape and unresolved radiation retain their own
//! accounts. The transport field is fixed over this observer's accepted time;
//! its Monte Carlo uncertainty is therefore fully correlated across substeps.
//! Energy standard errors scale with duration, not its square root. This is
//! deposited incident energy, not net radiative equilibrium or an LTE cooling law.

use super::{CheckedSum, OpticalHeatingError, finite, invalid, nonnegative};
use crate::volumes::transport::scattering::deposition::{PowerEstimate, ScatteringPowerDeposition};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

/// Explicit acceptance of roundoff, unresolved radiation and sampling error.
/// Statistical quantities are empirical fixed-sample estimates, not guaranteed
/// error bounds. These thresholds must be chosen for the intended thermal use.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatteringHeatingTolerance {
    /// Absolute per-lane power conservation and cell-reduction allowance, W.
    pub power_roundoff_w: f64,
    /// Maximum total mean unresolved power over all spectral lanes, W.
    pub max_unresolved_w: f64,
    /// Maximum sum of per-lane standard errors for either total absorption or
    /// total unresolved power, W. Summation conservatively avoids assuming
    /// independence between arbitrary publicly supplied report lanes.
    pub max_total_standard_error_w: f64,
    /// Maximum sum of per-lane, per-optical-cell absorption standard errors
    /// contributing to one receiving thermal cell, W. This conservative sum
    /// retains possible covariance when several optical cells map to one cell.
    pub max_cell_standard_error_w: f64,
}

/// Accepted incident-radiation accounts integrated over actual thermal time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatteringEnergyBudget<const N: usize> {
    pub duration_s: f64,
    pub incident_j: [f64; N],
    pub ballistic_escape_j: [f64; N],
    pub diffuse_escape_j: [f64; N],
    pub absorbed_j: [f64; N],
    pub unresolved_j: [f64; N],
    /// Estimated error of the fixed deposited-energy estimate, J.
    pub absorbed_standard_error_j: [f64; N],
    pub unresolved_standard_error_j: [f64; N],
    pub residual_j: [f64; N],
    /// Already included once in native `EnergyAudit::external_heat_j`.
    pub thermal_input_j: f64,
}

/// Fixed sampled heat sources; the native observer protocol owns transaction,
/// conduction, caloric storage, latent heat, events and rollback.
#[derive(Clone, Debug)]
pub struct ScatteringHeating<const N: usize> {
    // Candidate thermal advances share this immutable field, including its
    // bounded spatial uncertainty report, instead of cloning it at each trial.
    deposition: Arc<ScatteringPowerDeposition<N>>,
    sources_w: BTreeMap<CellId, f64>,
    thermal_input_w: f64,
    residual_w: [f64; N],
    time_s: f64,
    duration: CheckedSum,
    energy: ScatteringEnergyBudget<N>,
}

impl<const N: usize> ScatteringHeating<N> {
    /// Revalidate the report and map every positively absorbing cell exactly
    /// once to an active native thermal cell. Multiple grid cells may map into
    /// one thermal cell. The caller owns that geometric/control-volume mapping.
    /// Sources remain unchanged until a new observer is prepared at an accepted
    /// endpoint. No incident, escaped or unresolved power is deposited as heat.
    pub fn try_new(
        deposition: ScatteringPowerDeposition<N>,
        cell_bindings: &[(usize, CellId)],
        thermal: &DepositionSimulation<'_>,
        tolerance: ScatteringHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        checkpoint()?;
        if N == 0
            || thermal.cells().is_empty()
            || ![
                tolerance.power_roundoff_w,
                tolerance.max_unresolved_w,
                tolerance.max_total_standard_error_w,
                tolerance.max_cell_standard_error_w,
            ]
            .into_iter()
            .all(nonnegative)
        {
            return Err(invalid(
                "active thermal cells, nonempty packet and finite nonnegative scattering tolerances required",
            ));
        }
        let (lo, hi) = deposition.wavelength_bounds_nm;
        if !(lo.is_finite() && hi.is_finite() && lo > 0.0 && hi > lo)
            || cell_bindings.len() > deposition.cells.len()
        {
            return Err(invalid(
                "invalid spectral domain or excess thermal bindings",
            ));
        }
        let mut bindings = BTreeMap::new();
        for &(index, cell) in cell_bindings {
            checkpoint()?;
            if bindings.insert(index, cell).is_some() {
                return Err(OpticalHeatingError::DuplicateCell(index));
            }
            if !thermal.cells().contains_key(&cell) {
                return Err(OpticalHeatingError::InactiveCell(cell));
            }
        }
        for estimate in [
            &deposition.ballistic_escape,
            &deposition.diffuse_escape,
            &deposition.absorbed,
            &deposition.unresolved,
        ] {
            validate_estimate(estimate)?;
        }
        let mut unresolved = CheckedSum::default();
        let mut total_absorbed_error = CheckedSum::default();
        let mut total_unresolved_error = CheckedSum::default();
        let mut residual_w = [0.0; N];
        for lane in 0..N {
            checkpoint()?;
            let wavelength = deposition.wavelengths_nm[lane];
            let samples = deposition.sampled_histories[lane];
            if !wavelength.is_finite()
                || !(lo..=hi).contains(&wavelength)
                || !nonnegative(deposition.incident_power_w[lane])
                || (samples != 0
                    && (samples < 2 || samples != deposition.config.histories_per_lane))
                || deposition.precision_limited_histories[lane] > samples
                || (samples == 0
                    && (deposition.diffuse_escape.mean_w[lane] != 0.0
                        || deposition.unresolved.mean_w[lane] != 0.0
                        || [
                            &deposition.ballistic_escape,
                            &deposition.diffuse_escape,
                            &deposition.absorbed,
                            &deposition.unresolved,
                        ]
                        .iter()
                        .any(|e| e.standard_error_w[lane] != 0.0)))
            {
                return Err(invalid(
                    "invalid wavelength, source power or sampling account",
                ));
            }
            let mut balance = CheckedSum::default();
            balance.add(deposition.incident_power_w[lane])?;
            for power in [
                deposition.ballistic_escape.mean_w[lane],
                deposition.diffuse_escape.mean_w[lane],
                deposition.absorbed.mean_w[lane],
                deposition.unresolved.mean_w[lane],
            ] {
                balance.add(-power)?;
            }
            residual_w[lane] = balance.value()?;
            if residual_w[lane].abs() > tolerance.power_roundoff_w {
                return Err(invalid("sampled incident power is not conserved"));
            }
            unresolved.add(deposition.unresolved.mean_w[lane])?;
            total_absorbed_error.add(deposition.absorbed.standard_error_w[lane])?;
            total_unresolved_error.add(deposition.unresolved.standard_error_w[lane])?;
        }
        if unresolved.value()? > tolerance.max_unresolved_w
            || total_absorbed_error.value()? > tolerance.max_total_standard_error_w
            || total_unresolved_error.value()? > tolerance.max_total_standard_error_w
        {
            return Err(invalid(
                "scattering truncation or total sampling error exceeds the heating allowance",
            ));
        }
        let mut reported = BTreeSet::new();
        let mut cell_totals = [CheckedSum::default(); N];
        let mut sources = BTreeMap::<CellId, CheckedSum>::new();
        let mut receiving_errors = BTreeMap::<CellId, CheckedSum>::new();
        for cell in &deposition.cells {
            checkpoint()?;
            if !reported.insert(cell.cell_index) {
                return Err(invalid("duplicate cell in multiple-scattering deposition"));
            }
            validate_estimate(&cell.absorbed)?;
            let mut heat = CheckedSum::default();
            let mut error = CheckedSum::default();
            for lane in 0..N {
                if deposition.sampled_histories[lane] == 0
                    && cell.absorbed.standard_error_w[lane] != 0.0
                {
                    return Err(invalid(
                        "deterministic cell lane cannot claim a sampling error",
                    ));
                }
                cell_totals[lane].add(cell.absorbed.mean_w[lane])?;
                heat.add(cell.absorbed.mean_w[lane])?;
                error.add(cell.absorbed.standard_error_w[lane])?;
            }
            let heat = heat.value()?;
            if heat > 0.0 {
                let id = bindings
                    .get(&cell.cell_index)
                    .ok_or(OpticalHeatingError::UnmappedCell(cell.cell_index))?;
                sources.entry(*id).or_default().add(heat)?;
                receiving_errors
                    .entry(*id)
                    .or_default()
                    .add(error.value()?)?;
            }
        }
        if bindings.keys().any(|index| !reported.contains(index)) {
            return Err(invalid(
                "thermal binding names an unreported absorbing cell",
            ));
        }
        let mut absorbed = CheckedSum::default();
        for lane in 0..N {
            if (cell_totals[lane].value()? - deposition.absorbed.mean_w[lane]).abs()
                > tolerance.power_roundoff_w
            {
                return Err(invalid(
                    "spatial absorbed powers disagree with the sampled total",
                ));
            }
            absorbed.add(deposition.absorbed.mean_w[lane])?;
        }
        let mut sources_w = BTreeMap::new();
        let mut thermal_input = CheckedSum::default();
        for (id, sum) in sources {
            checkpoint()?;
            if receiving_errors[&id].value()? > tolerance.max_cell_standard_error_w {
                return Err(invalid(
                    "receiving thermal-cell sampling error exceeds the heating allowance",
                ));
            }
            let power = sum.value()?;
            thermal_input.add(power)?;
            sources_w.insert(id, power);
        }
        let thermal_input_w = thermal_input.value()?;
        if (thermal_input_w - absorbed.value()?).abs() > tolerance.power_roundoff_w {
            return Err(invalid(
                "mapped thermal source disagrees with sampled absorption",
            ));
        }
        let energy = integrate(&deposition, residual_w, thermal_input_w, 0.0)?;
        checkpoint()?;
        Ok(Self {
            deposition: Arc::new(deposition),
            sources_w,
            thermal_input_w,
            residual_w,
            time_s: thermal.time_s(),
            duration: CheckedSum::default(),
            energy,
        })
    }

    /// Original field, sampling errors, source stream and radiation accounts.
    pub fn deposition(&self) -> &ScatteringPowerDeposition<N> {
        &self.deposition
    }
    pub fn sources_w(&self) -> &BTreeMap<CellId, f64> {
        &self.sources_w
    }
    pub fn energy(&self) -> &ScatteringEnergyBudget<N> {
        &self.energy
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }

    fn check_clock(&self, thermal: &DepositionSimulation<'_>) -> Result<(), OpticalHeatingError> {
        if thermal.time_s() != self.time_s {
            Err(invalid(
                "scattering heating and native thermal clocks differ",
            ))
        } else {
            Ok(())
        }
    }
}

impl<const N: usize> DepositionObserver for ScatteringHeating<N> {
    type Error = OpticalHeatingError;

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.check_clock(thermal)?;
        for cell in self.sources_w.keys() {
            checkpoint()?;
            if !thermal.cells().contains_key(cell) {
                return Err(OpticalHeatingError::InactiveCell(*cell));
            }
        }
        Ok(self.sources_w.clone())
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.check_clock(thermal)
    }

    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        if start_s != self.time_s
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || thermal.time_s() <= start_s
        {
            return Err(invalid(
                "scattering energy requires the matching accepted thermal interval",
            ));
        }
        let mut duration = self.duration;
        duration.add(duration_s)?;
        let energy = integrate(
            &self.deposition,
            self.residual_w,
            self.thermal_input_w,
            duration.value()?,
        )?;
        checkpoint()?;
        self.duration = duration;
        self.energy = energy;
        self.time_s = thermal.time_s();
        Ok(())
    }
}

fn validate_estimate<const N: usize>(
    estimate: &PowerEstimate<N>,
) -> Result<(), OpticalHeatingError> {
    if estimate
        .mean_w
        .iter()
        .chain(estimate.standard_error_w.iter())
        .all(|v| nonnegative(*v))
    {
        Ok(())
    } else {
        Err(invalid(
            "sampled power and standard errors must be finite and nonnegative",
        ))
    }
}

fn integrate<const N: usize>(
    report: &ScatteringPowerDeposition<N>,
    residual: [f64; N],
    thermal_w: f64,
    duration_s: f64,
) -> Result<ScatteringEnergyBudget<N>, OpticalHeatingError> {
    let mut energy = ScatteringEnergyBudget {
        duration_s,
        incident_j: [0.0; N],
        ballistic_escape_j: [0.0; N],
        diffuse_escape_j: [0.0; N],
        absorbed_j: [0.0; N],
        unresolved_j: [0.0; N],
        absorbed_standard_error_j: [0.0; N],
        unresolved_standard_error_j: [0.0; N],
        residual_j: [0.0; N],
        thermal_input_j: finite(thermal_w * duration_s)?,
    };
    for lane in 0..N {
        energy.incident_j[lane] = finite(report.incident_power_w[lane] * duration_s)?;
        energy.ballistic_escape_j[lane] =
            finite(report.ballistic_escape.mean_w[lane] * duration_s)?;
        energy.diffuse_escape_j[lane] = finite(report.diffuse_escape.mean_w[lane] * duration_s)?;
        energy.absorbed_j[lane] = finite(report.absorbed.mean_w[lane] * duration_s)?;
        energy.unresolved_j[lane] = finite(report.unresolved.mean_w[lane] * duration_s)?;
        energy.absorbed_standard_error_j[lane] =
            finite(report.absorbed.standard_error_w[lane] * duration_s)?;
        energy.unresolved_standard_error_j[lane] =
            finite(report.unresolved.standard_error_w[lane] * duration_s)?;
        energy.residual_j[lane] = finite(residual[lane] * duration_s)?;
    }
    Ok(energy)
}
