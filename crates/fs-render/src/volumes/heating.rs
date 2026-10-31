//! Fixed incident-beam absorption on the existing conservative thermal timeline.
//!
//! This bridge consumes the exact cell powers from `deposit_beam_power`; it
//! introduces no absorption, heat-capacity, geometry or time-integration law.
//! Only true absorption heats the mapped cells. Directed transmission and
//! scattered power retain separate accounts and must be transported by their
//! original radiation owner if subsequent absorption is required.
//!
//! Optics, ray geometry and incident spectral powers are frozen for this
//! observer's lifetime. Recompute deposition and construct a new observer at an
//! accepted thermal endpoint when those inputs change. Visible-band beam
//! heating does not infer infrared emission/cooling, thermal equilibrium,
//! multiple scattering, fluorescence or a temperature-independent material law.

use super::transport::deposition::BeamPowerDeposition;
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN};
use fs_blake3::ContentHash;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
};
use std::collections::BTreeMap;

/// Retained SI power accounting for the original wavelength packet.
/// Each lane is already a quadrature-weighted power in W, not spectral radiance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OpticalPowerBudget<const N: usize> {
    /// Immutable source medium identity supplied by the transport owner.
    pub medium_identity: ContentHash,
    /// Original vacuum wavelengths, in nm.
    pub wavelengths_nm: [f64; N],
    /// Power supplied to the original directed segment, W.
    pub incident_w: [f64; N],
    /// Power continuing along the original direction, W.
    pub transmitted_w: [f64; N],
    /// Power removed by true absorption, W.
    pub absorbed_w: [f64; N],
    /// Power redirected out of the original beam, W; this is not heat.
    pub scattered_w: [f64; N],
    /// Incident minus transmitted, absorbed and scattered power, W.
    pub residual_w: [f64; N],
    /// Actual sum of prepared thermal-cell source rates, W.
    pub thermal_input_w: f64,
}

/// Cumulative physical energy over the native accepted thermal substeps.
/// The thermal input is already included in `EnergyAudit::external_heat_j`;
/// it must not be added to that audit again.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OpticalEnergyBudget<const N: usize> {
    /// Sum of the actual accepted thermal substep durations, s.
    pub duration_s: f64,
    /// Original source energy, J.
    pub incident_j: [f64; N],
    /// Still-directed outgoing energy, J.
    pub transmitted_j: [f64; N],
    /// True absorbed energy, J.
    pub absorbed_j: [f64; N],
    /// Redirected outgoing energy, J; no thermal credit is assigned.
    pub scattered_j: [f64; N],
    /// Original per-lane power residual integrated over this duration, J.
    pub residual_j: [f64; N],
    /// Integrated rate actually supplied to the mapped thermal cells, J.
    pub thermal_input_j: f64,
}

/// No refused construction or native advance publishes a partial heating state.
#[derive(Clone, Debug, PartialEq)]
pub enum OpticalHeatingError {
    /// Invalid report, mapping, clock, tolerance or arithmetic.
    Invalid(&'static str),
    /// A cell with positive true absorption has no receiving thermal binding.
    UnmappedCell(usize),
    /// A grid cell has more than one thermal assignment.
    DuplicateCell(usize),
    /// A binding names a thermal cell that is not active.
    InactiveCell(CellId),
    /// Source-owned thermal domain, conservation, work-budget or cancellation error.
    Thermal(DepositionError),
}

impl From<DepositionError> for OpticalHeatingError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl core::fmt::Display for OpticalHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "optical heating refused: {self:?}")
    }
}
impl std::error::Error for OpticalHeatingError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Thermal(error) => Some(error),
            _ => None,
        }
    }
}

/// Prepared absorption-only source on the original thermal cells.
///
/// Construction requires an already activated thermal snapshot: process due
/// births with the native `advance_to` first. Sparse bindings identify the
/// original DDA cell indices. Every positively absorbing cell needs exactly
/// one binding; several grid cells may contribute to the same thermal cell.
/// Geometry/control-volume compatibility is the caller's explicit assignment.
///
/// Pass this observer to `DepositionSimulation::advance_observed`, or compose
/// it with another physical owner using `fs_thermal::observers::CoupledObservers`.
/// Native candidate cloning makes heat, optical accounts and other observers
/// transactional together. Authored boundary sources are never overwritten.
/// Prepared source storage is O(receiving thermal cells); each substep needs no
/// optical traversal, spectral interpolation or retained grid-cell history.
#[derive(Clone, Debug)]
pub struct OpticalHeating<const N: usize> {
    power: OpticalPowerBudget<N>,
    sources_w: BTreeMap<CellId, f64>,
    time_s: f64,
    duration: CheckedSum,
    energy: OpticalEnergyBudget<N>,
}

impl<const N: usize> OpticalHeating<N> {
    /// Validate the transport report and prepare a fixed physical heat source.
    ///
    /// `power_tolerance_w` is an absolute allowance for each per-lane
    /// conservation and cell-reduction check, and for the total mapped heat.
    /// It is a roundoff admission bound, not a transport/model error estimate.
    /// Public report fields are revalidated against those allowances; negative,
    /// nonfinite, unbalanced or incomplete transport reports refuse.
    /// Bindings may cover zero-absorption cells, but must name reported cells.
    /// Cancellation is checked during admission before any source is published.
    pub fn try_new(
        deposition: BeamPowerDeposition<N>,
        cell_bindings: &[(usize, CellId)],
        thermal: &DepositionSimulation<'_>,
        power_tolerance_w: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        checkpoint()?;
        if N == 0 || !nonnegative(power_tolerance_w) || thermal.cells().is_empty() {
            return Err(invalid(
                "nonempty packet, active thermal cells and finite nonnegative power tolerance required",
            ));
        }
        if cell_bindings.len() > deposition.cells.len() {
            return Err(invalid(
                "bindings must refer only to cells in the deposition report",
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
        let mut residual_w = [0.0; N];
        for lane in 0..N {
            checkpoint()?;
            if !deposition.wavelengths_nm[lane].is_finite()
                || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&deposition.wavelengths_nm[lane])
                || ![
                    deposition.incident_power_w[lane],
                    deposition.transmitted_power_w[lane],
                    deposition.absorbed_power_w[lane],
                    deposition.scattered_power_w[lane],
                ]
                .into_iter()
                .all(nonnegative)
            {
                return Err(invalid(
                    "finite visible wavelengths and nonnegative finite beam powers required",
                ));
            }
            let mut balance = CheckedSum::default();
            balance.add(deposition.incident_power_w[lane])?;
            balance.add(-deposition.transmitted_power_w[lane])?;
            balance.add(-deposition.absorbed_power_w[lane])?;
            balance.add(-deposition.scattered_power_w[lane])?;
            residual_w[lane] = balance.value()?;
            if residual_w[lane].abs() > power_tolerance_w {
                return Err(invalid("incident beam power is not conserved"));
            }
        }
        let mut absorbed = [CheckedSum::default(); N];
        let mut scattered = [CheckedSum::default(); N];
        let mut sources = BTreeMap::<CellId, CheckedSum>::new();
        let mut reported = std::collections::BTreeSet::new();
        for cell in &deposition.cells {
            checkpoint()?;
            if !reported.insert(cell.cell_index) {
                return Err(invalid("duplicate cell in beam deposition report"));
            }
            let mut heat = CheckedSum::default();
            for lane in 0..N {
                let a = cell.absorbed_power_w[lane];
                let s = cell.scattered_power_w[lane];
                if !nonnegative(a) || !nonnegative(s) {
                    return Err(invalid(
                        "cell absorption and scattering powers must be finite and nonnegative",
                    ));
                }
                absorbed[lane].add(a)?;
                scattered[lane].add(s)?;
                heat.add(a)?;
            }
            let heat = heat.value()?;
            if heat > 0.0 {
                let id = bindings
                    .get(&cell.cell_index)
                    .ok_or(OpticalHeatingError::UnmappedCell(cell.cell_index))?;
                sources.entry(*id).or_default().add(heat)?;
            }
        }
        if bindings.keys().any(|index| !reported.contains(index)) {
            return Err(invalid("thermal binding names an unreported grid cell"));
        }
        let mut declared_heat = CheckedSum::default();
        for lane in 0..N {
            if (absorbed[lane].value()? - deposition.absorbed_power_w[lane]).abs()
                > power_tolerance_w
                || (scattered[lane].value()? - deposition.scattered_power_w[lane]).abs()
                    > power_tolerance_w
            {
                return Err(invalid(
                    "cell-resolved powers disagree with the beam totals",
                ));
            }
            declared_heat.add(deposition.absorbed_power_w[lane])?;
        }
        let mut sources_w = BTreeMap::new();
        let mut thermal_input = CheckedSum::default();
        for (cell, heat) in sources {
            checkpoint()?;
            let heat = heat.value()?;
            thermal_input.add(heat)?;
            sources_w.insert(cell, heat);
        }
        let thermal_input_w = thermal_input.value()?;
        if (thermal_input_w - declared_heat.value()?).abs() > power_tolerance_w {
            return Err(invalid(
                "mapped thermal input disagrees with true absorption",
            ));
        }
        let power = OpticalPowerBudget {
            medium_identity: deposition.medium_identity,
            wavelengths_nm: deposition.wavelengths_nm,
            incident_w: deposition.incident_power_w,
            transmitted_w: deposition.transmitted_power_w,
            absorbed_w: deposition.absorbed_power_w,
            scattered_w: deposition.scattered_power_w,
            residual_w,
            thermal_input_w,
        };
        let energy = power.integrate(0.0)?;
        checkpoint()?;
        Ok(Self {
            power,
            sources_w,
            time_s: thermal.time_s(),
            duration: CheckedSum::default(),
            energy,
        })
    }

    /// Immutable original beam budget and mapped heat rate.
    pub fn power(&self) -> &OpticalPowerBudget<N> {
        &self.power
    }

    /// Accepted cumulative accounts since this observer was constructed.
    pub fn energy(&self) -> &OpticalEnergyBudget<N> {
        &self.energy
    }

    /// Thermal clock to which this candidate's optical accounting is bound.
    pub fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Prepared true-absorption watts, keyed by stable thermal cell identity.
    pub fn sources_w(&self) -> &BTreeMap<CellId, f64> {
        &self.sources_w
    }

    fn check_clock(&self, thermal: &DepositionSimulation<'_>) -> Result<(), OpticalHeatingError> {
        if thermal.time_s() != self.time_s {
            return Err(invalid("optical and thermal accepted clocks differ"));
        }
        Ok(())
    }
}

impl<const N: usize> DepositionObserver for OpticalHeating<N> {
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
                "optical accounting requires the matching accepted thermal interval",
            ));
        }
        let mut duration = self.duration;
        duration.add(duration_s)?;
        let energy = self.power.integrate(duration.value()?)?;
        checkpoint()?;
        self.duration = duration;
        self.energy = energy;
        self.time_s = thermal.time_s();
        Ok(())
    }
}

impl<const N: usize> OpticalPowerBudget<N> {
    fn integrate(&self, duration_s: f64) -> Result<OpticalEnergyBudget<N>, OpticalHeatingError> {
        let mut energy = OpticalEnergyBudget {
            duration_s,
            incident_j: [0.0; N],
            transmitted_j: [0.0; N],
            absorbed_j: [0.0; N],
            scattered_j: [0.0; N],
            residual_j: [0.0; N],
            thermal_input_j: finite(self.thermal_input_w * duration_s)?,
        };
        for lane in 0..N {
            energy.incident_j[lane] = finite(self.incident_w[lane] * duration_s)?;
            energy.transmitted_j[lane] = finite(self.transmitted_w[lane] * duration_s)?;
            energy.absorbed_j[lane] = finite(self.absorbed_w[lane] * duration_s)?;
            energy.scattered_j[lane] = finite(self.scattered_w[lane] * duration_s)?;
            energy.residual_j[lane] = finite(self.residual_w[lane] * duration_s)?;
        }
        Ok(energy)
    }
}

fn invalid(reason: &'static str) -> OpticalHeatingError {
    OpticalHeatingError::Invalid(reason)
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn finite(value: f64) -> Result<f64, OpticalHeatingError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable optical heat account"))
    }
}

/// Compensated reduction, also retaining small accepted durations at long runs.
#[derive(Clone, Copy, Debug, Default)]
struct CheckedSum {
    sum: f64,
    correction: f64,
}
impl CheckedSum {
    fn add(&mut self, value: f64) -> Result<(), OpticalHeatingError> {
        let next = finite(self.sum + value)?;
        let error = if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.correction = finite(self.correction + error)?;
        self.sum = next;
        self.value()?;
        Ok(())
    }
    fn value(self) -> Result<f64, OpticalHeatingError> {
        finite(self.sum + self.correction)
    }
}
