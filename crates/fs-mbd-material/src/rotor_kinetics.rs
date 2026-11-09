//! Prescribed material temperatures on the original rotor Maxwell memories.
//!
//! Existing WLF/Arrhenius calibrations supply `tau(T)/tau(T_ref)` to the native
//! rotor cycle owner. Every branch retains its original stiffness, strain and
//! stress; temperature changes its dissipative clock without a storage jump.
//! Temperatures are held on the caller's coefficient segments. Approximating a
//! smooth thermal history requires independent segment refinement.
//!
//! This adapter changes Maxwell times only. The rotor's independently authored
//! C may include Kelvin viscosity, unrelated dampers or active bearings; it is
//! never implicitly scaled. No complete thermorheological material, thermal
//! feedback, expansion, heat capacity or temperature-dependent modulus law is
//! inferred. Native branch heat remains in original branch order for an explicit
//! thermal coupling owner to route.

use fs_material::visco::ViscoError;
use fs_material::visco::temperature::TemperatureShift;
use fs_mbd::rotordynamics::RotorError;
use fs_mbd::rotordynamics::periodic::{
    MAX_ROTOR_PERIODIC_SEGMENTS, PreparedRotorMaxwellCycle, PreparedRotorMaxwellPowerCycle,
    RotorFloquetBudget, RotorMaxwellFloquetResponse, RotorMaxwellKinetics, RotorPeriodicForcing,
    RotorPeriodicSegment,
};
use fs_mbd::rotordynamics::relaxation::{MAX_ROTOR_MAXWELL_BRANCHES, MaxwellDampingBranch};

/// Admission of an entire prescribed-temperature schedule; no partial output.
#[derive(Clone, Debug, PartialEq)]
pub enum RotorTemperatureError {
    Invalid(&'static str),
    /// The original calibration refused this exact segment/branch sample.
    Shift {
        segment: usize,
        branch: usize,
        source: ViscoError,
    },
    /// Native clock representability or cancellation refusal, retained intact.
    Rotor(RotorError),
}

impl core::fmt::Display for RotorTemperatureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(why) => write!(f, "rotor temperature schedule refused: {why}"),
            Self::Shift {
                segment,
                branch,
                source,
            } => {
                write!(
                    f,
                    "rotor temperature segment {segment}, branch {branch}: {source}"
                )
            }
            Self::Rotor(source) => source.fmt(f),
        }
    }
}
impl std::error::Error for RotorTemperatureError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Invalid(_) => None,
            Self::Shift { source, .. } => Some(source),
            Self::Rotor(source) => Some(source),
        }
    }
}
impl From<RotorError> for RotorTemperatureError {
    fn from(source: RotorError) -> Self {
        Self::Rotor(source)
    }
}

/// Immutable original branch bank and its prescribed absolute temperatures.
///
/// The caller binds one existing calibration to each original branch. Its
/// `relaxation_time_s` must be the time at that calibration's reference
/// temperature, including any original constant reference multiplier. This
/// association and calibration accuracy are caller declarations, not inferred
/// from material names. Different branches may have different calibrations.
///
/// All solver methods reuse this exact borrowed bank. Returned motion plans own
/// their maps; power plans borrow only the original models and storage laws.
/// Schedule storage is bounded by the native segment/branch caps, with two f64
/// values per sample, plus row metadata. Solver work uses its ordinary explicit
/// [`RotorFloquetBudget`].
#[derive(Clone, Debug)]
pub struct RotorTemperatureKinetics<'a> {
    branches: &'a [MaxwellDampingBranch],
    temperature_shifts: &'a [TemperatureShift],
    temperatures_k: Vec<Vec<f64>>,
    relaxation_time_scales: Vec<Vec<f64>>,
}

impl<'a> RotorTemperatureKinetics<'a> {
    /// Evaluate one temperature per original [segment][branch].
    ///
    /// Uses the existing calibrated law without extrapolation or refitting.
    /// Every ratio is evaluated against its original reference, never a prior
    /// segment's shifted time. Native admission checks effective times AND
    /// inverse times. Physical strain/stiffness and model admission remain with
    /// the downstream rotor owner. Cancellation is checked before allocation,
    /// per sample and before publication; failure leaves all inputs untouched.
    pub fn try_new(
        branches: &'a [MaxwellDampingBranch],
        temperature_shifts: &'a [TemperatureShift],
        temperatures_k: &[Vec<f64>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RotorTemperatureError> {
        poll(&mut cancelled)?;
        if temperatures_k.is_empty()
            || temperatures_k.len() > MAX_ROTOR_PERIODIC_SEGMENTS
            || branches.len() > MAX_ROTOR_MAXWELL_BRANCHES
            || temperature_shifts.len() != branches.len()
        {
            return Err(RotorTemperatureError::Invalid(
                "bounded segments and one original calibration per branch required",
            ));
        }
        for row in temperatures_k {
            poll(&mut cancelled)?;
            if row.len() != branches.len() {
                return Err(RotorTemperatureError::Invalid(
                    "one temperature per original branch required on every segment",
                ));
            }
        }
        let mut samples = Vec::with_capacity(temperatures_k.len());
        let mut scales = Vec::with_capacity(temperatures_k.len());
        for (segment, row) in temperatures_k.iter().enumerate() {
            poll(&mut cancelled)?;
            let mut ratios = Vec::with_capacity(branches.len());
            for (branch, (&temperature, shift)) in row.iter().zip(temperature_shifts).enumerate() {
                poll(&mut cancelled)?;
                ratios.push(
                    shift
                        .evaluate(temperature)
                        .map_err(|source| RotorTemperatureError::Shift {
                            segment,
                            branch,
                            source,
                        })?
                        .factor(),
                );
            }
            samples.push(row.clone());
            scales.push(ratios);
        }
        RotorMaxwellKinetics {
            relaxation_time_scales: &scales,
        }
        .validate(samples.len(), branches, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(Self {
            branches,
            temperature_shifts,
            temperatures_k: samples,
            relaxation_time_scales: scales,
        })
    }

    /// The unchanged physical branch identities and reference times.
    pub fn branches(&self) -> &'a [MaxwellDampingBranch] {
        self.branches
    }

    /// Original calibrations in that same branch order.
    pub fn temperature_shifts(&self) -> &'a [TemperatureShift] {
        self.temperature_shifts
    }

    /// Admitted absolute temperatures [K], in [segment][branch] order.
    pub fn temperatures_k(&self) -> &[Vec<f64>] {
        &self.temperatures_k
    }

    /// Borrow the evaluated schedule without allocating or evaluating again.
    pub fn as_kinetics(&self) -> RotorMaxwellKinetics<'_> {
        RotorMaxwellKinetics {
            relaxation_time_scales: &self.relaxation_time_scales,
        }
    }

    /// Complete native homogeneous Floquet map, retaining every material mode.
    pub fn floquet(
        &self,
        spin_rad_s: f64,
        segments: &[RotorPeriodicSegment<'_>],
        relative_generator_residual_tolerance: f64,
        budget: RotorFloquetBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<RotorMaxwellFloquetResponse, RotorError> {
        self.as_kinetics().floquet(
            spin_rad_s,
            segments,
            self.branches,
            relative_generator_residual_tolerance,
            budget,
            cancelled,
        )
    }

    /// Native repeated-cycle motion and periodic orbits on the original bank.
    /// Loading phase follows full cycle time, including temperature switches.
    pub fn prepare_cycle(
        &self,
        spin_rad_s: f64,
        segments: &[RotorPeriodicSegment<'_>],
        forcing: RotorPeriodicForcing<'_>,
        relative_generator_residual_tolerance: f64,
        budget: RotorFloquetBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedRotorMaxwellCycle, RotorError> {
        self.as_kinetics().prepare_cycle(
            spin_rad_s,
            segments,
            self.branches,
            forcing,
            relative_generator_residual_tolerance,
            budget,
            cancelled,
        )
    }

    /// Native source work, storage changes and individual branch heat using the
    /// same held temperatures as motion. Pure clock switches add no work.
    pub fn prepare_power_cycle<'model>(
        &self,
        spin_rad_s: f64,
        segments: &[RotorPeriodicSegment<'model>],
        forcing: RotorPeriodicForcing<'_>,
        relative_generator_residual_tolerance: f64,
        budget: RotorFloquetBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedRotorMaxwellPowerCycle<'model>, RotorError>
    where
        'a: 'model,
    {
        self.as_kinetics().prepare_power_cycle(
            spin_rad_s,
            segments,
            self.branches,
            forcing,
            relative_generator_residual_tolerance,
            budget,
            cancelled,
        )
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RotorTemperatureError> {
    if cancelled() {
        Err(RotorError::Cancelled.into())
    } else {
        Ok(())
    }
}
