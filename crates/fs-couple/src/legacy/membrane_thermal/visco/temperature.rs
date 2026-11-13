//! Actual-temperature feedback into native finite-strain Maxwell kinetics.
//!
//! Each original thermal receiver supplies its accepted temperature. The shared
//! WLF/Arrhenius owner shifts the ORIGINAL branch times once, then the native
//! membrane advances over the thermal-selected interval. Only its physical
//! relaxation loss returns to enthalpy. Accepted endpoints and thermal events
//! synchronize the kinetic law without evolving or reinitializing any history.
//!
//! This is first-order, lagged-temperature coupling. Refine the mechanical
//! timestep ceiling together with native thermal limits for varying temperature;
//! local return-map convergence is not a bound on this splitting error. The
//! equilibrium energy, branch energies and natural reference are fixed, so the
//! kinetic shift adds no reversible storage, entropy term or reference work.
//! Expansion and temperature-dependent elastic moduli require their own
//! thermodynamic closure and cannot be introduced through this interface.

use std::borrow::Cow;
use std::sync::Arc;

use fs_exec::CancelGate;
use fs_material::visco::temperature::{HyperelasticTemperatureError, TemperatureShift};
use fs_material::visco_hyper::ViscoHyperelastic;
use fs_solid::membrane::visco::dynamics::FiniteViscoMembraneMotion;
use fs_solid::membrane::visco::owned::OwnedFiniteViscoMembraneMotion;
use fs_solid::membrane::visco::{FiniteViscoMembraneError, ViscoelasticLawField};
use fs_solid::membrane::{MembraneError, MembranePatch};
use fs_thermal::enthalpy::{DepositionObserver, DepositionSimulation};

use super::super::{
    ElementThermalCell, MechanicalInterval, MembraneThermalConfig, MembraneThermalError,
    MembraneThermalReport, Result, ThermalMotion, advance_native, check,
};
use super::{ViscoMembraneThermalControls, admit_motion, mechanical_interval};

#[derive(Debug)]
pub enum ViscoMembraneTemperatureError {
    Invalid(&'static str),
    Shift {
        element: usize,
        source: HyperelasticTemperatureError,
    },
    Native(FiniteViscoMembraneError),
}

impl core::fmt::Display for ViscoMembraneTemperatureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(message) => write!(f, "membrane thermal kinetics: {message}"),
            Self::Shift { element, source } => write!(f, "membrane triangle {element}: {source}"),
            Self::Native(source) => source.fmt(f),
        }
    }
}

impl std::error::Error for ViscoMembraneTemperatureError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Invalid(_) => None,
            Self::Shift { source, .. } => Some(source),
            Self::Native(source) => Some(source),
        }
    }
}

fn poll(
    cancelled: &mut impl FnMut() -> bool,
) -> core::result::Result<(), ViscoMembraneTemperatureError> {
    if cancelled() {
        Err(ViscoMembraneTemperatureError::Native(
            MembraneError::Cancelled.into(),
        ))
    } else {
        Ok(())
    }
}

/// Immutable original spectra and explicit kinetic calibrations, in original
/// triangle order. Calibration is caller-declared; the sheet's elastic source
/// card does not by itself establish measured Maxwell or shift parameters.
#[derive(Clone, Debug)]
pub struct ViscoMembraneKinetics {
    original: Arc<[ViscoHyperelastic]>,
    shifts: Arc<[TemperatureShift]>,
}

impl ViscoMembraneKinetics {
    /// Declare each native spectrum at the corresponding shift's reference
    /// temperature. Capture the original unshifted laws once and retain this
    /// object across continuation calls; never construct a new calibration
    /// from a previously temperature-shifted output.
    ///
    /// One shift per original triangle is required, including for a uniform
    /// original law. Geometry, accepted displacement and branch age are not
    /// changed by this preparation. No source query or material solve runs.
    pub fn try_new(
        reference_motion: &FiniteViscoMembraneMotion<'_>,
        shifts: &[TemperatureShift],
        mut cancelled: impl FnMut() -> bool,
    ) -> core::result::Result<Self, ViscoMembraneTemperatureError> {
        poll(&mut cancelled)?;
        let count = reference_motion.membrane().triangle_states().len();
        if shifts.len() != count {
            return Err(ViscoMembraneTemperatureError::Invalid(
                "one kinetic calibration per original triangle required",
            ));
        }
        let field = reference_motion.membrane().material();
        let mut original = Vec::with_capacity(count);
        for (element, shift) in shifts.iter().enumerate() {
            poll(&mut cancelled)?;
            let law = match field {
                ViscoelasticLawField::Uniform(law) => law,
                ViscoelasticLawField::PerElement(laws) => &laws[element],
            };
            // Re-admit native public parameters and the positive time domain.
            // At the declared reference temperature the shift is exactly one.
            shift
                .shift_hyperelastic(law, shift.reference_temperature_k())
                .map_err(|source| ViscoMembraneTemperatureError::Shift { element, source })?;
            original.push(law.clone());
        }
        poll(&mut cancelled)?;
        Ok(Self {
            original: original.into(),
            shifts: shifts.into(),
        })
    }

    pub fn original_laws(&self) -> &[ViscoHyperelastic] {
        &self.original
    }

    pub fn shifts(&self) -> &[TemperatureShift] {
        &self.shifts
    }

    fn at_temperatures(
        &self,
        temperatures_k: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> core::result::Result<Vec<ViscoHyperelastic>, ViscoMembraneTemperatureError> {
        poll(&mut cancelled)?;
        if temperatures_k.len() != self.original.len() {
            return Err(ViscoMembraneTemperatureError::Invalid(
                "one actual temperature per original triangle required",
            ));
        }
        let mut laws = Vec::with_capacity(self.original.len());
        for (element, ((law, shift), &temperature)) in self
            .original
            .iter()
            .zip(self.shifts.iter())
            .zip(temperatures_k)
            .enumerate()
        {
            poll(&mut cancelled)?;
            laws.push(
                shift
                    .shift_hyperelastic(law, temperature)
                    .map_err(|source| ViscoMembraneTemperatureError::Shift { element, source })?
                    .into_law(),
            );
        }
        poll(&mut cancelled)?;
        Ok(laws)
    }
}

/// Complete proposal with an owned native motion and the original thermal owner.
pub struct TemperatureViscoMembraneThermalCandidate<'thermal, O> {
    motion: OwnedFiniteViscoMembraneMotion,
    thermal: DepositionSimulation<'thermal>,
    observer: O,
    report: MembraneThermalReport,
}

impl<'thermal, O> TemperatureViscoMembraneThermalCandidate<'thermal, O> {
    pub fn motion(&self) -> &OwnedFiniteViscoMembraneMotion {
        &self.motion
    }
    pub fn thermal(&self) -> &DepositionSimulation<'thermal> {
        &self.thermal
    }
    pub fn observer(&self) -> &O {
        &self.observer
    }
    pub fn report(&self) -> &MembraneThermalReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        OwnedFiniteViscoMembraneMotion,
        DepositionSimulation<'thermal>,
        O,
        MembraneThermalReport,
    ) {
        (self.motion, self.thermal, self.observer, self.report)
    }
}

#[derive(Clone)]
struct Adapter<'call, 'calibration, 'boundary> {
    motion: Cow<'call, OwnedFiniteViscoMembraneMotion>,
    kinetics: &'calibration ViscoMembraneKinetics,
    prescribed: &'boundary [(usize, f64)],
    controls: ViscoMembraneThermalControls,
    kinetic_j: f64,
    element_heat_j: Vec<f64>,
    source_temperature_k: Vec<Option<f64>>,
    last_temperatures_k: Vec<f64>,
}

impl ThermalMotion for Adapter<'_, '_, '_> {
    fn patch(&self) -> &MembranePatch {
        self.motion.patch()
    }
    fn time_s(&self) -> f64 {
        self.motion.as_native().membrane().time_s()
    }
    fn strain_energy_j(&self) -> f64 {
        self.motion.as_native().membrane().stored_energy_j()
    }
    fn kinetic_energy_j(&self) -> f64 {
        self.kinetic_j
    }
    fn max_step_s(&self) -> f64 {
        self.controls.maximum_step_s
    }
    fn source_temperature_k(&self, element: usize) -> Option<f64> {
        self.source_temperature_k[element]
    }
    fn element_dissipated_j(&self) -> &[f64] {
        &self.element_heat_j
    }
    fn uses_temperature_feedback(&self) -> bool {
        true
    }

    fn set_temperatures<E>(
        &mut self,
        temperatures_k: &[f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(), E> {
        // Constant-temperature runs do not rebuild laws after the first bind.
        if temperatures_k == self.last_temperatures_k.as_slice() {
            return Ok(());
        }
        let laws = self
            .kinetics
            .at_temperatures(temperatures_k, &mut *cancelled)
            .map_err(MembraneThermalError::Temperature)?;
        let motion = self
            .motion
            .with_kinetics(&laws, &mut *cancelled)
            .map_err(MembraneThermalError::Viscoelastic)?;
        // This native operation changes only tau: no transverse solve, history
        // reset, work increment, heat, displacement or velocity change occurs.
        self.motion = Cow::Owned(motion);
        self.last_temperatures_k = temperatures_k.to_vec();
        Ok(())
    }

    fn advance<E>(
        &mut self,
        dt_s: f64,
        dead_load_n: &[[f64; 3]],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<MechanicalInterval, E> {
        let step = self
            .motion
            .advance_with_pressure(
                dt_s,
                dead_load_n,
                self.controls.pressure_pa,
                self.prescribed,
                self.controls.solve,
                self.controls.motion,
                cancelled,
            )
            .map_err(MembraneThermalError::Viscoelastic)?;
        let interval = mechanical_interval(
            &step.material_energy,
            &step.inertia,
            &step.equilibrium,
            step.motion.as_native().membrane().time_s(),
            dt_s,
            self.prescribed,
        )?;
        self.kinetic_j = step.inertia.kinetic_after_j;
        self.element_heat_j = step.element_dissipated_j;
        self.motion = Cow::Owned(step.motion);
        Ok(interval)
    }
}

/// Advance relaxation, native momentum and enthalpy with actual-T feedback.
///
/// The original receiver/source/inventory, support, work and cancellation
/// admission is shared with the frozen-coefficient membrane composition. Each
/// interval uses its accepted starting temperatures; its physical branch heat
/// updates the caloric cells and therefore the kinetics of subsequent intervals.
/// The returned laws also match the accepted endpoint/event temperatures.
///
/// Reuse the SAME immutable kinetics declaration on continuation. Both its
/// temperature validity and every receiver's fully-solid applicability band are
/// enforced, including after the final step. Out-of-band heat, incompatible
/// energy/reference laws, a late observer refusal or cancellation publishes no
/// mechanical, thermal or observer half-result.
///
/// Temperature changes only the declared branch times. Fixed energy parameters
/// and natural strain imply no extra reversible thermal work; backward-Euler
/// kinetic damping and endpoint integration error remain separate from heat.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate<'thermal, O: DepositionObserver>(
    motion: &OwnedFiniteViscoMembraneMotion,
    thermal: &DepositionSimulation<'thermal>,
    observer: &O,
    target_s: f64,
    dead_load_n: &[[f64; 3]],
    prescribed: &[(usize, f64)],
    cells: &[ElementThermalCell<'_>],
    kinetics: &ViscoMembraneKinetics,
    mechanics: ViscoMembraneThermalControls,
    config: MembraneThermalConfig,
    gate: &CancelGate,
) -> Result<TemperatureViscoMembraneThermalCandidate<'thermal, O>, O::Error> {
    check(gate)?;
    let native = motion.as_native();
    let (kinetic_j, source_temperature_k) = admit_motion(&native, prescribed, mechanics, gate)?;
    if kinetics.original.len() != motion.patch().sheets().len() {
        return Err(MembraneThermalError::Temperature(
            ViscoMembraneTemperatureError::Invalid(
                "kinetic field must cover every original membrane triangle",
            ),
        ));
    }
    let adapter = Adapter {
        motion: Cow::Borrowed(motion),
        kinetics,
        prescribed,
        controls: mechanics,
        kinetic_j,
        element_heat_j: vec![0.0; kinetics.original.len()],
        source_temperature_k,
        last_temperatures_k: Vec::new(),
    };
    let candidate = advance_native(
        &adapter,
        thermal,
        observer,
        target_s,
        dead_load_n,
        cells,
        config,
        gate,
    )?;
    Ok(TemperatureViscoMembraneThermalCandidate {
        motion: candidate.motion.motion.into_owned(),
        thermal: candidate.thermal,
        observer: candidate.observer,
        report: candidate.report,
    })
}
