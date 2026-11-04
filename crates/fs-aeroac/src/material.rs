//! Material-bound flow sources and acoustic observations in physical units.
//!
//! The same admitted fluid supplies viscosity to the existing incompressible
//! jet solver, density to its surface-force conversion, and sound speed to
//! Curle radiation. Changing pressure, temperature or composition recompiles
//! those quantities together. No lattice sound speed is interpreted as a
//! physical acoustic speed. Gas, liquid and real-gas sources use their original
//! transport and acoustic laws without an ideal-gas conversion. This is a
//! frozen, homogeneous, low-Mach medium; chemistry, cavitation, compressible
//! feedback and atmospheric absorption are absent.
//! The existing two-dimensional shape/scaling scope still applies to Curle
//! observation. The three-dimensional adapter retains periodic-cell force;
//! it does not reinterpret that force as an isolated acoustic source.

mod frequency;
/// Same-material SI binding for the existing three-dimensional slot-jet rig.
pub mod slot_jet_3d;
pub use frequency::{ForceFrequencyRefinement, ForceToneFrequencyRefinement};
pub use slot_jet_3d::{
    MaterialSlotJet3d, MaterialSlotJet3dRun, SlotJet3dForceHarmonics, SlotJet3dForceLine,
    SlotJet3dForceSpectrum, SlotJet3dForceTones, SlotJet3dSetup,
};

use fs_fft::harmonic::{HarmonicFitError, HarmonicFitPlan, HarmonicFitSummary, ToneFitPlan};
use fs_lbm::core2::CollisionModel2;
use fs_lbm::material::{MaterialLbmError, MaterialScaling, plan_scaling_from_material};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::FluidTransport;
use fs_material::gas::GasState;
use fs_math::c64::C64;

use crate::curle2d::{DipoleField, PreparedDipoleObserver, dipole_pressure};
use crate::jetlab::{JetLabiumConfig, JetLabiumDiagnostics, run_jet_labium, validate_config};
use crate::{AeroacError, SCOPE_STATEMENT};

/// Physical applicability limit for the incompressible source approximation.
/// This is separate from the lattice solver's numerical Mach bound.
pub const MAX_PHYSICAL_MACH: f64 = 0.3;

/// Geometry and numerical run controls for the existing two-dimensional rig.
/// Lengths are in cells; [`JetPhysicalState::cell_size_m`] gives their SI size.
/// Flow velocity and relaxation are derived from the gas and physical inputs.
#[derive(Clone, Debug)]
pub struct JetLabiumSetup {
    /// Periodic domain length [cells].
    pub nx: usize,
    /// Domain height [cells].
    pub ny: usize,
    /// Slot half height [cells].
    pub slot_half: f64,
    /// Profile edge smoothing width [cells].
    pub slot_smoothing: f64,
    /// Distance to the plate's leading edge [cells].
    pub edge_distance: usize,
    /// Plate length [cells].
    pub plate_length: usize,
    /// Right-side fringe width [cells].
    pub fringe_width: usize,
    /// Numerical fringe strength.
    pub fringe_sigma: f64,
    /// Settling steps before the force record.
    pub steps_settle: usize,
    /// Recorded steps; power of two, at least 64.
    pub steps_record: usize,
    /// Transverse perturbation relative to jet speed.
    pub seed_amplitude: f64,
    /// Nozzle-wall thickness [cells]; zero disables the nozzle.
    pub nozzle_thickness: usize,
    /// Existing collision operator.
    pub collision: CollisionModel2,
}

/// Physical scenario and discretization, without independent gas coefficients.
#[derive(Clone, Copy, Debug)]
pub struct JetPhysicalState {
    /// Uniform cell size [m].
    pub cell_size_m: f64,
    /// Integration step [s].
    pub timestep_s: f64,
    /// Prescribed jet speed [m/s].
    pub jet_speed_m_s: f64,
    /// Declared speed envelope [m/s], including perturbations and transients.
    pub max_speed_m_s: f64,
    /// Maximum admitted grid cells before allocation.
    pub max_cells: usize,
}

/// Material/units admission and actual lattice execution share one error type.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialAeroacError {
    /// The common material-to-lattice owner refused its input.
    Scaling(MaterialLbmError),
    /// The original material acoustic owner refused its state.
    Material(MaterialError),
    /// Geometry, applicability, radiation or solver refusal.
    Aeroac(AeroacError),
}

impl core::fmt::Display for MaterialAeroacError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Scaling(error) => write!(f, "{error}"),
            Self::Material(error) => write!(f, "{error}"),
            Self::Aeroac(error) => write!(f, "{error}"),
        }
    }
}
impl std::error::Error for MaterialAeroacError {}
impl From<MaterialLbmError> for MaterialAeroacError {
    fn from(error: MaterialLbmError) -> Self {
        Self::Scaling(error)
    }
}
impl From<AeroacError> for MaterialAeroacError {
    fn from(error: AeroacError) -> Self {
        Self::Aeroac(error)
    }
}

fn invalid(what: &'static str) -> MaterialAeroacError {
    AeroacError::InvalidParameter { what }.into()
}

// Both rigs retain one frozen acoustic snapshot beside the actual transport
// binding. A caller implementing the two views cannot mix different densities.
fn acoustic_from_scaling<F: FluidTransport + FluidAcoustics + ?Sized>(
    scaling: &MaterialScaling<'_, F>,
) -> Result<FluidAcousticState, MaterialAeroacError> {
    let acoustic = scaling
        .material()
        .acoustic_state()
        .map_err(MaterialAeroacError::Material)?;
    if acoustic.density_kg_m3().to_bits() != scaling.transport().density_kg_m3().to_bits() {
        return Err(invalid(
            "flow and acoustic material views must share the same density",
        ));
    }
    Ok(acoustic)
}

/// One compiled medium/geometry binding. Derived controls are read-only.
pub struct MaterialJetLabium<'a, F: FluidTransport + FluidAcoustics + ?Sized = GasState> {
    scaling: MaterialScaling<'a, F>,
    acoustic: FluidAcousticState,
    config: JetLabiumConfig,
    max_speed_m_s: f64,
    force_per_span_unit_n_m: f64,
}

impl<'a, F: FluidTransport + FluidAcoustics + ?Sized> MaterialJetLabium<'a, F> {
    /// Bind one original gas, liquid or real-gas material state to the real rig.
    /// Transport and acoustic views must have identical density. Source receipts
    /// and phase policy remain on the borrowed material; no gas surrogate is made.
    ///
    /// # Errors
    /// Refuses inconsistent gas fields, invalid geometry/budgets, unrepresentable
    /// units, an insufficient speed envelope or excessive physical/lattice Mach.
    pub fn try_new(
        material: &'a F,
        setup: &JetLabiumSetup,
        physical: JetPhysicalState,
    ) -> Result<Self, MaterialAeroacError> {
        if setup
            .nx
            .checked_mul(setup.ny)
            .is_none_or(|n| n == 0 || n > physical.max_cells)
        {
            return Err(invalid("jet grid exceeds the declared cell budget"));
        }
        if !physical.jet_speed_m_s.is_finite() || physical.jet_speed_m_s <= 0.0 {
            return Err(invalid("physical jet speed must be finite and positive"));
        }
        let scaling = plan_scaling_from_material(
            material,
            physical.cell_size_m,
            physical.timestep_s,
            physical.max_speed_m_s,
        )?;
        let acoustic = acoustic_from_scaling(&scaling)?;
        if physical.max_speed_m_s / acoustic.sound_speed_m_s() >= MAX_PHYSICAL_MACH {
            return Err(invalid(
                "fluid speed envelope exceeds the incompressible physical Mach limit",
            ));
        }
        let config = JetLabiumConfig {
            nx: setup.nx,
            ny: setup.ny,
            slot_half: setup.slot_half,
            slot_smoothing: setup.slot_smoothing,
            u_jet: scaling.velocity_lattice(physical.jet_speed_m_s)?,
            tau: scaling.flow().tau,
            edge_distance: setup.edge_distance,
            plate_length: setup.plate_length,
            fringe_width: setup.fringe_width,
            fringe_sigma: setup.fringe_sigma,
            steps_settle: setup.steps_settle,
            steps_record: setup.steps_record,
            seed_amplitude: setup.seed_amplitude,
            nozzle_thickness: setup.nozzle_thickness,
            collision: setup.collision,
        };
        validate_config(&config)?;
        if !config.steps_record.is_power_of_two() || config.steps_record < 64 {
            return Err(invalid("steps_record must be a power of two >= 64"));
        }
        let initial_speed = physical.jet_speed_m_s
            * fs_math::det::sqrt(1.0 + setup.seed_amplitude * setup.seed_amplitude);
        if initial_speed > physical.max_speed_m_s {
            return Err(invalid(
                "seeded jet exceeds the declared physical speed envelope",
            ));
        }
        // A 2-D cell has mass rho*dx^2 per unit span. Momentum exchange
        // per step therefore scales as rho*dx^3/dt^2, i.e. N/m, not N.
        let velocity_unit = scaling.velocity_si(1.0)?;
        let force_unit = scaling.transport().density_kg_m3()
            * velocity_unit
            * velocity_unit
            * physical.cell_size_m;
        if !force_unit.is_finite() || force_unit <= 0.0 {
            return Err(invalid("force per unit span conversion is unrepresentable"));
        }
        Ok(Self {
            scaling,
            acoustic,
            config,
            max_speed_m_s: physical.max_speed_m_s,
            force_per_span_unit_n_m: force_unit,
        })
    }

    /// Shared material/transport/units binding used by the flow solver.
    #[must_use]
    pub const fn scaling(&self) -> &MaterialScaling<'a, F> {
        &self.scaling
    }

    /// Frozen acoustic coefficients shared by Mach admission and observations.
    #[must_use]
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Compiled lattice inputs; inspectable without independently mutating fluid
    /// viscosity, jet speed, density or acoustic speed.
    #[must_use]
    pub const fn lattice_config(&self) -> &JetLabiumConfig {
        &self.config
    }

    /// Execute the existing jet-labium solver and return its physical force.
    /// The diagnostic bound is sampled on two planes every 16 steps; passing
    /// it is not a proof of a global speed bound or grid convergence.
    ///
    /// # Errors
    /// Existing rig refusals, an observed speed-envelope violation, or
    /// unrepresentable physical force refuse before any result is returned.
    pub fn run(&self) -> Result<MaterialJetLabiumRun<'_, F>, MaterialAeroacError> {
        let mut run = run_jet_labium(&self.config)?;
        let observed_speed = self.scaling.velocity_si(run.diagnostics.mach_max_lattice)?;
        if observed_speed > self.max_speed_m_s {
            return Err(invalid(
                "sampled jet speed exceeds the declared physical envelope",
            ));
        }
        for force in &mut run.force_series {
            for component in force {
                *component *= self.force_per_span_unit_n_m;
                if !component.is_finite() {
                    return Err(invalid("physical force per unit span is unrepresentable"));
                }
            }
        }
        Ok(MaterialJetLabiumRun {
            material: self.scaling.material(),
            acoustic: self.acoustic,
            timestep_s: self.scaling.timestep_s(),
            force_series_n_m: run.force_series,
            sampled_physical_mach_max: observed_speed / self.acoustic.sound_speed_m_s(),
            lattice_diagnostics: run.diagnostics,
        })
    }
}

/// The actual force history and observations in the same frozen fluid state.
pub struct MaterialJetLabiumRun<'a, F: FluidTransport + FluidAcoustics + ?Sized = GasState> {
    material: &'a F,
    acoustic: FluidAcousticState,
    timestep_s: f64,
    force_series_n_m: Vec<[f64; 2]>,
    /// Maximum sampled physical speed divided by this fluid's acoustic speed.
    /// Samples are restricted to the existing rig's diagnostic planes/ticks.
    pub sampled_physical_mach_max: f64,
    /// Existing Reynolds and numerical flow diagnostics, with their original
    /// lattice-unit conventions retained.
    pub lattice_diagnostics: JetLabiumDiagnostics,
}

/// One observed frequency of the actual surface-force record.
#[derive(Clone, Debug)]
pub struct ObservedJetLine {
    /// Physical frequency [Hz].
    pub frequency_hz: f64,
    /// Complex peak force per span [N/m], in exp(-i omega t) convention.
    /// FFT observations use Hann coherent gain; fitted tones use the
    /// joint least-squares amplitudes without a window-gain approximation.
    pub force_n_m: [C64; 2],
    /// Same-medium Curle pressure and its existing scope boundary.
    pub field: DipoleField,
}

/// Prepared Hann-windowed spectrum of the actual 2-D source record.
///
/// Both force components are transformed once. Subsequent bin/receiver queries
/// allocate no FFT workspace and preserve the original material sound speed.
/// Finite-window leakage remains; this does not fit off-bin tones.
pub struct JetForceSpectrum {
    acoustic: FluidAcousticState,
    spectrum: ForceSpectrum,
}

/// Joint Welch spectrum of the actual physical x/y force record. Retains the
/// original acoustic material and measured force-component interference.
pub struct JetForceCrossSpectrum {
    acoustic: FluidAcousticState,
    spectrum: fs_fft::welch::WelchCrossSpectrum,
}

impl JetForceCrossSpectrum {
    /// Recorded window, sample clock, coverage and /Hz cross-spectrum factors.
    pub fn spectrum(&self) -> &fs_fft::welch::WelchCrossSpectrum {
        &self.spectrum
    }

    /// Same frozen fluid used to admit and execute the original flow.
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Observe the power in one resolved interior bin, with its original
    /// source-force units and shared material. This does not fit an isolated
    /// tone or promote the 2D flow source into an absolute-SPL authority.
    pub fn observe_bin(
        &self,
        bin: usize,
        observer_m: [f64; 2],
        source_m: [f64; 2],
        budget: crate::curle2d::CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<crate::curle2d::DipoleFieldStatistics, AeroacError> {
        if cancelled() {
            return Err(AeroacError::Cancelled {
                what: "recorded flow-source observation",
            });
        }
        let frequency =
            self.spectrum
                .frequency_hz(bin)
                .map_err(|_| AeroacError::InvalidParameter {
                    what: "flow-source Welch bin",
                })?;
        let observer =
            PreparedDipoleObserver::new(self.acoustic, frequency, observer_m, &[source_m], 1)?;
        observer.evaluate_welch_bin(&self.spectrum, bin, budget, cancelled)
    }
}

impl JetForceSpectrum {
    /// Original record length; interior positive bins are 1..sample_count/2.
    #[must_use]
    pub const fn sample_count(&self) -> usize {
        self.spectrum.sample_count
    }

    /// Original sample interval [s].
    #[must_use]
    pub const fn timestep_s(&self) -> f64 {
        self.spectrum.timestep_s
    }

    /// Same frozen fluid snapshot used by the actual flow.
    #[must_use]
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Observe one cached force bin at an actual source/receiver geometry.
    /// DC, Nyquist, invalid bins, geometry and nonfinite fields refuse.
    pub fn observe_spectral_line(
        &self,
        bin: usize,
        observer_m: [f64; 2],
        source_m: [f64; 2],
    ) -> Result<ObservedJetLine, AeroacError> {
        let (force_n_m, frequency_hz) = self.spectrum.line(bin)?;
        observe_force(self.acoustic, force_n_m, frequency_hz, observer_m, source_m)
    }
}

/// Joint DC and specified-harmonic fit of one actual 2-D force record.
///
/// No source harmonic is invented: the common L1 fitter separates the declared
/// frequencies in the recorded samples, including a non-bin fundamental.
/// Residual RMS measures unexplained samples, not uncertainty or omitted physics.
pub struct JetForceHarmonics {
    acoustic: FluidAcousticState,
    harmonics: Vec<usize>,
    fundamental_hz: f64,
    force_n_m: Vec<[C64; 2]>,
    /// Simultaneously fitted constant `(Fx, Fy)` forces [N/m].
    pub dc_force_n_m: [f64; 2],
    /// Unweighted sample residual RMS for each force component [N/m].
    pub residual_rms_n_m: [f64; 2],
}

impl JetForceHarmonics {
    /// Requested harmonic orders, in the shared plan's original order.
    #[must_use]
    pub fn harmonics(&self) -> &[usize] {
        &self.harmonics
    }

    /// Fundamental used by the joint fit [Hz], supplied or locally refined.
    #[must_use]
    pub const fn fundamental_hz(&self) -> f64 {
        self.fundamental_hz
    }

    /// Complex peak `(Fx, Fy)` amplitudes [N/m], aligned with `harmonics()`.
    /// Phase is referenced to the first recorded sample, in exp(-i omega t).
    #[must_use]
    pub fn force_n_m(&self) -> &[[C64; 2]] {
        &self.force_n_m
    }

    /// Exact source medium retained for all harmonic radiation paths.
    #[must_use]
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Radiate a fitted harmonic using its actual frequency and original fluid.
    /// An undeclared harmonic order refuses rather than interpolating a tone.
    pub fn observe_harmonic(
        &self,
        harmonic: usize,
        observer_m: [f64; 2],
        source_m: [f64; 2],
    ) -> Result<ObservedJetLine, AeroacError> {
        let index =
            self.harmonics
                .binary_search(&harmonic)
                .map_err(|_| AeroacError::InvalidParameter {
                    what: "requested force harmonic was not fitted",
                })?;
        observe_force(
            self.acoustic,
            self.force_n_m[index],
            self.fundamental_hz * harmonic as f64,
            observer_m,
            source_m,
        )
    }
}

/// Joint DC and arbitrary-frequency fit of the actual 2-D surface forces.
///
/// The declared frequencies may be noncommensurate, such as independent
/// shedding tones. They are measured from the same force record, with no
/// integer-harmonic rounding, source synthesis or invented phase locking.
pub struct JetForceTones {
    acoustic: FluidAcousticState,
    frequencies_hz: Vec<f64>,
    force_n_m: Vec<[C64; 2]>,
    /// Simultaneously fitted constant `(Fx, Fy)` forces [N/m].
    pub dc_force_n_m: [f64; 2],
    /// Unweighted sample residual RMS for each force component [N/m].
    /// Unmodeled tones or frequency errors can also bias fitted amplitudes.
    pub residual_rms_n_m: [f64; 2],
}

impl JetForceTones {
    /// Positive fitted frequencies [Hz], in the admitted increasing order.
    #[must_use]
    pub fn frequencies_hz(&self) -> &[f64] {
        &self.frequencies_hz
    }

    /// Complex peak `(Fx, Fy)` amplitudes [N/m], aligned with `frequencies_hz()`.
    /// Phase is referenced to the first recorded sample, in exp(-i omega t).
    #[must_use]
    pub fn force_n_m(&self) -> &[[C64; 2]] {
        &self.force_n_m
    }

    /// Exact source medium retained for every physical radiation path.
    #[must_use]
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Radiate a fitted tone at its true frequency in the original fluid.
    /// `index` is zero-based in `frequencies_hz()`; an absent index refuses.
    /// No FFT-bin lookup, spectral interpolation or harmonic-order inference.
    pub fn observe_tone(
        &self,
        index: usize,
        observer_m: [f64; 2],
        source_m: [f64; 2],
    ) -> Result<ObservedJetLine, AeroacError> {
        let frequency_hz =
            *self
                .frequencies_hz
                .get(index)
                .ok_or(AeroacError::InvalidParameter {
                    what: "requested force tone was not fitted",
                })?;
        observe_force(
            self.acoustic,
            self.force_n_m[index],
            frequency_hz,
            observer_m,
            source_m,
        )
    }
}

fn observe_force(
    acoustic: FluidAcousticState,
    force_n_m: [C64; 2],
    frequency_hz: f64,
    observer_m: [f64; 2],
    source_m: [f64; 2],
) -> Result<ObservedJetLine, AeroacError> {
    let wavenumber = core::f64::consts::TAU * frequency_hz / acoustic.sound_speed_m_s();
    let field = dipole_pressure(force_n_m, wavenumber, observer_m, source_m)?;
    Ok(ObservedJetLine {
        frequency_hz,
        force_n_m,
        field,
    })
}

impl<F: FluidTransport + FluidAcoustics + ?Sized> MaterialJetLabiumRun<'_, F> {
    /// Surface-force history `(Fx, Fy)` per unit span [N/m], sampled every dt.
    #[must_use]
    pub fn force_series_n_m(&self) -> &[[f64; 2]] {
        &self.force_series_n_m
    }

    /// Physical sample period [s].
    #[must_use]
    pub const fn timestep_s(&self) -> f64 {
        self.timestep_s
    }

    /// Original material shared by flow and observation, retaining source data.
    #[must_use]
    pub const fn material(&self) -> &F {
        self.material
    }

    /// Exact frozen acoustic state used by the flow's physical Mach admission.
    #[must_use]
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Estimate the joint force cross spectrum from the accepted physical
    /// history, without copying/transposing it or inventing x/y independence.
    /// The original material and sample clock follow every bin observation.
    pub fn prepare_force_cross_spectrum(
        &self,
        config: fs_fft::welch::WelchConfig,
        budget: fs_fft::welch::WelchBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<JetForceCrossSpectrum, fs_fft::welch::WelchError> {
        Ok(JetForceCrossSpectrum {
            acoustic: self.acoustic,
            spectrum: fs_fft::welch::WelchCrossSpectrum::estimate_interleaved(
                &self.force_series_n_m,
                self.timestep_s,
                config,
                budget,
                cancelled,
            )?,
        })
    }

    /// Authority boundary retained from the existing two-dimensional rig.
    #[must_use]
    pub const fn scope(&self) -> &'static str {
        SCOPE_STATEMENT
    }

    /// Observe a bin directly from the accepted physical force history.
    ///
    /// Each component uses the existing jet spectral instrument's mean removal
    /// and symmetric Hann window. The shared FFT is normalized by twice the
    /// inverse window sum, so a resolved sinusoid retains its force amplitude.
    /// Its coefficient is conjugated into `exp(-i omega t)` convention, with
    /// phase measured from the first recorded force sample. Finite-window
    /// leakage remains; this does not certify a stationary or isolated tone.
    ///
    /// # Errors
    /// DC, Nyquist and out-of-range bins refuse, as do invalid observations.
    pub fn observe_spectral_line(
        &self,
        bin: usize,
        observer_m: [f64; 2],
        source_m: [f64; 2],
    ) -> Result<ObservedJetLine, AeroacError> {
        check_force_bin(self.force_series_n_m.len(), bin)?;
        self.prepare_force_spectrum(self.force_series_n_m.len())?
            .observe_spectral_line(bin, observer_m, source_m)
    }

    /// Transform both recorded forces once for repeated harmonic/bin queries.
    /// The caller's sample-count bound is checked before allocating O(N)
    /// spectrum and FFT storage. Preparation is O(N log N); each lookup is O(1)
    /// plus the existing acoustic transfer. No physical source law changes.
    pub fn prepare_force_spectrum(
        &self,
        max_samples: usize,
    ) -> Result<JetForceSpectrum, AeroacError> {
        Ok(JetForceSpectrum {
            acoustic: self.acoustic,
            spectrum: ForceSpectrum::new(&self.force_series_n_m, self.timestep_s, max_samples)?,
        })
    }

    /// Fit DC and all declared harmonics jointly on the original sample clock.
    /// Reuse one L1 QR plan for both axes; frequencies need not be FFT bins.
    /// An incorrect record length or sample period refuses. The plan controls
    /// memory and conditioning; one additional N-sample axis buffer is reused.
    /// This estimates resolved periodic content, not stationarity or source
    /// compactness, and does not recover harmonics already aliased by sampling.
    pub fn fit_force_harmonics(
        &self,
        plan: &HarmonicFitPlan,
    ) -> Result<JetForceHarmonics, AeroacError> {
        let (force_n_m, summary) =
            fit_force_harmonics(&self.force_series_n_m, self.timestep_s, plan)?;
        Ok(JetForceHarmonics {
            acoustic: self.acoustic,
            harmonics: plan.harmonics().to_vec(),
            fundamental_hz: plan.fundamental_hz(),
            force_n_m,
            dc_force_n_m: summary.map(|s| s.dc),
            residual_rms_n_m: summary.map(|s| s.residual_rms),
        })
    }

    /// Refine one common fundamental from both recorded physical force axes.
    /// The original clock, relative component strengths, DC, phase origin and
    /// N/m units are retained. Every candidate shares one QR across both axes;
    /// no preferred coordinate direction or separate pitch average is chosen.
    ///
    /// The seed and local interval come from the caller. This does not discover
    /// a global shedding frequency, remove aliasing, or validate compactness.
    /// The total memory budget includes transposition of the accepted record.
    pub fn fit_force_harmonics_refined(
        &self,
        initial_fundamental_hz: f64,
        harmonics: &[usize],
        budget: fs_fft::harmonic::HarmonicFrequencyBudget,
    ) -> Result<
        (JetForceHarmonics, ForceFrequencyRefinement),
        fs_fft::harmonic::HarmonicFrequencyError,
    > {
        let (force_n_m, summary, refinement) = frequency::refine_force_harmonics(
            &self.force_series_n_m,
            self.timestep_s,
            initial_fundamental_hz,
            harmonics,
            budget,
        )?;
        Ok((
            JetForceHarmonics {
                acoustic: self.acoustic,
                harmonics: harmonics.to_vec(),
                fundamental_hz: refinement.fundamental_hz,
                force_n_m,
                dc_force_n_m: summary.map(|s| s.dc),
                residual_rms_n_m: summary.map(|s| s.residual_rms),
            },
            refinement,
        ))
    }

    /// Fit the actual force record at declared, possibly noncommensurate tones.
    /// The common L1 plan jointly resolves DC and every tone; both force axes
    /// reuse its factorization and scratch on the original count and clock.
    /// Known frequencies are inputs, not estimated pitches. This does not
    /// remove aliasing, or establish stationarity or source compactness.
    pub fn fit_force_tones(&self, plan: &ToneFitPlan) -> Result<JetForceTones, AeroacError> {
        let (force_n_m, summary) = fit_force_tones(&self.force_series_n_m, self.timestep_s, plan)?;
        Ok(JetForceTones {
            acoustic: self.acoustic,
            frequencies_hz: plan.frequencies_hz().to_vec(),
            force_n_m,
            dc_force_n_m: summary.map(|s| s.dc),
            residual_rms_n_m: summary.map(|s| s.residual_rms),
        })
    }

    /// Jointly refine independent, possibly inharmonic force frequencies.
    /// Every candidate reuses one L1 fit across both recorded axes, retaining
    /// their physical relative amplitudes, independent DC and sample-zero phase.
    /// The caller supplies a local seed and nonoverlapping interval per tone;
    /// resolved lines continue through same-fluid Curle radiation in N/m.
    ///
    /// This is local stationary estimation, not unknown-source discovery,
    /// alias removal or a statistical frequency certificate. Ill-conditioned
    /// frequency separation refuses, including a requested but absent tone.
    /// The total memory budget includes transposition of the accepted record.
    pub fn fit_force_tones_refined(
        &self,
        initial_frequencies_hz: &[f64],
        frequency_bounds_hz: &[[f64; 2]],
        budget: fs_fft::harmonic::ToneFrequencyBudget,
    ) -> Result<(JetForceTones, ForceToneFrequencyRefinement), fs_fft::harmonic::ToneFrequencyError>
    {
        let fit = frequency::refine_force_tones(
            &self.force_series_n_m,
            self.timestep_s,
            initial_frequencies_hz,
            frequency_bounds_hz,
            budget,
        )?;
        Ok((
            JetForceTones {
                acoustic: self.acoustic,
                frequencies_hz: fit.frequencies_hz,
                force_n_m: fit.force,
                dc_force_n_m: fit.summary.map(|s| s.dc),
                residual_rms_n_m: fit.summary.map(|s| s.residual_rms),
            },
            fit.refinement,
        ))
    }

    /// Radiate a caller-windowed force line [N/m] in the workspace's
    /// `exp(-i omega t)` convention. Frequency is in Hz and coordinates in m;
    /// wavenumber is always derived from the very fluid that generated the flow.
    /// No source-amplitude calibration or absorption model is introduced.
    ///
    /// # Errors
    /// Refuses nonpositive/nonfinite frequency, frequency at or above this
    /// record's Nyquist limit, and any existing Curle-kernel refusal.
    pub fn radiate_force_line(
        &self,
        force_n_m: [C64; 2],
        frequency_hz: f64,
        observer_m: [f64; 2],
        source_m: [f64; 2],
    ) -> Result<DipoleField, AeroacError> {
        if !frequency_hz.is_finite() || frequency_hz <= 0.0 || frequency_hz * self.timestep_s >= 0.5
        {
            return Err(AeroacError::InvalidParameter {
                what: "frequency must be positive, finite and below the force record's Nyquist limit",
            });
        }
        let wavenumber =
            (2.0 * core::f64::consts::PI * frequency_hz) / self.acoustic.sound_speed_m_s();
        dipole_pressure(force_n_m, wavenumber, observer_m, source_m)
    }

    /// Prepare interference-aware observation of several compact source forces
    /// in this run's exact fluid, frequency band and physical coordinate frame.
    /// Reuse the returned operator for many same-frequency source excitations.
    /// Unlike collapsing their net force, this retains each propagation path
    /// and the near-field particle velocity. Force phasors must share one time
    /// origin; independently started flow records do not imply phase locking.
    /// The returned operator also accepts explicitly independent sources or
    /// caller-supplied coherence factors, retaining the same material snapshot.
    ///
    /// # Errors
    /// Refuses frequencies outside the recorded band, invalid geometry, or a
    /// source count above `max_sources`; see [`PreparedDipoleObserver::new`].
    pub fn prepare_coherent_observer(
        &self,
        frequency_hz: f64,
        observer_m: [f64; 2],
        source_positions_m: &[[f64; 2]],
        max_sources: usize,
    ) -> Result<PreparedDipoleObserver, AeroacError> {
        if !frequency_hz.is_finite() || frequency_hz <= 0.0 || frequency_hz * self.timestep_s >= 0.5
        {
            return Err(AeroacError::InvalidParameter {
                what: "frequency must be positive, finite and below the force record's Nyquist limit",
            });
        }
        PreparedDipoleObserver::new(
            self.acoustic,
            frequency_hz,
            observer_m,
            source_positions_m,
            max_sources,
        )
    }
}

impl MaterialJetLabiumRun<'_, GasState> {
    /// Original gas state for gas-bound runs, including mixture properties.
    #[must_use]
    pub const fn gas(&self) -> &GasState {
        self.material
    }
}

fn check_force_bin(n: usize, bin: usize) -> Result<(), AeroacError> {
    if bin == 0 || bin >= n / 2 {
        return Err(AeroacError::InvalidParameter {
            what: "observed force bin must be strictly between DC and Nyquist",
        });
    }
    Ok(())
}

// Both spatial rigs use one amplitude, window and phase convention. Units
// stay those of the input history: the 2-D run supplies N/m, the 3-D run N.
struct ForceSpectrum {
    sample_count: usize,
    timestep_s: f64,
    lines: Vec<[C64; 2]>,
}

impl ForceSpectrum {
    fn new(forces: &[[f64; 2]], timestep_s: f64, max_samples: usize) -> Result<Self, AeroacError> {
        let n = forces.len();
        if n < 4 || !n.is_power_of_two() || n > max_samples {
            return Err(AeroacError::InvalidParameter {
                what: "force spectrum requires a power-of-two record within the sample budget",
            });
        }
        if !timestep_s.is_finite() || timestep_s <= 0.0 {
            return Err(AeroacError::InvalidParameter {
                what: "force spectrum requires a positive finite sample interval",
            });
        }
        let n_f = n as f64;
        let window: Vec<f64> = (0..n)
            .map(|i| 0.5 - 0.5 * fs_math::det::cos(core::f64::consts::TAU * i as f64 / (n_f - 1.0)))
            .collect();
        let normalization = 2.0 / window.iter().sum::<f64>();
        let fft = fs_fft::Fft::new(n);
        let mut buffer = vec![fs_fft::C64::default(); n];
        let mut scratch = vec![fs_fft::C64::default(); n];
        let mut lines = vec![[C64::ZERO; 2]; n / 2];
        for axis in 0..2 {
            let mean = forces.iter().map(|f| f[axis]).sum::<f64>() / n_f;
            for (i, value) in forces.iter().enumerate() {
                buffer[i] = fs_fft::C64::new((value[axis] - mean) * window[i], 0.0);
            }
            fft.forward(&mut buffer, &mut scratch);
            for (bin, line) in lines.iter_mut().enumerate().skip(1) {
                line[axis] = C64::new(buffer[bin].re, -buffer[bin].im).scale(normalization);
            }
        }
        Ok(Self {
            sample_count: n,
            timestep_s,
            lines,
        })
    }

    fn line(&self, bin: usize) -> Result<([C64; 2], f64), AeroacError> {
        check_force_bin(self.sample_count, bin)?;
        let line = self.lines[bin];
        let frequency_hz = bin as f64 / (self.sample_count as f64 * self.timestep_s);
        if !frequency_hz.is_finite()
            || frequency_hz <= 0.0
            || line
                .iter()
                .any(|value| !value.re.is_finite() || !value.im.is_finite())
        {
            return Err(AeroacError::NonFinite {
                what: "physical force spectrum or frequency is unrepresentable",
            });
        }
        Ok((line, frequency_hz))
    }
}

fn fit_force_harmonics(
    forces: &[[f64; 2]],
    timestep_s: f64,
    plan: &HarmonicFitPlan,
) -> Result<(Vec<[C64; 2]>, [HarmonicFitSummary; 2]), AeroacError> {
    fit_force_components(
        forces,
        timestep_s,
        plan.sample_count(),
        plan.timestep_s(),
        plan.harmonics().len(),
        |samples, scratch, phasors| plan.fit_into(samples, scratch, phasors),
    )
}

fn fit_force_tones(
    forces: &[[f64; 2]],
    timestep_s: f64,
    plan: &ToneFitPlan,
) -> Result<(Vec<[C64; 2]>, [HarmonicFitSummary; 2]), AeroacError> {
    fit_force_components(
        forces,
        timestep_s,
        plan.sample_count(),
        plan.timestep_s(),
        plan.frequencies_hz().len(),
        |samples, scratch, phasors| plan.fit_into(samples, scratch, phasors),
    )
}

fn fit_force_components(
    forces: &[[f64; 2]],
    timestep_s: f64,
    sample_count: usize,
    plan_timestep_s: f64,
    tone_count: usize,
    fit: impl Fn(&[f64], &mut [f64], &mut [fs_fft::C64]) -> Result<HarmonicFitSummary, HarmonicFitError>,
) -> Result<(Vec<[C64; 2]>, [HarmonicFitSummary; 2]), AeroacError> {
    if sample_count != forces.len() || plan_timestep_s.to_bits() != timestep_s.to_bits() {
        return Err(AeroacError::InvalidParameter {
            what: "tone or harmonic fit plan must match the original force record and sample clock",
        });
    }
    let mut axis_samples = vec![0.0; forces.len()];
    let mut scratch = vec![0.0; forces.len()];
    let mut phasors = vec![fs_fft::C64::default(); tone_count];
    let mut forces_fitted = vec![[C64::ZERO; 2]; tone_count];
    let fit_axis = |axis: usize,
                    samples: &mut [f64],
                    scratch: &mut [f64],
                    phasors: &mut [fs_fft::C64],
                    fitted: &mut [[C64; 2]]|
     -> Result<HarmonicFitSummary, AeroacError> {
        for (sample, force) in samples.iter_mut().zip(forces) {
            *sample = force[axis];
        }
        let summary = fit(samples, scratch, phasors).map_err(AeroacError::HarmonicFit)?;
        for (force, phasor) in fitted.iter_mut().zip(phasors) {
            force[axis] = C64::new(phasor.re, phasor.im);
        }
        Ok(summary)
    };
    let x = fit_axis(
        0,
        &mut axis_samples,
        &mut scratch,
        &mut phasors,
        &mut forces_fitted,
    )?;
    let y = fit_axis(
        1,
        &mut axis_samples,
        &mut scratch,
        &mut phasors,
        &mut forces_fitted,
    )?;
    Ok((forces_fitted, [x, y]))
}

#[cfg(test)]
#[path = "material/spectrum_tests.rs"]
mod spectrum_tests;
