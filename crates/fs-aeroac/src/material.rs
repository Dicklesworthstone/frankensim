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

/// Same-material SI binding for the existing three-dimensional slot-jet rig.
pub mod slot_jet_3d;
pub use slot_jet_3d::{
    MaterialSlotJet3d, MaterialSlotJet3dRun, SlotJet3dForceLine, SlotJet3dSetup,
};

use fs_lbm::core2::CollisionModel2;
use fs_lbm::material::{MaterialLbmError, MaterialScaling, plan_scaling_from_material};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::FluidTransport;
use fs_material::gas::GasState;
use fs_math::c64::C64;

use crate::curle2d::{DipoleField, dipole_pressure};
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

/// One observed bin of the actual surface-force record.
#[derive(Clone, Debug)]
pub struct ObservedJetLine {
    /// Physical bin frequency [Hz].
    pub frequency_hz: f64,
    /// Complex force per span [N/m], normalized by the Hann coherent gain.
    pub force_n_m: [C64; 2],
    /// Same-medium Curle pressure and its existing scope boundary.
    pub field: DipoleField,
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
        let (force_n_m, frequency_hz) =
            force_spectral_line(&self.force_series_n_m, self.timestep_s, bin)?;
        let field = self.radiate_force_line(force_n_m, frequency_hz, observer_m, source_m)?;
        Ok(ObservedJetLine {
            frequency_hz,
            force_n_m,
            field,
        })
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
}

impl MaterialJetLabiumRun<'_, GasState> {
    /// Original gas state for gas-bound runs, including mixture properties.
    #[must_use]
    pub const fn gas(&self) -> &GasState {
        self.material
    }
}

// Both spatial rigs use one amplitude, window and phase convention. Units
// stay those of the input history: the 2-D run supplies N/m, the 3-D run N.
fn force_spectral_line(
    forces: &[[f64; 2]],
    timestep_s: f64,
    bin: usize,
) -> Result<([C64; 2], f64), AeroacError> {
    let n = forces.len();
    if bin == 0 || bin >= n / 2 {
        return Err(AeroacError::InvalidParameter {
            what: "observed force bin must be strictly between DC and Nyquist",
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
    let mut line = [C64::new(0.0, 0.0); 2];
    for axis in 0..2 {
        let mean = forces.iter().map(|f| f[axis]).sum::<f64>() / n_f;
        for (i, value) in forces.iter().enumerate() {
            buffer[i] = fs_fft::C64::new((value[axis] - mean) * window[i], 0.0);
        }
        fft.forward(&mut buffer, &mut scratch);
        line[axis] = C64::new(buffer[bin].re, -buffer[bin].im).scale(normalization);
    }
    let frequency_hz = bin as f64 / (n_f * timestep_s);
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
