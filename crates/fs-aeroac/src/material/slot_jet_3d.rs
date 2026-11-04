//! One physical fluid state drives the original three-dimensional slot-jet rig.
//!
//! The central-moment viscosity, SI force and acoustic wavenumber share the
//! same gas, liquid or real-gas source. The aggregate record integrates one
//! periodic span in N; [`distributed`] retains spatial force patches for the
//! existing three-dimensional acoustic observer. These are finite prescribed
//! force distributions, without an isolated-radiator qualification.

use super::{
    ForceSpectrum, JetPhysicalState, MAX_PHYSICAL_MACH, MaterialAeroacError, acoustic_from_scaling,
    check_force_bin, fit_force_harmonics, fit_force_tones, invalid,
};
use crate::slot_jet_3d::{
    SlotJet3dConfig, SlotJet3dDiagnostics, SlotJet3dRun, SlotJet3dSgsDiagnostics,
    SlotJet3dSgsError, run_slot_jet_3d, run_slot_jet_3d_with_sgs, validate_config,
};
use crate::{AeroacError, SCOPE_STATEMENT};
use fs_exec::CancelGate;
use fs_fft::harmonic::{HarmonicFitPlan, ToneFitPlan};
use fs_lbm::d3q19::CollisionModel3;
use fs_lbm::material::{MaterialScaling, plan_scaling_from_material};
use fs_lbm::sgs3::SgsConfig;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::FluidTransport;
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};

/// Spatial plate forces and coherent 3-D observation in this material's units.
pub mod distributed;

/// Geometry in cells and numerical controls, with no independent viscosity.
#[derive(Clone, Debug)]
pub struct SlotJet3dSetup {
    /// Periodic domain length [cells].
    pub nx: usize,
    /// Cross-stream domain height [cells].
    pub ny: usize,
    /// Periodic span [cells].
    pub nz: usize,
    /// Binary slit half-width [cells].
    pub slot_half: f64,
    /// Central-moment higher-order relaxation, independent of physical viscosity.
    pub higher_order_rate: f64,
    /// Nozzle-wall thickness [cells].
    pub nozzle_thickness: usize,
    /// Nozzle exit to splitter leading edge [cells].
    pub edge_distance: usize,
    /// Splitter length [cells].
    pub plate_length: usize,
    /// Right-side fringe width [cells].
    pub fringe_width: usize,
    /// Numerical fringe strength.
    pub fringe_sigma: f64,
    /// Transverse seed relative to jet speed, before the rig's spanwise factor.
    pub seed_amplitude: f64,
    /// Settling steps before recording.
    pub steps_settle: usize,
    /// Recorded physical steps, at least 64; arbitrary lengths are retained.
    pub steps_record: usize,
}

/// Immutable fluid/geometry/units binding for the existing D3Q19 rig.
pub struct MaterialSlotJet3d<'a, F: FluidTransport + FluidAcoustics + ?Sized = GasState> {
    scaling: MaterialScaling<'a, F>,
    acoustic: FluidAcousticState,
    config: SlotJet3dConfig,
    max_speed_m_s: f64,
    force_unit_n: f64,
    span_m: f64,
}

/// Native SGS/source failures and the physical material binding stay distinct.
#[derive(Clone, Debug)]
pub enum MaterialSlotJet3dSgsError {
    Source(SlotJet3dSgsError),
    Material(MaterialAeroacError),
}

impl core::fmt::Display for MaterialSlotJet3dSgsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(error) => write!(f, "{error}"),
            Self::Material(error) => write!(f, "{error}"),
        }
    }
}
impl std::error::Error for MaterialSlotJet3dSgsError {}

impl<'a, F: FluidTransport + FluidAcoustics + ?Sized> MaterialSlotJet3d<'a, F> {
    /// Compile an original gas, liquid or real-gas source into the original rig.
    /// The shared material lattice owner supplies the second-order rate; its
    /// lattice sound speed never substitutes for this fluid's acoustic speed.
    /// Both source views must retain identical density and original phase policy.
    ///
    /// # Errors
    /// Refuses incoherent gas, unsupported rig geometry, allocation overflow,
    /// invalid numerical rates, unrepresentable SI scales, and either physical
    /// or numerical Mach violations before allocating a flow grid.
    pub fn try_new(
        material: &'a F,
        setup: &SlotJet3dSetup,
        physical: JetPhysicalState,
    ) -> Result<Self, MaterialAeroacError> {
        if setup
            .nx
            .checked_mul(setup.ny)
            .and_then(|n| n.checked_mul(setup.nz))
            .is_none_or(|n| n == 0 || n > physical.max_cells)
        {
            return Err(invalid("3-D jet grid exceeds the declared cell budget"));
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
        let config = SlotJet3dConfig {
            nx: setup.nx,
            ny: setup.ny,
            nz: setup.nz,
            slot_half: setup.slot_half,
            u_jet: scaling.velocity_lattice(physical.jet_speed_m_s)?,
            collision: CollisionModel3::CentralMoment {
                second_order_rate: 1.0 / scaling.flow().tau,
                higher_order_rate: setup.higher_order_rate,
            },
            nozzle_thickness: setup.nozzle_thickness,
            edge_distance: setup.edge_distance,
            plate_length: setup.plate_length,
            fringe_width: setup.fringe_width,
            fringe_sigma: setup.fringe_sigma,
            seed_amplitude: setup.seed_amplitude,
            steps_settle: setup.steps_settle,
            steps_record: setup.steps_record,
        };
        validate_config(&config)?;
        // The unchanged rig multiplies its y seed by 1+0.5*sin(2*pi*z/nz).
        let transverse_fraction = 1.5 * setup.seed_amplitude;
        let initial_speed =
            physical.jet_speed_m_s * det::sqrt(1.0 + transverse_fraction * transverse_fraction);
        if initial_speed > physical.max_speed_m_s {
            return Err(invalid(
                "seeded 3-D jet exceeds the declared physical speed envelope",
            ));
        }
        // 3-D cell mass rho*dx^3 times dx/dt, per timestep: rho*dx^4/dt^2 N.
        let velocity_unit = scaling.velocity_si(1.0)?;
        let force_unit_n = scaling.transport().density_kg_m3()
            * velocity_unit
            * velocity_unit
            * physical.cell_size_m
            * physical.cell_size_m;
        let span_m = setup.nz as f64 * physical.cell_size_m;
        if !force_unit_n.is_finite() || force_unit_n <= 0.0 || !span_m.is_finite() || span_m <= 0.0
        {
            return Err(invalid(
                "3-D jet force or span conversion is unrepresentable",
            ));
        }
        Ok(Self {
            scaling,
            acoustic,
            config,
            max_speed_m_s: physical.max_speed_m_s,
            force_unit_n,
            span_m,
        })
    }

    /// Original material and its shared transport/lattice conversion.
    #[must_use]
    pub const fn scaling(&self) -> &MaterialScaling<'a, F> {
        &self.scaling
    }

    /// Frozen acoustic coefficients shared by Mach admission and spectral data.
    #[must_use]
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Derived lattice controls, without mutable fluid-coefficient overrides.
    #[must_use]
    pub const fn lattice_config(&self) -> &SlotJet3dConfig {
        &self.config
    }

    /// Execute the original three-dimensional solver and convert its forces.
    /// The speed check uses the final two diagnostic planes. It does not prove
    /// a full-field or all-times bound, flow convergence, or acoustic accuracy.
    ///
    /// # Errors
    /// Original rig instability, observed envelope violations and nonfinite
    /// physical force refuse without returning a partial physical result.
    pub fn run(&self) -> Result<MaterialSlotJet3dRun<'_, F>, MaterialAeroacError> {
        self.physical_run(run_slot_jet_3d(&self.config)?)
    }

    /// Execute the same material-bound source with an explicitly authored SGS
    /// closure. The filter width is in cells; molecular viscosity, density and
    /// force conversion remain those of this original material. Returned SGS
    /// diagnostics retain their lattice units and exploration-only scope.
    ///
    /// # Errors
    /// Refuses native closure/step/cancellation failures and the same physical
    /// speed-envelope or force-conversion failures as [`Self::run`].
    pub fn run_with_sgs(
        &self,
        closure: SgsConfig,
        gate: &CancelGate,
    ) -> Result<(MaterialSlotJet3dRun<'_, F>, SlotJet3dSgsDiagnostics), MaterialSlotJet3dSgsError>
    {
        let (run, diagnostics) = run_slot_jet_3d_with_sgs(&self.config, closure, gate)
            .map_err(MaterialSlotJet3dSgsError::Source)?;
        let physical = self
            .physical_run(run)
            .map_err(MaterialSlotJet3dSgsError::Material)?;
        Ok((physical, diagnostics))
    }

    fn physical_run(
        &self,
        mut run: SlotJet3dRun,
    ) -> Result<MaterialSlotJet3dRun<'_, F>, MaterialAeroacError> {
        let speed = self.scaling.velocity_si(run.diagnostics.mach_max_lattice)?;
        if speed > self.max_speed_m_s {
            return Err(invalid(
                "final-plane 3-D jet speed exceeds the declared physical envelope",
            ));
        }
        for force in &mut run.force_series {
            for component in force {
                *component *= self.force_unit_n;
                if !component.is_finite() {
                    return Err(invalid("physical 3-D jet force is unrepresentable"));
                }
            }
        }
        Ok(MaterialSlotJet3dRun {
            material: self.scaling.material(),
            acoustic: self.acoustic,
            timestep_s: self.scaling.timestep_s(),
            span_m: self.span_m,
            force_series_n: run.force_series,
            final_plane_physical_mach_max: speed / self.acoustic.sound_speed_m_s(),
            lattice_diagnostics: run.diagnostics,
        })
    }
}

/// Actual physical force history from one periodic three-dimensional cell.
pub struct MaterialSlotJet3dRun<'a, F: FluidTransport + FluidAcoustics + ?Sized = GasState> {
    material: &'a F,
    acoustic: FluidAcousticState,
    timestep_s: f64,
    span_m: f64,
    force_series_n: Vec<[f64; 2]>,
    /// Maximum final-plane speed divided by this fluid's sound speed.
    pub final_plane_physical_mach_max: f64,
    /// Original numerical diagnostics, without relabeling them as SI quantities.
    pub lattice_diagnostics: SlotJet3dDiagnostics,
}

/// Peak force amplitude at one physical frequency, not sound pressure.
#[derive(Clone, Debug)]
pub struct SlotJet3dForceLine {
    /// Physical frequency [Hz].
    pub frequency_hz: f64,
    /// Complex `(Fx, Fy)` force [N] across the recorded periodic span.
    /// Convention `exp(-i omega t)`, phase from the first recorded sample.
    pub force_n: [C64; 2],
    /// Acoustic wavenumber [rad/m] derived from the same material fluid state.
    /// This does not identify a radiation geometry or propagation model.
    pub acoustic_wavenumber_rad_m: f64,
}

/// Cached Hann spectrum of forces integrated over the original periodic span.
/// It retains N, not N/m, and does not declare an isolated acoustic radiator.
pub struct SlotJet3dForceSpectrum {
    acoustic: FluidAcousticState,
    span_m: f64,
    spectrum: ForceSpectrum,
}

impl SlotJet3dForceSpectrum {
    /// Original number of uniformly recorded samples.
    #[must_use]
    pub const fn sample_count(&self) -> usize {
        self.spectrum.sample_count
    }
    /// Original physical sampling interval [s].
    #[must_use]
    pub const fn timestep_s(&self) -> f64 {
        self.spectrum.timestep_s
    }
    /// Original periodic span represented by the integrated force [m].
    #[must_use]
    pub const fn span_m(&self) -> f64 {
        self.span_m
    }
    /// Same frozen fluid snapshot used by the actual flow.
    #[must_use]
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Read a cached interior FFT bin; no transform or allocation is repeated.
    /// DC, Nyquist, invalid bins and unrepresentable wave scales refuse.
    pub fn force_spectral_line(&self, bin: usize) -> Result<SlotJet3dForceLine, AeroacError> {
        let (force_n, frequency_hz) = self.spectrum.line(bin)?;
        force_line(self.acoustic, force_n, frequency_hz)
    }
}

/// Joint DC and harmonic fit of the accepted periodic-cell forces, in N.
///
/// A non-bin fundamental is allowed. This estimates the declared frequencies,
/// not a new source model; no isolated-source radiation or SPL is inferred.
pub struct SlotJet3dForceHarmonics {
    acoustic: FluidAcousticState,
    span_m: f64,
    harmonics: Vec<usize>,
    fundamental_hz: f64,
    force_n: Vec<[C64; 2]>,
    /// Simultaneously fitted constant `(Fx, Fy)` forces [N].
    pub dc_force_n: [f64; 2],
    /// Unweighted sample residual RMS [N]; not a confidence interval.
    pub residual_rms_n: [f64; 2],
}

impl SlotJet3dForceHarmonics {
    /// Harmonic orders, in the original fit plan's order.
    #[must_use]
    pub fn harmonics(&self) -> &[usize] {
        &self.harmonics
    }
    /// Fundamental used by the joint fit [Hz], supplied or locally refined.
    #[must_use]
    pub const fn fundamental_hz(&self) -> f64 {
        self.fundamental_hz
    }
    /// Peak force phasors in N, aligned with `harmonics()`, exp(-i omega t).
    #[must_use]
    pub fn force_n(&self) -> &[[C64; 2]] {
        &self.force_n
    }
    /// Original periodic span represented by the integrated force [m].
    #[must_use]
    pub const fn span_m(&self) -> f64 {
        self.span_m
    }
    /// The source's actual frozen fluid, used for acoustic wavenumbers.
    #[must_use]
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Return a requested harmonic with its same-fluid acoustic wavenumber.
    /// An undeclared harmonic refuses; no spectral interpolation is invented.
    pub fn force_harmonic(&self, harmonic: usize) -> Result<SlotJet3dForceLine, AeroacError> {
        let index =
            self.harmonics
                .binary_search(&harmonic)
                .map_err(|_| AeroacError::InvalidParameter {
                    what: "requested force harmonic was not fitted",
                })?;
        force_line(
            self.acoustic,
            self.force_n[index],
            self.fundamental_hz * harmonic as f64,
        )
    }
}

/// Joint DC and arbitrary-tone fit of the original periodic-cell forces [N].
/// Known frequencies may be inharmonic or belong to independent sources;
/// they are not rounded to harmonic orders or to FFT bins.
pub struct SlotJet3dForceTones {
    acoustic: FluidAcousticState,
    span_m: f64,
    frequencies_hz: Vec<f64>,
    force_n: Vec<[C64; 2]>,
    /// Simultaneously fitted constant `(Fx, Fy)` forces [N].
    pub dc_force_n: [f64; 2],
    /// Unweighted sample residual RMS [N]; not a confidence interval.
    /// Unmodeled tones or frequency errors can also bias fitted amplitudes.
    pub residual_rms_n: [f64; 2],
}

impl SlotJet3dForceTones {
    /// Positive fitted frequencies [Hz], in the admitted increasing order.
    #[must_use]
    pub fn frequencies_hz(&self) -> &[f64] {
        &self.frequencies_hz
    }
    /// Peak force phasors [N], aligned with `frequencies_hz()`, exp(-i omega t).
    /// Phase is measured from the first recorded force sample.
    #[must_use]
    pub fn force_n(&self) -> &[[C64; 2]] {
        &self.force_n
    }
    /// Original periodic span represented by the integrated forces [m].
    #[must_use]
    pub const fn span_m(&self) -> f64 {
        self.span_m
    }
    /// Frozen fluid used by the actual source and acoustic wavenumbers.
    #[must_use]
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Read a fitted tone and its actual same-fluid acoustic wavenumber.
    /// `index` is zero-based in `frequencies_hz()`; an absent index refuses.
    /// The force stays in N; no isolated 3-D radiator or span conversion.
    pub fn force_tone(&self, index: usize) -> Result<SlotJet3dForceLine, AeroacError> {
        let frequency_hz =
            *self
                .frequencies_hz
                .get(index)
                .ok_or(AeroacError::InvalidParameter {
                    what: "requested force tone was not fitted",
                })?;
        force_line(self.acoustic, self.force_n[index], frequency_hz)
    }
}

fn force_line(
    acoustic: FluidAcousticState,
    force_n: [C64; 2],
    frequency_hz: f64,
) -> Result<SlotJet3dForceLine, AeroacError> {
    let acoustic_wavenumber_rad_m =
        core::f64::consts::TAU * frequency_hz / acoustic.sound_speed_m_s();
    if !acoustic_wavenumber_rad_m.is_finite() || acoustic_wavenumber_rad_m <= 0.0 {
        return Err(AeroacError::NonFinite {
            what: "physical acoustic wavenumber is unrepresentable",
        });
    }
    Ok(SlotJet3dForceLine {
        frequency_hz,
        force_n,
        acoustic_wavenumber_rad_m,
    })
}

impl<F: FluidTransport + FluidAcoustics + ?Sized> MaterialSlotJet3dRun<'_, F> {
    /// Splitter `(Fx, Fy)` forces [N] across one periodic span at every dt.
    /// The existing rig does not record a spanwise force component.
    #[must_use]
    pub fn force_series_n(&self) -> &[[f64; 2]] {
        &self.force_series_n
    }
    /// Physical sampling period [s].
    #[must_use]
    pub const fn timestep_s(&self) -> f64 {
        self.timestep_s
    }
    /// Periodic span represented by each integrated force [m].
    #[must_use]
    pub const fn span_m(&self) -> f64 {
        self.span_m
    }
    /// Original source used by the flow, retaining its material receipts/data.
    #[must_use]
    pub const fn material(&self) -> &F {
        self.material
    }
    /// Exact frozen acoustic state used by the flow's physical Mach admission.
    #[must_use]
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }
    /// Original rig's spectral-shape/scaling boundary; no SPL is inferred.
    #[must_use]
    pub const fn scope(&self) -> &'static str {
        SCOPE_STATEMENT
    }

    /// Extract a force bin with the same mean removal, symmetric Hann window,
    /// coherent-gain normalization and phase as the two-dimensional observer.
    /// The integrated force remains a periodic source datum; neither a 2-D
    /// Green function nor an isolated 3-D dipole is applied to it.
    ///
    /// # Errors
    /// DC, Nyquist, invalid bins or unrepresentable frequency/force refuse.
    pub fn force_spectral_line(&self, bin: usize) -> Result<SlotJet3dForceLine, AeroacError> {
        check_force_bin(self.force_series_n.len(), bin)?;
        self.prepare_force_spectrum(self.force_series_n.len())?
            .force_spectral_line(bin)
    }

    /// Transform the accepted record once for many frequency-bin queries.
    /// Checks the caller's sample bound before O(N) storage/O(N log N) work.
    pub fn prepare_force_spectrum(
        &self,
        max_samples: usize,
    ) -> Result<SlotJet3dForceSpectrum, AeroacError> {
        Ok(SlotJet3dForceSpectrum {
            acoustic: self.acoustic,
            span_m: self.span_m,
            spectrum: ForceSpectrum::new(&self.force_series_n, self.timestep_s, max_samples)?,
        })
    }

    /// Jointly fit DC and requested harmonics with the shared L1 QR plan.
    /// The plan must match the actual sample count and physical clock exactly;
    /// both recorded force axes share its factorization and scratch. No force
    /// per span, acoustic pressure or absent spanwise component is invented.
    pub fn fit_force_harmonics(
        &self,
        plan: &HarmonicFitPlan,
    ) -> Result<SlotJet3dForceHarmonics, AeroacError> {
        let (force_n, summary) = fit_force_harmonics(&self.force_series_n, self.timestep_s, plan)?;
        Ok(SlotJet3dForceHarmonics {
            acoustic: self.acoustic,
            span_m: self.span_m,
            harmonics: plan.harmonics().to_vec(),
            fundamental_hz: plan.fundamental_hz(),
            force_n,
            dc_force_n: summary.map(|s| s.dc),
            residual_rms_n: summary.map(|s| s.residual_rms),
        })
    }

    /// Refine one common fundamental from both recorded integrated-force axes.
    /// Original N units, sample clock, phase origin and periodic span survive.
    /// No absent spanwise force or isolated-source radiation is manufactured.
    /// The shared local vector fit retains physical component weights and
    /// explicitly budgets its transposed record and numerical workspace.
    pub fn fit_force_harmonics_refined(
        &self,
        initial_fundamental_hz: f64,
        harmonics: &[usize],
        budget: fs_fft::harmonic::HarmonicFrequencyBudget,
    ) -> Result<
        (SlotJet3dForceHarmonics, super::ForceFrequencyRefinement),
        fs_fft::harmonic::HarmonicFrequencyError,
    > {
        let (force_n, summary, refinement) = super::frequency::refine_force_harmonics(
            &self.force_series_n,
            self.timestep_s,
            initial_fundamental_hz,
            harmonics,
            budget,
        )?;
        Ok((
            SlotJet3dForceHarmonics {
                acoustic: self.acoustic,
                span_m: self.span_m,
                harmonics: harmonics.to_vec(),
                fundamental_hz: refinement.fundamental_hz,
                force_n,
                dc_force_n: summary.map(|s| s.dc),
                residual_rms_n: summary.map(|s| s.residual_rms),
            },
            refinement,
        ))
    }

    /// Jointly fit arbitrary declared tones on the original physical clock.
    /// One shared L1 factorization separates DC and all tones for both force
    /// axes. Noncommensurate lines retain their true frequencies, phase from
    /// sample zero, and forces integrated over the actual periodic span.
    /// No noise, stationarity or absent spatial-radiation model is inferred.
    pub fn fit_force_tones(&self, plan: &ToneFitPlan) -> Result<SlotJet3dForceTones, AeroacError> {
        let (force_n, summary) = fit_force_tones(&self.force_series_n, self.timestep_s, plan)?;
        Ok(SlotJet3dForceTones {
            acoustic: self.acoustic,
            span_m: self.span_m,
            frequencies_hz: plan.frequencies_hz().to_vec(),
            force_n,
            dc_force_n: summary.map(|s| s.dc),
            residual_rms_n: summary.map(|s| s.residual_rms),
        })
    }

    /// Jointly refine independent frequencies from the integrated force axes.
    /// Local nonoverlapping intervals preserve each tone's identity. The
    /// common L1 owner uses both physical components without renormalizing
    /// their relative strengths or imposing harmonic spacing. Original N
    /// units, sample-zero phase, material and periodic span survive unchanged.
    /// No absent spanwise force or isolated-source radiation is inferred.
    ///
    /// The complete workspace includes the transposed input. Unresolved
    /// frequencies, Nyquist violations and incomplete convergence refuse.
    pub fn fit_force_tones_refined(
        &self,
        initial_frequencies_hz: &[f64],
        frequency_bounds_hz: &[[f64; 2]],
        budget: fs_fft::harmonic::ToneFrequencyBudget,
    ) -> Result<
        (SlotJet3dForceTones, super::ForceToneFrequencyRefinement),
        fs_fft::harmonic::ToneFrequencyError,
    > {
        let fit = super::frequency::refine_force_tones(
            &self.force_series_n,
            self.timestep_s,
            initial_frequencies_hz,
            frequency_bounds_hz,
            budget,
        )?;
        Ok((
            SlotJet3dForceTones {
                acoustic: self.acoustic,
                span_m: self.span_m,
                frequencies_hz: fit.frequencies_hz,
                force_n: fit.force,
                dc_force_n: fit.summary.map(|s| s.dc),
                residual_rms_n: fit.summary.map(|s| s.residual_rms),
            },
            fit.refinement,
        ))
    }
}

impl MaterialSlotJet3dRun<'_, GasState> {
    /// Original gas state for gas-bound runs, including mixture properties.
    #[must_use]
    pub const fn gas(&self) -> &GasState {
        self.material
    }
}

#[cfg(test)]
mod tone_tests {
    use super::*;
    use fs_fft::harmonic::HarmonicFitBudget;

    #[test]
    fn g1_refined_independent_tones_preserve_periodic_cell_force_and_span() {
        let gas = GasState::try_new_moist_air(300.0, 93_000.0, 0.2).unwrap();
        let dt = 0.00025;
        let frequencies = [137.3, 211.7];
        let expected = [
            [C64::new(0.02, -0.01), C64::new(-0.004, 0.008)],
            [C64::new(0.003, 0.009), C64::new(0.014, -0.002)],
        ];
        let record = (0..1024)
            .map(|i| {
                let mut force = [0.05, -0.03];
                for (frequency, phasors) in frequencies.iter().zip(expected) {
                    let (s, c) = (core::f64::consts::TAU * frequency * i as f64 * dt).sin_cos();
                    for axis in 0..2 {
                        force[axis] += phasors[axis].re * c + phasors[axis].im * s;
                    }
                }
                force
            })
            .collect();
        // A manufactured force history exercises the material-run observation
        // seam without making a flow-convergence or isolated-radiator claim.
        let run = MaterialSlotJet3dRun {
            material: &gas,
            acoustic: gas.acoustic_state().unwrap(),
            timestep_s: dt,
            span_m: 0.013,
            force_series_n: record,
            final_plane_physical_mach_max: 0.0,
            lattice_diagnostics: SlotJet3dDiagnostics {
                mach_max_lattice: 0.0,
                flux_plate_plane: 0.0,
                flux_fringe_plane: 0.0,
                reynolds: 0.0,
                record_len: 1024,
                strouhal_bin_width: 0.0,
            },
        };
        let (fitted, diagnostics) = run
            .fit_force_tones_refined(
                &[137.7, 211.3],
                &[[136.0, 139.0], [210.0, 213.0]],
                fs_fft::harmonic::ToneFrequencyBudget {
                    absolute_tolerance_hz: 1e-8,
                    max_evaluations: 32,
                    min_frequency_reciprocal_condition: 1e-9,
                    fit: HarmonicFitBudget {
                        max_workspace_bytes: 1_000_000,
                        min_reciprocal_condition: 1e-10,
                    },
                },
            )
            .unwrap();
        assert_eq!(fitted.span_m(), run.span_m());
        assert_eq!(fitted.acoustic(), run.acoustic());
        assert!(diagnostics.evaluations > 1 && diagnostics.evaluations <= 32);
        assert!((fitted.dc_force_n[0] - 0.05).abs() < 1e-10);
        assert!((fitted.dc_force_n[1] + 0.03).abs() < 1e-10);
        assert!(fitted.residual_rms_n.iter().all(|x| *x < 1e-10));
        for (index, frequency) in frequencies.iter().enumerate() {
            let line = fitted.force_tone(index).unwrap();
            assert!((line.frequency_hz - frequency).abs() < 1e-8);
            assert_eq!(
                line.acoustic_wavenumber_rad_m,
                core::f64::consts::TAU * line.frequency_hz / gas.sound_speed
            );
            for axis in 0..2 {
                assert!((line.force_n[axis] - expected[index][axis]).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn g1_refined_harmonics_retain_periodic_span_and_integrated_newtons() {
        let dt = 0.00025;
        let f0 = 137.3;
        let expected = [C64::new(0.02, -0.01), C64::new(-0.004, 0.008)];
        let record: Vec<_> = (0..1021)
            .map(|i| {
                let (s, c) = (core::f64::consts::TAU * f0 * i as f64 * dt).sin_cos();
                [
                    0.05 + expected[0].re * c + expected[0].im * s,
                    -0.03 + expected[1].re * c + expected[1].im * s,
                ]
            })
            .collect();
        let (force_n, summary, refinement) = super::super::frequency::refine_force_harmonics(
            &record,
            dt,
            137.8,
            &[1],
            fs_fft::harmonic::HarmonicFrequencyBudget {
                frequency_bounds_hz: [136.0, 139.0],
                absolute_tolerance_hz: 1e-8,
                max_evaluations: 24,
                fit: HarmonicFitBudget {
                    max_workspace_bytes: 1_000_000,
                    min_reciprocal_condition: 1e-10,
                },
            },
        )
        .unwrap();
        let fitted = SlotJet3dForceHarmonics {
            acoustic: FluidAcousticState::try_new(1.2, 340.0).unwrap(),
            span_m: 0.013,
            harmonics: vec![1],
            fundamental_hz: refinement.fundamental_hz,
            force_n,
            dc_force_n: summary.map(|s| s.dc),
            residual_rms_n: summary.map(|s| s.residual_rms),
        };
        assert_eq!(fitted.span_m(), 0.013);
        let line = fitted.force_harmonic(1).unwrap();
        assert!((line.frequency_hz - f0).abs() < 1e-8);
        assert!(
            (line.acoustic_wavenumber_rad_m - core::f64::consts::TAU * f0 / 340.0).abs() < 1e-9
        );
        for axis in 0..2 {
            assert!((line.force_n[axis] - expected[axis]).abs() < 1e-10);
            assert!(fitted.residual_rms_n[axis] < 1e-10);
        }
        assert!((fitted.dc_force_n[0] - 0.05).abs() < 1e-10);
        assert!((fitted.dc_force_n[1] + 0.03).abs() < 1e-10);
    }

    #[test]
    fn g1_arbitrary_tones_retain_integrated_force_and_periodic_span() {
        let n = 512;
        let dt = 0.00025;
        let frequencies = [137.3, 211.7];
        let expected = [
            [C64::new(0.02, -0.01), C64::new(-0.004, 0.008)],
            [C64::new(-0.003, 0.005), C64::new(0.006, 0.002)],
        ];
        let record: Vec<_> = (0..n)
            .map(|i| {
                let mut force = [0.05, -0.03];
                for (frequency, phasors) in frequencies.iter().zip(expected) {
                    let phase = core::f64::consts::TAU * frequency * i as f64 * dt;
                    for axis in 0..2 {
                        force[axis] +=
                            phasors[axis].re * phase.cos() + phasors[axis].im * phase.sin();
                    }
                }
                force
            })
            .collect();
        let plan = ToneFitPlan::new(
            n,
            dt,
            &frequencies,
            HarmonicFitBudget {
                max_workspace_bytes: 1_000_000,
                min_reciprocal_condition: 1e-10,
            },
        )
        .unwrap();
        let (force_n, summary) = fit_force_tones(&record, dt, &plan).unwrap();
        let fitted = SlotJet3dForceTones {
            acoustic: FluidAcousticState::try_new(1.2, 340.0).unwrap(),
            span_m: 0.013,
            frequencies_hz: frequencies.to_vec(),
            force_n,
            dc_force_n: summary.map(|s| s.dc),
            residual_rms_n: summary.map(|s| s.residual_rms),
        };
        assert_eq!(fitted.span_m(), 0.013);
        assert_eq!(fitted.frequencies_hz(), frequencies);
        for (index, amplitude) in expected.iter().enumerate() {
            let line = fitted.force_tone(index).unwrap();
            assert_eq!(line.frequency_hz, frequencies[index]);
            assert_eq!(
                line.acoustic_wavenumber_rad_m,
                core::f64::consts::TAU * frequencies[index] / 340.0
            );
            for axis in 0..2 {
                assert!((line.force_n[axis] - amplitude[axis]).abs() < 1e-12);
                assert!(fitted.residual_rms_n[axis] < 1e-12);
            }
        }
        assert!((fitted.dc_force_n[0] - 0.05).abs() < 1e-12);
        assert!((fitted.dc_force_n[1] + 0.03).abs() < 1e-12);
        assert!(fitted.force_tone(2).is_err());
        assert!(fitted.force_tone(usize::MAX).is_err());
    }
}
