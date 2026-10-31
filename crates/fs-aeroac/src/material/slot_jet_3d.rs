//! One physical fluid state drives the original three-dimensional slot-jet rig.
//!
//! The central-moment viscosity, SI force and acoustic wavenumber share the
//! same gas, liquid or real-gas source. Forces are integrated across one periodic span, in N; they are
//! not the N/m line sources of the two-dimensional Curle observer. This module
//! deliberately leaves spatial radiation to an appropriate acoustic boundary.

use super::{
    JetPhysicalState, MAX_PHYSICAL_MACH, MaterialAeroacError, acoustic_from_scaling,
    force_spectral_line, invalid,
};
use crate::slot_jet_3d::{SlotJet3dConfig, SlotJet3dDiagnostics, run_slot_jet_3d, validate_config};
use crate::{AeroacError, SCOPE_STATEMENT};
use fs_lbm::d3q19::CollisionModel3;
use fs_lbm::material::{MaterialScaling, plan_scaling_from_material};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::FluidTransport;
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};

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
    /// Recorded steps; a power of two, at least 64.
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
        let mut run = run_slot_jet_3d(&self.config)?;
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

/// Hann-windowed force amplitude at one physical frequency, not sound pressure.
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
        let (force_n, frequency_hz) =
            force_spectral_line(&self.force_series_n, self.timestep_s, bin)?;
        let acoustic_wavenumber_rad_m =
            core::f64::consts::TAU * frequency_hz / self.acoustic.sound_speed_m_s();
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
}

impl MaterialSlotJet3dRun<'_, GasState> {
    /// Original gas state for gas-bound runs, including mixture properties.
    #[must_use]
    pub const fn gas(&self) -> &GasState {
        self.material
    }
}
