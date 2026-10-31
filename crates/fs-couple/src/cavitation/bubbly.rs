//! Small-amplitude coherent waves in a dilute monodisperse bubbly liquid.
//!
//! The original Keller--Miksis owner supplies its reference-equilibrium
//! linearization. For exp(-i omega t), its fractional-radius susceptibility is
//! `chi = -(1-i omega tau)/(rho R0^2 (omega0^2-a omega^2-i b omega))`.
//! The point-monopole scattering amplitude is `f = -rho R0^3 omega^2 chi`.
//! Foldy's coherent dispersion law then gives
//! `k^2 = omega^2/c_l^2 + 4 pi n f`, where `n = 3 alpha/(4 pi R0^3)`.
//! See Foldy (1945), equations 15 and 54:
//! <https://doi.org/10.1103/PhysRev.67.107>.
//!
//! Radiation back-reaction, viscosity, gas stiffness and surface tension come
//! from that native bubble. In particular, pressure-derivative forcing is
//! retained; adding a fitted damping coefficient would change the model.
//! The bubble mass is negligible and the mixture momentum density is rho_l
//! to leading dilute order. This is not the finite-mass frozen Wood closure.
//! Thermal gas relaxation, vapor mass transfer, slip, nonlinear oscillations,
//! cloud collapse and spatial correlations are outside this model.
//!
//! Small alpha alone does not justify propagation close to resonance. Besides
//! compactness and long-wavelength homogenization, this owner checks the
//! neighboring scattered-pressure ratio |f|/d and the multiple-scattering
//! parameter n sigma_s/k_l, with d=n^(-1/3), sigma_s=4 pi |f|^2. These follow
//! Commander and Prosperetti (1989), section VI, equations 64--66:
//! <https://pages.jh.edu/aprospe1/publications/PapersPublished/BubblyLiquids/Commander.pdf>.
//! The declared ceilings are conservative admission choices, not error bounds.
//! Coherent attenuation includes scattering out of the mean wave; it must not
//! be deposited wholesale into a thermal body's heat ledger.

use fs_duct::compose::{ComposeOptions, PreparedChain};
use fs_duct::fluid::{self, FluidDuctError};
use fs_duct::{ConeSlices, Duct, SegmentWave, WallModel};
use fs_flux::cavitation::{BubbleLinearization, CavitationError, KellerMiksisBubble};
use fs_material::MaterialError;
use fs_material::cavitation::{CavitationProperties, CavitationState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::{c64::C64, det};

/// Declared physical bounds; no value is a damping fit or a numerical clamp.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubblyLiquidLimits {
    /// Largest equilibrium bubble volume fraction, in (0, 0.01].
    pub maximum_void_fraction: f64,
    /// Largest max(k_l, |k|) R0, in (0, 0.1].
    pub maximum_wavenumber_radius: f64,
    /// Largest max(k_l, |k|) d, in (0, 1]. The characteristic wavelength and
    /// attenuation length must exceed the mean interbubble spacing.
    pub maximum_wavenumber_spacing: f64,
    /// Largest |delta R|/R0 under the declared LOCAL pressure, in (0, 0.05].
    pub maximum_fractional_radius: f64,
    /// Largest |f|/d, in (0, 0.1].
    pub maximum_neighbor_scattered_pressure_ratio: f64,
    /// Largest n sigma_s/k_l, in (0, 0.1].
    pub maximum_multiple_scattering_parameter: f64,
}

/// Original source/solver errors or a named physical admission failure.
#[derive(Clone, Debug, PartialEq)]
pub enum BubblyLiquidError {
    Material(MaterialError),
    Bubble(CavitationError),
    Duct(FluidDuctError),
    MaterialMismatch,
    InvalidInput(&'static str),
    OutsideModel {
        quantity: &'static str,
        value: f64,
        maximum: f64,
    },
}

impl core::fmt::Display for BubblyLiquidError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "bubbly liquid material: {error}"),
            Self::Bubble(error) => write!(f, "bubbly liquid dynamics: {error}"),
            Self::Duct(error) => write!(f, "bubbly liquid propagation: {error}"),
            Self::MaterialMismatch => write!(f, "bubble and source must use the identical liquid"),
            Self::InvalidInput(what) => write!(f, "invalid bubbly liquid input: {what}"),
            Self::OutsideModel {
                quantity,
                value,
                maximum,
            } => {
                write!(f, "bubbly liquid {quantity}={value} exceeds {maximum}")
            }
        }
    }
}
impl std::error::Error for BubblyLiquidError {}
impl From<MaterialError> for BubblyLiquidError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<CavitationError> for BubblyLiquidError {
    fn from(error: CavitationError) -> Self {
        Self::Bubble(error)
    }
}
impl From<FluidDuctError> for BubblyLiquidError {
    fn from(error: FluidDuctError) -> Self {
        Self::Duct(error)
    }
}

fn bound(quantity: &'static str, value: f64, maximum: f64) -> Result<(), BubblyLiquidError> {
    if !value.is_finite() || value < 0.0 || value > maximum {
        return Err(BubblyLiquidError::OutsideModel {
            quantity,
            value,
            maximum,
        });
    }
    Ok(())
}
fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

impl BubblyLiquidLimits {
    fn validate(self) -> Result<(), BubblyLiquidError> {
        for (name, value, ceiling) in [
            ("void-fraction limit", self.maximum_void_fraction, 0.01),
            ("compactness limit", self.maximum_wavenumber_radius, 0.1),
            ("homogenization limit", self.maximum_wavenumber_spacing, 1.0),
            ("linear-radius limit", self.maximum_fractional_radius, 0.05),
            (
                "neighbor-scattering limit",
                self.maximum_neighbor_scattered_pressure_ratio,
                0.1,
            ),
            (
                "multiple-scattering limit",
                self.maximum_multiple_scattering_parameter,
                0.1,
            ),
        ] {
            if !value.is_finite() || value <= 0.0 || value > ceiling {
                return Err(BubblyLiquidError::InvalidInput(name));
            }
        }
        Ok(())
    }
}

/// One original liquid source, native compressible bubble and bubble population.
/// The coefficients and population are frozen for this owner's lifetime.
pub struct DiluteBubblyLiquid<'a, F: CavitationProperties + ?Sized> {
    source: &'a F,
    material: CavitationState,
    bubble: KellerMiksisBubble,
    linearization: Option<BubbleLinearization>,
    void_fraction: f64,
    limits: BubblyLiquidLimits,
}

impl<'a, F: CavitationProperties + ?Sized> DiluteBubblyLiquid<'a, F> {
    pub fn new(
        source: &'a F,
        bubble: KellerMiksisBubble,
        void_fraction: f64,
        limits: BubblyLiquidLimits,
    ) -> Result<Self, BubblyLiquidError> {
        limits.validate()?;
        bound("void fraction", void_fraction, limits.maximum_void_fraction)?;
        let material = source.cavitation_state()?;
        if material != bubble.liquid() {
            return Err(BubblyLiquidError::MaterialMismatch);
        }
        let linearization = if void_fraction == 0.0 {
            None
        } else {
            Some(bubble.linearization()?)
        };
        Ok(Self {
            source,
            material,
            bubble,
            linearization,
            void_fraction,
            limits,
        })
    }

    pub const fn source(&self) -> &'a F {
        self.source
    }
    pub const fn material(&self) -> CavitationState {
        self.material
    }
    pub const fn bubble(&self) -> &KellerMiksisBubble {
        &self.bubble
    }
    pub const fn void_fraction(&self) -> f64 {
        self.void_fraction
    }
    pub const fn linearization(&self) -> Option<BubbleLinearization> {
        self.linearization
    }

    /// Admit one harmonic and its maximum pressure amplitude ANYWHERE in the
    /// bubbly region. This is not an inlet-pressure bound: duct resonances and
    /// reflections can amplify local pressure. The downstream field owner must
    /// verify that its solved local pressures respect this declared envelope.
    pub fn response(
        &self,
        omega_rad_s: f64,
        maximum_local_pressure_amplitude_pa: f64,
    ) -> Result<BubblyLiquidResponse, BubblyLiquidError> {
        if !omega_rad_s.is_finite()
            || omega_rad_s <= 0.0
            || !maximum_local_pressure_amplitude_pa.is_finite()
            || maximum_local_pressure_amplitude_pa < 0.0
        {
            return Err(BubblyLiquidError::InvalidInput(
                "frequency and local pressure amplitude",
            ));
        }
        let medium = self.material.acoustic_state()?;
        let omega = omega_rad_s;
        let k_l = omega / medium.sound_speed_m_s();
        if !k_l.is_finite() || k_l <= 0.0 {
            return Err(BubblyLiquidError::InvalidInput(
                "unrepresentable liquid wavenumber",
            ));
        }
        let mut diagnostics = BubblyResponseDiagnostics {
            void_fraction: self.void_fraction,
            maximum_local_pressure_amplitude_pa,
            fractional_radius_per_pa: C64::ZERO,
            fractional_radius_amplitude: 0.0,
            wall_mach_amplitude: 0.0,
            mean_bubble_spacing_m: None,
            wavenumber_radius: 0.0,
            wavenumber_spacing: 0.0,
            neighbor_scattered_pressure_ratio: 0.0,
            multiple_scattering_parameter: 0.0,
        };
        let Some(linear) = self.linearization else {
            return Ok(BubblyLiquidResponse {
                medium,
                omega_rad_s,
                wavenumber: C64::from_re(k_l),
                diagnostics,
            });
        };
        let radius = linear.equilibrium_radius_m;
        let density = medium.density_kg_m3();
        let inertia = density * radius * radius;
        let numerator = C64::new(1.0, -omega * linear.pressure_derivative_time_s);
        let denominator = C64::new(
            linear.stiffness_coefficient_per_s2 - linear.inertia_multiplier * omega * omega,
            -linear.damping_coefficient_per_s * omega,
        );
        let susceptibility = -(numerator / denominator).scale(1.0 / inertia);
        if !finite(susceptibility) || susceptibility.abs() == 0.0 || susceptibility.im >= 0.0 {
            return Err(BubblyLiquidError::InvalidInput(
                "unresolved passive bubble susceptibility",
            ));
        }
        diagnostics.fractional_radius_per_pa = susceptibility;
        diagnostics.fractional_radius_amplitude =
            susceptibility.abs() * maximum_local_pressure_amplitude_pa;
        bound(
            "fractional radius",
            diagnostics.fractional_radius_amplitude,
            self.limits.maximum_fractional_radius,
        )?;
        let radial_limits = self.bubble.limits();
        let radius_envelope = (1.0 - radial_limits.minimum_radius_m / radius)
            .min(radial_limits.maximum_radius_m / radius - 1.0);
        bound(
            "native radius interval",
            diagnostics.fractional_radius_amplitude,
            radius_envelope,
        )?;
        diagnostics.wall_mach_amplitude = diagnostics.fractional_radius_amplitude * k_l * radius;
        bound(
            "wall Mach amplitude",
            diagnostics.wall_mach_amplitude,
            radial_limits.maximum_wall_mach_number,
        )?;
        // This is the native linearized wall-minus-drive pressure, not a
        // second radiation-power subtraction or a heuristic damping term.
        let wall_stiffness = C64::new(
            inertia * linear.stiffness_coefficient_per_s2,
            -4.0 * self.material.dynamic_viscosity_pa_s() * omega,
        );
        let pressure_ratio = (C64::ONE + wall_stiffness * susceptibility).abs()
            * (maximum_local_pressure_amplitude_pa / medium.adiabatic_bulk_modulus_pa());
        bound(
            "native liquid pressure ratio",
            pressure_ratio,
            self.bubble.compressibility_limits().maximum_pressure_ratio,
        )?;
        let inverse_sound_speed = 1.0 / medium.sound_speed_m_s();
        let slowness_squared = C64::from_re(inverse_sound_speed * inverse_sound_speed)
            - susceptibility.scale(3.0 * self.void_fraction * density);
        let wavenumber = slowness_squared.sqrt().scale(omega);
        if !finite(wavenumber) || wavenumber.re <= 0.0 || wavenumber.im <= 0.0 {
            return Err(BubblyLiquidError::InvalidInput(
                "unresolved forward passive mixture wave",
            ));
        }
        let largest_wavenumber = k_l.max(wavenumber.abs());
        diagnostics.wavenumber_radius = largest_wavenumber * radius;
        bound(
            "bubble compactness",
            diagnostics.wavenumber_radius,
            self.limits.maximum_wavenumber_radius,
        )?;
        let spacing_over_radius = det::pow(
            4.0 * core::f64::consts::PI / (3.0 * self.void_fraction),
            1.0 / 3.0,
        );
        let spacing = radius * spacing_over_radius;
        if !spacing.is_finite() || spacing <= 0.0 {
            return Err(BubblyLiquidError::InvalidInput(
                "unrepresentable interbubble spacing",
            ));
        }
        diagnostics.mean_bubble_spacing_m = Some(spacing);
        diagnostics.wavenumber_spacing = largest_wavenumber * spacing;
        bound(
            "homogenization",
            diagnostics.wavenumber_spacing,
            self.limits.maximum_wavenumber_spacing,
        )?;
        let scattering_over_radius = (numerator / denominator).scale(omega * omega).abs();
        diagnostics.neighbor_scattered_pressure_ratio =
            scattering_over_radius / spacing_over_radius;
        bound(
            "neighbor scattered pressure",
            diagnostics.neighbor_scattered_pressure_ratio,
            self.limits.maximum_neighbor_scattered_pressure_ratio,
        )?;
        diagnostics.multiple_scattering_parameter =
            3.0 * self.void_fraction * scattering_over_radius * scattering_over_radius
                / (k_l * radius);
        bound(
            "multiple scattering",
            diagnostics.multiple_scattering_parameter,
            self.limits.maximum_multiple_scattering_parameter,
        )?;
        Ok(BubblyLiquidResponse {
            medium,
            omega_rad_s,
            wavenumber,
            diagnostics,
        })
    }
}

/// Physical envelope and response under the admitted local pressure amplitude.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubblyResponseDiagnostics {
    pub void_fraction: f64,
    pub maximum_local_pressure_amplitude_pa: f64,
    /// Complex (delta R/R0)/delta p in Pa^-1, under exp(-i omega t).
    pub fractional_radius_per_pa: C64,
    pub fractional_radius_amplitude: f64,
    pub wall_mach_amplitude: f64,
    /// None at the exact pure-liquid endpoint.
    pub mean_bubble_spacing_m: Option<f64>,
    pub wavenumber_radius: f64,
    pub wavenumber_spacing: f64,
    pub neighbor_scattered_pressure_ratio: f64,
    pub multiple_scattering_parameter: f64,
}

/// An immutable, admitted harmonic response; density and k remain inseparable.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubblyLiquidResponse {
    medium: FluidAcousticState,
    omega_rad_s: f64,
    wavenumber: C64,
    diagnostics: BubblyResponseDiagnostics,
}

impl BubblyLiquidResponse {
    pub const fn liquid(&self) -> FluidAcousticState {
        self.medium
    }
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }
    pub const fn wavenumber_per_m(&self) -> C64 {
        self.wavenumber
    }
    pub const fn diagnostics(&self) -> BubblyResponseDiagnostics {
        self.diagnostics
    }

    /// The native p/volume-flow wave; coherent attenuation is not all heat.
    pub fn segment_wave(&self, duct_radius_m: f64) -> Result<SegmentWave, BubblyLiquidError> {
        Ok(fluid::segment_wave_with_wavenumber(
            &self.medium,
            duct_radius_m,
            self.omega_rad_s,
            self.wavenumber,
        )?)
    }

    /// Reuse native analytic/Magnus propagation and stable scattering assembly.
    ///
    /// The geometry must describe a representative bulk bubbly region in the
    /// duct's plane-wave regime; wall-near bubbles are not independently solved.
    /// Local field amplitudes must obey the pressure envelope in diagnostics.
    pub fn duct_chain(
        &self,
        duct: &Duct,
        wall: WallModel<'_>,
        cone_slices: ConeSlices,
        options: ComposeOptions,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedChain, BubblyLiquidError> {
        Ok(fluid::duct_chain_with_wavenumber(
            duct,
            &self.medium,
            self.omega_rad_s,
            self.wavenumber,
            wall,
            cone_slices,
            options,
            cancelled,
        )?)
    }
}

#[cfg(test)]
#[path = "bubbly/tests.rs"]
mod tests;
