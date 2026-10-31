//! Coherent waves from the complete native equilibrium bubble distribution.
//!
//! The occupied-volume susceptibility is 3 sum_j alpha_j chi_j, so
//! k^2 = omega^2/c_l^2 - 3 rho_l omega^2 sum_j alpha_j chi_j.
//! Number weights and volume weights are never interchanged. The scattering
//! cross section is summed as sum_j n_j 4 pi |f_j|^2; unlike coherent forward
//! amplitudes, scattered power cannot cancel between radius classes.
//!
//! Every occupied class uses the same single-bubble response and admission
//! kernel as DiluteBubblyLiquid. The conservative neighbor check uses the
//! largest |f_j| divided by the total mean spacing. Splitting one histogram
//! class into several identical classes therefore cannot relax physical limits.
//! This remains a frozen, dilute, uncorrelated linear ensemble, not interacting
//! cloud dynamics, a distribution-fitting routine or a thermal model.

use super::*;
use fs_flux::cavitation::BubblePopulation;

/// One source and its original, number-conserving native bubble distribution.
pub struct PolydisperseBubblyLiquid<'a, F: CavitationProperties + ?Sized> {
    source: &'a F,
    material: CavitationState,
    population: BubblePopulation,
    limits: BubblyLiquidLimits,
}

impl<'a, F: CavitationProperties + ?Sized> PolydisperseBubblyLiquid<'a, F> {
    pub fn new(
        source: &'a F,
        population: BubblePopulation,
        limits: BubblyLiquidLimits,
    ) -> Result<Self, BubblyLiquidError> {
        limits.validate()?;
        let material = source.cavitation_state()?;
        if material != population.liquid() {
            return Err(BubblyLiquidError::MaterialMismatch);
        }
        bound(
            "void fraction",
            population.void_fraction(),
            limits.maximum_void_fraction,
        )?;
        Ok(Self {
            source,
            material,
            population,
            limits,
        })
    }

    pub const fn source(&self) -> &'a F {
        self.source
    }
    pub const fn material(&self) -> CavitationState {
        self.material
    }
    pub const fn population(&self) -> &BubblePopulation {
        &self.population
    }

    /// Admit the complete distribution under one maximum LOCAL pressure
    /// envelope. This bound applies everywhere in the eventual solved field.
    /// Preparation allocates no per-frequency storage and never changes weights.
    pub fn response(
        &self,
        omega_rad_s: f64,
        maximum_local_pressure_amplitude_pa: f64,
    ) -> Result<PolydisperseBubblyLiquidResponse, BubblyLiquidError> {
        self.response_controlled(
            omega_rad_s,
            maximum_local_pressure_amplitude_pa,
            &mut || false,
        )
    }

    /// Fixed-order compensated accumulation with a cancellation check per class.
    /// A response is returned only after all classes and aggregate limits admit.
    pub fn response_controlled(
        &self,
        omega_rad_s: f64,
        maximum_local_pressure_amplitude_pa: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PolydisperseBubblyLiquidResponse, BubblyLiquidError> {
        poll(cancelled)?;
        let (medium, k_l) = admit_harmonic(
            self.material,
            omega_rad_s,
            maximum_local_pressure_amplitude_pa,
        )?;
        let mut diagnostics = PopulationResponseDiagnostics {
            void_fraction: self.population.void_fraction(),
            number_density_per_m3: self.population.number_density_per_m3(),
            occupied_classes: self.population.occupied_classes(),
            maximum_local_pressure_amplitude_pa,
            void_fraction_per_pa: C64::ZERO,
            maximum_fractional_radius_amplitude: 0.0,
            maximum_wall_mach_amplitude: 0.0,
            maximum_liquid_pressure_ratio: 0.0,
            mean_bubble_spacing_m: None,
            wavenumber_radius: 0.0,
            wavenumber_spacing: 0.0,
            neighbor_scattered_pressure_ratio: 0.0,
            multiple_scattering_parameter: 0.0,
        };
        if self.population.occupied_classes() == 0 {
            poll(cancelled)?;
            return Ok(PolydisperseBubblyLiquidResponse {
                medium,
                omega_rad_s,
                wavenumber: C64::from_re(k_l),
                diagnostics,
            });
        }
        let mut compliance_real = Sum::default();
        let mut compliance_imaginary = Sum::default();
        let mut scattering = Sum::default();
        let mut maximum_radius: f64 = 0.0;
        let mut maximum_scattering_amplitude: f64 = 0.0;
        for class in self.population.classes() {
            poll(cancelled)?;
            let Some(linear) = class.linearization() else {
                continue;
            };
            let response = single_bubble_response(
                class.bubble(),
                linear,
                medium,
                omega_rad_s,
                maximum_local_pressure_amplitude_pa,
                self.limits,
            )?;
            let radius = linear.equilibrium_radius_m;
            let volume_response = response.susceptibility.scale(3.0 * class.void_fraction());
            compliance_real.add(volume_response.re);
            compliance_imaginary.add(volume_response.im);
            // n_j 4 pi |f_j|^2 = 3 alpha_j |f_j/R_j|^2/R_j.
            // This form avoids cubing radii again and uses the exact admitted
            // class volume; it does not square the coherent amplitude sum.
            scattering.add(
                3.0 * class.void_fraction()
                    * response.scattering_over_radius
                    * response.scattering_over_radius
                    / radius,
            );
            maximum_radius = maximum_radius.max(radius);
            maximum_scattering_amplitude =
                maximum_scattering_amplitude.max(response.scattering_over_radius * radius);
            diagnostics.maximum_fractional_radius_amplitude = diagnostics
                .maximum_fractional_radius_amplitude
                .max(response.fractional_radius_amplitude);
            diagnostics.maximum_wall_mach_amplitude = diagnostics
                .maximum_wall_mach_amplitude
                .max(response.wall_mach_amplitude);
            diagnostics.maximum_liquid_pressure_ratio = diagnostics
                .maximum_liquid_pressure_ratio
                .max(response.liquid_pressure_ratio);
        }
        diagnostics.void_fraction_per_pa =
            C64::new(compliance_real.value(), compliance_imaginary.value());
        let wavenumber =
            coherent_wavenumber(medium, omega_rad_s, diagnostics.void_fraction_per_pa)?;
        let spacing = det::pow(self.population.number_density_per_m3(), -1.0 / 3.0);
        let collective = collective_admission(
            k_l.max(wavenumber.abs()),
            maximum_radius,
            spacing,
            maximum_scattering_amplitude,
            scattering.value(),
            k_l,
            self.limits,
        )?;
        diagnostics.mean_bubble_spacing_m = Some(spacing);
        diagnostics.wavenumber_radius = collective.wavenumber_radius;
        diagnostics.wavenumber_spacing = collective.wavenumber_spacing;
        diagnostics.neighbor_scattered_pressure_ratio =
            collective.neighbor_scattered_pressure_ratio;
        diagnostics.multiple_scattering_parameter = collective.multiple_scattering_parameter;
        poll(cancelled)?;
        Ok(PolydisperseBubblyLiquidResponse {
            medium,
            omega_rad_s,
            wavenumber,
            diagnostics,
        })
    }
}

/// Aggregate physical response of the original distribution, with per-class
/// amplitude envelopes represented by their maxima, never a mean-radius proxy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PopulationResponseDiagnostics {
    pub void_fraction: f64,
    pub number_density_per_m3: f64,
    pub occupied_classes: usize,
    pub maximum_local_pressure_amplitude_pa: f64,
    /// Complex delta(alpha)/delta(p) [1/Pa] under exp(-i omega t).
    pub void_fraction_per_pa: C64,
    pub maximum_fractional_radius_amplitude: f64,
    pub maximum_wall_mach_amplitude: f64,
    pub maximum_liquid_pressure_ratio: f64,
    pub mean_bubble_spacing_m: Option<f64>,
    /// Largest max(k_l, |k|) R_j among occupied classes.
    pub wavenumber_radius: f64,
    pub wavenumber_spacing: f64,
    pub neighbor_scattered_pressure_ratio: f64,
    /// sum_j n_j 4 pi |f_j|^2/k_l.
    pub multiple_scattering_parameter: f64,
}

/// Immutable admitted harmonic propagation through the original duct owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolydisperseBubblyLiquidResponse {
    medium: FluidAcousticState,
    omega_rad_s: f64,
    wavenumber: C64,
    diagnostics: PopulationResponseDiagnostics,
}
impl PolydisperseBubblyLiquidResponse {
    pub const fn liquid(&self) -> FluidAcousticState {
        self.medium
    }
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }
    pub const fn wavenumber_per_m(&self) -> C64 {
        self.wavenumber
    }
    pub const fn diagnostics(&self) -> PopulationResponseDiagnostics {
        self.diagnostics
    }

    pub fn segment_wave(&self, duct_radius_m: f64) -> Result<SegmentWave, BubblyLiquidError> {
        Ok(fluid::segment_wave_with_wavenumber(
            &self.medium,
            duct_radius_m,
            self.omega_rad_s,
            self.wavenumber,
        )?)
    }

    /// The declared local pressure bound applies to the resulting complete
    /// field. Coherent scattering attenuation is not deposited as latent heat.
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

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), BubblyLiquidError> {
    if cancelled() {
        Err(CavitationError::Cancelled.into())
    } else {
        Ok(())
    }
}
#[derive(Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }
    fn value(&self) -> f64 {
        self.sum + self.correction
    }
}

#[cfg(test)]
#[path = "population/tests.rs"]
mod tests;
