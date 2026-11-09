//! Coherent thin-body material observations through the original 3D memories.

use super::super::super::{
    TensorPronyHarmonicSpectrum, TensorPronyHarmonicWindowPoint, TensorPronySpectrumWindow,
};
use super::super::projection::{
    Sum,
    window::{MAX_SPECTRUM_LINES, StrainParts},
};
use super::{MaterialPlane, PlaneStressHarmonic, TensorProny, ViscoError};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_math::{c64::C64, dd::Dd};

/// Settled plane-stress lines with their original frequency-dependent transverse strains.
///
/// The engineering coordinates are `[epsilon_x, epsilon_y, gamma_xy]` in the
/// ordered material plane. Each distinct frequency eliminates the complete
/// complex tensor, including Kelvin viscosity and every Maxwell branch. The
/// resulting six-component strains feed the original 3D spring, dashpot and
/// stress-work contractions; no reduced Prony branches or constant Poisson
/// ratio are substituted. Peak phasors share the `exp(-i omega t)` clock.
pub struct PlaneStressHarmonicSpectrum<'a> {
    source: TensorPronyHarmonicSpectrum<'a>,
    responses: Vec<PlaneStressHarmonic>,
    // Cache the original DC constraint or its refusal without requiring DC
    // applicability from callers observing an exactly zero static strain.
    relaxed_response: Result<PlaneStressHarmonic, ViscoError>,
    plane: MaterialPlane,
}

/// Exact finite-interval material work and heat for coherent plane-stress lines.
///
/// Temporal carrier and beat moments are prepared once. Point evaluation uses
/// fixed stack scratch with no allocation or constitutive solve, and retains
/// the same O(F²*(branches+1)*6) bound as the full-tensor spectrum observer.
pub struct PlaneStressSpectrumWindow<'p, 'a> {
    source: &'p PlaneStressHarmonicSpectrum<'a>,
    response: TensorPronySpectrumWindow<'p, 'a>,
}

impl TensorProny {
    /// Prepare coherent plane-stress observation for up to 64 authored lines.
    ///
    /// Frequencies are positive, finite, source-admitted and need no common
    /// period. Exact duplicates combine before transverse recovery and strain
    /// admission. Empty input is the zero AC field. DC/prestrain and causal
    /// startup are not inferred from settled harmonic strains.
    ///
    /// Retained numeric payload is at most size_of::<PlaneStressHarmonicSpectrum>()
    /// + 128*F*(branches+1) + F*size_of::<PlaneStressHarmonic>() bytes, excluding
    /// the borrowed original material and allocator metadata.
    pub fn prepare_plane_stress_harmonic_spectrum(
        &self,
        plane: MaterialPlane,
        angular_frequencies_rad_s: &[f64],
    ) -> Result<PlaneStressHarmonicSpectrum<'_>, ViscoError> {
        let source = self.prepare_harmonic_spectrum(angular_frequencies_rad_s)?;
        let mut responses = Vec::with_capacity(source.distinct_frequencies().len());
        for &omega in source.distinct_frequencies() {
            responses.push(self.prepare_plane_stress_harmonic(plane, omega)?);
        }
        Ok(PlaneStressHarmonicSpectrum {
            source,
            responses,
            relaxed_response: self.prepare_plane_stress_harmonic(plane, 0.0),
            plane,
        })
    }
}

impl<'a> PlaneStressHarmonicSpectrum<'a> {
    /// Original 3D constitutive source, including its frequency and strain domain.
    pub const fn law(&self) -> &'a TensorProny {
        self.source.law()
    }

    pub const fn plane(&self) -> MaterialPlane {
        self.plane
    }

    /// Original caller order, including repeated frequencies [rad/s].
    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        self.source.angular_frequencies_rad_s()
    }

    /// Retained numeric payload, excluding the separately borrowed source law.
    pub fn numerical_payload_bytes(&self) -> usize {
        core::mem::size_of::<Self>() + self.source.numerical_payload_bytes()
            - core::mem::size_of::<TensorPronyHarmonicSpectrum<'_>>()
            + self.responses.capacity() * core::mem::size_of::<PlaneStressHarmonic>()
    }

    /// Admit the complete recovered 3D waveform against the original strain bound.
    ///
    /// Distinct-frequency ellipse peak norms are summed conservatively for all
    /// time, including transverse normal and shear strains. Equal-frequency
    /// interference is resolved first, using retained compensated coordinates.
    /// This can refuse a valid waveform whose distinct peaks cannot coincide;
    /// it is not a sampled-window or exact multitone supremum certificate.
    pub fn validate_strains(&self, in_plane: &[[C64; 3]]) -> Result<(), ViscoError> {
        let strains = self.recover_grouped_strains(in_plane)?;
        self.source
            .admit_grouped_strains(&strains[..self.responses.len()])
    }

    /// Admit a relaxed in-plane DC bias together with the complete AC waveform.
    /// The bias uses the original zero-frequency transverse constraint, while
    /// each AC line retains its own dispersive constraint. Nonzero bias requires
    /// original DC applicability; an exactly zero bias leaves AC admission intact.
    pub fn validate_strains_with_bias(
        &self,
        in_plane: &[[C64; 3]],
        bias: &[f64; 3],
    ) -> Result<(), ViscoError> {
        let strains = self.recover_grouped_strains(in_plane)?;
        self.source.admit_grouped_strains_with_bias(
            &strains[..self.responses.len()],
            &self.recover_bias(bias)?,
        )
    }

    fn recover_bias(&self, bias: &[f64; 3]) -> Result<[Dd; 6], ViscoError> {
        if bias.iter().any(|value| !value.is_finite()) {
            return Err(ViscoError::Parameters {
                what: "relaxed plane-stress bias must be a finite engineering strain",
            });
        }
        if bias.iter().all(|&value| value == 0.0) {
            return Ok([Dd::ZERO; 6]);
        }
        let response = self.relaxed_response.as_ref().map_err(Clone::clone)?;
        Ok(response.strain_map.map(|row| {
            let mut value = Sum::default();
            for (map, &strain) in row.iter().zip(bias) {
                value.add([map.re, strain, 1.0, 1.0]);
            }
            value.retained_value()
        }))
    }

    fn recover_grouped_strains(
        &self,
        in_plane: &[[C64; 3]],
    ) -> Result<[StrainParts; MAX_SPECTRUM_LINES], ViscoError> {
        let grouped = self.source.group_strains(in_plane)?;
        let mut strains = [[[Dd::ZERO; 2]; 6]; MAX_SPECTRUM_LINES];
        for ((result, coordinates), response) in
            strains.iter_mut().zip(&grouped).zip(&self.responses)
        {
            for (value, row) in result.iter_mut().zip(&response.strain_map) {
                let (mut re, mut im) = (Sum::default(), Sum::default());
                for (map, z) in row.iter().zip(coordinates) {
                    // Restore physical units only after the full map contracts.
                    // Both words survive equal-tone cancellation and recovery.
                    for (real, imaginary) in [(z[0].hi, z[1].hi), (z[0].lo, z[1].lo)] {
                        re.add([map.re, real, 1.0, 1.0]);
                        re.add([-map.im, imaginary, 1.0, 1.0]);
                        im.add([map.re, imaginary, 1.0, 1.0]);
                        im.add([map.im, real, 1.0, 1.0]);
                    }
                }
                *value = [re.retained_value(), im.retained_value()];
            }
        }
        Ok(strains)
    }

    /// Cache a uniform physical interval once for repeated thin-body points.
    /// Additional temporal payload is bounded by size_of::<PlaneStressSpectrumWindow>()
    /// + 256*F + 64*F² bytes. The source cache is borrowed separately.
    pub fn prepare_window(
        &self,
        window: HarmonicWindow,
    ) -> Result<PlaneStressSpectrumWindow<'_, 'a>, ViscoError> {
        Ok(PlaneStressSpectrumWindow {
            source: self,
            response: self.source.prepare_window(window)?,
        })
    }
}

impl<'p, 'a> PlaneStressSpectrumWindow<'p, 'a> {
    pub const fn source(&self) -> &'p PlaneStressHarmonicSpectrum<'a> {
        self.source
    }

    pub const fn window(&self) -> HarmonicWindow {
        self.response.window()
    }

    /// Additional temporal payload, excluding the borrowed material source cache.
    pub fn numerical_payload_bytes(&self) -> usize {
        core::mem::size_of::<Self>() + self.response.numerical_payload_bytes()
            - core::mem::size_of::<TensorPronySpectrumWindow<'_, '_>>()
    }

    /// Original material stored energy, Maxwell/Kelvin heat and independent work.
    ///
    /// Every quantity is per original material volume [J/m³]. Shell/plate
    /// callers supply their actual through-thickness quadrature volumes; these
    /// are neither section resultants nor cycle-averaged heating powers.
    pub fn evaluate(
        &self,
        in_plane: &[[C64; 3]],
    ) -> Result<TensorPronyHarmonicWindowPoint, ViscoError> {
        self.evaluate_with_checkpoint(in_plane, || Ok(()))
    }

    /// Observe the real thin-body waveform with an explicitly relaxed DC bias.
    ///
    /// Only the original equilibrium spring carries the static strain. Its
    /// transverse recovery uses omega=0, not an AC Poisson ratio or branchwise
    /// condensation. DC/AC storage and stress-work cross terms are retained;
    /// the original Maxwell/Kelvin heat is unchanged. Outputs remain energy
    /// densities per material volume. No startup memory, finite-strain or
    /// geometric-stiffness/prestress effects are inferred.
    pub fn evaluate_with_bias(
        &self,
        in_plane: &[[C64; 3]],
        bias: &[f64; 3],
    ) -> Result<TensorPronyHarmonicWindowPoint, ViscoError> {
        self.evaluate_with_bias_and_checkpoint(in_plane, bias, || Ok(()))
    }

    /// Biased observation with the same bounded caller-owned cancellation.
    pub fn evaluate_with_bias_and_checkpoint<E: From<ViscoError>>(
        &self,
        in_plane: &[[C64; 3]],
        bias: &[f64; 3],
        mut checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<TensorPronyHarmonicWindowPoint, E> {
        checkpoint()?;
        let strains = self.source.recover_grouped_strains(in_plane)?;
        self.response.evaluate_grouped_with_bias_and_checkpoint(
            &strains[..self.source.responses.len()],
            &self.source.recover_bias(bias)?,
            checkpoint,
        )
    }

    /// Caller-owned cancellation before recovery and between original branches.
    /// The bounded recovery takes at most 64*6*3 complex two-word contractions.
    /// No accepted material history or partially returned field is mutated.
    pub fn evaluate_with_checkpoint<E: From<ViscoError>>(
        &self,
        in_plane: &[[C64; 3]],
        mut checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<TensorPronyHarmonicWindowPoint, E> {
        checkpoint()?;
        let strains = self.source.recover_grouped_strains(in_plane)?;
        self.response
            .evaluate_grouped_with_checkpoint(&strains[..self.source.responses.len()], checkpoint)
    }
}
