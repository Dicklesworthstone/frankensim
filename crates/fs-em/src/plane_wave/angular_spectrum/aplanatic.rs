//! Vectorial aplanatic focusing on the original physical Fourier grid.
//!
//! The reference-sphere field is sqrt(cos(theta)) times the rotated pupil
//! Jones vector. Cartesian angular quadrature has dOmega=dsx*dsy/cos(theta),
//! so the Fourier coefficient carries 1/sqrt(cos(theta)), not sqrt(cos(theta)).
//! The common quadrature/focal-length factor is fixed by the requested incident
//! power. Every mode keeps its complex pupil phase, including aberrations.
//!
//! Reference: Liu et al., "Revisiting PSF models: unifying framework and
//! high-performance implementation", equations (14), (24):
//! <https://arxiv.org/html/2502.03170v1>. The exp(+i omega t) convention here is
//! conjugate to that reference. The common pupil-to-focus phase is chosen zero;
//! an authored Jones phase remains the coherent phase relative to other beams.

use super::modes::{basis, cartesian, finite, incident_cosine, polarized};
use super::{AngularSpectrumError, AngularSpectrumGrid, bad, poll, transverse_number};
use fs_math::{c64::C64, dd::Dd, det};

/// Ideal lossless aplanatic lens illumination, expressed in the incident medium.
///
/// `sine_aperture` is sin(theta_max)=NA/n, strictly between zero and one. The
/// pupil is centered on +z. Its Jones shape can describe apodization, stops,
/// vortices, aberrations and spatially varying polarization. Incident watts
/// mean the total cycle-mean +z power of the resulting periodic Fourier beam,
/// after the authored pupil shape; they are not power before a lossy stop.
///
/// This is the convergent Cartesian quadrature of the Richards–Wolf focal
/// model, not a solution of a finite material lens or a subwavelength screen.
/// Increase the physical cell to refine angular quadrature, and decrease its
/// sample spacing to resolve finer output detail. Neither change is automatic.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AplanaticPupil {
    sine_aperture: f64,
    focus_m: [f64; 3],
    incident_power_w: f64,
}

impl AplanaticPupil {
    pub fn try_new(
        sine_aperture: f64,
        focus_m: [f64; 3],
        incident_power_w: f64,
    ) -> Result<Self, AngularSpectrumError> {
        if !sine_aperture.is_finite() || sine_aperture <= 0.0 || sine_aperture >= 1.0 {
            return Err(bad(
                "aplanatic sine aperture must lie strictly between zero and one",
            ));
        }
        if !focus_m.iter().all(|x| x.is_finite())
            || !incident_power_w.is_finite()
            || incident_power_w <= 0.0
        {
            return Err(bad(
                "aplanatic focus must be finite and incident watts positive and finite",
            ));
        }
        Ok(Self {
            sine_aperture,
            focus_m,
            incident_power_w,
        })
    }

    pub fn sine_aperture(self) -> f64 {
        self.sine_aperture
    }
    pub fn focus_m(self) -> [f64; 3] {
        self.focus_m
    }
    pub fn incident_power_w(self) -> f64 {
        self.incident_power_w
    }

    // The parent has admitted the grid, its area, and these N complex-pair
    // samples before this function. No source callback runs before admission.
    pub(super) fn sample_modes(
        self,
        grid: AngularSpectrumGrid,
        area_m2: f64,
        k: f64,
        admittance_s: f64,
        mut pupil_jones: impl FnMut(AplanaticPupilSample) -> [C64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[C64; 2]>, AngularSpectrumError> {
        poll(cancelled)?;
        if !k.is_finite() || k <= 0.0 || !admittance_s.is_finite() || admittance_s <= 0.0 {
            return Err(bad(
                "aplanatic incident wave number or admittance is unrepresentable",
            ));
        }
        let radius = k * self.sine_aperture;
        if !radius.is_finite()
            || radius <= 0.0
            || grid
                .spacing_m
                .iter()
                .any(|&dx| radius >= core::f64::consts::PI / dx)
        {
            return Err(bad(
                "the complete aplanatic pupil must lie strictly inside both grid Nyquist limits",
            ));
        }
        let n = grid.width * grid.height;
        let mut modes = Vec::new();
        modes
            .try_reserve_exact(n)
            .map_err(|_| bad("aplanatic source allocation failed"))?;
        let mut largest = 0.0_f64;
        for row in 0..grid.height {
            for column in 0..grid.width {
                if (row * grid.width + column) % 32 == 0 {
                    poll(cancelled)?;
                }
                let value = match self.sample(grid, k, radius, column, row)? {
                    Some(sample) => pupil_jones(sample),
                    None => [C64::ZERO; 2],
                };
                if !value.iter().all(|&z| finite(z)) {
                    return Err(bad("aplanatic Jones samples must be finite"));
                }
                for z in value {
                    largest = largest.max(z.re.abs()).max(z.im.abs());
                }
                modes.push(value);
            }
        }
        if largest == 0.0 {
            return Err(bad(
                "a positive-power aplanatic pupil must transmit a nonzero sample",
            ));
        }

        // Normalize the authored shape before squaring or rotating, so even
        // subnormal/very large Jones magnitudes describe the same finite shape.
        // Rotation plus 1/sqrt(cos) leaves cos*|E|^2 = |Jones|^2 mode by mode.
        let mut shape_power = Dd::ZERO;
        for (index, value) in modes.iter_mut().enumerate() {
            if index % 32 == 0 {
                poll(cancelled)?;
            }
            if *value == [C64::ZERO; 2] {
                continue;
            }
            let jones = value.map(|z| C64::new(z.re / largest, z.im / largest));
            for z in jones {
                let re = Dd::from_f64(z.re);
                let im = Dd::from_f64(z.im);
                shape_power = shape_power + re * re + im * im;
            }
            let column = index % grid.width;
            let row = index / grid.width;
            let kx = transverse_number(column, grid.width, grid.spacing_m[0]);
            let ky = transverse_number(row, grid.height, grid.spacing_m[1]);
            let cosine = incident_cosine(k, kx, ky)
                .ok_or_else(|| bad("aplanatic mode left the incident light cone"))?;
            let p = basis(kx, ky);
            let a = polarized(jones, p);
            let root = det::sqrt(cosine);
            let transverse = cartesian([a[0].scale(1.0 / root), a[1].scale(root)], p);
            let phase = Dd::from_f64(kx) * Dd::from_f64(self.focus_m[0])
                + Dd::from_f64(ky) * Dd::from_f64(self.focus_m[1])
                + Dd::from_f64(k * cosine) * Dd::from_f64(self.focus_m[2]);
            if !phase.to_f64().is_finite() {
                return Err(bad("aplanatic focus phase is unrepresentable"));
            }
            let (sin, cos) = det::sin_cos_dd(phase);
            let phasor = C64::new(cos.to_f64(), sin.to_f64());
            *value = transverse.map(|z| z * phasor);
        }
        // Separate square roots avoid overflowing P/area, admittance*area, or
        // a squared field solely while computing a representable amplitude.
        let amplitude = (det::sqrt(self.incident_power_w) / det::sqrt(area_m2))
            / det::sqrt(admittance_s)
            * core::f64::consts::SQRT_2
            / det::sqrt(shape_power.to_f64());
        if !amplitude.is_finite() || amplitude <= 0.0 {
            return Err(bad("aplanatic power normalization is unrepresentable"));
        }
        for (index, value) in modes.iter_mut().enumerate() {
            if index % 256 == 0 {
                poll(cancelled)?;
            }
            *value = value.map(|z| z.scale(amplitude));
            if !value.iter().all(|&z| finite(z)) {
                return Err(bad("aplanatic focused amplitudes are unrepresentable"));
            }
        }
        poll(cancelled)?;
        Ok(modes)
    }

    fn sample(
        self,
        grid: AngularSpectrumGrid,
        k: f64,
        radius: f64,
        column: usize,
        row: usize,
    ) -> Result<Option<AplanaticPupilSample>, AngularSpectrumError> {
        let kx = transverse_number(column, grid.width, grid.spacing_m[0]);
        let ky = transverse_number(row, grid.height, grid.spacing_m[1]);
        if !kx.is_finite() || !ky.is_finite() {
            return Err(bad("aplanatic transverse wave number is unrepresentable"));
        }
        let x = kx / radius;
        let y = ky / radius;
        if x.abs() > 1.0
            || y.abs() > 1.0
            || (Dd::from_f64(x) * Dd::from_f64(x) + Dd::from_f64(y) * Dd::from_f64(y) - Dd::ONE)
                .to_f64()
                > 0.0
        {
            return Ok(None);
        }
        let cosine = incident_cosine(k, kx, ky)
            .ok_or_else(|| bad("aplanatic angular support reaches grazing incidence"))?;
        Ok(Some(AplanaticPupilSample {
            column,
            row,
            pupil_coordinate: [x, y],
            direction: [kx / k, ky / k, cosine],
        }))
    }
}

/// One Cartesian angular-quadrature point in the lens pupil.
///
/// `pupil_coordinate` is the transverse direction divided by the sine aperture,
/// so the circular pupil has radius one. `direction` is the outgoing unit ray
/// in the incident medium. The Jones callback returns its transverse input
/// polarization in the original laboratory x/y basis, before ray rotation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AplanaticPupilSample {
    pub column: usize,
    pub row: usize,
    pub pupil_coordinate: [f64; 2],
    pub direction: [f64; 3],
}

#[cfg(test)]
mod tests;
