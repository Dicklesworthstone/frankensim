//! Source-specific stress birefringence of Corning 7980 synthetic fused silica.
//!
//! Richard Priestley (Corning), "Birefringence dispersion in fused silica for
//! DUV lithography", Proc. SPIE 4346, 1300--1305 (2001), equation (2), Table 1,
//! <https://doi.org/10.1117/12.435666>. The three numerical coefficient facts
//! below are transcribed with attribution; no article, figure or fit is copied.
//!
//! The experiment used uniaxial compression. Table 1 reports specimen averages
//! plus/minus one standard deviation, not confidence bounds or a per-ray random
//! law. Temperature was not specified. The load-cell capacity was 1.24 MPa;
//! wavelength-specific plotted fit ranges differ. That instrument capacity is
//! NOT a certified multi-axial stress envelope. No temperature law or continuous
//! spectral fit is inferred. Only the reported 193, 248 and 633 nm lines exist.
//!
//! We use tension-positive Cauchy stress in the existing tensor-shear Voigt
//! order. The general transverse-tensor calculation is an explicitly adopted
//! isotropic, linear stress-optic constitutive extension of the measured law
//! n_x - n_y = R (sigma_x - sigma_y). It is not a new multi-axial calibration.
//! Callers own applicability to their stress, temperature and specimen state.
//! Only differential indices and retardation are supplied: the common index
//! shift, hydrostatic response, mode walk-off and absorption remain unknown.

use crate::Voigt;
use fs_blake3::{ContentHash, DomainHasher};
use fs_math::det;

pub const SOURCE_RELEASE: &str = "Priestley 2001, Corning 7980, equation (2), Table 1";
pub const SOURCE_URL: &str = "https://doi.org/10.1117/12.435666";
pub const SOURCE_PDF_URL: &str =
    "https://www.phys.hawaii.edu/~idlab/taskAndSchedule/iTOP/PSI001300.pdf";
/// Nominal source laser lines; these points are not interpolation endpoints.
pub const SOURCE_WAVELENGTHS_NM: [f64; 3] = [193.0, 248.0, 633.0];

// (nominal wavelength nm, mean R, specimen standard deviation), nm/cm/MPa.
const SOURCE_COEFFICIENTS: [(f64, f64, f64); 3] =
    [(193.0, 51.0, 3.0), (248.0, 42.0, 1.0), (633.0, 35.0, 1.0)];
const NM_PER_CM_PER_MPA_TO_PER_PA: f64 = 1.0e-13;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SilicaStressOpticError {
    UnsupportedSourceWavelength,
    InvalidStress,
    InvalidTransverseFrame,
    InvalidDistance,
    NumericalRange,
}

impl core::fmt::Display for SilicaStressOpticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Corning 7980 reference stress optics: {self:?}")
    }
}
impl core::error::Error for SilicaStressOpticError {}

/// One original source line, without an unstated temperature or spectral model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Corning7980StressOptics {
    wavelength_nm: f64,
    coefficient_per_pa: f64,
    specimen_stddev_per_pa: f64,
}

impl Corning7980StressOptics {
    pub fn at_source_wavelength(wavelength_nm: f64) -> Result<Self, SilicaStressOpticError> {
        let &(_, mean, stddev) = SOURCE_COEFFICIENTS
            .iter()
            .find(|&&(wavelength, _, _)| wavelength == wavelength_nm)
            .ok_or(SilicaStressOpticError::UnsupportedSourceWavelength)?;
        Ok(Self {
            wavelength_nm,
            coefficient_per_pa: mean * NM_PER_CM_PER_MPA_TO_PER_PA,
            specimen_stddev_per_pa: stddev * NM_PER_CM_PER_MPA_TO_PER_PA,
        })
    }

    pub const fn wavelength_nm(self) -> f64 {
        self.wavelength_nm
    }
    /// Positive R in n_x - n_y = R (sigma_x - sigma_y), in inverse Pa.
    pub const fn coefficient_per_pa(self) -> f64 {
        self.coefficient_per_pa
    }
    /// Source's between-specimen one-sigma spread; no confidence level inferred.
    pub const fn specimen_standard_deviation_per_pa(self) -> f64 {
        self.specimen_stddev_per_pa
    }

    pub fn source_identity() -> ContentHash {
        let mut hash = DomainHasher::new("fs-material:corning7980-stress-optics:v1");
        hash.update(SOURCE_RELEASE.as_bytes());
        hash.update(SOURCE_URL.as_bytes());
        hash.update(b"nominal_nm;tension_positive_Pa;nm/cm/MPa;specimen_mean_plus_minus_1sigma;temperature_unspecified;isotropic_linear_extension");
        hash.update(&NM_PER_CM_PER_MPA_TO_PER_PA.to_le_bytes());
        for &(wavelength, mean, stddev) in &SOURCE_COEFFICIENTS {
            for value in [wavelength, mean, stddev] {
                hash.update(&value.to_le_bytes());
            }
        }
        hash.finalize()
    }

    /// Project the ORIGINAL Cauchy tensor onto an orthonormal transverse frame.
    ///
    /// Borrowed Voigt stress uses [xx, yy, zz, xy, yz, zx] in Pa, with tensor
    /// shear, not engineering shear. The reported pair is R*(sigma_hh-sigma_vv)
    /// and 2*R*sigma_hv. It represents only the trace-free index difference.
    /// A common pressure cancels without creating artificial birefringence.
    pub fn project(
        self,
        stress_pa: &Voigt,
        horizontal: [f64; 3],
        vertical: [f64; 3],
    ) -> Result<TransverseBirefringence, SilicaStressOpticError> {
        if !stress_pa.iter().all(|v| v.is_finite()) {
            return Err(SilicaStressOpticError::InvalidStress);
        }
        if !horizontal.into_iter().chain(vertical).all(f64::is_finite)
            || (dot(horizontal, horizontal) - 1.0).abs() > 1.0e-12
            || (dot(vertical, vertical) - 1.0).abs() > 1.0e-12
            || dot(horizontal, vertical).abs() > 1.0e-12
        {
            return Err(SilicaStressOpticError::InvalidTransverseFrame);
        }
        // Remove a common pressure before projection. This makes a hydrostatic
        // tensor exactly null even for rounded oblique orthonormal frames.
        let [xx, yy, zz, xy, yz, zx] = *stress_pa;
        let tensor = [[xx - zz, xy, zx], [xy, yy - zz, yz], [zx, yz, 0.0]];
        let h = tensor.map(|row| dot(row, horizontal));
        let v = tensor.map(|row| dot(row, vertical));
        let q = self.coefficient_per_pa * (dot(horizontal, h) - dot(vertical, v));
        let u = self.coefficient_per_pa * (2.0 * dot(horizontal, v));
        if !q.is_finite() || !u.is_finite() {
            return Err(SilicaStressOpticError::NumericalRange);
        }
        Ok(TransverseBirefringence {
            horizontal_minus_vertical: q,
            twice_cross: u,
            wavelength_nm: self.wavelength_nm,
        })
    }
}

/// Trace-free differential index in one transverse frame, not two absolute n's.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransverseBirefringence {
    horizontal_minus_vertical: f64,
    twice_cross: f64,
    wavelength_nm: f64,
}

impl TransverseBirefringence {
    pub const fn horizontal_minus_vertical(self) -> f64 {
        self.horizontal_minus_vertical
    }
    pub const fn twice_cross(self) -> f64 {
        self.twice_cross
    }

    /// Retardance and twice-angle slow-axis components for a metric path.
    /// The higher-n (slow) principal axis has (cos(2*theta), sin(2*theta)).
    pub fn retardance(self, distance_m: f64) -> Result<LinearRetardance, SilicaStressOpticError> {
        if !distance_m.is_finite() || distance_m < 0.0 {
            return Err(SilicaStressOpticError::InvalidDistance);
        }
        let q = self.horizontal_minus_vertical;
        let u = self.twice_cross;
        let scale = q.abs().max(u.abs());
        if scale == 0.0 {
            return Ok(LinearRetardance {
                phase_rad: 0.0,
                cos_double_slow_axis: 1.0,
                sin_double_slow_axis: 0.0,
            });
        }
        let gap = scale * det::sqrt((q / scale).powi(2) + (u / scale).powi(2));
        let phase = core::f64::consts::TAU * (gap * distance_m / (self.wavelength_nm * 1.0e-9));
        if !phase.is_finite() {
            return Err(SilicaStressOpticError::NumericalRange);
        }
        Ok(LinearRetardance {
            phase_rad: phase,
            cos_double_slow_axis: q / gap,
            sin_double_slow_axis: u / gap,
        })
    }
}

/// Nonnegative slow-minus-fast phase, retaining the signed principal direction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearRetardance {
    pub phase_rad: f64,
    pub cos_double_slow_axis: f64,
    pub sin_double_slow_axis: f64,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[cfg(test)]
mod tests;
