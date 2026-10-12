//! Complete modal loading of a planar radiator in an infinite rigid baffle.
//!
//! This projects the same half-space Rayleigh kernel and equivalent-disc self
//! integral as `baffled_piston_impedance`. It does not solve a closed-body BEM
//! problem. Positive area cells carry signed, work-conjugate velocity shapes;
//! every retained field contributes to the full reciprocal impedance matrix.
//! Fields use `exp(-i omega t)`: a mass-like load has negative imaginary part.
use super::{HelmholtzError, Medium, MAX_DENSE_PANELS, MIN_PANELS_PER_WAVELENGTH, green};
use fs_math::{c64::C64, det};
use core::f64::consts::PI;

/// Complete structural velocity fields retained by the loading integral.
/// Downstream passive realizations may project these fields into fewer acoustic
/// ports; that reduction must not remove coordinates from this physical matrix.
pub const MAX_MODAL_RAYLEIGH_FIELDS: usize = 128;

/// Bound on `cells^2 * modes + cells * modes^2` per frequency. This admits every
/// previous 32-field/8192-cell case while bounding the cost of larger modal
/// spaces. A caller preparing many frequencies must additionally bound its
/// complete sweep. No cells or fields are silently discarded to meet this cap.
pub const MAX_MODAL_RAYLEIGH_WORK: u64 = 2_200_000_000;

/// One positive-area quadrature cell on the radiator's plane, in SI units.
#[derive(Clone, Debug)]
pub struct AreaSample {
    /// Coordinates in the common baffle plane [m].
    pub position_m: [f64; 2],
    /// Physical integration weight [m^2], included once in each work integral.
    pub area_m2: f64,
    /// Normal velocity per generalized modal velocity; all cells cover all modes.
    pub mode_shape: Vec<f64>,
}

/// Geometry prepared once for a bounded frequency sweep; no spatial truncation.
pub struct PreparedSurface {
    count: usize,
    modes: usize,
    areas: Vec<f64>,
    radii: Vec<f64>,
    weighted_shape: Vec<f64>,
    /// Strict upper triangle, in row order. At 8192 cells this is under 256 MiB.
    distances: Vec<f64>,
    maximum_radius_m: f64,
}

fn invalid(what: &'static str) -> HelmholtzError {
    HelmholtzError::BadParameter { what }
}

fn reserved<T>(capacity: usize, what: &'static str) -> Result<Vec<T>, HelmholtzError> {
    let mut values = Vec::new();
    values.try_reserve_exact(capacity)
        .map_err(|_| HelmholtzError::AllocationFailed { what })?;
    Ok(values)
}

impl PreparedSurface {
    /// Retain every supplied cell and mode. Coincident cells refuse: their
    /// distinct physical self integrals cannot be inferred from area alone.
    /// The complete geometry and one impedance evaluation use less than
    /// 281 MiB at the independent cell/field caps; every reservation is fallible.
    /// The modal work cap is checked before any geometry storage is allocated.
    pub fn new(samples: &[AreaSample]) -> Result<Self, HelmholtzError> {
        let count = samples.len();
        if count == 0 || count > MAX_DENSE_PANELS {
            return Err(HelmholtzError::WorkCap { panels: count });
        }
        let modes = samples[0].mode_shape.len();
        if !(1..=MAX_MODAL_RAYLEIGH_FIELDS).contains(&modes) {
            return Err(invalid("Rayleigh loading requires 1..128 complete modal fields"));
        }
        let cells = count as u64;
        let fields = modes as u64;
        if cells * fields * (cells + fields) > MAX_MODAL_RAYLEIGH_WORK {
            return Err(invalid("Rayleigh loading exceeds the per-frequency modal work limit"));
        }
        let mut areas = reserved(count, "complete Rayleigh cell areas")?;
        let mut radii = reserved(count, "complete Rayleigh cell radii")?;
        let mut weighted_shape = reserved(count * modes, "complete Rayleigh weighted shapes")?;
        let mut maximum_radius_m = 0.0_f64;
        for cell in samples {
            if !cell.area_m2.is_finite() || cell.area_m2 <= 0.0
                || cell.position_m.iter().any(|v| !v.is_finite())
                || cell.mode_shape.len() != modes
                || cell.mode_shape.iter().any(|v| !v.is_finite()) {
                return Err(invalid("Rayleigh cells need finite coordinates, positive areas and complete finite shapes"));
            }
            let radius = det::sqrt(cell.area_m2 / PI);
            if !radius.is_finite() || radius <= 0.0 {
                return Err(invalid("Rayleigh cell radius exceeds finite numerical range"));
            }
            areas.push(cell.area_m2);
            radii.push(radius);
            maximum_radius_m = maximum_radius_m.max(radius);
            for &shape in &cell.mode_shape {
                let weighted = cell.area_m2 * shape;
                if !weighted.is_finite() || (shape != 0.0 && weighted == 0.0) {
                    return Err(invalid("Rayleigh area-weighted shape exceeds finite numerical range"));
                }
                weighted_shape.push(weighted);
            }
        }
        let pairs = count * (count - 1) / 2;
        let mut distances = reserved(pairs, "complete Rayleigh pair distances")?;
        for i in 0..count {
            for j in i + 1..count {
                let a = samples[i].position_m;
                let b = samples[j].position_m;
                let r = (a[0] - b[0]).hypot(a[1] - b[1]);
                if !r.is_finite() || r <= 0.0 {
                    return Err(invalid("Rayleigh area cells must have distinct, finitely separated centres"));
                }
                distances.push(r);
            }
        }
        Ok(Self { count, modes, areas, radii, weighted_shape, distances, maximum_radius_m })
    }

    /// Wavelength divided by the largest equivalent-disc diameter. This is an
    /// area-resolution diagnostic, not a bound on shape, aspect-ratio or
    /// quadrature error. The same six-cells-per-wavelength floor as the BEM
    /// owner is enforced before a sampled impedance is formed.
    pub fn equivalent_disc_ppw(&self, omega: f64, medium: Medium) -> Result<f64, HelmholtzError> {
        if !omega.is_finite() || omega <= 0.0 || !medium.density.is_finite()
            || medium.density <= 0.0 || !medium.sound_speed.is_finite() || medium.sound_speed <= 0.0 {
            return Err(invalid("Rayleigh frequency, density and sound speed must be positive and finite"));
        }
        let k = omega / medium.sound_speed;
        let diameter_in_wavelengths = k * self.maximum_radius_m / PI;
        let ppw = 1.0 / diameter_in_wavelengths;
        if !k.is_finite() || k <= 0.0 || !ppw.is_finite() || ppw <= 0.0 {
            return Err(invalid("Rayleigh wavelength exceeds finite numerical range"));
        }
        Ok(ppw)
    }

    /// Generalized force per generalized velocity for ONE radiating half-space.
    /// This is `Phi^T A K A Phi`, with the Rayleigh `2G` kernel off the diagonal
    /// and the existing equivalent-disc source integral on each diagonal.
    /// Both off-diagonal directions share one kernel evaluation. Projecting
    /// the pressure fields costs O(cells^2*modes + cells*modes^2), not one
    /// spatial double integral per matrix entry. Each reciprocal entry is
    /// evaluated once; data are never averaged or clipped to appear passive.
    pub fn impedance(&self, omega: f64, medium: Medium) -> Result<Vec<C64>, HelmholtzError> {
        let ppw = self.equivalent_disc_ppw(omega, medium)?;
        if ppw < MIN_PANELS_PER_WAVELENGTH {
            return Err(HelmholtzError::TooCoarse { panels_per_wavelength: ppw });
        }
        let k = omega / medium.sound_speed;
        let prefactor = C64::new(0.0, -2.0 * omega * medium.density);
        if !prefactor.im.is_finite() {
            return Err(invalid("Rayleigh pressure scale exceeds finite numerical range"));
        }
        let mut pressure = reserved(self.count * self.modes, "complete Rayleigh pressure fields")?;
        pressure.resize(self.count * self.modes, C64::ZERO);
        for i in 0..self.count {
            let phase = k * self.radii[i];
            // rho*c*(1-exp(i*k*a))/A is the existing disc integral after
            // multiplying by -i*omega*rho/(2*pi) and dividing by source area.
            // 1-cos(x)=2*sin(x/2)^2 retains weak low-frequency resistance.
            let scale = medium.density * medium.sound_speed / self.areas[i];
            let self_kernel = C64::new(2.0 * det::sin(0.5 * phase).powi(2), -det::sin(phase)).scale(scale);
            for mode in 0..self.modes {
                pressure[i * self.modes + mode] = self_kernel.scale(self.weighted_shape[i * self.modes + mode]);
            }
        }
        let mut pair = 0;
        for i in 0..self.count {
            for j in i + 1..self.count {
                let r = self.distances[pair];
                pair += 1;
                if !(k * r).is_finite() {
                    return Err(invalid("Rayleigh travel phase exceeds finite numerical range"));
                }
                let kernel = green(k, r) * prefactor;
                for mode in 0..self.modes {
                    let a = i * self.modes + mode;
                    let b = j * self.modes + mode;
                    pressure[a] = pressure[a] + kernel.scale(self.weighted_shape[b]);
                    pressure[b] = pressure[b] + kernel.scale(self.weighted_shape[a]);
                }
            }
        }
        let mut result = reserved(self.modes * self.modes, "complete Rayleigh modal impedance")?;
        result.resize(self.modes * self.modes, C64::ZERO);
        for row in 0..self.modes {
            for column in row..self.modes {
                let mut value = C64::ZERO;
                for cell in 0..self.count {
                    value = value + pressure[cell * self.modes + column].scale(self.weighted_shape[cell * self.modes + row]);
                }
                if !value.re.is_finite() || !value.im.is_finite() {
                    return Err(invalid("Rayleigh modal impedance exceeds finite numerical range"));
                }
                result[row * self.modes + column] = value;
                result[column * self.modes + row] = value;
            }
        }
        Ok(result)
    }
}

#[cfg(test)]
#[path = "rayleigh_modal_tests.rs"]
mod tests;

