//! SI Newtonian no-slip, shear-free falling film on a locally planar wall.
//! q = g_t h^3/(3 nu), u(y) = g_t (h y - y^2/2)/nu, Re = q/nu (mean speed).
//! For exp(s t + i(k x - omega t)), the linear Benney approximation is
//! s = q k^2 [(6/5) Re - g_n/g_t - sigma k^2/(rho g_t)], omega = 3 (q/h) k.
//! The normal acceleration points into the wall; neutral Re = (5/6) g_n/g_t.
//! This near-onset long-wave surface mode is not turbulent transition or a
//! free-surface Orr-Sommerfeld solve. Wavy films can remain laminar. Sampled
//! decay proves neither global stability nor unsampled wavelengths. Side walls,
//! curvature, gas shear, non-Newtonian response and contact lines are omitted.
//! Source: Scheid et al., Phys. Fluids 14 (2002), DOI 10.1063/1.1515270,
//! equation (7), zero thermal terms; their R = 3 Re. Author manuscript:
//! <https://www.uni-muenster.de/Physik.TP/pauli/~thiele/Paper/SOCT2002pf.pdf>.
pub mod material;

/// Prescribed uniform Newtonian liquid properties; no material calibration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilmMaterial {
    /// Mass density in kg/m^3.
    pub density_kg_m3: f64,
    /// Kinematic viscosity in m^2/s.
    pub kinematic_viscosity_m2_s: f64,
    /// Surface tension in N/m; zero explicitly omits capillarity.
    pub surface_tension_n_m: f64,
}
/// Caller-selected asymptotic budgets, each strictly between zero and one.
/// These limit kh and Re*kh; they are not calibrated error bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LongWaveLimits {
    /// Maximum wavenumber times film thickness.
    pub max_kh: f64,
    /// Maximum Reynolds number times wavenumber times thickness.
    pub max_re_kh: f64,
}
/// One admitted infinitesimal surface mode; positive growth means amplification.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilmMode {
    /// Temporal amplitude growth rate in 1/s.
    pub growth_per_s: f64,
    /// Downstream angular frequency in rad/s.
    pub angular_frequency_rad_s: f64,
    /// Wavenumber times film thickness.
    pub kh: f64,
    /// Reynolds number times kh.
    pub re_kh: f64,
}
/// A physical or numerical admission failure; no clipped substitute is returned.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FilmError {
    /// A named SI input is outside the documented domain.
    InvalidInput(&'static str),
    /// Zero flow supplies no finite-thickness film to analyze.
    DryFilm,
    /// A required nonzero derived value underflowed, or a value overflowed.
    NonRepresentable,
    /// The authored mode exceeds at least one long-wave budget.
    OutsideLongWaveRegime,
}

impl core::fmt::Display for FilmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput(name) => write!(f, "invalid falling-film input: {name}"),
            Self::DryFilm => f.write_str("zero flow has no finite-thickness film"),
            Self::NonRepresentable => f.write_str("falling-film state is not representable"),
            Self::OutsideLongWaveRegime => f.write_str("mode exceeds kh or Re*kh budget"),
        }
    }
}
impl std::error::Error for FilmError {}

/// Admitted Nusselt base flow. Construction is allocation-free; profile and
/// mode queries use fixed arithmetic and never refit a material or geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NusseltFilm {
    material: FilmMaterial,
    flow: f64,
    tangent: f64,
    normal_ratio: f64,
    thickness: f64,
    mean_velocity: f64,
    reynolds: f64,
    wall_shear: f64,
}

impl NusseltFilm {
    /// Positive flow per unit width, positive tangential acceleration, and
    /// nonnegative wall-directed normal acceleration define the local film.
    /// Dry flow refuses explicitly. The free surface has zero shear stress.
    pub fn new(
        material: FilmMaterial,
        flow_per_width_m2_s: f64,
        tangential_acceleration_m_s2: f64,
        normal_acceleration_m_s2: f64,
    ) -> Result<Self, FilmError> {
        let q = flow_per_width_m2_s;
        let gt = tangential_acceleration_m_s2;
        let gn = normal_acceleration_m_s2;
        let rho = material.density_kg_m3;
        let nu = material.kinematic_viscosity_m2_s;
        let sigma = material.surface_tension_n_m;
        for (name, value, zero_allowed) in [
            ("density", rho, false),
            ("kinematic viscosity", nu, false),
            ("surface tension", sigma, true),
            ("flow per width", q, true),
            ("tangential acceleration", gt, false),
            ("normal acceleration", gn, true),
        ] {
            if !value.is_finite() || value < 0.0 || (!zero_allowed && value == 0.0) {
                return Err(FilmError::InvalidInput(name));
            }
        }
        if q == 0.0 {
            return Err(FilmError::DryFilm);
        }
        // Log scaling avoids underflow in nu*q and overflow in 3*nu*q.
        // Use the existing deterministic owner, never platform cbrt/powf.
        let h = fs_math::det::exp(
            (fs_math::det::ln(3.0) + fs_math::det::ln(nu)
                + fs_math::det::ln(q)
                - fs_math::det::ln(gt))
                / 3.0,
        );
        let mean_velocity = q / h;
        let reynolds = q / nu;
        let wall_shear = (rho * h) * gt;
        for value in [h, mean_velocity, 1.5 * mean_velocity, reynolds, wall_shear] {
            represented(value, true)?;
        }
        let normal_ratio = gn / gt;
        represented(normal_ratio, gn > 0.0)?;
        Ok(Self {
            material,
            flow: q,
            tangent: gt,
            normal_ratio,
            thickness: h,
            mean_velocity,
            reynolds,
            wall_shear,
        })
    }

    /// Equilibrium film thickness in metres.
    #[must_use]
    pub fn thickness_m(self) -> f64 {
        self.thickness
    }
    /// Mean downstream velocity in m/s.
    #[must_use]
    pub fn mean_velocity_m_s(self) -> f64 {
        self.mean_velocity
    }
    /// Shear-free surface velocity in m/s.
    #[must_use]
    pub fn surface_velocity_m_s(self) -> f64 {
        1.5 * self.mean_velocity
    }
    /// Reynolds number on mean speed and film thickness: q/nu.
    #[must_use]
    pub fn reynolds(self) -> f64 {
        self.reynolds
    }
    /// Tangential wall-shear magnitude in Pa.
    #[must_use]
    pub fn wall_shear_pa(self) -> f64 {
        self.wall_shear
    }

    /// Velocity at wall-normal distance y in [0,h], in m/s.
    pub fn velocity_m_s(self, y_m: f64) -> Result<f64, FilmError> {
        let z = self.coordinate(y_m)?;
        let value = self.surface_velocity_m_s() * z * (2.0 - z);
        represented(value, y_m > 0.0)?;
        Ok(value)
    }

    /// Shear stress magnitude at wall-normal distance y in [0,h], in Pa.
    pub fn shear_pa(self, y_m: f64) -> Result<f64, FilmError> {
        self.coordinate(y_m)?;
        let value = self.wall_shear * ((self.thickness - y_m) / self.thickness);
        represented(value, y_m < self.thickness)?;
        Ok(value)
    }

    /// Screen one nonnegative physical wavenumber. k=0 is exactly neutral by
    /// mass conservation. Nonzero modes must pass both caller-supplied budgets;
    /// no extrapolation or asymptotic-accuracy certificate is supplied.
    pub fn long_wave_mode(
        self,
        wavenumber_rad_m: f64,
        limits: LongWaveLimits,
    ) -> Result<FilmMode, FilmError> {
        let k = wavenumber_rad_m;
        if !k.is_finite() || k < 0.0 {
            return Err(FilmError::InvalidInput("wavenumber"));
        }
        for limit in [limits.max_kh, limits.max_re_kh] {
            if !limit.is_finite() || limit <= 0.0 || limit >= 1.0 {
                return Err(FilmError::InvalidInput("long-wave budget"));
            }
        }
        let kh = k * self.thickness;
        let re_kh = self.reynolds * kh;
        if kh > limits.max_kh || re_kh > limits.max_re_kh {
            return Err(FilmError::OutsideLongWaveRegime);
        }
        if k == 0.0 {
            return Ok(FilmMode {
                growth_per_s: 0.0,
                angular_frequency_rad_s: 0.0,
                kh,
                re_kh,
            });
        }
        represented(kh, true)?;
        represented(re_kh, true)?;
        let rate = (self.flow * k) * k;
        let frequency = 3.0 * (self.mean_velocity * k);
        represented(rate, true)?;
        represented(frequency, true)?;
        let sigma = self.material.surface_tension_n_m;
        let capillary = if sigma == 0.0 {
            0.0
        } else {
            (((sigma / self.material.density_kg_m3) / self.tangent) * k) * k
        };
        represented(capillary, sigma > 0.0)?;
        let driving = (6.0_f64 / 5.0).mul_add(self.reynolds, -self.normal_ratio) - capillary;
        let growth = rate * driving;
        represented(growth.abs(), driving != 0.0)?;
        Ok(FilmMode {
            growth_per_s: growth,
            angular_frequency_rad_s: frequency,
            kh,
            re_kh,
        })
    }

    fn coordinate(self, y: f64) -> Result<f64, FilmError> {
        if !y.is_finite() || y < 0.0 || y > self.thickness {
            return Err(FilmError::InvalidInput("wall-normal coordinate"));
        }
        Ok(y / self.thickness)
    }
}

fn represented(value: f64, nonzero: bool) -> Result<(), FilmError> {
    if !value.is_finite() || (nonzero && value == 0.0) {
        Err(FilmError::NonRepresentable)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const WATER: FilmMaterial = FilmMaterial {
        density_kg_m3: 1000.0,
        kinematic_viscosity_m2_s: 1e-6,
        surface_tension_n_m: 0.072,
    };
    const LIMITS: LongWaveLimits = LongWaveLimits {
        max_kh: 0.2,
        max_re_kh: 0.2,
    };
    fn near(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 2e-12 * a.abs().max(b.abs()).max(f64::MIN_POSITIVE),
            "{a:e} != {b:e}"
        );
    }
    #[test]
    fn profile_closes_flow_weight_and_dissipation() {
        let q = 2e-6;
        let gt = 4.0;
        let film = NusseltFilm::new(WATER, q, gt, 3.0).unwrap();
        let h = film.thickness_m();
        near(h * h * h, 3.0 * WATER.kinematic_viscosity_m2_s * q / gt);
        assert_eq!(film.velocity_m_s(0.0), Ok(0.0));
        assert_eq!(film.shear_pa(h), Ok(0.0));
        near(film.wall_shear_pa(), WATER.density_kg_m3 * gt * h);
        // Three-point Simpson quadrature is exact for these quadratic fields.
        let mut flow = 0.0;
        let mut dissipation = 0.0;
        for (z, w) in [(0.0, 1.0), (0.5, 4.0), (1.0, 1.0)] {
            let u = film.velocity_m_s(z * h).unwrap();
            let tau = film.shear_pa(z * h).unwrap();
            flow += w * u * h / 6.0;
            dissipation += w * tau * tau * h
                / (6.0 * WATER.density_kg_m3 * WATER.kinematic_viscosity_m2_s);
        }
        near(flow, q);
        near(dissipation, WATER.density_kg_m3 * gt * q);
        near(film.surface_velocity_m_s(), 1.5 * q / h);
    }
    #[test]
    fn neutral_threshold_and_capillary_cutoff_follow_independent_coefficients() {
        let dry_capillary = FilmMaterial { surface_tension_n_m: 0.0, ..WATER };
        for factor in [0.99, 1.01] {
            let q = (5.0 / 6.0) * WATER.kinematic_viscosity_m2_s * factor;
            let film = NusseltFilm::new(dry_capillary, q, 4.0, 4.0).unwrap();
            let mode = film.long_wave_mode(10.0, LIMITS).unwrap();
            assert_eq!(mode.growth_per_s > 0.0, factor > 1.0);
        }
        let film = NusseltFilm::new(WATER, 1e-6, 4.0, 4.0).unwrap();
        let h = film.thickness_m();
        let inertial = 2.0 * 4.0 * 4.0 * h.powi(6)
            / (15.0 * WATER.kinematic_viscosity_m2_s.powi(3));
        let hydrostatic = 4.0 * h.powi(3) / (3.0 * WATER.kinematic_viscosity_m2_s);
        let capillary = WATER.surface_tension_n_m * h.powi(3)
            / (3.0 * WATER.density_kg_m3 * WATER.kinematic_viscosity_m2_s);
        let cutoff = ((inertial - hydrostatic) / capillary).sqrt();
        for factor in [0.5, 2.0] {
            let k = cutoff * factor;
            let mode = film.long_wave_mode(k, LIMITS).unwrap();
            near(mode.growth_per_s, (inertial - hydrostatic) * k * k - capillary * k.powi(4));
            assert_eq!(mode.growth_per_s > 0.0, factor < 1.0);
            near(mode.angular_frequency_rad_s, 3.0 * film.mean_velocity_m_s() * k);
        }
        assert_eq!(film.long_wave_mode(0.0, LIMITS).unwrap().growth_per_s, 0.0);
    }
    #[test]
    fn material_scaling_and_small_flow_do_not_invent_a_reynolds_floor() {
        let first = NusseltFilm::new(WATER, 1e-6, 4.0, 3.0).unwrap();
        let viscous = FilmMaterial { kinematic_viscosity_m2_s: 8e-6, ..WATER };
        let second = NusseltFilm::new(viscous, 1e-6, 4.0, 3.0).unwrap();
        near(second.thickness_m(), 2.0 * first.thickness_m());
        near(second.mean_velocity_m_s(), 0.5 * first.mean_velocity_m_s());
        near(second.reynolds(), first.reynolds() / 8.0);
        let tiny = NusseltFilm::new(WATER, 1e-320, 4.0, 3.0).unwrap();
        assert!(tiny.thickness_m() > 0.0 && tiny.reynolds() > 0.0);
        assert!(tiny.reynolds() < 1e-300);
    }
    #[test]
    fn invalid_or_unresolved_modes_refuse_without_clamping() {
        assert_eq!(NusseltFilm::new(WATER, 0.0, 4.0, 3.0), Err(FilmError::DryFilm));
        assert!(NusseltFilm::new(WATER, 1e-6, 0.0, 3.0).is_err());
        assert!(NusseltFilm::new(WATER, f64::NAN, 4.0, 3.0).is_err());
        let film = NusseltFilm::new(WATER, 2e-6, 4.0, 3.0).unwrap();
        assert!(film.velocity_m_s(-1.0).is_err());
        assert!(film.shear_pa(film.thickness_m() * 2.0).is_err());
        assert_eq!(film.long_wave_mode(0.15 / film.thickness_m(), LIMITS),
            Err(FilmError::OutsideLongWaveRegime)); // kh passes; Re*kh fails.
        assert_eq!(film.long_wave_mode(0.25 / film.thickness_m(), LIMITS),
            Err(FilmError::OutsideLongWaveRegime));
        assert!(film.long_wave_mode(f64::NAN, LIMITS).is_err());
        assert!(film.long_wave_mode(1.0, LongWaveLimits { max_kh: 1.0, ..LIMITS }).is_err());
    }
}
