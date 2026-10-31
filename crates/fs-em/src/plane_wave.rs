//! Exact monochromatic scattering by one planar isotropic dielectric interface.
//!
//! This explicit plane-wave rung is separate from the mesh MQS solvers.
//! Sources declare either a nondispersive lossless law or the existing passive
//! Debye law; selecting either is not experimental qualification at a frequency.
//! Both sides retain their complete source states, including absolute SI mu.
//!
//! Convention: exp(+i omega t), forward exp(-i k_z z), interface at z=0,
//! incidence from z<0. Amplitudes use the SAME tangential electric-field axis
//! for incident, reflected and transmitted waves. TM is therefore not the
//! signed rotating p-polarization basis commonly used by optical renderers.
//! Continuity gives r=(Y_i-Y_t)/(Y_i+Y_t), t=1+r, with
//! Y_TE=k_z/(omega mu), Y_TM=omega epsilon/k_z.
//!
//! Reference: MIT 6.974, "Electromagnetic Waves and Interfaces", §§2.2.1/2.2.3:
//! <https://ocw.mit.edu/courses/6-974-fundamentals-of-photonics-quantum-electronics-spring-2006/fafcb9e36b7749f4eeb099d666bdf70a_electr_wav_inter.pdf>.
//! No finite layers, roughness, anisotropy, conducting leakage, magnetic loss,
//! transient wave evolution or measured-material accuracy claim is supplied.

use fs_material::{
    field_storage::{
        FieldStorageError, dielectric::ResolvedDebyeCapacitor,
        linear_electric_permittivity, linear_magnetic_permeability,
    },
    state_point::ResolvedMaterialStatePoint,
};
use fs_math::{c64::C64, det};

/// Electric field transverse to the incidence plane, or magnetic field transverse.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Polarization {
    Te,
    Tm,
}

/// Explicit existing material law for the transmitted half-space.
#[derive(Clone, Copy, Debug)]
pub enum DielectricResponse<'a> {
    /// Positive, isotropic, nondispersive lossless epsilon and mu.
    Lossless(&'a ResolvedMaterialStatePoint),
    /// Original polarization response, with a lossless sourced mu in the SAME
    /// retained bundle. The capacitor geometry cancels from bulk permittivity.
    Debye(&'a ResolvedDebyeCapacitor),
}
impl<'a> DielectricResponse<'a> {
    pub fn material(self) -> &'a ResolvedMaterialStatePoint {
        match self {
            Self::Lossless(material) => material,
            Self::Debye(model) => model.material(),
        }
    }

    fn permittivity(self, omega: f64) -> Result<C64, ReflectionError> {
        let epsilon = match self {
            Self::Lossless(material) => C64::from_re(linear_electric_permittivity(material)?),
            Self::Debye(model) => {
                // Reuse the original frequency response: Y=i*omega*epsilon*A/d.
                // No second pole law, loss tangent or guessed conductivity.
                let y = model.admittance(omega)?;
                let geometry = positive(model.gap_m() / model.area_m2())?;
                C64::new(
                    (y.susceptance_s / omega) * geometry,
                    -(y.conductance_s / omega) * geometry,
                )
            }
        };
        if !finite(epsilon) || epsilon.re <= 0.0 || epsilon.im > 0.0 {
            return Err(ReflectionError::Invalid("dielectric response must be finite and passive"));
        }
        Ok(epsilon)
    }
}

/// Source/model-domain or unrepresentable scattering refusal.
#[derive(Debug)]
pub enum ReflectionError {
    Material(FieldStorageError),
    Invalid(&'static str),
}
impl From<FieldStorageError> for ReflectionError {
    fn from(value: FieldStorageError) -> Self {
        Self::Material(value)
    }
}
impl core::fmt::Display for ReflectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar dielectric reflection: {self:?}")
    }
}
impl std::error::Error for ReflectionError {}

/// Sources and derived constant incident data are retained together.
/// Incident loss must be absent so incoming and reflected normal powers have
/// an unambiguous additive decomposition. The transmitted side may absorb.
pub struct PlanarDielectricInterface<'a> {
    incident: &'a ResolvedMaterialStatePoint,
    transmitted: DielectricResponse<'a>,
    incident_epsilon_f_m: f64,
    incident_mu_h_m: f64,
    transmitted_mu_h_m: f64,
}

/// Complex field amplitudes and time-averaged NORMAL Poynting power fractions.
/// Transmitted power is flux entering the half-space at z=0, not power surviving
/// at infinity. For an absorbing half-space it is ultimately dissipated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InterfaceReflection {
    pub reflected_tangential_e: C64,
    pub transmitted_tangential_e: C64,
    pub reflected_power_fraction: f64,
    pub transmitted_power_fraction: f64,
    /// Positive real incident normal wavenumber [rad/m].
    pub incident_normal_wavenumber_rad_m: f64,
    /// Re >= 0, Im <= 0: propagation/decay into z>0, including TIR.
    pub transmitted_normal_wavenumber_rad_m: C64,
    /// Actual bulk constitutive value, using exp(+i omega t).
    pub transmitted_permittivity_f_m: C64,
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn positive(x: f64) -> Result<f64, ReflectionError> {
    if x.is_finite() && x > 0.0 {
        Ok(x)
    } else {
        Err(ReflectionError::Invalid("positive representable wave data required"))
    }
}

impl<'a> PlanarDielectricInterface<'a> {
    /// Uses the original field-storage source admission for lossless constants.
    /// Missing mu, nonlinear/field-dependent coefficients or unsupported source
    /// coordinates refuse; vacuum or nonmagnetic defaults are never inserted.
    pub fn try_new(
        incident: &'a ResolvedMaterialStatePoint,
        transmitted: DielectricResponse<'a>,
    ) -> Result<Self, ReflectionError> {
        let incident_epsilon_f_m = linear_electric_permittivity(incident)?;
        let incident_mu_h_m = linear_magnetic_permeability(incident)?;
        let transmitted_mu_h_m = linear_magnetic_permeability(transmitted.material())?;
        if let DielectricResponse::Lossless(material) = transmitted {
            linear_electric_permittivity(material)?;
        }
        Ok(Self {
            incident,
            transmitted,
            incident_epsilon_f_m,
            incident_mu_h_m,
            transmitted_mu_h_m,
        })
    }

    pub const fn incident_material(&self) -> &'a ResolvedMaterialStatePoint {
        self.incident
    }
    pub fn transmitted_material(&self) -> &'a ResolvedMaterialStatePoint {
        self.transmitted.material()
    }

    /// Allocation-free after material preparation; bounded by the original
    /// MAX_DIELECTRIC_POLES. Exact grazing (zero incident normal power) and DC
    /// refuse; arbitrarily close admitted positive incidence remains meaningful.
    pub fn response(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        polarization: Polarization,
    ) -> Result<InterfaceReflection, ReflectionError> {
        let omega = positive(omega_rad_s)?;
        let c = positive(incidence_cosine)?;
        if c > 1.0 {
            return Err(ReflectionError::Invalid("incidence cosine must lie in (0,1]"));
        }
        let epsilon = self.transmitted.permittivity(omega)?;
        let epsilon_ratio = C64::new(
            epsilon.re / self.incident_epsilon_f_m,
            epsilon.im / self.incident_epsilon_f_m,
        );
        let mu_ratio = positive(self.transmitted_mu_h_m / self.incident_mu_h_m)?;
        let index_squared_ratio = epsilon_ratio.scale(mu_ratio);
        if !finite(index_squared_ratio)
            || !finite(epsilon_ratio)
            || epsilon_ratio.re <= 0.0
            || index_squared_ratio.re <= 0.0
        {
            return Err(ReflectionError::Invalid("material contrast is unrepresentable"));
        }
        // k_z,t / k_i. This form retains c^2 for identical media even when
        // 1-c^2 rounds to 1 close to grazing incidence.
        let normal_real = if c > 0.5 {
            // Preserve a small positive index contrast at normal incidence.
            index_squared_ratio.re - (1.0 - c) * (1.0 + c)
        } else {
            (index_squared_ratio.re - 1.0) + c * c
        };
        let normal_squared = C64::new(normal_real, index_squared_ratio.im);
        let mut normal = if index_squared_ratio == C64::ONE {
            C64::from_re(c)
        } else {
            normal_squared.sqrt()
        };
        // The shared principal root chooses +i for a negative real argument.
        // With exp(+i omega t), the physical evanescent root is instead -i.
        if normal.im > 0.0 {
            normal = -normal;
        }
        let wave_number = positive(
            omega * det::sqrt(self.incident_epsilon_f_m) * det::sqrt(self.incident_mu_h_m),
        )?;
        let incident_normal = positive(wave_number * c)?;
        let transmitted_normal = normal.scale(wave_number);
        if !finite(normal) || !finite(transmitted_normal) || normal.re < 0.0 {
            return Err(ReflectionError::Invalid("outgoing wave number is unrepresentable"));
        }

        // Multiply through by k_z for TM: no singular admittance evaluation
        // at the critical angle. a/b = Y_i/Y_t for both polarizations.
        let (a, b) = match polarization {
            Polarization::Te => (C64::from_re(mu_ratio * c), normal),
            Polarization::Tm => (normal, epsilon_ratio.scale(c)),
        };
        let scale = positive(a.re.abs().max(a.im.abs()).max(b.re.abs()).max(b.im.abs()))?;
        let a = C64::new(a.re / scale, a.im / scale);
        let b = C64::new(b.re / scale, b.im / scale);
        let denominator = a + b;
        let norm = positive(denominator.norm_sq())?;
        let reflected = (a - b) / denominator;
        let transmitted = a.scale(2.0) / denominator;
        let reflected_power = reflected.norm_sq();
        // Stable Poynting flux, not |t|^2 alone and not a clipped 1-|r|^2.
        let transmitted_power = 4.0 * (a.re * b.re + a.im * b.im) / norm;
        if !finite(reflected)
            || !finite(transmitted)
            || !transmitted_power.is_finite()
            || transmitted_power < 0.0
            || (reflected_power + transmitted_power - 1.0).abs() > 128.0 * f64::EPSILON
        {
            return Err(ReflectionError::Invalid("passive interface power balance failed"));
        }
        Ok(InterfaceReflection {
            reflected_tangential_e: reflected,
            transmitted_tangential_e: transmitted,
            reflected_power_fraction: reflected_power,
            transmitted_power_fraction: transmitted_power,
            incident_normal_wavenumber_rad_m: incident_normal,
            transmitted_normal_wavenumber_rad_m: transmitted_normal,
            transmitted_permittivity_f_m: epsilon,
        })
    }
}
