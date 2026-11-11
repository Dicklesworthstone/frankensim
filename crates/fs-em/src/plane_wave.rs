//! Exact monochromatic scattering by one planar isotropic dielectric interface.
//!
//! This explicit plane-wave rung is separate from the mesh MQS solvers.
//! Sources declare a nondispersive lossless law, the existing passive Debye
//! law, source-bound inertial Lorentz polarization, an explicitly selected
//! isotropic Ohmic conductor with displacement current, an original spectral
//! complex-index metal response, or explicitly
//! prepared collisional plasma. Selecting a
//! law is not experimental qualification at a frequency.
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
//! Ohmic displacement current: MIT 6.013, Electromagnetics and Applications,
//! §9.2.4, equations 9.2.37–9.2.39:
//! <https://ocw.mit.edu/courses/6-013-electromagnetics-and-applications-spring-2009/d3be4ea78b036a6362230fb41780cf54_MIT6_013S09_notes.pdf>.
//! See [`layers`] for finite coherent isotropic films. The separate [`uniaxial`]
//! owner supplies transparent-crystal refraction in a principal incidence plane.
//! No roughness, magnetic loss, transient wave evolution or measured-material
//! accuracy claim is supplied by these interfaces.

use fs_material::{
    conductor::{ResolvedConductorState, planar::PlanarConductorMaterial},
    field_storage::{
        FieldStorageError,
        dielectric::{ResolvedDebyeCapacitor, lorentz::ResolvedLorentzMaterial},
        linear_electric_permittivity, linear_magnetic_permeability,
    },
    state_point::ResolvedMaterialStatePoint,
};
use fs_math::{c64::C64, dd::Dd, det};

/// Coherent finite dielectric films, with the same source laws and conventions.
pub mod layers;

/// Local interface fields and signed power when the incident medium absorbs.
pub mod absorbing_incidence;

/// Source-bound ordinary/extraordinary refraction in a principal crystal section.
pub mod uniaxial;

/// Coherent vector diffraction of sampled beams through the original films.
pub mod angular_spectrum;

/// Analytic frequency derivatives and propagation delays of original bulk laws.
pub mod dispersion;

/// Original spectral metal/dielectric optics with explicit same-state permeability.
pub mod optical;
pub use optical::{
    MolecularGasOpticalResponse, MolecularGasOpticalSource, OpticalConductorResponse,
    OpticalMaterialResponse, RefreshedSpectralOpticalSource, SilicaOpticalResponse,
    SpectralOpticalResponse, SpectralOpticalSource, ThermoOpticOpticalResponse,
    UniaxialOpticalResponse, WaterOpticalResponse,
};

/// Original collisional carrier transport in an admitted outgoing half-space wave.
pub mod plasma;
pub use plasma::{IonizedGasWavePoint, IonizedGasWaveResponse, PreparedIonizedGasWave};

// Shared power-to-momentum contractions. Source/geometry admission belongs
// to the original film and beam owners, not to an arbitrary power input.
mod momentum;

/// Electric field transverse to the incidence plane, or magnetic field transverse.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Polarization {
    Te,
    Tm,
}

/// Explicit material law for a homogeneous isotropic electromagnetic medium.
#[derive(Clone, Copy, Debug)]
pub enum DielectricResponse<'a> {
    /// Positive, isotropic, nondispersive lossless epsilon and mu.
    Lossless(&'a ResolvedMaterialStatePoint),
    /// Original polarization response, with a lossless sourced mu in the SAME
    /// retained bundle. The capacitor geometry cancels from bulk permittivity.
    Debye(&'a ResolvedDebyeCapacitor),
    /// Original bound-charge oscillators, including inertial storage and
    /// physical damping. Real permittivity may be negative above resonance.
    /// This steady response does not replace the retained transient history.
    Lorentz(&'a ResolvedLorentzMaterial),
    /// Original scalar Ohmic law J=sigma E plus lossless displacement current,
    /// with epsilon and mu from the SAME complete conductor material bundle.
    /// For exp(+i omega t), epsilon_eff = epsilon - i sigma/omega. This admits
    /// skin-depth attenuation and Joule loss, without treating a DC source as
    /// measured optical data or inventing carrier dispersion or magnetic loss.
    Ohmic(&'a ResolvedConductorState),
    /// Original Debye polarization and Ohmic conduction acting in parallel.
    /// Both owners must retain the identical complete material snapshot; a
    /// conductivity from another source state cannot be attached to this law.
    /// Polarization loss and Joule loss add once in the effective permittivity.
    DebyeWithOhmic {
        dielectric: &'a ResolvedDebyeCapacitor,
        conductor: &'a ResolvedConductorState,
    },
    /// Independent bound-charge resonance and DC Ohmic transport from the
    /// identical complete material snapshot. Each loss channel enters once.
    LorentzWithOhmic {
        dielectric: &'a ResolvedLorentzMaterial,
        conductor: &'a ResolvedConductorState,
    },
    /// Visible complex-index source with its original same-state permeability.
    /// Includes the source's electric absorption once; real epsilon may be
    /// negative. Frequencies outside the source's spectral band refuse.
    OpticalConductor(OpticalConductorResponse<'a>),
    /// Original metal, Cauchy, tabulated or Sellmeier optical family. Dielectric
    /// intensity absorption supplies k=alpha*lambda_vac/(4*pi); RGB-only input
    /// refuses. The same complete source supplies real scalar permeability.
    OpticalMaterial(OpticalMaterialResponse<'a>),
    /// Original finite-band complex-index or phase-index/absorption source.
    /// The authored interpolation-gap limit and same-state permeability remain
    /// attached; no visible-band restriction or extrapolated tail is supplied.
    SpectralOptical(SpectralOpticalResponse<'a>),
    /// Original IAPWS ordinary-water phase dispersion at retained T/rho/mu,
    /// with explicitly zero loss or independently sourced current absorption.
    /// This index release supplies no causal transient reservoir.
    WaterOptical(WaterOpticalResponse<'a>),
    /// Original NASA CHARMS Corning 7980 phase dispersion at retained T/mu,
    /// with explicitly zero loss or independently sourced current absorption.
    /// Original thermal/group derivatives remain; no causal reservoir is added.
    SilicaOptical(SilicaOpticalResponse<'a>),
    /// Original dielectric reference with source-owned SCHOTT temperature
    /// correction and actual-state mu. Explicitly zero loss or retained original
    /// spectral absorption; no inferred thermal loss or causal reservoir.
    ThermoOpticOptical(ThermoOpticOpticalResponse<'a>),
    /// Original molecular phase/group law and an explicitly selected magnetic
    /// source, prepared at one frequency. Lossless coherent approximation;
    /// molecular scattering is not misrepresented as absorbed heat.
    MolecularGasOptical(MolecularGasOpticalResponse<'a>),
    /// Original collisional free-carrier response plus an explicitly sourced
    /// lossless bound-charge background, prepared at one angular frequency.
    /// Normalized half-space scattering does not admit a physical amplitude;
    /// prepare_ionized_gas_wave checks the complete actual surface drift.
    /// Finite films refuse until their whole-field drift admission is supplied.
    IonizedGas(IonizedGasWaveResponse<'a>),
}
impl<'a> DielectricResponse<'a> {
    /// Original scalar-property bundle. Ionized gas retains its bound-charge
    /// background; molecular optics retains its separately selected magnetic
    /// source. The enum also retains the actual carrier or molecular owner.
    pub fn material(self) -> &'a ResolvedMaterialStatePoint {
        match self {
            Self::Lossless(material) => material,
            Self::Debye(model) => model.material(),
            Self::Lorentz(model) => model.material(),
            Self::Ohmic(conductor) => conductor.material(),
            Self::DebyeWithOhmic { dielectric, .. } => dielectric.material(),
            Self::LorentzWithOhmic { dielectric, .. } => dielectric.material(),
            Self::OpticalConductor(response) => response.material(),
            Self::OpticalMaterial(response) => response.material(),
            Self::SpectralOptical(response) => response.material(),
            Self::WaterOptical(response) => response.material(),
            Self::SilicaOptical(response) => response.material(),
            Self::ThermoOpticOptical(response) => response.material(),
            Self::MolecularGasOptical(response) => response.material(),
            Self::IonizedGas(response) => response.material(),
        }
    }

    /// The complex permittivity (F/m) this response presents at angular
    /// frequency `omega`: the value the interface, film and slab solvers
    /// use, with the same refusals.
    pub fn permittivity_f_m(self, omega: f64) -> Result<C64, ReflectionError> {
        self.permittivity(omega)
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
            Self::Lorentz(model) => model.permittivity_f_m(omega)?,
            Self::Ohmic(conductor) => {
                let epsilon = linear_electric_permittivity(conductor.material())?;
                let sigma = ohmic_conductivity(conductor)?;
                // Maxwell-Ampere: curl H = (sigma + i omega epsilon) E.
                // Reuse the original scalar source and its independent basis;
                // the film/field owners then account for sigma |E|^2 / 2 once.
                C64::new(epsilon, -positive(sigma / omega)?)
            }
            Self::DebyeWithOhmic {
                dielectric,
                conductor,
            } => {
                let sigma = shared_ohmic_conductivity(dielectric, conductor)?;
                let epsilon = Self::Debye(dielectric).permittivity(omega)?;
                C64::new(epsilon.re, epsilon.im - positive(sigma / omega)?)
            }
            Self::LorentzWithOhmic {
                dielectric,
                conductor,
            } => {
                let sigma = shared_lorentz_ohmic_conductivity(dielectric, conductor)?;
                let epsilon = dielectric.permittivity_f_m(omega)?;
                C64::new(epsilon.re, epsilon.im - positive(sigma / omega)?)
            }
            Self::OpticalConductor(response) => response.permittivity_f_m(omega)?,
            Self::OpticalMaterial(response) => response.permittivity_f_m(omega)?,
            Self::SpectralOptical(response) => response.permittivity_f_m(omega)?,
            Self::WaterOptical(response) => response.permittivity_f_m(omega)?,
            Self::SilicaOptical(response) => response.permittivity_f_m(omega)?,
            Self::ThermoOpticOptical(response) => response.permittivity_f_m(omega)?,
            Self::MolecularGasOptical(response) => response.permittivity_f_m(omega)?,
            Self::IonizedGas(response) => response.permittivity(omega)?,
        };
        if !finite(epsilon) || epsilon == C64::ZERO || epsilon.im > 0.0 {
            return Err(ReflectionError::Invalid(
                "dielectric response must be finite and passive",
            ));
        }
        Ok(epsilon)
    }

    // Incoming/reflected normal powers have the original additive port
    // interpretation only in a lossless incident medium. A dispersive source
    // is still valid at a monochromatic frequency; zero loss is checked there,
    // without replacing its spectral law by a separately authored constant.
    fn lossless_permittivity(self, omega: f64) -> Result<f64, ReflectionError> {
        let epsilon = self.permittivity(omega)?;
        if epsilon.im != 0.0 || epsilon.re <= 0.0 {
            return Err(ReflectionError::Invalid(
                "incident wave response must be positive and lossless at the selected frequency",
            ));
        }
        Ok(epsilon.re)
    }
}

// Same static source admission as the native spatial Ohmic field owner. A
// field-, frequency- or strain-dependent transport snapshot is not a constant
// instantaneous J=sigma E law. Retain explicitly owned stationary phase data.
fn ohmic_conductivity(conductor: &ResolvedConductorState) -> Result<f64, ReflectionError> {
    if !conductor.has_admitted_static_sources() {
        return Err(ReflectionError::Invalid(
            "Ohmic waves require the original stationary scalar conductivity law",
        ));
    }
    positive(conductor.resistivity_ohm_m().recip())
}

fn shared_ohmic_conductivity(
    dielectric: &ResolvedDebyeCapacitor,
    conductor: &ResolvedConductorState,
) -> Result<f64, ReflectionError> {
    if dielectric.material().identity() != conductor.material().identity() {
        return Err(ReflectionError::Invalid(
            "Debye and Ohmic wave laws require the same complete material snapshot",
        ));
    }
    ohmic_conductivity(conductor)
}

fn shared_lorentz_ohmic_conductivity(
    dielectric: &ResolvedLorentzMaterial,
    conductor: &ResolvedConductorState,
) -> Result<f64, ReflectionError> {
    if dielectric.material().identity() != conductor.material().identity() {
        return Err(ReflectionError::Invalid(
            "Lorentz and Ohmic wave laws require the same complete material snapshot",
        ));
    }
    ohmic_conductivity(conductor)
}

/// Source/model-domain or unrepresentable scattering refusal.
#[derive(Debug)]
pub enum ReflectionError {
    Material(FieldStorageError),
    IonizedGas(fs_material::ionized_gas::IonizedGasError),
    Invalid(&'static str),
}
impl From<FieldStorageError> for ReflectionError {
    fn from(value: FieldStorageError) -> Self {
        Self::Material(value)
    }
}
impl From<fs_material::ionized_gas::IonizedGasError> for ReflectionError {
    fn from(value: fs_material::ionized_gas::IonizedGasError) -> Self {
        Self::IonizedGas(value)
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
    incident_response: DielectricResponse<'a>,
    transmitted: DielectricResponse<'a>,
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
    /// Transmitted power-conjugate H divided by incident H. Its independent
    /// Fresnel numerator preserves a small H when reflection rounds to one,
    /// including critical TM fields; it is not formed by subtracting 1-r.
    pub transmitted_conjugate_h: C64,
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
        Err(ReflectionError::Invalid(
            "positive representable wave data required",
        ))
    }
}

/// Frequency-resolved material contrast shared by half-spaces and finite films.
#[derive(Clone, Copy)]
struct WaveMedium {
    epsilon_f_m: C64,
    epsilon_ratio: C64,
    mu_ratio: f64,
    index_squared_ratio: C64,
    // Preserve the constitutive product before subtracting the conserved
    // tangential wave number. A rounded n^2 can otherwise cross the critical
    // angle or falsely identify unequal media as a grazing null interface.
    index_squared_real_residual: f64,
}
impl WaveMedium {
    fn new(
        epsilon: C64,
        mu: f64,
        incident_epsilon: f64,
        incident_mu: f64,
    ) -> Result<Self, ReflectionError> {
        let epsilon_ratio = C64::new(epsilon.re / incident_epsilon, epsilon.im / incident_epsilon);
        let mu_ratio = positive(mu / incident_mu)?;
        let index_squared_ratio = epsilon_ratio.scale(mu_ratio);
        if !finite(index_squared_ratio)
            || !finite(epsilon_ratio)
            || epsilon_ratio == C64::ZERO
            || index_squared_ratio == C64::ZERO
            || epsilon_ratio.im > 0.0
        {
            return Err(ReflectionError::Invalid(
                "material contrast is unrepresentable",
            ));
        }
        Ok(Self {
            epsilon_f_m: epsilon,
            epsilon_ratio,
            mu_ratio,
            index_squared_ratio,
            index_squared_real_residual: epsilon_ratio
                .re
                .mul_add(mu_ratio, -index_squared_ratio.re),
        })
    }

    fn normal_squared(self, cosine: f64) -> C64 {
        // Preserve small positive contrast at normal incidence, and c^2 for
        // equal indices close to grazing where 1-c^2 rounds to one.
        let normal_real = if cosine > 0.5 {
            self.index_squared_ratio.re - (1.0 - cosine) * (1.0 + cosine)
        } else {
            (self.index_squared_ratio.re - 1.0) + cosine * cosine
        };
        // Only the cancellation neighbourhood needs the precision ladder.
        // Away from it the ordinary scalar path stays unchanged. Inside it,
        // retain both the material-product residual and the exact cosine
        // square; rounding either first can erase a real TIR phase or weak
        // positive transmission. This is numerical accuracy for the admitted
        // material contrast, not additional accuracy of its source data.
        let normal_real = if normal_real.abs()
            <= 8.0 * f64::EPSILON * self.index_squared_ratio.re.abs().max(1.0)
        {
            let index = Dd::from_pair(
                self.index_squared_ratio.re,
                self.index_squared_real_residual,
            );
            let c = Dd::from_f64(cosine);
            if cosine > 0.5 {
                (index - (Dd::ONE - c) * (Dd::ONE + c)).to_f64()
            } else {
                ((index - Dd::ONE) + c * c).to_f64()
            }
        } else {
            normal_real
        };
        C64::new(normal_real, self.index_squared_ratio.im)
    }

    fn normal(self, cosine: f64) -> Result<C64, ReflectionError> {
        let normal =
            if self.index_squared_ratio == C64::ONE && self.index_squared_real_residual == 0.0 {
                C64::from_re(cosine)
            } else {
                self.normal_squared(cosine).sqrt()
            };
        checked_outgoing_normal(normal)
    }
}

fn checked_outgoing_normal(mut normal: C64) -> Result<C64, ReflectionError> {
    // The principal square root has +i on the negative real axis; the
    // outgoing exp(+i omega t) evanescent branch instead has -i.
    if normal.im > 0.0 {
        normal = -normal;
    }
    if !finite(normal) || normal.re < 0.0 {
        return Err(ReflectionError::Invalid(
            "outgoing wave number is unrepresentable",
        ));
    }
    Ok(normal)
}

// Shared bounded boundary pair. The caller gives a/b=Y_i/Y_t without
// evaluating a singular transmitted TM admittance at critical incidence.
fn normalized_interface_pair(a: C64, b: C64) -> Result<(C64, C64, C64, f64), ReflectionError> {
    let scale = positive(a.re.abs().max(a.im.abs()).max(b.re.abs()).max(b.im.abs()))?;
    let a = C64::new(a.re / scale, a.im / scale);
    let b = C64::new(b.re / scale, b.im / scale);
    let denominator = a + b;
    let norm = positive(denominator.norm_sq())?;
    Ok((a, b, denominator, norm))
}

impl<'a> PlanarDielectricInterface<'a> {
    /// Uses the original field-storage source admission for lossless constants.
    /// Missing mu, nonlinear/field-dependent coefficients or unsupported source
    /// coordinates refuse; vacuum or nonmagnetic defaults are never inserted.
    pub fn try_new(
        incident: &'a ResolvedMaterialStatePoint,
        transmitted: DielectricResponse<'a>,
    ) -> Result<Self, ReflectionError> {
        Self::try_new_with_incident_response(DielectricResponse::Lossless(incident), transmitted)
    }

    /// Retain both original constitutive laws, including spectral incident
    /// glass. The angle-based `response` requires positive, exactly lossless
    /// incident permittivity at each queried frequency. The explicit transverse
    /// wave-number query separately admits absorbing incidence with interference
    /// flux. Out-of-band optical queries refuse; an optical source does not
    /// require a second static electric-permittivity property.
    pub fn try_new_with_incident_response(
        incident: DielectricResponse<'a>,
        transmitted: DielectricResponse<'a>,
    ) -> Result<Self, ReflectionError> {
        let incident_mu_h_m = linear_magnetic_permeability(incident.material())?;
        let transmitted_mu_h_m = linear_magnetic_permeability(transmitted.material())?;
        for source in [incident, transmitted] {
            match source {
                DielectricResponse::Lossless(material) => {
                    linear_electric_permittivity(material)?;
                }
                DielectricResponse::Ohmic(conductor) => {
                    linear_electric_permittivity(conductor.material())?;
                    ohmic_conductivity(conductor)?;
                }
                DielectricResponse::DebyeWithOhmic {
                    dielectric,
                    conductor,
                } => {
                    shared_ohmic_conductivity(dielectric, conductor)?;
                }
                DielectricResponse::LorentzWithOhmic {
                    dielectric,
                    conductor,
                } => {
                    shared_lorentz_ohmic_conductivity(dielectric, conductor)?;
                }
                DielectricResponse::Debye(_)
                | DielectricResponse::Lorentz(_)
                | DielectricResponse::OpticalConductor(_)
                | DielectricResponse::OpticalMaterial(_)
                | DielectricResponse::SpectralOptical(_)
                | DielectricResponse::WaterOptical(_)
                | DielectricResponse::SilicaOptical(_)
                | DielectricResponse::ThermoOpticOptical(_)
                | DielectricResponse::MolecularGasOptical(_)
                | DielectricResponse::IonizedGas(_) => {}
            }
        }
        Ok(Self {
            incident: incident.material(),
            incident_response: incident,
            transmitted,
            incident_mu_h_m,
            transmitted_mu_h_m,
        })
    }

    pub const fn incident_material(&self) -> &'a ResolvedMaterialStatePoint {
        self.incident
    }
    /// Original incident law; its spectral parameterization is retained.
    pub const fn incident_response(&self) -> DielectricResponse<'a> {
        self.incident_response
    }
    /// Full transmitted constitutive owner, including separately sourced
    /// plasma carriers when present.
    pub const fn transmitted_response(&self) -> DielectricResponse<'a> {
        self.transmitted
    }
    /// Original scalar-property bundle; plasma returns its explicit background.
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
            return Err(ReflectionError::Invalid(
                "incidence cosine must lie in (0,1]",
            ));
        }
        let incident_epsilon_f_m = self.incident_response.lossless_permittivity(omega)?;
        let epsilon = self.transmitted.permittivity(omega)?;
        let medium = WaveMedium::new(
            epsilon,
            self.transmitted_mu_h_m,
            incident_epsilon_f_m,
            self.incident_mu_h_m,
        )?;
        let wave_number =
            positive(omega * det::sqrt(incident_epsilon_f_m) * det::sqrt(self.incident_mu_h_m))?;
        interface_response(medium, wave_number, c, polarization)
    }
}

fn interface_response(
    medium: WaveMedium,
    wave_number: f64,
    cosine: f64,
    polarization: Polarization,
) -> Result<InterfaceReflection, ReflectionError> {
    let normal = medium.normal(cosine)?;
    let mut response =
        interface_response_from_normal(medium, wave_number, cosine, polarization, normal)?;
    if response.reflected_power_fraction < 64.0 * f64::EPSILON {
        // Near Brewster or impedance matching, subtracting the rounded
        // square root from the other boundary coefficient can erase a finite
        // reflected field, including its sign/phase. Rationalize with the
        // original constitutive and angle data only in this cancellation
        // neighbourhood. The ordinary reflection path is unchanged.
        let reflected = near_zero_reflection(medium, cosine, polarization, normal)?;
        response.reflected_tangential_e = reflected;
        response.reflected_power_fraction = reflected.norm_sq();
    }
    Ok(response)
}

fn near_zero_reflection(
    medium: WaveMedium,
    cosine: f64,
    polarization: Polarization,
    normal: C64,
) -> Result<C64, ReflectionError> {
    if (medium.epsilon_ratio == C64::ONE && medium.mu_ratio == 1.0)
        || (cosine == 1.0 && medium.epsilon_ratio == C64::from_re(medium.mu_ratio))
    {
        // Identical media are transparent at every angle. Equal impedances
        // also match at normal incidence even when refractive indices differ.
        return Ok(C64::ZERO);
    }
    let e = Dd::from_f64(medium.epsilon_ratio.re);
    let f = Dd::from_f64(medium.epsilon_ratio.im);
    let m = Dd::from_f64(medium.mu_ratio);
    let c = Dd::from_f64(cosine);
    let contrast = e * m - Dd::ONE;
    // r=(a²-b²)/(a+b)², with the same tangential-E convention as the
    // native Fresnel pair. These factorizations also keep identical media
    // exactly reflection-free at arbitrarily small incidence cosines.
    let (real, imaginary, a, b) = match polarization {
        Polarization::Te => (
            ((m - Dd::ONE) * c) * ((m + Dd::ONE) * c) - contrast,
            -(f * m),
            C64::from_re(medium.mu_ratio * cosine),
            normal,
        ),
        Polarization::Tm => (
            contrast - ((e - Dd::ONE) * c) * ((e + Dd::ONE) * c) + (f * c) * (f * c),
            f * m - ((e * c) * (f * c)) * Dd::from_f64(2.0),
            normal,
            medium.epsilon_ratio.scale(cosine),
        ),
    };
    let scale = Dd::from_f64(a.re.abs().max(a.im.abs()).max(b.re.abs()).max(b.im.abs()));
    let numerator = C64::new(
        ((real / scale) / scale).to_f64(),
        ((imaginary / scale) / scale).to_f64(),
    );
    let (_, _, denominator, _) = normalized_interface_pair(a, b)?;
    let reflected = (numerator / denominator) / denominator;
    if !finite(reflected) {
        return Err(ReflectionError::Invalid(
            "near-zero reflected field is unrepresentable",
        ));
    }
    Ok(reflected)
}

// Also used by a film expressed in a real numerical reference port. Its
// resolved normal wave number may come from an explicit conserved kx.
fn interface_response_from_normal(
    medium: WaveMedium,
    wave_number: f64,
    cosine: f64,
    polarization: Polarization,
    normal: C64,
) -> Result<InterfaceReflection, ReflectionError> {
    let incident_normal = positive(wave_number * cosine)?;
    let transmitted_normal = normal.scale(wave_number);
    if !finite(transmitted_normal) {
        return Err(ReflectionError::Invalid(
            "outgoing wave number is unrepresentable",
        ));
    }

    // Multiply through by k_z for TM: no singular admittance evaluation
    // at the critical angle. a/b = Y_i/Y_t for both polarizations.
    let (a, b) = match polarization {
        Polarization::Te => (C64::from_re(medium.mu_ratio * cosine), normal),
        Polarization::Tm => (normal, medium.epsilon_ratio.scale(cosine)),
    };
    let (a, b, denominator, norm) = normalized_interface_pair(a, b)?;
    let reflected = (a - b) / denominator;
    let transmitted = a.scale(2.0) / denominator;
    let transmitted_h = b.scale(2.0) / denominator;
    let reflected_power = reflected.norm_sq();
    // Stable Poynting flux, not |t|^2 alone and not a clipped 1-|r|^2.
    let transmitted_power = 4.0 * (a.re * b.re + a.im * b.im) / norm;
    if !finite(reflected)
        || !finite(transmitted)
        || !finite(transmitted_h)
        || !transmitted_power.is_finite()
        || transmitted_power < 0.0
        || (reflected_power + transmitted_power - 1.0).abs() > 128.0 * f64::EPSILON
    {
        return Err(ReflectionError::Invalid(
            "passive interface power balance failed",
        ));
    }
    Ok(InterfaceReflection {
        reflected_tangential_e: reflected,
        transmitted_tangential_e: transmitted,
        transmitted_conjugate_h: transmitted_h,
        reflected_power_fraction: reflected_power,
        transmitted_power_fraction: transmitted_power,
        incident_normal_wavenumber_rad_m: incident_normal,
        transmitted_normal_wavenumber_rad_m: transmitted_normal,
        transmitted_permittivity_f_m: medium.epsilon_f_m,
    })
}

#[cfg(test)]
#[path = "plane_wave/reflection_tests.rs"]
mod reflection_tests;
