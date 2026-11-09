//! Analytic bulk dispersion of the original passive electric material laws.
//!
//! The outgoing branch is `exp(+i omega t - i k z)`, with `Im(k) <= 0`.
//! The native constitutive owner supplies epsilon; original Debye/Lorentz
//! poles or finite-band optical laws supply analytic derivatives. No frequency
//! differencing, pole fitting, material query, or allocation occurs. Ohmic
//! leakage is differentiated once; optical intensity absorption keeps its
//! original wavelength parameterization.
//!
//! `Re(dk/domega)` is group slowness and `Re(d²k/domega²)` is GVD. For an
//! absorbing/anomalous band, the inverse phase slope is NOT a signal-front or
//! energy velocity, and may be negative or exceed the vacuum light speed.
//! Complex derivatives also retain attenuation slope and curvature. A local
//! frequency jet alone neither propagates a broadband pulse nor bounds its
//! higher-order distortion. Native angular-spectrum exposure already sums
//! original frequency-resolved fields when full spectral propagation is needed.
//!
//! Definitions and limits: F. X. Kärtner, MIT 6.977, Ultrafast Optics, §2.6,
//! <https://ocw.mit.edu/courses/6-977-ultrafast-optics-spring-2005/d99af77c9cdca438616ee745ba41e548_chapter2.pdf>.

use super::{
    DielectricResponse, ReflectionError, finite, linear_magnetic_permeability,
    ohmic_conductivity, positive, shared_lorentz_ohmic_conductivity,
    shared_ohmic_conductivity,
};
use fs_math::{c64::C64, dd::Dd};

mod interface;
pub use interface::{
    FresnelAmplitudeDispersion, FresnelLogarithmicDispersion, IncidenceFrequencyConstraint,
    InterfaceFrequencyDispersion,
};

type Result<T> = core::result::Result<T, ReflectionError>;

/// Frequency-local propagation data with the complete original source attached.
/// Derivatives hold the material snapshot and its scalar permeability fixed.
/// This is a homogeneous bulk wave: no transverse waveguide constraint,
/// interface delay, diffraction, or spatially varying medium is inferred.
#[derive(Clone, Copy, Debug)]
pub struct BulkWaveDispersion<'a> {
    source: DielectricResponse<'a>,
    omega: f64,
    epsilon: C64,
    epsilon_derivatives: [C64; 2],
    permeability: f64,
    wave_number: C64,
    first: C64,
    second: C64,
}

/// Exact carrier transfer plus local derivatives of that same propagation.
/// These are field-amplitude quantities; power transmission is `|transfer|²`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BulkPropagation {
    /// `exp(-i k distance)` from the existing native propagation owner.
    pub transfer: C64,
    /// `Re(k) distance / omega` [s]; distinct from envelope delay.
    pub phase_delay_s: f64,
    /// `Re(dk/domega) distance` [s]; signed in anomalous/absorbing media.
    pub group_delay_s: f64,
    /// `Re(d²k/domega²) distance` [s²], the group-delay dispersion.
    pub group_delay_dispersion_s2: f64,
    /// `d ln|transfer| / domega = Im(dk/domega) distance` [s].
    pub log_amplitude_derivative_s: f64,
    /// `d² ln|transfer| / domega² = Im(d²k/domega²) distance` [s²].
    pub log_amplitude_second_derivative_s2: f64,
}

impl<'a> DielectricResponse<'a> {
    /// Prepare phase, group and second-order dispersion at a positive angular
    /// frequency [rad/s], using the original passive nondispersive, Debye,
    /// Lorentz, Ohmic, collisional plasma and finite-band optical laws.
    /// Source/material semantics match native interfaces. Plasma derivatives
    /// retain the original frozen carrier spectrum, including mutual drag;
    /// the cached native frequency and response must already be admitted.
    /// Optical derivatives retain Cauchy, Sellmeier, linear wavelength tables
    /// or the original IAPWS water and CHARMS silica laws; endpoints and
    /// nonsmooth knots refuse.
    /// Frequency-prepared gases retain their original derivatives at the SAME
    /// admitted interior frequency: fixed partial densities for molecular
    /// Lorentz--Lorenz mixtures, or fixed T/P/composition for Ciddor air.
    pub fn bulk_dispersion(self, omega_rad_s: f64) -> Result<BulkWaveDispersion<'a>> {
        let omega = positive(omega_rad_s)?;
        let optical = match self {
            Self::OpticalConductor(source) => Some(source.frequency_jet(omega)?),
            Self::OpticalMaterial(source) => Some(source.frequency_jet(omega)?),
            Self::SpectralOptical(source) => Some(source.frequency_jet(omega)?),
            Self::WaterOptical(source) => Some(source.frequency_jet(omega)?),
            Self::SilicaOptical(source) => Some(source.frequency_jet(omega)?),
            Self::ThermoOpticOptical(source) => Some(source.frequency_jet(omega)?),
            Self::MolecularGasOptical(source) => Some(source.frequency_jet(omega)?),
            _ => None,
        };
        if let Some(jet) = optical {
            return optical_bulk(self, omega, jet);
        }
        let epsilon = self.permittivity(omega)?;
        let permeability = linear_magnetic_permeability(self.material())?;
        let (epsilon_derivatives, loss_moments) = derivatives(self, omega)?;
        // Use the SAME binary-scaled outgoing bulk root as native absorbing
        // incidence. This preserves real-negative/evanescent branch selection
        // and avoids forming omega² or the dimensional epsilon*mu product.
        let wave_number = super::absorbing_incidence::kernel::Medium::new(
            omega, 0.0, epsilon, permeability,
        )?
        .normal;
        if wave_number == C64::ZERO {
            return Err(invalid("zero bulk wave number has no admitted frequency jet"));
        }
        let [slope, curvature] = epsilon_derivatives.map(|value| value / epsilon);
        let slowness = C64::new(wave_number.re / omega, wave_number.im / omega);
        // k = omega sqrt(mu epsilon). Store omega*epsilon' and
        // omega²*epsilon'' until the last dimensional divisions. The exact
        // constant-law limit gives second == 0 without subtracting k/omega².
        let mut first = slowness * (C64::ONE + slope.scale(0.5));
        let second_scaled =
            slowness * (slope + curvature.scale(0.5) - slope * slope.scale(0.25));
        let mut second = C64::new(second_scaled.re / omega, second_scaled.im / omega);
        if epsilon.re > 0.0 && epsilon.im < 0.0 && slowness.re > 0.0 && slowness.im < 0.0 {
            // Weak absorption makes Im(h) + Im(omega*h') cancel to much
            // smaller order than either term. Differentiate 2*Re(h)*Im(h)
            // = mu*Im(epsilon), using original per-law loss moments:
            // m1 = Im(epsilon + omega*epsilon'),
            // m2 = Im(2*omega*epsilon' + omega²*epsilon'').
            // Ohmic m1/m2 are exactly zero; Debye tails remain explicit.
            // A positive real epsilon keeps Re(h) the dominant root part.
            let lambda = (slowness * slope.scale(0.5)).re / slowness.re;
            let eta = (slowness
                * (curvature.scale(0.5) - slope * slope.scale(0.25))).re
                / slowness.re;
            let [loss_first, loss_second] = loss_moments.map(|value| value / epsilon.im);
            let attenuation_first = slowness.im * (loss_first - lambda);
            let attenuation_second = (slowness.im
                * (loss_second - 2.0 * lambda * loss_first + 2.0 * lambda * lambda - eta))
                / omega;
            // An extreme weak pole can have finite complex derivatives while
            // normalization by its tiny loss overflows. In that case retain
            // the already evaluated general jet, including evanescent limits.
            if [lambda, eta, loss_first, loss_second, attenuation_first, attenuation_second]
                .into_iter().all(f64::is_finite)
            {
                first.im = attenuation_first;
                second.im = attenuation_second;
            }
        }
        if ![slope, curvature, first, second].into_iter().all(finite) {
            return Err(invalid("bulk dispersion derivatives are unrepresentable"));
        }
        Ok(BulkWaveDispersion {
            source: self,
            omega,
            epsilon,
            epsilon_derivatives,
            permeability,
            wave_number,
            first,
            second,
        })
    }
}

impl<'a> BulkWaveDispersion<'a> {
    pub const fn source(self) -> DielectricResponse<'a> {
        self.source
    }
    pub const fn angular_frequency_rad_s(self) -> f64 {
        self.omega
    }
    pub const fn permittivity_f_m(self) -> C64 {
        self.epsilon
    }
    /// `[omega * epsilon', omega² * epsilon'']`, both in F/m. The second
    /// entry is a scaled ordinary second derivative, not a log derivative.
    pub const fn scaled_permittivity_derivatives_f_m(self) -> [C64; 2] {
        self.epsilon_derivatives
    }
    pub const fn permeability_h_m(self) -> f64 {
        self.permeability
    }
    pub const fn wavenumber_rad_m(self) -> C64 {
        self.wave_number
    }
    /// Complex `dk/domega` [s/m]. The imaginary part is loss dispersion.
    pub const fn wavenumber_derivative_s_m(self) -> C64 {
        self.first
    }
    /// Complex `d²k/domega²` [s²/m], including attenuation curvature.
    pub const fn wavenumber_second_derivative_s2_m(self) -> C64 {
        self.second
    }
    pub const fn group_slowness_s_m(self) -> f64 {
        self.first.re
    }
    pub const fn group_velocity_dispersion_s2_m(self) -> f64 {
        self.second.re
    }
    /// Phase-crest velocity, absent for a purely evanescent branch or an
    /// unrepresentable quotient. No group/energy/front interpretation.
    pub fn phase_velocity_m_s(self) -> Option<f64> {
        if self.wave_number.re == 0.0 {
            return None;
        }
        let speed = self.omega / self.wave_number.re;
        speed.is_finite().then_some(speed)
    }
    /// Inverse real phase slope. A vanishing slope returns None; negative
    /// and superluminal values remain signed, without interpreting them as
    /// information transport. Attenuation can reshape a pulse in such bands.
    pub fn group_velocity_m_s(self) -> Option<f64> {
        let speed = self.first.re.recip();
        speed.is_finite().then_some(speed)
    }
    /// Propagate the carrier through a nonnegative physical distance [m] with
    /// the existing bounded complex exponential and double-double phase
    /// reduction. Delays and GDD describe this original branch locally, not
    /// an invented constant-index or second-order replacement medium.
    pub fn propagation(self, distance_m: f64) -> Result<BulkPropagation> {
        if !distance_m.is_finite() || distance_m < 0.0 {
            return Err(invalid("bulk propagation requires finite nonnegative distance"));
        }
        let transfer = super::angular_spectrum::modes::propagation(self.wave_number, distance_m)
            .ok_or_else(|| invalid("bulk propagation is unrepresentable"))?;
        let phase_delay_s = (self.wave_number.re / self.omega) * distance_m;
        let group_delay_s = self.first.re * distance_m;
        let group_delay_dispersion_s2 = self.second.re * distance_m;
        let log_amplitude_derivative_s = self.first.im * distance_m;
        let log_amplitude_second_derivative_s2 = self.second.im * distance_m;
        if ![
            phase_delay_s, group_delay_s, group_delay_dispersion_s2,
            log_amplitude_derivative_s, log_amplitude_second_derivative_s2,
        ].iter().all(|value| value.is_finite()) {
            return Err(invalid("bulk propagation delay is unrepresentable"));
        }
        Ok(BulkPropagation {
            transfer,
            phase_delay_s,
            group_delay_s,
            group_delay_dispersion_s2,
            log_amplitude_derivative_s,
            log_amplitude_second_derivative_s2,
        })
    }
}

// Optical phase and attenuation are already the original vacuum-referenced N.
// Keep the native carrier branch, but differentiate k=omega*N/c0 directly:
// reconstructing optical group/GVD by squaring N and taking another root loses
// exact constant-alpha and linear-index cancellations.
fn optical_bulk<'a>(
    source: DielectricResponse<'a>,
    omega: f64,
    jet: super::optical::frequency::OpticalFrequencyJet,
) -> Result<BulkWaveDispersion<'a>> {
    let permeability = linear_magnetic_permeability(source.material())?;
    let wave_number = super::absorbing_incidence::kernel::Medium::new(
        omega, 0.0, jet.epsilon, permeability,
    )?.normal;
    let second = C64::new(
        jet.second_scaled_s_m.re / omega,
        jet.second_scaled_s_m.im / omega,
    );
    if wave_number == C64::ZERO || !finite(second) {
        return Err(invalid("optical bulk frequency derivatives are unrepresentable"));
    }
    Ok(BulkWaveDispersion {
        source,
        omega,
        epsilon: jet.epsilon,
        epsilon_derivatives: jet.epsilon_derivatives,
        permeability,
        wave_number,
        first: jet.first_s_m,
        second,
    })
}

fn invalid(message: &'static str) -> ReflectionError {
    ReflectionError::Invalid(message)
}

// Compensate opposing resonant slopes without allocating per-pole storage.
#[derive(Default)]
struct DerivativeSum {
    real: [Dd; 2],
    imaginary: [Dd; 2],
    loss_moments: [Dd; 2],
}
impl DerivativeSum {
    fn add(&mut self, values: [C64; 2], loss_moments: [f64; 2]) -> Result<()> {
        for (index, value) in values.into_iter().enumerate() {
            if !finite(value) || !loss_moments[index].is_finite() {
                return Err(invalid("material dispersion derivative is unrepresentable"));
            }
            self.real[index] = self.real[index] + Dd::from_f64(value.re);
            self.imaginary[index] = self.imaginary[index] + Dd::from_f64(value.im);
            self.loss_moments[index] =
                self.loss_moments[index] + Dd::from_f64(loss_moments[index]);
        }
        Ok(())
    }
    fn finish(self) -> Result<([C64; 2], [f64; 2])> {
        let values = core::array::from_fn(|index| {
            C64::new(self.real[index].to_f64(), self.imaginary[index].to_f64())
        });
        let loss_moments = self.loss_moments.map(Dd::to_f64);
        if values.into_iter().all(finite)
            && loss_moments.into_iter().all(f64::is_finite)
        {
            Ok((values, loss_moments))
        } else {
            Err(invalid("material dispersion derivative sum is unrepresentable"))
        }
    }
}

fn derivatives(
    source: DielectricResponse<'_>,
    omega: f64,
) -> Result<([C64; 2], [f64; 2])> {
    let mut sum = DerivativeSum::default();
    let (debye, lorentz, conductivity) = match source {
        DielectricResponse::IonizedGas(response) => {
            // bulk_dispersion already admits the original cached frequency and
            // native value. Differentiate the retained coupled carrier modes.
            let jet = response.transport()
                .permittivity_frequency_derivatives(omega, || false)?;
            return Ok((
                jet.scaled_permittivity_derivatives_f_m,
                jet.permittivity_loss_moments_f_m,
            ));
        }
        DielectricResponse::Lossless(_) => (None, None, 0.0),
        DielectricResponse::Debye(model) => (Some(model), None, 0.0),
        DielectricResponse::Lorentz(model) => (None, Some(model), 0.0),
        DielectricResponse::Ohmic(conductor) => (None, None, ohmic_conductivity(conductor)?),
        DielectricResponse::DebyeWithOhmic { dielectric, conductor } =>
            (Some(dielectric), None, shared_ohmic_conductivity(dielectric, conductor)?),
        DielectricResponse::LorentzWithOhmic { dielectric, conductor } =>
            (None, Some(dielectric), shared_lorentz_ohmic_conductivity(dielectric, conductor)?),
        _ => return Err(invalid("this original wave law has no admitted analytic bulk frequency jet")),
    };
    if let Some(model) = debye {
        for pole in model.poles() {
            let x = omega * pole.relaxation_time_s;
            // h=1/(1+ix), q=ix/(1+ix). Form both directly so a weak
            // relaxation is not lost by subtracting a rounded h from one.
            let (h, q) = if x <= 1.0 {
                let denominator = 1.0 + x * x;
                (C64::new(1.0 / denominator, -x / denominator),
                 C64::new((x * x) / denominator, x / denominator))
            } else {
                let inverse = (1.0 / omega) / pole.relaxation_time_s;
                let denominator = 1.0 + inverse * inverse;
                (C64::new((inverse * inverse) / denominator, -inverse / denominator),
                 C64::new(1.0 / denominator, inverse / denominator))
            };
            let first = -(h.scale(pole.permittivity_increment_f_m) * q);
            let loss = h.im * pole.permittivity_increment_f_m;
            let loss_first = (2.0 * loss) * h.re;
            sum.add(
                [first, (first * q).scale(-2.0)],
                [loss_first, loss_first * (4.0 * h.re - 3.0)],
            )?;
        }
    }
    if let Some(model) = lorentz {
        for &pole in model.poles() {
            let susceptibility = pole.susceptibility_f_m(omega)
                .map_err(fs_material::field_storage::FieldStorageError::from)?;
            let natural = pole.angular_frequency_rad_s();
            let damping = pole.damping_rate_per_s();
            let scale = natural.max(omega).max(damping);
            let x = omega / scale;
            let denominator = C64::new(
                ((natural - omega) / scale) * (natural / scale + x),
                (damping / scale) * x,
            );
            // p=omega D'/D and r=-omega² D''/D, D=omega0²-omega²+i gamma omega.
            let p = C64::new(-2.0 * x * x, (damping / scale) * x) / denominator;
            let r = C64::from_re(2.0 * x * x) / denominator;
            let first = -(susceptibility * p);
            let second = (first * p).scale(-2.0) + susceptibility * r;
            // Im(chi + omega*chi') = 2*Im(chi)*Re((omega0²+omega²)/D).
            // Keep the original real-zero at resonance instead of subtracting
            // two almost equal loss contributions.
            let natural_scaled = natural / scale;
            let loss_factor = (C64::from_re(natural_scaled * natural_scaled + x * x)
                / denominator).re;
            sum.add(
                [first, second],
                [susceptibility.im * (2.0 * loss_factor), 2.0 * first.im + second.im],
            )?;
        }
    }
    if conductivity != 0.0 {
        let loss = positive(conductivity / omega)?;
        // The loss of omega*epsilon is constant for the original Ohmic law.
        sum.add([C64::new(0.0, loss), C64::new(0.0, -2.0 * loss)], [0.0; 2])?;
    }
    sum.finish()
}
