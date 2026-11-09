//! Analytic frequency derivatives of the original Fresnel boundary solve.
//!
//! With exp(+i omega t), a transfer H has delay -Im(H'/H). Its logarithmic
//! curvature is H''/H-(H'/H)^2. These are local narrowband phase slopes, not
//! signal-front or energy velocities; a zero amplitude has no phase or delay.
//! The authored transverse-frequency constraint is part of the result.
//!
//! Fresnel fields use the reference and tangential-axis convention of the
//! parent plane-wave owner. Phase/group-delay definitions and their limits:
//! <https://ocw.mit.edu/courses/6-977-ultrafast-optics-spring-2005/d99af77c9cdca438616ee745ba41e548_chapter2.pdf>.

use super::{BulkWaveDispersion, Result, invalid};
use crate::plane_wave::{
    InterfaceReflection, PlanarDielectricInterface, Polarization, WaveMedium, finite,
    normalized_interface_pair,
};
use fs_math::c64::C64;

/// What stays fixed as frequency changes near the admitted carrier.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IncidenceFrequencyConstraint {
    /// One incident direction: k_parallel changes with the original incident
    /// material's dispersive bulk wavenumber. A nondispersive interface then
    /// has no frequency-dependent Fresnel phase, including total reflection.
    FixedAngle,
    /// One transverse spatial Fourier component: k_parallel stays fixed and
    /// the incidence angle changes with frequency. This is the path used by
    /// a fixed transverse mode of the original angular-spectrum field owner.
    FixedTransverseWavenumber,
}

/// Original complex carrier amplitude and its analytic frequency derivatives.
/// Scaled derivatives avoid unnecessary powers of optical angular frequency.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FresnelAmplitudeDispersion {
    amplitude: C64,
    omega: f64,
    scaled_derivatives: [C64; 2],
}

/// Local logarithmic amplitude and phase dispersion of a nonzero transfer.
/// Signed delays in absorbing or anomalous bands do not imply signal speed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FresnelLogarithmicDispersion {
    /// -d arg(H)/domega [s], at the physical interface z=0.
    pub group_delay_s: f64,
    /// -d² arg(H)/domega² [s²].
    pub group_delay_dispersion_s2: f64,
    /// d ln|H|/domega [s].
    pub log_amplitude_derivative_s: f64,
    /// d² ln|H|/domega² [s²].
    pub log_amplitude_second_derivative_s2: f64,
}

impl FresnelAmplitudeDispersion {
    pub const fn amplitude(self) -> C64 {
        self.amplitude
    }
    /// [omega H', omega² H'']; the second entry is an ordinary derivative.
    pub const fn scaled_derivatives(self) -> [C64; 2] {
        self.scaled_derivatives
    }
    pub fn derivative_s(self) -> C64 {
        divide_real(self.scaled_derivatives[0], self.omega)
    }
    pub fn second_derivative_s2(self) -> C64 {
        divide_real(
            divide_real(self.scaled_derivatives[1], self.omega),
            self.omega,
        )
    }
    /// None at an exact Fresnel/Brewster zero or an unrepresentable logarithmic
    /// quotient. Complex amplitude derivatives remain available at a zero;
    /// assigning it an arbitrary phase or zero delay would be incorrect.
    pub fn logarithmic_dispersion(self) -> Option<FresnelLogarithmicDispersion> {
        if self.amplitude == C64::ZERO {
            return None;
        }
        let first = self.scaled_derivatives[0] / self.amplitude;
        let second = self.scaled_derivatives[1] / self.amplitude - first * first;
        let first = divide_real(first, self.omega);
        let second = divide_real(divide_real(second, self.omega), self.omega);
        if !finite(first) || !finite(second) {
            return None;
        }
        Some(FresnelLogarithmicDispersion {
            group_delay_s: -first.im,
            group_delay_dispersion_s2: -second.im,
            log_amplitude_derivative_s: first.re,
            log_amplitude_second_derivative_s2: second.re,
        })
    }
}

/// Source-bound Fresnel jets, using the unchanged native passive carrier solve.
#[derive(Clone, Copy, Debug)]
pub struct InterfaceFrequencyDispersion<'a> {
    incident: BulkWaveDispersion<'a>,
    transmitted: BulkWaveDispersion<'a>,
    cosine: f64,
    polarization: Polarization,
    constraint: IncidenceFrequencyConstraint,
    carrier: InterfaceReflection,
    amplitudes: [FresnelAmplitudeDispersion; 3],
}

impl<'a> InterfaceFrequencyDispersion<'a> {
    pub const fn incident_dispersion(self) -> BulkWaveDispersion<'a> {
        self.incident
    }
    pub const fn transmitted_dispersion(self) -> BulkWaveDispersion<'a> {
        self.transmitted
    }
    pub const fn incidence_cosine(self) -> f64 {
        self.cosine
    }
    pub const fn polarization(self) -> Polarization {
        self.polarization
    }
    pub const fn frequency_constraint(self) -> IncidenceFrequencyConstraint {
        self.constraint
    }
    pub const fn carrier(self) -> InterfaceReflection {
        self.carrier
    }
    pub const fn reflected_tangential_e(self) -> FresnelAmplitudeDispersion {
        self.amplitudes[0]
    }
    pub const fn transmitted_tangential_e(self) -> FresnelAmplitudeDispersion {
        self.amplitudes[1]
    }
    pub const fn transmitted_conjugate_h(self) -> FresnelAmplitudeDispersion {
        self.amplitudes[2]
    }
}

impl<'a> PlanarDielectricInterface<'a> {
    /// Differentiate the native TE/TM boundary amplitudes at a positive carrier
    /// frequency, retaining the original material jets and explicitly chosen
    /// direction constraint. Preparation and queries allocate no storage.
    ///
    /// Incident epsilon and its frequency jet must be lossless. The transmitted
    /// source may absorb or be evanescent. Original finite-band/derivative-domain
    /// refusals propagate; exact critical transmission has a singular square-root
    /// jet and refuses here, while the carrier-only `response` remains available.
    /// No finite-band pulse, finite-film delay or additional material law is fitted.
    pub fn frequency_dispersion(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        polarization: Polarization,
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<InterfaceFrequencyDispersion<'a>> {
        let carrier = self.response(omega_rad_s, incidence_cosine, polarization)?;
        let incident = self.incident_response.bulk_dispersion(omega_rad_s)?;
        let transmitted = self.transmitted.bulk_dispersion(omega_rad_s)?;
        let epsilon = incident.permittivity_f_m().re;
        let incident_derivatives = incident.scaled_permittivity_derivatives_f_m();
        if incident_derivatives.iter().any(|value| value.im != 0.0) {
            return Err(invalid("incident frequency jet must remain lossless"));
        }
        let medium = WaveMedium::new(
            carrier.transmitted_permittivity_f_m,
            transmitted.permeability_h_m(),
            epsilon,
            incident.permeability_h_m(),
        )?;
        let ei = incident_derivatives.map(|value| divide_real(value, epsilon));
        let et = transmitted
            .scaled_permittivity_derivatives_f_m()
            .map(|value| divide_real(value, epsilon));
        let derivatives =
            boundary_derivatives(medium, incidence_cosine, polarization, constraint, ei, et)?;
        let values = [
            carrier.reflected_tangential_e,
            carrier.transmitted_tangential_e,
            carrier.transmitted_conjugate_h,
        ];
        let amplitudes = core::array::from_fn(|index| FresnelAmplitudeDispersion {
            amplitude: values[index],
            omega: omega_rad_s,
            scaled_derivatives: if index == 2 {
                derivatives.map(core::ops::Neg::neg)
            } else {
                derivatives
            },
        });
        if amplitudes
            .iter()
            .any(|jet| !finite(jet.derivative_s()) || !finite(jet.second_derivative_s2()))
        {
            return Err(invalid("Fresnel frequency derivatives are unrepresentable"));
        }
        Ok(InterfaceFrequencyDispersion {
            incident,
            transmitted,
            cosine: incidence_cosine,
            polarization,
            constraint,
            carrier,
            amplitudes,
        })
    }
}

fn divide_real(value: C64, denominator: f64) -> C64 {
    C64::new(value.re / denominator, value.im / denominator)
}

// Derivatives of k_z/omega remove the common omega factor before the boundary
// quotient. Constant material/angle laws therefore have exactly zero jets,
// without subtracting almost equal propagation terms. All derivatives below
// are scaled by the carrier omega (first) or omega² (ordinary second).
fn boundary_derivatives(
    medium: WaveMedium,
    cosine: f64,
    polarization: Polarization,
    constraint: IncidenceFrequencyConstraint,
    ei: [C64; 2],
    et: [C64; 2],
) -> Result<[C64; 2]> {
    let qi = C64::from_re(cosine);
    let qt = medium.normal(cosine)?;
    if qt == C64::ZERO {
        return Err(invalid("critical transmission has no finite frequency jet"));
    }
    let sine_squared = (1.0 - cosine) * (1.0 + cosine);
    let (di, dt) = match constraint {
        IncidenceFrequencyConstraint::FixedAngle => (
            [ei[0].scale(cosine * cosine), ei[1].scale(cosine * cosine)],
            core::array::from_fn(|index| {
                et[index].scale(medium.mu_ratio) - ei[index].scale(sine_squared)
            }),
        ),
        IncidenceFrequencyConstraint::FixedTransverseWavenumber => {
            // k_parallel²/omega² contributes +2 s² and -6 s².
            let transverse = [2.0 * sine_squared, -6.0 * sine_squared];
            (
                core::array::from_fn(|index| ei[index] + C64::from_re(transverse[index])),
                core::array::from_fn(|index| {
                    et[index].scale(medium.mu_ratio) + C64::from_re(transverse[index])
                }),
            )
        }
    };
    let di = root_derivatives(qi, di);
    let dt = root_derivatives(qt, dt);
    let (a, b, da, db) = match polarization {
        Polarization::Te => (
            qi.scale(medium.mu_ratio),
            qt,
            di.map(|value| value.scale(medium.mu_ratio)),
            dt,
        ),
        Polarization::Tm => (
            qt,
            medium.epsilon_ratio * qi,
            [
                ei[0] * qt + dt[0],
                ei[1] * qt + (ei[0] * dt[0]).scale(2.0) + dt[1],
            ],
            [
                et[0] * qi + medium.epsilon_ratio * di[0],
                et[1] * qi + (et[0] * di[0]).scale(2.0) + medium.epsilon_ratio * di[1],
            ],
        ),
    };
    let scale = a.re.abs().max(a.im.abs()).max(b.re.abs()).max(b.im.abs());
    let (a, b, denominator, _) = normalized_interface_pair(a, b)?;
    let da = da.map(|value| divide_real(value, scale));
    let db = db.map(|value| divide_real(value, scale));
    // Differentiate r=(a-b)/(a+b) by its cross product, so a Brewster zero
    // never divides by r and weak transmitted E/H are not lost by subtraction.
    let first = ((da[0] * b - a * db[0]).scale(2.0) / denominator) / denominator;
    let second = ((da[1] * b - a * db[1]).scale(2.0) / denominator
        - (first * (da[0] + db[0])).scale(2.0))
        / denominator;
    if !finite(first) || !finite(second) {
        return Err(invalid("Fresnel frequency derivatives are unrepresentable"));
    }
    // t_E=1+r and t_H=1-r have the same/opposite derivatives, but each
    // carrier value above still comes from its independent native numerator.
    Ok([first, second])
}

fn root_derivatives(root: C64, squared_derivatives: [C64; 2]) -> [C64; 2] {
    let first = (squared_derivatives[0] / root).scale(0.5);
    [
        first,
        (squared_derivatives[1].scale(0.5) - first * first) / root,
    ]
}
