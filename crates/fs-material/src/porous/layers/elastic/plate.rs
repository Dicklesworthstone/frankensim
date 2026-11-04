//! Exact P/S refraction through a finite, fluid-loaded isotropic elastic plate.
//!
//! Both plate faces enforce continuous normal velocity/traction and zero shear
//! traction. The finite plate retains coherent longitudinal/SV conversion,
//! multiple reflection, Lamb resonances and evanescent tunneling. It is the
//! original frozen, unstressed, lossless elastic state: no fitted plate theory,
//! artificial attenuation, anisotropic projection or new material authority.
//! Existing equivalent-fluid entrance layers keep their original preparation.
//!
//! Even/odd elastic potentials give two real pressure/velocity boundary pairs.
//! Each normal-wave sine/cosine is scaled by exp(-|Im(k_n)| d/2); their sinc
//! limits retain exact critical refraction without dividing by a zero cosine.
//! The receiving fluid load is eliminated projectively, including pressure or
//! velocity zeros. A factored Wronskian keeps exponentially small transmitted
//! fields without subtracting almost equal growing transfer-matrix entries.
//!
//! Exact fluid/solid boundary conditions and the potential model are described
//! in Almeida et al., *Frustrated total internal reflection of ultrasonic waves
//! at a fluid-coupled elastic plate*, <https://arxiv.org/abs/2302.08826>, and
//! Johansson et al., JASA 118, 3436 (2005), <https://doi.org/10.1121/1.2126927>.
//! The exp(+i omega t), outgoing exp(-i k_n z) convention matches the parent.

use super::super::transmission::{log_scaled, positive_product_ratio, transmit_state};
use super::super::{
    Layer, PreparedLayer, RigidBackedLayers, bad, divided, finite, magnitude_scale,
    normal_wavenumber_ratio, retained_fluid_snapshot,
};
use super::ElasticTransmittingLayers;
use crate::MaterialError;
use crate::elastic::acoustic::ElasticAcousticState;
use crate::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::{c64::C64, dd::Dd, det};

/// Original equivalent-fluid entrance layers, one finite isotropic elastic
/// plate, and an outgoing lossless fluid. Material states and Christoffel
/// speeds are evaluated once; incidence queries allocate nothing.
pub struct ElasticPlateLayers<'a> {
    solid: ElasticTransmittingLayers<'a>,
    receiver_source: &'a dyn FluidAcoustics,
    receiver: FluidAcousticState,
    longitudinal_fluid: FluidAcousticState,
    thickness: f64,
    speed_ratio: f64,
    density_ratio: f64,
}

/// Pressure R at the entrance, pressure T at the receiving face, and normal
/// power fractions. P/S wavenumbers retain the actual refracted branches, even
/// when both modes are evanescent inside a transmitting finite plate.
///
/// The solid and receiving fluid are lossless. Absorption belongs only to the
/// finite porous entrance layers. Evanescent fields are never counted as heat;
/// a propagating receiver can carry power tunneled across an evanescent plate.
#[derive(Clone, Copy, Debug)]
pub struct LayeredElasticPlateTransmission {
    pub pressure_reflection: C64,
    pub pressure_transmission: C64,
    pub plate_front_pressure: C64,
    pub longitudinal_normal_wavenumber_rad_m: C64,
    pub shear_normal_wavenumber_rad_m: C64,
    pub receiving_normal_wavenumber_rad_m: C64,
    pub reflected_fraction: f64,
    pub transmitted_fraction: f64,
    pub absorbed_fraction: f64,
}

impl<'a> ElasticPlateLayers<'a> {
    /// The original solid must be explicitly isotropic, as for the half-space
    /// owner. Thickness is finite and nonnegative. Zero thickness is the exact
    /// fluid-interface limit, while material admission remains source-backed.
    pub fn new(
        incident: &'a dyn ThermoviscousFluid,
        plate: &'a ElasticAcousticState,
        receiver_source: &'a dyn FluidAcoustics,
        omega: f64,
        thickness_m: f64,
        layers: &[Layer<'a>],
    ) -> Result<Self, MaterialError> {
        if !thickness_m.is_finite() || thickness_m < 0.0 {
            return Err(bad(
                "elastic plate thickness must be finite and nonnegative",
            ));
        }
        let solid = ElasticTransmittingLayers::new(incident, plate, omega, layers)?;
        let prepared = &solid.prepared;
        let incident_state = prepared.fluid().acoustic();
        let original_incident: &dyn FluidAcoustics = incident;
        let receiver = match retained_fluid_snapshot(
            receiver_source,
            original_incident,
            incident_state,
            &prepared.layers,
        ) {
            Some(state) => state,
            None => receiver_source.acoustic_state()?,
        };
        let speed_ratio = incident_state.sound_speed_m_s() / receiver.sound_speed_m_s();
        let density_ratio = incident_state.density_kg_m3() / receiver.density_kg_m3();
        if [speed_ratio, density_ratio]
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(bad("elastic plate receiving contrast is unrepresentable"));
        }
        // The exact normal-incidence P limit uses the existing finite-fluid
        // kernel with the original solid's rho and c_P. It introduces no new
        // source query or oblique fluid approximation to the elastic plate.
        let longitudinal_fluid =
            FluidAcousticState::try_new(plate.density_kg_m3(), solid.longitudinal_speed)?;
        Ok(Self {
            solid,
            receiver_source,
            receiver,
            longitudinal_fluid,
            thickness: thickness_m,
            speed_ratio,
            density_ratio,
        })
    }

    pub const fn prepared_layers(&self) -> &RigidBackedLayers<'a> {
        &self.solid.prepared
    }

    pub const fn plate_state(&self) -> &'a ElasticAcousticState {
        self.solid.receiver
    }

    /// Original longitudinal and shear Christoffel speeds, in that order.
    pub const fn plate_phase_speeds_m_s(&self) -> [f64; 2] {
        [self.solid.longitudinal_speed, self.solid.shear_speed]
    }

    pub const fn thickness_m(&self) -> f64 {
        self.thickness
    }

    pub const fn receiving_source(&self) -> &'a dyn FluidAcoustics {
        self.receiver_source
    }

    pub const fn receiving_state(&self) -> &FluidAcousticState {
        &self.receiver
    }

    /// Exact coherent plate scattering at a propagating incident-fluid angle.
    /// The P, SV and receiving-fluid critical angles all retain finite limits.
    /// No finite-thickness or frequency approximation is used.
    pub fn scattering(
        &self,
        cosine: f64,
    ) -> Result<LayeredElasticPlateTransmission, MaterialError> {
        if !cosine.is_finite() || cosine <= 0.0 || cosine > 1.0 {
            return Err(bad("elastic plate needs incidence cosine in (0,1]"));
        }
        let prepared = &self.solid.prepared;
        let ci = prepared.fluid().acoustic().sound_speed_m_s();
        let qp = normal_wavenumber_ratio(
            self.solid.longitudinal_ratio,
            cosine,
            ci,
            self.solid.longitudinal_speed,
        );
        let qs =
            normal_wavenumber_ratio(self.solid.shear_ratio, cosine, ci, self.solid.shear_speed);
        let qr = normal_wavenumber_ratio(
            self.speed_ratio,
            cosine,
            ci,
            self.receiver.sound_speed_m_s(),
        );
        let p_kn = qp.scale(prepared.k0);
        let s_kn = qs.scale(prepared.k0);
        let r_kn = qr.scale(prepared.k0);
        for (q, kn) in [(qp, p_kn), (qs, s_kn), (qr, r_kn)] {
            if !finite(kn) || (q != C64::ZERO && kn == C64::ZERO) {
                return Err(bad("elastic plate normal wavenumber is unrepresentable"));
            }
        }
        let load = if qr == C64::ZERO {
            C64::ZERO
        } else {
            let value = positive_product_ratio(self.density_ratio, magnitude_scale(qr), cosine)?;
            if qr.re > 0.0 {
                C64::from_re(value)
            } else {
                C64::new(0.0, -value)
            }
        };
        let boundary = if self.thickness == 0.0 {
            PlateBoundary {
                pressure: C64::ONE,
                velocity: load,
                transmission: 1.0,
                transmission_log_scale: 0.0,
            }
        } else if cosine == 1.0 {
            self.normal_boundary(load)?
        } else {
            let cp = divided(qp, self.solid.longitudinal_ratio);
            let cs = divided(qs, self.solid.shear_ratio);
            let sine = det::sqrt((1.0 - cosine) * (1.0 + cosine));
            let ss = sine / self.solid.shear_ratio;
            let zeta = self.solid.impedance_ratio / cosine;
            if !finite(cp) || !finite(cs) || !ss.is_finite() || !zeta.is_finite() || zeta <= 0.0 {
                return Err(bad("elastic plate refraction contrast is unrepresentable"));
            }
            parity_boundary(
                bounded_wave(
                    cp,
                    p_kn,
                    prepared.omega / self.solid.longitudinal_speed,
                    self.thickness,
                )?,
                bounded_wave(
                    cs,
                    s_kn,
                    prepared.omega / self.solid.shear_speed,
                    self.thickness,
                )?,
                ss,
                self.solid.shear_speed / self.solid.longitudinal_speed,
                zeta,
                load,
            )?
        };
        let state = transmit_state(prepared, cosine, boundary.pressure, boundary.velocity)?;
        let transmission = C64::from_re(boundary.transmission);
        let pressure_transmission =
            state.scaled_coefficient(transmission, boundary.transmission_log_scale)?;
        let transmitted = if load.re == 0.0 {
            0.0
        } else {
            state.power(
                transmission,
                det::ln(load.re) + 2.0 * boundary.transmission_log_scale,
            )
        };
        let reflected = state.reflection.norm_sq();
        let absorbed = if self.solid.lossless {
            0.0
        } else {
            state.entering - transmitted
        };
        let tolerance = 2048.0 * f64::EPSILON;
        if !transmitted.is_finite()
            || transmitted < 0.0
            || !absorbed.is_finite()
            || absorbed < -tolerance
            || (reflected + transmitted + absorbed - 1.0).abs() > tolerance
        {
            return Err(bad("elastic plate failed passive normal-power balance"));
        }
        Ok(LayeredElasticPlateTransmission {
            pressure_reflection: state.reflection,
            pressure_transmission,
            plate_front_pressure: state.coefficient(boundary.pressure)?,
            longitudinal_normal_wavenumber_rad_m: p_kn,
            shear_normal_wavenumber_rad_m: s_kn,
            receiving_normal_wavenumber_rad_m: r_kn,
            reflected_fraction: reflected,
            transmitted_fraction: transmitted,
            absorbed_fraction: absorbed,
        })
    }

    fn normal_boundary(&self, load: C64) -> Result<PlateBoundary, MaterialError> {
        let initial_scale = magnitude_scale(load).max(1.0);
        let layer = PreparedLayer {
            source: Layer::FluidLayer {
                fluid: &self.longitudinal_fluid,
                thickness_m: self.thickness,
            },
            porous: None,
            fluid: Some(self.longitudinal_fluid),
        };
        let (pressure, velocity, normalization, _) = self.solid.prepared.propagate_layer(
            &layer,
            1.0,
            C64::from_re(1.0 / initial_scale),
            divided(load, initial_scale),
        )?;
        Ok(PlateBoundary {
            pressure,
            velocity,
            transmission: 1.0,
            transmission_log_scale: -det::ln(initial_scale) - log_scaled(16.0, normalization),
        })
    }
}

pub(super) struct PlateBoundary {
    pub(super) pressure: C64,
    pub(super) velocity: C64,
    pub(super) transmission: f64,
    pub(super) transmission_log_scale: f64,
}

// C=cos(k_n h), H=sin(k_n h)/cos(theta), L=cos(theta)*sin(k_n h),
// each multiplied by exp(-attenuation). All three are real on either outgoing
// lossless branch, including the critical cosine=0 limit H=k*h.
pub(super) struct BoundedWave {
    c: f64,
    h: f64,
    l: f64,
    attenuation: f64,
}

pub(super) fn bounded_wave(
    cosine: C64,
    kn: C64,
    k: f64,
    thickness: f64,
) -> Result<BoundedWave, MaterialError> {
    let half = thickness * 0.5;
    let kh = k * half;
    let wave = if kn.im == 0.0 {
        let phase = kn.re * half;
        if !phase.is_finite() {
            return Err(bad("elastic plate phase is unrepresentable"));
        }
        let (s, c) = det::sin_cos(phase);
        let h = if phase.abs() < 1.0e-4 {
            kh * if phase == 0.0 { 1.0 } else { s / phase }
        } else {
            s / cosine.re
        };
        BoundedWave {
            c,
            h,
            l: cosine.re * s,
            attenuation: 0.0,
        }
    } else {
        let attenuation = -kn.im * half;
        let difference = -det::expm1(-2.0 * attenuation);
        let h = if attenuation < 1.0e-4 {
            kh * if attenuation == 0.0 {
                1.0
            } else {
                difference / (2.0 * attenuation)
            }
        } else {
            (0.5 * difference) / -cosine.im
        };
        BoundedWave {
            c: 1.0 - 0.5 * difference,
            h,
            l: 0.5 * cosine.im * difference,
            attenuation,
        }
    };
    if [wave.c, wave.h, wave.l, wave.attenuation]
        .iter()
        .any(|value| !value.is_finite())
    {
        return Err(bad("elastic plate bounded wave is unrepresentable"));
    }
    Ok(wave)
}

fn product_sum(a: f64, b: f64, c: f64, d: f64) -> f64 {
    (Dd::from_f64(a) * Dd::from_f64(b) + Dd::from_f64(c) * Dd::from_f64(d)).to_f64()
}

fn parity_boundary(
    p: BoundedWave,
    s: BoundedWave,
    ss: f64,
    speed_ratio: f64,
    zeta: f64,
    load: C64,
) -> Result<PlateBoundary, MaterialError> {
    parity_boundary_pair(p, s, ss, speed_ratio, zeta, C64::ONE, load)
}

// The same shear-free parity reduction accepts any homogeneous downstream
// pressure/velocity pair. Neither a pressure/velocity zero nor a plate
// transmission zero is divided out when fluid-coupled plates are composed.
pub(super) fn parity_boundary_pair(
    p: BoundedWave,
    s: BoundedWave,
    ss: f64,
    speed_ratio: f64,
    zeta: f64,
    exit_pressure: C64,
    exit_velocity: C64,
) -> Result<PlateBoundary, MaterialError> {
    let a = (-2.0 * ss).mul_add(ss, 1.0);
    let aa = a * a;
    let b = 4.0 * speed_ratio * ss * ss;
    // Shear-free even channel: [e_p, i e_q] at the front, [e_p, -i e_q]
    // at the back. Odd channel: [i o_p, o_q] and [-i o_p, o_q].
    let ep = product_sum(aa, p.c * s.h, b, p.l * s.c);
    let eq = zeta * p.l * s.h;
    let op = product_sum(aa, s.c * p.h, b, s.l * p.c);
    let oq = zeta * p.c * s.c;
    let even_scale = ep.abs().max(eq.abs());
    let odd_scale = op.abs().max(oq.abs());
    if [ep, eq, op, oq].iter().any(|value| !value.is_finite())
        || even_scale == 0.0
        || odd_scale == 0.0
    {
        return Err(bad("elastic plate parity pair is unrepresentable"));
    }
    let (ep, eq) = (ep / even_scale, eq / even_scale);
    let (op, oq) = (op / odd_scale, oq / odd_scale);
    let difference = product_sum(ep, oq, -op, eq);
    let pressure_product = 2.0 * ep * op;
    let velocity_product = 2.0 * eq * oq;
    let load_scale = magnitude_scale(exit_pressure).max(magnitude_scale(exit_velocity));
    if !finite(exit_pressure) || !finite(exit_velocity) || load_scale == 0.0 {
        return Err(bad(
            "elastic plate receiving boundary pair is unrepresentable",
        ));
    }
    let scaled_load = divided(exit_velocity, load_scale);
    let (pressure, velocity) = if exit_pressure == C64::ONE {
        // Retain the original single-plate arithmetic exactly.
        (
            C64::from_re(difference / load_scale)
                + C64::new(-scaled_load.im, scaled_load.re).scale(pressure_product),
            scaled_load.scale(difference) + C64::new(0.0, velocity_product / load_scale),
        )
    } else {
        let scaled_pressure = divided(exit_pressure, load_scale);
        // A downstream cavity may cancel the two complex contributions at a
        // resonance. Preserve the product tails before rounding each sum.
        (
            C64::new(
                product_sum(
                    difference,
                    scaled_pressure.re,
                    -pressure_product,
                    scaled_load.im,
                ),
                product_sum(
                    difference,
                    scaled_pressure.im,
                    pressure_product,
                    scaled_load.re,
                ),
            ),
            C64::new(
                product_sum(
                    difference,
                    scaled_load.re,
                    -velocity_product,
                    scaled_pressure.im,
                ),
                product_sum(
                    difference,
                    scaled_load.im,
                    velocity_product,
                    scaled_pressure.re,
                ),
            ),
        )
    };

    let direct = product_sum(ep, oq, op, eq);
    let direct_scale = (ep * oq).abs() + (op * eq).abs();
    let (transmission, transmission_log_scale) = if direct.abs() > 1.0e-4 * direct_scale {
        // Well-conditioned determinant from the same normalized boundary
        // pairs keeps lossless flux correlated at narrow plate resonances.
        (direct / load_scale, 0.0)
    } else {
        // W = zeta*(a² Hs Cs exp(-2 Ap) + b Lp Cp exp(-2 As)).
        // Removing the common attenuation BEFORE the sum retains tunneling
        // after the two direct determinant products become equal in f64.
        let attenuation = p.attenuation.min(s.attenuation);
        let first = aa * s.h * s.c * det::exp(-2.0 * (p.attenuation - attenuation));
        let second = b * p.l * p.c * det::exp(-2.0 * (s.attenuation - attenuation));
        let log_scale = det::ln(zeta)
            - 2.0 * attenuation
            - det::ln(even_scale)
            - det::ln(odd_scale)
            - det::ln(load_scale);
        (first + second, log_scale)
    };
    if !transmission.is_finite() || !transmission_log_scale.is_finite() {
        return Err(bad(
            "elastic plate transmission determinant is unrepresentable",
        ));
    }
    Ok(PlateBoundary {
        pressure,
        velocity,
        transmission,
        transmission_log_scale,
    })
}

#[cfg(test)]
mod pair_tests {
    use super::*;

    #[test]
    fn pressure_and_velocity_nodes_keep_projective_linearity() {
        let solve = |pressure, velocity| {
            let cp = (1.0_f64 - (0.4_f64 / 0.6).powi(2)).sqrt();
            let cs = (1.0_f64 - 0.4_f64.powi(2)).sqrt();
            parity_boundary_pair(
                bounded_wave(C64::from_re(cp), C64::from_re(3.0 * cp), 3.0, 0.7).unwrap(),
                bounded_wave(C64::from_re(cs), C64::from_re(5.0 * cs), 5.0, 0.7).unwrap(),
                0.4,
                0.6,
                0.2,
                pressure,
                velocity,
            )
            .unwrap()
        };
        let physical = |pair: PlateBoundary| {
            let scale = det::exp(-pair.transmission_log_scale) / pair.transmission;
            [pair.pressure.scale(scale), pair.velocity.scale(scale)]
        };
        let factor = C64::new(-0.7, 1.3);
        for (pressure, velocity) in [(C64::ZERO, C64::ONE), (C64::ONE, C64::ZERO)] {
            let original = physical(solve(pressure, velocity));
            let scaled = physical(solve(pressure * factor, velocity * factor));
            for component in 0..2 {
                assert!((scaled[component] - original[component] * factor).abs() < 3e-14);
            }
        }
    }

    #[test]
    fn exact_plate_transmission_zero_retains_a_nonzero_pressure_pair() {
        // Exact phase limit k_Pn*h=pi, k_Sn*h=3pi/2, c_S/c_P=0.8.
        // The shared tangential slowness gives sin²(theta_S)=0.352. Both
        // parity channels have zero normal velocity, so T=0 and R=1; the
        // homogeneous pressure pair must survive without dividing by T.
        let ss = 0.352_f64.sqrt();
        let cs = (1.0_f64 - 0.352).sqrt();
        let pair = parity_boundary_pair(
            BoundedWave {
                c: -1.0,
                h: 0.0,
                l: 0.0,
                attenuation: 0.0,
            },
            BoundedWave {
                c: 0.0,
                h: -1.0 / cs,
                l: -cs,
                attenuation: 0.0,
            },
            ss,
            0.8,
            0.5,
            C64::ONE,
            C64::ONE,
        )
        .unwrap();
        assert_eq!(pair.transmission, 0.0);
        assert_eq!(pair.velocity, C64::ZERO);
        assert!(finite(pair.pressure) && pair.pressure != C64::ZERO);
        assert!(pair.transmission_log_scale.is_finite());
    }
}
