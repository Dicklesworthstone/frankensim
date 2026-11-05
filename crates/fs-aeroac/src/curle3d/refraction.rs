//! The transmitted geometrical-acoustics arrival through a stationary fluid plane.
//!
//! Original fluids supply density and sound speed independently. Fermat's path
//! preserves `sin(theta_i)/c_i`; pressure and normal-velocity continuity give
//! `Tp = 2 Z2 cos(theta1)/(Z2 cos(theta1) + Z1 cos(theta2))`. See C. C. Mei,
//! MIT Wave Propagation, chapter 3, equations (1.21), (1.26):
//! <https://web.mit.edu/1.138j/www/material/chap-3.pdf>.
//!
//! A spherical source also needs refracted ray-tube spreading, not inverse
//! total path length. For leg lengths L1,L2, cosines u1,u2 and m=c2/c1,
//! `D² = (L1+m L2)(L1+m L2 (u1/u2)²)`. The perpendicular ray-tube area per
//! source solid angle is `(u2/u1) D²`; transmitted power and receiving
//! impedance give `p/Q = -i omega rho1 Tp exp(i omega T)/(4 pi D)`.
//! This is also the stationary-phase term of the transmitted Weyl integral.
//! Ray-tube reference: Hudimac, JASA 29 (1957), 916–917,
//! <https://doi.org/10.1121/1.1909097>.
//!
//! This is a high-frequency approximation, not the exact spherical-wave field.
//! Only the selected transmitted ray is returned. Leading-order velocity is
//! `p n2/Z2`; near-field induction, head/lateral waves, evanescent transmission,
//! critical-angle diffraction, mean flow, absorption and source feedback are
//! absent. Caller-selected clearance and grazing limits are admission controls,
//! not certified wave-error bounds. Full-wave material/FE owners remain
//! necessary outside that envelope. No material or source is inferred.

use super::CoherentCompactField;
use super::statistics::{finite, invalid, nonfinite};
use crate::AeroacError;
use crate::observation::HarmonicFieldLine;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd, det};

/// Applicability retained on every returned field.
pub const REFRACTED_MONOPOLE_SCOPE: &str = "Prescribed stationary compact 3D volume source; \
leading-order high-frequency transmitted ray across one infinite planar boundary between \
two stationary homogeneous lossless fluids. No near field, head/lateral waves, evanescent \
or critical diffraction, mean flow, absorption, source feedback or measured-SPL qualification.";

/// Infinite plane separating two original, stationary, lossless fluid snapshots.
/// The normal points from the incident (negative) into the receiving half-space.
#[derive(Clone, Copy, Debug)]
pub struct PlanarFluidInterface {
    incident: FluidAcousticState,
    receiving: FluidAcousticState,
    point_m: [f64; 3],
    normal: [f64; 3],
}

impl PlanarFluidInterface {
    /// Retain both states and normalize a finite, nonzero receiving normal.
    pub fn new(
        incident: FluidAcousticState,
        receiving: FluidAcousticState,
        point_m: [f64; 3],
        receiving_normal: [f64; 3],
    ) -> Result<Self, AeroacError> {
        if point_m
            .iter()
            .chain(&receiving_normal)
            .any(|v| !v.is_finite())
        {
            return Err(nonfinite("acoustic transmission plane geometry"));
        }
        let scale = receiving_normal.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
        if scale == 0.0 {
            return Err(invalid(
                "acoustic transmission plane needs a nonzero normal",
            ));
        }
        let scaled = receiving_normal.map(|v| v / scale);
        let length = norm(scaled);
        Ok(Self {
            incident,
            receiving,
            point_m,
            normal: scaled.map(|v| v / length),
        })
    }

    /// Original incident state.
    pub const fn incident_state(self) -> FluidAcousticState {
        self.incident
    }
    /// Original receiving state.
    pub const fn receiving_state(self) -> FluidAcousticState {
        self.receiving
    }
    /// Original plane origin [m].
    pub const fn point_m(self) -> [f64; 3] {
        self.point_m
    }
    /// Unit normal from the incident into the receiving fluid.
    pub const fn receiving_normal(self) -> [f64; 3] {
        self.normal
    }

    fn distance(self, position: [f64; 3]) -> Dd {
        let mut result = Dd::ZERO;
        for axis in 0..3 {
            result = result
                + (Dd::from_f64(position[axis]) - Dd::from_f64(self.point_m[axis]))
                    * Dd::from_f64(self.normal[axis]);
        }
        result
    }
}

/// Explicit physical applicability and numerical-work controls. No defaults.
#[derive(Clone, Copy, Debug)]
pub struct RefractedMonopoleLimits {
    /// Minimum k_i*h_i on EACH side, where h_i is perpendicular clearance [m].
    /// Positive and finite; choosing it does not certify an approximation error.
    pub minimum_kh: f64,
    /// Minimum cosine on BOTH ray legs, in (0,1]; bounds grazing/critical approach.
    pub minimum_direction_cosine: f64,
    /// Maximum bisection half-width along the plane [m], excluding roundoff.
    pub path_tolerance_m: f64,
    /// Maximum bisection iterations, in 1..=128. Exhaustion refuses the result.
    pub maximum_iterations: usize,
}

/// Original geometry and the selected transmitted geometrical arrival.
#[derive(Clone, Copy, Debug)]
pub struct RefractedMonopolePath {
    /// Crossing position on the original plane [m].
    pub interface_position_m: [f64; 3],
    /// Source-to-plane ray direction.
    pub incident_direction: [f64; 3],
    /// Plane-to-microphone ray direction.
    pub receiving_direction: [f64; 3],
    /// Incident and receiving propagation distances [m].
    pub leg_lengths_m: [f64; 2],
    /// Incident and receiving cosines to the interface normal.
    pub direction_cosines: [f64; 2],
    /// Physical travel time L1/c1 + L2/c2 [s], already in the transfer phase.
    pub delay_s: f64,
    /// Pressure amplitude T at the plane; may exceed one without creating power.
    pub pressure_transmission: f64,
    /// Ratio of transmitted to incident normal acoustic power at the plane.
    pub normal_power_transmission: f64,
    /// Refracted spreading D [m], generally different from L1+L2.
    pub spreading_distance_m: f64,
    /// Original incident/receiving k_i*h_i, checked against the caller's envelope.
    pub clearance_kh: [f64; 2],
    /// Final root-bracket half-width [m]; not an interval-arithmetic certificate.
    pub crossing_half_width_m: f64,
}

/// Cached source-to-microphone transfer on one original physical source clock.
/// Preparation and application allocate nothing. Preparation uses a bounded
/// monotone root solve; application uses four complex multiplies and no trig.
#[derive(Clone, Debug)]
pub struct PreparedRefractedMonopole {
    interface: PlanarFluidInterface,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    path: RefractedMonopolePath,
    transfer: [C64; 4],
}

impl PreparedRefractedMonopole {
    /// Prepare one physical volume source [m³/s] and one microphone, strictly
    /// on opposite sides. No reflected path or other source type is inserted.
    pub fn new(
        interface: PlanarFluidInterface,
        frequency_hz: f64,
        source_m: [f64; 3],
        observer_m: [f64; 3],
        limits: RefractedMonopoleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        if source_m.iter().chain(&observer_m).any(|v| !v.is_finite()) {
            return Err(nonfinite("refracted monopole positions"));
        }
        let omega = core::f64::consts::TAU * frequency_hz;
        if !positive(omega)
            || !positive(limits.minimum_kh)
            || !positive(limits.path_tolerance_m)
            || !positive(limits.minimum_direction_cosine)
            || limits.minimum_direction_cosine > 1.0
            || !(1..=128).contains(&limits.maximum_iterations)
        {
            return Err(invalid(
                "refracted monopole frequency, envelope or work limits",
            ));
        }
        let a_dd = -interface.distance(source_m);
        let b_dd = interface.distance(observer_m);
        let a = a_dd.to_f64();
        let b = b_dd.to_f64();
        if !positive(a) || !positive(b) {
            return Err(invalid(
                "source and receiver must lie strictly on their fluid sides",
            ));
        }
        let c1 = interface.incident.sound_speed_m_s();
        let c2 = interface.receiving.sound_speed_m_s();
        let clearance_kh = [(omega / c1) * a, (omega / c2) * b];
        if clearance_kh
            .iter()
            .any(|v| !v.is_finite() || *v < limits.minimum_kh)
        {
            return Err(invalid(
                "refracted monopole is outside the admitted high-frequency clearance",
            ));
        }
        let lateral = core::array::from_fn(|axis| {
            (Dd::from_f64(observer_m[axis])
                - Dd::from_f64(source_m[axis])
                - (a_dd + b_dd) * Dd::from_f64(interface.normal[axis]))
            .to_f64()
        });
        let range = norm(lateral);
        if !range.is_finite() {
            return Err(nonfinite("refracted monopole lateral separation"));
        }
        let tangent = if range == 0.0 {
            [0.0; 3]
        } else {
            lateral.map(|v| v / range)
        };
        let (crossing, half_width) = crossing(a, b, range, c1, c2, limits, cancelled)?;
        let l1 = det::hypot(a, crossing);
        let l2 = det::hypot(b, range - crossing);
        let u1 = a / l1;
        let u2 = b / l2;
        if [u1, u2]
            .iter()
            .any(|u| !u.is_finite() || *u < limits.minimum_direction_cosine)
        {
            return Err(invalid(
                "refracted monopole is outside the admitted grazing envelope",
            ));
        }
        let incident_direction =
            core::array::from_fn(|i| tangent[i] * (crossing / l1) + interface.normal[i] * u1);
        let receiving_direction = core::array::from_fn(|i| {
            tangent[i] * ((range - crossing) / l2) + interface.normal[i] * u2
        });
        let interface_position_m = core::array::from_fn(|i| {
            (Dd::from_f64(source_m[i])
                + a_dd * Dd::from_f64(interface.normal[i])
                + Dd::from_f64(crossing) * Dd::from_f64(tangent[i]))
            .to_f64()
        });
        let z1 = interface.incident.characteristic_impedance_pa_s_m();
        let z2 = interface.receiving.characteristic_impedance_pa_s_m();
        let impedance_scale = z1.max(z2);
        let za = (z2 / impedance_scale) * u1;
        let zb = (z1 / impedance_scale) * u2;
        let tp = 2.0 * (za / (za + zb));
        let power = tp * (2.0 * (zb / (za + zb)));
        // The two principal ray-tube widths remain different after refraction.
        // Multiplying square roots avoids overflowing their squared area.
        let m_l2 = (c2 / c1) * l2;
        let cosine_ratio = u1 / u2;
        let spreading = det::sqrt(l1 + m_l2) * det::sqrt(l1 + (m_l2 * cosine_ratio) * cosine_ratio);
        let delay = Dd::from_f64(l1) / Dd::from_f64(c1) + Dd::from_f64(l2) / Dd::from_f64(c2);
        let phase = delay * Dd::from_f64(omega);
        let scale = (omega / (4.0 * core::f64::consts::PI))
            * (interface.incident.density_kg_m3() / spreading)
            * tp;
        if !positive(spreading)
            || !positive(delay.to_f64())
            || !positive(scale)
            || !phase.hi.is_finite()
            || !phase.lo.is_finite()
            || !tp.is_finite()
            || !power.is_finite()
            || interface_position_m.iter().any(|v| !v.is_finite())
        {
            return Err(nonfinite("refracted monopole path or transfer scales"));
        }
        let (sin, cos) = det::sin_cos_dd(phase);
        let pressure = C64::new(sin.to_f64(), -cos.to_f64()).scale(scale);
        let transfer = [
            pressure,
            pressure.scale(receiving_direction[0] / z2),
            pressure.scale(receiving_direction[1] / z2),
            pressure.scale(receiving_direction[2] / z2),
        ];
        if transfer.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("refracted monopole field transfer"));
        }
        poll(cancelled)?;
        Ok(Self {
            interface,
            frequency_hz,
            source_m,
            observer_m,
            transfer,
            path: RefractedMonopolePath {
                interface_position_m,
                incident_direction,
                receiving_direction,
                leg_lengths_m: [l1, l2],
                direction_cosines: [u1, u2],
                delay_s: delay.to_f64(),
                pressure_transmission: tp,
                normal_power_transmission: power,
                spreading_distance_m: spreading,
                clearance_kh,
                crossing_half_width_m: half_width,
            },
        })
    }

    /// Selected transmitted ray and its physical scales.
    pub const fn path(&self) -> &RefractedMonopolePath {
        &self.path
    }
    /// Original plane and fluid snapshots.
    pub const fn interface(&self) -> PlanarFluidInterface {
        self.interface
    }
    /// Original positive frequency [Hz].
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    /// Original source position [m].
    pub const fn source_position_m(&self) -> [f64; 3] {
        self.source_m
    }
    /// Original microphone position [m].
    pub const fn observer_position_m(&self) -> [f64; 3] {
        self.observer_m
    }

    fn amplitudes(&self, volume_flow_m3_s: C64) -> Result<[C64; 4], AeroacError> {
        if !finite(volume_flow_m3_s) {
            return Err(nonfinite("refracted monopole volume-flow phasor"));
        }
        let amplitudes = self.transfer.map(|v| v * volume_flow_m3_s);
        if amplitudes.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("refracted monopole field amplitude"));
        }
        Ok(amplitudes)
    }

    /// Coherent peak phasor Q [m³/s], exp(-i omega t), on the source clock.
    /// Propagation phase is already applied. Velocity/intensity are leading
    /// geometrical-acoustics terms; reactive near-field energy is not modeled.
    pub fn evaluate(&self, volume_flow_m3_s: C64) -> Result<CoherentCompactField, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(volume_flow_m3_s)?;
        let velocity = [ux, uy, uz];
        let intensity = velocity.map(|v| pressure.scale(0.5) * v.conj());
        if intensity.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("refracted monopole intensity"));
        }
        Ok(CoherentCompactField {
            pressure_pa: pressure,
            particle_velocity_m_s: velocity,
            active_intensity_w_m2: intensity.map(|v| v.re),
            reactive_intensity_w_m2: intensity.map(|v| v.im),
            scope: REFRACTED_MONOPOLE_SCOPE,
        })
    }

    /// Feed the existing finite-window/receiver observer without a second delay.
    /// The original high-frequency scope still applies to every supplied line.
    pub fn evaluate_harmonic_line(
        &self,
        volume_flow_m3_s: C64,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(volume_flow_m3_s)?;
        Ok(HarmonicFieldLine {
            frequency_hz: self.frequency_hz,
            pressure_pa: pressure,
            particle_velocity_m_s: [ux, uy, uz],
        })
    }
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn norm(x: [f64; 3]) -> f64 {
    det::hypot(det::hypot(x[0], x[1]), x[2])
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AeroacError> {
    if cancelled() {
        Err(AeroacError::Cancelled {
            what: "planar acoustic refraction",
        })
    } else {
        Ok(())
    }
}

// Monotone Fermat derivative: sin(theta1)/c1 - sin(theta2)/c2.
// Normalize speeds rather than multiplying unbounded slownesses/distances.
fn crossing(
    a: f64,
    b: f64,
    range: f64,
    c1: f64,
    c2: f64,
    limits: RefractedMonopoleLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64), AeroacError> {
    if range == 0.0 {
        return Ok((0.0, 0.0));
    }
    let speed = c1.max(c2);
    let mut low = 0.0;
    let mut high = range;
    for _ in 0..limits.maximum_iterations {
        poll(cancelled)?;
        let half_width = (high - low) * 0.5;
        let x = low + half_width;
        if half_width <= limits.path_tolerance_m {
            return Ok((x, half_width));
        }
        if x == low || x == high {
            break;
        }
        let slope = (c2 / speed) * (x / det::hypot(a, x))
            - (c1 / speed) * ((range - x) / det::hypot(b, range - x));
        if slope > 0.0 {
            high = x;
        } else {
            low = x;
        }
    }
    let half_width = (high - low) * 0.5;
    if half_width <= limits.path_tolerance_m {
        return Ok((low + half_width, half_width));
    }
    Err(AeroacError::NotConverged {
        what: "planar acoustic crossing half-width [m]",
        residual: half_width,
    })
}
