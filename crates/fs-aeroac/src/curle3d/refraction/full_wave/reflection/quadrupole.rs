//! Original symmetric volume stress reflected by the two-fluid plane.
//!
//! For the parent's W=q*dq*R*exp(i*zeta1*(a+b))/zeta1, the scalar reflected
//! Green function is i/(4*pi) integral W*J0. Its original source covector is
//! (q*m,zeta1), while its receiver covector is (q*m,-zeta1). Thus its source
//! Hessian has the transmitted owner's four angular channels; the receiver
//! normal changes sign in every third-derivative velocity channel. Source
//! trace and both physical shear contractions remain present. Equal original
//! speeds give the exact image stress S*M*S, S=I-2*n*n, in the incident fluid.
//!
//! Angular moments, cylindrical tensor reconstruction and numerical controls
//! belong to the original transmitted quadrupole owner. The original reflected
//! Green owner supplies branch maps, grazing cancellation and scaled weights.

use super::super::force::j1_over_argument;
use super::super::quadrupole::{
    FullWaveQuadrupoleDiagnostics, FullWaveQuadrupoleLimits, cartesian, j2_over_argument,
};
use super::super::{
    Dd, PlanarFluidInterface, QuadratureLimits, Spectrum, finite, integrate_spectrum, j0, j1,
    nonfinite, poll, positive,
};
use super::{
    ReflectionCoordinates, ReflectionKernel, ReflectionSpectralWeight, exact_reflection, field,
    line, singular_total,
};
use crate::AeroacError;
use crate::curle3d::CoherentCompactField;
use crate::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use crate::curle3d::statistics::ComplexSum;
use crate::observation::HarmonicFieldLine;
use core::f64::consts::{PI, SQRT_2};
use fs_math::{c64::C64, det};

/// Reflected field only; original source clock and complete stress tensor.
pub const FULL_WAVE_REFLECTED_QUADRUPOLE_SCOPE: &str = "Prescribed stationary symmetric \
3D Lighthill stress moment [N m]; reflected-only spherical-wave field in the incident \
fluid at one infinite plane between original homogeneous stationary lossless fluids, \
including trace, physical shear, reactive, transverse, critical and evanescent fields. \
Estimated quadrature error and analytic tail bounds are not interval certificates. \
No singular direct self-field, source-size validation, unresolved SGS sound, source \
feedback, absorption, mean flow, finite edges or measured-SPL qualification.";

/// Coherent pressure/velocity sum, including direct/reflected energy interference.
pub const FULL_WAVE_TOTAL_QUADRUPOLE_SCOPE: &str = "Prescribed stationary symmetric \
3D Lighthill stress moment [N m]; coherent direct-plus-reflected spherical-wave field \
in the original incident fluid at one infinite lossless fluid interface. Numerical \
relative control applies to the reflected transfer; only its absolute error estimates \
carry to total fields near interference nulls. No interval certificate, singular \
direct self-field, source-size validation, unresolved SGS sound, source feedback, \
absorption, mean flow, finite edges or measured-SPL qualification.";

type Transfer = [[C64; 6]; 4];
const ENTRIES: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

/// Cached original stress-to-pressure/velocity maps on the incident side.
/// Warm application is allocation-free and uses DD reductions of fixed 4-by-6
/// Cartesian transfers; no quadrature, derivatives or transcendental work.
/// A coincident source/receiver retains its finite reflected response, while
/// total requests refuse the singular prescribed point-source direct field.
#[derive(Clone, Debug)]
pub struct PreparedReflectedFullWaveQuadrupole {
    interface: PlanarFluidInterface,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    reflected: Transfer,
    total: Option<Transfer>,
    diagnostics: FullWaveQuadrupoleDiagnostics,
}

impl PreparedReflectedFullWaveQuadrupole {
    /// Source and receiver must lie strictly in the original incident fluid.
    /// Limits control the complete reflected stress transfer in the original
    /// fourteen-channel Frobenius norms, before a particular moment is applied.
    pub fn new(
        interface: PlanarFluidInterface,
        frequency_hz: f64,
        source_m: [f64; 3],
        observer_m: [f64; 3],
        limits: FullWaveQuadrupoleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        QuadratureLimits::from(limits).admit()?;
        let coordinates =
            ReflectionCoordinates::new(interface, frequency_hz, source_m, observer_m)?;
        let incident = interface.incident_state();
        if !positive(coordinates.omega / incident.sound_speed_m_s()) {
            return Err(nonfinite("reflected quadrupole wavenumber"));
        }
        let (reflected, diagnostics) = if incident.sound_speed_m_s()
            == interface.receiving_state().sound_speed_m_s()
        {
            (
                equal_speed_image(interface, frequency_hz, source_m, observer_m, cancelled)?,
                FullWaveQuadrupoleDiagnostics {
                    evaluations: 0,
                    maximum_depth_reached: 0,
                    retained_panels: 0,
                    transverse_wavenumber_cutoff_rad_m: 0.0,
                    maximum_bessel_argument: 0.0,
                    pressure_quadrature_error_pa_per_n_m: [0.0; 4],
                    velocity_quadrature_error_m_s_per_n_m: [0.0; 10],
                    pressure_tail_bound_pa_per_n_m: [0.0; 4],
                    velocity_tail_bound_m_s_per_n_m: [0.0; 10],
                    exact_equal_speed: true,
                },
            )
        } else {
            let log_pressure = -det::ln(4.0 * PI);
            let log_velocity =
                log_pressure - det::ln(coordinates.omega) - det::ln(incident.density_kg_m3());
            let mut common = ReflectionKernel::new(interface, coordinates)?;
            // Each angular moment is bounded by one, |R|<=1 above both
            // critical circles, and q,|zeta1|<=t+kmax. The effective counts
            // include squared sqrt(2) weights for physical symmetric shear.
            let [p, v] = common.prepare_tail(
                limits.into(),
                [2, 3],
                [log_pressure, log_velocity],
                [5, 14],
                cancelled,
            )?;
            let odd_p = if coordinates.range == 0.0 { 0.0 } else { p };
            let odd_v = if coordinates.range == 0.0 { 0.0 } else { v };
            let tail = [
                p,
                p,
                SQRT_2 * odd_p,
                p,
                odd_v,
                odd_v,
                SQRT_2 * v,
                odd_v,
                SQRT_2 * odd_v,
                SQRT_2 * v,
                v,
                v,
                SQRT_2 * odd_v,
                v,
            ];
            let kernel = ReflectionQuadrupoleKernel {
                common,
                log_pressure,
                log_velocity,
            };
            let result = integrate_spectrum::<14, 4>(&kernel, tail, limits.into(), cancelled)?;
            (
                cartesian(
                    result.value,
                    coordinates.tangent,
                    interface.receiving_normal(),
                )?,
                FullWaveQuadrupoleDiagnostics {
                    evaluations: result.evaluations,
                    maximum_depth_reached: result.maximum_depth_reached,
                    retained_panels: result.retained_panels,
                    transverse_wavenumber_cutoff_rad_m: result.cutoff,
                    maximum_bessel_argument: result.maximum_argument,
                    pressure_quadrature_error_pa_per_n_m: core::array::from_fn(|i| result.error[i]),
                    velocity_quadrature_error_m_s_per_n_m: core::array::from_fn(|i| {
                        result.error[i + 4]
                    }),
                    pressure_tail_bound_pa_per_n_m: core::array::from_fn(|i| tail[i]),
                    velocity_tail_bound_m_s_per_n_m: core::array::from_fn(|i| tail[i + 4]),
                    exact_equal_speed: false,
                },
            )
        };
        if reflected.iter().flatten().any(|v| !finite(*v)) {
            return Err(nonfinite("reflected quadrupole Cartesian transfer"));
        }
        let total = if source_m == observer_m {
            None
        } else {
            let direct = PreparedQuadrupoleObserver::new_cancellable(
                incident,
                frequency_hz,
                observer_m,
                &[source_m],
                1,
                &mut *cancelled,
            )?;
            let mut total = reflected;
            for column in 0..6 {
                poll(cancelled)?;
                let mut source = QuadrupolePhasor::default();
                source.stress_moment_n_m[column] = C64::ONE;
                let original = direct.evaluate_harmonic_line(&[source])?;
                let values = [
                    original.pressure_pa,
                    original.particle_velocity_m_s[0],
                    original.particle_velocity_m_s[1],
                    original.particle_velocity_m_s[2],
                ];
                for row in 0..4 {
                    let mut sum = ComplexSum::default();
                    sum.add(values[row])?;
                    sum.add(reflected[row][column])?;
                    total[row][column] = sum.value();
                }
            }
            if total.iter().flatten().any(|v| !finite(*v)) {
                return Err(nonfinite("direct-plus-reflected quadrupole transfer"));
            }
            Some(total)
        };
        poll(cancelled)?;
        Ok(Self {
            interface,
            frequency_hz,
            source_m,
            observer_m,
            reflected,
            total,
            diagnostics,
        })
    }

    pub const fn interface(&self) -> PlanarFluidInterface {
        self.interface
    }
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    pub const fn source_position_m(&self) -> [f64; 3] {
        self.source_m
    }
    pub const fn observer_position_m(&self) -> [f64; 3] {
        self.observer_m
    }
    /// Numerical evidence for reflected transfer only; original fourteen-channel
    /// ordering and physical Frobenius weights match the transmitted owner.
    pub const fn diagnostics(&self) -> &FullWaveQuadrupoleDiagnostics {
        &self.diagnostics
    }
    pub const fn has_total_field(&self) -> bool {
        self.total.is_some()
    }

    /// Physical symmetric stress [N m] on its original exp(-i omega t) clock.
    pub fn evaluate_reflected(
        &self,
        source: QuadrupolePhasor,
    ) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply(&self.reflected, source)?,
            FULL_WAVE_REFLECTED_QUADRUPOLE_SCOPE,
        )
    }

    /// Add full direct/reflected pressure and velocity before forming intensity.
    pub fn evaluate_total(
        &self,
        source: QuadrupolePhasor,
    ) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply(self.total.as_ref().ok_or_else(singular_total)?, source)?,
            FULL_WAVE_TOTAL_QUADRUPOLE_SCOPE,
        )
    }

    pub fn evaluate_reflected_harmonic_line(
        &self,
        source: QuadrupolePhasor,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(self.frequency_hz, apply(&self.reflected, source)?))
    }

    pub fn evaluate_total_harmonic_line(
        &self,
        source: QuadrupolePhasor,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            apply(self.total.as_ref().ok_or_else(singular_total)?, source)?,
        ))
    }
}

pub(in crate::curle3d::refraction::full_wave) fn apply(
    transfer: &Transfer,
    source: QuadrupolePhasor,
) -> Result<[C64; 4], AeroacError> {
    if source.stress_moment_n_m.iter().any(|v| !finite(*v)) {
        return Err(nonfinite("reflected quadrupole stress moment"));
    }
    let mut sums = [ComplexSum::default(); 4];
    for (column, moment) in source.stress_moment_n_m.into_iter().enumerate() {
        if moment != C64::ZERO {
            for (row, sum) in transfer.iter().zip(&mut sums) {
                sum.add(row[column] * moment)?;
            }
        }
    }
    let values = sums.map(ComplexSum::value);
    if values.iter().any(|v| !finite(*v)) {
        return Err(nonfinite("reflected quadrupole field amplitude"));
    }
    Ok(values)
}

fn equal_speed_image(
    interface: PlanarFluidInterface,
    frequency_hz: f64,
    source: [f64; 3],
    observer: [f64; 3],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Transfer, AeroacError> {
    let reflection = exact_reflection(interface);
    if reflection == 0.0 {
        return Ok([[C64::ZERO; 6]; 4]);
    }
    let n = interface.receiving_normal();
    let twice_distance = interface.distance(source) * Dd::from_f64(2.0);
    // Work with original DD separation rather than materializing a distant
    // Cartesian image location that can lose a small physical clearance.
    let displacement = core::array::from_fn(|i| {
        (Dd::from_f64(observer[i]) - Dd::from_f64(source[i]) + twice_distance * Dd::from_f64(n[i]))
            .to_f64()
    });
    let original = PreparedQuadrupoleObserver::new_cancellable(
        interface.incident_state(),
        frequency_hz,
        displacement,
        &[[0.0; 3]],
        1,
        &mut *cancelled,
    )?;
    let mirror = |i: usize, j: usize| (if i == j { 1.0 } else { 0.0 }) - 2.0 * n[i] * n[j];
    let mut transfer = [[C64::ZERO; 6]; 4];
    for (column, (a, b)) in ENTRIES.into_iter().enumerate() {
        poll(cancelled)?;
        let image = QuadrupolePhasor {
            stress_moment_n_m: ENTRIES.map(|(i, j)| {
                C64::from_re(
                    mirror(i, a) * mirror(j, b)
                        + if a != b {
                            mirror(i, b) * mirror(j, a)
                        } else {
                            0.0
                        },
                )
            }),
        };
        let value = original.evaluate_harmonic_line(&[image])?;
        transfer[0][column] = value.pressure_pa.scale(reflection);
        for axis in 0..3 {
            transfer[axis + 1][column] = value.particle_velocity_m_s[axis].scale(reflection);
        }
    }
    Ok(transfer)
}

struct ReflectionQuadrupoleKernel {
    common: ReflectionKernel,
    log_pressure: f64,
    log_velocity: f64,
}

impl Spectrum<14> for ReflectionQuadrupoleKernel {
    fn regions(&self) -> usize {
        <ReflectionKernel as Spectrum<3>>::regions(&self.common)
    }
    fn extent(&self) -> (f64, f64) {
        <ReflectionKernel as Spectrum<3>>::extent(&self.common)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        <ReflectionKernel as Spectrum<3>>::resolved(&self.common, region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 14], AeroacError> {
        quadrupole_channels(
            self.common.sample_weight(region, x)?,
            self.common.coordinates.range,
            self.log_pressure,
            self.log_velocity,
        )
    }
}

/// Shared original source Hessian and returning receiver derivative. A planar
/// layer changes only the reflected scalar weight, never these angular signs,
/// physical shear factors or logarithmically balanced derivative products.
pub(in crate::curle3d::refraction::full_wave) fn quadrupole_channels(
    weight: ReflectionSpectralWeight,
    range: f64,
    log_pressure: f64,
    log_velocity: f64,
) -> Result<[C64; 14], AeroacError> {
    let argument = weight.q * range;
    if !argument.is_finite() {
        return Err(nonfinite("reflected quadrupole Bessel argument"));
    }
    let zero = j0(argument);
    let one = j1(argument);
    let ratio = j1_over_argument(argument, one);
    let derivative = zero - ratio;
    let second_ratio = j2_over_argument(argument, zero, ratio);
    let third = one - second_ratio;
    let q_log = det::ln(weight.q);
    let magnitude = weight.zeta.abs();
    let z_log = det::ln(magnitude);
    let z = if magnitude == 0.0 {
        C64::ZERO
    } else {
        C64::new(weight.zeta.re / magnitude, weight.zeta.im / magnitude)
    };
    let zz = z * z;
    let i = C64::new(0.0, 1.0);
    let p = log_pressure;
    let v = log_velocity;
    // Scale physical derivatives BEFORE exponentiation. Even q³ or zeta³
    // can restore a representable field after the bare decay underflows.
    // Preserve the original four pressure and TEN velocity diagnostics.
    let values = [
        weight.scaled(-i, p + 2.0 * q_log, derivative),
        weight.scaled(-i, p + 2.0 * q_log, ratio),
        weight.scaled(z, p + q_log + z_log, SQRT_2 * one),
        weight.scaled(-i * zz, p + 2.0 * z_log, zero),
        weight.scaled(C64::ONE, v + 3.0 * q_log, third),
        weight.scaled(C64::ONE, v + 3.0 * q_log, second_ratio),
        weight.scaled(-i * z, v + 2.0 * q_log + z_log, SQRT_2 * derivative),
        weight.scaled(zz, v + q_log + 2.0 * z_log, one),
        weight.scaled(C64::ONE, v + 3.0 * q_log, SQRT_2 * second_ratio),
        weight.scaled(-i * z, v + 2.0 * q_log + z_log, SQRT_2 * ratio),
        weight.scaled(i * z, v + 2.0 * q_log + z_log, derivative),
        weight.scaled(i * z, v + 2.0 * q_log + z_log, ratio),
        weight.scaled(-zz, v + q_log + 2.0 * z_log, SQRT_2 * one),
        weight.scaled(i * zz * z, v + 3.0 * z_log, zero),
    ];
    if values.iter().any(|value| !finite(*value)) {
        return Err(nonfinite("reflected quadrupole spectrum"));
    }
    Ok(values)
}
