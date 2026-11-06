//! Analytic transmission of the original symmetric volume-stress moment.
//!
//! For the parent's outgoing scalar Green function H, p=M_ij d_si d_sj H
//! and u_l=d_xl p/(i omega rho2). Source-normal derivatives use zeta1;
//! receiver-normal derivatives use zeta2. All six physical stress entries,
//! their trace, and transverse receiving motion are retained.
//!
//! In the receiver radial/transverse/normal frame (t,s,n), put f=J1(x)/x,
//! g=J2(x)/x=-f', d=J0-f=J1', and e=J1-g=-J1'', x=qR. Angular integration
//! gives pressure (tt,ss,tn,nn) proportional to
//! [-i q²d,-i q²f,q zeta1 J1,-i zeta1²J0]. Receiving differentiation gives
//! velocity (ttt,tss,ttn,tnn,sts,ssn,ntt,nss,ntn,nnn) proportional to
//! [q³e,q³g,-i q²zeta1 d,q zeta1²J1,q³g,-i q²zeta1 f,
//!  -i q²zeta2 d,-i q²zeta2 f,q zeta1 zeta2 J1,-i zeta1²zeta2 J0].
//! The pressure/velocity prefactors are rho2/(2 pi), 1/(2 pi omega).
//! Their shared weight is q*dq*E/(rho2*zeta1+rho1*zeta2).
//! DLMF 10.9.2 gives these angular moments and their absolute bounds;
//! 10.6.2 gives the recurrences: <https://dlmf.nist.gov/10.9>,
//! <https://dlmf.nist.gov/10.6>. The small-x continuation avoids subtraction
//! of two nearly equal order-one Bessel values when forming J2/x.

use super::force::j1_over_argument;
use super::reflection::{field, line, prepare_reflected_tail};
use super::{
    ComplexSum, Dd, Kernel, PlanarFluidInterface, QuadratureLimits, SpectralPoint, Spectrum,
    finite, integrate_spectrum, invalid, j0, j1, log_weighted, minus_i, nonfinite, norm, poll,
    positive,
};
use crate::AeroacError;
use crate::curle3d::CoherentCompactField;
use crate::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use crate::observation::HarmonicFieldLine;
use core::f64::consts::{FRAC_1_SQRT_2, SQRT_2, TAU};
use fs_math::{c64::C64, det};

/// Coherent original volume cells and geometry-bound recorded stress spectra.
pub mod array;

/// Physical scope of the transmitted stress field, including its original clock.
pub const FULL_WAVE_QUADRUPOLE_SCOPE: &str = "Prescribed stationary symmetric 3D \
Lighthill stress moment [N m] in the original incident fluid; outgoing spherical \
full-wave transmission across one infinite planar interface between stationary \
homogeneous lossless fluids, including trace, shear, near, transverse, critical \
and evanescent fields. Estimated quadrature error and analytic tail bounds are \
not interval certificates. No source-size validation, unresolved SGS sound, \
mean flow, absorption, finite boundaries, source feedback or measured-SPL qualification.";

/// Accuracy of the unit symmetric stress-to-field map, before any particular
/// source amplitude is applied. No numerical defaults or ray fallback.
#[derive(Clone, Copy, Debug)]
pub struct FullWaveQuadrupoleLimits {
    /// Absolute Frobenius pressure-Hessian tolerance [Pa/(N m)].
    pub absolute_pressure_tolerance_pa_per_n_m: f64,
    /// Absolute Frobenius velocity-tensor tolerance [(m/s)/(N m)].
    pub absolute_velocity_tolerance_m_s_per_n_m: f64,
    /// Relative transfer-norm tolerance in [0,1), combined with absolute by max.
    /// It does not promise relative field accuracy at a stress-pattern null.
    pub relative_tolerance: f64,
    /// Shared fourteen-channel work cap in 1..=MAX_FULL_WAVE_EVALUATIONS.
    /// Exactly equal original speeds require no spectral evaluations.
    pub maximum_evaluations: usize,
    /// Transformed-region bisection depth in 1..=48; root depth is zero.
    pub maximum_depth: usize,
}

impl From<FullWaveQuadrupoleLimits> for QuadratureLimits {
    fn from(limits: FullWaveQuadrupoleLimits) -> Self {
        Self {
            absolute_pressure: limits.absolute_pressure_tolerance_pa_per_n_m,
            absolute_velocity: limits.absolute_velocity_tolerance_m_s_per_n_m,
            relative: limits.relative_tolerance,
            maximum_evaluations: limits.maximum_evaluations,
            maximum_depth: limits.maximum_depth,
        }
    }
}

/// Pressure channels are (tt,ss,sqrt(2)tn,nn). Velocity channels are
/// (ttt,tss,sqrt(2)ttn,tnn,sqrt(2)sts,sqrt(2)ssn,ntt,nss,sqrt(2)ntn,nnn),
/// where the first index is the receiving direction and the last two are
/// symmetric source-stress indices. Their Euclidean norms are the complete
/// pressure/velocity transfer Frobenius norms, including physical shear.
/// Multiply those errors by ||M||_F, with off-diagonal source entries counted
/// twice, to bound a particular stress field. Geometry, native Bessel and other
/// floating-point errors are outside the estimates and analytic tail bounds.
#[derive(Clone, Copy, Debug)]
pub struct FullWaveQuadrupoleDiagnostics {
    pub evaluations: usize,
    pub maximum_depth_reached: usize,
    pub retained_panels: usize,
    pub transverse_wavenumber_cutoff_rad_m: f64,
    pub maximum_bessel_argument: f64,
    pub pressure_quadrature_error_pa_per_n_m: [f64; 4],
    pub velocity_quadrature_error_m_s_per_n_m: [f64; 10],
    pub pressure_tail_bound_pa_per_n_m: [f64; 4],
    pub velocity_tail_bound_m_s_per_n_m: [f64; 10],
    /// True only when the original sound speeds match exactly.
    pub exact_equal_speed: bool,
}

/// Cached Cartesian pressure/velocity response to one original stress moment.
/// Warm application uses fixed-order DD reductions and no allocations,
/// quadrature, finite differences, transcendental calls or extra travel delay.
#[derive(Clone, Debug)]
pub struct PreparedFullWaveQuadrupole {
    interface: PlanarFluidInterface,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    // Rows p,ux,uy,uz; physical columns Mxx,Myy,Mzz,Mxy,Mxz,Myz [N m].
    pub(super) transfer: [[C64; 6]; 4],
    diagnostics: FullWaveQuadrupoleDiagnostics,
}

impl PreparedFullWaveQuadrupole {
    /// Prepare source and receiver strictly on opposite original fluid sides.
    /// Every source cell is a prescribed compact stress moment; the owner does
    /// not infer unresolved stresses, source quadrature accuracy or calibration.
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
        let omega = TAU * frequency_hz;
        if !positive(omega) {
            return Err(invalid("full-wave quadrupole frequency"));
        }
        if source_m.iter().chain(&observer_m).any(|v| !v.is_finite()) {
            return Err(nonfinite("full-wave quadrupole positions"));
        }
        let a = -interface.distance(source_m);
        let b = interface.distance(observer_m);
        if !positive(a.to_f64()) || !positive(b.to_f64()) || !positive((a + b).to_f64()) {
            return Err(invalid(
                "full-wave stress source and observer must lie strictly on their fluid sides",
            ));
        }
        let lateral = core::array::from_fn(|i| {
            (Dd::from_f64(observer_m[i])
                - Dd::from_f64(source_m[i])
                - (a + b) * Dd::from_f64(interface.normal[i]))
            .to_f64()
        });
        let range = norm(lateral);
        if !range.is_finite() {
            return Err(nonfinite("full-wave quadrupole lateral separation"));
        }
        let tangent = if range == 0.0 {
            [0.0; 3]
        } else {
            lateral.map(|v| v / range)
        };
        let c1 = interface.incident.sound_speed_m_s();
        let c2 = interface.receiving.sound_speed_m_s();
        let rho1 = interface.incident.density_kg_m3();
        let rho2 = interface.receiving.density_kg_m3();
        let density_scale = rho1.max(rho2);
        let density = [rho1 / density_scale, rho2 / density_scale];
        if density.iter().any(|v| !positive(*v)) {
            return Err(nonfinite("full-wave quadrupole density contrast"));
        }
        let (transfer, diagnostics) = if c1 == c2 {
            // The native homogeneous tensor carries the original source fluid.
            // Its pressure scales by Tp, while its velocity scales by Tp*rho1/rho2.
            // This retains the exact full near field for a pure density jump.
            let original = PreparedQuadrupoleObserver::new_cancellable(
                interface.incident,
                frequency_hz,
                observer_m,
                &[source_m],
                1,
                &mut *cancelled,
            )?;
            let scales = [
                2.0 * density[1] / (density[0] + density[1]),
                2.0 * density[0] / (density[0] + density[1]),
            ];
            let mut transfer = [[C64::ZERO; 6]; 4];
            for column in 0..6 {
                poll(cancelled)?;
                let mut source = QuadrupolePhasor::default();
                source.stress_moment_n_m[column] = C64::ONE;
                let original = original.evaluate_harmonic_line(&[source])?;
                transfer[0][column] = original.pressure_pa.scale(scales[0]);
                for axis in 0..3 {
                    transfer[axis + 1][column] =
                        original.particle_velocity_m_s[axis].scale(scales[1]);
                }
            }
            (
                transfer,
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
            let k = [omega / c1, omega / c2];
            if k.iter().any(|v| !positive(*v)) || k[0] == k[1] {
                return Err(nonfinite("full-wave quadrupole distinct fluid wavenumbers"));
            }
            let kmin = k[0].min(k[1]);
            let kmax = k[0].max(k[1]);
            let gap = det::sqrt(kmax - kmin) * det::sqrt(kmax) * det::sqrt(1.0 + kmin / kmax);
            let pressure_scale = density[1] / TAU;
            let velocity_scale = ((1.0 / TAU) / omega) / density_scale;
            if !positive(gap) || !positive(pressure_scale) || !positive(velocity_scale) {
                return Err(nonfinite("full-wave quadrupole spectral scales"));
            }
            let mut kernel = QuadrupoleKernel(Kernel {
                a,
                b,
                range,
                density,
                kmin,
                kmax,
                gap,
                incident_has_larger_k: k[0] > k[1],
                pressure_scale,
                velocity_scale,
                tail_extent: 0.0,
            });
            let tail = kernel.prepare_tail(limits, cancelled)?;
            let result = integrate_spectrum::<14, 4>(&kernel, tail, limits.into(), cancelled)?;
            (
                cartesian(result.value, tangent, interface.normal)?,
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
        if transfer.iter().flatten().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave quadrupole Cartesian transfer"));
        }
        poll(cancelled)?;
        Ok(Self {
            interface,
            frequency_hz,
            source_m,
            observer_m,
            transfer,
            diagnostics,
        })
    }

    pub const fn diagnostics(&self) -> &FullWaveQuadrupoleDiagnostics {
        &self.diagnostics
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

    fn amplitudes(&self, source: QuadrupolePhasor) -> Result<[C64; 4], AeroacError> {
        if source.stress_moment_n_m.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave quadrupole stress moment"));
        }
        let mut values = [C64::ZERO; 4];
        for (row, coefficients) in self.transfer.iter().enumerate() {
            let mut sum = ComplexSum::default();
            for (coefficient, moment) in coefficients.iter().zip(source.stress_moment_n_m) {
                sum.add(*coefficient * moment)?;
            }
            values[row] = sum.value();
        }
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave quadrupole field amplitude"));
        }
        Ok(values)
    }

    /// Physical symmetric stress [N m], retaining the original exp(-i omega t)
    /// clock. Intensity uses the full coherent pressure and vector velocity.
    pub fn evaluate(&self, source: QuadrupolePhasor) -> Result<CoherentCompactField, AeroacError> {
        field(self.amplitudes(source)?, FULL_WAVE_QUADRUPOLE_SCOPE)
    }

    /// Native receiver line with the full propagation phase already applied.
    pub fn evaluate_harmonic_line(
        &self,
        source: QuadrupolePhasor,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(self.frequency_hz, self.amplitudes(source)?))
    }
}

/// Rotate the normalized cylindrical tensor without choosing an arbitrary
/// transverse axis. P=I-nn and S=P-tt retain the exact axis limit t=0.
pub(super) fn cartesian(
    a: [C64; 14],
    t: [f64; 3],
    n: [f64; 3],
) -> Result<[[C64; 6]; 4], AeroacError> {
    let [
        ptt,
        pss,
        ptn,
        pnn,
        vttt,
        vtss,
        vttn,
        vtnn,
        vsts,
        vssn,
        vntt,
        vnss,
        vntn,
        vnnn,
    ] = a;
    let ptn = ptn.scale(FRAC_1_SQRT_2);
    let vttn = vttn.scale(FRAC_1_SQRT_2);
    let vsts = vsts.scale(FRAC_1_SQRT_2);
    let vssn = vssn.scale(FRAC_1_SQRT_2);
    let vntn = vntn.scale(FRAC_1_SQRT_2);
    let p = |i: usize, j: usize| (if i == j { 1.0 } else { 0.0 }) - n[i] * n[j];
    let s = |i: usize, j: usize| p(i, j) - t[i] * t[j];
    let mut transfer = [[C64::ZERO; 6]; 4];
    for (entry, (i, j)) in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
        .into_iter()
        .enumerate()
    {
        let multiplicity = if i == j { 1.0 } else { 2.0 };
        let tn = t[i] * n[j] + n[i] * t[j];
        let mut pressure = ComplexSum::default();
        pressure.add(pss.scale(p(i, j)))?;
        pressure.add((ptt - pss).scale(t[i] * t[j]))?;
        pressure.add(ptn.scale(tn))?;
        pressure.add(pnn.scale(n[i] * n[j]))?;
        transfer[0][entry] = pressure.value().scale(multiplicity);
        for l in 0..3 {
            let mut velocity = ComplexSum::default();
            velocity.add(vtss.scale(t[l] * p(i, j)))?;
            velocity.add((vttt - vtss).scale(t[l] * t[i] * t[j]))?;
            velocity.add(vttn.scale(t[l] * tn))?;
            velocity.add(vtnn.scale(t[l] * n[i] * n[j]))?;
            velocity.add(vsts.scale(t[i] * s(l, j) + t[j] * s(l, i)))?;
            velocity.add(vssn.scale(n[i] * s(l, j) + n[j] * s(l, i)))?;
            velocity.add(vnss.scale(n[l] * p(i, j)))?;
            velocity.add((vntt - vnss).scale(n[l] * t[i] * t[j]))?;
            velocity.add(vntn.scale(n[l] * tn))?;
            velocity.add(vnnn.scale(n[l] * n[i] * n[j]))?;
            transfer[l + 1][entry] = velocity.value().scale(multiplicity);
        }
    }
    Ok(transfer)
}

struct QuadrupoleKernel(Kernel);

impl Spectrum<14> for QuadrupoleKernel {
    fn extent(&self) -> (f64, f64) {
        Spectrum::<3>::extent(&self.0)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        self.0.resolved(region, lower, upper)
    }

    fn sample(&self, region: usize, x: f64) -> Result<[C64; 14], AeroacError> {
        let kernel = &self.0;
        let SpectralPoint { q, zeta, measure } = kernel.point(region, x);
        let phase = kernel.a * Dd::from_f64(zeta[0].re) + kernel.b * Dd::from_f64(zeta[1].re);
        let attenuation =
            (kernel.a * Dd::from_f64(zeta[0].im) + kernel.b * Dd::from_f64(zeta[1].im)).to_f64();
        let argument = q * kernel.range;
        if !phase.hi.is_finite()
            || !phase.lo.is_finite()
            || !attenuation.is_finite()
            || !argument.is_finite()
        {
            return Err(nonfinite("full-wave quadrupole spectral phase"));
        }
        let (sin, cos) = det::sin_cos_dd(phase);
        let phasor = C64::new(cos.to_f64(), sin.to_f64());
        let denominator = zeta[0].scale(kernel.density[1]) + zeta[1].scale(kernel.density[0]);
        let zero = j0(argument);
        let one = j1(argument);
        let ratio = j1_over_argument(argument, one);
        let derivative = zero - ratio;
        let second_ratio = j2_over_argument(argument, zero, ratio);
        let third = one - second_ratio;
        let pressure = kernel.pressure_scale;
        let velocity = kernel.velocity_scale;
        let angular = [
            derivative,
            ratio,
            one,
            zero,
            third,
            second_ratio,
            derivative,
            one,
            second_ratio,
            ratio,
            derivative,
            ratio,
            one,
            zero,
        ];
        if attenuation <= 700.0 {
            let weight = phasor.scale(det::exp(-attenuation))
                * (C64::from_re(measure[0]) / denominator).scale(measure[1]);
            let first = weight * zeta[0];
            let second = weight * zeta[1];
            let first_squared = first * zeta[0];
            let values = [
                minus_i(weight).scale(((pressure * q) * q) * derivative),
                minus_i(weight).scale(((pressure * q) * q) * ratio),
                first.scale((pressure * q) * (SQRT_2 * one)),
                minus_i(first_squared).scale(pressure * zero),
                weight.scale((((velocity * q) * q) * q) * third),
                weight.scale((((velocity * q) * q) * q) * second_ratio),
                minus_i(first).scale(((velocity * q) * q) * (SQRT_2 * derivative)),
                first_squared.scale((velocity * q) * one),
                weight.scale((((velocity * q) * q) * q) * (SQRT_2 * second_ratio)),
                minus_i(first).scale(((velocity * q) * q) * (SQRT_2 * ratio)),
                minus_i(second).scale(((velocity * q) * q) * derivative),
                minus_i(second).scale(((velocity * q) * q) * ratio),
                (first * zeta[1]).scale((velocity * q) * (SQRT_2 * one)),
                minus_i(first_squared * zeta[1]).scale(velocity * zero),
            ];
            if values
                .iter()
                .zip(angular)
                .all(|(v, factor)| finite(*v) && (factor == 0.0 || *v != C64::ZERO))
            {
                return Ok(values);
            }
        }
        // Balance physical derivatives before exponentiation: q³ or zeta1²zeta2
        // can restore a representable field after raw decay underflows, or can
        // overflow separately even though the complete attenuated field is finite.
        let scale = denominator.re.abs().max(denominator.im.abs());
        if !positive(scale) {
            return Err(nonfinite("full-wave quadrupole interface denominator"));
        }
        let normalized = C64::new(denominator.re / scale, denominator.im / scale);
        let magnitude = normalized.abs();
        let unit = phasor * normalized.conj().scale(1.0 / magnitude);
        let log_weight = det::ln(measure[0]) + det::ln(measure[1])
            - det::ln(scale)
            - det::ln(magnitude)
            - attenuation;
        let zmag = zeta.map(C64::abs);
        let zunit = core::array::from_fn::<_, 2, _>(|i| {
            if zmag[i] == 0.0 {
                C64::ZERO
            } else {
                C64::new(zeta[i].re / zmag[i], zeta[i].im / zmag[i])
            }
        });
        let lp = log_weight + det::ln(pressure);
        let lv = log_weight + det::ln(velocity);
        let lq = det::ln(q);
        let lz = zmag.map(det::ln);
        let first = unit * zunit[0];
        let second = unit * zunit[1];
        let first_squared = first * zunit[0];
        let values = [
            log_weighted(minus_i(unit), lp + 2.0 * lq, derivative),
            log_weighted(minus_i(unit), lp + 2.0 * lq, ratio),
            log_weighted(first, lp + lq + lz[0], SQRT_2 * one),
            log_weighted(minus_i(first_squared), lp + 2.0 * lz[0], zero),
            log_weighted(unit, lv + 3.0 * lq, third),
            log_weighted(unit, lv + 3.0 * lq, second_ratio),
            log_weighted(minus_i(first), lv + 2.0 * lq + lz[0], SQRT_2 * derivative),
            log_weighted(first_squared, lv + lq + 2.0 * lz[0], one),
            log_weighted(unit, lv + 3.0 * lq, SQRT_2 * second_ratio),
            log_weighted(minus_i(first), lv + 2.0 * lq + lz[0], SQRT_2 * ratio),
            log_weighted(minus_i(second), lv + 2.0 * lq + lz[1], derivative),
            log_weighted(minus_i(second), lv + 2.0 * lq + lz[1], ratio),
            log_weighted(first * zunit[1], lv + lq + lz[0] + lz[1], SQRT_2 * one),
            log_weighted(
                minus_i(first_squared * zunit[1]),
                lv + 2.0 * lz[0] + lz[1],
                zero,
            ),
        ];
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave quadrupole spectral field"));
        }
        Ok(values)
    }
}

pub(super) fn j2_over_argument(argument: f64, zero: f64, ratio: f64) -> f64 {
    if argument.abs() < 0.1 {
        // DLMF 10.2.2: the omitted x^11 term is <4.1e-21 in this branch.
        let square = argument * argument;
        argument
            * (1.0 / 8.0
                - square
                    * (1.0 / 96.0
                        - square
                            * (1.0 / 3072.0 - square * (1.0 / 184_320.0 - square / 17_694_720.0))))
    } else {
        (2.0 * ratio - zero) / argument
    }
}

impl QuadrupoleKernel {
    fn prepare_tail(
        &mut self,
        limits: FullWaveQuadrupoleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 14], AeroacError> {
        let kernel = &mut self.0;
        // q²=kmax²+t² gives q*dq/|D| <= dt/(d1+d2). Each q or |zeta|
        // is <=t+kmax and |E|<=exp(-h*t). Thus pressure has the degree-two
        // positive envelope and velocity degree three. Angular moments have
        // magnitude <=1; effective channel counts include squared sqrt(2)
        // shear weights, so their norms control every physical tensor entry.
        let log_density = det::ln(kernel.density[0] + kernel.density[1]);
        let (extent, [p, v]) = prepare_reflected_tail(
            (kernel.a + kernel.b).to_f64(),
            kernel.kmax,
            limits.into(),
            [2, 3],
            [
                det::ln(kernel.pressure_scale) - log_density,
                det::ln(kernel.velocity_scale) - log_density,
            ],
            [5, 14],
            cancelled,
        )?;
        kernel.tail_extent = extent;
        let odd_p = if kernel.range == 0.0 { 0.0 } else { p };
        let odd_v = if kernel.range == 0.0 { 0.0 } else { v };
        Ok([
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
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_third_derivatives_survive_raw_decay_underflow_against_angular_plane_waves() {
        let scale = 1e100;
        let kernel = QuadrupoleKernel(Kernel {
            a: Dd::from_f64(0.5 / scale),
            b: Dd::from_f64(0.5 / scale),
            range: 1e-103,
            density: [1e-200, 1.0],
            kmin: scale,
            kmax: 2.0 * scale,
            gap: det::sqrt(3.0) * scale,
            incident_has_larger_k: true,
            pressure_scale: 1.0 / TAU,
            velocity_scale: 1e-150,
            tail_extent: 1000.0 * scale,
        });
        let point = kernel.0.point(2, 0.8);
        let attenuation = (kernel.0.a * Dd::from_f64(point.zeta[0].im)
            + kernel.0.b * Dd::from_f64(point.zeta[1].im))
        .to_f64();
        assert_eq!(det::exp(-attenuation), 0.0);
        let half = det::exp(-0.5 * attenuation);
        let denominator =
            kernel.0.density[1] * point.zeta[0].im + kernel.0.density[0] * point.zeta[1].im;
        let measure = (point.measure[0] / denominator) * (point.measure[1] * half);
        let mut expected = [C64::ZERO; 14];
        // Independent angular plane-wave quadrature, with no Bessel function,
        // recurrence, logarithmic reconstruction or cylindrical derivative
        // formula. Two half-decays balance the original dimensional products.
        for angle in 0..128 {
            let phi = TAU * (angle as f64 + 0.5) / 128.0;
            let (sine, cosine) = phi.sin_cos();
            let source_wave = [
                C64::from_re(point.q * cosine),
                C64::from_re(point.q * sine),
                point.zeta[0],
            ];
            let receiver_wave = [source_wave[0], source_wave[1], point.zeta[1]];
            let (sine, cosine) = (point.q * kernel.0.range * cosine).sin_cos();
            let angular = C64::new(cosine, sine).scale(1.0 / 128.0);
            for (channel, (i, j, normalization)) in
                [(0, 0, 1.0), (1, 1, 1.0), (0, 2, SQRT_2), (2, 2, 1.0)]
                    .into_iter()
                    .enumerate()
            {
                let value = (-source_wave[i].scale(half) * source_wave[j])
                    .scale(kernel.0.pressure_scale * measure)
                    * angular;
                expected[channel] = expected[channel] + value.scale(normalization);
            }
            for (channel, (l, i, j, normalization)) in [
                (0, 0, 0, 1.0),
                (0, 1, 1, 1.0),
                (0, 0, 2, SQRT_2),
                (0, 2, 2, 1.0),
                (1, 0, 1, SQRT_2),
                (1, 1, 2, SQRT_2),
                (2, 0, 0, 1.0),
                (2, 1, 1, 1.0),
                (2, 0, 2, SQRT_2),
                (2, 2, 2, 1.0),
            ]
            .into_iter()
            .enumerate()
            {
                let value = ((-source_wave[i].scale(half)
                    * source_wave[j].scale(kernel.0.velocity_scale))
                    * receiver_wave[l])
                    .scale(measure)
                    * angular;
                expected[channel + 4] = expected[channel + 4] + value.scale(normalization);
            }
        }
        for (actual, reference) in kernel.sample(2, 0.8).unwrap().into_iter().zip(expected) {
            assert!(actual.abs() > 0.0 && reference.abs() > 0.0);
            assert!(
                (actual - reference).abs() <= 5e-11 * reference.abs(),
                "actual={actual:?}, angular-plane-wave reference={reference:?}"
            );
        }
    }
}
