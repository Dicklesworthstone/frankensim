//! Analytic source and receiver derivatives of the transmitted Green function.
//!
//! For a prescribed force ON the incident fluid, `p_F = F . grad_source H`,
//! where `p_Q = -i omega rho1 H Q`. Thus source-normal differentiation brings
//! in zeta1, whereas receiver-normal differentiation brings in zeta2. Replacing
//! either by a refracted direction or p/Z loses the near-interface tensor.
//!
//! In the radial/transverse/normal frame, seven Sommerfeld integrals determine
//! the pressure vector and full velocity matrix. With w=q*dq*E/D, their
//! integrands are rho2/(2*pi)*[i*q*J1, zeta1*J0] for pressure and
//! 1/(2*pi*omega)*[q²*J1', q²*J1/x, i*q*zeta1*J1,
//! i*q*zeta2*J1, zeta1*zeta2*J0] for velocity (tt,ss,tn,nt,nn).
//! `x=qR`, `J1'=J0-J1/x`, and `J1/x -> 1/2` on axis. See DLMF
//! 10.6.2 and 10.9.2/10.9.4: <https://dlmf.nist.gov/10.6>,
//! <https://dlmf.nist.gov/10.9>. The latter integral forms also give the
//! real-argument bounds |J0|,|J1|,|J1'|<=1 and |J1/x|<=1/2 used below.

use super::{
    ComplexSum, Dd, Kernel, PlanarFluidInterface, QuadratureLimits, SpectralPoint, Spectrum,
    finite, integrate_spectrum, invalid, j0, j1, log_weighted, minus_i, nonfinite, norm, poll,
    positive,
};
use crate::AeroacError;
use crate::curle3d::CoherentCompactField;
use crate::observation::HarmonicFieldLine;
use core::f64::consts::{PI, TAU};
use fs_math::{c64::C64, det};

/// Applicability of the physical force-to-field transfer.
pub const FULL_WAVE_FORCE_SCOPE: &str = "Prescribed stationary compact 3D force ON fluid [N]; \
outgoing spherical-wave transmission across one infinite planar boundary between two \
stationary homogeneous lossless fluids, including transverse, reactive and evanescent fields. \
Estimated finite quadrature error and a separate analytic tail bound; no interval \
certificate, flow, absorption, finite boundaries, source feedback or measured-SPL qualification.";

/// Numerical accuracy of the force-to-pressure vector and force-to-velocity
/// tensor, before multiplying by a particular complex force vector [N].
#[derive(Clone, Copy, Debug)]
pub struct FullWaveForceLimits {
    /// Positive absolute tolerance on the Euclidean pressure-transfer norm [Pa/N].
    pub absolute_pressure_tolerance_pa_per_n: f64,
    /// Positive absolute tolerance on the Frobenius velocity-transfer norm [m/s/N].
    /// This retains both transverse diagonal terms, including on the normal axis.
    pub absolute_velocity_tolerance_m_s_per_n: f64,
    /// Relative tolerance in [0,1), combined with each absolute tolerance by max.
    /// This controls the transfer norm, not the relative error at a force null.
    pub relative_tolerance: f64,
    /// Integrand evaluations in 1..=MAX_FULL_WAVE_EVALUATIONS, shared by all
    /// seven complex channels. The exact equal-speed expression needs none.
    pub maximum_evaluations: usize,
    /// Transformed-region bisection depth in 1..=48; root depth is zero.
    pub maximum_depth: usize,
}

impl From<FullWaveForceLimits> for QuadratureLimits {
    fn from(limits: FullWaveForceLimits) -> Self {
        Self {
            absolute_pressure: limits.absolute_pressure_tolerance_pa_per_n,
            absolute_velocity: limits.absolute_velocity_tolerance_m_s_per_n,
            relative: limits.relative_tolerance,
            maximum_evaluations: limits.maximum_evaluations,
            maximum_depth: limits.maximum_depth,
        }
    }
}

pub(super) fn admit(limits: FullWaveForceLimits) -> Result<(), AeroacError> {
    QuadratureLimits::from(limits).admit()
}

/// Numerical estimates and analytic tail bounds for unit force transfers.
/// Pressure entries are (radial,normal); velocity entries are (tt,ss,tn,nt,nn),
/// where rows are velocity directions and columns are force directions.
/// Euclidean norms of these arrays are respectively pressure-vector and
/// velocity-tensor errors; multiply by ||F||_2 for a particular force [N].
/// Quadrature estimates are not certificates. Neither estimate nor tail bound
/// includes geometry, Bessel, coordinate transformation or floating-point error.
#[derive(Clone, Copy, Debug)]
pub struct FullWaveForceDiagnostics {
    pub evaluations: usize,
    pub maximum_depth_reached: usize,
    pub retained_panels: usize,
    /// Retained upper transverse wavenumber [rad/m], zero for the exact case.
    pub transverse_wavenumber_cutoff_rad_m: f64,
    /// Largest q*R; the native Bessel accuracy scope is stated by the parent.
    pub maximum_bessel_argument: f64,
    pub pressure_quadrature_error_pa_per_n: [f64; 2],
    pub velocity_quadrature_error_m_s_per_n: [f64; 5],
    pub pressure_tail_bound_pa_per_n: [f64; 2],
    pub velocity_tail_bound_m_s_per_n: [f64; 5],
    /// True only for exactly equal ORIGINAL sound speeds, including density jumps.
    pub exact_equal_speed: bool,
}

/// Complete transmitted pressure and velocity for any prescribed complex force
/// at one source position and frequency. Preparation caches twelve complex
/// coefficients. Application has fixed work, DD coherent reductions, no
/// allocation, quadrature, finite differences or transcendental evaluation.
#[derive(Clone, Debug)]
pub struct PreparedFullWaveForce {
    interface: PlanarFluidInterface,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    // Rows: p,ux,uy,uz. Columns: original Cartesian Fx,Fy,Fz [N].
    pub(super) transfer: [[C64; 3]; 4],
    diagnostics: FullWaveForceDiagnostics,
}

impl PreparedFullWaveForce {
    /// Prepare a stationary force on the incident fluid and a receiver strictly
    /// in the receiving fluid. Every propagating, critical and evanescent
    /// component is retained up to explicit quadrature and tail controls.
    pub fn new(
        interface: PlanarFluidInterface,
        frequency_hz: f64,
        source_m: [f64; 3],
        observer_m: [f64; 3],
        limits: FullWaveForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        admit(limits)?;
        let omega = TAU * frequency_hz;
        if !positive(omega) {
            return Err(invalid("full-wave force frequency"));
        }
        if source_m.iter().chain(&observer_m).any(|v| !v.is_finite()) {
            return Err(nonfinite("full-wave force positions"));
        }
        let a = -interface.distance(source_m);
        let b = interface.distance(observer_m);
        let h = (a + b).to_f64();
        if !positive(a.to_f64()) || !positive(b.to_f64()) || !positive(h) {
            return Err(invalid(
                "full-wave force and observer must lie strictly on their fluid sides",
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
            return Err(nonfinite("full-wave force lateral separation"));
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
            return Err(nonfinite("full-wave force density contrast"));
        }
        let (field, diagnostics) = if c1 == c2 {
            (
                equal_speed(omega, c1, range, h, density_scale, density)?,
                FullWaveForceDiagnostics {
                    evaluations: 0,
                    maximum_depth_reached: 0,
                    retained_panels: 0,
                    transverse_wavenumber_cutoff_rad_m: 0.0,
                    maximum_bessel_argument: 0.0,
                    pressure_quadrature_error_pa_per_n: [0.0; 2],
                    velocity_quadrature_error_m_s_per_n: [0.0; 5],
                    pressure_tail_bound_pa_per_n: [0.0; 2],
                    velocity_tail_bound_m_s_per_n: [0.0; 5],
                    exact_equal_speed: true,
                },
            )
        } else {
            let k = [omega / c1, omega / c2];
            if k.iter().any(|v| !positive(*v)) || k[0] == k[1] {
                return Err(nonfinite("full-wave force distinct fluid wavenumbers"));
            }
            let kmin = k[0].min(k[1]);
            let kmax = k[0].max(k[1]);
            let gap = det::sqrt(kmax - kmin) * det::sqrt(kmax) * det::sqrt(1.0 + kmin / kmax);
            let pressure_scale = density[1] / TAU;
            let velocity_scale = ((1.0 / TAU) / omega) / density_scale;
            if !positive(gap) || !positive(pressure_scale) || !positive(velocity_scale) {
                return Err(nonfinite("full-wave force spectral scales"));
            }
            let mut kernel = ForceKernel(Kernel {
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
            let result = integrate_spectrum::<7, 2>(&kernel, tail, limits.into(), cancelled)?;
            (
                result.value,
                FullWaveForceDiagnostics {
                    evaluations: result.evaluations,
                    maximum_depth_reached: result.maximum_depth_reached,
                    retained_panels: result.retained_panels,
                    transverse_wavenumber_cutoff_rad_m: result.cutoff,
                    maximum_bessel_argument: result.maximum_argument,
                    pressure_quadrature_error_pa_per_n: [result.error[0], result.error[1]],
                    velocity_quadrature_error_m_s_per_n: core::array::from_fn(|i| {
                        result.error[i + 2]
                    }),
                    pressure_tail_bound_pa_per_n: [tail[0], tail[1]],
                    velocity_tail_bound_m_s_per_n: core::array::from_fn(|i| tail[i + 2]),
                    exact_equal_speed: false,
                },
            )
        };
        let [pt, pn, vtt, vss, vtn, vnt, vnn] = field;
        let normal = interface.normal;
        let mut transfer = [[C64::ZERO; 3]; 4];
        for column in 0..3 {
            transfer[0][column] = pt.scale(tangent[column]) + pn.scale(normal[column]);
            for row in 0..3 {
                let nn = normal[row] * normal[column];
                let identity = if row == column { 1.0 } else { 0.0 };
                let mut sum = ComplexSum::default();
                sum.add(vss.scale(identity - nn))?;
                sum.add((vtt - vss).scale(tangent[row] * tangent[column]))?;
                sum.add(vtn.scale(tangent[row] * normal[column]))?;
                sum.add(vnt.scale(normal[row] * tangent[column]))?;
                sum.add(vnn.scale(nn))?;
                transfer[row + 1][column] = sum.value();
            }
        }
        if transfer.iter().flatten().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave force transfer"));
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

    pub const fn diagnostics(&self) -> &FullWaveForceDiagnostics {
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

    fn amplitudes(&self, force_on_fluid_n: [C64; 3]) -> Result<[C64; 4], AeroacError> {
        if force_on_fluid_n.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave force phasor"));
        }
        let mut field = [C64::ZERO; 4];
        for (row, coefficients) in self.transfer.iter().enumerate() {
            let mut sum = ComplexSum::default();
            for column in 0..3 {
                sum.add(coefficients[column] * force_on_fluid_n[column])?;
            }
            field[row] = sum.value();
        }
        if field.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave force field amplitude"));
        }
        Ok(field)
    }

    /// Complex peak force ON the original fluid [N], with the source clock
    /// exp(-i omega t). Intensity follows the coherent full vector field.
    pub fn evaluate(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<CoherentCompactField, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(force_on_fluid_n)?;
        let velocity = [ux, uy, uz];
        let intensity = velocity.map(|v| pressure.scale(0.5) * v.conj());
        if intensity.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave force intensity"));
        }
        Ok(CoherentCompactField {
            pressure_pa: pressure,
            particle_velocity_m_s: velocity,
            active_intensity_w_m2: intensity.map(|v| v.re),
            reactive_intensity_w_m2: intensity.map(|v| v.im),
            scope: FULL_WAVE_FORCE_SCOPE,
        })
    }

    /// Native pressure/vector-velocity observer line; all propagation phase is
    /// already included, so callers must not apply a second travel delay.
    pub fn evaluate_harmonic_line(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(force_on_fluid_n)?;
        Ok(HarmonicFieldLine {
            frequency_hz: self.frequency_hz,
            pressure_pa: pressure,
            particle_velocity_m_s: [ux, uy, uz],
        })
    }
}

fn equal_speed(
    omega: f64,
    speed: f64,
    range: f64,
    h: f64,
    density_scale: f64,
    density: [f64; 2],
) -> Result<[C64; 7], AeroacError> {
    let distance = det::hypot(range, h);
    let inverse = 1.0 / distance;
    let k = omega / speed;
    let phase = Dd::from_f64(omega) * Dd::from_f64(distance) / Dd::from_f64(speed);
    let pressure_scale = 2.0 * density[1] / (density[0] + density[1]);
    let velocity_scale = ((2.0 / (density[0] + density[1])) / density_scale) / omega;
    if !positive(distance)
        || !positive(inverse)
        || !positive(k)
        || !positive(velocity_scale)
        || !phase.hi.is_finite()
        || !phase.lo.is_finite()
    {
        return Err(nonfinite("equal-speed force distance, phase or scale"));
    }
    let (sin, cos) = det::sin_cos_dd(phase);
    let green = C64::new(cos.to_f64(), sin.to_f64()).scale((0.25 / PI) / distance);
    // a=-G'; retain BOTH its radial derivative and the angular derivative a/r.
    let radial = green * C64::new(inverse, -k);
    let derivative = green * C64::new(k * k - 2.0 * inverse * inverse, 2.0 * k * inverse);
    let tangential = radial.scale(inverse);
    let r = range / distance;
    let n = h / distance;
    let velocity = |value: C64| minus_i(value).scale(velocity_scale);
    let cross = velocity((derivative - tangential).scale(r * n));
    Ok([
        radial.scale(pressure_scale * r),
        radial.scale(pressure_scale * n),
        velocity(derivative.scale(r * r) + tangential.scale(1.0 - r * r)),
        velocity(tangential),
        cross,
        cross,
        velocity(derivative.scale(n * n) + tangential.scale(1.0 - n * n)),
    ])
}

struct ForceKernel(Kernel);

impl Spectrum<7> for ForceKernel {
    fn geometry(&self) -> &Kernel {
        &self.0
    }

    fn sample(&self, region: usize, x: f64) -> Result<[C64; 7], AeroacError> {
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
            return Err(nonfinite("full-wave force spectral phase"));
        }
        let (sin, cos) = det::sin_cos_dd(phase);
        let phasor = C64::new(cos.to_f64(), sin.to_f64());
        let denominator = zeta[0].scale(kernel.density[1]) + zeta[1].scale(kernel.density[0]);
        let jzero = j0(argument);
        let jone = j1(argument);
        let ratio = j1_over_argument(argument, jone);
        let derivative = jzero - ratio;
        let pressure = kernel.pressure_scale;
        let velocity = kernel.velocity_scale;
        if attenuation <= 700.0 {
            let weight = phasor.scale(det::exp(-attenuation))
                * (C64::from_re(measure[0]) / denominator).scale(measure[1]);
            let values = [
                times_i(weight).scale((pressure * q) * jone),
                (weight * zeta[0]).scale(pressure * jzero),
                weight.scale(((velocity * q) * q) * derivative),
                weight.scale(((velocity * q) * q) * ratio),
                times_i(weight * zeta[0]).scale((velocity * q) * jone),
                times_i(weight * zeta[1]).scale((velocity * q) * jone),
                ((weight * zeta[0]) * zeta[1]).scale(velocity * jzero),
            ];
            let nonzero = [
                jone != 0.0,
                zeta[0] != C64::ZERO && jzero != 0.0,
                derivative != 0.0,
                ratio != 0.0,
                zeta[0] != C64::ZERO && jone != 0.0,
                zeta[1] != C64::ZERO && jone != 0.0,
                zeta[0] != C64::ZERO && zeta[1] != C64::ZERO && jzero != 0.0,
            ];
            if values
                .iter()
                .zip(nonzero)
                .all(|(v, expected_nonzero)| finite(*v) && (!expected_nonzero || *v != C64::ZERO))
            {
                return Ok(values);
            }
        }
        // Preserve physically representable derivatives even when the raw
        // exponential or a q² product alone would underflow/overflow.
        let scale = denominator.re.abs().max(denominator.im.abs());
        if !positive(scale) {
            return Err(nonfinite("full-wave force interface denominator"));
        }
        let normalized = C64::new(denominator.re / scale, denominator.im / scale);
        let magnitude = normalized.abs();
        let unit = phasor * normalized.conj().scale(1.0 / magnitude);
        let log_weight = det::ln(measure[0]) + det::ln(measure[1])
            - det::ln(scale)
            - det::ln(magnitude)
            - attenuation;
        let zeta_magnitude = zeta.map(C64::abs);
        let zeta_unit = core::array::from_fn::<_, 2, _>(|i| {
            if zeta_magnitude[i] == 0.0 {
                C64::ZERO
            } else {
                C64::new(
                    zeta[i].re / zeta_magnitude[i],
                    zeta[i].im / zeta_magnitude[i],
                )
            }
        });
        let lp = log_weight + det::ln(pressure);
        let lv = log_weight + det::ln(velocity);
        let lq = det::ln(q);
        let lz = zeta_magnitude.map(det::ln);
        let values = [
            log_weighted(times_i(unit), lp + lq, jone),
            log_weighted(unit * zeta_unit[0], lp + lz[0], jzero),
            log_weighted(unit, lv + 2.0 * lq, derivative),
            log_weighted(unit, lv + 2.0 * lq, ratio),
            log_weighted(times_i(unit * zeta_unit[0]), lv + lq + lz[0], jone),
            log_weighted(times_i(unit * zeta_unit[1]), lv + lq + lz[1], jone),
            log_weighted(
                unit * zeta_unit[0] * zeta_unit[1],
                lv + lz[0] + lz[1],
                jzero,
            ),
        ];
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave force spectral field"));
        }
        Ok(values)
    }
}

fn times_i(value: C64) -> C64 {
    C64::new(-value.im, value.re)
}

fn j1_over_argument(argument: f64, jone: f64) -> f64 {
    if argument.abs() < 1e-3 {
        // Analytic continuation at zero; the omitted x^8 term is <2e-32.
        let square = argument * argument;
        0.5 - square * (1.0 / 16.0 - square * (1.0 / 384.0 - square / 18_432.0))
    } else {
        jone / argument
    }
}

impl ForceKernel {
    fn prepare_tail(
        &mut self,
        limits: FullWaveForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 7], AeroacError> {
        let kernel = &mut self.0;
        let h = (kernel.a + kernel.b).to_f64();
        let hk = h * kernel.kmax;
        if !hk.is_finite() {
            return Err(nonfinite("full-wave force clearance wavenumber"));
        }
        // For q²=kmax²+t² the base measure is <=e^(-ht)/(d1+d2) dt.
        // Every q and |zeta_i| is <=t+kmax. Integrating its degree-one and
        // degree-two powers gives the exact positive tail envelopes below.
        let log_density = det::ln(kernel.density[0] + kernel.density[1]);
        let log_h = det::ln(h);
        let log_pressure = det::ln(kernel.pressure_scale) - log_density - 2.0 * log_h;
        let log_velocity = det::ln(kernel.velocity_scale) - log_density - 3.0 * log_h;
        let goals = [
            det::ln(limits.absolute_pressure_tolerance_pa_per_n) - det::ln(4.0 * det::sqrt(2.0)),
            det::ln(limits.absolute_velocity_tolerance_m_s_per_n) - det::ln(4.0 * det::sqrt(5.0)),
        ];
        let logs = |s: f64| {
            let shifted = s + hk + 1.0;
            let inverse = 1.0 / shifted;
            [
                log_pressure - s + det::ln(shifted),
                log_velocity - s + 2.0 * det::ln(shifted) + det::ln(1.0 + inverse * inverse),
            ]
        };
        let mut low = 0.0;
        let mut high = 1.0;
        let mut admitted = false;
        for _ in 0..32 {
            poll(cancelled)?;
            let bound = logs(high);
            if bound[0] <= goals[0] && bound[1] <= goals[1] {
                admitted = true;
                break;
            }
            low = high;
            high *= 2.0;
        }
        if !admitted {
            return Err(invalid(
                "full-wave force tail cannot meet the requested tolerance",
            ));
        }
        for _ in 0..24 {
            poll(cancelled)?;
            let middle = low + (high - low) * 0.5;
            let bound = logs(middle);
            if bound[0] <= goals[0] && bound[1] <= goals[1] {
                high = middle;
            } else {
                low = middle;
            }
        }
        kernel.tail_extent = high / h;
        if !positive(kernel.tail_extent) {
            return Err(nonfinite("full-wave force evanescent cutoff"));
        }
        let bound = logs(high).map(det::exp);
        let off_axis = kernel.range != 0.0;
        Ok([
            if off_axis { bound[0] } else { 0.0 },
            bound[0],
            if off_axis { bound[1] } else { bound[1] * 0.5 },
            bound[1] * 0.5,
            if off_axis { bound[1] } else { 0.0 },
            if off_axis { bound[1] } else { 0.0 },
            bound[1],
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn force_derivatives_survive_underflow_of_the_unscaled_evanescent_exponential() {
        let scale = 1e160;
        let kernel = ForceKernel(Kernel {
            a: Dd::from_f64(0.5 / scale),
            b: Dd::from_f64(0.5 / scale),
            range: 1e-162,
            density: [1e-200, 1.0],
            kmin: scale,
            kmax: 2.0 * scale,
            gap: det::sqrt(3.0) * scale,
            incident_has_larger_k: true,
            pressure_scale: 1.0 / TAU,
            velocity_scale: 1e-250,
            tail_extent: 1000.0 * scale,
        });
        let point = kernel.0.point(2, 0.8);
        let attenuation = (kernel.0.a * Dd::from_f64(point.zeta[0].im)
            + kernel.0.b * Dd::from_f64(point.zeta[1].im))
        .to_f64();
        assert_eq!(det::exp(-attenuation), 0.0);
        // Independent balanced products use two half-decays and no logs.
        // This keeps every expected physical derivative representable.
        let half = det::exp(-0.5 * attenuation);
        let denominator =
            kernel.0.density[1] * point.zeta[0].im + kernel.0.density[0] * point.zeta[1].im;
        let measure = (point.measure[0] / denominator) * (point.measure[1] * half);
        let q = point.q;
        let [z1, z2] = point.zeta.map(|z| z.im);
        let argument = q * kernel.0.range;
        let zero = j0(argument);
        let one = j1(argument);
        let ratio = one / argument;
        let pressure = kernel.0.pressure_scale;
        let velocity = kernel.0.velocity_scale;
        let expected = [
            C64::from_re(measure * ((pressure * q) * half) * one),
            C64::from_re(measure * ((pressure * z1) * half) * zero),
            C64::new(
                0.0,
                -measure * (((velocity * q) * q) * half) * (zero - ratio),
            ),
            C64::new(0.0, -measure * (((velocity * q) * q) * half) * ratio),
            C64::new(0.0, measure * (((velocity * q) * z1) * half) * one),
            C64::new(0.0, measure * (((velocity * q) * z2) * half) * one),
            C64::new(0.0, measure * (((velocity * z1) * z2) * half) * zero),
        ];
        for (actual, reference) in kernel.sample(2, 0.8).unwrap().into_iter().zip(expected) {
            assert!(actual.abs() > 0.0 && reference.abs() > 0.0);
            assert!((actual - reference).abs() <= 3e-12 * reference.abs());
        }
    }
}
