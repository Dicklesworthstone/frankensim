//! Passive thin Jones layers on the original physical entrance pupil.
//!
//! J(u) = R(theta) diag(t0,t1) R(theta)^T, theta = theta0 + g.u,
//! with |tj|² in [0,1]. Writing p=(t0+t1)/2 and d=(t0-t1)/2 gives
//! J = p I + d [[cos(2 theta), sin(2 theta)], [sin(2 theta), -cos(2 theta)]].
//! Thus each patterned layer has at most three spatial Fourier modes. Products
//! of coplanar layers convolve those modes, retaining their matrix order.
//! A mode exp(i q.u) uses the ORIGINAL pupil at theta_obs-theta_source-q/k.
//! No pupil raster, separate diffraction integrator or efficiency fitting is used.
//!
//! This is a bounded thin-screen Jones model, not inter-layer propagation,
//! multiple reflection, arbitrary spatial textures or high-NA vector optics.
//! Principal transmission and phase are authored at EACH admitted frequency;
//! their common phase is retained for spectral interference. See Doelman et al.,
//! *Minimizing the polarization leakage of geometric-phase coronagraphs with
//! multiple grating pattern combinations*: <https://arxiv.org/abs/2002.05622>.

use super::super::super::super::PreparedDiffractionAmplitude;
use super::{
    CoherentDiffractionError, CoherentJonesPupilWave, Dd, PreparedCoherentJonesDiffractionExposure,
    add_count, admit_coordinate, admit_count, limit, poll, rounded,
};
use fs_math::det;

type Complex = [Dd; 2];
type Matrix = [[Complex; 2]; 2];
const ZERO: Complex = [Dd::ZERO; 2];
const ONE: Complex = [Dd::ONE, Dd::ZERO];
const ZERO_MATRIX: Matrix = [[ZERO; 2]; 2];
const IDENTITY: Matrix = [[ONE, ZERO], [ZERO, ONE]];

/// Hard bound on each intermediate frequency's distinct Jones-pupil modes.
pub const MAX_COHERENT_JONES_PUPIL_MODES: usize = 64;

/// A passive linear retarder/diattenuator with spatially rotating principal axes.
/// All layers in a response occupy the SAME pupil plane; list order is optical order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentJonesPupilLayer {
    power_transmission: [f64; 2],
    phase_rad: [f64; 2],
    axis_angle_rad: f64,
    axis_gradient_rad_m: [f64; 2],
}

impl CoherentJonesPupilLayer {
    /// Each principal power transmission lies in [0,1], with field response
    /// sqrt(Tj) exp(i phase_j). All parameters must be finite. The first axis
    /// is at `axis_angle_rad + dot(axis_gradient_rad_m, pupil_position_m)` from
    /// pupil x. A gradient pi/period produces a polarization grating of that
    /// physical period. Zero gradient is a uniform waveplate or diattenuator.
    pub fn try_new(
        power_transmission: [f64; 2],
        phase_rad: [f64; 2],
        axis_angle_rad: f64,
        axis_gradient_rad_m: [f64; 2],
    ) -> Result<Self, CoherentDiffractionError> {
        if power_transmission
            .iter()
            .any(|&t| !t.is_finite() || !(0.0..=1.0).contains(&t))
            || phase_rad
                .iter()
                .chain(&axis_gradient_rad_m)
                .any(|v| !v.is_finite())
            || !axis_angle_rad.is_finite()
        {
            return Err(CoherentDiffractionError::InvalidJonesPupil);
        }
        Ok(Self {
            power_transmission,
            phase_rad,
            axis_angle_rad,
            axis_gradient_rad_m,
        })
    }

    fn harmonics(self) -> Result<([Harmonic; 3], usize), CoherentDiffractionError> {
        let t: [Complex; 2] = std::array::from_fn(|axis| {
            let amplitude = Dd::from_f64(self.power_transmission[axis]).sqrt();
            let (s, c) = det::sin_cos_dd(Dd::from_f64(self.phase_rad[axis]));
            [amplitude * c, amplitude * s]
        });
        let (s, c) = det::sin_cos_dd(Dd::from_f64(self.axis_angle_rad));
        let mut modes = [Harmonic {
            gradient: [Dd::ZERO; 2],
            matrix: ZERO_MATRIX,
        }; 3];
        if self.axis_gradient_rad_m == [0.0; 2] {
            // Direct principal-channel products retain a very weak aligned
            // transmission instead of subtracting rounded mean/difference terms.
            let cross = scale(sub(t[0], t[1]), s * c);
            modes[0].matrix = [
                [add(scale(t[0], c * c), scale(t[1], s * s)), cross],
                [cross, add(scale(t[0], s * s), scale(t[1], c * c))],
            ];
            return Ok((modes, 1));
        }
        let half = Dd::from_f64(0.5);
        let mean = scale(add(t[0], t[1]), half);
        let difference = scale(sub(t[0], t[1]), half);
        modes[0].matrix = [[mean, ZERO], [ZERO, mean]];
        // Double-angle identities avoid overflowing 2*axis_angle for a finite angle.
        let phase = [c * c - s * s, Dd::from_f64(2.0) * s * c];
        let plus = scale(mul(difference, phase), half);
        let minus = scale(mul(difference, [phase[0], -phase[1]]), half);
        let neg_i_plus = [plus[1], -plus[0]];
        let i_minus = [-minus[1], minus[0]];
        let gradient = self
            .axis_gradient_rad_m
            .map(|g| Dd::from_f64(g) * Dd::from_f64(2.0));
        for q in gradient {
            rounded(q)?;
        }
        modes[1] = Harmonic {
            gradient,
            matrix: [[plus, neg_i_plus], [neg_i_plus, neg(plus)]],
        };
        modes[2] = Harmonic {
            gradient: gradient.map(|q| -q),
            matrix: [[minus, i_minus], [i_minus, neg(minus)]],
        };
        Ok((modes, 3))
    }
}

/// Explicit pupil response for one exact vacuum frequency, without interpolation.
/// Repeated source tones at that frequency share this SAME physical response.
#[derive(Clone, Copy, Debug)]
pub struct CoherentJonesPupilResponse<'a> {
    pub frequency_hz: f64,
    /// Empty layers mean the original clear pupil. First layer acts first.
    pub layers: &'a [CoherentJonesPupilLayer],
}

/// Global preparation bounds, checked before publishing a changed pupil.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoherentJonesPupilLimits {
    /// Total authored layers over all distinct frequencies.
    pub maximum_layers: usize,
    /// Total modes in completed frequencies plus the active intermediate
    /// expansion, including its identity seed. Exact duplicate gradients merge;
    /// no small coefficients are thresholded away to satisfy this budget.
    pub maximum_modes: usize,
    /// Total original pupil evaluations for ONE point, over all waves/factors.
    /// Pixel limits additionally bound this count times the spatial sample count.
    pub maximum_pupil_evaluations: usize,
}

#[derive(Clone, Copy)]
struct Harmonic {
    gradient: [Dd; 2], // spatial phase gradient in rad/m, BEFORE medium conversion
    matrix: Matrix,
}

struct PupilMode {
    tangent_shift: [Dd; 2],
    matrix: Matrix,
}

pub(super) struct PreparedJonesPupil {
    modes: Vec<PupilMode>,
}

impl PreparedCoherentJonesDiffractionExposure<'_, '_> {
    /// Install a complete frequency-indexed set of passive, coplanar pupil layers.
    ///
    /// Each distinct source frequency requires exactly one response. This
    /// prevents equal-frequency source entries from acquiring different optics.
    /// Spectral transmission/retardance and any material validity remain authored;
    /// no missing frequency is extrapolated. Full common phase is retained.
    /// Layers compile to exact finite Fourier/Jones modes of the original pupil;
    /// no truncated spatial series or separate pupil quadrature is introduced.
    ///
    /// Source amplitudes keep their original CLEAR-pupil normalization. Passive
    /// loss is not renormalized away. Shifted incident directions must stay in
    /// the original paraxial domain. Cancellation, bad data or any exceeded
    /// bound preserves the previously installed pupil and its work counts.
    pub fn set_pupil_responses(
        &mut self,
        responses: &[CoherentJonesPupilResponse<'_>],
        limits: CoherentJonesPupilLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), CoherentDiffractionError> {
        poll(&mut cancelled)?;
        if responses.len() != self.frequency_count() {
            return Err(CoherentDiffractionError::InconsistentJonesPupilSpectrum);
        }
        let mut layer_count = 0;
        for (index, response) in responses.iter().enumerate() {
            poll(&mut cancelled)?;
            if !response.frequency_hz.is_finite() || response.frequency_hz <= 0.0 {
                return Err(CoherentDiffractionError::InvalidFrequency);
            }
            #[allow(clippy::float_cmp)]
            if responses[..index]
                .iter()
                .any(|r| r.frequency_hz == response.frequency_hz)
            {
                return Err(CoherentDiffractionError::InconsistentJonesPupilSpectrum);
            }
            layer_count = add_count(
                layer_count,
                response.layers.len(),
                "Jones pupil layers",
                limits.maximum_layers,
            )?;
        }
        let mut staged = Vec::with_capacity(self.frequency_count());
        let mut total_modes = 0;
        let mut evaluations = 0;
        for (group_index, group) in self.optics.groups.iter().enumerate() {
            poll(&mut cancelled)?;
            #[allow(clippy::float_cmp)]
            let response = responses
                .iter()
                .find(|r| r.frequency_hz == group.frequency_hz)
                .ok_or(CoherentDiffractionError::InconsistentJonesPupilSpectrum)?;
            let maximum = limits
                .maximum_modes
                .saturating_sub(total_modes)
                .min(MAX_COHERENT_JONES_PUPIL_MODES);
            let harmonics = compile(response.layers, maximum, &mut cancelled)?;
            total_modes = add_count(
                total_modes,
                harmonics.len(),
                "Jones pupil modes",
                limits.maximum_modes,
            )?;
            let first = self
                .optics
                .tone_groups
                .iter()
                .position(|&g| g == group_index)
                .ok_or(CoherentDiffractionError::InconsistentJonesPupilSpectrum)?;
            // Original scale is lambda_medium/(2*pi*radius), so scale*radius
            // is 1/k_medium. Never substitute vacuum wavelength in a medium.
            let inverse_k = Dd::from_f64(self.optics.tones[first].pupil.scale)
                * Dd::from_f64(self.pupil_radius_m);
            let mut modes = Vec::with_capacity(harmonics.len());
            for harmonic in harmonics {
                poll(&mut cancelled)?;
                let tangent_shift = harmonic.gradient.map(|q| q * inverse_k);
                for shift in tangent_shift {
                    rounded(shift)?;
                }
                modes.push(PupilMode {
                    tangent_shift,
                    matrix: harmonic.matrix,
                });
            }
            for (source, &source_group) in self.sources.iter().zip(&self.optics.tone_groups) {
                if source_group != group_index {
                    continue;
                }
                for bank in
                    core::iter::once(source.mean).chain(source.factors.iter().map(|f| f.waves))
                {
                    poll(&mut cancelled)?;
                    let requested = bank.len().checked_mul(modes.len()).ok_or_else(|| {
                        limit("Jones pupil evaluations", limits.maximum_pupil_evaluations)
                    })?;
                    evaluations = add_count(
                        evaluations,
                        requested,
                        "Jones pupil evaluations",
                        limits.maximum_pupil_evaluations,
                    )?;
                    for wave in bank {
                        for mode in &modes {
                            poll(&mut cancelled)?;
                            admit_coordinate([
                                rounded(
                                    Dd::from_f64(wave.tangent_direction[0]) + mode.tangent_shift[0],
                                )?,
                                rounded(
                                    Dd::from_f64(wave.tangent_direction[1]) + mode.tangent_shift[1],
                                )?,
                            ])?;
                        }
                    }
                }
            }
            staged.push(PreparedJonesPupil { modes });
        }
        poll(&mut cancelled)?;
        self.pupils = Some(staged);
        self.pupil_evaluations = evaluations;
        Ok(())
    }
}

impl PreparedJonesPupil {
    pub(super) fn field(
        &self,
        pupil: &PreparedDiffractionAmplitude<'_>,
        wave: &CoherentJonesPupilWave,
        at: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Complex; 2], CoherentDiffractionError> {
        let input = wave
            .amplitude_sqrt_w
            .map(|a| [Dd::from_f64(a.re), Dd::from_f64(a.im)]);
        let offset = std::array::from_fn::<_, 2, _>(|axis| {
            Dd::from_f64(at[axis]) - Dd::from_f64(wave.tangent_direction[axis])
        });
        let mut sum = [ZERO; 2];
        for mode in &self.modes {
            poll(cancelled)?;
            let h = pupil.amplitude([
                rounded(offset[0] - mode.tangent_shift[0])?,
                rounded(offset[1] - mode.tangent_shift[1])?,
            ])?;
            let transfer = [Dd::from_f64(h.re), Dd::from_f64(h.im)];
            for (target, row) in sum.iter_mut().zip(&mode.matrix) {
                let field = add(
                    checked_mul(row[0], input[0])?,
                    checked_mul(row[1], input[1])?,
                );
                *target = add(*target, checked_mul(transfer, field)?);
            }
        }
        for value in sum.iter().flatten() {
            rounded(*value)?;
        }
        Ok(sum)
    }
}

fn compile(
    layers: &[CoherentJonesPupilLayer],
    maximum_modes: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Harmonic>, CoherentDiffractionError> {
    admit_count("Jones pupil modes per frequency", 1, maximum_modes)?;
    let mut modes = vec![Harmonic {
        gradient: [Dd::ZERO; 2],
        matrix: IDENTITY,
    }];
    for layer in layers {
        poll(cancelled)?;
        let (local, count) = layer.harmonics()?;
        let mut next: Vec<Harmonic> = Vec::with_capacity((modes.len() * count).min(maximum_modes));
        for old in &modes {
            for added in &local[..count] {
                poll(cancelled)?;
                let matrix = matrix_product(added.matrix, old.matrix)?;
                if matrix.iter().flatten().flatten().all(|v| *v == Dd::ZERO) {
                    continue;
                }
                let gradient =
                    std::array::from_fn(|axis| old.gradient[axis] + added.gradient[axis]);
                for q in gradient {
                    rounded(q)?;
                }
                if let Some(existing) = next.iter_mut().find(|term| term.gradient == gradient) {
                    for (a, b) in existing
                        .matrix
                        .iter_mut()
                        .flatten()
                        .zip(matrix.into_iter().flatten())
                    {
                        *a = add(*a, b);
                    }
                } else {
                    admit_count(
                        "Jones pupil modes per frequency",
                        next.len() + 1,
                        maximum_modes,
                    )?;
                    next.push(Harmonic { gradient, matrix });
                }
            }
        }
        // Remove only exact cancellation. Small physical leakage is retained.
        next.retain(|mode| {
            mode.matrix
                .iter()
                .flatten()
                .flatten()
                .any(|v| *v != Dd::ZERO)
        });
        modes = next;
    }
    Ok(modes)
}

fn add(a: Complex, b: Complex) -> Complex {
    [a[0] + b[0], a[1] + b[1]]
}
fn sub(a: Complex, b: Complex) -> Complex {
    [a[0] - b[0], a[1] - b[1]]
}
fn neg(a: Complex) -> Complex {
    [-a[0], -a[1]]
}
fn scale(a: Complex, b: Dd) -> Complex {
    [a[0] * b, a[1] * b]
}
fn mul(a: Complex, b: Complex) -> Complex {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}
fn checked_mul(a: Complex, b: Complex) -> Result<Complex, CoherentDiffractionError> {
    let product = mul(a, b);
    if product == ZERO && a != ZERO && b != ZERO {
        return Err(CoherentDiffractionError::Unrepresentable);
    }
    Ok(product)
}
fn matrix_product(a: Matrix, b: Matrix) -> Result<Matrix, CoherentDiffractionError> {
    let mut product = ZERO_MATRIX;
    for r in 0..2 {
        for c in 0..2 {
            product[r][c] = add(
                checked_mul(a[r][0], b[0][c])?,
                checked_mul(a[r][1], b[1][c])?,
            );
        }
    }
    Ok(product)
}
