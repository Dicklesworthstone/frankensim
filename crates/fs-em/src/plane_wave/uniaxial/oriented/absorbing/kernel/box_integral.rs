//! Original Maxwell-mode products on a rectangular control volume.
//!
//! For the native exp(+i omega t-i k.r) convention, conj(E_i) E_j is
//! weighted by exp(i (Re(k_i)-Re(k_j)).r) and positive decay. Transverse
//! integrals use the existing DD harmonic-window kernel; depth uses the
//! original complex absorption integral. Six face fluxes are integrated
//! independently of dielectric work. No flux residual is deposited as heat.
//! Finite films retain entrance/exit mode anchors and select a decaying
//! endpoint for every forward/backward product, including opaque films.
//!
//! Short boxes use positive tensor Gauss integration of the COHERENT field,
//! avoiding cancellation between nearly equal diagonal/cross integrals.
//! Equal normal propagation is integrated analytically even for opaque depth.
//! This is a bounded numerical evaluation, not an interval certificate.

use super::{
    AbsorbingCrystalResponse, IlluminationError, absorption, poll_illumination,
};
use super::super::{AbsorbingCrystalBox, AbsorbingCrystalBoxLimits, slab::AbsorbingSlabResponse};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_math::{c64::C64, dd::Dd, decay, det};

type Result<T> = core::result::Result<T, IlluminationError>;
type Complex = [Dd; 2];
type Fields = [Complex; 6];
const ZERO: Complex = [Dd::ZERO; 2];
const ONE: Complex = [Dd::ONE, Dd::ZERO];

#[derive(Clone, Copy)]
struct Wave {
    k: [Dd; 3],
    alpha: f64,
    // Positive attenuation away from the original entrance (+1) or exit (-1).
    orientation: f64,
    anchor: Dd,
    fields: Fields,
}

#[derive(Clone, Copy)]
struct Sample {
    fields: Fields,
    projected: [Complex; 4],
}

#[derive(Clone, Copy)]
struct Observation {
    origin: [Dd; 3],
    extent: [f64; 3],
    integrated: [bool; 3],
}

#[derive(Clone, Copy)]
enum Quantity {
    Heat([f64; 4]),
    Flux(usize),
}

struct Loss {
    axis: [Dd; 3],
    norm: Dd,
    weight: [f64; 4],
    omega: f64,
}

pub(in super::super) struct BoxIntegral {
    pub(in super::super) absorbed_power_w: f64,
    pub(in super::super) outward_face_power_w: [f64; 6],
    pub(in super::super) residual_w: f64,
    pub(in super::super) mode_pair_count: usize,
    pub(in super::super) workspace_bytes: usize,
}

/// Shared logical admission for an original coherent half-space factor.
pub(in super::super) fn coherent_box_requirements(
    component_count: usize,
    result_bytes: usize,
) -> Result<(usize, usize)> {
    let modes = component_count.checked_mul(2).ok_or(bad("crystal box mode count overflows"))?;
    box_requirements(modes, result_bytes)
}

/// A finite-film factor retains four incoming-face-anchored waves per beam.
pub(in super::super) fn slab_box_requirements(
    component_count: usize,
    result_bytes: usize,
) -> Result<(usize, usize)> {
    let modes = component_count.checked_mul(4).ok_or(bad("film box mode count overflows"))?;
    box_requirements(modes, result_bytes)
}

fn box_requirements(mode_count: usize, result_bytes: usize) -> Result<(usize, usize)> {
    let mode_pair_count = mode_count.checked_mul(mode_count).ok_or(bad("crystal box pair count overflows"))?;
    let workspace_bytes = mode_count
        .checked_mul(core::mem::size_of::<Wave>() + core::mem::size_of::<Sample>())
        .and_then(|n| n.checked_add(result_bytes))
        .and_then(|n| n.checked_add(core::mem::size_of::<BoxIntegral>()))
        .and_then(|n| n.checked_add(core::mem::size_of::<[Fields; 2]>()))
        .and_then(|n| n.checked_add(core::mem::size_of::<[[Dd; 8]; 2]>()))
        .ok_or(bad("crystal box workspace overflows"))?;
    Ok((mode_pair_count, workspace_bytes))
}

/// A rounded upper coordinate must not admit support outside the original film.
pub(in super::super) fn admit_slab_box(region: AbsorbingCrystalBox, thickness: f64) -> Result<()> {
    let upper = Dd::from_f64(region.origin_m()[2]) + Dd::from_f64(region.extent_m()[2]);
    if !thickness.is_finite() || thickness <= 0.0 || !finite(upper)
        || Dd::from_f64(thickness).lt(upper)
    {
        return Err(bad("film box extends beyond the original crystal thickness"));
    }
    Ok(())
}

fn original_loss(axis: [f64; 3], weights: [f64; 2], omega: f64) -> Loss {
    let axis = axis.map(Dd::from_f64);
    let norm = axis.iter().fold(Dd::ZERO, |sum, &v| sum + v * v);
    Loss {
        axis,
        norm,
        weight: [weights[0], weights[0], weights[0], weights[1] / norm.to_f64()],
        omega,
    }
}

pub(in super::super) fn coherent_box_absorption<'a>(
    responses: impl ExactSizeIterator<Item = &'a AbsorbingCrystalResponse> + Clone,
    region: AbsorbingCrystalBox,
    limits: AbsorbingCrystalBoxLimits,
    result_bytes: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<BoxIntegral> {
    poll_illumination(&mut cancelled)?;
    let component_count = responses.len();
    let (mode_pair_count, workspace_bytes) =
        coherent_box_requirements(component_count, result_bytes)?;
    // The shared requirements check admitted this exact multiplication.
    let mode_count = component_count * 2;
    if mode_pair_count > limits.maximum_mode_pairs || workspace_bytes > limits.maximum_workspace_bytes {
        return Err(bad("crystal box exceeds its pair work or workspace budget"));
    }
    let reference = responses.clone().next();
    let loss = original_loss(
        reference.map_or([0.0, 0.0, 1.0], |r| r.axis),
        reference.map_or([0.0; 2], |r| r.loss_f_m),
        reference.map_or(1.0, |r| r.omega),
    );
    let waves = coalesce(responses, mode_count, &mut cancelled)?;
    integrate_box(&waves, mode_count, loss, region, mode_pair_count, workspace_bytes, &mut cancelled)
}

pub(in super::super) fn coherent_slab_box_absorption<'a>(
    responses: impl ExactSizeIterator<Item = &'a AbsorbingSlabResponse> + Clone,
    region: AbsorbingCrystalBox,
    thickness: f64,
    limits: AbsorbingCrystalBoxLimits,
    result_bytes: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<BoxIntegral> {
    poll_illumination(&mut cancelled)?;
    admit_slab_box(region, thickness)?;
    let component_count = responses.len();
    let (mode_pair_count, workspace_bytes) = slab_box_requirements(component_count, result_bytes)?;
    if mode_pair_count > limits.maximum_mode_pairs || workspace_bytes > limits.maximum_workspace_bytes {
        return Err(bad("film box exceeds its pair work or workspace budget"));
    }
    let mode_count = component_count * 4; // checked by the shared admission
    let reference = responses.clone().next();
    let loss = original_loss(
        reference.map_or([0.0, 0.0, 1.0], |r| r.axis),
        reference.map_or([0.0; 2], |r| r.principal_loss_permittivities_f_m()),
        reference.map_or(1.0, |r| r.omega_rad_s()),
    );
    let waves = coalesce_slab(responses, mode_count, &mut cancelled)?;
    integrate_box(&waves, mode_count, loss, region, mode_pair_count, workspace_bytes, &mut cancelled)
}

fn integrate_box(
    waves: &[Wave],
    mode_count: usize,
    loss: Loss,
    region: AbsorbingCrystalBox,
    mode_pair_count: usize,
    workspace_bytes: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<BoxIntegral> {
    let mut samples = Vec::new();
    samples.try_reserve_exact(mode_count).map_err(|_| bad("crystal box field allocation failed"))?;
    let volume = Observation {
        origin: region.origin_m().map(Dd::from_f64),
        extent: region.extent_m(),
        integrated: [true; 3],
    };
    let absorbed_power_w = if loss.weight == [0.0; 4] {
        0.0
    } else {
        measure(waves, &mut samples, volume, Quantity::Heat(loss.weight), &loss, cancelled)?
    };
    let mut outward_face_power_w = [0.0; 6];
    for axis in 0..3 {
        for upper in 0..2 {
            poll_illumination(cancelled)?;
            let mut face = volume;
            face.integrated[axis] = false;
            if upper == 1 {
                face.origin[axis] = face.origin[axis] + Dd::from_f64(face.extent[axis]);
            }
            let watts = measure(
                waves, &mut samples, face, Quantity::Flux(axis), &loss, cancelled,
            )?;
            outward_face_power_w[2 * axis + upper] = if upper == 0 { -watts } else { watts };
        }
    }
    let outward = outward_face_power_w.iter().fold(Dd::ZERO, |sum, &v| sum + Dd::from_f64(v));
    let residual_w = (-outward - Dd::from_f64(absorbed_power_w)).to_f64();
    if !residual_w.is_finite() {
        return Err(bad("crystal box power residual is unrepresentable"));
    }
    poll_illumination(cancelled)?;
    Ok(BoxIntegral {
        absorbed_power_w, outward_face_power_w, residual_w, mode_pair_count, workspace_bytes,
    })
}

fn bad(message: &'static str) -> IlluminationError {
    IlluminationError::Electromagnetic(message)
}

fn finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}

fn coalesce<'a>(
    responses: impl Iterator<Item = &'a AbsorbingCrystalResponse>,
    capacity: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Wave>> {
    let mut waves = Vec::new();
    waves.try_reserve_exact(capacity).map_err(|_| bad("crystal box mode allocation failed"))?;
    for response in responses {
        poll_illumination(cancelled)?;
        for mode in response.transmitted {
            let wave = Wave {
                k: [
                    Dd::from_f64(response.tangential_wavenumber_rad_m)
                        * Dd::from_f64(response.incidence_plane_axis[0]),
                    Dd::from_f64(response.tangential_wavenumber_rad_m)
                        * Dd::from_f64(response.incidence_plane_axis[1]),
                    Dd::from_f64(mode.normal_wavenumber_rad_m.re),
                ],
                alpha: -mode.normal_wavenumber_rad_m.im,
                orientation: 1.0,
                anchor: Dd::ZERO,
                fields: original_fields(mode.electric_v_m, mode.magnetic_a_m),
            };
            merge_wave(&mut waves, wave, cancelled)?;
        }
    }
    retain_active(waves, cancelled)
}

fn coalesce_slab<'a>(
    responses: impl Iterator<Item = &'a AbsorbingSlabResponse>,
    capacity: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Wave>> {
    let mut waves = Vec::new();
    waves.try_reserve_exact(capacity).map_err(|_| bad("film box mode allocation failed"))?;
    for response in responses {
        poll_illumination(cancelled)?;
        let plane = response.incidence_plane_axis();
        for (index, mode) in response.forward.iter().chain(&response.backward).enumerate() {
            let orientation = if index < 2 { 1.0 } else { -1.0 };
            let wave = Wave {
                k: [
                    Dd::from_f64(response.tangential_wavenumber_rad_m) * Dd::from_f64(plane[0]),
                    Dd::from_f64(response.tangential_wavenumber_rad_m) * Dd::from_f64(plane[1]),
                    Dd::from_f64(mode.normal_wavenumber_rad_m.re),
                ],
                alpha: -orientation * mode.normal_wavenumber_rad_m.im,
                orientation,
                anchor: if index < 2 { Dd::ZERO } else { Dd::from_f64(response.thickness_m) },
                fields: original_fields(mode.electric_v_m, mode.magnetic_a_m),
            };
            merge_wave(&mut waves, wave, cancelled)?;
        }
    }
    retain_active(waves, cancelled)
}

fn original_fields(electric: [C64; 3], magnetic: [C64; 3]) -> Fields {
    core::array::from_fn(|j| {
        let value = if j < 3 { electric[j] } else { magnetic[j - 3] };
        [Dd::from_f64(value.re), Dd::from_f64(value.im)]
    })
}

fn merge_wave(waves: &mut Vec<Wave>, wave: Wave, cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if !wave.alpha.is_finite() || wave.alpha < 0.0 || !finite(wave.anchor)
        || wave.k.iter().any(|&value| !finite(value))
    {
        return Err(bad("crystal box requires finite passive anchored eigenwaves"));
    }
    for (index, old) in waves.iter_mut().enumerate() {
        if index % 32 == 0 { poll_illumination(cancelled)?; }
        if old.k == wave.k && old.alpha == wave.alpha
            && old.anchor == wave.anchor && old.orientation == wave.orientation
        {
            for (left, right) in old.fields.iter_mut().zip(wave.fields) {
                *left = [left[0] + right[0], left[1] + right[1]];
            }
            return Ok(());
        }
    }
    waves.push(wave);
    Ok(())
}

fn retain_active(mut waves: Vec<Wave>, cancelled: &mut impl FnMut() -> bool) -> Result<Vec<Wave>> {
    let mut active = 0;
    for index in 0..waves.len() {
        poll_illumination(cancelled)?;
        let wave = waves[index];
        if wave.fields.iter().flatten().any(|&v| !finite(v)) {
            return Err(bad("coalesced crystal box field is unrepresentable"));
        }
        if wave.fields.iter().any(|v| *v != ZERO) {
            waves[active] = wave;
            active += 1;
        }
    }
    waves.truncate(active);
    Ok(waves)
}

fn multiply(a: Complex, b: Complex) -> Complex {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}
fn conjugate_product(a: Complex, b: Complex) -> Complex {
    [a[0] * b[0] + a[1] * b[1], a[0] * b[1] - a[1] * b[0]]
}
fn squared(a: Complex) -> Dd {
    a[0] * a[0] + a[1] * a[1]
}
fn exposure(alpha: f64, depth: Dd) -> f64 {
    let leading = alpha * depth.hi;
    if leading.is_infinite() { leading } else { (Dd::from_f64(alpha) * depth).to_f64() }
}
fn attenuation(value: Complex, x: f64) -> Complex {
    let remaining = decay::remaining_fraction(x);
    if remaining >= f64::MIN_POSITIVE {
        let factor = if x < 0.5 {
            Dd::ONE + Dd::from_f64(det::expm1(-x))
        } else {
            Dd::from_f64(remaining)
        };
        value.map(|v| v * factor)
    } else {
        value.map(|v| {
            Dd::from_f64(decay::attenuate(v.hi, x)) + Dd::from_f64(decay::attenuate(v.lo, x))
        })
    }
}
fn rotation(phase: Dd) -> Result<Complex> {
    if !finite(phase) {
        return Err(bad("crystal box relative phase is unrepresentable"));
    }
    let (sine, cosine) = det::sin_cos_dd(phase);
    Ok([cosine, -sine])
}

fn prepare_samples(
    waves: &[Wave],
    samples: &mut Vec<Sample>,
    observation: Observation,
    loss: &Loss,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 2]> {
    samples.clear();
    let mut reference = None;
    for (index, wave) in waves.iter().enumerate() {
        poll_illumination(cancelled)?;
        // A fixed z face is ONE absolute plane. Its unused extent must not
        // move a backward wave to the opposite, unobserved box endpoint.
        let near = if wave.orientation < 0.0 && observation.integrated[2] {
            observation.origin[2] + Dd::from_f64(observation.extent[2])
        } else {
            observation.origin[2]
        };
        let distance = (near - wave.anchor) * Dd::from_f64(wave.orientation);
        let x = exposure(wave.alpha, distance);
        let fields = wave.fields.map(|v| attenuation(v, x));
        if reference.is_none() && fields.iter().any(|v| *v != ZERO) {
            reference = Some(index);
        }
        samples.push(Sample { fields, projected: [ZERO; 4] });
    }
    let Some(reference) = reference else { return Ok([0.0; 2]); };
    let mut scale = [0.0_f64; 2];
    for (wave, sample) in waves.iter().zip(samples.iter_mut()) {
        poll_illumination(cancelled)?;
        if sample.fields.iter().all(|v| *v == ZERO) { continue; }
        let reference = &waves[reference];
        let transverse = (0..2).fold(Dd::ZERO, |sum, j| {
            sum + (wave.k[j] - reference.k[j]) * observation.origin[j]
        });
        let da = observation.origin[2] - wave.anchor;
        let db = observation.origin[2] - reference.anchor;
        let beat = wave.k[2] - reference.k[2];
        let normal = if wave.anchor == reference.anchor {
            beat * da
        } else if da.to_f64().abs() <= db.to_f64().abs() {
            beat * da + reference.k[2] * (reference.anchor - wave.anchor)
        } else {
            beat * db + wave.k[2] * (reference.anchor - wave.anchor)
        };
        let phase = transverse + normal;
        let phase = rotation(phase)?;
        for (j, value) in sample.fields.iter_mut().enumerate() {
            *value = multiply(*value, phase);
            for part in *value {
                if !finite(part) { return Err(bad("crystal box field is unrepresentable")); }
                scale[j / 3] = scale[j / 3].max(part.hi.abs()).max(part.lo.abs());
            }
        }
    }
    for sample in samples.iter_mut() {
        poll_illumination(cancelled)?;
        for (j, value) in sample.fields.iter_mut().enumerate() {
            if scale[j / 3] != 0.0 {
                let old = *value;
                *value = value.map(|v| v / Dd::from_f64(scale[j / 3]));
                if old != ZERO && *value == ZERO {
                    return Err(bad("crystal box field normalization underflows"));
                }
            }
        }
        let dot = [0, 1].map(|part| {
            (0..3).fold(Dd::ZERO, |sum, j| sum + sample.fields[j][part] * loss.axis[j])
        });
        sample.projected = core::array::from_fn(|j| {
            if j == 3 { dot } else {
                [0, 1].map(|part| sample.fields[j][part] - dot[part] * loss.axis[j] / loss.norm)
            }
        });
    }
    Ok(scale)
}

fn depth_integral(a: f64, b: f64, beat: Dd, width: f64) -> Result<Complex> {
    let length = Dd::from_f64(width);
    if a == b && beat == Dd::ZERO {
        return Ok([absorption::decay_integral(a, length), Dd::ZERO]);
    }
    let rate = [
        -Dd::from_f64(a) * length - Dd::from_f64(b) * length,
        beat * length,
    ];
    if rate.iter().all(|&v| finite(v)) && rate[0].to_f64().abs() + rate[1].to_f64().abs() < 0.5 {
        // Integral_0^1 exp(rate*u) du. Retain the interval even when opacity
        // underflows; 32 coefficients bound truncation below DD roundoff.
        let mut term = ONE;
        let mut mean = ONE;
        for n in 1..=32 {
            term = multiply(term, rate).map(|v| v / Dd::from_f64((n + 1) as f64));
            mean = [mean[0] + term[0], mean[1] + term[1]];
        }
        return Ok(mean.map(|v| v * length));
    }
    absorption::cross_integral([a, b], beat, length).map_err(IlluminationError::from)
}

// Samples retain phase at the observation origin but attenuation at their
// nearest physical endpoint. For a backward-selected pair, u=L-v contributes
// exp(+i beat*L) as well as reversing the beat inside the decaying integral.
fn anchored_depth_integral(a: &Wave, b: &Wave, width: f64) -> Result<Complex> {
    let beat = a.k[2] - b.k[2];
    if a.orientation > 0.0 && b.orientation > 0.0 {
        return depth_integral(a.alpha, b.alpha, beat, width);
    }
    let same = a.orientation == b.orientation;
    let orientation = if same || a.alpha > b.alpha {
        a.orientation
    } else if b.alpha > a.alpha {
        b.orientation
    } else {
        1.0 // Equal opposing rates use the same endpoint in BOTH pair orders.
    };
    let rates = if same { [a.alpha, b.alpha] } else { [(a.alpha - b.alpha).abs(), 0.0] };
    let mut integral = depth_integral(
        rates[0], rates[1], if orientation > 0.0 { beat } else { -beat }, width,
    )?;
    if orientation < 0.0 {
        integral = multiply(integral, rotation(-beat * Dd::from_f64(width))?);
    }
    if !same {
        integral = attenuation(integral, exposure(a.alpha.min(b.alpha), Dd::from_f64(width)));
    }
    Ok(integral)
}

fn depth_scale(waves: &[Wave], observation: Observation) -> f64 {
    if !observation.integrated[2] { return 1.0; }
    let alpha = waves.iter().map(|wave| wave.alpha).fold(f64::INFINITY, f64::min);
    if alpha == 0.0 {
        observation.extent[2]
    } else {
        observation.extent[2].min((1.0 / alpha) * 0.5)
    }
}

fn pair_kernel(
    a: &Wave,
    b: &Wave,
    windows: &[Option<HarmonicWindow>; 2],
    observation: Observation,
    depth: f64,
) -> Result<Complex> {
    let mut weight = ONE;
    for (axis, window) in windows.iter().enumerate() {
        if let Some(window) = window {
            // Reverse the time-owner order: conj(E_i) E_j carries the
            // positive spatial difference phase in the native convention.
            let kernel = window.proper_product_kernel_dd(b.k[axis], a.k[axis])
                .map_err(|_| bad("crystal box transverse integral is unrepresentable"))?;
            let re = kernel.mean_real_product_dd(C64::ONE, C64::ONE)
                .map_err(|_| bad("crystal box transverse kernel is unrepresentable"))?;
            let negative_im = kernel.mean_real_product_dd(C64::new(0.0, 1.0), C64::ONE)
                .map_err(|_| bad("crystal box transverse kernel is unrepresentable"))?;
            weight = multiply(weight, [re * Dd::from_f64(2.0), negative_im * Dd::from_f64(-2.0)]);
        }
    }
    if observation.integrated[2] {
        let axial = anchored_depth_integral(a, b, observation.extent[2])?
            .map(|v| v / Dd::from_f64(depth));
        weight = multiply(weight, axial);
    }
    Ok(weight)
}

fn analytic(
    waves: &[Wave],
    samples: &[Sample],
    observation: Observation,
    quantity: Quantity,
    depth: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Dd; 4]> {
    let window = |axis: usize| -> Result<Option<HarmonicWindow>> {
        if observation.integrated[axis] {
            Ok(Some(HarmonicWindow::new(0.0, observation.extent[axis])
                .map_err(|_| bad("crystal box transverse extent is invalid"))?))
        } else { Ok(None) }
    };
    let windows = [window(0)?, window(1)?];
    let mut sum = [Dd::ZERO; 4];
    let mut magnitude = [Dd::ZERO; 4];
    for (i, a) in waves.iter().enumerate() {
        for (j, b) in waves.iter().enumerate() {
            poll_illumination(cancelled)?;
            let weight = pair_kernel(a, b, &windows, observation, depth)?;
            match quantity {
                Quantity::Heat(weights) => {
                    for channel in 0..4 {
                        if weights[channel] == 0.0 { continue; }
                        let overlap = conjugate_product(samples[i].projected[channel], samples[j].projected[channel]);
                        let term = multiply(overlap, weight)[0];
                        sum[channel] = sum[channel] + term;
                        magnitude[channel] = magnitude[channel] + Dd::from_f64(term.to_f64().abs());
                    }
                }
                Quantity::Flux(axis) => {
                    let p = (axis + 1) % 3;
                    let q = (axis + 2) % 3;
                    let first = conjugate_product(samples[i].fields[p], samples[j].fields[q + 3]);
                    let second = conjugate_product(samples[i].fields[q], samples[j].fields[p + 3]);
                    let overlap = [first[0] - second[0], first[1] - second[1]];
                    sum[0] = sum[0] + multiply(overlap, weight)[0];
                }
            }
        }
    }
    if let Quantity::Heat(weights) = quantity {
        for (channel, (value, scale)) in sum.iter().zip(magnitude).enumerate() {
            if weights[channel] == 0.0 { continue; }
            if value.to_f64() < 0.0 || (*value != Dd::ZERO && value.to_f64() <= 64.0 * f64::EPSILON * scale.to_f64())
                || (*value == Dd::ZERO && scale != Dd::ZERO)
            {
                return Err(bad("crystal box loss cancellation is unresolved; subdivide the box"));
            }
        }
    }
    Ok(sum)
}

// Positive 8-point Gauss-Legendre rule on [0,1]. We normalize the retained
// weights in DD so a constant integrand has exactly unit quadrature measure.
const NODES: [f64; 8] = [
    0.019_855_071_751_231_884, 0.101_666_761_293_186_63,
    0.237_233_795_041_835_5, 0.408_282_678_752_175_1,
    0.591_717_321_247_824_9, 0.762_766_204_958_164_5,
    0.898_333_238_706_813_4, 0.980_144_928_248_768_1,
];
const WEIGHTS: [f64; 8] = [
    0.050_614_268_145_188_13, 0.111_190_517_226_687_24,
    0.156_853_322_938_943_65, 0.181_341_891_689_180_98,
    0.181_341_891_689_180_98, 0.156_853_322_938_943_65,
    0.111_190_517_226_687_24, 0.050_614_268_145_188_13,
];

struct ShortBox {
    varying: [bool; 3],
    common_normal: bool,
    derivative_scale: [f64; 3],
}

fn short_box(waves: &[Wave], observation: Observation) -> Option<ShortBox> {
    let first = waves.first()?;
    let mut phase = [0.0_f64; 3];
    let mut alpha = 0.0_f64;
    let common_normal = waves.iter().all(|wave| {
        wave.k[2] == first.k[2] && wave.alpha == first.alpha && wave.orientation == first.orientation
    });
    for wave in waves {
        for axis in 0..3 {
            if observation.integrated[axis] {
                let width_phase = ((wave.k[axis] - first.k[axis]) * Dd::from_f64(observation.extent[axis])).to_f64().abs();
                if !width_phase.is_finite() { return None; }
                phase[axis] = phase[axis].max(width_phase);
            }
        }
        if observation.integrated[2] && !common_normal {
            alpha = alpha.max(wave.alpha * observation.extent[2]);
        }
    }
    if phase.iter().any(|&v| v > 0.125) || !alpha.is_finite() || alpha > 0.125 {
        return None;
    }
    Some(ShortBox {
        varying: [
            observation.integrated[0] && phase[0] != 0.0,
            observation.integrated[1] && phase[1] != 0.0,
            observation.integrated[2] && !common_normal,
        ],
        common_normal,
        derivative_scale: [2.0 * phase[0], 2.0 * phase[1], 2.0 * (phase[2] + alpha)],
    })
}

fn short_integral(
    waves: &[Wave],
    samples: &[Sample],
    observation: Observation,
    quantity: Quantity,
    depth: f64,
    short: ShortBox,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Dd; 4]> {
    let weight_sum = WEIGHTS.iter().fold(Dd::ZERO, |sum, &v| sum + Dd::from_f64(v));
    let weight = |axis: usize, index: usize| if short.varying[axis] {
        Dd::from_f64(WEIGHTS[index]) / weight_sum
    } else { Dd::ONE };
    let coordinate = |axis: usize, index: usize| if short.varying[axis] {
        Dd::from_f64(observation.extent[axis]) * Dd::from_f64(NODES[index])
    } else { Dd::ZERO };
    let count = short.varying.map(|v| if v { 8 } else { 1 });
    let mut sum = [Dd::ZERO; 4];
    for x in 0..count[0] {
        for y in 0..count[1] {
            for z in 0..count[2] {
                poll_illumination(cancelled)?;
                let point = [coordinate(0, x), coordinate(1, y), coordinate(2, z)];
                let w = weight(0, x) * weight(1, y) * weight(2, z);
                let mut fields = [ZERO; 6];
                let mut projected = [ZERO; 4];
                for (index, (wave, sample)) in waves.iter().zip(samples).enumerate() {
                    if index % 32 == 0 { poll_illumination(cancelled)?; }
                    let phase = (0..3).fold(Dd::ZERO, |value, axis| {
                        value + (wave.k[axis] - waves[0].k[axis]) * point[axis]
                    });
                    let phase = rotation(phase)?;
                    let distance = if wave.orientation < 0.0 && observation.integrated[2] {
                        Dd::from_f64(observation.extent[2]) - point[2]
                    } else {
                        point[2]
                    };
                    let x = if short.common_normal { 0.0 } else { exposure(wave.alpha, distance) };
                    match quantity {
                        Quantity::Heat(_) => {
                            for (sum, value) in projected.iter_mut().zip(sample.projected) {
                                let v = multiply(attenuation(value, x), phase);
                                *sum = [sum[0] + v[0], sum[1] + v[1]];
                            }
                        }
                        Quantity::Flux(_) => {
                            for (sum, value) in fields.iter_mut().zip(sample.fields) {
                                let v = multiply(attenuation(value, x), phase);
                                *sum = [sum[0] + v[0], sum[1] + v[1]];
                            }
                        }
                    }
                }
                match quantity {
                    Quantity::Heat(_) => {
                        for (sum, value) in sum.iter_mut().zip(projected) {
                            *sum = *sum + w * squared(value);
                        }
                    }
                    Quantity::Flux(axis) => {
                        let p = (axis + 1) % 3;
                        let q = (axis + 2) % 3;
                        sum[0] = sum[0] + w * (conjugate_product(fields[p], fields[q + 3])[0]
                            - conjugate_product(fields[q], fields[p + 3])[0]);
                    }
                }
            }
        }
    }
    if let Quantity::Heat(weights) = quantity {
        // For the monic shifted Legendre polynomial, the Gauss remainder
        // coefficient is 8!^4 / (17 * 16!^3) < 1.7e-23. Each exponential
        // pair has |16th derivative| <= amplitude_bound * scale^16.
        // Tensor errors add because both exact and positive quadrature
        // measures have norm one. This screens analytic truncation ONLY.
        let remainder = 1.7e-23 * short.derivative_scale.iter().map(|v| {
            let square = v * v;
            let fourth = square * square;
            let eighth = fourth * fourth;
            eighth * eighth
        }).sum::<f64>();
        for channel in 0..4 {
            if weights[channel] == 0.0 { continue; }
            let amplitude = samples.iter().fold(Dd::ZERO, |a, sample| {
                a + Dd::from_f64(sample.projected[channel][0].to_f64().abs())
                    + Dd::from_f64(sample.projected[channel][1].to_f64().abs())
            });
            let bound = amplitude * amplitude * Dd::from_f64(remainder);
            if sum[channel] == Dd::ZERO && amplitude != Dd::ZERO
                || bound.to_f64() > 64.0 * f64::EPSILON * sum[channel].to_f64()
            {
                return Err(bad("short crystal box loss integral is unresolved"));
            }
        }
    }
    let axial = if observation.integrated[2] {
        if short.common_normal {
            absorption::decay_integral(waves[0].alpha, Dd::from_f64(observation.extent[2]))
                / Dd::from_f64(depth)
        } else {
            Dd::from_f64(observation.extent[2]) / Dd::from_f64(depth)
        }
    } else { Dd::ONE };
    Ok(sum.map(|value| value * axial))
}

fn restore(value: Dd, factors: &[f64]) -> Result<f64> {
    if !finite(value) { return Err(bad("crystal box integral is unrepresentable")); }
    if value == Dd::ZERO { return Ok(0.0); }
    let mut product = value;
    for &factor in factors { product = product * Dd::from_f64(factor); }
    let direct = product.to_f64();
    if direct.is_finite() && direct != 0.0 { return Ok(direct); }
    let original = value.to_f64();
    if !original.is_finite() || original == 0.0 {
        return Err(bad("crystal box normalized power is unrepresentable"));
    }
    let logarithm = factors.iter().fold(det::ln(original.abs()), |sum, &factor| sum + det::ln(factor));
    let magnitude = det::exp(logarithm);
    if !magnitude.is_finite() || magnitude == 0.0 {
        return Err(bad("crystal box physical power is unrepresentable"));
    }
    Ok(if original < 0.0 { -magnitude } else { magnitude })
}

fn measure(
    waves: &[Wave],
    samples: &mut Vec<Sample>,
    observation: Observation,
    quantity: Quantity,
    loss: &Loss,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    poll_illumination(cancelled)?;
    if waves.is_empty() { return Ok(0.0); }
    let [electric, magnetic] = prepare_samples(waves, samples, observation, loss, cancelled)?;
    if electric == 0.0 || matches!(quantity, Quantity::Flux(_)) && magnetic == 0.0 {
        return Ok(0.0);
    }
    let depth = depth_scale(waves, observation);
    if !depth.is_finite() || depth <= 0.0 {
        return Err(bad("crystal box integration depth scale is unrepresentable"));
    }
    let integral = if let Some(short) = short_box(waves, observation) {
        short_integral(waves, samples, observation, quantity, depth, short, cancelled)?
    } else {
        analytic(waves, samples, observation, quantity, depth, cancelled)?
    };
    let transverse = (0..2).filter(|&axis| observation.integrated[axis])
        .fold(Dd::ONE, |area, axis| area * Dd::from_f64(observation.extent[axis]))
        .to_f64();
    let power = match quantity {
        Quantity::Heat(_) => {
            let mut sum = Dd::ZERO;
            for (value, &weight) in integral.iter().zip(&loss.weight) {
                if weight == 0.0 { continue; }
                if value.to_f64() < 0.0 { return Err(bad("crystal box dielectric work is negative")); }
                sum = sum + Dd::from_f64(restore(*value, &[depth, transverse, weight, 0.5, loss.omega, electric, electric])?);
            }
            sum.to_f64()
        }
        Quantity::Flux(_) => restore(integral[0], &[depth, transverse, 0.5, electric, magnetic])?,
    };
    if !power.is_finite() { return Err(bad("crystal box power is unrepresentable")); }
    poll_illumination(cancelled)?;
    Ok(power)
}

#[cfg(test)]
#[path = "box_integral/tests.rs"]
mod tests;
