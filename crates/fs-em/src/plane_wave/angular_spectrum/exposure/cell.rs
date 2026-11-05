//! Continuous real-field products for a common nonzero transverse carrier.
//!
//! Center each mode in the physical cell before forming products. Reflection
//! symmetry leaves cosine/cosine and sine/sine Gram matrices on each axis.
//! Their tensor product is applied in two one-dimensional contractions, never
//! by an N-by-N spatial pair loop. Native harmonic kernels retain the common
//! carrier and each envelope bin separately through their DD sums/differences.

use super::*;
use fs_fft::harmonic_window::WindowError;

type Parts = [Dd; 2];
type PairPlane = [[Parts; 2]; 6];

#[derive(Clone, Copy)]
pub(super) struct Coefficient {
    source: Parts,
    cosine: Parts,
    sine: Parts,
}
impl Coefficient {
    const ZERO: Self = Self {
        source: [Dd::ZERO; 2],
        cosine: [Dd::ZERO; 2],
        sine: [Dd::ZERO; 2],
    };
}

pub(super) struct CarrierCell {
    grid: AngularSpectrumGrid,
    x: Vec<[Dd; 2]>,
    y: Vec<[Dd; 2]>,
    phase: Vec<Parts>,
}

pub(super) struct CarrierFields {
    fields: Vec<Vec<[Coefficient; 6]>>,
}

impl CarrierFields {
    pub(super) fn at(&self, frequency: usize, point: usize) -> [Coefficient; 6] {
        self.fields[frequency][point]
    }
}

/// Additional axis contractions and peak logical numerical bytes. The
/// original bank separately admits its C64 modal planes and temporal kernels.
/// One dual-Gram scratch plane is reused across all frequency groups.
pub(super) fn admission(
    grid: AngularSpectrumGrid,
    frequencies: usize,
    carrier: [f64; 2],
) -> Result<(usize, usize), AngularSpectrumError> {
    if carrier == [0.0; 2] || frequencies == 0 {
        return Ok((0, 0));
    }
    let n = grid
        .width
        .checked_mul(grid.height)
        .ok_or_else(|| bad("carrier cell point count overflow"))?;
    let axes = grid
        .width
        .checked_add(grid.height)
        .ok_or_else(|| bad("carrier cell axis work overflow"))?;
    let work = n
        .checked_mul(frequencies)
        .and_then(|v| v.checked_mul(axes))
        .ok_or_else(|| bad("carrier cell contraction work overflow"))?;
    let axis_pairs = grid
        .width
        .checked_mul(grid.width)
        .and_then(|x| {
            grid.height
                .checked_mul(grid.height)
                .and_then(|y| x.checked_add(y))
        })
        .ok_or_else(|| bad("carrier cell kernel count overflow"))?;
    let bytes = axis_pairs
        .checked_mul(core::mem::size_of::<[Dd; 2]>())
        .and_then(|v| v.checked_add(n.checked_mul(core::mem::size_of::<Parts>())?))
        .and_then(|v| v.checked_add(n.checked_mul(core::mem::size_of::<PairPlane>())?))
        .and_then(|v| {
            v.checked_add(
                n.checked_mul(frequencies)?
                    .checked_mul(core::mem::size_of::<[Coefficient; 6]>())?,
            )
        })
        .and_then(|v| {
            v.checked_add(frequencies.checked_mul(core::mem::size_of::<Vec<[Coefficient; 6]>>())?)
        })
        .and_then(|v| {
            v.checked_add(
                core::mem::size_of::<CarrierCell>()
                    + core::mem::size_of::<CarrierFields>()
                    + core::mem::size_of::<Vec<PairPlane>>(),
            )
        })
        .ok_or_else(|| bad("carrier cell workspace overflow"))?;
    Ok((work, bytes))
}

impl CarrierCell {
    pub(super) fn prepare(
        beam: &PlanarAngularSpectrum<'_, '_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        poll(cancelled)?;
        let grid = beam.grid();
        let extent = [
            grid.width as f64 * grid.spacing_m[0],
            grid.height as f64 * grid.spacing_m[1],
        ];
        let x = axis(
            grid.width,
            extent[0],
            beam.carrier_rad_m[0],
            |i| beam.envelope_wavenumber(0, i),
            cancelled,
        )?;
        let y = axis(
            grid.height,
            extent[1],
            beam.carrier_rad_m[1],
            |i| beam.envelope_wavenumber(1, i),
            cancelled,
        )?;
        let n = beam.modes.len();
        let mut phase = Vec::new();
        phase
            .try_reserve_exact(n)
            .map_err(|_| bad("carrier cell phase allocation failed"))?;
        for point in 0..n {
            if point % 32 == 0 {
                poll(cancelled)?;
            }
            let kx = Dd::from_f64(beam.carrier_rad_m[0])
                + Dd::from_f64(beam.envelope_wavenumber(0, point % grid.width));
            let ky = Dd::from_f64(beam.carrier_rad_m[1])
                + Dd::from_f64(beam.envelope_wavenumber(1, point / grid.width));
            let half = Dd::from_f64(0.5);
            let angle = half_product(kx, extent[0], half) + half_product(ky, extent[1], half);
            if !finite_dd(angle) {
                return Err(bad("carrier cell center phase is unrepresentable"));
            }
            let (sine, cosine) = det::sin_cos_dd(angle);
            if !finite_dd(sine) || !finite_dd(cosine) {
                return Err(bad("carrier cell center rotation is unrepresentable"));
            }
            phase.push([cosine, sine]);
        }
        poll(cancelled)?;
        Ok(Self { grid, x, y, phase })
    }

    pub(super) fn fields(
        &self,
        fields: &[Vec<[C64; 6]>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CarrierFields, AngularSpectrumError> {
        let n = self.phase.len();
        let (width, height) = (self.grid.width, self.grid.height);
        let mut result = Vec::new();
        result
            .try_reserve_exact(fields.len())
            .map_err(|_| bad("carrier cell frequency allocation failed"))?;
        let mut scratch = Vec::new();
        scratch
            .try_reserve_exact(n)
            .map_err(|_| bad("carrier cell contraction allocation failed"))?;
        for point in 0..n {
            if point % 64 == 0 {
                poll(cancelled)?;
            }
            scratch.push([[[Dd::ZERO; 2]; 2]; 6]);
        }
        for field in fields {
            poll(cancelled)?;
            if field.len() != n {
                return Err(bad("carrier cell field differs from the original grid"));
            }
            let mut plane = Vec::new();
            plane
                .try_reserve_exact(n)
                .map_err(|_| bad("carrier cell field allocation failed"))?;
            for (point, value) in field.iter().enumerate() {
                if point % 32 == 0 {
                    poll(cancelled)?;
                }
                let [cosine, sine] = self.phase[point];
                // Conjugate the EM amplitude into the temporal owner's
                // negative-time convention, including its spatial center phase.
                plane.push(value.map(|z| Coefficient {
                    source: [
                        Dd::from_f64(z.re) * cosine + Dd::from_f64(z.im) * sine,
                        Dd::from_f64(z.re) * sine - Dd::from_f64(z.im) * cosine,
                    ],
                    ..Coefficient::ZERO
                }));
            }
            for row in 0..height {
                for column in 0..width {
                    poll(cancelled)?;
                    let target = row * width + column;
                    scratch[target] = [[[Dd::ZERO; 2]; 2]; 6];
                    for other in 0..width {
                        if other % 32 == 0 {
                            poll(cancelled)?;
                        }
                        let weights = self.x[column * width + other];
                        let source = &plane[row * width + other];
                        for component in 0..6 {
                            for quadrature in 0..2 {
                                for part in 0..2 {
                                    scratch[target][component][quadrature][part] = scratch[target]
                                        [component][quadrature][part]
                                        + weights[quadrature] * source[component].source[part];
                                }
                            }
                        }
                    }
                }
            }
            for row in 0..height {
                for column in 0..width {
                    poll(cancelled)?;
                    let target = &mut plane[row * width + column];
                    for other in 0..height {
                        if other % 32 == 0 {
                            poll(cancelled)?;
                        }
                        let [c, s] = self.y[row * height + other];
                        let source = &scratch[other * width + column];
                        for component in 0..6 {
                            for part in 0..2 {
                                target[component].cosine[part] = target[component].cosine[part]
                                    + c * source[component][0][part]
                                    + s * source[component][1][part];
                                target[component].sine[part] = target[component].sine[part]
                                    + c * source[component][1][part]
                                    + s * source[component][0][part];
                            }
                        }
                    }
                    if target.iter().any(|value| {
                        value
                            .source
                            .iter()
                            .chain(&value.cosine)
                            .chain(&value.sine)
                            .any(|&part| !finite_dd(part))
                    }) {
                        return Err(bad("carrier cell contraction is unrepresentable"));
                    }
                }
            }
            result.push(plane);
        }
        poll(cancelled)?;
        Ok(CarrierFields { fields: result })
    }
}

fn axis(
    count: usize,
    extent: f64,
    carrier: f64,
    envelope: impl Fn(usize) -> f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<[Dd; 2]>, AngularSpectrumError> {
    let window = HarmonicWindow::new(0.0, extent)
        .map_err(|_| bad("carrier cell axis extent is unrepresentable"))?;
    let mut result = Vec::new();
    result
        .try_reserve_exact(count * count)
        .map_err(|_| bad("carrier cell kernel allocation failed"))?;
    for i in 0..count {
        for j in 0..count {
            if j % 32 == 0 {
                poll(cancelled)?;
            }
            result.push(
                window
                    .product_kernel_shifted(envelope(i), envelope(j), carrier)
                    .map_err(|_| bad("carrier cell spatial kernel is unrepresentable"))?
                    .centered_gram_dd(),
            );
        }
    }
    Ok(result)
}

fn half_product(number: Dd, extent: f64, half: Dd) -> Dd {
    if number.hi.abs() >= extent {
        (number * half) * Dd::from_f64(extent)
    } else {
        number * Dd::from_f64(0.5 * extent)
    }
}

pub(super) fn product<P: Copy>(
    kernel: P,
    left: Coefficient,
    right: Coefficient,
    temporal: impl Fn(P, Parts, Parts) -> Result<Dd, WindowError>,
) -> Result<Dd, WindowError> {
    let rotate = |[real, imaginary]: Parts| [-imaginary, real];
    let result = temporal(kernel, left.source, right.cosine)?
        + temporal(kernel, rotate(left.source), rotate(right.sine))?;
    if !finite_dd(result) {
        return Err(WindowError::Unrepresentable);
    }
    Ok(result)
}

// Temporal differentiation of the already conjugated amplitudes rotates -i.
// Rotate all retained representations so it commutes with spatial integration.
pub(super) fn rotate_minus_i(value: Coefficient) -> Coefficient {
    let rotate = |[real, imaginary]: Parts| [imaginary, -real];
    Coefficient {
        source: rotate(value.source),
        cosine: rotate(value.cosine),
        sine: rotate(value.sine),
    }
}

fn finite_dd(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}

pub(super) fn power<P: Copy>(
    grid: AngularSpectrumGrid,
    frequencies: usize,
    fields: &CarrierFields,
    kernels: &[P],
    temporal: impl Fn(P, Parts, Parts) -> Result<Dd, WindowError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Dd, AngularSpectrumError> {
    let mut sum = Dd::ZERO;
    for point in 0..grid.width * grid.height {
        poll(cancelled)?;
        for i in 0..frequencies {
            let a = fields.at(i, point);
            for j in 0..frequencies {
                if j % 32 == 0 {
                    poll(cancelled)?;
                }
                let b = fields.at(j, point);
                let kernel = kernels[i * frequencies + j];
                sum = sum
                    + product(kernel, a[0], b[4], &temporal)
                        .map_err(|_| bad("carrier cell power product is unrepresentable"))?
                    - product(kernel, a[1], b[3], &temporal)
                        .map_err(|_| bad("carrier cell power product is unrepresentable"))?;
            }
        }
    }
    let area = Dd::from_f64(grid.width as f64 * grid.spacing_m[0])
        * Dd::from_f64(grid.height as f64 * grid.spacing_m[1]);
    let power = sum * area;
    if !finite_dd(power) {
        return Err(bad("integrated carrier cell power is unrepresentable"));
    }
    poll(cancelled)?;
    Ok(power)
}
