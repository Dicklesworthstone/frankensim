//! Continuous rectangular support for the original irreversible reservoirs.
//!
//! Centered spatial Gram weights supply the normalized difference/sum Fourier
//! integrals. Original grouped fields receive one joint x/y center rotation
//! before products, retaining their source time origin. Exponent-scaled DD
//! contractions preserve the original heat owner's representable dynamic range.

use super::*;

pub(super) type Electric = [[ScaledDd; 2]; 3];

pub(super) struct Plane {
    modes: Vec<Electric>,
    center: Electric,
}

#[derive(Clone, Copy)]
struct Weight {
    difference_decrement: Dd,
    sum_decrement: Dd,
}

struct Axis {
    wavenumbers: Vec<f64>,
    weights: Vec<Weight>,
}

pub(super) struct PreparedRectangle {
    grid: AngularSpectrumGrid,
    aperture: AngularSpectrumAperture,
    x: Axis,
    y: Axis,
}

impl PreparedRectangle {
    /// Logical retained cache storage; the caller separately admits fields,
    /// depth/frequency/modal work, original reconstruction and spectral output.
    pub(super) fn workspace_bytes(
        grid: AngularSpectrumGrid,
    ) -> Result<usize, AngularSpectrumError> {
        if grid.width == 0 || grid.height == 0 {
            return Err(bad("rectangular heat requires a nonempty original grid"));
        }
        let axes = grid
            .width
            .checked_add(grid.height)
            .and_then(|v| v.checked_mul(core::mem::size_of::<f64>()))
            .ok_or_else(|| bad("rectangular heat axis storage overflows"))?;
        grid.width
            .checked_mul(grid.width)
            .and_then(|x| {
                grid.height
                    .checked_mul(grid.height)
                    .and_then(|y| x.checked_add(y))
            })
            .and_then(|v| v.checked_mul(core::mem::size_of::<Weight>()))
            .and_then(|v| v.checked_add(axes))
            .and_then(|v| v.checked_add(core::mem::size_of::<Self>()))
            .ok_or_else(|| bad("rectangular heat kernel storage overflows"))
    }

    pub(super) fn prepare(
        grid: AngularSpectrumGrid,
        aperture: AngularSpectrumAperture,
        kx: impl Fn(usize) -> f64,
        ky: impl Fn(usize) -> f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        poll(cancelled)?;
        Self::workspace_bytes(grid)?;
        let extent = aperture.extent_m();
        Ok(Self {
            grid,
            aperture,
            x: Axis::prepare(grid.width, extent[0], kx, cancelled)?,
            y: Axis::prepare(grid.height, extent[1], ky, cancelled)?,
        })
    }

    // Keep start and half-width phases separate, as in the original detector
    // owner, so a tiny patch is not rounded into a distant spatial coordinate.
    // Sum x/y phases before trig: translating a spatial dark fringe must not
    // leave a false residual from independently rounded axis rotations.
    fn electric(
        &self,
        field: Vec<[[Dd; 2]; 6]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Plane, AngularSpectrumError> {
        if field.len() != self.grid.width * self.grid.height {
            return Err(bad(
                "rectangular heat field must cover the original modal grid",
            ));
        }
        let mut result = Vec::new();
        result
            .try_reserve_exact(field.len())
            .map_err(|_| bad("rectangular heat field allocation failed"))?;
        let origin = self.aperture.origin_m();
        let extent = self.aperture.extent_m();
        let mut center = [[ScaledDd::ZERO; 2]; 3];
        for (i, field) in field.into_iter().enumerate() {
            if i % 32 == 0 {
                poll(cancelled)?;
            }
            let k = [
                self.x.wavenumbers[i % self.grid.width],
                self.y.wavenumbers[i / self.grid.width],
            ];
            let (mut start, mut half_width) = (Dd::ZERO, Dd::ZERO);
            for axis in 0..2 {
                let number = Dd::from_f64(k[axis]);
                start = start + number * Dd::from_f64(origin[axis]);
                half_width = half_width
                    + if k[axis].abs() >= extent[axis] {
                        (number * Dd::from_f64(0.5)) * Dd::from_f64(extent[axis])
                    } else {
                        number * Dd::from_f64(0.5 * extent[axis])
                    };
            }
            if !finite_dd(start) || !finite_dd(half_width) {
                return Err(bad("rectangular heat center phase is unrepresentable"));
            }
            let (s, c) = det::sin_cos_dd(start);
            let (t, d) = det::sin_cos_dd(half_width);
            let rotation = [
                ScaledDd::from_dd(c * d - s * t),
                ScaledDd::from_dd(-(s * d + c * t)),
            ];
            let electric: Electric =
                core::array::from_fn(|axis| multiply(field[axis].map(ScaledDd::from_dd), rotation));
            if electric.iter().any(|&v| !finite_parts(v)) {
                return Err(bad(
                    "rectangular heat centered electric field is unrepresentable",
                ));
            }
            for axis in 0..3 {
                add_parts(&mut center[axis], electric[axis]);
            }
            result.push(electric);
        }
        Ok(Plane {
            modes: result,
            center,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn contract_bank<F: FnMut() -> bool>(
        &self,
        frequencies: &[f64],
        proper: bool,
        quadrature: AngularSpectrumHeatVolumeQuadrature<'_>,
        volume: VolumeAdmission,
        reservoir: impl FnMut(usize, usize, &mut F) -> Result<PairIntegral, AngularSpectrumError>,
        mut field: impl FnMut(usize, f64, &mut F) -> Result<Vec<[[Dd; 2]; 6]>, AngularSpectrumError>,
        terms: &mut Vec<AngularHarmonicSpectrumTerm>,
        cancelled: &mut F,
    ) -> Result<(), AngularSpectrumError> {
        poll(cancelled)?;
        if frequencies.is_empty() {
            return Ok(());
        }
        let mut pairs = prepare_pairs(frequencies, reservoir, cancelled)?;
        let mut fields = Vec::new();
        fields
            .try_reserve_exact(frequencies.len())
            .map_err(|_| bad("rectangular heat frequency plane allocation failed"))?;
        for sample in quadrature.depth_samples {
            poll(cancelled)?;
            fields.clear();
            for group in 0..frequencies.len() {
                let original = field(group, sample.depth_m, cancelled)?;
                fields.push(self.electric(original, cancelled)?);
            }
            let area_depth =
                ScaledDd::from_dd(volume.area_m2) * ScaledDd::from_f64(sample.weight_m);
            // Exact identity: W = 1 + (W-1). Sum the actual center fields
            // BEFORE products, then add spatial-variation corrections. A dark
            // patch or weak grouped residual must not be recovered by
            // subtracting large modal powers after its low part was lost.
            for (p, pair) in pairs.iter_mut().enumerate() {
                if p % 32 == 0 {
                    poll(cancelled)?;
                }
                let dot = electric_product(
                    &fields[pair.left].center,
                    &fields[pair.right].center,
                    proper,
                );
                accumulate(&mut pair.spatial, dot, area_depth, area_depth, proper);
            }
            for i in 0..volume.points {
                poll(cancelled)?;
                for j in 0..volume.points {
                    if j % 32 == 0 {
                        poll(cancelled)?;
                    }
                    let x = self.x.weights
                        [(i % self.grid.width) * self.grid.width + j % self.grid.width];
                    let y = self.y.weights
                        [(i / self.grid.width) * self.grid.height + j / self.grid.width];
                    // (1+dx)(1+dy)-1 = dx+dy+dx*dy. Form these decrements
                    // without subtracting two nearly equal transforms near one.
                    let difference =
                        area_depth * decrement(x.difference_decrement, y.difference_decrement);
                    let sum = if proper {
                        ScaledDd::ZERO
                    } else {
                        area_depth * decrement(x.sum_decrement, y.sum_decrement)
                    };
                    if difference.is_zero() && (proper || sum.is_zero()) {
                        continue;
                    }
                    for (p, pair) in pairs.iter_mut().enumerate() {
                        if p % 32 == 0 {
                            poll(cancelled)?;
                        }
                        let dot = electric_product(
                            &fields[pair.left].modes[i],
                            &fields[pair.right].modes[j],
                            proper,
                        );
                        accumulate(&mut pair.spatial, dot, difference, sum, proper);
                    }
                }
            }
        }
        emit_pairs(frequencies, proper, pairs, terms, cancelled)
    }
}

impl Axis {
    fn prepare(
        count: usize,
        extent: f64,
        number: impl Fn(usize) -> f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        let pairs = count
            .checked_mul(count)
            .ok_or_else(|| bad("rectangular heat axis work overflows"))?;
        let window = HarmonicWindow::new(0.0, extent)
            .map_err(|_| bad("rectangular heat extent is unrepresentable"))?;
        let mut wavenumbers = Vec::new();
        wavenumbers
            .try_reserve_exact(count)
            .map_err(|_| bad("rectangular heat axis allocation failed"))?;
        for i in 0..count {
            if i % 32 == 0 {
                poll(cancelled)?;
            }
            let value = number(i);
            if !value.is_finite() {
                return Err(bad("rectangular heat wavenumber is unrepresentable"));
            }
            wavenumbers.push(value);
        }
        let mut weights = Vec::new();
        weights
            .try_reserve_exact(pairs)
            .map_err(|_| bad("rectangular heat kernel allocation failed"))?;
        for &a in &wavenumbers {
            for (j, &b) in wavenumbers.iter().enumerate() {
                if j % 32 == 0 {
                    poll(cancelled)?;
                }
                let [c, s] = window
                    .product_kernel(a, b)
                    .map_err(|_| bad("rectangular heat spatial phase is unrepresentable"))?
                    .centered_gram_dd();
                weights.push(Weight {
                    difference_decrement: (c - Dd::ONE) + s,
                    sum_decrement: (c - Dd::ONE) - s,
                });
            }
        }
        Ok(Self {
            wavenumbers,
            weights,
        })
    }
}

fn decrement(x: Dd, y: Dd) -> ScaledDd {
    let (x, y) = (ScaledDd::from_dd(x), ScaledDd::from_dd(y));
    x + y + x * y
}

fn electric_product(a: &Electric, b: &Electric, proper: bool) -> PairIntegral {
    let mut dot = PairIntegral::default();
    for axis in 0..3 {
        add_parts(&mut dot.difference, multiply_conjugate(a[axis], b[axis]));
        if !proper {
            add_parts(&mut dot.sum, multiply(a[axis], b[axis]));
        }
    }
    dot
}

fn accumulate(
    target: &mut PairIntegral,
    value: PairIntegral,
    difference: ScaledDd,
    sum: ScaledDd,
    proper: bool,
) {
    for part in 0..2 {
        target.difference[part] = target.difference[part] + value.difference[part] * difference;
        if !proper {
            target.sum[part] = target.sum[part] + value.sum[part] * sum;
        }
    }
}

#[cfg(test)]
mod tests;
