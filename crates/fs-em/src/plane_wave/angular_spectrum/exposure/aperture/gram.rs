//! Shared separable aperture/shutter Gram application on original E/H fields.
use super::*;

pub(super) struct GramWorkspace {
    pub(super) residual: Vec<FieldParts>,
    pub(super) covectors: Vec<FieldParts>,
    odd_x: Vec<FieldParts>,
    odd_y: Vec<FieldParts>,
    modes: usize,
}
impl GramWorkspace {
    pub(super) fn new(
        modes: usize,
        frequencies: usize,
        real_frequencies: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        let n = modes
            .checked_mul(frequencies)
            .ok_or_else(|| bad("aperture Gram extent overflow"))?;
        let real = modes
            .checked_mul(real_frequencies)
            .ok_or_else(|| bad("aperture Gram extent overflow"))?;
        Ok(Self {
            residual: zeros(n, cancelled)?,
            covectors: zeros(n, cancelled)?,
            odd_x: zeros(n, cancelled)?,
            odd_y: zeros(real, cancelled)?,
            modes,
        })
    }

    // Original real even/odd operation graph: preserve small positive spatial
    // and shutter variances without subtracting rounded near-unit kernels.
    pub(super) fn real(
        &mut self,
        spatial: &Spatial<RealKernel>,
        temporal: &[RealKernel],
        f: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AngularSpectrumError> {
        let n = self.modes * f;
        clear(&mut self.covectors[..n], cancelled)?;
        clear(&mut self.odd_x[..n], cancelled)?;
        clear(&mut self.odd_y[..n], cancelled)?;
        let (nx, ny) = (spatial.grid.width, spatial.grid.height);
        // x: seed both even and odd spatial parity for each real and
        // imaginary field component. No difference of near-unit kernels.
        for frequency in 0..f {
            for y in 0..ny {
                for x in 0..nx {
                    poll(cancelled)?;
                    let target = frequency * self.modes + y * nx + x;
                    for q in 0..nx {
                        if q % 32 == 0 {
                            poll(cancelled)?;
                        }
                        let [c, s] = spatial.x[x * nx + q].0;
                        let input = &self.residual[frequency * self.modes + y * nx + q];
                        for component in 0..6 {
                            for part in 0..2 {
                                self.covectors[target][component][part] = self.covectors[target]
                                    [component][part]
                                    + c * input[component][part];
                                self.odd_x[target][component][part] = self.odd_x[target][component]
                                    [part]
                                    + s * input[component][part];
                            }
                        }
                    }
                }
            }
        }
        // y: compose parity. Initial residual storage is now reusable;
        // the original source fields remain immutable for the pullback.
        clear(&mut self.residual[..n], cancelled)?;
        for frequency in 0..f {
            for y in 0..ny {
                for x in 0..nx {
                    poll(cancelled)?;
                    let target = frequency * self.modes + y * nx + x;
                    for q in 0..ny {
                        if q % 32 == 0 {
                            poll(cancelled)?;
                        }
                        let [c, s] = spatial.y[y * ny + q].0;
                        let input = frequency * self.modes + q * nx + x;
                        for component in 0..6 {
                            for part in 0..2 {
                                let even = self.covectors[input][component][part];
                                let odd = self.odd_x[input][component][part];
                                self.residual[target][component][part] =
                                    self.residual[target][component][part] + c * even + s * odd;
                                self.odd_y[target][component][part] =
                                    self.odd_y[target][component][part] + s * even + c * odd;
                            }
                        }
                    }
                }
            }
        }
        // time: real fields take total even parity; imaginary fields
        // take total odd parity. The original kernel owns all cycle halves.
        clear(&mut self.covectors[..n], cancelled)?;
        for a in 0..f {
            for mode in 0..self.modes {
                poll(cancelled)?;
                let target = a * self.modes + mode;
                for b in 0..f {
                    if b % 32 == 0 {
                        poll(cancelled)?;
                    }
                    let [c, s] = temporal[a * f + b].0;
                    let input = b * self.modes + mode;
                    for component in 0..6 {
                        self.covectors[target][component][0] = self.covectors[target][component][0]
                            + c * self.residual[input][component][0]
                            + s * self.odd_y[input][component][0];
                        self.covectors[target][component][1] = self.covectors[target][component][1]
                            + s * self.residual[input][component][1]
                            + c * self.odd_y[input][component][1];
                    }
                }
            }
        }

        self.check(n, cancelled)
    }

    // Proper source products use conjugate Gram actions on their complex
    // amplitudes. Each spatial/time axis is separable; the real-field factor
    // 1/2 is applied ONCE after all three passes.
    pub(super) fn proper(
        &mut self,
        spatial: &Spatial<ProperKernel>,
        temporal: &[ProperKernel],
        f: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AngularSpectrumError> {
        let (nx, ny) = (spatial.grid.width, spatial.grid.height);
        for frequency in 0..f {
            for y in 0..ny {
                for x in 0..nx {
                    poll(cancelled)?;
                    let mut sum = [[Dd::ZERO; 2]; 6];
                    for q in 0..nx {
                        if q % 32 == 0 {
                            poll(cancelled)?;
                        }
                        proper_add(
                            &mut sum,
                            spatial.x[x * nx + q],
                            &self.residual[frequency * self.modes + y * nx + q],
                        );
                    }
                    self.covectors[frequency * self.modes + y * nx + x] = sum;
                }
            }
        }
        for frequency in 0..f {
            for y in 0..ny {
                for x in 0..nx {
                    poll(cancelled)?;
                    let mut sum = [[Dd::ZERO; 2]; 6];
                    for q in 0..ny {
                        if q % 32 == 0 {
                            poll(cancelled)?;
                        }
                        proper_add(
                            &mut sum,
                            spatial.y[y * ny + q],
                            &self.covectors[frequency * self.modes + q * nx + x],
                        );
                    }
                    self.odd_x[frequency * self.modes + y * nx + x] = sum;
                }
            }
        }
        for a in 0..f {
            for mode in 0..self.modes {
                poll(cancelled)?;
                let mut sum = [[Dd::ZERO; 2]; 6];
                for b in 0..f {
                    if b % 32 == 0 {
                        poll(cancelled)?;
                    }
                    proper_add(
                        &mut sum,
                        temporal[a * f + b],
                        &self.odd_x[b * self.modes + mode],
                    );
                }
                self.covectors[a * self.modes + mode] =
                    sum.map(|value| value.map(|part| part * Dd::from_f64(0.5)));
            }
        }
        self.check(self.modes * f, cancelled)
    }

    fn check(
        &self,
        n: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AngularSpectrumError> {
        for values in self.covectors[..n].chunks(128) {
            poll(cancelled)?;
            if values.iter().flatten().flatten().any(|&v| !finite_dd(v)) {
                return Err(bad("aperture Gram field is unrepresentable"));
            }
        }
        Ok(())
    }
}
fn proper_add(output: &mut FieldParts, kernel: ProperKernel, input: &FieldParts) {
    let [real, imaginary] = kernel.0;
    for (sum, &[r, i]) in output.iter_mut().zip(input) {
        sum[0] = sum[0] + real * r + imaginary * i;
        sum[1] = sum[1] + real * i - imaginary * r;
    }
}
fn zeros(
    count: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<FieldParts>, AngularSpectrumError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| bad("aperture Gram allocation failed"))?;
    for i in 0..count {
        if i % 256 == 0 {
            poll(cancelled)?;
        }
        result.push([[Dd::ZERO; 2]; 6]);
    }
    Ok(result)
}
pub(super) fn clear(
    fields: &mut [FieldParts],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AngularSpectrumError> {
    for tile in fields.chunks_mut(128) {
        poll(cancelled)?;
        tile.fill([[Dd::ZERO; 2]; 6]);
    }
    Ok(())
}
