//! Dynamic Favre gradient-diffusion transport on the original periodic filter.
//!
//! For resolved scalar `phi`, the modeled flux is
//! `q_i = -rho C_phi h^2 |S| grad_i(phi)`. The scalar Germano identity gives
//! `L_i = hat(rho u_i phi) - hat(rho) uhat_i phihat` and
//! `M_i = hat(rho |S| grad_i(phi)) - 9 hat(rho) |Shat| grad_i(phihat)`
//! in cell coordinates. We fit `C_phi = <L.M>/<M.M>`; unlike the stress
//! identity there is no factor two. The physical diffusivity is
//! `alpha_t = C_phi h^2 |S|`. It is fitted independently of momentum viscosity,
//! avoiding a singular turbulent-Prandtl ratio when that viscosity vanishes.
//!
//! Reference: Squires, CTR Annual Research Briefs (1991), equations
//! (2.13)--(2.19), https://web.stanford.edu/group/ctr/ResBriefs/1991/19_SQUIRES.pdf;
//! Moin, Squires, Cabot & Lee (1991), doi:10.1063/1.858164. The existing box
//! filter has test width 3h (squared ratio nine). Coefficient averaging and
//! nonnegative diffusion are explicit policies. This supplies no wall model,
//! reaction-rate model, variable-property thermodynamic closure, or DNS claim.

use super::*;

/// Explicit scalar-fit controls. The consuming transport owner enforces the
/// diffusivity/molecular bound using its actual material and face coefficients.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DynamicScalarConfig {
    pub averaging: DynamicAveraging,
    /// Admission bound on C_phi itself, not on its square or square root.
    pub maximum_coefficient: f64,
    pub maximum_diffusivity_to_molecular: f64,
    pub usage: SgsUse,
}

impl DynamicScalarConfig {
    pub fn validate(self) -> Result<(), DynamicSgsError> {
        if self.usage == SgsUse::AcceptanceRender {
            return Err(DynamicSgsError::Invalid(
                "scalar acceptance use is not supported",
            ));
        }
        if !self.maximum_coefficient.is_finite()
            || self.maximum_coefficient < 0.0
            || !self.maximum_diffusivity_to_molecular.is_finite()
            || self.maximum_diffusivity_to_molecular < 0.0
        {
            return Err(DynamicSgsError::Invalid(
                "finite nonnegative scalar bounds required",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DynamicScalarReport {
    pub maximum_coefficient: f64,
    pub minimum_raw_coefficient: f64,
    /// Cells where negative modeled diffusivity was explicitly excluded.
    pub backscatter_limited_cells: usize,
    pub undetermined_cells: usize,
    /// Same length^2/time units as the returned modeled momentum viscosity.
    pub maximum_diffusivity: f64,
}

/// Reused scalar arrays: eight f64 values per cell, excluding fixed metadata.
/// Cell order and topology match [`PeriodicDynamicSgs`]. Scratch/outputs may
/// be consumed only after the combined preparation returns success.
pub struct PeriodicDynamicScalar {
    shape: [usize; 3],
    scalar: Vec<f64>,
    filtered: Vec<f64>,
    basis: Vec<[f64; 3]>,
    contraction: Vec<[f64; 2]>,
    diffusivity: Vec<f64>,
}

fn gradient(shape: [usize; 3], cell: usize, field: &[f64]) -> [f64; 3] {
    let p = coordinates(shape, cell);
    core::array::from_fn(|axis| {
        let sides = neighbors(shape[axis], p[axis]);
        let mut low = p;
        let mut high = p;
        low[axis] = sides[0];
        high[axis] = sides[2];
        0.5 * (field[index(shape, high)] - field[index(shape, low)])
    })
}

impl PeriodicDynamicScalar {
    pub fn required_bytes(cells: usize) -> Result<usize, DynamicSgsError> {
        cells
            .checked_mul(8 * core::mem::size_of::<f64>())
            .ok_or(DynamicSgsError::Invalid("scalar workspace size overflow"))
    }

    pub fn new(shape: [usize; 3], maximum_bytes: usize) -> Result<Self, DynamicSgsError> {
        if shape.iter().any(|&n| n < 3) {
            return Err(DynamicSgsError::Invalid(
                "at least three cells per scalar periodic axis required",
            ));
        }
        let n = shape
            .into_iter()
            .try_fold(1usize, |n, d| n.checked_mul(d))
            .ok_or(DynamicSgsError::Invalid("scalar cell count overflow"))?;
        let required = Self::required_bytes(n)?;
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            });
        }
        let result = Self {
            shape,
            scalar: allocated(n, 0.0)?,
            filtered: allocated(n, 0.0)?,
            basis: allocated(n, [0.0; 3])?,
            contraction: allocated(n, [0.0; 2])?,
            diffusivity: allocated(n, 0.0)?,
        };
        let required = result.workspace_bytes();
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            });
        }
        Ok(result)
    }

    pub fn workspace_bytes(&self) -> usize {
        (self.scalar.capacity()
            + self.filtered.capacity()
            + 3 * self.basis.capacity()
            + 2 * self.contraction.capacity()
            + self.diffusivity.capacity())
            * core::mem::size_of::<f64>()
    }

    pub fn diffusivities(&self) -> &[f64] {
        &self.diffusivity
    }

    fn normalize(
        &mut self,
        scalar: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, DynamicSgsError> {
        let mut scale: f64 = 0.0;
        let reference = scalar[0];
        for (cell, (&value, normalized)) in scalar.iter().zip(&mut self.scalar).enumerate() {
            poll(cancelled, cell)?;
            let difference = value - reference;
            if !difference.is_finite() {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "finite representable scalar differences required",
                });
            }
            *normalized = difference;
            scale = scale.max(difference.abs());
        }
        if scale != 0.0 {
            for (cell, value) in self.scalar.iter_mut().enumerate() {
                poll(cancelled, cell)?;
                *value /= scale;
            }
        }
        Ok(scale)
    }

    fn prepare_from(
        &mut self,
        flow: &PeriodicDynamicSgs,
        config: DynamicScalarConfig,
        density: &[f64],
        scales: InputScales,
        scalar_scale: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DynamicScalarReport, DynamicSgsError> {
        let n = self.scalar.len();
        // Momentum deliberately skips all filtered scratch for uniform flow.
        // Scalar input has already been validated even on this exact-zero path.
        if scales.velocity == 0.0 || scalar_scale == 0.0 {
            for (tile, values) in self.diffusivity.chunks_mut(64).enumerate() {
                poll(cancelled, tile * 64)?;
                values.fill(0.0);
            }
            poll(cancelled, 0)?;
            return Ok(DynamicScalarReport {
                undetermined_cells: n,
                ..DynamicScalarReport::default()
            });
        }
        for cell in 0..n {
            poll(cancelled, cell)?;
            self.basis[cell] = gradient(self.shape, cell, &self.scalar)
                .map(|v| (density[cell] / scales.density) * flow.strain_rate[cell] * v);
            let mut value = 0.0;
            each_box(self.shape, cell, |j| {
                value += (density[j] / scales.density / 27.0) * self.scalar[j];
            });
            self.filtered[cell] = value / flow.filtered[cell][3];
        }
        let mut numerator = Sum::default();
        let mut denominator = Sum::default();
        for cell in 0..n {
            poll(cancelled, cell)?;
            let mean = flow.filtered[cell];
            let scalar_mean = self.filtered[cell];
            let mut a = [0.0; 3];
            let mut l = [0.0; 3];
            each_box(self.shape, cell, |j| {
                let rho = density[j] / scales.density / 27.0;
                let difference = self.scalar[j] - scalar_mean;
                for axis in 0..3 {
                    l[axis] += rho * (flow.velocity[j][axis] - mean[axis]) * difference;
                    a[axis] += self.basis[j][axis] / 27.0;
                }
            });
            let (rate, _) = strain(self.shape, cell, &flow.filtered);
            let g = gradient(self.shape, cell, &self.filtered);
            let mut pair = [0.0; 2];
            for axis in 0..3 {
                let m = a[axis] - 9.0 * mean[3] * rate * g[axis];
                pair[0] += l[axis] * m;
                pair[1] += m * m;
            }
            self.contraction[cell] = pair;
            numerator.add(pair[0]);
            denominator.add(pair[1]);
        }
        let domain = [numerator.value, denominator.value];
        let mut report = DynamicScalarReport {
            minimum_raw_coefficient: f64::INFINITY,
            ..DynamicScalarReport::default()
        };
        for cell in 0..n {
            poll(cancelled, cell)?;
            let pair = match config.averaging {
                DynamicAveraging::Domain => domain,
                DynamicAveraging::LocalBox => {
                    let mut pair = [0.0; 2];
                    each_box(self.shape, cell, |j| {
                        for k in 0..2 {
                            pair[k] += self.contraction[j][k] / 27.0;
                        }
                    });
                    pair
                }
            };
            let raw = if pair[1] == 0.0 {
                report.undetermined_cells += 1;
                0.0
            } else {
                pair[0] / pair[1]
            };
            if !raw.is_finite() {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "unrepresentable scalar least-squares coefficient",
                });
            }
            report.minimum_raw_coefficient = report.minimum_raw_coefficient.min(raw);
            if raw < 0.0 {
                report.backscatter_limited_cells += 1;
            }
            let coefficient = raw.max(0.0);
            if coefficient > config.maximum_coefficient {
                return Err(DynamicSgsError::Coefficient {
                    cell,
                    observed: coefficient,
                    maximum: config.maximum_coefficient,
                });
            }
            let diffusivity = coefficient * flow.strain_rate[cell] * scales.velocity * flow.spacing;
            if !diffusivity.is_finite() {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "unrepresentable modeled scalar diffusivity",
                });
            }
            self.diffusivity[cell] = diffusivity;
            report.maximum_coefficient = report.maximum_coefficient.max(coefficient);
            report.maximum_diffusivity = report.maximum_diffusivity.max(diffusivity);
        }
        poll(cancelled, 0)?;
        Ok(report)
    }
}

impl PeriodicDynamicSgs {
    /// Fit momentum and scalar diffusivity from the SAME current fields. This
    /// reuses normalized velocity, Favre velocity/density and grid strain;
    /// no independent velocity sweep or public stale-scratch API is exposed.
    /// Neither field is consumable if either fit or cancellation refuses.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_with_scalar(
        &mut self,
        momentum: DynamicSgsConfig,
        scalar_config: DynamicScalarConfig,
        density: &[f64],
        velocity: &[[f64; 3]],
        scalar: &[f64],
        scalar_workspace: &mut PeriodicDynamicScalar,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(DynamicSgsReport, DynamicScalarReport), DynamicSgsError> {
        poll(&mut cancelled, 0)?;
        scalar_config.validate()?;
        if scalar_workspace.shape != self.shape || scalar.len() != self.velocity.len() {
            return Err(DynamicSgsError::Invalid(
                "scalar field and flow workspace shapes differ",
            ));
        }
        let scalar_scale = scalar_workspace.normalize(scalar, &mut cancelled)?;
        let (momentum_report, scales) =
            self.prepare_scaled(momentum, density, velocity, &mut cancelled)?;
        let scalar_report = scalar_workspace.prepare_from(
            self,
            scalar_config,
            density,
            scales,
            scalar_scale,
            &mut cancelled,
        )?;
        Ok((momentum_report, scalar_report))
    }
}

#[cfg(test)]
#[path = "scalar/tests.rs"]
mod tests;
