//! Dynamic Favre gradient-diffusion transport on the original periodic filter.
//!
//! For resolved scalar `phi`, the modeled flux is
//! `q_i = -rho C_phi Delta^2 |S| grad_i(phi)`. The scalar Germano identity gives
//! `L_i = hat(rho u_i phi) - hat(rho) uhat_i phihat` and
//! `M_i = hat(rho |S| grad_i(phi)) - 9 hat(rho) |Shat| grad_i(phihat)`
//! with derivatives scaled by the grid filter width `Delta`, each using its
//! physical axis spacing. We fit `C_phi = <L.M>/<M.M>`; unlike the stress
//! identity there is no factor two. The physical diffusivity is
//! `alpha_t = C_phi Delta^2 |S|`. It is fitted independently of momentum viscosity,
//! avoiding a singular turbulent-Prandtl ratio when that viscosity vanishes.
//!
//! Reference: Squires, CTR Annual Research Briefs (1991), equations
//! (2.13)--(2.19), https://web.stanford.edu/group/ctr/ResBriefs/1991/19_SQUIRES.pdf;
//! Moin, Squires, Cabot & Lee (1991), doi:10.1063/1.858164. The existing box
//! filter has each test-box width three times its grid-box width (squared
//! geometric-mean ratio nine). Coefficient averaging and
//! nonnegative diffusion are explicit policies. This supplies no wall model,
//! reaction-rate model, variable-property thermodynamic closure, or DNS claim.
//! The sibling [`super::channel`] wrappers reuse this fit with a tangential
//! test filter, plane averaging, and explicit homogeneous scalar wall data.

use super::*;

/// One common diffusivity fitted to several same-unit scalar components.
#[path = "scalar/components.rs"]
pub mod components;

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
    pub(super) shape: [usize; 3],
    walls: [ScalarWall; 2],
    pub(super) boundary_gradient: ScalarBoundaryGradient,
    scalar: Vec<f64>,
    filtered: Vec<f64>,
    pub(super) basis: Vec<[f64; 3]>,
    pub(super) contraction: Vec<[f64; 2]>,
    pub(super) diffusivity: Vec<f64>,
}

/// A flux owner may constrain a coupled flux without prescribing the scalar's
/// normal gradient. Such components use interior reconstruction at both scales.
#[derive(Clone, Copy)]
pub(super) enum ScalarBoundaryGradient {
    Declared,
    Interior,
    /// Primitive values at full-spacing exterior reservoir ghosts.
    Reservoir,
}

#[allow(clippy::too_many_arguments)]
fn gradient(
    shape: [usize; 3],
    metric: [f64; 3],
    geometry: FilterGeometry,
    walls: [ScalarWall; 2],
    boundary_gradient: ScalarBoundaryGradient,
    filtered: bool,
    cell: usize,
    field: &[f64],
) -> Result<[f64; 3], DynamicSgsError> {
    let at = coordinates(shape, cell);
    let mut gradient = [0.0; 3];
    for axis in 0..3 {
        let stencil = geometry.stencil(shape, metric, cell, at, axis);
        gradient[axis] = if matches!(boundary_gradient, ScalarBoundaryGradient::Reservoir) {
            stencil.ghost_derivative(|j| field[j], walls)?
        } else if matches!(boundary_gradient, ScalarBoundaryGradient::Interior) {
            stencil.interior_scalar_derivative(field)?
        } else if filtered {
            stencil.filtered_scalar_derivative(field, walls)?
        } else {
            stencil.derivative(|j| field[j], walls)?
        };
    }
    Ok(gradient)
}

impl PeriodicDynamicScalar {
    pub fn required_bytes(cells: usize) -> Result<usize, DynamicSgsError> {
        cells
            .checked_mul(8 * core::mem::size_of::<f64>())
            .ok_or(DynamicSgsError::Invalid("scalar workspace size overflow"))
    }

    pub fn new(shape: [usize; 3], maximum_bytes: usize) -> Result<Self, DynamicSgsError> {
        Self::new_geometry(shape, FilterGeometry::Periodic, maximum_bytes)
    }

    pub(super) fn new_geometry(
        shape: [usize; 3],
        geometry: FilterGeometry,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        geometry.validate_shape(shape)?;
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
            walls: [ScalarWall::ZeroGradient; 2],
            boundary_gradient: ScalarBoundaryGradient::Declared,
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

    /// Fitted nonnegative C_phi, independently of centered cell strain.
    /// Conservative face consumers must combine it with their current face
    /// strain and the flow's squared filter width. Valid only after successful
    /// combined preparation; out-of-grid indices return None.
    pub fn coefficient(&self, cell: usize) -> Option<f64> {
        self.basis.get(cell).map(|value| value[0])
    }

    /// Allocation-free traversal of successfully fitted C_phi values.
    pub fn coefficients(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.basis.iter().map(|value| value[0])
    }

    pub(super) fn normalize(
        &mut self,
        scalar: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, DynamicSgsError> {
        self.normalize_with_walls(scalar, None, cancelled)
    }

    pub(super) fn normalize_with_walls(
        &mut self,
        scalar: &[f64],
        walls: Option<[ScalarWall; 2]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, DynamicSgsError> {
        self.normalize_boundary(scalar, walls, 0.0, cancelled)
    }

    pub(super) fn normalize_boundary(
        &mut self,
        scalar: &[f64],
        walls: Option<[ScalarWall; 2]>,
        external_extent: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, DynamicSgsError> {
        let mut scale: f64 = external_extent;
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
        if let Some(walls) = walls {
            self.walls = walls;
            for wall in &mut self.walls {
                if let ScalarWall::Value(value) = wall {
                    *value -= reference;
                    if !value.is_finite() {
                        return Err(DynamicSgsError::Invalid(
                            "finite representable scalar wall differences required",
                        ));
                    }
                    scale = scale.max(value.abs());
                }
            }
        }
        if scale != 0.0 {
            for (cell, value) in self.scalar.iter_mut().enumerate() {
                poll(cancelled, cell)?;
                *value /= scale;
            }
            if walls.is_some() {
                for wall in &mut self.walls {
                    if let ScalarWall::Value(value) = wall {
                        *value /= scale;
                    }
                }
            }
        }
        Ok(scale)
    }

    pub(super) fn prepare_from(
        &mut self,
        flow: &PeriodicDynamicSgs,
        config: DynamicScalarConfig,
        density: &[f64],
        scales: InputScales,
        scalar_scale: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DynamicScalarReport, DynamicSgsError> {
        self.prepare_from_profile(flow, config, density, scales, scalar_scale, None, cancelled)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare_from_profile(
        &mut self,
        flow: &PeriodicDynamicSgs,
        config: DynamicScalarConfig,
        density: &[f64],
        scales: InputScales,
        scalar_scale: f64,
        profile: Option<&reservoir::ProfileBoundary<'_>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DynamicScalarReport, DynamicSgsError> {
        let n = self.scalar.len();
        // Momentum deliberately skips all filtered scratch for uniform flow.
        // Scalar input has already been validated even on this exact-zero path.
        if scales.velocity == 0.0 || scalar_scale == 0.0 {
            for (cell, (diffusivity, basis)) in
                self.diffusivity.iter_mut().zip(&mut self.basis).enumerate()
            {
                poll(cancelled, cell)?;
                *diffusivity = 0.0;
                basis[0] = 0.0;
            }
            poll(cancelled, 0)?;
            return Ok(DynamicScalarReport {
                undetermined_cells: n,
                ..DynamicScalarReport::default()
            });
        }
        let domain = self.prepare_contractions_with_profile(
            flow, density, scales, scalar_scale, profile, cancelled,
        )?;
        if !flow.normal_ghosts || config.averaging == DynamicAveraging::Domain {
            flow.geometry
                .average_planes(self.shape, &mut self.contraction, cancelled)?;
        }
        self.fit_from(flow, config, scales, domain, cancelled)
    }

    // Both the single-scalar and common-component fits use this exact kernel.
    // `scalar` has already been centered and normalized by the caller.
    pub(super) fn prepare_contractions(
        &mut self,
        flow: &PeriodicDynamicSgs,
        density: &[f64],
        scales: InputScales,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2], DynamicSgsError> {
        self.prepare_contractions_with_profile(flow, density, scales, 0.0, None, cancelled)
    }

    fn prepare_contractions_with_profile(
        &mut self,
        flow: &PeriodicDynamicSgs,
        density: &[f64],
        scales: InputScales,
        scalar_scale: f64,
        profile: Option<&reservoir::ProfileBoundary<'_>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2], DynamicSgsError> {
        let n = self.scalar.len();
        let samples = flow.geometry.samples();
        let ratio_squared = flow.geometry.ratio_squared();
        for cell in 0..n {
            poll(cancelled, cell)?;
            let walls = profile.map_or(Ok(self.walls), |profile|
                profile.scalar(cell, false, scalar_scale))?;
            let g = gradient(
                self.shape,
                flow.derivative_metric,
                flow.geometry,
                walls,
                self.boundary_gradient,
                false,
                cell,
                &self.scalar,
            )?;
            let rho = density[cell] / scales.density;
            let rate = flow.strain_rate[cell];
            self.basis[cell] = g.map(|v| rho * rate * v);
            if g.iter().zip(self.basis[cell]).any(|(g, b)| {
                !b.is_finite() || (b == 0.0 && rho != 0.0 && rate != 0.0 && *g != 0.0)
            }) {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "unrepresentable normalized scalar basis",
                });
            }
            let mut value = 0.0;
            flow.geometry.each_box(self.shape, cell, |j| {
                value += (density[j] / scales.density / samples) * self.scalar[j];
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
            flow.geometry.each_box(self.shape, cell, |j| {
                let rho = density[j] / scales.density / samples;
                let difference = self.scalar[j] - scalar_mean;
                for axis in 0..3 {
                    l[axis] += rho * (flow.velocity[j][axis] - mean[axis]) * difference;
                    a[axis] += self.basis[j][axis] / samples;
                }
            });
            let wall_velocity = profile.map_or(Ok(flow.wall_velocity), |profile|
                profile.velocity(cell, true, scales.velocity))?;
            let (rate, _) = if flow.normal_ghosts {
                strain_and_gradient(
                    self.shape,
                    flow.derivative_metric,
                    flow.geometry,
                    wall_velocity,
                    true,
                    cell,
                    &flow.filtered,
                )
                .map(|(rate, strain, _)| (rate, strain))
            } else {
                strain(
                    self.shape,
                    flow.derivative_metric,
                    flow.geometry,
                    flow.wall_velocity,
                    cell,
                    &flow.filtered,
                )
            }?;
            let walls = profile.map_or(Ok(self.walls), |profile|
                profile.scalar(cell, true, scalar_scale))?;
            let g = gradient(
                self.shape,
                flow.derivative_metric,
                flow.geometry,
                walls,
                self.boundary_gradient,
                true,
                cell,
                &self.filtered,
            )?;
            let mut pair = [0.0; 2];
            let mut modeled = false;
            for axis in 0..3 {
                let m = a[axis] - ratio_squared * mean[3] * rate * g[axis];
                modeled |= m != 0.0;
                pair[0] += l[axis] * m;
                pair[1] += m * m;
            }
            if pair.iter().any(|v| !v.is_finite()) || (pair[1] == 0.0 && modeled) {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "unrepresentable scalar least-squares contraction",
                });
            }
            self.contraction[cell] = pair;
            numerator.add(pair[0]);
            denominator.add(pair[1]);
        }
        Ok([numerator.value, denominator.value])
    }

    fn fit_from(
        &mut self,
        flow: &PeriodicDynamicSgs,
        config: DynamicScalarConfig,
        scales: InputScales,
        domain: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DynamicScalarReport, DynamicSgsError> {
        let n = self.scalar.len();
        let mut report = DynamicScalarReport {
            minimum_raw_coefficient: f64::INFINITY,
            ..DynamicScalarReport::default()
        };
        for cell in 0..n {
            poll(cancelled, cell)?;
            let pair = match config.averaging {
                DynamicAveraging::Domain => match flow.geometry {
                    FilterGeometry::Periodic => domain,
                    FilterGeometry::Channel(_) => self.contraction[cell],
                },
                DynamicAveraging::LocalBox => {
                    let mut pair = [0.0; 2];
                    let mut determined = false;
                    let samples = flow.averaging_samples();
                    flow.each_averaging_box(cell, |j| {
                        determined |= self.contraction[j][1] != 0.0;
                        for k in 0..2 {
                            pair[k] += self.contraction[j][k] / samples;
                        }
                    });
                    if determined && pair[1] == 0.0 {
                        return Err(DynamicSgsError::Cell {
                            cell,
                            reason: "underflow of local scalar least-squares denominator",
                        });
                    }
                    pair
                }
            };
            if pair.iter().any(|v| !v.is_finite()) || (pair[1] == 0.0 && pair[0] != 0.0) {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "unrepresentable averaged scalar least-squares contraction",
                });
            }
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
            let diffusivity =
                coefficient * flow.strain_rate[cell] * scales.velocity * flow.filter_width;
            if !diffusivity.is_finite() {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "unrepresentable modeled scalar diffusivity",
                });
            }
            self.diffusivity[cell] = diffusivity;
            // Contractions have consumed the vector basis. Keep C_phi in its
            // first lane even when centered strain makes cell diffusion zero.
            self.basis[cell][0] = coefficient;
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
