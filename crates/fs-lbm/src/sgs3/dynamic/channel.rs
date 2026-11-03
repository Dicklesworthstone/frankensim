//! Dynamic momentum and scalar transport for a resolved planar channel.
//!
//! The active cell centers occupy a uniform Cartesian box between two
//! impermeable, no-slip walls half a normal cell outside the first/last center.
//! The other two axes are periodic and statistically homogeneous. Both wall
//! velocities are spatially uniform; tangential motion is allowed. Scalar wall
//! data are independently uniform fixed values or zero normal derivatives.
//!
//! Only the two tangential directions are test filtered: nine equal cell
//! averages, with widths `(3hp, 3hq, hn)`. As on the periodic path the model
//! length is the volume width `Delta = (hp hq hn)^(1/3)`, so the squared test
//! ratio is `9^(2/3) = 3^(4/3)`, NOT nine. The numerator and denominator are
//! averaged independently within each entire homogeneous plane. No averaging
//! joins different wall distances, and no test sample crosses a wall.
//!
//! Cabot, CTR Annual Research Briefs 1991, pp.191--205, §2.1, explicitly
//! compares this volume-width convention with a two-dimensional width for
//! homogeneous-direction filtering and plane-averaged dynamic scalar fits:
//! <https://web.stanford.edu/group/ctr/ResBriefs/1991/18_CABOT.pdf>.
//! Our finite-volume box and half-cell derivatives are a declared discretization,
//! not a reproduction of that paper's spectral channel calibration.
//!
//! Constant tangential weights avoid the variable-width/truncated-filter
//! commutation terms of a wall-normal test filter. Uniform Dirichlet wall
//! values also survive Favre filtering. A physical Neumann condition generally
//! does NOT survive density-weighted filtering: the filtered scalar therefore
//! uses a second-order, three-fluid-point normal derivative at such walls,
//! without imposing a wall condition on the Favre quotient. Its physical
//! counterpart uses the declared zero wall derivative. This retains derivative
//! truncation error; it does not claim exact discrete Favre commutation or
//! remove implicit grid-filter effects. Nor does it model curved walls, corners, variable grid
//! widths, spatially varying wall data, wall functions, or near-wall resolution
//! error. Scalar eddy diffusion and isotropic eddy viscosity retain their
//! scale-similarity and aspect-ratio limitations; this is exploratory LES.

use super::scalar::{DynamicScalarConfig, DynamicScalarReport, PeriodicDynamicScalar};
use super::*;

/// Common component diffusion without an invented composition wall gradient.
#[path = "channel/components.rs"]
pub mod components;

/// Spatially uniform scalar data at one physical half-cell wall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ScalarWall {
    /// Prescribed value in the same units as the cell scalar.
    Value(f64),
    /// Zero physical normal derivative (an insulated/impermeable scalar wall).
    /// Requires at least three normal fluid cells for the filtered quotient's
    /// independent one-sided derivative; it is not assigned a Neumann datum.
    ZeroGradient,
}

/// Momentum controls; fitting over each homogeneous plane is intrinsic.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChannelDynamicSgsConfig {
    pub maximum_coefficient_squared: f64,
    pub maximum_eddy_to_molecular: f64,
    pub usage: SgsUse,
}

impl ChannelDynamicSgsConfig {
    pub fn validate(self) -> Result<(), DynamicSgsError> {
        self.fit_config().validate()
    }
    fn fit_config(self) -> DynamicSgsConfig {
        DynamicSgsConfig {
            averaging: DynamicAveraging::Domain,
            maximum_coefficient_squared: self.maximum_coefficient_squared,
            maximum_eddy_to_molecular: self.maximum_eddy_to_molecular,
            usage: self.usage,
        }
    }
}

/// Independent scalar controls; the native owner enforces the molecular cap.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChannelDynamicScalarConfig {
    pub maximum_coefficient: f64,
    pub maximum_diffusivity_to_molecular: f64,
    pub usage: SgsUse,
}

impl ChannelDynamicScalarConfig {
    pub fn validate(self) -> Result<(), DynamicSgsError> {
        self.fit_config().validate()
    }
    fn fit_config(self) -> DynamicScalarConfig {
        DynamicScalarConfig {
            averaging: DynamicAveraging::Domain,
            maximum_coefficient: self.maximum_coefficient,
            maximum_diffusivity_to_molecular: self.maximum_diffusivity_to_molecular,
            usage: self.usage,
        }
    }
}

/// The same 136 numerical bytes/cell as the periodic core. Only active fluid
/// cells are stored, x-fastest; wall planes are not part of the input fields.
/// Preparation allocates nothing, and outputs are valid only after success.
pub struct ChannelDynamicSgs {
    inner: PeriodicDynamicSgs,
    normal_axis: usize,
}

/// Independent scalar fit with 64 numerical bytes/cell, excluding metadata.
pub struct ChannelDynamicScalar {
    inner: PeriodicDynamicScalar,
    normal_axis: usize,
}

impl ChannelDynamicSgs {
    pub fn required_bytes(cells: usize) -> Result<usize, DynamicSgsError> {
        PeriodicDynamicSgs::required_bytes(cells)
    }

    pub fn new(
        shape: [usize; 3],
        spacing: f64,
        normal_axis: usize,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        Self::new_cartesian(shape, [spacing; 3], normal_axis, maximum_bytes)
    }

    /// At least two fluid cells normal to the walls and three in each periodic
    /// direction. The grid and wall normal remain fixed for this workspace.
    pub fn new_cartesian(
        shape: [usize; 3],
        spacing: [f64; 3],
        normal_axis: usize,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        let geometry = FilterGeometry::Channel(normal_axis);
        geometry.validate_shape(shape)?;
        Ok(Self {
            inner: PeriodicDynamicSgs::new_cartesian_geometry(
                shape,
                spacing,
                geometry,
                maximum_bytes,
            )?,
            normal_axis,
        })
    }

    pub fn workspace_bytes(&self) -> usize {
        self.inner.workspace_bytes()
    }
    pub const fn shape(&self) -> [usize; 3] {
        self.inner.shape
    }
    pub const fn spacing(&self) -> [f64; 3] {
        self.inner.spacing
    }
    pub const fn normal_axis(&self) -> usize {
        self.normal_axis
    }
    pub const fn filter_width(&self) -> f64 {
        self.inner.filter_width
    }
    pub fn viscosities(&self) -> &[f64] {
        self.inner.viscosities()
    }

    /// Fit from the current physical flow and both current wall velocities.
    /// Each wall's normal velocity must be zero. Wall values participate in
    /// reference subtraction and normalization, including uniform fluid fields.
    pub fn prepare(
        &mut self,
        config: ChannelDynamicSgsConfig,
        density: &[f64],
        velocity: &[[f64; 3]],
        wall_velocity: [[f64; 3]; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<DynamicSgsReport, DynamicSgsError> {
        self.inner
            .prepare_scaled_with_walls(
                config.fit_config(),
                density,
                velocity,
                wall_velocity,
                cancelled,
            )
            .map(|(report, _)| report)
    }

    /// Fit momentum and scalar diffusion from one same-time field. No public
    /// path can consume stale normalized flow scratch. Both outputs require
    /// success, and any refusal leaves the caller's source fields unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_with_scalar(
        &mut self,
        momentum: ChannelDynamicSgsConfig,
        heat: ChannelDynamicScalarConfig,
        density: &[f64],
        velocity: &[[f64; 3]],
        wall_velocity: [[f64; 3]; 2],
        scalar: &[f64],
        scalar_walls: [ScalarWall; 2],
        scalar_workspace: &mut ChannelDynamicScalar,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(DynamicSgsReport, DynamicScalarReport), DynamicSgsError> {
        poll(&mut cancelled, 0)?;
        heat.validate()?;
        if scalar_workspace.normal_axis != self.normal_axis
            || scalar_workspace.inner.shape != self.inner.shape
            || scalar.len() != self.inner.velocity.len()
        {
            return Err(DynamicSgsError::Invalid(
                "channel scalar and flow geometry differ",
            ));
        }
        if scalar_walls.contains(&ScalarWall::ZeroGradient) && self.shape()[self.normal_axis] < 3 {
            return Err(DynamicSgsError::Invalid(
                "Neumann scalar filtering requires three normal fluid cells",
            ));
        }
        let scalar_scale = scalar_workspace.inner.normalize_with_walls(
            scalar,
            Some(scalar_walls),
            &mut cancelled,
        )?;
        let (momentum_report, scales) = self.inner.prepare_scaled_with_walls(
            momentum.fit_config(),
            density,
            velocity,
            wall_velocity,
            &mut cancelled,
        )?;
        let scalar_report = scalar_workspace.inner.prepare_from(
            &self.inner,
            heat.fit_config(),
            density,
            scales,
            scalar_scale,
            &mut cancelled,
        )?;
        Ok((momentum_report, scalar_report))
    }
}

impl ChannelDynamicScalar {
    pub fn required_bytes(cells: usize) -> Result<usize, DynamicSgsError> {
        PeriodicDynamicScalar::required_bytes(cells)
    }
    pub fn new(
        shape: [usize; 3],
        normal_axis: usize,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        let geometry = FilterGeometry::Channel(normal_axis);
        geometry.validate_shape(shape)?;
        Ok(Self {
            inner: PeriodicDynamicScalar::new_geometry(shape, geometry, maximum_bytes)?,
            normal_axis,
        })
    }
    pub fn workspace_bytes(&self) -> usize {
        self.inner.workspace_bytes()
    }
    pub fn diffusivities(&self) -> &[f64] {
        self.inner.diffusivities()
    }
}

#[derive(Clone, Copy)]
pub(super) enum FilterGeometry {
    Periodic,
    Channel(usize),
}

// One stack-only stencil per derivative axis, reused across velocity components.
// Address calculation and periodic integer divisions stay outside component loops.
#[derive(Clone, Copy)]
pub(super) struct AxisStencil {
    cell: usize,
    metric: f64,
    span: StencilSpan,
}

#[derive(Clone, Copy)]
enum StencilSpan {
    Centered {
        low: usize,
        high: usize,
    },
    Wall {
        near: usize,
        far: Option<usize>,
        side: usize,
    },
}

impl AxisStencil {
    fn checked(self, unscaled: f64, derivative: f64) -> Result<f64, DynamicSgsError> {
        if !derivative.is_finite() || (derivative == 0.0 && unscaled != 0.0) {
            return Err(DynamicSgsError::Cell {
                cell: self.cell,
                reason: "unrepresentable normalized derivative",
            });
        }
        Ok(derivative)
    }

    pub(super) fn derivative(
        self,
        value: impl Fn(usize) -> f64,
        walls: [ScalarWall; 2],
    ) -> Result<f64, DynamicSgsError> {
        match self.span {
            StencilSpan::Centered { low, high } => {
                let difference = value(high) - value(low);
                self.checked(difference, 0.5 * difference * self.metric)
            }
            StencilSpan::Wall { near, side, .. } => {
                let center = value(self.cell);
                let neighbor = value(near);
                let sign = if side == 0 { 1.0 } else { -1.0 };
                let unscaled = match walls[side] {
                    ScalarWall::Value(wall) => sign * ((center - wall) + (neighbor - wall) / 3.0),
                    ScalarWall::ZeroGradient => sign * 0.5 * (neighbor - center),
                };
                self.checked(unscaled, unscaled * self.metric)
            }
        }
    }

    pub(super) fn filtered_scalar_derivative(
        self,
        field: &[f64],
        walls: [ScalarWall; 2],
    ) -> Result<f64, DynamicSgsError> {
        if let StencilSpan::Wall { side, .. } = self.span {
            if walls[side] == ScalarWall::ZeroGradient {
                return self.interior_scalar_derivative(field);
            }
        }
        self.derivative(|j| field[j], walls)
    }

    // No scalar boundary datum is inferred. This also handles physical
    // composition gradients when an impermeable wall constrains a coupled
    // molecular flux (including thermal diffusion), rather than grad(Y).
    pub(super) fn interior_scalar_derivative(self, field: &[f64]) -> Result<f64, DynamicSgsError> {
        if let StencilSpan::Wall { near, far, side } = self.span {
            let far = far.ok_or(DynamicSgsError::Invalid(
                "interior scalar wall derivatives require three normal fluid cells",
            ))?;
            let d1 = field[near] - field[self.cell];
            let d2 = field[far] - field[self.cell];
            let difference = if side == 0 {
                4.0 * d1 - d2
            } else {
                d2 - 4.0 * d1
            };
            return self.checked(difference, 0.5 * difference * self.metric);
        }
        self.derivative(|j| field[j], [ScalarWall::ZeroGradient; 2])
    }
}

impl FilterGeometry {
    pub(super) fn validate_shape(self, shape: [usize; 3]) -> Result<(), DynamicSgsError> {
        let normal = match self {
            Self::Periodic => None,
            Self::Channel(axis) if axis < 3 => Some(axis),
            Self::Channel(_) => {
                return Err(DynamicSgsError::Invalid(
                    "channel normal axis must be 0, 1 or 2",
                ));
            }
        };
        if shape
            .iter()
            .enumerate()
            .any(|(axis, n)| *n < if Some(axis) == normal { 2 } else { 3 })
        {
            return Err(DynamicSgsError::Invalid(
                "three cells per periodic axis and two per wall-normal axis required",
            ));
        }
        Ok(())
    }

    pub(super) fn samples(self) -> f64 {
        match self {
            Self::Periodic => 27.0,
            Self::Channel(_) => 9.0,
        }
    }

    pub(super) fn ratio_squared(self) -> f64 {
        match self {
            Self::Periodic => 9.0,
            Self::Channel(_) => det::exp(det::ln(9.0) * (2.0 / 3.0)),
        }
    }

    pub(super) fn each_box(self, shape: [usize; 3], cell: usize, mut visit: impl FnMut(usize)) {
        match self {
            Self::Periodic => each_box(shape, cell, visit),
            Self::Channel(normal) => {
                let p = coordinates(shape, cell);
                let mut sides =
                    core::array::from_fn::<_, 3, _>(|axis| neighbors(shape[axis], p[axis]));
                sides[normal] = [p[normal]; 3];
                let lengths =
                    core::array::from_fn::<_, 3, _>(|axis| if axis == normal { 1 } else { 3 });
                for z in 0..lengths[2] {
                    for y in 0..lengths[1] {
                        for x in 0..lengths[0] {
                            visit(index(shape, [sides[0][x], sides[1][y], sides[2][z]]));
                        }
                    }
                }
            }
        }
    }

    // Exact derivative of the quadratic through a half-cell Dirichlet wall,
    // this center and its interior neighbor. ZeroGradient instead supplies
    // the derivative at the wall; no periodic ghost is ever used normally.
    pub(super) fn stencil(
        self,
        shape: [usize; 3],
        metric: [f64; 3],
        cell: usize,
        p: [usize; 3],
        axis: usize,
    ) -> AxisStencil {
        let boundary = matches!(self, Self::Channel(normal) if normal == axis)
            && (p[axis] == 0 || p[axis] + 1 == shape[axis]);
        let span = if boundary {
            let low = p[axis] == 0;
            let side = if low { 0 } else { 1 };
            let mut next = p;
            next[axis] = if low { 1 } else { shape[axis] - 2 };
            let far = if shape[axis] >= 3 {
                let mut far = p;
                far[axis] = if low { 2 } else { shape[axis] - 3 };
                Some(index(shape, far))
            } else {
                None
            };
            StencilSpan::Wall {
                near: index(shape, next),
                far,
                side,
            }
        } else {
            let sides = neighbors(shape[axis], p[axis]);
            let mut low = p;
            let mut high = p;
            low[axis] = sides[0];
            high[axis] = sides[2];
            StencilSpan::Centered {
                low: index(shape, low),
                high: index(shape, high),
            }
        };
        AxisStencil {
            cell,
            metric: metric[axis],
            span,
        }
    }

    // Reuse contractions as plane means: O(cells), no additional allocation or
    // wall-normal-length-dependent memory. Sum before division, not coefficients.
    pub(super) fn average_planes(
        self,
        shape: [usize; 3],
        pairs: &mut [[f64; 2]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DynamicSgsError> {
        let Self::Channel(normal) = self else {
            return Ok(());
        };
        let axes = match normal {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        };
        let plane_cells = (shape[axes[0]] * shape[axes[1]]) as f64;
        let mut visited = 0;
        for layer in 0..shape[normal] {
            let mut sum = [Sum::default(), Sum::default()];
            let mut determined = false;
            for j in 0..shape[axes[1]] {
                for i in 0..shape[axes[0]] {
                    poll(cancelled, visited)?;
                    visited += 1;
                    let mut p = [0; 3];
                    p[normal] = layer;
                    p[axes[0]] = i;
                    p[axes[1]] = j;
                    let pair = pairs[index(shape, p)];
                    determined |= pair[1] != 0.0;
                    sum[0].add(pair[0] / plane_cells);
                    sum[1].add(pair[1] / plane_cells);
                }
            }
            if determined && sum[1].value == 0.0 {
                return Err(DynamicSgsError::Invalid(
                    "underflow of plane least-squares denominator",
                ));
            }
            for j in 0..shape[axes[1]] {
                for i in 0..shape[axes[0]] {
                    poll(cancelled, visited)?;
                    visited += 1;
                    let mut p = [0; 3];
                    p[normal] = layer;
                    p[axes[0]] = i;
                    p[axes[1]] = j;
                    pairs[index(shape, p)] = [sum[0].value, sum[1].value];
                }
            }
        }
        poll(cancelled, 0)
    }
}

#[cfg(test)]
#[path = "channel/tests.rs"]
mod tests;
