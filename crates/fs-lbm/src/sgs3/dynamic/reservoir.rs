//! Dynamic momentum and scalar mixing between two prescribed reservoirs.
//!
//! Two transverse axes are periodic. At each end of the remaining axis the
//! caller supplies a primitive ghost one FULL cell spacing
//! from the boundary center. Its normal velocity can carry flow in or out.
//! Uniform reservoirs use the same exterior primitive at both filter levels.
//! Spatial profiles use each physical ghost and the corresponding transverse
//! density-weighted test-filtered ghost. This preserves exterior covariance
//! instead of silently treating a varying inlet as a uniform reservoir.
//!
//! The original equal-volume nine-cell transverse filter is retained, with
//! test widths `(hn, 3hp, 3hq)` and squared volume-width ratio `9^(2/3)`.
//! No filter sample crosses either open plane. This is the declared transverse
//! filter convention discussed by Cabot (1991), §2.1, not an isotropic 27-cell
//! filter near an opening: <https://web.stanford.edu/group/ctr/ResBriefs/1991/18_CABOT.pdf>.
//! Domain averaging means one coefficient per transverse plane; LocalBox
//! averages contractions only over the transverse nine-cell neighborhood.
//! Neither operation wraps the nonperiodic coordinate or mixes the two ends.
//!
//! This owner fits Smagorinsky momentum and an independent gradient-diffusion
//! scalar model using the original Germano kernels. The native solver owns
//! reservoir density, pressure, mass flow, energy and the actual boundary flux.
//! No characteristic boundary, localized opening, unresolved inflow spectrum,
//! radiation condition, SGS-energy storage or DNS calibration is supplied.

use super::scalar::{
    DynamicScalarConfig, DynamicScalarReport, PeriodicDynamicScalar, ScalarBoundaryGradient,
};
use super::*;

#[path = "reservoir/profile.rs"]
mod profile;
pub use profile::ReservoirPrimitive;
pub(super) use profile::ProfileBoundary;

/// Reused 136 numerical bytes per fluid cell; no halo arrays are constructed.
pub struct ReservoirDynamicSgs {
    inner: PeriodicDynamicSgs,
    normal_axis: usize,
}

/// Independent scalar transport with 64 numerical bytes per fluid cell.
pub struct ReservoirDynamicScalar {
    inner: PeriodicDynamicScalar,
    normal_axis: usize,
}

impl ReservoirDynamicSgs {
    pub fn required_bytes(cells: usize) -> Result<usize, DynamicSgsError> {
        PeriodicDynamicSgs::required_bytes(cells)
    }

    /// At least two fluid cells along the reservoir normal and three on each
    /// periodic transverse axis. Length units must match the supplied velocity.
    pub fn new(
        shape: [usize; 3],
        spacing: [f64; 3],
        normal_axis: usize,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        let mut inner = PeriodicDynamicSgs::new_cartesian_geometry(
            shape,
            spacing,
            FilterGeometry::Channel(normal_axis),
            maximum_bytes,
        )?;
        inner.normal_ghosts = true;
        Ok(Self { inner, normal_axis })
    }

    pub fn workspace_bytes(&self) -> usize {
        self.inner.workspace_bytes()
    }
    pub const fn shape(&self) -> [usize; 3] {
        self.inner.shape
    }
    pub const fn normal_axis(&self) -> usize {
        self.normal_axis
    }
    pub const fn spacing(&self) -> [f64; 3] {
        self.inner.spacing
    }
    pub const fn filter_width(&self) -> f64 {
        self.inner.filter_width
    }
    pub fn viscosities(&self) -> &[f64] {
        self.inner.viscosities()
    }
    pub fn coefficient_squared(&self, cell: usize) -> Option<f64> {
        self.inner.coefficient_squared(cell)
    }
    pub fn coefficients_squared(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.inner.coefficients_squared()
    }

    /// Fit from the current physical field and low/high full-spacing ghosts.
    /// Scratch and outputs are consumable only after success. No source field
    /// is mutated; preparation allocates no numerical arrays.
    pub fn prepare(
        &mut self,
        config: DynamicSgsConfig,
        density: &[f64],
        velocity: &[[f64; 3]],
        ghost_velocity: [[f64; 3]; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<DynamicSgsReport, DynamicSgsError> {
        self.inner
            .prepare_scaled_with_walls(config, density, velocity, ghost_velocity, cancelled)
            .map(|(report, _)| report)
    }

    /// Same-state independent scalar fit. The two scalar ghosts use the same
    /// full-spacing planes as the velocity ghosts, not halfway-wall values.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_with_scalar(
        &mut self,
        momentum: DynamicSgsConfig,
        heat: DynamicScalarConfig,
        density: &[f64],
        velocity: &[[f64; 3]],
        ghost_velocity: [[f64; 3]; 2],
        scalar: &[f64],
        ghost_scalar: [f64; 2],
        scalar_workspace: &mut ReservoirDynamicScalar,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(DynamicSgsReport, DynamicScalarReport), DynamicSgsError> {
        poll(&mut cancelled, 0)?;
        heat.validate()?;
        if scalar_workspace.normal_axis != self.normal_axis
            || scalar_workspace.inner.shape != self.inner.shape
            || scalar.len() != self.inner.velocity.len()
        {
            return Err(DynamicSgsError::Invalid(
                "reservoir scalar and momentum geometry differ",
            ));
        }
        let scalar_scale = scalar_workspace.inner.normalize_with_walls(
            scalar,
            Some(ghost_scalar.map(ScalarWall::Value)),
            &mut cancelled,
        )?;
        let (momentum_report, scales) = self.inner.prepare_scaled_with_walls(
            momentum,
            density,
            velocity,
            ghost_velocity,
            &mut cancelled,
        )?;
        let scalar_report = scalar_workspace.inner.prepare_from(
            &self.inner,
            heat,
            density,
            scales,
            scalar_scale,
            &mut cancelled,
        )?;
        Ok((momentum_report, scalar_report))
    }
}

impl ReservoirDynamicScalar {
    pub fn required_bytes(cells: usize) -> Result<usize, DynamicSgsError> {
        PeriodicDynamicScalar::required_bytes(cells)
    }

    pub fn new(
        shape: [usize; 3],
        normal_axis: usize,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        let mut inner = PeriodicDynamicScalar::new_geometry(
            shape,
            FilterGeometry::Channel(normal_axis),
            maximum_bytes,
        )?;
        inner.boundary_gradient = ScalarBoundaryGradient::Reservoir;
        Ok(Self { inner, normal_axis })
    }

    pub fn workspace_bytes(&self) -> usize {
        self.inner.workspace_bytes()
    }
    pub fn diffusivities(&self) -> &[f64] {
        self.inner.diffusivities()
    }
    pub fn coefficient(&self, cell: usize) -> Option<f64> {
        self.inner.coefficient(cell)
    }
    pub fn coefficients(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.inner.coefficients()
    }
}

#[cfg(test)]
#[path = "reservoir/tests.rs"]
mod tests;
