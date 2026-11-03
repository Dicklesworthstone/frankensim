//! One common dynamic diffusivity for same-unit components in a plane channel.
//!
//! The joint Germano objective, common component normalization and bounded
//! storage are shared with [`super::super::scalar::components`]. Raw component
//! contractions are summed before one homogeneous-plane reduction. This fits
//! all species, including the closure species, without giving trace components
//! artificial weight or allocating a component-by-cell copy.
//!
//! No composition wall datum is invented. An impermeable coupled molecular
//! flux can balance composition and temperature gradients (Soret diffusion),
//! so zero species flux does not imply a zero mass-fraction gradient. Both the
//! physical and Favre-filtered component gradients use a second-order one-sided
//! reconstruction from three fluid centres at the walls. The native transport
//! owner still imposes its actual wall flux and discrete zero-sum mass flux.
//! This does not prescribe wall composition, replace Maxwell--Stefan transport,
//! or model heterogeneous reactions or turbulent reaction-rate correlations.

use super::super::scalar::components::{DynamicComponentsConfig, PeriodicDynamicComponents};
use super::*;

/// Plane averaging is intrinsic; the native owner bounds physical diffusivity
/// against its coupled transport operator without an arbitrary Schmidt ratio.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChannelDynamicComponentsConfig {
    pub maximum_coefficient: f64,
    pub usage: SgsUse,
}

impl ChannelDynamicComponentsConfig {
    pub fn validate(self) -> Result<(), DynamicSgsError> {
        self.fit_config().validate()
    }

    fn fit_config(self) -> DynamicComponentsConfig {
        DynamicComponentsConfig {
            averaging: DynamicAveraging::Domain,
            maximum_coefficient: self.maximum_coefficient,
            usage: self.usage,
        }
    }
}

/// Eighty numerical bytes per fluid cell, independent of component count.
/// Shape, wall normal and component count are fixed for the workspace. Scratch
/// and diffusivity are consumable only after the complete preparation succeeds.
pub struct ChannelDynamicComponents {
    inner: PeriodicDynamicComponents,
    normal_axis: usize,
}

impl ChannelDynamicComponents {
    pub fn required_bytes(cells: usize) -> Result<usize, DynamicSgsError> {
        PeriodicDynamicComponents::required_bytes(cells)
    }

    /// At least three cells along every axis, including the wall normal.
    /// All components must use the same units; include every mass fraction.
    pub fn new(
        shape: [usize; 3],
        normal_axis: usize,
        component_count: usize,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        Ok(Self {
            inner: PeriodicDynamicComponents::new_geometry(
                shape,
                FilterGeometry::Channel(normal_axis),
                component_count,
                maximum_bytes,
            )?,
            normal_axis,
        })
    }

    pub fn workspace_bytes(&self) -> usize {
        self.inner.workspace_bytes()
    }

    pub fn component_count(&self) -> usize {
        self.inner.component_count()
    }

    pub fn diffusivities(&self) -> &[f64] {
        self.inner.diffusivities()
    }
}

impl ChannelDynamicSgs {
    /// Fit momentum, heat and common component diffusion from one current flow
    /// preparation. Heat uses the supplied physical wall data; component wall
    /// gradients are reconstructed from interior values at both filter scales.
    /// The component accessor must return the same immutable snapshot on every
    /// visit, including every closure or absent species.
    ///
    /// All three modeled transports are independent. No output is consumable
    /// if any fit refuses. Preparation allocates nothing and polls cancellation
    /// at most every 64 cells in each component sweep. Physical inputs are read
    /// only; the native owner admits and publishes the complete transport step.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_with_scalar_and_components(
        &mut self,
        momentum: ChannelDynamicSgsConfig,
        heat_config: ChannelDynamicScalarConfig,
        component_config: ChannelDynamicComponentsConfig,
        density: &[f64],
        velocity: &[[f64; 3]],
        wall_velocity: [[f64; 3]; 2],
        heat: &[f64],
        heat_walls: [ScalarWall; 2],
        heat_workspace: &mut ChannelDynamicScalar,
        components: impl Fn(usize, usize) -> f64,
        component_workspace: &mut ChannelDynamicComponents,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(DynamicSgsReport, DynamicScalarReport, DynamicScalarReport), DynamicSgsError> {
        poll(&mut cancelled, 0)?;
        heat_config.validate()?;
        component_config.validate()?;
        if heat_workspace.normal_axis != self.normal_axis
            || component_workspace.normal_axis != self.normal_axis
            || heat_workspace.inner.shape != self.shape()
            || component_workspace.inner.shape() != self.shape()
            || heat.len() != self.inner.velocity.len()
        {
            return Err(DynamicSgsError::Invalid(
                "channel heat, components and flow geometry differ",
            ));
        }
        let common_scale = component_workspace
            .inner
            .measure_scale(&components, &mut cancelled)?;
        let heat_scale =
            heat_workspace
                .inner
                .normalize_with_walls(heat, Some(heat_walls), &mut cancelled)?;
        let (momentum_report, scales) = self.inner.prepare_scaled_with_walls(
            momentum.fit_config(),
            density,
            velocity,
            wall_velocity,
            &mut cancelled,
        )?;
        let heat_report = heat_workspace.inner.prepare_from(
            &self.inner,
            heat_config.fit_config(),
            density,
            scales,
            heat_scale,
            &mut cancelled,
        )?;
        let component_report = component_workspace.inner.prepare_from(
            &self.inner,
            component_config.fit_config(),
            density,
            scales,
            common_scale,
            &components,
            &mut cancelled,
        )?;
        Ok((momentum_report, heat_report, component_report))
    }
}

#[cfg(test)]
#[path = "components/tests.rs"]
mod tests;
