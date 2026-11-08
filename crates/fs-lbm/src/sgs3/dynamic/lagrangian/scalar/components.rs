//! Independent pathline memory for stress, heat and a common component flux.
//!
//! The existing component kernel sums every same-unit scalar's Germano
//! contractions before fitting one coefficient. One common normalization
//! preserves trace-component weights. Its physical history scale is
//! rho_scale^2 u_scale^2 component_scale^2, independently of the heat scalar
//! and the momentum scale. No species-by-cell field is allocated.
//!
//! All three fits share the current Favre flow and sample time. A trial
//! publishes all three histories together only after its consumer accepts.
//! Zero current sources still advect and decay the accepted signed history.
//! The periodic or resolved-channel geometry, authored timescales and
//! exploratory limits of the parent scalar and momentum histories apply.
//! Channel component gradients use three interior fluid centers at the walls;
//! impermeability never invents a composition wall value or zero gradient.

use super::super::super::scalar::components::PeriodicDynamicComponents;
use super::*;

/// Memory controls for a common same-unit component diffusivity.
/// The native consumer separately bounds physical species diffusion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LagrangianComponentsConfig {
    pub memory_time: f64,
    pub maximum_departure_cells: f64,
    /// Bound on C_phi, not on its square or square root.
    pub maximum_coefficient: f64,
    pub usage: SgsUse,
}

impl LagrangianComponentsConfig {
    pub fn validate(self) -> Result<(), DynamicSgsError> {
        LagrangianScalarConfig {
            memory_time: self.memory_time,
            maximum_departure_cells: self.maximum_departure_cells,
            maximum_coefficient: self.maximum_coefficient,
            // A coupled molecular species operator has no unique scalar ratio.
            maximum_diffusivity_to_molecular: 0.0,
            usage: self.usage,
        }
        .validate()
    }
}

struct ComponentMemory {
    workspace: PeriodicDynamicComponents,
    history: Vec<[f64; 2]>,
    candidate: Vec<[f64; 2]>,
    history_log_scale: Option<f64>,
}

/// One trajectory's shared momentum, heat and common-component histories.
/// Numerical payload is 376 bytes/cell, independent of component count.
/// Preparation and commit reuse the bounded allocations.
pub struct LagrangianDynamicComponents {
    flow: LagrangianDynamicScalar,
    components: ComponentMemory,
}

/// All three fits are admitted, but no accepted history has advanced.
#[must_use = "commit after a successful flow step, or drop to reject all three histories"]
pub struct LagrangianComponentsTrial<'a> {
    joint: LagrangianScalarTrial<'a>,
    components: &'a mut ComponentMemory,
    log_scale: Option<f64>,
    report: DynamicScalarReport,
}

impl LagrangianComponentsTrial<'_> {
    pub const fn momentum_report(&self) -> DynamicSgsReport {
        self.joint.momentum_report()
    }

    pub const fn momentum_memory_time_range(&self) -> [f64; 2] {
        self.joint.momentum_memory_time_range()
    }

    pub const fn scalar_report(&self) -> DynamicScalarReport {
        self.joint.scalar_report()
    }

    pub const fn components_report(&self) -> DynamicScalarReport {
        self.report
    }

    pub fn viscosities(&self) -> &[f64] {
        self.joint.viscosities()
    }

    pub fn scalar_diffusivities(&self) -> &[f64] {
        self.joint.diffusivities()
    }

    pub fn component_diffusivities(&self) -> &[f64] {
        self.components.workspace.diffusivities()
    }

    pub fn coefficients_squared(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.joint.coefficients_squared()
    }

    pub fn coefficient_squared(&self, cell: usize) -> Option<f64> {
        self.joint.coefficient_squared(cell)
    }

    pub fn scalar_coefficients(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.joint.scalar_coefficients()
    }

    pub fn scalar_coefficient(&self, cell: usize) -> Option<f64> {
        self.joint.scalar_coefficient(cell)
    }

    /// Common C_phi independent of centered strain. Conservative faces
    /// combine it with their own current strain and the squared filter width.
    pub fn component_coefficients(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.components.workspace.coefficients()
    }

    pub fn component_coefficient(&self, cell: usize) -> Option<f64> {
        self.components.workspace.coefficient(cell)
    }

    /// Infallible, allocation-free publication of all histories at one time.
    pub fn commit(self) {
        self.joint.commit();
        core::mem::swap(&mut self.components.history, &mut self.components.candidate);
        self.components.history_log_scale = self.log_scale;
    }
}

impl LagrangianDynamicComponents {
    /// Select adaptive momentum memory before the first joint fit. Heat and
    /// component relaxation keep their independently authored memory times.
    pub fn with_adaptive_momentum_memory(
        mut self,
        config: super::super::AdaptiveLagrangianMemory,
    ) -> Result<Self, DynamicSgsError> {
        self.flow = self.flow.with_adaptive_momentum_memory(config)?;
        Ok(self)
    }

    pub const fn adaptive_momentum_memory(&self) -> Option<super::super::AdaptiveLagrangianMemory> {
        self.flow.adaptive_momentum_memory()
    }

    pub const fn accepted_momentum_memory_time_range(&self) -> Option<[f64; 2]> {
        self.flow.accepted_momentum_memory_time_range()
    }

    pub fn required_bytes(cells: usize) -> Result<usize, DynamicSgsError> {
        let history =
            cells
                .checked_mul(4 * core::mem::size_of::<f64>())
                .ok_or(DynamicSgsError::Invalid(
                    "Lagrangian component history size overflow",
                ))?;
        LagrangianDynamicScalar::required_bytes(cells)?
            .checked_add(PeriodicDynamicComponents::required_bytes(cells)?)
            .and_then(|bytes| bytes.checked_add(history))
            .ok_or(DynamicSgsError::Invalid(
                "Lagrangian component workspace size overflow",
            ))
    }

    pub fn new(
        shape: [usize; 3],
        spacing: [f64; 3],
        component_count: usize,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        Self::new_geometry(
            shape,
            spacing,
            FilterGeometry::Periodic,
            component_count,
            maximum_bytes,
        )
    }

    /// Joint momentum, heat and component history in a resolved plane channel.
    /// Each dimension needs at least three fluid cells. Numerical storage is
    /// unchanged and remains independent of the number of components.
    pub fn new_channel(
        shape: [usize; 3],
        spacing: [f64; 3],
        normal_axis: usize,
        component_count: usize,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        Self::new_geometry(
            shape,
            spacing,
            FilterGeometry::Channel(normal_axis),
            component_count,
            maximum_bytes,
        )
    }

    fn new_geometry(
        shape: [usize; 3],
        spacing: [f64; 3],
        geometry: FilterGeometry,
        component_count: usize,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        geometry.validate_shape(shape)?;
        if let FilterGeometry::Channel(normal) = geometry {
            if shape[normal] < 3 {
                return Err(DynamicSgsError::Invalid(
                    "component wall derivatives require three normal fluid cells",
                ));
            }
        }
        if component_count == 0 {
            return Err(DynamicSgsError::Invalid("at least one component required"));
        }
        let cells = shape
            .into_iter()
            .try_fold(1usize, |n, d| n.checked_mul(d))
            .ok_or(DynamicSgsError::Invalid(
                "Lagrangian component cell count overflow",
            ))?;
        cells
            .checked_mul(component_count)
            .ok_or(DynamicSgsError::Invalid(
                "Lagrangian component work count overflow",
            ))?;
        let required = Self::required_bytes(cells)?;
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            });
        }
        let result = Self {
            flow: LagrangianDynamicScalar::new_geometry(shape, spacing, geometry, maximum_bytes)?,
            components: ComponentMemory {
                workspace: PeriodicDynamicComponents::new_geometry(
                    shape,
                    geometry,
                    component_count,
                    maximum_bytes,
                )?,
                history: allocated(cells, [0.0; 2])?,
                candidate: allocated(cells, [0.0; 2])?,
                history_log_scale: None,
            },
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
        self.flow.workspace_bytes()
            + self.components.workspace.workspace_bytes()
            + (self.components.history.capacity() + self.components.candidate.capacity())
                * core::mem::size_of::<[f64; 2]>()
    }

    pub const fn shape(&self) -> [usize; 3] {
        self.flow.shape()
    }

    pub const fn filter_width(&self) -> f64 {
        self.flow.filter_width()
    }

    pub fn component_count(&self) -> usize {
        self.components.workspace.component_count()
    }

    pub const fn normal_axis(&self) -> Option<usize> {
        self.flow.normal_axis()
    }

    pub const fn accepted_time(&self) -> Option<f64> {
        self.flow.accepted_time()
    }

    /// Prepare three independent fits from the same immutable physical snapshot.
    /// Include every component, including any reconstructed closure species.
    /// The accessor must return the same value on every visit; work is linear
    /// in cells times components with bounded cancellation checks.
    /// Errors, cancellation and dropped trials preserve every accepted bank,
    /// logarithmic scale and sample time. Only private scratch may change.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare(
        &mut self,
        momentum: LagrangianSgsConfig,
        heat: LagrangianScalarConfig,
        config: LagrangianComponentsConfig,
        time: f64,
        density: &[f64],
        velocity: &[[f64; 3]],
        scalar: &[f64],
        components: impl Fn(usize, usize) -> f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<LagrangianComponentsTrial<'_>, DynamicSgsError> {
        self.prepare_impl(
            momentum, heat, config, time, density, velocity, None, scalar, None, components,
            cancelled,
        )
    }

    /// Prepare all three channel fits from one physical sample. Wall velocities
    /// and heat laws are lower then upper; composition wall derivatives are
    /// reconstructed from the original fluid values without a prescribed datum.
    /// A refusal or dropped trial leaves every accepted bank and clock unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_with_walls(
        &mut self,
        momentum: LagrangianSgsConfig,
        heat: LagrangianScalarConfig,
        config: LagrangianComponentsConfig,
        time: f64,
        density: &[f64],
        velocity: &[[f64; 3]],
        wall_velocity: [[f64; 3]; 2],
        scalar: &[f64],
        scalar_walls: [ScalarWall; 2],
        components: impl Fn(usize, usize) -> f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<LagrangianComponentsTrial<'_>, DynamicSgsError> {
        self.prepare_impl(
            momentum,
            heat,
            config,
            time,
            density,
            velocity,
            Some(wall_velocity),
            scalar,
            Some(scalar_walls),
            components,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare_impl(
        &mut self,
        momentum: LagrangianSgsConfig,
        heat: LagrangianScalarConfig,
        config: LagrangianComponentsConfig,
        time: f64,
        density: &[f64],
        velocity: &[[f64; 3]],
        wall_velocity: Option<[[f64; 3]; 2]>,
        scalar: &[f64],
        scalar_walls: Option<[ScalarWall; 2]>,
        components: impl Fn(usize, usize) -> f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LagrangianComponentsTrial<'_>, DynamicSgsError> {
        poll(&mut cancelled, 0)?;
        config.validate()?;
        let common_scale = self
            .components
            .workspace
            .measure_scale(&components, &mut cancelled)?;
        let previous = self.flow.accepted_time();
        let elapsed = previous.map_or(0.0, |accepted| time - accepted);
        let (joint, scales) = self.flow.prepare_scaled_with_walls(
            momentum,
            heat,
            time,
            density,
            velocity,
            wall_velocity,
            scalar,
            scalar_walls,
            &mut cancelled,
        )?;
        let flow = &joint.momentum.owner.flow;
        self.components.workspace.prepare_contractions(
            flow,
            density,
            scales,
            common_scale,
            &components,
            &mut cancelled,
        )?;
        let current_log = (scales.velocity > 0.0 && common_scale > 0.0).then(|| {
            2.0 * det::ln(scales.density)
                + 2.0 * det::ln(scales.velocity)
                + 2.0 * det::ln(common_scale)
        });
        let update = HistoryUpdate::new(
            flow,
            previous.is_some(),
            self.components.history_log_scale,
            current_log,
            elapsed,
            config.memory_time,
            config.maximum_departure_cells,
        )?;
        let (report, log_scale) = prepare_history(
            &mut self.components.workspace.scalar,
            &self.components.history,
            &mut self.components.candidate,
            &update,
            flow,
            scales,
            velocity,
            config.maximum_coefficient,
            &mut cancelled,
        )?;
        Ok(LagrangianComponentsTrial {
            joint,
            components: &mut self.components,
            log_scale,
            report,
        })
    }
}

#[cfg(test)]
#[path = "components/tests.rs"]
mod tests;

#[cfg(test)]
#[path = "components/channel_tests.rs"]
mod channel_tests;
