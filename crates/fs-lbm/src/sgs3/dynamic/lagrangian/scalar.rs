//! Joint Lagrangian histories for momentum and one resolved scalar.
//!
//! The scalar flux is `-rho C_phi Delta^2 |S| grad(phi)`. Its Germano
//! numerator `L_i M_i` and denominator `M_i M_i` use the existing scalar
//! filter kernel, independently of the momentum contraction. Each pair is
//! transported with the same physical velocity by the parent pathline update;
//! the caller supplies separate, positive memory times. A fixed turbulent
//! Prandtl/Schmidt number and division by momentum viscosity are unnecessary.
//!
//! Scalar contractions scale as `rho_scale^2 u_scale^2 phi_scale^2`, whereas
//! stress contractions scale as `rho_scale^2 u_scale^4`. Their accepted
//! logarithmic scales and history banks are therefore separate. Both fits
//! share one current Favre velocity/density/strain preparation and one sample
//! time. Only a joint trial can publish history, preventing a scalar from
//! advancing alone or being paired with another momentum trajectory.
//!
//! This combines the scalar identity of Squires (1991), equations (2.13)--
//! (2.19), <https://web.stanford.edu/group/ctr/ResBriefs/1991/19_SQUIRES.pdf>,
//! with the parent first-order pathline averaging (Meneveau, Lund & Cabot,
//! 1996, doi:10.1017/S0022112096007379). It is a declared discretization with
//! authored memory times, not the papers' adaptive/calibrated timescale.
//! Periodic Cartesian domains and resolved plane channels are admitted.
//! Channel walls use the original scalar derivative policy at both filter
//! scales, and the parent's positive, wall-bounded history interpolation.
//! Backscatter is excluded from emitted diffusion while signed history is
//! retained. There is no thermal wall-function, reaction-rate, scalar-variance,
//! isotropic SGS-energy or DNS qualification.

use super::super::channel::ScalarWall;
use super::super::scalar::{DynamicScalarConfig, DynamicScalarReport, PeriodicDynamicScalar};
use super::*;

/// Joint stress, heat and common-component pathline histories.
#[path = "scalar/components.rs"]
pub mod components;

/// Explicit scalar-memory controls, independent of momentum-memory controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LagrangianScalarConfig {
    /// Positive relaxation time in the same units as the physical sample time.
    pub memory_time: f64,
    /// Maximum departure along each axis in cells, in (0, 1].
    pub maximum_departure_cells: f64,
    /// Admission bound on C_phi, not on its square or square root.
    pub maximum_coefficient: f64,
    /// The consuming transport owner enforces this against actual molecular
    /// diffusivity at its face/material state. This is never a clamp.
    pub maximum_diffusivity_to_molecular: f64,
    pub usage: SgsUse,
}

impl LagrangianScalarConfig {
    pub fn validate(self) -> Result<(), DynamicSgsError> {
        DynamicScalarConfig {
            averaging: DynamicAveraging::LocalBox,
            maximum_coefficient: self.maximum_coefficient,
            maximum_diffusivity_to_molecular: self.maximum_diffusivity_to_molecular,
            usage: self.usage,
        }
        .validate()?;
        if !self.memory_time.is_finite() || self.memory_time <= 0.0 {
            return Err(DynamicSgsError::Invalid(
                "positive finite scalar Lagrangian memory time required",
            ));
        }
        if !self.maximum_departure_cells.is_finite()
            || self.maximum_departure_cells <= 0.0
            || self.maximum_departure_cells > 1.0
        {
            return Err(DynamicSgsError::Invalid(
                "scalar Lagrangian departure limit must be in (0, 1]",
            ));
        }
        Ok(())
    }
}

struct ScalarMemory {
    workspace: PeriodicDynamicScalar,
    history: Vec<[f64; 2]>,
    candidate: Vec<[f64; 2]>,
    history_log_scale: Option<f64>,
}

/// One trajectory's shared momentum/scalar filters and accepted histories.
/// The additional numerical payload is 264 bytes/cell, excluding metadata.
/// Preparation and publication reuse all storage and allocate no field arrays.
pub struct LagrangianDynamicScalar {
    flow: LagrangianDynamicSgs,
    scalar: ScalarMemory,
}

/// Both fits have been admitted, but neither history has advanced. Dropping
/// the trial rejects both; commit only after the consuming flow step accepts.
#[must_use = "commit after a successful flow step, or drop to reject both histories"]
pub struct LagrangianScalarTrial<'a> {
    momentum: LagrangianSgsTrial<'a>,
    scalar: &'a mut ScalarMemory,
    log_scale: Option<f64>,
    report: DynamicScalarReport,
}

impl LagrangianScalarTrial<'_> {
    pub const fn momentum_report(&self) -> DynamicSgsReport {
        self.momentum.report()
    }

    pub const fn scalar_report(&self) -> DynamicScalarReport {
        self.report
    }

    pub fn viscosities(&self) -> &[f64] {
        self.momentum.viscosities()
    }

    pub fn diffusivities(&self) -> &[f64] {
        self.scalar.workspace.diffusivities()
    }

    /// Fitted C_s squared; conservative faces supply their own strain.
    pub fn coefficients_squared(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.momentum.coefficients_squared()
    }

    pub fn coefficient_squared(&self, cell: usize) -> Option<f64> {
        self.momentum.coefficient_squared(cell)
    }

    /// Fitted nonnegative C_phi, independent of centered cell strain and of
    /// C_s squared. Conservative scalar faces supply their own current strain.
    pub fn scalar_coefficients(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.scalar.workspace.basis.iter().map(|value| value[0])
    }

    pub fn scalar_coefficient(&self, cell: usize) -> Option<f64> {
        self.scalar.workspace.basis.get(cell).map(|value| value[0])
    }

    /// Infallible, allocation-free publication of both histories at one time.
    pub fn commit(self) {
        self.momentum.commit();
        core::mem::swap(&mut self.scalar.history, &mut self.scalar.candidate);
        self.scalar.history_log_scale = self.log_scale;
    }
}

impl LagrangianDynamicScalar {
    pub fn required_bytes(cells: usize) -> Result<usize, DynamicSgsError> {
        let history =
            cells
                .checked_mul(4 * core::mem::size_of::<f64>())
                .ok_or(DynamicSgsError::Invalid(
                    "Lagrangian scalar history size overflow",
                ))?;
        LagrangianDynamicSgs::required_bytes(cells)?
            .checked_add(PeriodicDynamicScalar::required_bytes(cells)?)
            .and_then(|bytes| bytes.checked_add(history))
            .ok_or(DynamicSgsError::Invalid(
                "Lagrangian scalar workspace size overflow",
            ))
    }

    pub fn new(
        shape: [usize; 3],
        spacing: [f64; 3],
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        Self::new_geometry(shape, spacing, FilterGeometry::Periodic, maximum_bytes)
    }

    /// Joint history on a resolved plane channel, with the same byte budget
    /// as the periodic owner. Supply both physical wall laws at each sample.
    pub fn new_channel(
        shape: [usize; 3],
        spacing: [f64; 3],
        normal_axis: usize,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        Self::new_geometry(
            shape,
            spacing,
            FilterGeometry::Channel(normal_axis),
            maximum_bytes,
        )
    }

    fn new_geometry(
        shape: [usize; 3],
        spacing: [f64; 3],
        geometry: FilterGeometry,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        geometry.validate_shape(shape)?;
        let cells = shape
            .into_iter()
            .try_fold(1usize, |n, d| n.checked_mul(d))
            .ok_or(DynamicSgsError::Invalid(
                "Lagrangian scalar cell count overflow",
            ))?;
        let required = Self::required_bytes(cells)?;
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            });
        }
        let result = Self {
            flow: LagrangianDynamicSgs::new_geometry(shape, spacing, geometry, maximum_bytes)?,
            scalar: ScalarMemory {
                workspace: PeriodicDynamicScalar::new_geometry(shape, geometry, maximum_bytes)?,
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

    /// Actual Vec capacities, including both filters and all four history banks.
    pub fn workspace_bytes(&self) -> usize {
        self.flow.workspace_bytes()
            + self.scalar.workspace.workspace_bytes()
            + (self.scalar.history.capacity() + self.scalar.candidate.capacity())
                * core::mem::size_of::<[f64; 2]>()
    }

    pub const fn shape(&self) -> [usize; 3] {
        self.flow.shape()
    }

    pub const fn filter_width(&self) -> f64 {
        self.flow.filter_width()
    }

    pub const fn normal_axis(&self) -> Option<usize> {
        self.flow.normal_axis()
    }

    /// Common last committed sample time, not the consuming flow's endpoint.
    pub const fn accepted_time(&self) -> Option<f64> {
        self.flow.accepted_time()
    }

    /// Prepare both fits from the same physical fields at one increasing time.
    /// The first sample uses each original local-box fit. Subsequent samples
    /// transport and relax each signed contraction pair independently.
    ///
    /// No output is exposed if either fit refuses. All errors, cancellation and
    /// dropped trials leave both accepted histories unchanged; only reusable
    /// private scratch may change. A new trajectory needs a new owner.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare(
        &mut self,
        momentum: LagrangianSgsConfig,
        config: LagrangianScalarConfig,
        time: f64,
        density: &[f64],
        velocity: &[[f64; 3]],
        scalar: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<LagrangianScalarTrial<'_>, DynamicSgsError> {
        self.prepare_scaled(momentum, config, time, density, velocity, scalar, cancelled)
            .map(|(trial, _)| trial)
    }

    /// Prepare both channel fits with the current lower/upper physical walls.
    /// Fixed scalar values include their wall-to-fluid variation in the common
    /// normalization. ZeroGradient and Interior policies require three normal
    /// fluid cells even if the current scalar field is uniform.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_with_walls(
        &mut self,
        momentum: LagrangianSgsConfig,
        config: LagrangianScalarConfig,
        time: f64,
        density: &[f64],
        velocity: &[[f64; 3]],
        wall_velocity: [[f64; 3]; 2],
        scalar: &[f64],
        scalar_walls: [ScalarWall; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<LagrangianScalarTrial<'_>, DynamicSgsError> {
        self.prepare_scaled_with_walls(
            momentum,
            config,
            time,
            density,
            velocity,
            Some(wall_velocity),
            scalar,
            Some(scalar_walls),
            cancelled,
        )
        .map(|(trial, _)| trial)
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare_scaled(
        &mut self,
        momentum: LagrangianSgsConfig,
        config: LagrangianScalarConfig,
        time: f64,
        density: &[f64],
        velocity: &[[f64; 3]],
        scalar: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<(LagrangianScalarTrial<'_>, InputScales), DynamicSgsError> {
        self.prepare_scaled_with_walls(
            momentum, config, time, density, velocity, None, scalar, None, cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare_scaled_with_walls(
        &mut self,
        momentum: LagrangianSgsConfig,
        config: LagrangianScalarConfig,
        time: f64,
        density: &[f64],
        velocity: &[[f64; 3]],
        wall_velocity: Option<[[f64; 3]; 2]>,
        scalar: &[f64],
        scalar_walls: Option<[ScalarWall; 2]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(LagrangianScalarTrial<'_>, InputScales), DynamicSgsError> {
        poll(&mut cancelled, 0)?;
        config.validate()?;
        if scalar.len() != self.scalar.history.len() {
            return Err(DynamicSgsError::Invalid(
                "scalar field and Lagrangian workspace lengths differ",
            ));
        }
        match (self.normal_axis(), scalar_walls) {
            (Some(normal), Some(walls)) => {
                if self.shape()[normal] < 3
                    && walls
                        .iter()
                        .any(|wall| matches!(wall, ScalarWall::ZeroGradient | ScalarWall::Interior))
                {
                    return Err(DynamicSgsError::Invalid(
                        "scalar wall reconstruction requires three normal fluid cells",
                    ));
                }
            }
            (None, None) => {}
            _ => {
                return Err(DynamicSgsError::Invalid(
                    "Lagrangian scalar wall data must match workspace geometry",
                ));
            }
        }
        let scalar_scale =
            self.scalar
                .workspace
                .normalize_with_walls(scalar, scalar_walls, &mut cancelled)?;
        let previous = self.flow.accepted_time();
        let elapsed = previous.map_or(0.0, |t| time - t);
        let (momentum, scales) = self.flow.prepare_scaled_with_walls(
            momentum,
            time,
            density,
            velocity,
            wall_velocity,
            &mut cancelled,
        )?;
        let flow = &momentum.owner.flow;
        if scales.velocity == 0.0 || scalar_scale == 0.0 {
            // The current source is zero, but old pathline history must still
            // advect and decay. Do not take the instantaneous scalar zero exit.
            for (tile, pairs) in self.scalar.workspace.contraction.chunks_mut(64).enumerate() {
                poll(&mut cancelled, tile * 64)?;
                pairs.fill([0.0; 2]);
            }
        } else {
            self.scalar
                .workspace
                .prepare_contractions(flow, density, scales, &mut cancelled)?;
        }
        let current_log = (scales.velocity > 0.0 && scalar_scale > 0.0).then(|| {
            2.0 * det::ln(scales.density)
                + 2.0 * det::ln(scales.velocity)
                + 2.0 * det::ln(scalar_scale)
        });
        let update = HistoryUpdate::new(
            flow,
            previous.is_some(),
            self.scalar.history_log_scale,
            current_log,
            elapsed,
            config.memory_time,
            config.maximum_departure_cells,
        )?;
        let (report, log_scale) = prepare_history(
            &mut self.scalar.workspace,
            &self.scalar.history,
            &mut self.scalar.candidate,
            &update,
            flow,
            scales,
            velocity,
            config.maximum_coefficient,
            &mut cancelled,
        )?;
        Ok((
            LagrangianScalarTrial {
                momentum,
                scalar: &mut self.scalar,
                log_scale,
                report,
            },
            scales,
        ))
    }
}

// Heat and a same-unit component bank share emission and history rescaling,
// while keeping their sources, physical scales and accepted banks separate.
#[allow(clippy::too_many_arguments)]
fn prepare_history(
    workspace: &mut PeriodicDynamicScalar,
    history: &[[f64; 2]],
    candidate: &mut [[f64; 2]],
    update: &HistoryUpdate,
    flow: &PeriodicDynamicSgs,
    scales: InputScales,
    velocity: &[[f64; 3]],
    maximum_coefficient: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(DynamicScalarReport, Option<f64>), DynamicSgsError> {
    let mut report = DynamicScalarReport {
        minimum_raw_coefficient: f64::INFINITY,
        ..DynamicScalarReport::default()
    };
    let mut amplitude: f64 = 0.0;
    for cell in 0..history.len() {
        poll(cancelled, cell)?;
        let pair = update.pair(cell, velocity[cell], history, &workspace.contraction)?;
        let raw = if pair[1] == 0.0 {
            report.undetermined_cells += 1;
            0.0
        } else {
            pair[0] / pair[1]
        };
        if !raw.is_finite() {
            return Err(DynamicSgsError::Cell {
                cell,
                reason: "unrepresentable Lagrangian scalar coefficient",
            });
        }
        report.minimum_raw_coefficient = report.minimum_raw_coefficient.min(raw);
        report.backscatter_limited_cells += usize::from(raw < 0.0);
        let coefficient = raw.max(0.0);
        if coefficient > maximum_coefficient {
            return Err(DynamicSgsError::Coefficient {
                cell,
                observed: coefficient,
                maximum: maximum_coefficient,
            });
        }
        let diffusivity =
            coefficient * flow.strain_rate[cell] * scales.velocity * flow.filter_width;
        if !diffusivity.is_finite() {
            return Err(DynamicSgsError::Cell {
                cell,
                reason: "unrepresentable Lagrangian scalar diffusivity",
            });
        }
        workspace.diffusivity[cell] = diffusivity;
        // The vector basis has been fully consumed by contractions. Reuse
        // its first lane for C_phi, even when centered strain vanishes.
        workspace.basis[cell][0] = coefficient;
        candidate[cell] = pair;
        amplitude = amplitude.max(pair[0].abs()).max(pair[1]);
        report.maximum_coefficient = report.maximum_coefficient.max(coefficient);
        report.maximum_diffusivity = report.maximum_diffusivity.max(diffusivity);
    }
    let log_scale = update.normalize(candidate, amplitude, cancelled)?;
    poll(cancelled, 0)?;
    Ok((report, log_scale))
}

#[cfg(test)]
#[path = "scalar/tests.rs"]
mod tests;
