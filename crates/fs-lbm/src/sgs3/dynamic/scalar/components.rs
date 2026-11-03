//! A constrained dynamic fit for several scalars in the same units.
//!
//! For mass fractions `Y_s`, one common diffusivity preserves the continuum
//! constraint `sum_s j_s = -rho D_t grad(sum_s Y_s) = 0`. The original face
//! transport owner remains responsible for enforcing that constraint in its
//! discrete flux, preserving absent species, and carrying partial enthalpies.
//! This module fits the common coefficient by minimizing
//! `sum_s <|L_s - C M_s|^2>`, hence
//! `C = sum_s <L_s.M_s> / sum_s <M_s.M_s>`.
//!
//! This is the explicit equal-weight vector extension of the scalar Germano
//! identity documented in the parent module. Every component uses the SAME
//! normalization; independently normalized component fits would give trace
//! species spurious weight. Include every species, including any species
//! reconstructed from the mass constraint. The exact-arithmetic fit is
//! invariant to component permutation and orthogonal changes of same-unit
//! component basis. Floating reductions retain the caller's component order.
//! Molecular differential diffusion and reaction-rate closure are not modeled.

use super::*;

/// Controls for one common same-unit scalar diffusivity. The native consumer
/// separately bounds the physical diffusivity using its transport operator;
/// a coupled multispecies molecular operator has no unique scalar denominator.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DynamicComponentsConfig {
    /// Spatial averaging of the summed contractions, before division.
    pub averaging: DynamicAveraging,
    /// Admission bound on the fitted coefficient, never a positive-value clamp.
    pub maximum_coefficient: f64,
    /// Dynamic transport remains an explicit exploratory choice.
    pub usage: SgsUse,
}

impl DynamicComponentsConfig {
    /// Validate the fit policy before any field work.
    pub fn validate(self) -> Result<(), DynamicSgsError> {
        self.fit_config().validate()
    }

    fn fit_config(self) -> DynamicScalarConfig {
        DynamicScalarConfig {
            averaging: self.averaging,
            maximum_coefficient: self.maximum_coefficient,
            // The private fit kernel does not consume this native-owner bound.
            maximum_diffusivity_to_molecular: 0.0,
            usage: self.usage,
        }
    }
}

/// One scalar workspace and two joint contractions per cell: 80 numerical
/// bytes/cell, independent of component count. Preparation allocates nothing.
/// Output and scratch may only be consumed after the entire combined fit
/// succeeds. A refusal may overwrite scratch but never modifies input fields.
pub struct PeriodicDynamicComponents {
    scalar: PeriodicDynamicScalar,
    contractions: Vec<[f64; 2]>,
    components: usize,
}

impl PeriodicDynamicComponents {
    /// Numerical payload, excluding fixed owner and allocator metadata.
    pub fn required_bytes(cells: usize) -> Result<usize, DynamicSgsError> {
        cells
            .checked_mul(10 * core::mem::size_of::<f64>())
            .ok_or(DynamicSgsError::Invalid(
                "component workspace size overflow",
            ))
    }

    /// Bind the periodic shape and exact component count. Work is linear in
    /// their product; cancellation is checked at most every 64 cells of each
    /// component sweep. No transposed component-by-cell allocation is needed.
    pub fn new(
        shape: [usize; 3],
        component_count: usize,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        if shape.iter().any(|&n| n < 3) || component_count == 0 {
            return Err(DynamicSgsError::Invalid(
                "at least three cells per periodic axis and one component required",
            ));
        }
        Self::new_geometry(
            shape,
            FilterGeometry::Periodic,
            component_count,
            maximum_bytes,
        )
    }

    pub(in super::super) fn new_geometry(
        shape: [usize; 3],
        geometry: FilterGeometry,
        component_count: usize,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        geometry.validate_shape(shape)?;
        if component_count == 0 {
            return Err(DynamicSgsError::Invalid("at least one component required"));
        }
        if let FilterGeometry::Channel(normal) = geometry {
            if shape[normal] < 3 {
                return Err(DynamicSgsError::Invalid(
                    "component wall derivatives require three normal fluid cells",
                ));
            }
        }
        let cells = shape
            .into_iter()
            .try_fold(1usize, |n, d| n.checked_mul(d))
            .ok_or(DynamicSgsError::Invalid("component cell count overflow"))?;
        cells
            .checked_mul(component_count)
            .ok_or(DynamicSgsError::Invalid("component work count overflow"))?;
        let required = Self::required_bytes(cells)?;
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            });
        }
        let mut scalar = PeriodicDynamicScalar::new_geometry(shape, geometry, maximum_bytes)?;
        if matches!(geometry, FilterGeometry::Channel(_)) {
            scalar.boundary_gradient = ScalarBoundaryGradient::Interior;
        }
        let result = Self {
            scalar,
            contractions: allocated(cells, [0.0; 2])?,
            components: component_count,
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

    /// Retained numerical allocation, including vector spare capacity.
    pub fn workspace_bytes(&self) -> usize {
        self.scalar.workspace_bytes()
            + self.contractions.capacity() * core::mem::size_of::<[f64; 2]>()
    }

    /// All components included in the joint least-squares objective.
    pub fn component_count(&self) -> usize {
        self.components
    }

    /// Common modeled diffusivity, with the flow's length-squared/time units.
    pub fn diffusivities(&self) -> &[f64] {
        self.scalar.diffusivities()
    }

    pub(in super::super) fn shape(&self) -> [usize; 3] {
        self.scalar.shape
    }

    pub(in super::super) fn measure_scale(
        &self,
        values: &impl Fn(usize, usize) -> f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, DynamicSgsError> {
        let mut scale: f64 = 0.0;
        for component in 0..self.components {
            let reference = values(component, 0);
            for cell in 0..self.contractions.len() {
                poll(cancelled, cell)?;
                let difference = values(component, cell) - reference;
                if !difference.is_finite() {
                    return Err(DynamicSgsError::Cell {
                        cell,
                        reason: "finite representable component differences required",
                    });
                }
                scale = scale.max(difference.abs());
            }
        }
        Ok(scale)
    }

    #[allow(clippy::too_many_arguments)]
    pub(in super::super) fn prepare_from(
        &mut self,
        flow: &PeriodicDynamicSgs,
        config: DynamicComponentsConfig,
        density: &[f64],
        scales: InputScales,
        common_scale: f64,
        values: &impl Fn(usize, usize) -> f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DynamicScalarReport, DynamicSgsError> {
        if common_scale == 0.0 || scales.velocity == 0.0 {
            return self.scalar.prepare_from(
                flow,
                config.fit_config(),
                density,
                scales,
                0.0,
                cancelled,
            );
        }
        for (tile, pairs) in self.contractions.chunks_mut(64).enumerate() {
            poll(cancelled, tile * 64)?;
            pairs.fill([0.0; 2]);
        }
        for component in 0..self.components {
            let reference = values(component, 0);
            let mut varies = false;
            for (cell, scalar) in self.scalar.scalar.iter_mut().enumerate() {
                poll(cancelled, cell)?;
                let difference = values(component, cell) - reference;
                *scalar = difference / common_scale;
                if !scalar.is_finite() || (*scalar == 0.0 && difference != 0.0) {
                    return Err(DynamicSgsError::Cell {
                        cell,
                        reason: "unrepresentable common-normalized component",
                    });
                }
                varies |= difference != 0.0;
            }
            // A constant or absent component contributes exactly zero.
            if !varies {
                continue;
            }
            self.scalar
                .prepare_contractions(flow, density, scales, cancelled)?;
            for (cell, (total, pair)) in self
                .contractions
                .iter_mut()
                .zip(&self.scalar.contraction)
                .enumerate()
            {
                poll(cancelled, cell)?;
                for i in 0..2 {
                    total[i] += pair[i];
                }
                if total.iter().any(|v| !v.is_finite()) {
                    return Err(DynamicSgsError::Cell {
                        cell,
                        reason: "unrepresentable joint component contraction",
                    });
                }
            }
        }
        let mut numerator = Sum::default();
        let mut denominator = Sum::default();
        for (cell, (target, pair)) in self
            .scalar
            .contraction
            .iter_mut()
            .zip(&self.contractions)
            .enumerate()
        {
            poll(cancelled, cell)?;
            *target = *pair;
            numerator.add(pair[0]);
            denominator.add(pair[1]);
        }
        // Sum raw component contractions first, then average each plane once.
        // Averaging every component separately adds sweeps and can discard a
        // small component's weight before the common objective is formed.
        flow.geometry
            .average_planes(self.scalar.shape, &mut self.scalar.contraction, cancelled)?;
        self.scalar.fit_from(
            flow,
            config.fit_config(),
            scales,
            [numerator.value, denominator.value],
            cancelled,
        )
    }
}

impl PeriodicDynamicSgs {
    /// Fit momentum, one independent heat scalar, and a common same-unit
    /// component diffusivity from one preparation of the current flow.
    /// `components(component, cell)` must return the same immutable current
    /// scalar snapshot on every visit; an AoS mass-fraction accessor avoids
    /// materializing a second transposed composition field. Every component,
    /// including a closure species, must be represented in the workspace.
    ///
    /// Heat and component diffusivity are independent of each other and of
    /// momentum viscosity. Neither output is consumable if any fit refuses.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_with_scalar_and_components(
        &mut self,
        momentum: DynamicSgsConfig,
        heat_config: DynamicScalarConfig,
        component_config: DynamicComponentsConfig,
        density: &[f64],
        velocity: &[[f64; 3]],
        heat: &[f64],
        heat_workspace: &mut PeriodicDynamicScalar,
        components: impl Fn(usize, usize) -> f64,
        component_workspace: &mut PeriodicDynamicComponents,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(DynamicSgsReport, DynamicScalarReport, DynamicScalarReport), DynamicSgsError> {
        poll(&mut cancelled, 0)?;
        heat_config.validate()?;
        component_config.validate()?;
        if heat_workspace.shape != self.shape
            || component_workspace.scalar.shape != self.shape
            || heat.len() != self.velocity.len()
        {
            return Err(DynamicSgsError::Invalid(
                "heat, component and flow workspace shapes differ",
            ));
        }
        let common_scale = component_workspace.measure_scale(&components, &mut cancelled)?;
        let heat_scale = heat_workspace.normalize(heat, &mut cancelled)?;
        let (momentum_report, scales) =
            self.prepare_scaled(momentum, density, velocity, &mut cancelled)?;
        let heat_report = heat_workspace.prepare_from(
            self,
            heat_config,
            density,
            scales,
            heat_scale,
            &mut cancelled,
        )?;
        let component_report = component_workspace.prepare_from(
            self,
            component_config,
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
