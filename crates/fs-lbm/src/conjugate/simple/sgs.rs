//! Explicit SGS on the unsteady staggered finite-volume owner.
//!
//! This is a resolved-velocity closure, separate from equilibrium LVEL.
//! Momentum consumes the existing complete symmetric effective-viscosity
//! stress. Heat uses gradient diffusion with Pr_t=0.9; neither dissipation
//! heating nor a turbulent-DNS or thermal-wall-resolution claim is inferred.

use super::{FvBoundary, Solver, Turbulence};
use crate::conjugate::{ChtError, poll};
use crate::sgs3::{SgsConfig, SgsModel, SgsUse, eddy_viscosity, vreman_eddy_viscosity};
use fs_exec::CancelGate;

/// Explicit exploration controls for transient finite-volume flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FvSgsConfig {
    /// Shared Smagorinsky, WALE or Vreman model and its authored coefficient.
    pub model: SgsModel,
    /// Positive dimensionless filter factor. Smagorinsky/WALE multiply the
    /// cube root of volume; Vreman multiplies EACH directional cell width.
    pub filter_width_factor: f64,
    /// Refuse, rather than clip, a cell exceeding this nu_t/nu ratio.
    pub maximum_eddy_to_molecular: f64,
    /// Acceptance renders refuse; no existing transient enables this model.
    pub usage: SgsUse,
}

impl FvSgsConfig {
    /// Admit coefficients and use policy before allocating or advancing flow.
    ///
    /// # Errors
    /// Invalid coefficients, filter factor or viscosity cap, or acceptance use.
    pub fn validate(self) -> Result<(), ChtError> {
        SgsConfig {
            model: self.model,
            filter_width_lu: self.filter_width_factor,
            maximum_eddy_to_molecular: self.maximum_eddy_to_molecular,
            usage: self.usage,
        }
        .validate()
        .map_err(|error| refuse(error.to_string()))
    }

    fn disabled(self) -> bool {
        self.model.coefficient() == 0.0
    }
}

pub(super) struct State {
    config: FvSgsConfig,
    velocity: Vec<[f64; 3]>,
    /// Scalar-filter cache for Smagorinsky/WALE. Vreman borrows directional
    /// geometry instead, avoiding a redundant three-vector per cell.
    filter_width: Vec<f64>,
}

fn refuse(reason: impl Into<String>) -> ChtError {
    ChtError::InvalidInput {
        field: "simple.sgs",
        reason: reason.into(),
    }
}

/// A sample is either a prescribed value or a prescribed zero derivative
/// at its actual distance from the cell centre.
#[derive(Clone, Copy)]
enum Sample {
    Value(f64),
    ZeroGradient,
}

/// Derivative of the quadratic through the centre and two boundary/sample
/// conditions. This is exact for quadratics on arbitrary nonuniform spacing;
/// zero derivative applies at the boundary face, not at the adjacent centre.
fn derivative(center: f64, minus: (f64, Sample), plus: (f64, Sample)) -> f64 {
    let (hm, hp) = (minus.0, plus.0);
    match (minus.1, plus.1) {
        (Sample::Value(low), Sample::Value(high)) => {
            ((hm / hp) * (high - center) + (hp / hm) * (center - low)) / (hm + hp)
        }
        (Sample::ZeroGradient, Sample::Value(high)) => {
            ((high - center) / hp) * (2.0 * hm / (2.0 * hm + hp))
        }
        (Sample::Value(low), Sample::ZeroGradient) => {
            ((center - low) / hm) * (2.0 * hp / (2.0 * hp + hm))
        }
        (Sample::ZeroGradient, Sample::ZeroGradient) => 0.0,
    }
}

impl Solver<'_> {
    /// Borrow the accepted tangent closure, without LVEL's legacy thermal
    /// secant substitution. Empty is the exact laminar representation. This
    /// is read-only: a thermal coupler must not recompute or accumulate nu_t.
    pub(in crate::conjugate) fn tangent_eddy_viscosity(&self) -> &[f64] {
        &self.eddy_viscosity
    }

    pub(in crate::conjugate) fn enable_sgs(
        &mut self,
        config: FvSgsConfig,
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        poll(gate)?;
        config.validate()?;
        if self.time.is_none() {
            return Err(refuse(
                "SGS requires unsteady flow, not a steady turbulence closure",
            ));
        }
        if self.config.turbulence != Turbulence::Laminar {
            return Err(refuse("SGS and equilibrium LVEL cannot be combined"));
        }
        if !self.molecular_viscosity.is_empty() {
            return Err(refuse(
                "this transient SGS entry requires the declared constant molecular viscosity",
            ));
        }
        // Preserve ordinary laminar arithmetic, including its existing stress
        // shortcut, for the exact zero-coefficient comparison/control.
        if config.disabled() {
            self.sgs = None;
            self.eddy_viscosity.clear();
            return Ok(());
        }
        let directional = matches!(config.model, SgsModel::Vreman { .. });
        let mut filter_width = Vec::with_capacity(if directional { 0 } else { self.domain.cell_count() });
        for cell in 0..self.domain.cell_count() {
            if cell % 64 == 0 {
                poll(gate)?;
            }
            if directional {
                for width in self.domain.widths(cell) {
                    let width = config.filter_width_factor * width;
                    if !(width.is_finite() && width > 0.0) {
                        return Err(refuse("directional filter width is not representable on this cell geometry"));
                    }
                }
            } else {
                let width = config.filter_width_factor
                    * fs_math::det::pow(self.domain.volume(cell), 1.0 / 3.0);
                if !(width.is_finite() && width > 0.0) {
                    return Err(refuse("filter width is not representable on this cell geometry"));
                }
                filter_width.push(width);
            }
        }
        self.eddy_viscosity = vec![0.0; self.domain.cell_count()];
        self.sgs = Some(State {
            config,
            filter_width,
            velocity: vec![[0.0; 3]; self.domain.cell_count()],
        });
        Ok(())
    }

    pub(in crate::conjugate) fn has_sgs(&self) -> bool {
        self.sgs.is_some()
    }

    pub(in crate::conjugate) fn update_sgs(&mut self, gate: &CancelGate) -> Result<(), ChtError> {
        let Some(mut state) = self.sgs.take() else {
            return Ok(());
        };
        // Restore the configuration/cache even on refusal; public marches
        // return an error rather than publishing a partially advanced result.
        let result = (|| {
            for cell in 0..self.domain.cell_count() {
                if cell % 64 == 0 {
                    poll(gate)?;
                }
                state.velocity[cell] = self.cell_velocity(cell);
            }
            let nu = self.mu / self.rho;
            for cell in 0..self.domain.cell_count() {
                if cell % 64 == 0 {
                    poll(gate)?;
                }
                if !self.domain.is_fluid(cell) {
                    self.eddy_viscosity[cell] = 0.0;
                    continue;
                }
                let gradient = self.sgs_gradient(cell, &state.velocity);
                let eddy = match state.config.model {
                    SgsModel::Vreman { cv } => vreman_eddy_viscosity(
                        cv,
                        self.domain.widths(cell).map(|d| d * state.config.filter_width_factor),
                        gradient,
                    ),
                    _ => eddy_viscosity(state.config.model, state.filter_width[cell], gradient),
                }
                .map_err(|error| refuse(format!("cell {cell}: {error}")))?;
                let ratio = eddy / nu;
                if !ratio.is_finite() || ratio > state.config.maximum_eddy_to_molecular {
                    return Err(refuse(format!(
                        "cell {cell}: eddy/molecular ratio {ratio} exceeds {}",
                        state.config.maximum_eddy_to_molecular
                    )));
                }
                self.eddy_viscosity[cell] = eddy;
            }
            Ok(())
        })();
        self.sgs = Some(state);
        result
    }

    fn sgs_sample(
        &self,
        cell: usize,
        direction: usize,
        component: usize,
        velocity: &[[f64; 3]],
    ) -> (f64, Sample) {
        let axis = direction / 2;
        let half = 0.5 * self.domain.widths(cell)[axis];
        match self.domain.neighbor(cell, direction) {
            Some(other) if self.domain.is_fluid(other) => (
                half + 0.5 * self.domain.widths(other)[axis],
                Sample::Value(velocity[other][component]),
            ),
            Some(_) => (half, Sample::Value(0.0)),
            None => (
                half,
                match self.faces[direction] {
                    FvBoundary::Wall { velocity } | FvBoundary::Inlet { velocity } => {
                        Sample::Value(velocity[component])
                    }
                    FvBoundary::Symmetry | FvBoundary::Outlet | FvBoundary::Pressure { .. } => {
                        Sample::ZeroGradient
                    }
                },
            ),
        }
    }

    fn sgs_gradient(&self, cell: usize, velocity: &[[f64; 3]]) -> [[f64; 3]; 3] {
        let at = self.domain.coords(cell);
        let mut gradient = [[0.0; 3]; 3];
        for (component, row) in gradient.iter_mut().enumerate() {
            for (axis, value) in row.iter_mut().enumerate() {
                *value = if component == axis {
                    (self.vel[axis][self.cell_face(axis, at, true)]
                        - self.vel[axis][self.cell_face(axis, at, false)])
                        / self.domain.width(axis, at[axis])
                } else {
                    derivative(
                        velocity[cell][component],
                        self.sgs_sample(cell, 2 * axis, component, velocity),
                        self.sgs_sample(cell, 2 * axis + 1, component, velocity),
                    )
                };
            }
        }
        gradient
    }
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod vreman_tests;
