//! Original primitive gas fields in the continuous multimodal wave chart.
//!
//! At uniform static pressure and zero mean flow, the local pressure/volume
//! flow equations use the original rho(x), c(x), viscosity and thermal law.
//! No independently interpolated acoustic coefficients or additional fitted
//! attenuation are needed. The transverse basis and its geometric coupling
//! are unchanged. This is the same finite-mode, locally cylindrical loss
//! approximation as the uniform continuous modal owner, not a full 3-D
//! thermoviscous or convected-flow model.

use super::*;
use crate::profile::MAX_RICCATI_STEPS;
use crate::profile::gas_field::GasFieldDuctError;
use fs_material::gas_field::AxialGasField;
use std::sync::Arc;

mod molecular;
pub use molecular::MolecularModalGasFieldProfile;

#[derive(Clone, Debug)]
struct QuadraturePoint {
    x_m: f64,
    gas: GasState,
    radius: f64,
    flare: f64,
}

#[derive(Clone, Debug)]
struct FieldStep {
    h: f64,
    left_m: f64,
    points: [QuadraturePoint; 2],
}

/// A smooth axisymmetric profile bound to its complete original gas field.
///
/// Material states and geometry are evaluated once at the actual two-point
/// Magnus quadrature of every span. Frequency sweeps reuse these evaluations;
/// source changes reuse the resulting [`PreparedModalProfile`]. Every
/// retained radial mode and signed geometric coupling participates.
///
/// Nonuniform fields split exactly at primitive material stations and spline
/// knots. `steps_per_span` refines each smooth interval, independently of the
/// retained mode count. Uniform fields delegate the original modal path in
/// its original operation order. The original field and source receipts stay
/// available after preparation; no averaged gas is substituted.
#[derive(Clone, Debug)]
pub struct ModalGasFieldProfile {
    profile: RadiusProfile,
    field: Arc<AxialGasField>,
    inlet: GasState,
    outlet: GasState,
    steps_per_span: RiccatiSteps,
    /// Ordered from mouth to inlet, as in the original wave-chart solve.
    quadrature: Vec<FieldStep>,
}

impl ModalGasFieldProfile {
    /// Admit geometry and all local material states before publishing a binding.
    /// The field must cover the exact profile length. Source-domain, finite
    /// coordinate-resolution and total integration-work refusals propagate.
    pub fn new(
        profile: &RadiusProfile,
        field: &AxialGasField,
        steps_per_span: RiccatiSteps,
    ) -> Result<Self, GasFieldDuctError> {
        Self::new_controlled(profile, field, steps_per_span, &mut || false)
    }

    /// As [`Self::new`], with cancellation before each pair of original
    /// material evaluations and before publication. No partial binding escapes.
    pub fn new_controlled(
        profile: &RadiusProfile,
        field: &AxialGasField,
        steps_per_span: RiccatiSteps,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, GasFieldDuctError> {
        poll(cancelled, 0)?;
        steps_per_span.validate()?;
        let length = profile.length();
        if length != field.length_m() {
            return Err(bad("modal gas field must cover the exact profile length").into());
        }
        let inlet = *field.state_at(0.0)?.state();
        let outlet = *field.state_at(length)?.state();
        let mut result = Self {
            profile: profile.clone(),
            field: Arc::new(field.clone()),
            inlet,
            outlet,
            steps_per_span,
            quadrature: Vec::new(),
        };
        if field.uniform_state().is_some() {
            poll(cancelled, 0)?;
            return Ok(result);
        }
        let mut boundaries = field.stations_m().to_vec();
        if let Some(stations) = profile.stations() {
            let origin = stations[0].0;
            boundaries.extend(
                stations[1..stations.len() - 1]
                    .iter()
                    .map(|&(x, _)| x - origin),
            );
        }
        boundaries.sort_by(f64::total_cmp);
        boundaries.dedup();
        let count = (boundaries.len() - 1)
            .checked_mul(steps_per_span.0)
            .filter(|&count| count <= MAX_RICCATI_STEPS)
            .ok_or_else(|| bad("modal gas-field integration exceeds MAX_RICCATI_STEPS"))?;
        result.quadrature.reserve(count);
        let root3 = det::sqrt(3.0);
        let (c1, c2) = (0.5 - root3 / 6.0, 0.5 + root3 / 6.0);
        for span in boundaries.windows(2).rev() {
            for j in (0..steps_per_span.0).rev() {
                poll(cancelled, result.quadrature.len())?;
                let left = span[0] + (span[1] - span[0]) * (j as f64 / steps_per_span.0 as f64);
                let right = if j + 1 == steps_per_span.0 {
                    span[1]
                } else {
                    span[0] + (span[1] - span[0]) * ((j + 1) as f64 / steps_per_span.0 as f64)
                };
                let h = left - right;
                let positions = [right + c1 * h, right + c2 * h];
                if !(h.is_finite()
                    && h < 0.0
                    && positions[1] > left
                    && positions[0] < right
                    && positions[1] < positions[0])
                {
                    return Err(bad("modal gas-field quadrature is numerically unresolved").into());
                }
                let point = |x| -> Result<QuadraturePoint, GasFieldDuctError> {
                    let gas = *field.state_at(x)?.state();
                    let radius = profile.radius(x)?;
                    let flare = profile.slope(x)? / radius;
                    if !flare.is_finite() {
                        return Err(bad("modal gas-field geometric coupling is not finite").into());
                    }
                    Ok(QuadraturePoint { x_m: x, gas, radius, flare })
                };
                result.quadrature.push(FieldStep {
                    h,
                    left_m: left,
                    points: [point(positions[0])?, point(positions[1])?],
                });
            }
        }
        poll(cancelled, result.quadrature.len())?;
        Ok(result)
    }

    /// Original smooth geometry, with every authored spline station.
    pub const fn profile(&self) -> &RadiusProfile {
        &self.profile
    }
    /// Original primitive field and its complete material source owner.
    pub fn field(&self) -> &AxialGasField {
        &self.field
    }
    /// Original gas at the inlet, used for the positive real wave reference.
    pub const fn inlet_gas(&self) -> &GasState {
        &self.inlet
    }
    /// Original gas at the mouth, used for local radiation and continuation.
    pub const fn outlet_gas(&self) -> &GasState {
        &self.outlet
    }
    /// Actual number of Magnus steps. A source or frequency query never
    /// silently refines, averages or discards these spans.
    pub fn steps_used(&self) -> usize {
        if self.quadrature.is_empty() {
            self.steps_per_span.0
        } else {
            self.quadrature.len()
        }
    }

    /// Prepare an analytic, scalar tabulated or matched outlet in the original
    /// outlet gas. Ideal closed/release boundaries act on every mode; scalar
    /// loads retain the original matched higher-mode convention.
    pub fn prepare(
        &self,
        omega: f64,
        loss: LossModel,
        load: MouthLoad<'_>,
        modes: ModeCount,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedModalProfile, DuctError> {
        self.prepare_load(omega, loss, load.into(), modes, cancelled)
    }

    /// Prepare the complete flanged radiation matrix in its explicit exterior
    /// gas. Exterior admission compares against the local mouth gas.
    pub fn prepare_with_flanged_mouth(
        &self,
        exterior: &GasState,
        omega: f64,
        loss: LossModel,
        modes: ModeCount,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedModalProfile, DuctError> {
        self.prepare_load(
            omega,
            loss,
            ProfileLoad::Flanged(exterior),
            modes,
            cancelled,
        )
    }

    /// Prepare an original full mouth matrix with all signed mutual loads.
    /// Its frequency, radius, mode count and source gas must match the mouth.
    pub fn prepare_with_mouth_matrix(
        &self,
        omega: f64,
        loss: LossModel,
        load: &ModalMouthImpedance,
        modes: ModeCount,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedModalProfile, DuctError> {
        poll(cancelled, 0)?;
        load.admit(&self.outlet, omega, self.profile.outlet_radius(), modes.0)?;
        self.prepare_load(omega, loss, ProfileLoad::Matrix(load), modes, cancelled)
    }

    fn prepare_load(
        &self,
        omega: f64,
        loss: LossModel,
        load: ProfileLoad<'_>,
        modes: ModeCount,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedModalProfile, DuctError> {
        poll(cancelled, 0)?;
        modes.validate()?;
        let run = if self.quadrature.is_empty() {
            propagate_controlled(
                &self.profile,
                &self.inlet,
                omega,
                loss,
                load,
                modes.0,
                self.steps_per_span,
                true,
                cancelled,
            )?
        } else {
            propagate_resolved(
                &self.profile,
                &self.inlet,
                &self.outlet,
                omega,
                loss,
                None,
                load,
                modes.0,
                self.quadrature.len(),
                true,
                cancelled,
                &mut |j, basis, zref| {
                    let step = &self.quadrature[j];
                    let local = |p: &QuadraturePoint| {
                        generator_at_radius(&p.gas, omega, loss, basis, zref, p.radius, p.flare)
                    };
                    Ok(PropagationStep {
                        h: step.h,
                        left_m: step.left_m,
                        local: [local(&step.points[0])?, local(&step.points[1])?],
                    })
                },
            )?
        };
        PreparedModalProfile::from_propagated(
            &self.inlet,
            &self.outlet,
            Some(Arc::clone(&self.field)),
            omega,
            modes,
            run,
        )
    }
}
