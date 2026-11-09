//! Native source-temperature mechanics through the original density pipeline.
//!
//! The CutFEM owner differentiates C(T), integrated free expansion, physical
//! stress, B-bar, ghost work and embedded support motion. This adapter composes
//! that scale derivative with the original filter/projection/SIMP map and adds
//! the physical occupancy derivative of the original reference-density force.
//! No mechanical derivative matrix, extra adjoint or heat solve is introduced.

use std::ops::ControlFlow;

use super::{CutDensityStudy3, Sdf3Elasticity};
use crate::{EvaluationStop, SolveControl, SolveWork};
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::thermal::{
    ThermoelasticAdjointSettings3, ThermoelasticDesignObjective3, ThermoelasticDesignSensitivity3,
    ThermoelasticMotion3, ThermoelasticResponse3,
};
use fs_cutfem::elastic3::{CutElasticity3, ElasticityError3};

/// Scalar objective partials at one accepted native thermoelastic state.
#[derive(Debug, Clone)]
pub struct CutThermoelasticObjective3 {
    /// Native partials with respect to u, physical stress, T and stiffness scale.
    pub native: ThermoelasticDesignObjective3,
    /// Explicit derivative with respect to PROJECTED physical occupancy while
    /// u, stress, scale, T and effective external force stay fixed. None is zero.
    pub explicit_occupancy_derivative: Option<Vec<f64>>,
    /// Explicit derivative with respect to the callback's effective external
    /// force while u, stress, occupancy, scale and T stay fixed. For J=f.u this
    /// is u, in addition to the native displacement partial f. None is zero.
    /// Do not also place its body-load chain in the occupancy partial.
    pub explicit_force_derivative: Option<Vec<f64>>,
}

/// Complete physical evaluation and derivatives in the original coordinates.
#[derive(Debug, Clone)]
pub struct CutThermoelasticEvaluation3 {
    /// Native state, scale/T partials and equilibrium-only force derivative.
    /// Native primal and adjoint iteration/residual reports remain unchanged.
    pub sensitivity: ThermoelasticDesignSensitivity3,
    /// Total derivative with respect to the caller's raw filter controls.
    pub gradient: Vec<f64>,
    /// Total derivative with respect to the authored external-force vector,
    /// including its explicit objective partial. Strong-clamp entries are zero.
    pub force_derivative: Vec<f64>,
    pub projected_rho: Vec<f64>,
    pub scales: Vec<f64>,
    pub volume_fraction: f64,
    pub volume_gradient: Vec<f64>,
    /// Cumulative filter/pullback work from the supplied SolveControl.
    /// Native mechanical CG work is separately reported in sensitivity.
    pub work: SolveWork,
}

#[derive(Debug)]
pub enum CutThermoelasticError3 {
    Evaluation(EvaluationStop),
    Native(ElasticityError3),
    Invalid(&'static str),
}
impl From<EvaluationStop> for CutThermoelasticError3 {
    fn from(error: EvaluationStop) -> Self {
        Self::Evaluation(error)
    }
}
impl From<ElasticityError3> for CutThermoelasticError3 {
    fn from(error: ElasticityError3) -> Self {
        Self::Native(error)
    }
}
impl core::fmt::Display for CutThermoelasticError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "native thermoelastic density evaluation refused: {self:?}"
        )
    }
}
impl std::error::Error for CutThermoelasticError3 {}
type Result<T> = core::result::Result<T, CutThermoelasticError3>;

/// Native Cartesian/adaptive thermoelastic dispatch on the sealed SDF owners.
/// The existing native primal/adjoint solver and its explicit budget are used.
pub trait ThermoelasticSdf3Elasticity: Sdf3Elasticity {
    #[doc(hidden)]
    fn native_thermoelastic_design(
        &self,
        force: &[f64],
        motion: Option<ThermoelasticMotion3<'_>>,
        settings: ThermoelasticAdjointSettings3,
        objective: impl FnOnce(
            &[f64],
            &[ThermoelasticResponse3],
        )
            -> core::result::Result<ThermoelasticDesignObjective3, ElasticityError3>,
        checkpoint: impl FnMut(usize) -> ControlFlow<()>,
    ) -> core::result::Result<ThermoelasticDesignSensitivity3, ElasticityError3>;
}
macro_rules! native_backend {
    ($owner:ty) => {
        impl ThermoelasticSdf3Elasticity for $owner {
            fn native_thermoelastic_design(
                &self,
                force: &[f64],
                motion: Option<ThermoelasticMotion3<'_>>,
                settings: ThermoelasticAdjointSettings3,
                objective: impl FnOnce(
                    &[f64],
                    &[ThermoelasticResponse3],
                ) -> core::result::Result<
                    ThermoelasticDesignObjective3,
                    ElasticityError3,
                >,
                checkpoint: impl FnMut(usize) -> ControlFlow<()>,
            ) -> core::result::Result<ThermoelasticDesignSensitivity3, ElasticityError3> {
                self.solve_thermoelastic_design_adjoint(
                    force, motion, settings, objective, checkpoint,
                )
            }
        }
    };
}
native_backend!(CutElasticity3);
native_backend!(AdaptiveElasticity3);

// The source, geometry and number of cells remain immutable during evaluation.
// Restore on every return, including unwind from a caller's objective/checkpoint.
struct RestoreScales<'a, O: Sdf3Elasticity> {
    study: &'a mut CutDensityStudy3<O>,
    previous: Vec<f64>,
}
impl<O: Sdf3Elasticity> Drop for RestoreScales<'_, O> {
    fn drop(&mut self) {
        self.study
            .operator
            .set_scales(&self.previous)
            .expect("previous admitted density scales remain valid");
    }
}

fn partial(
    values: Option<Vec<f64>>,
    expected: usize,
    control: &mut SolveControl<'_>,
) -> Result<Vec<f64>> {
    let values = values.unwrap_or_else(|| vec![0.0; expected]);
    if values.len() != expected {
        return Err(CutThermoelasticError3::Invalid(
            "objective partial shape mismatch",
        ));
    }
    for chunk in values.chunks(192) {
        control.checkpoint("sdf3-thermoelastic-partials")?;
        if chunk.iter().any(|value| !value.is_finite()) {
            return Err(CutThermoelasticError3::Invalid(
                "nonfinite objective partial",
            ));
        }
    }
    Ok(values)
}

impl<O: ThermoelasticSdf3Elasticity> CutDensityStudy3<O> {
    /// Evaluate one native source-temperature experiment and its complete
    /// raw-density derivative on this exact Cartesian or adaptive operator.
    ///
    /// The pure objective callback runs exactly once AFTER physical primal
    /// admission. Arguments are accepted u, native thermal bulk responses,
    /// projected physical occupancy, and the actual external+body force after
    /// strong-clamp masking. Optional partials hold the other arguments and
    /// native scale/T coordinates fixed. The native force adjoint and explicit
    /// force partial BOTH contribute to the body's occupancy transpose.
    ///
    /// Source C(T), integrated expansion and thermal motion derivatives remain
    /// native. Source temperature is independent of raw density; the returned
    /// T derivative holds reference inventory and external force fixed. Changing
    /// source assignments, thermal density, geometry, acceleration or a solved
    /// temperature field requires an additional caller-owned chain rule.
    ///
    /// settings bounds the native primal+adjoint CG work PER EVALUATION and its
    /// physical/source point domain. The original control bounds the filter and
    /// both density pullbacks, retaining its cumulative spent work. It also
    /// supplies every native cancellation checkpoint; no fresh control is made.
    /// Native iterations are not charged again as filter iterations.
    ///
    /// Success, refusal, cancellation and unwinding callbacks restore incoming
    /// scales. Returned fields carry their own evaluated scales. This is a
    /// discrete evaluator, not an optimizer, heat adjoint or continuum bound.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_thermoelastic(
        &mut self,
        rho: &[f64],
        force: &[f64],
        motion: Option<ThermoelasticMotion3<'_>>,
        settings: ThermoelasticAdjointSettings3,
        objective: impl FnOnce(
            &[f64],
            &[ThermoelasticResponse3],
            &[f64],
            &[f64],
        )
            -> core::result::Result<CutThermoelasticObjective3, ElasticityError3>,
        control: &mut SolveControl<'_>,
    ) -> Result<CutThermoelasticEvaluation3> {
        control.checkpoint("sdf3-thermoelastic-start")?;
        if rho.len() != self.cells() || force.len() != self.operator.n() {
            return Err(CutThermoelasticError3::Invalid(
                "complete raw design and native force required",
            ));
        }
        for chunk in rho.chunks(192) {
            control.checkpoint("sdf3-thermoelastic-design")?;
            if chunk
                .iter()
                .any(|r| !r.is_finite() || !(0.0..=1.0).contains(r))
            {
                return Err(CutThermoelasticError3::Invalid(
                    "raw densities must lie in [0,1]",
                ));
            }
        }
        for chunk in force.chunks(192) {
            control.checkpoint("sdf3-thermoelastic-force")?;
            if chunk.iter().any(|f| !f.is_finite()) {
                return Err(CutThermoelasticError3::Invalid(
                    "finite external forces required",
                ));
            }
        }
        let design = self.design(rho, control)?;
        let previous = self.operator.scales().to_vec();
        let restore = RestoreScales {
            study: self,
            previous,
        };
        restore.study.operator.set_scales(&design.scales)?;
        let study = &*restore.study;
        let body = study.material_body_force(&design.projected, control)?;
        let external = study.external_force(force, body.as_deref(), control)?;
        let mut partials = None;
        let mut stopped = None;
        let native = study.operator.native_thermoelastic_design(
            &external,
            motion,
            settings,
            |u, response| {
                let value = objective(u, response, &design.projected, &external)?;
                partials = Some((
                    value.explicit_occupancy_derivative,
                    value.explicit_force_derivative,
                ));
                Ok(value.native)
            },
            |_| match control.checkpoint("sdf3-thermoelastic-native") {
                Ok(()) => ControlFlow::Continue(()),
                Err(stop) => {
                    stopped = Some(stop);
                    ControlFlow::Break(())
                }
            },
        );
        if let Some(stop) = stopped {
            return Err(stop.into());
        }
        let sensitivity = native?;
        let (occupancy, force_partial) = partials.ok_or(CutThermoelasticError3::Invalid(
            "native objective was not evaluated",
        ))?;
        let mut occupancy = partial(occupancy, study.cells(), control)?;
        let mut force_derivative = partial(force_partial, external.len(), control)?;
        for (i, (total, &implicit)) in force_derivative
            .iter_mut()
            .zip(&sensitivity.thermoelastic.force_derivative)
            .enumerate()
        {
            if i % 192 == 0 {
                control.checkpoint("sdf3-thermoelastic-force-pullback")?;
            }
            *total = if study.operator.fixed()[i / 3] {
                0.0
            } else {
                *total + implicit
            };
            if !total.is_finite() {
                return Err(CutThermoelasticError3::Invalid("force derivative overflow"));
            }
        }
        study.add_material_body_gradient(&force_derivative, 1.0, &mut occupancy, control)?;
        let gradient =
            study.physical_gradient(&design, &sensitivity.scale_derivative, &occupancy, control)?;
        let local_volume: Vec<_> = study
            .mass
            .iter()
            .zip(&design.slope)
            .map(|(mass, slope)| mass * slope)
            .collect();
        let volume_gradient = study.pullback(&local_volume, control)?;
        if !design.volume.is_finite() || volume_gradient.iter().any(|v| !v.is_finite()) {
            return Err(CutThermoelasticError3::Invalid(
                "volume derivative overflow",
            ));
        }
        control.checkpoint("sdf3-thermoelastic-publish")?;
        Ok(CutThermoelasticEvaluation3 {
            sensitivity,
            gradient,
            force_derivative,
            projected_rho: design.projected,
            scales: design.scales,
            volume_fraction: design.volume,
            volume_gradient,
            work: control.work(),
        })
    }
}
