//! Complete objective derivatives through the original thermoelastic equilibrium.
//!
//! The native solver owns both primal and adjoint equilibrium. The same source
//! potential supplies stress and its temperature partial; the native residual
//! transpose supplies C'(T), integrated expansion, B-bar, ghost and Nitsche work.
//! No Krylov iteration is differentiated and no derivative matrix is assembled.

use std::ops::ControlFlow;

use fs_solver::krylov::ResidualClaim;
use fs_solver::op::LinearOp;

use super::ThermoelasticResponse3;
use crate::elastic3::adaptive::AdaptiveElasticity3;
use crate::elastic3::{CutElasticity3, ElasticityError3, ElasticitySolution3};

/// Bounds shared by the physical solve and its algebraic adjoint.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermoelasticAdjointSettings3 {
    pub relative_residual_tolerance: f64,
    /// Aggregate native CG iterations across BOTH solves.
    pub maximum_iterations: usize,
    pub poll_iterations: usize,
    /// Actual and constitutive total-strain Frobenius ceiling in (0,1).
    /// Original source domains independently constrain elastic strain.
    pub maximum_kinematic_strain: f64,
    /// Native point allowance, including selected supports and ghost work
    /// where required by the residual temperature transpose.
    pub maximum_points: usize,
}

/// The same pure embedded motion used by equilibrium and its derivative.
/// The slope is metres per kelvin of the owning cell's temperature.
/// Strong box clamps remain homogeneous.
#[derive(Clone, Copy)]
pub struct ThermoelasticMotion3<'a> {
    pub displacement: &'a dyn Fn([f64; 3], [f64; 3]) -> [f64; 3],
    pub temperature_derivative: Option<&'a dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
}

/// A scalar objective and its partials at the accepted physical state.
///
/// The callback defines J(u, sigma, T). Displacement entries use the original
/// independent coordinates. Optional stress cotangents pair DIRECTLY with the
/// physical [xx, yy, zz, xy, yz, zx] components at each returned bulk point;
/// the caller includes all desired quadrature weights. They are derivatives
/// with respect to the scaled physical stress, not the reference stress.
///
/// An explicit temperature partial holds u and sigma fixed and has one entry
/// per original active cell. Do not include C'(T), expansion or equilibrium
/// effects here: this owner adds those exactly once. Missing optional partials
/// mean exact zero. Nonlinear objectives may compute these partials from the
/// supplied accepted displacement and native response. A sampled stress
/// objective is not a continuum peak-stress bound.
#[derive(Debug, Clone)]
pub struct ThermoelasticObjective3 {
    pub value: f64,
    pub displacement_derivative: Vec<f64>,
    pub stress_derivatives: Option<Vec<[f64; 6]>>,
    pub explicit_temperature_derivative: Option<Vec<f64>>,
}

/// Accepted primal state and total equilibrium objective derivatives.
///
/// Geometry, original reference density, cell scales, support assignment and
/// independent external force stay fixed during temperature differentiation.
/// force_derivative holds the authored objective and temperatures fixed.
/// Explicit dependence of the objective on force, changing density/body force,
/// geometry, source-selection or a heat-equation adjoint belongs to the caller.
#[derive(Debug, Clone)]
pub struct ThermoelasticSensitivity3 {
    pub solution: ElasticitySolution3,
    pub value: f64,
    pub temperature_derivative: Vec<f64>,
    pub force_derivative: Vec<f64>,
    pub adjoint_iterations: usize,
    pub adjoint_residual: ResidualClaim,
}

#[derive(Clone, Copy)]
enum Owner<'a> {
    Cartesian(&'a CutElasticity3),
    Adaptive(&'a AdaptiveElasticity3),
}

// Both public entries compose the same native operations; no second assembly
// or independent hanging-node reconstruction is introduced.
macro_rules! native {
    ($owner:expr, $method:ident $(, $argument:expr)* $(,)?) => {
        match $owner {
            Owner::Cartesian(operator) => operator.$method($($argument),*),
            Owner::Adaptive(operator) => operator.$method($($argument),*),
        }
    };
}

fn invalid(message: &'static str) -> ElasticityError3 {
    ElasticityError3::Invalid(message)
}

fn poll(
    completed: usize,
    checkpoint: &mut impl FnMut(usize) -> ControlFlow<()>,
) -> Result<(), ElasticityError3> {
    if checkpoint(completed).is_break() {
        Err(ElasticityError3::Cancelled)
    } else {
        Ok(())
    }
}

fn finite(
    values: &[f64],
    completed: usize,
    checkpoint: &mut impl FnMut(usize) -> ControlFlow<()>,
) -> Result<(), ElasticityError3> {
    for chunk in values.chunks(192) {
        poll(completed, checkpoint)?;
        if chunk.iter().any(|value| !value.is_finite()) {
            return Err(invalid(
                "nonfinite thermoelastic objective data or derivative",
            ));
        }
    }
    Ok(())
}

fn evaluate(
    owner: Owner<'_>,
    external_force: &[f64],
    motion: Option<ThermoelasticMotion3<'_>>,
    settings: ThermoelasticAdjointSettings3,
    objective: impl FnOnce(
        &[f64],
        &[ThermoelasticResponse3],
    ) -> Result<ThermoelasticObjective3, ElasticityError3>,
    mut checkpoint: impl FnMut(usize) -> ControlFlow<()>,
) -> Result<ThermoelasticSensitivity3, ElasticityError3> {
    poll(0, &mut checkpoint)?;
    let n = native!(owner, n);
    let cells = native!(owner, cells);
    if external_force.len() != n
        || !settings.relative_residual_tolerance.is_finite()
        || settings.relative_residual_tolerance <= 0.0
        || settings.relative_residual_tolerance >= 1.0
        || settings.poll_iterations == 0
        || !settings.maximum_kinematic_strain.is_finite()
        || settings.maximum_kinematic_strain <= 0.0
        || settings.maximum_kinematic_strain >= 1.0
    {
        return Err(invalid(
            "invalid thermoelastic adjoint controls or force shape",
        ));
    }
    finite(external_force, 0, &mut checkpoint)?;
    // Admit the complete derivative point work before preparation, solves or
    // user callbacks, using the residual owner's identical preflight.
    let raw = match owner {
        Owner::Cartesian(operator) => operator,
        Owner::Adaptive(operator) => &operator.raw,
    };
    raw.admit_thermal_temperature_points(settings.maximum_points, || checkpoint(0))?;

    let prepared = native!(
        owner,
        prepare_thermoelastic,
        settings.maximum_points,
        || checkpoint(0)
    )?;
    let mut force = native!(owner, eigenstrain_load, &prepared, || checkpoint(0))?;
    let motion_force = match motion {
        Some(motion) => Some(native!(
            owner,
            prescribed_displacement_load,
            motion.displacement,
            || checkpoint(0)
        )?),
        None => None,
    };
    for i in 0..n {
        if i % 192 == 0 {
            poll(0, &mut checkpoint)?;
        }
        force[i] += external_force[i];
        if let Some(motion_force) = &motion_force {
            force[i] += motion_force[i];
        }
        if !force[i].is_finite() {
            return Err(invalid("thermoelastic equilibrium force overflow"));
        }
    }
    let solution = native!(
        owner,
        solve_controlled,
        &force,
        settings.relative_residual_tolerance,
        settings.maximum_iterations,
        settings.poll_iterations,
        &mut checkpoint
    )?;
    let primal_iterations = solution.iterations();
    let response = native!(
        owner,
        thermoelastic_response,
        solution.coefficients(),
        settings.maximum_kinematic_strain,
        settings.maximum_points,
        || checkpoint(primal_iterations)
    )?;
    poll(primal_iterations, &mut checkpoint)?;
    let objective = objective(solution.coefficients(), &response)?;
    poll(primal_iterations, &mut checkpoint)?;
    if !objective.value.is_finite()
        || objective.displacement_derivative.len() != n
        || objective
            .stress_derivatives
            .as_ref()
            .is_some_and(|values| values.len() != response.len())
        || objective
            .explicit_temperature_derivative
            .as_ref()
            .is_some_and(|values| values.len() != cells)
    {
        return Err(invalid(
            "thermoelastic objective value or partial shape is invalid",
        ));
    }
    let mut dual_force = objective.displacement_derivative;
    finite(&dual_force, primal_iterations, &mut checkpoint)?;
    let mut temperature = objective
        .explicit_temperature_derivative
        .unwrap_or_else(|| vec![0.0; cells]);
    finite(&temperature, primal_iterations, &mut checkpoint)?;
    if let Some(mut stress_derivatives) = objective.stress_derivatives {
        // Differentiate physical scaled stress. Native reference stress's
        // transpose supplies the same material and B-bar basis on either grid.
        let scales = native!(owner, scales);
        for (point, derivative) in response.iter().zip(&mut stress_derivatives) {
            poll(primal_iterations, &mut checkpoint)?;
            if derivative.iter().any(|value| !value.is_finite()) {
                return Err(invalid("nonfinite thermoelastic stress cotangent"));
            }
            temperature[point.mechanical.cell] += derivative
                .iter()
                .zip(point.stress_temperature_derivative)
                .map(|(weight, slope)| weight * slope)
                .sum::<f64>();
            for value in derivative {
                *value *= scales[point.mechanical.cell];
                if !value.is_finite() {
                    return Err(invalid("scaled thermoelastic stress cotangent overflow"));
                }
            }
        }
        let stress_force = native!(
            owner,
            reference_bulk_stress_pullback,
            &stress_derivatives,
            || checkpoint(primal_iterations)
        )?;
        for (i, (force, contribution)) in dual_force.iter_mut().zip(stress_force).enumerate() {
            if i % 192 == 0 {
                poll(primal_iterations, &mut checkpoint)?;
            }
            *force += contribution;
        }
        finite(&dual_force, primal_iterations, &mut checkpoint)?;
        finite(&temperature, primal_iterations, &mut checkpoint)?;
    }

    let remaining = settings
        .maximum_iterations
        .checked_sub(primal_iterations)
        .ok_or_else(|| invalid("native primal exceeded thermoelastic iteration budget"))?;
    let adjoint = native!(
        owner,
        solve_controlled,
        &dual_force,
        settings.relative_residual_tolerance,
        remaining,
        settings.poll_iterations,
        |iterations| checkpoint(primal_iterations + iterations)
    )?;
    let completed = primal_iterations + adjoint.iterations();
    let residual_derivative = match motion {
        Some(motion) => native!(
            owner,
            thermal_residual_temperature_pullback_with_motion,
            solution.coefficients(),
            adjoint.coefficients(),
            motion.displacement,
            motion.temperature_derivative,
            settings.maximum_kinematic_strain,
            settings.maximum_points,
            || checkpoint(completed)
        ),
        None => native!(
            owner,
            thermal_residual_temperature_pullback,
            solution.coefficients(),
            adjoint.coefficients(),
            settings.maximum_kinematic_strain,
            settings.maximum_points,
            || checkpoint(completed)
        ),
    }?;
    // R=K(T)u-b(T); K^T z=J_u. The implicit derivative is J_T-z^T R_T.
    for (i, (derivative, residual)) in temperature.iter_mut().zip(residual_derivative).enumerate() {
        if i % 64 == 0 {
            poll(completed, &mut checkpoint)?;
        }
        *derivative -= residual;
    }
    finite(&temperature, completed, &mut checkpoint)?;
    poll(completed, &mut checkpoint)?;
    Ok(ThermoelasticSensitivity3 {
        value: objective.value,
        temperature_derivative: temperature,
        force_derivative: adjoint.coefficients().to_vec(),
        adjoint_iterations: adjoint.iterations(),
        adjoint_residual: adjoint.residual_claim(),
        solution,
    })
}

impl CutElasticity3 {
    /// Solve and differentiate an objective of native displacement, stress and T.
    ///
    /// Source-temperature changes include the full equilibrium response, not
    /// only the local constitutive partial. The objective callback executes once
    /// after physical primal admission. The dual is algebraic and is never
    /// subjected to a physical strain-amplitude ceiling. Both solves must pass
    /// their original true-residual gates within one aggregate iteration budget.
    /// Successful or refused evaluation leaves this operator unchanged.
    ///
    /// This assumes the original supported elastic problem is well posed; a
    /// small residual alone does not establish uniqueness or continuum accuracy.
    pub fn solve_thermoelastic_adjoint(
        &self,
        external_force: &[f64],
        motion: Option<ThermoelasticMotion3<'_>>,
        settings: ThermoelasticAdjointSettings3,
        objective: impl FnOnce(
            &[f64],
            &[ThermoelasticResponse3],
        ) -> Result<ThermoelasticObjective3, ElasticityError3>,
        checkpoint: impl FnMut(usize) -> ControlFlow<()>,
    ) -> Result<ThermoelasticSensitivity3, ElasticityError3> {
        evaluate(
            Owner::Cartesian(self),
            external_force,
            motion,
            settings,
            objective,
            checkpoint,
        )
    }
}

impl AdaptiveElasticity3 {
    /// The same complete thermoelastic objective adjoint on original master
    /// coordinates. Native hanging-node reconstruction and its exact transpose
    /// govern primal, stress and temperature work; no alternate projection is
    /// introduced. Semantics and bounds match
    /// [CutElasticity3::solve_thermoelastic_adjoint].
    pub fn solve_thermoelastic_adjoint(
        &self,
        external_force: &[f64],
        motion: Option<ThermoelasticMotion3<'_>>,
        settings: ThermoelasticAdjointSettings3,
        objective: impl FnOnce(
            &[f64],
            &[ThermoelasticResponse3],
        ) -> Result<ThermoelasticObjective3, ElasticityError3>,
        checkpoint: impl FnMut(usize) -> ControlFlow<()>,
    ) -> Result<ThermoelasticSensitivity3, ElasticityError3> {
        evaluate(
            Owner::Adaptive(self),
            external_force,
            motion,
            settings,
            objective,
            checkpoint,
        )
    }
}
