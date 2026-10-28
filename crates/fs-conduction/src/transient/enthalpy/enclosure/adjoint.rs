//! Complete reflected-radiation Jacobian in specific-enthalpy coordinates.
//!
//! For patch mean W*T(h), flux q and P1 flux distribution G, the additional
//! residual tangent is dt*G*(dq/dT_patch)*W*diag(dT/dh). It is retained as one
//! nodal factor pair per patch. Latent plateaus have exact zero temperature
//! slope; direct enthalpy history is never divided by that slope.

use fs_exec::Cx;

use super::{
    bind, convection, initial_target, means, prescribed_neumann_heat, prepared_radiosity,
    shifted_boundary,
};
use crate::{
    AmbientRadiationConfig, ConductionProblem, ThermalInterfaces,
    assemble::ASSEMBLY_TILE,
    radiation::{GrayDiffuseEnclosure, PreparedRadiosity, STEFAN_BOLTZMANN_W_M2_K4},
    transient::enthalpy::{
        EnthalpyBackwardEuler, EnthalpyStepConfig, EnthalpyStepSolution, StepContext,
        adjoint::{EnthalpyAdjointError, EnthalpyStepLinearization},
        poll,
    },
};

impl<'m> EnthalpyBackwardEuler<'m, '_> {
    /// Recheck an accepted enclosure endpoint and bind its complete h tangent.
    /// The supplied h field is authoritative; temperatures, source/boundary
    /// accounting and residual telemetry are independently reconstructed.
    /// `problem.boundary` contains the ORIGINAL convection and prescribed flux,
    /// not shifted radiation rows. Vacuum traces retain the same reflected
    /// radiation tangent; they do not create a convective control port.
    ///
    /// Geometry, view factors, emissivities, phase charts and reference masses
    /// are fixed. Source/history/convection derivatives retain all reflected
    /// radiation. `max_feedback_entries` caps 2*vertices*patches scalars.
    /// No nonlinear or Krylov iteration is differentiated.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn linearize_accepted_with_enclosure_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        step_config: EnthalpyStepConfig,
        enclosure: &GrayDiffuseEnclosure,
        radiation_config: AmbientRadiationConfig,
        accepted: EnthalpyStepSolution,
        max_feedback_entries: usize,
    ) -> Result<EnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        self.admit_step(cx, problem, old_h, dt_s, step_config)?;
        radiation_config.validate()?;
        let n = self.masses.len();
        let p = enclosure.surfaces().len();
        let entries = n.checked_mul(2).and_then(|v| v.checked_mul(p));
        if p == 0
            || p > 64
            || entries.is_none_or(|v| v > max_feedback_entries)
            || entries
                .and_then(|v| v.checked_mul(std::mem::size_of::<f64>()))
                .is_none_or(|v| isize::try_from(v).is_err())
        {
            return Err(EnthalpyAdjointError::InvalidInput(
                "enclosure h tangent exceeds its two-vector cap or 1..=64 patches",
            ));
        }
        if accepted.specific_enthalpy_j_kg.len() != n {
            return Err(EnthalpyAdjointError::InvalidInput(
                "enclosure endpoint length mismatch",
            ));
        }
        let context = StepContext {
            storage: self,
            cx,
            problem,
            interfaces,
            old: old_h,
            dt: dt_s,
            pointwise_radiation: None,
        };
        let temperature = context.temperatures(&accepted.specific_enthalpy_j_kg)?;
        let patches = bind(cx, problem, enclosure, &temperature)?;
        let prepared = enclosure.prepare_radiosity(cx)?;
        let target = initial_target(&context, step_config, &prepared, &patches, radiation_config)?;
        let means = means(cx, problem, enclosure, &temperature)?;
        let fresh = prepared_radiosity(cx, &prepared, &means, radiation_config)?;
        let boundary = shifted_boundary(problem, &patches, &fresh)?;
        let mut transport = self.linearize_accepted_with_target(
            cx,
            ConductionProblem {
                boundary: &boundary,
                ..problem
            },
            interfaces,
            old_h,
            dt_s,
            step_config,
            accepted,
            Some(target),
        )?;
        // Reuse the primal factorization for every derivative column while
        // independently checking the original derivative equations.
        let slopes = flux_temperature_jacobian(
            cx,
            &prepared,
            &means,
            step_config.newton.relative_tolerance,
        )?;
        for column in 0..p {
            poll(cx, column)?;
            let mut left = vec![0.0; n];
            let mut weights = vec![0.0; n];
            for (row, surface) in enclosure.surfaces().iter().enumerate() {
                let coefficient = checked(dt_s * slopes[row][column])?;
                for (ordinal, &slot) in surface.face_slots().iter().enumerate() {
                    if ordinal % ASSEMBLY_TILE == 0 {
                        poll(cx, ordinal)?;
                    }
                    let face = &problem.mesh.boundary()[slot];
                    for &vertex in &face.vertices {
                        let vertex = vertex as usize;
                        left[vertex] = checked(left[vertex] + coefficient * face.area / 3.0)?;
                        if row == column {
                            weights[vertex] =
                                checked(weights[vertex] + face.area / surface.area_m2() / 3.0)?;
                        }
                    }
                }
            }
            let right = transport.temperature_pullback(cx, &weights)?;
            transport.add_radiation_feedback(cx, left, right)?;
        }
        let (_, convective_out) = convection(cx, problem, transport.primal())?;
        let prescribed_neumann_out = prescribed_neumann_heat(cx, problem)?;
        let primal = transport.primal();
        let energy = checked(
            primal.stored_energy_change_j
                + dt_s
                    * (prescribed_neumann_out + convective_out + fresh.enclosure_energy_closure_w
                        - primal.source_w - primal.dirichlet_in_w),
        )?;
        if energy.abs() > step_config.energy_tolerance_j {
            return Err(super::super::EnthalpyError::EnergyBalance {
                residual_j: energy,
                tolerance_j: step_config.energy_tolerance_j,
            }
            .into());
        }
        transport.retain_convection_boundary(problem.boundary);
        poll(cx, n)?;
        Ok(transport)
    }
}

fn flux_temperature_jacobian(
    cx: &Cx<'_>,
    prepared: &PreparedRadiosity<'_>,
    temperatures: &[f64],
    relative: f64,
) -> Result<Vec<Vec<f64>>, EnthalpyAdjointError> {
    let enclosure = prepared.enclosure();
    let n = enclosure.surfaces().len();
    let factors = enclosure.view_factors().factors();
    let mut matrix = vec![vec![0.0; n]; n];
    let mut slopes = Vec::with_capacity(n);
    for (row, surface) in enclosure.surfaces().iter().enumerate() {
        poll(cx, row)?;
        let t = temperatures[row];
        if let Some((lower, upper)) = surface.emissivity().temperature_validity_k()
            && (t <= lower || t >= upper)
        {
            return Err(EnthalpyAdjointError::InvalidInput(
                "enclosure derivative requires interior emissivity temperature validity",
            ));
        }
        let epsilon = surface.emissivity().value();
        slopes.push(checked(
            4.0 * epsilon * STEFAN_BOLTZMANN_W_M2_K4 * t * t * t,
        )?);
        for (column, &factor) in factors[row].iter().enumerate() {
            matrix[row][column] =
                checked((if row == column { 1.0 } else { 0.0 }) - (1.0 - epsilon) * factor)?;
        }
    }
    let mut result = vec![vec![0.0; n]; n];
    for column in 0..n {
        poll(cx, column)?;
        let mut rhs = vec![0.0; n];
        rhs[column] = slopes[column];
        let response = prepared.solve_rhs(cx, &rhs)?;
        let scale = slopes[column].abs().max(f64::MIN_POSITIVE);
        for row in 0..n {
            let mut applied = 0.0;
            let mut irradiation = 0.0;
            for k in 0..n {
                applied = checked(matrix[row][k].mul_add(response[k], applied))?;
                irradiation = checked(factors[row][k].mul_add(response[k], irradiation))?;
            }
            if checked((applied - rhs[row]) / scale)?.abs() > relative {
                return Err(EnthalpyAdjointError::InvalidInput(
                    "enclosure derivative system failed its original equation residual",
                ));
            }
            result[row][column] = checked(response[row] - irradiation)?;
        }
    }
    Ok(result)
}

fn checked(value: f64) -> Result<f64, EnthalpyAdjointError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(EnthalpyAdjointError::NonFiniteArithmetic)
    }
}
