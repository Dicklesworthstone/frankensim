//! Two-field comparison through the producer's card-backed radiative secant.
//! Convection references may depend on other walls; radiation never heats air.

use super::{ConductionError, ConductionProblem, Cx, DiscreteGoalComparison,
    LinearConfig, RobinResponse, ThermalInterfaces, checked, invalid, poll, vector};
use crate::radiation::{AmbientRadiationGradient, AmbientRadiationLinearization, AmbientRadiationPatch};
use crate::{ScalarField, ThermalBc};
use std::collections::BTreeSet;

impl AmbientRadiationGradient {
    /// Compare a linear nodal-temperature goal through convection and radiation.
    ///
    /// `problem.boundary` is the ORIGINAL convection/zero-flux partition. In `regions`
    /// order, `reference_convection` supplies [h_c, r_c, dh_c/dmean, dr_c/dmean]
    /// at the reference field and `approximate_convection` supplies [h_c, r_c]
    /// at the approximate field. Units are W/(m^2 K), K, W/(m^2 K^2), and K/K.
    /// The additional row-major m-by-m `reference_feedback` maps all wall means
    /// to convective references. Do not repeat its diagonal in the local slope.
    /// Supply zeros for fixed references. The caller owns the convection law,
    /// evaluates BOTH physical states, and supplies its exact derivative.
    ///
    /// A radiation-only trace retains its original uniform zero-flux boundary.
    /// Its coefficients, coefficient/reference slopes and feedback ROW must
    /// all be zero; positive reference inputs are inert placeholders. Its own
    /// patch supplies the radiative reservoir without manufacturing convection.
    ///
    /// Every radiation patch must occur exactly once in `regions`. Its original
    /// emissivity card, query, domain and reservoir are retained. The actual
    /// secant is re-evaluated at each field's own area-mean wall temperature.
    /// Both the total coefficient slope and weighted-reference slope enter the
    /// genuine nonsymmetric transpose; cross-wall convection feedback is scaled
    /// by h_c/(h_c+h_rad), never by one. Consistent P1 trace mass, smooth k(T),
    /// matching contact, and prescribed temperatures remain in the core owner.
    /// No extra primal solve or differentiated solver iteration is performed.
    ///
    /// The physical reference residual is checked even for a zero goal. The
    /// observed nonlinear remainder is retained separately from signed nodal
    /// residual contributions. This is an Estimated discrete comparison, not a
    /// continuum, uncertainty, stability or physical-validation certificate.
    /// It is the mean-temperature secant model, not pointwise T(x)^4 quadrature.
    /// At most 64 regions and patches are admitted; four nodal factor vectors
    /// per region and the m-by-m feedback must fit `max_feedback_entries`.
    /// That cap is not total allocator/RSS. Failure never mutates input fields.
    /// Only the reference needs a unique emissivity tangent. The approximate
    /// field may reach a supported validity endpoint because it needs physical
    /// secant values, not a derivative or neighboring source-temperature query.
    #[allow(clippy::too_many_arguments)]
    pub fn compare_goal_at(
        cx: &Cx<'_>, problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
        linear: LinearConfig, reference: &[f64], approximate: &[f64], regions: &[&str],
        reference_convection: &[[f64; 4]], approximate_convection: &[[f64; 2]],
        reference_feedback: &[f64], patches: &[AmbientRadiationPatch],
        weights: &[f64], max_feedback_entries: usize,
    ) -> Result<DiscreteGoalComparison, ConductionError> {
        poll(cx, 0)?;
        let m = regions.len();
        let n = problem.mesh.vertex_count();
        if m == 0 || m > 64 || patches.is_empty() || patches.len() > 64
            || reference_convection.len() != m || approximate_convection.len() != m
            || reference_feedback.len() != m*m || m*m > max_feedback_entries
            || n.checked_mul(4).and_then(|v| v.checked_mul(m))
                .is_none_or(|v| v > max_feedback_entries)
        {
            return Err(invalid("radiative goal requires matching law rows, radiation patches and bounded feedback"));
        }
        for field in [reference, approximate, weights] { vector(cx, field, n)?; }
        vector(cx, reference_feedback, m*m)?;
        let mut seen = BTreeSet::new();
        for &name in regions {
            if !seen.insert(name) { return Err(invalid("duplicate radiative-goal region")); }
        }
        seen.clear();
        for patch in patches {
            if !seen.insert(patch.region()) || !regions.contains(&patch.region()) {
                return Err(invalid("every radiation patch needs one unique selected thermal region"));
            }
        }
        let mut selected = [0; 64];
        let mut reference_rows = [(0, 0.0, 0.0); 64];
        let mut approximate_rows = [(0, 0.0, 0.0); 64];
        let mut slopes = [0.0; 64];
        let mut convection_rows = 0;
        let mut feedback = reference_feedback.to_vec();
        for (i, &name) in regions.iter().enumerate() {
            poll(cx, i)?;
            vector(cx, &reference_convection[i], 4)?;
            vector(cx, &approximate_convection[i], 2)?;
            let [hc, tc, dhc, dtc] = reference_convection[i];
            let [ha, ta] = approximate_convection[i];
            if tc <= 0.0 || ta <= 0.0 {
                return Err(invalid("radiative goal needs positive absolute reference values"));
            }
            let region = problem.boundary.region_names().iter().position(|r| r == name)
                .ok_or_else(|| invalid("radiative-goal region is absent from the original partition"))?;
            match problem.boundary.conditions()[region] {
                ThermalBc::Robin { htc: ScalarField::Uniform(_), t_ref: ScalarField::Uniform(_) }
                    if hc > 0.0 && ha > 0.0 => {
                    reference_rows[convection_rows] = (region, hc, tc);
                    approximate_rows[convection_rows] = (region, ha, ta);
                    convection_rows += 1;
                }
                ThermalBc::Neumann { outward_flux: ScalarField::Uniform(flux) }
                    if flux == 0.0 && hc == 0.0 && ha == 0.0 && dhc == 0.0 && dtc == 0.0
                        && patches.iter().any(|patch| patch.region() == name)
                        && reference_feedback[i*m..(i+1)*m].iter().all(|&v| v == 0.0) => {}
                _ => return Err(invalid(
                    "retain uniform convection or radiation-only zero-flux ownership without invented convection or feedback")),
            }
            selected[i] = region;
            slopes[i] = dhc;
            // Radiation weights the entire original reference derivative.
            // Its own local weighted-reference slope is added separately.
            feedback[i*m+i] = checked(reference_feedback[i*m+i] + dtc)?;
        }
        let reference_boundary = problem.boundary
            .with_uniform_robin_replacements(&reference_rows[..convection_rows])?;
        let approximate_boundary = problem.boundary
            .with_uniform_robin_replacements(&approximate_rows[..convection_rows])?;
        let bound = AmbientRadiationLinearization::bind(cx,
            ConductionProblem { boundary: &reference_boundary, ..problem }, reference,
            regions, &slopes[..m], Some(&feedback), patches, max_feedback_entries)?;
        let other_boundary = AmbientRadiationLinearization::physical_boundary_at(cx,
            ConductionProblem { boundary: &approximate_boundary, ..problem }, approximate,
            regions, patches, max_feedback_entries)?;
        let mut points = [[0.0; 4]; 64];
        let mut other = [[0.0; 2]; 64];
        for (i, &region) in selected[..m].iter().enumerate() {
            poll(cx, i)?;
            let ThermalBc::Robin { htc: ScalarField::Uniform(h), t_ref: ScalarField::Uniform(r) }
                = bound.boundary().conditions()[region]
            else { return Err(invalid("radiation binding lost its combined reference Robin row")); };
            let ThermalBc::Robin { htc: ScalarField::Uniform(ha), t_ref: ScalarField::Uniform(ra) }
                = other_boundary.conditions()[region]
            else { return Err(invalid("radiation binding lost its combined approximate Robin row")); };
            points[i] = [h, r, bound.coefficient_mean_slopes()[i], bound.local_reference_slopes()[i]];
            other[i] = [ha, ra];
        }
        let total_feedback = bound.reference_feedback()
            .ok_or_else(|| invalid("radiation binding lost its original reference-feedback map"))?;
        RobinResponse::compare_mean_robin_goal_with_reference_feedback_at(
            cx, ConductionProblem { boundary: bound.boundary(), ..problem }, interfaces,
            linear, reference, approximate, regions, &points[..m], &other[..m],
            total_feedback, weights, max_feedback_entries,
        )
    }
}

#[cfg(test)]
mod tests;
