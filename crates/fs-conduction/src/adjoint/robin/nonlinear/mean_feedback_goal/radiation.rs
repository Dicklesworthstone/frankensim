//! Two-field comparison through the producer's card-backed radiative secant.
//! Convection references may depend on other walls; radiation never heats air.

use super::{ConductionError, ConductionProblem, Cx, DiscreteGoalComparison,
    LinearConfig, RobinResponse, ThermalInterfaces, checked, invalid, poll, vector};
use crate::radiation::{AmbientRadiationGradient, AmbientRadiationPatch};
use std::collections::BTreeSet;

impl AmbientRadiationGradient {
    /// Compare a linear nodal-temperature goal through convection and radiation.
    ///
    /// `problem.boundary` is the ORIGINAL convective partition. In `regions`
    /// order, `reference_convection` supplies [h_c, r_c, dh_c/dmean, dr_c/dmean]
    /// at the reference field and `approximate_convection` supplies [h_c, r_c]
    /// at the approximate field. Units are W/(m^2 K), K, W/(m^2 K^2), and K/K.
    /// The additional row-major m-by-m `reference_feedback` maps all wall means
    /// to convective references. Do not repeat its diagonal in the local slope.
    /// Supply zeros for fixed references. The caller owns the convection law,
    /// evaluates BOTH physical states, and supplies its exact derivative.
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
                return Err(invalid("every radiation patch needs one unique selected convective region"));
            }
        }
        let mut points = reference_convection.to_vec();
        let mut other = approximate_convection.to_vec();
        let mut feedback = reference_feedback.to_vec();
        for (i, &name) in regions.iter().enumerate() {
            poll(cx, i)?;
            vector(cx, &points[i], 4)?;
            vector(cx, &other[i], 2)?;
            let [hc, tc, dhc, dtc] = points[i];
            let [ha, ta] = other[i];
            if hc <= 0.0 || ha <= 0.0 || tc <= 0.0 || ta <= 0.0 {
                return Err(invalid("radiative goal needs positive convection coefficients and absolute references"));
            }
            let Some(patch) = patches.iter().find(|p| p.region() == name) else { continue; };
            let region = problem.boundary.region_names().iter().position(|r| r == name)
                .ok_or_else(|| invalid("radiative-goal region is absent from the original partition"))?;
            let mut area = 0.0;
            let mut sums = [0.0; 2];
            for (slot, face) in problem.mesh.boundary().iter().enumerate() {
                if slot % 512 == 0 { poll(cx, slot)?; }
                if problem.boundary.region_for(slot) != Some(region) { continue; }
                area = checked(area + face.area)?;
                for &v in &face.vertices {
                    sums[0] = checked(sums[0] + (face.area/3.0)*reference[v as usize])?;
                    sums[1] = checked(sums[1] + (face.area/3.0)*approximate[v as usize])?;
                }
            }
            if area <= 0.0 { return Err(invalid("radiative-goal patch has no positive trace area")); }
            let wall = checked(sums[0]/area)?;
            let other_wall = checked(sums[1]/area)?;
            let hr = patch.secant_coefficient_w_m2_k(wall)?;
            let hra = patch.secant_coefficient_w_m2_k(other_wall)?;
            let dhr = patch.secant_partials_w_m2_k2(wall)?[0];
            let reservoir = patch.ambient_temperature_k();
            let h = checked(hc + hr)?;
            let other_h = checked(ha + hra)?;
            let r = checked((hc/h)*tc + (hr/h)*reservoir)?;
            points[i] = [h, r, checked(dhc + dhr)?, checked(
                (dhc/h)*(tc-r) + (dhr/h)*(reservoir-r) + (hc/h)*dtc)?];
            other[i] = [other_h, checked((ha/other_h)*ta + (hra/other_h)*reservoir)?];
            for value in &mut feedback[i*m..(i+1)*m] {
                *value = checked((hc/h)*(*value))?;
            }
        }
        RobinResponse::compare_mean_robin_goal_with_reference_feedback_at(
            cx, problem, interfaces, linear, reference, approximate, regions,
            &points, &other, &feedback, weights, max_feedback_entries,
        )
    }
}

#[cfg(test)]
mod tests;
