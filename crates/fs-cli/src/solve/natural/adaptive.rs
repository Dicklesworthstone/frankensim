//! Actual Churchill-Chu laws in an enriched thermal-goal comparison.
//! Own the same unmodified card, multiplier and pressure as the primal.

use super::{Coefficient, CorrelationId, NaturalLaw, SolveRefusal, coefficient, conduction_error};
use fs_conduction::adjoint::{DiscreteGoalComparison, RobinResponse};
use fs_conduction::{ConductionError, ConductionProblem, LinearConfig, ScalarField, ThermalBc,
    ThermalInterfaces};
use fs_exec::Cx;

fn invalid(what: impl Into<String>) -> SolveRefusal {
    conduction_error("cli-solve-conduction-natural-adaptive", what,
        "use admitted Churchill-Chu vertical-plate laws and fund the complete enriched comparison; no frozen-coefficient fallback is available")
}
fn finite(value: f64) -> Result<f64, SolveRefusal> {
    if value.is_finite() { Ok(value) } else { Err(invalid("nonfinite natural-goal arithmetic")) }
}
fn poll(cx: &Cx<'_>) -> Result<(), SolveRefusal> {
    cx.checkpoint().map_err(|_| conduction_error("cli-solve-cancelled",
        "natural-convection goal comparison cancelled", "resume the run"))
}
fn lower(error: ConductionError) -> SolveRefusal {
    if matches!(error, ConductionError::Cancelled { .. }) {
        conduction_error("cli-solve-cancelled", "natural-convection goal comparison cancelled",
            "resume the run")
    } else {
        invalid(format!("complete natural-convection goal comparison refused: {error}"))
    }
}

/// A rung retains its actual constitutive inputs, not a rounded receipt or a
/// Robin coefficient evaluated at an earlier fixed-point iterate.
pub(in crate::solve) struct NaturalGoal {
    laws: Vec<NaturalLaw>,
    pressure_pa: f64,
    max_feedback_entries: usize,
}

impl NaturalGoal {
    pub(in crate::solve) fn new(
        laws: Vec<NaturalLaw>, pressure_pa: f64, memory_bytes: u64,
    ) -> Result<Self, SolveRefusal> {
        if laws.is_empty() || laws.len() > 64 || !(pressure_pa.is_finite() && pressure_pa > 0.0)
            || laws.iter().any(|law| law.card != CorrelationId::ChurchillChuVerticalPlate)
        {
            return Err(invalid("a natural goal needs one to 64 differentiated vertical-plate laws and positive pressure"));
        }
        let mut names = std::collections::BTreeSet::new();
        for law in &laws {
            if !names.insert(law.target.as_str()) {
                return Err(invalid("duplicate natural-goal target"));
            }
        }
        Ok(Self { laws, pressure_pa,
            max_feedback_entries: usize::try_from(memory_bytes / 32).unwrap_or(usize::MAX) })
    }

    /// Check both fields against the real law. Smooth k(T), prescribed values,
    /// fixed contact, consistent face mass and the total transpose remain in
    /// the conduction owner. The finite nonlinear remainder stays observable.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::solve) fn compare(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, linear: LinearConfig,
        reference: &[f64], approximate: &[f64], weights: &[f64],
    ) -> Result<DiscreteGoalComparison, SolveRefusal> {
        poll(cx)?;
        let n = problem.mesh.vertex_count();
        if n.checked_mul(2).and_then(|v| v.checked_mul(self.laws.len()))
            .is_none_or(|v| v > self.max_feedback_entries)
        {
            return Err(invalid("natural-goal feedback exceeds the declared memory allowance"));
        }
        for field in [reference, approximate, weights] {
            if field.len() != n { return Err(invalid("natural-goal field length differs from the mesh")); }
            for (i, value) in field.iter().enumerate() {
                if i % 512 == 0 { poll(cx)?; }
                if !value.is_finite() { return Err(invalid("nonfinite natural-goal field")); }
            }
        }
        let mut names = Vec::with_capacity(self.laws.len());
        let mut points = Vec::with_capacity(self.laws.len());
        let mut approximate_values = Vec::with_capacity(self.laws.len());
        for law in &self.laws {
            poll(cx)?;
            let region = problem.boundary.region_names().iter().position(|name| name == &law.target)
                .ok_or_else(|| invalid("natural-goal target has no retained boundary region"))?;
            let ThermalBc::Robin { htc: ScalarField::Uniform(_), t_ref: ScalarField::Uniform(ambient) }
                = &problem.boundary.conditions()[region]
            else { return Err(invalid("natural goals require their original uniform Robin rows")); };
            if *ambient != law.ambient_k { return Err(invalid("natural-goal reference changed the declared ambient")); }
            let mut area = 0.0;
            let mut integral = 0.0;
            let mut approximate_integral = 0.0;
            for (slot, face) in problem.mesh.boundary().iter().enumerate() {
                if slot % 512 == 0 { poll(cx)?; }
                if problem.boundary.region_for(slot) != Some(region) { continue; }
                area = finite(area + face.area)?;
                for &v in &face.vertices {
                    integral = finite(integral + face.area / 3.0 * reference[v as usize])?;
                    approximate_integral = finite(approximate_integral
                        + face.area / 3.0 * approximate[v as usize])?;
                }
            }
            if area <= 0.0 { return Err(invalid("natural-goal target has no positive area")); }
            let mean = finite(integral / area)?;
            let approximate_mean = finite(approximate_integral / area)?;
            let point = coefficient(law, finite(mean - law.ambient_k)?, self.pressure_pa)?;
            let other = coefficient(law, finite(approximate_mean - law.ambient_k)?, self.pressure_pa)?;
            let slope = mean_slope(law, mean, point)?;
            names.push(law.target.as_str());
            points.push([point.htc_w_m2_k, law.ambient_k, slope, 0.0]);
            approximate_values.push([other.htc_w_m2_k, law.ambient_k]);
        }
        RobinResponse::compare_mean_robin_goal_at(cx, problem, interfaces, linear,
            reference, approximate, &names, &points, &approximate_values, weights,
            self.max_feedback_entries).map_err(lower)
    }
}

/// At fixed pressure, Ra ~ |Tw-Ta|/Tfilm^3 and Nu=(0.825+x)^2,
/// x ~ Ra^(1/6). Use the card's UNSCALED Nu and the APPLIED h so a model
/// discrepancy multiplier affects the derivative without changing the card.
/// The signed delta denominator is necessary on the cooled-wall branch.
fn mean_slope(law: &NaturalLaw, wall: f64, point: Coefficient) -> Result<f64, SolveRefusal> {
    let delta = finite(wall - law.ambient_k)?;
    let film = finite(law.ambient_k + 0.5 * delta)?;
    if law.card != CorrelationId::ChurchillChuVerticalPlate
        || delta == 0.0 || wall <= 0.0 || film <= 0.0
        || !(point.nusselt.is_finite() && point.nusselt > 0.0)
        || !(point.htc_w_m2_k.is_finite() && point.htc_w_m2_k > 0.0)
    {
        return Err(invalid("natural-goal tangent lies outside its admitted smooth physical branch"));
    }
    let exponent = finite((1.0 - 0.825 / fs_math::det::sqrt(point.nusselt)) / 3.0)?;
    if !(exponent > 0.0 && exponent < 1.0 / 3.0) {
        return Err(invalid("natural-goal Nusselt value has no admitted Churchill-Chu slope"));
    }
    finite(point.htc_w_m2_k * exponent * (1.0 / delta - 1.5 / film))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn law(scale: f64) -> NaturalLaw {
        NaturalLaw { target: "wall".into(), length_m: 0.06, ambient_k: 300.0,
            card: CorrelationId::ChurchillChuVerticalPlate, htc_multiplier: scale }
    }
    #[test]
    fn nominal_adjoint_natural_goal_slopes_follow_the_actual_scaled_card_on_both_branches() {
        for scale in [0.75, 1.0, 1.25] { for pressure in [80000.0, 101325.0] {
            let input = law(scale);
            for delta in [-80.0, -25.0, -1.0, 1.0, 25.0, 80.0] {
                let point = coefficient(&input, delta, pressure).unwrap();
                let slope = mean_slope(&input, input.ambient_k + delta, point).unwrap();
                let step = 1e-3;
                let expected = (coefficient(&input, delta + step, pressure).unwrap().htc_w_m2_k
                    - coefficient(&input, delta - step, pressure).unwrap().htc_w_m2_k) / (2.0 * step);
                assert!((slope - expected).abs() < 2e-6 * expected.abs().max(1e-3),
                    "scale {scale}, delta {delta}: {slope:e} vs {expected:e}");
                let nominal = law(1.0);
                let unscaled = mean_slope(&nominal, nominal.ambient_k + delta,
                    coefficient(&nominal, delta, pressure).unwrap()).unwrap();
                assert!((slope - scale * unscaled).abs() < 1e-13 * slope.abs().max(1.0));
            }
        } }
    }

    #[test]
    fn nominal_adjoint_natural_goal_rejects_invalid_law_sets_and_singular_tangents() {
        assert!(NaturalGoal::new(Vec::new(), 101325.0, 1 << 20).is_err());
        assert!(NaturalGoal::new(vec![law(1.0), law(1.0)], 101325.0, 1 << 20).is_err());
        assert!(NaturalGoal::new(vec![law(1.0)], f64::NAN, 1 << 20).is_err());
        let point = coefficient(&law(1.0), 10.0, 101325.0).unwrap();
        for wall in [0.0, 300.0, f64::NAN] {
            assert!(mean_slope(&law(1.0), wall, point).is_err());
        }
    }
}
