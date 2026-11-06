//! Two physical ambient-radiation states and the complete live-contact dual.

use super::{TemperatureContacts, compare_goal_fields, finite, invalid, poll};
use crate::adjoint::DiscreteGoalComparison;
use crate::radiation::{AmbientRadiationLinearization, AmbientRadiationPatch};
use crate::{ConductionError, ConductionProblem, LinearConfig, ScalarField, ThermalBc};
use fs_exec::Cx;

impl TemperatureContacts<'_, '_> {
    /// Compare live-contact fields under the original mean ambient-radiation law.
    ///
    /// `problem.boundary` is the ORIGINAL convection/zero-flux partition. In
    /// `regions` order, reference points are `[h_c, r_c, dh_c/dmean, dr_c/dmean]`
    /// and approximate values are `[h_c, r_c]`, evaluated by their convection
    /// owner at each field. The row-major `reference_feedback[i*m+j]` adds
    /// dr_c_i/dmean_j; do not repeat its diagonal in the local reference slope.
    /// Every radiating and mean-dependent convection trace must be selected.
    ///
    /// A radiation-only trace remains an original uniform zero-flux region:
    /// its convection coefficients, slopes and feedback ROW must be zero.
    /// Its finite positive reference input is an inert placeholder; the
    /// radiation patch retains the actual black-reservoir temperature. No
    /// convection is manufactured on that trace or counted as air heat.
    ///
    /// Both residuals re-query the original contact source and evaluate the
    /// radiative secant at their own area-mean wall temperatures. The reference
    /// uses the complete bulk/contact tangent plus coefficient, weighted-
    /// reference and cross-wall feedback from the original radiation owner.
    /// Convection feedback is scaled by h_c/(h_c+h_rad), exactly once. The
    /// reference must pass the original watt residual and heat gates, even
    /// for a zero goal; no additional primal solve or tolerance change occurs.
    ///
    /// Only the reference needs a unique material/contact/emissivity tangent.
    /// The approximation may lie on a supported emissivity interval endpoint
    /// or cross a contact source segment: it needs physical values only. The
    /// observed nonlinear remainder is retained with the signed nodal terms.
    /// This is the mean-temperature secant model and a discrete comparison,
    /// not pointwise T(x)^4 quadrature, enclosure exchange or a continuum bound.
    /// At most 64 regions/patches, four nodal factors per region and the m*m
    /// feedback map must fit their explicit allowance, alongside contact work.
    #[allow(clippy::too_many_arguments)]
    pub fn compare_steady_goal_with_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        reference_temperature: &[f64],
        approximate_temperature: &[f64],
        linear: LinearConfig,
        residual_tolerance_w: f64,
        energy_tolerance_w: f64,
        regions: &[&str],
        reference_points: &[[f64; 4]],
        approximate_values: &[[f64; 2]],
        reference_feedback: &[f64],
        patches: &[AmbientRadiationPatch],
        full_nodal_weights: &[f64],
        max_feedback_entries: usize,
        max_face_pairs: usize,
    ) -> Result<DiscreteGoalComparison, ConductionError> {
        self.admit_goal_comparison(
            cx,
            problem,
            reference_temperature,
            approximate_temperature,
            linear,
            residual_tolerance_w,
            energy_tolerance_w,
            full_nodal_weights,
            max_face_pairs,
        )?;
        let m = regions.len();
        let entries = problem
            .mesh
            .vertex_count()
            .checked_mul(4)
            .and_then(|n| n.checked_mul(m));
        if m == 0
            || m > 64
            || patches.is_empty()
            || patches.len() > 64
            || reference_points.len() != m
            || approximate_values.len() != m
            || reference_feedback.len() != m * m
            || m * m > max_feedback_entries
            || entries.is_none_or(|n| {
                n > max_feedback_entries
                    || n.checked_mul(8)
                        .is_none_or(|bytes| bytes > isize::MAX as usize)
            })
        {
            return Err(invalid(
                "<steady-contact-radiation-goal>",
                "matching radiation/convection rows, a complete reference map and bounded factor storage are required",
            ));
        }
        for (i, patch) in patches.iter().enumerate() {
            poll(cx, i)?;
            if !regions.contains(&patch.region())
                || patches[..i]
                    .iter()
                    .any(|other| other.region() == patch.region())
            {
                return Err(invalid(
                    "<steady-contact-radiation-goal>",
                    "each radiation patch must occur once in the selected original regions",
                ));
            }
        }
        for (i, &value) in reference_feedback.iter().enumerate() {
            if i % 256 == 0 {
                poll(cx, i)?;
            }
            finite(value)?;
        }
        let mut reference_rows = [(0, 0.0, 0.0); 64];
        let mut approximate_rows = [(0, 0.0, 0.0); 64];
        let mut slopes = [0.0; 64];
        let mut diagonal = [0.0; 64];
        let mut convection_rows = 0;
        for (i, &name) in regions.iter().enumerate() {
            poll(cx, i)?;
            for &value in reference_points[i].iter().chain(&approximate_values[i]) {
                finite(value)?;
            }
            let region = problem
                .boundary
                .region_names()
                .iter()
                .position(|r| r == name)
                .ok_or_else(|| {
                    invalid(
                        "<steady-contact-radiation-goal>",
                        "unknown original thermal region",
                    )
                })?;
            let [h, r, dh, dr] = reference_points[i];
            let [ha, ra] = approximate_values[i];
            if regions[..i].contains(&name) || r <= 0.0 || ra <= 0.0 {
                return Err(invalid(
                    "<steady-contact-radiation-goal>",
                    "selected regions must be unique with positive absolute reference values",
                ));
            }
            match problem.boundary.conditions()[region] {
                ThermalBc::Robin {
                    htc: ScalarField::Uniform(_),
                    t_ref: ScalarField::Uniform(_),
                } if h > 0.0 && ha > 0.0 => {
                    reference_rows[convection_rows] = (region, h, r);
                    approximate_rows[convection_rows] = (region, ha, ra);
                    convection_rows += 1;
                }
                ThermalBc::Neumann {
                    outward_flux: ScalarField::Uniform(flux),
                } if flux == 0.0
                    && h == 0.0
                    && ha == 0.0
                    && dh == 0.0
                    && dr == 0.0
                    && patches.iter().any(|patch| patch.region() == name)
                    && reference_feedback[i * m..(i + 1) * m]
                        .iter()
                        .all(|&v| v == 0.0) => {}
                _ => {
                    return Err(invalid(
                        "<steady-contact-radiation-goal>",
                        "retain uniform convection or radiation-only zero-flux ownership with no invented convective feedback",
                    ));
                }
            }
            slopes[i] = dh;
            diagonal[i] = finite(reference_feedback[i * m + i] + dr)?;
        }
        let mut feedback = reference_feedback.to_vec();
        for i in 0..m {
            feedback[i * m + i] = diagonal[i];
        }
        let reference_convection = problem
            .boundary
            .with_uniform_robin_replacements(&reference_rows[..convection_rows])?;
        let approximate_convection = problem
            .boundary
            .with_uniform_robin_replacements(&approximate_rows[..convection_rows])?;
        let bound = AmbientRadiationLinearization::bind(
            cx,
            ConductionProblem {
                boundary: &reference_convection,
                ..problem
            },
            reference_temperature,
            regions,
            &slopes[..m],
            Some(&feedback),
            patches,
            max_feedback_entries,
        )?;
        let reference_problem = ConductionProblem {
            boundary: bound.boundary(),
            ..problem
        };
        let endpoint = self.linearize_steady_endpoint(
            cx,
            reference_problem,
            reference_temperature,
            linear,
            residual_tolerance_w,
            energy_tolerance_w,
            regions,
            max_face_pairs,
        )?;
        let approximate_boundary = AmbientRadiationLinearization::physical_boundary_at(
            cx,
            ConductionProblem {
                boundary: &approximate_convection,
                ..problem
            },
            approximate_temperature,
            regions,
            patches,
            max_feedback_entries,
        )?;
        let approximate_contacts =
            self.evaluate(cx, &approximate_boundary, approximate_temperature)?;
        let total_feedback = bound.reference_feedback().ok_or_else(|| {
            invalid(
                "<steady-contact-radiation-goal>",
                "radiation binding lost its supplied reference map",
            )
        })?;
        compare_goal_fields(
            cx,
            reference_problem,
            &approximate_boundary,
            approximate_temperature,
            full_nodal_weights,
            &endpoint,
            &approximate_contacts,
            || {
                endpoint.pullback_endpoint_with_mean_feedback(
                    cx,
                    bound.coefficient_mean_slopes(),
                    bound.local_reference_slopes(),
                    total_feedback,
                    full_nodal_weights,
                    max_feedback_entries,
                )
            },
        )
    }
}
