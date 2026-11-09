//! Resolution convergence of the original nonlinear contact problem.
//!
//! Every rung closes the complete original material law before comparison.
//! Differences are numerical estimates, never physical heat or certified bounds.

use super::ContactResponseError;

/// Automatic refinement of the original contact quadrature.
///
/// A refinement doubles the subdivisions and repeats the complete constitutive
/// solve from the same externally supplied physical state. At least two solved
/// resolutions are compared. Every component uses its dimensional absolute
/// tolerance plus `relative_tolerance * max(abs(coarse), abs(fine))`.
///
/// Before each rung, the work owner reserves its COMPLETE worst-case nonlinear
/// solve: `maximum_iterations * solves_per_iteration` field trials and twice
/// `maximum_iterations * (bulk_sources + contact_points)` original-source queries.
/// The second query allowance covers relaxed-coordinate evaluation. Each rung's
/// iteration cap is narrowed to the remaining work before reservation. Admission
/// can therefore refuse a rung that might have converged more cheaply. Successful
/// reports retain actual work, and unused reservations are released after a rung.
/// External callbacks retain their own bounds for work hidden inside one trial.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactQuadratureRefinementControls {
    /// Maximum subdivision doublings after the initial solved resolution.
    pub maximum_refinements: usize,
    pub maximum_field_solves: usize,
    pub maximum_source_evaluations: usize,
    pub absolute_voltage_tolerance_v: f64,
    pub absolute_current_tolerance_a: f64,
    pub absolute_power_tolerance_w: f64,
    pub absolute_current_density_tolerance_a_m2: f64,
    pub relative_tolerance: f64,
}

impl ContactQuadratureRefinementControls {
    pub fn validate(self) -> Result<(), ContactResponseError> {
        let absolute = [
            self.absolute_voltage_tolerance_v,
            self.absolute_current_tolerance_a,
            self.absolute_power_tolerance_w,
            self.absolute_current_density_tolerance_a_m2,
        ];
        if self.maximum_refinements == 0
            || self.maximum_refinements >= usize::BITS as usize
            || self.maximum_field_solves == 0
            || self.maximum_source_evaluations == 0
            || !self.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_tolerance)
            || absolute.iter().any(|v| {
                !v.is_finite() || *v < 0.0 || (*v == 0.0 && self.relative_tolerance == 0.0)
            })
        {
            return Err(ContactResponseError::Invalid(
                "positive bounded contact refinement and dimensional tolerances required",
            ));
        }
        Ok(())
    }

    pub(crate) fn compare(
        self,
        coarse: &ContactQuadratureSnapshot,
        fine: &ContactQuadratureSnapshot,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ContactQuadratureDifferences, ContactResponseError> {
        self.validate()?;
        let mut result = ContactQuadratureDifferences {
            converged: true,
            ..ContactQuadratureDifferences::default()
        };
        for (a, b, absolute, maximum) in [
            (
                &coarse.voltage_differences_v,
                &fine.voltage_differences_v,
                self.absolute_voltage_tolerance_v,
                &mut result.maximum_voltage_difference_v,
            ),
            (
                &coarse.current_a,
                &fine.current_a,
                self.absolute_current_tolerance_a,
                &mut result.maximum_current_difference_a,
            ),
            (
                &coarse.power_w,
                &fine.power_w,
                self.absolute_power_tolerance_w,
                &mut result.maximum_power_difference_w,
            ),
            (
                &coarse.current_density_a_m2,
                &fine.current_density_a_m2,
                self.absolute_current_density_tolerance_a_m2,
                &mut result.maximum_current_density_difference_a_m2,
            ),
        ] {
            if a.is_empty() || a.len() != b.len() {
                return Err(ContactResponseError::Invalid(
                    "matching original contact refinement observations required",
                ));
            }
            for (&old, &new) in a.iter().zip(b) {
                if cancelled() {
                    return Err(ContactResponseError::Cancelled);
                }
                let difference = (new - old).abs();
                let tolerance = absolute + self.relative_tolerance * old.abs().max(new.abs());
                if !old.is_finite()
                    || !new.is_finite()
                    || !difference.is_finite()
                    || !tolerance.is_finite()
                {
                    return Err(ContactResponseError::Invalid(
                        "finite contact refinement observations and differences required",
                    ));
                }
                *maximum = maximum.max(difference);
                result.converged &= difference <= tolerance;
            }
        }
        Ok(result)
    }
}

/// Componentwise differences between the last two fully converged resolutions.
/// Gauge-invariant edge drops and original contact endpoint jumps supply volts;
/// individual contact currents, bulk/contact heat and actual current densities
/// are retained separately, so cancellation cannot hide an unresolved contact.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ContactQuadratureDifferences {
    pub maximum_voltage_difference_v: f64,
    pub maximum_current_difference_a: f64,
    pub maximum_power_difference_w: f64,
    pub maximum_current_density_difference_a_m2: f64,
    pub converged: bool,
}

/// Actual resolution and work of an accepted contact-refinement solve.
/// This is an observed resolution-convergence estimate, not an error certificate.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactQuadratureRefinementReport {
    pub initial_subdivisions: usize,
    pub final_subdivisions: usize,
    pub refinements: usize,
    pub field_solves: usize,
    pub source_evaluations: usize,
    pub differences: ContactQuadratureDifferences,
}

/// Geometry adapters supply observations in the SAME original ordering at
/// every resolution; quadrature-point arrays themselves change length.
pub(crate) struct ContactQuadratureSnapshot {
    pub voltage_differences_v: Vec<f64>,
    pub current_a: Vec<f64>,
    pub power_w: Vec<f64>,
    pub current_density_a_m2: Vec<f64>,
}

pub(crate) struct RefinementWork {
    controls: ContactQuadratureRefinementControls,
    field_solves: usize,
    source_evaluations: usize,
    reserved_fields: usize,
    reserved_sources: usize,
}

impl RefinementWork {
    pub(crate) fn new(
        controls: ContactQuadratureRefinementControls,
    ) -> Result<Self, ContactResponseError> {
        controls.validate()?;
        Ok(Self {
            controls,
            field_solves: 0,
            source_evaluations: 0,
            reserved_fields: 0,
            reserved_sources: 0,
        })
    }

    /// Admit a nonlinear iteration cap that fits both remaining global budgets.
    /// Source work uses a conservative complete relaxed-iteration reservation.
    pub(crate) fn iteration_allowance(
        &self,
        maximum_iterations: usize,
        solves_per_iteration: usize,
        bulk_sources: usize,
        total_contact_points: usize,
    ) -> Result<usize, ContactResponseError> {
        let sources_per_iteration = bulk_sources
            .checked_add(total_contact_points)
            .and_then(|count| count.checked_mul(2))
            .filter(|&count| count > 0);
        if maximum_iterations == 0 || solves_per_iteration == 0 || sources_per_iteration.is_none() {
            return Err(ContactResponseError::Invalid(
                "representable positive contact refinement iteration work required",
            ));
        }
        let sources_per_iteration = sources_per_iteration.unwrap();
        let mut iterations = maximum_iterations;
        for (work, completed, required, maximum) in [
            (
                "field solves",
                self.field_solves,
                solves_per_iteration,
                self.controls.maximum_field_solves,
            ),
            (
                "source evaluations",
                self.source_evaluations,
                sources_per_iteration,
                self.controls.maximum_source_evaluations,
            ),
        ] {
            let available = maximum - completed;
            if available < required {
                return Err(ContactResponseError::QuadratureBudget {
                    work,
                    completed,
                    required,
                    maximum,
                });
            }
            iterations = iterations.min(available / required);
        }
        Ok(iterations)
    }

    pub(crate) fn reserve(
        &mut self,
        maximum_iterations: usize,
        solves_per_iteration: usize,
        bulk_sources: usize,
        total_contact_points: usize,
    ) -> Result<(), ContactResponseError> {
        let fields = maximum_iterations.checked_mul(solves_per_iteration);
        let sources = bulk_sources
            .checked_add(total_contact_points)
            .and_then(|count| count.checked_mul(maximum_iterations))
            .and_then(|count| count.checked_mul(2));
        if maximum_iterations == 0
            || solves_per_iteration == 0
            || fields.is_none()
            || sources.is_none()
        {
            return Err(ContactResponseError::Invalid(
                "representable complete contact refinement work reservation required",
            ));
        }
        let (fields, sources) = (fields.unwrap(), sources.unwrap());
        for (work, completed, required, maximum) in [
            (
                "field solves",
                self.field_solves,
                fields,
                self.controls.maximum_field_solves,
            ),
            (
                "source evaluations",
                self.source_evaluations,
                sources,
                self.controls.maximum_source_evaluations,
            ),
        ] {
            if completed
                .checked_add(required)
                .is_none_or(|total| total > maximum)
            {
                return Err(ContactResponseError::QuadratureBudget {
                    work,
                    completed,
                    required,
                    maximum,
                });
            }
        }
        self.reserved_fields = fields;
        self.reserved_sources = sources;
        Ok(())
    }

    pub(crate) fn record(
        &mut self,
        field_solves: usize,
        source_evaluations: usize,
    ) -> Result<(), ContactResponseError> {
        if field_solves > self.reserved_fields || source_evaluations > self.reserved_sources {
            return Err(ContactResponseError::Invalid(
                "native contact refinement work exceeded its admitted reservation",
            ));
        }
        // Each addition is bounded by the checked reservation above.
        self.field_solves += field_solves;
        self.source_evaluations += source_evaluations;
        self.reserved_fields = 0;
        self.reserved_sources = 0;
        Ok(())
    }

    pub(crate) fn report(
        &self,
        initial_subdivisions: usize,
        final_subdivisions: usize,
        refinements: usize,
        differences: ContactQuadratureDifferences,
    ) -> ContactQuadratureRefinementReport {
        ContactQuadratureRefinementReport {
            initial_subdivisions,
            final_subdivisions,
            refinements,
            field_solves: self.field_solves,
            source_evaluations: self.source_evaluations,
            differences,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn controls() -> ContactQuadratureRefinementControls {
        ContactQuadratureRefinementControls {
            maximum_refinements: 3,
            maximum_field_solves: 40,
            maximum_source_evaluations: 2000,
            absolute_voltage_tolerance_v: 1e-6,
            absolute_current_tolerance_a: 1e-6,
            absolute_power_tolerance_w: 1e-6,
            absolute_current_density_tolerance_a_m2: 1e-6,
            relative_tolerance: 1e-6,
        }
    }

    fn observations(values: Vec<f64>) -> ContactQuadratureSnapshot {
        ContactQuadratureSnapshot {
            voltage_differences_v: values.clone(),
            current_a: values.clone(),
            power_w: values.clone(),
            current_density_a_m2: values,
        }
    }

    #[test]
    fn independent_contact_changes_cannot_cancel_in_the_acceptance_gate() {
        let a = observations(vec![1.0, 2.0]);
        let mut b = observations(vec![1.0, 2.0]);
        b.current_a = vec![1.01, 1.99];
        b.power_w = vec![1.01, 1.99];
        let diff = controls().compare(&a, &b, &mut || false).unwrap();
        assert!(!diff.converged);
        assert!((diff.maximum_current_difference_a - 0.01).abs() < 1e-15);
        assert!(controls().compare(&a, &a, &mut || false).unwrap().converged);
        b = observations(vec![1.0 + 1e-7, 2.0 - 1e-7]);
        assert!(controls().compare(&a, &b, &mut || false).unwrap().converged);
    }

    #[test]
    fn reservation_bounds_full_rungs_but_charges_only_actual_work() {
        let mut work = RefinementWork::new(controls()).unwrap();
        work.reserve(20, 1, 2, 5).unwrap();
        work.record(3, 20).unwrap();
        assert_eq!(work.iteration_allowance(40, 1, 2, 9).unwrap(), 37);
        work.reserve(20, 1, 2, 9).unwrap();
        work.record(4, 42).unwrap();
        let report = work.report(1, 2, 1, ContactQuadratureDifferences::default());
        assert_eq!((report.field_solves, report.source_evaluations), (7, 62));
        assert!(matches!(
            work.reserve(34, 1, 2, 9),
            Err(ContactResponseError::QuadratureBudget {
                work: "field solves",
                ..
            })
        ));
        assert!(matches!(
            work.reserve(1, 1, 2, 1000),
            Err(ContactResponseError::QuadratureBudget {
                work: "source evaluations",
                ..
            })
        ));
        assert!(work.reserve(usize::MAX, 2, 2, 9).is_err());
    }

    #[test]
    fn invalid_observations_and_late_cancellation_refuse_without_acceptance() {
        let a = observations(vec![1.0, 2.0]);
        let mut b = observations(vec![1.0, 2.0]);
        b.power_w[1] = f64::NAN;
        assert!(controls().compare(&a, &b, &mut || false).is_err());
        let mut polls = 0;
        assert!(matches!(
            controls().compare(&a, &a, &mut || {
                polls += 1;
                polls == 5
            }),
            Err(ContactResponseError::Cancelled)
        ));
        assert!(controls().compare(&a, &a, &mut || false).unwrap().converged);
        let mut bad = controls();
        bad.maximum_refinements = 0;
        assert!(bad.validate().is_err());
    }
}
