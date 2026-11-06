//! Nonlinear free flight composed from the native periodic slender-jet owner.
//!
//! The authored cell-centred radii and staggered laboratory axial velocities
//! remain the native initial condition. No inlet radius, flow rate, growth mode,
//! or velocity perturbation is inferred. Transport uses the initial conserved
//! volume-weighted speed U0, and a requested flight length L admits at most L/U0
//! seconds of the native temporal evolution.
//! This is a bulk temporal-to-distance mapping; individual fluid parcels can
//! have different trajectories in the nonuniform evolving velocity field.
//!
//! This is a periodic constant-property Newtonian segment, not a nozzle inlet
//! boundary-value problem. Neglected gravity-induced stretching, gas stresses,
//! and end effects remain assumptions of the native slender approximation.
//! An admitted neck or slope limit is a stopping condition, never a resolved
//! pinch-off, detached droplet, or measured breakup certificate.
//!
//! The owner keeps the original native or material-bound jet. Native stepping
//! supplies all mass, momentum, curvature and viscous dynamics; this module
//! only sets the flight duration and retains the accepted state and stop reason.

use fs_exec::CancelGate;
use fs_lbm::jet::nonlinear::{
    NonlinearJetError, PeriodicSlenderJet, SlenderJetAdvance, SlenderJetDiagnostics,
    SlenderJetEvolution, SlenderJetOutcome, SlenderJetStepBudget, SlenderJetViolation,
};

/// Progress of this flight, distinct from completion of one native time chunk.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NonlinearFreeJetOutcome {
    /// An admitted chunk completed before the requested flight time.
    InFlight,
    /// The native clock reached the requested residence-time target.
    FlightComplete,
    /// The native model stopped before accepting an out-of-domain substep.
    DomainLimit(SlenderJetViolation),
}

/// Accepted flight progress and the original native physical result.
#[derive(Debug)]
pub struct NonlinearFreeJetReport {
    /// Completion or the original physical approximation-limit stop.
    pub outcome: NonlinearFreeJetOutcome,
    /// Laboratory transport velocity from the initial momentum and volume [m/s].
    pub mean_velocity_m_s: f64,
    /// Authored downstream travel distance [m].
    pub requested_flight_length_m: f64,
    /// Nominal residence time L/U0 [s].
    pub residence_time_s: f64,
    /// Native accepted time since this flight started [s].
    pub elapsed_time_s: f64,
    /// U0 times accepted elapsed time [m], including native clock rounding.
    pub travelled_distance_m: f64,
    /// Physical diagnostics of the same accepted, evolving native state.
    pub state: SlenderJetDiagnostics,
    /// The actual advance, with its before/after state, accepted and attempted
    /// substeps, viscous/upwind dissipation and energy residual. None means this
    /// call only observed the accepted state without a native time advance.
    pub native: Option<SlenderJetAdvance>,
}

/// A native or source-bound periodic jet travelling through a specified flight.
///
/// Both admitted initial profiles are retained for comparison with the evolving
/// state. Velocities are laboratory-frame face velocities, not unspecified mode
/// amplitudes or velocity deviations. No mean velocity is silently subtracted.
/// The original J remains owned here, so a material binding retains its exact
/// immutable source throughout every accepted advance.
pub struct NonlinearFreeJet<J = PeriodicSlenderJet> {
    jet: J,
    initial_radii_m: Vec<f64>,
    initial_velocities_m_s: Vec<f64>,
    initial_time_s: f64,
    target_time_s: f64,
    flight_length_m: f64,
    mean_velocity_m_s: f64,
    residence_time_s: f64,
    stop: Option<SlenderJetViolation>,
}

impl<J: SlenderJetEvolution> NonlinearFreeJet<J> {
    /// Start a flight from the jet's current admitted physical state.
    ///
    /// U0 is the native mean axial velocity, equivalent to initial momentum
    /// divided by density times volume. Internal periodic capillary/viscous
    /// fluxes conserve that transport momentum; no second flow-rate input is
    /// accepted. A new flight can start from an already evolved native state.
    ///
    /// # Errors
    /// Refuses nonfinite/negative length, nonpositive downstream transport, or
    /// an unrepresentable residence time, native clock target or travel distance.
    pub fn new(jet: J, flight_length_m: f64) -> Result<Self, NonlinearJetError> {
        if !flight_length_m.is_finite() || flight_length_m < 0.0 {
            return Err(NonlinearJetError::InvalidInput("jet flight length"));
        }
        let state = jet.state();
        let mean_velocity_m_s = state.diagnostics().mean_velocity_m_s;
        if !mean_velocity_m_s.is_finite() || mean_velocity_m_s <= 0.0 {
            return Err(NonlinearJetError::InvalidInput(
                "a downstream free flight requires positive initial mean axial velocity",
            ));
        }
        let residence_time_s = flight_length_m / mean_velocity_m_s;
        let initial_time_s = state.time_s();
        let target_time_s = initial_time_s + residence_time_s;
        let represented_distance_m = mean_velocity_m_s * (target_time_s - initial_time_s);
        if !residence_time_s.is_finite()
            || !target_time_s.is_finite()
            || !represented_distance_m.is_finite()
            || (flight_length_m > 0.0
                && (residence_time_s == 0.0
                    || target_time_s <= initial_time_s
                    || represented_distance_m == 0.0))
        {
            return Err(NonlinearJetError::NonRepresentable);
        }
        // These are snapshots of the admitted native values, including its
        // h^2 and staggered-flux representation, not a regenerated eigenmode.
        let initial_radii_m = state.radii_m().collect();
        let initial_velocities_m_s = state.velocities_m_s().collect();
        Ok(Self {
            jet,
            initial_radii_m,
            initial_velocities_m_s,
            initial_time_s,
            target_time_s,
            flight_length_m,
            mean_velocity_m_s,
            residence_time_s,
            stop: None,
        })
    }

    /// Original native or material-bound owner, with its current accepted state.
    #[must_use]
    pub const fn jet(&self) -> &J {
        &self.jet
    }

    /// Admitted cell-centred radii at the start of this flight [m].
    #[must_use]
    pub fn initial_radii_m(&self) -> &[f64] {
        &self.initial_radii_m
    }

    /// Admitted laboratory face velocities at the start of this flight [m/s].
    #[must_use]
    pub fn initial_velocities_m_s(&self) -> &[f64] {
        &self.initial_velocities_m_s
    }

    /// Remaining native time until the residence target [s].
    ///
    /// A physical stop can leave a positive remainder: it is not converted into
    /// completion by shortening the authored flight.
    #[must_use]
    pub fn remaining_time_s(&self) -> f64 {
        (self.target_time_s - self.jet.state().time_s()).max(0.0)
    }

    /// Observe the accepted state without advancing it.
    #[must_use]
    pub fn report(&self) -> NonlinearFreeJetReport {
        self.make_report(None)
    }

    /// Advance at most one authored time chunk, capped by the remaining flight.
    ///
    /// Pass remaining_time_s() to request the whole remaining flight, or a
    /// shorter duration to inspect/resume actual nonlinear evolution. Native
    /// step-size and work budgets remain explicit. An ordinary chunk completing
    /// before the residence target returns InFlight.
    ///
    /// Native domain limits commit only the last admitted substep and terminate
    /// this flight. Repeated calls then observe the same stopped state. A zero-
    /// length, zero-duration or completed flight performs no native step.
    ///
    /// # Errors
    /// The chunk duration must be finite and nonnegative. Cancellation, numerical
    /// failure and exhausted substep budgets retain the entry state of this
    /// entire call, as guaranteed by the native evolution owner.
    pub fn advance(
        &mut self,
        maximum_duration_s: f64,
        budget: SlenderJetStepBudget,
        gate: &CancelGate,
    ) -> Result<NonlinearFreeJetReport, NonlinearJetError> {
        if !maximum_duration_s.is_finite() || maximum_duration_s < 0.0 {
            return Err(NonlinearJetError::InvalidInput(
                "free-flight chunk duration must be finite and nonnegative",
            ));
        }
        if maximum_duration_s == 0.0 || self.stop.is_some() || self.remaining_time_s() == 0.0 {
            return Ok(self.report());
        }
        let duration_s = maximum_duration_s.min(self.remaining_time_s());
        let native = self.jet.advance(duration_s, budget, gate)?;
        if let SlenderJetOutcome::DomainLimit(violation) = native.outcome {
            self.stop = Some(violation);
        }
        Ok(self.make_report(Some(native)))
    }

    /// Recover the original owner with its final accepted state and source.
    #[must_use]
    pub fn into_jet(self) -> J {
        self.jet
    }

    fn make_report(&self, native: Option<SlenderJetAdvance>) -> NonlinearFreeJetReport {
        let elapsed_time_s = self.jet.state().time_s() - self.initial_time_s;
        let state = native
            .as_ref()
            .map_or_else(|| self.jet.state().diagnostics(), |advance| advance.after);
        let outcome = if let Some(violation) = self.stop {
            NonlinearFreeJetOutcome::DomainLimit(violation)
        } else if self.remaining_time_s() == 0.0 {
            NonlinearFreeJetOutcome::FlightComplete
        } else {
            NonlinearFreeJetOutcome::InFlight
        };
        NonlinearFreeJetReport {
            outcome,
            mean_velocity_m_s: self.mean_velocity_m_s,
            requested_flight_length_m: self.flight_length_m,
            residence_time_s: self.residence_time_s,
            elapsed_time_s,
            travelled_distance_m: self.mean_velocity_m_s * elapsed_time_s,
            state,
            native,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_lbm::jet::nonlinear::{SlenderJetLimitKind, SlenderJetLimits, SlenderJetParameters};

    fn budget() -> SlenderJetStepBudget {
        SlenderJetStepBudget {
            maximum_step_s: 0.01,
            cfl: 0.2,
            maximum_substeps: 20_000,
        }
    }

    fn seeded_jet(mean_velocity_m_s: f64, minimum_radius_m: f64) -> PeriodicSlenderJet {
        let cells = 128;
        let k = 0.3;
        let length_m = core::f64::consts::TAU / k;
        let radii_m: Vec<_> = (0..cells)
            .map(|i| 1.0 + 0.01 * (core::f64::consts::TAU * i as f64 / cells as f64).cos())
            .collect();
        let velocities_m_s = vec![mean_velocity_m_s; cells];
        PeriodicSlenderJet::new(
            SlenderJetParameters {
                density_kg_m3: 1.0,
                surface_tension_n_m: 1.0,
                kinematic_viscosity_m2_s: 0.01,
            },
            length_m,
            &radii_m,
            &velocities_m_s,
            SlenderJetLimits {
                minimum_radius_m,
                maximum_abs_slope: 0.2,
                minimum_radius_cells: 2.0,
            },
            &CancelGate::new_clock_free(),
        )
        .unwrap()
    }

    fn near(actual: f64, expected: f64, relative: f64) {
        assert!(
            (actual - expected).abs() <= relative * expected.abs().max(1e-15),
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn g1_nonlinear_flight_evolves_the_original_shape_and_preserves_inventory() {
        let gate = CancelGate::new_clock_free();
        let jet = seeded_jet(0.05, 0.5);
        let before = jet.diagnostics();
        let mut flight = NonlinearFreeJet::new(jet, 0.05).unwrap();
        let initial_radii = flight.initial_radii_m().to_vec();
        let initial_velocities = flight.initial_velocities_m_s().to_vec();

        let first = flight.advance(0.4, budget(), &gate).unwrap();
        assert_eq!(first.outcome, NonlinearFreeJetOutcome::InFlight);
        near(first.elapsed_time_s, 0.4, 2e-14);
        near(first.travelled_distance_m, 0.02, 2e-14);
        assert!(first.native.as_ref().unwrap().accepted_substeps > 0);

        let report = flight.advance(10.0, budget(), &gate).unwrap();
        assert_eq!(report.outcome, NonlinearFreeJetOutcome::FlightComplete);
        near(report.residence_time_s, 1.0, 2e-14);
        near(report.elapsed_time_s, 1.0, 2e-14);
        near(report.travelled_distance_m, 0.05, 2e-14);
        near(report.state.volume_m3, before.volume_m3, 2e-12);
        near(
            report.state.axial_momentum_kg_m_s,
            before.axial_momentum_kg_m_s,
            2e-12,
        );
        assert!(report.state.minimum_radius_m < before.minimum_radius_m - 5e-5);
        assert!(
            flight
                .jet()
                .velocities_m_s()
                .any(|u| (u - 0.05).abs() > 1e-5)
        );
        assert_eq!(flight.initial_radii_m(), initial_radii);
        assert_eq!(flight.initial_velocities_m_s(), initial_velocities);
        let native = report.native.unwrap();
        assert!(native.viscous_dissipation_j >= 0.0);
        assert!(native.upwind_dissipation_j >= 0.0);
        assert!(native.energy_balance_residual_j.is_finite());

        let completed = flight.advance(0.1, budget(), &gate).unwrap();
        assert!(completed.native.is_none());
        assert_eq!(completed.state, report.state);
    }

    #[test]
    fn g3_authored_velocity_controls_residence_and_nonlinear_necking() {
        let gate = CancelGate::new_clock_free();
        let mut slow = NonlinearFreeJet::new(seeded_jet(0.05, 0.5), 0.05).unwrap();
        let mut fast = NonlinearFreeJet::new(seeded_jet(0.1, 0.5), 0.05).unwrap();
        let slow_report = slow.advance(2.0, budget(), &gate).unwrap();
        let fast_report = fast.advance(2.0, budget(), &gate).unwrap();
        near(
            fast_report.residence_time_s,
            0.5 * slow_report.residence_time_s,
            2e-14,
        );
        near(
            fast_report.travelled_distance_m,
            slow_report.travelled_distance_m,
            2e-14,
        );
        assert!(slow_report.state.minimum_radius_m < fast_report.state.minimum_radius_m);
        assert_eq!(slow_report.outcome, NonlinearFreeJetOutcome::FlightComplete);
        assert_eq!(fast_report.outcome, NonlinearFreeJetOutcome::FlightComplete);
    }

    #[test]
    fn g1_neck_limit_retains_partial_flight_and_the_last_admitted_shape() {
        let gate = CancelGate::new_clock_free();
        let mut flight = NonlinearFreeJet::new(seeded_jet(0.05, 0.9899), 0.1).unwrap();
        let report = flight.advance(5.0, budget(), &gate).unwrap();
        let NonlinearFreeJetOutcome::DomainLimit(violation) = report.outcome else {
            panic!("expected a physical neck limit, got {:?}", report.outcome);
        };
        assert_eq!(violation.kind, SlenderJetLimitKind::MinimumRadius);
        assert!(report.elapsed_time_s > 0.0);
        assert!(report.elapsed_time_s < report.residence_time_s);
        assert!(report.travelled_distance_m < report.requested_flight_length_m);
        assert!(report.state.minimum_radius_m >= 0.9899);
        assert!(flight.remaining_time_s() > 0.0);
        let radii = flight.jet().radius_squared_m2().to_vec();
        let stopped = flight.advance(5.0, budget(), &gate).unwrap();
        assert_eq!(stopped.outcome, report.outcome);
        assert!(stopped.native.is_none());
        assert_eq!(stopped.state, report.state);
        assert_eq!(flight.jet().radius_squared_m2(), radii);
    }

    #[test]
    fn g0_zero_flight_invalid_transport_and_budget_refusal_do_not_advance() {
        let gate = CancelGate::new_clock_free();
        let jet = seeded_jet(0.05, 0.5);
        let initial = jet.radius_squared_m2().to_vec();
        let mut zero = NonlinearFreeJet::new(jet, 0.0).unwrap();
        let report = zero.advance(1.0, budget(), &gate).unwrap();
        assert_eq!(report.outcome, NonlinearFreeJetOutcome::FlightComplete);
        assert_eq!(report.elapsed_time_s, 0.0);
        assert!(report.native.is_none());
        assert_eq!(zero.jet().radius_squared_m2(), initial);
        assert!(NonlinearFreeJet::new(seeded_jet(0.0, 0.5), 0.1).is_err());
        assert!(NonlinearFreeJet::new(seeded_jet(-0.05, 0.5), 0.1).is_err());
        assert!(NonlinearFreeJet::new(seeded_jet(0.05, 0.5), f64::NAN).is_err());

        let mut flight = NonlinearFreeJet::new(seeded_jet(0.05, 0.5), 0.05).unwrap();
        let radii = flight.jet().radius_squared_m2().to_vec();
        let flux = flight.jet().radius_squared_flux_m3_s().to_vec();
        let tiny_budget = SlenderJetStepBudget {
            maximum_substeps: 1,
            ..budget()
        };
        assert!(matches!(
            flight.advance(1.0, tiny_budget, &gate),
            Err(NonlinearJetError::SubstepBudget)
        ));
        assert_eq!(flight.report().elapsed_time_s, 0.0);
        assert_eq!(flight.jet().radius_squared_m2(), radii);
        assert_eq!(flight.jet().radius_squared_flux_m3_s(), flux);
        let resumed = flight.advance(2.0, budget(), &gate).unwrap();
        assert_eq!(resumed.outcome, NonlinearFreeJetOutcome::FlightComplete);
    }
    #[test]
    fn g1_implicit_viscosity_completes_high_viscosity_flight_with_bounded_work() {
        use fs_lbm::jet::nonlinear::SlenderJetIntegrationMethod;
        let gate = CancelGate::new_clock_free();
        let initial = seeded_jet(0.05, 0.5);
        let radii: Vec<_> = initial.radii_m().collect();
        let velocities: Vec<_> = initial.velocities_m_s().collect();
        let mut parameters = initial.parameters();
        parameters.kinematic_viscosity_m2_s = 50.0;
        let mut implicit = PeriodicSlenderJet::new(
            parameters,
            initial.length_m(),
            &radii,
            &velocities,
            initial.limits(),
            &gate,
        )
        .unwrap();
        let before = implicit.diagnostics();
        let work = SlenderJetStepBudget {
            maximum_substeps: 1000,
            ..budget()
        };
        // The original explicit viscous restriction cannot reach this same
        // physical residence within the authored work budget.
        let mut explicit = implicit.clone();
        assert!(matches!(
            explicit.advance(0.1, work, &gate),
            Err(NonlinearJetError::SubstepBudget),
        ));
        assert_eq!(explicit.time_s(), 0.0);
        implicit
            .set_integration_method(SlenderJetIntegrationMethod::ImplicitViscosity, &gate)
            .unwrap();
        let mut flight = NonlinearFreeJet::new(implicit, 0.005).unwrap();
        let first = flight.advance(0.04, work, &gate).unwrap();
        assert_eq!(first.outcome, NonlinearFreeJetOutcome::InFlight);
        let last = flight.advance(0.1, work, &gate).unwrap();
        assert_eq!(last.outcome, NonlinearFreeJetOutcome::FlightComplete);
        near(last.elapsed_time_s, 0.1, 2e-14);
        near(last.travelled_distance_m, 0.005, 2e-14);
        near(last.state.volume_m3, before.volume_m3, 3e-12);
        near(
            last.state.axial_momentum_kg_m_s,
            before.axial_momentum_kg_m_s,
            3e-12,
        );
        let native = last.native.unwrap();
        assert_eq!(
            native.integration_method,
            SlenderJetIntegrationMethod::ImplicitViscosity
        );
        assert!(native.accepted_substeps < work.maximum_substeps);
        let error = native.temporal_error.unwrap();
        assert!(error.maximum_radius_error_m.is_finite());
        assert!(error.maximum_velocity_error_m_s.is_finite());
        assert!(native.viscous_dissipation_j >= 0.0);
        assert_eq!(flight.initial_radii_m(), radii);
        assert_eq!(flight.initial_velocities_m_s(), velocities);
    }
}
