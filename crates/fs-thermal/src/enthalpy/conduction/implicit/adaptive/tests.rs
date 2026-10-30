use super::super::tests::{near, setup};
use super::*;

struct Law {
    datum: f64,
    latent: bool,
    phase: bool,
}
impl ThermalLaw for Law {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let x = h - self.datum;
        let maximum = if self.latent { 550.0 } else { 1000.0 };
        if !x.is_finite() || !(100.0..=maximum).contains(&x) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if !self.latent || x < 300.0 {
            x
        } else if x <= 500.0 {
            300.0
        } else {
            x - 200.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let primitive = |h| {
            let x = h - self.datum;
            let t = self.temperature_k(h)?;
            Ok::<_, DepositionError>(if !self.latent || x < 300.0 {
                t.ln()
            } else if x <= 500.0 {
                300.0_f64.ln() + (x - 300.0) / 300.0
            } else {
                t.ln() + 200.0 / 300.0
            })
        };
        Ok(Some(primitive(b)? - primitive(a)?))
    }
}
impl ConductivityLaw for Law {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.temperature_k(h)?;
        Ok(if self.latent {
            (600.0 - (h - self.datum)) / 200.0
        } else {
            1.0
        })
    }
}
impl ImplicitConductionLaw for Law {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [
            self.datum + 100.0,
            self.datum + if self.latent { 550.0 } else { 1000.0 },
        ]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: self.temperature_k(h)?,
            temperature_slope_kg_k_j: if self.latent && (300.0..=500.0).contains(&(h - self.datum))
            {
                0.0
            } else {
                1.0
            },
            conductivity_w_m_k: self.conductivity_w_m_k(h)?,
            conductivity_slope_kg_w_m_k_j: if self.latent { -0.005 } else { 0.0 },
        })
    }
    fn time_observables(&self, h: f64) -> Result<ImplicitConductionObservables, DepositionError> {
        Ok(ImplicitConductionObservables {
            temperature_k: self.temperature_k(h)?,
            phase_mass_fraction: self.phase.then(|| {
                if self.latent {
                    ((h - self.datum - 300.0) / 200.0).clamp(0.0, 1.0)
                } else {
                    0.0
                }
            }),
        })
    }
}
fn policy(dt: f64, absolute: f64) -> AdaptiveImplicitConductionControls {
    AdaptiveImplicitConductionControls {
        step: ImplicitConductionControls {
            maximum_dt_s: dt,
            cell_tolerance_j: 1e-10,
            total_tolerance_j: 1e-9,
            ..ImplicitConductionControls::default()
        },
        initial_dt_s: dt,
        maximum_trials: 4096,
        enthalpy_absolute_tolerance_j_kg: absolute,
        temperature_absolute_tolerance_k: absolute,
        phase_fraction_absolute_tolerance: None,
        relative_tolerance: 0.0,
    }
}

#[test]
fn adaptive_cooling_reduces_actual_error_without_enthalpy_datum_bias() {
    let law = Law {
        datum: 0.0,
        latent: false,
        phase: false,
    };
    let (thermal, observer) = setup(&law, &[400.0], Some(300.0), vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let mut loose = policy(1.0, 0.5);
    loose.relative_tolerance = 0.01;
    let mut tight = policy(1.0, 0.005);
    tight.relative_tolerance = 0.0001;
    let (a, ra) = chain
        .advance_adaptive_to(&thermal, 1.0, loose, &mut workspace, || Ok(()))
        .unwrap();
    let (b, rb) = chain
        .advance_adaptive_to(&thermal, 1.0, tight, &mut workspace, || Ok(()))
        .unwrap();
    let exact = 300.0 + 100.0 * (-1.0_f64).exp();
    let error_a = (a.temperature_k(CellId(0)).unwrap() - exact).abs();
    let error_b = (b.temperature_k(CellId(0)).unwrap() - exact).abs();
    assert!(error_b < 0.25 * error_a, "loose={error_a}, tight={error_b}");
    assert!(rb.accepted_comparisons > ra.accepted_comparisons);
    assert!(rb.maximum_accepted_normalized_error <= 1.0 && rb.accuracy_rejections > 0);
    assert_eq!(rb.solve.energy.substeps, 2 * rb.accepted_comparisons);
    near(
        rb.solve.energy.external_heat_j,
        b.cells()[&CellId(0)].specific_h_j_kg - 400.0,
        1e-9,
    );

    let shifted = Law {
        datum: 8192.0,
        latent: false,
        phase: false,
    };
    let (old, observer) = setup(&shifted, &[8592.0], Some(300.0), vec![]);
    let chain =
        ImplicitConductionChain::new(&old, &observer, &[(MaterialId(0), &shifted)], || Ok(()))
            .unwrap();
    let (same, shifted_report) = chain
        .advance_adaptive_to(&old, 1.0, tight, &mut workspace, || Ok(()))
        .unwrap();
    near(
        same.temperature_k(CellId(0)).unwrap(),
        b.temperature_k(CellId(0)).unwrap(),
        2e-8,
    );
    assert_eq!(shifted_report.accepted_comparisons, rb.accepted_comparisons);
    assert_eq!(shifted_report.accuracy_rejections, rb.accuracy_rejections);
}

#[test]
fn phase_observable_refines_latent_inventory_with_identical_temperatures() {
    let law = Law {
        datum: 0.0,
        latent: true,
        phase: true,
    };
    let (thermal, observer) = setup(&law, &[400.0], Some(250.0), vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let mut loose = policy(1.0, 1000.0);
    loose.phase_fraction_absolute_tolerance = Some(0.01);
    let mut tight = loose;
    tight.phase_fraction_absolute_tolerance = Some(0.0001);
    let (a, ra) = chain
        .advance_adaptive_to(&thermal, 1.0, loose, &mut workspace, || Ok(()))
        .unwrap();
    let (b, rb) = chain
        .advance_adaptive_to(&thermal, 1.0, tight, &mut workspace, || Ok(()))
        .unwrap();
    let exact = 600.0 - 200.0 * 0.25_f64.exp();
    let error_a = (a.cells()[&CellId(0)].specific_h_j_kg - exact).abs();
    let error_b = (b.cells()[&CellId(0)].specific_h_j_kg - exact).abs();
    assert!(error_b < 0.25 * error_a, "loose={error_a}, tight={error_b}");
    assert!(rb.accepted_comparisons > ra.accepted_comparisons && rb.accuracy_rejections > 0);
    assert_eq!(a.temperature_k(CellId(0)).unwrap(), 300.0);
    assert_eq!(b.temperature_k(CellId(0)).unwrap(), 300.0);
    assert_eq!(rb.last_error.unwrap().maximum_temperature_difference_k, 0.0);
    assert!(rb.maximum_accepted_normalized_error <= 1.0);
    let heat = b.cells()[&CellId(0)].specific_h_j_kg - 400.0;
    near(rb.solve.energy.external_heat_j, heat, 1e-9);
    near(
        rb.solve.entropy.thermal_storage_change_j_k,
        heat / 300.0,
        1e-12,
    );
    near(
        rb.solve.entropy.reservoir_entropy_input_j_k,
        heat / 250.0,
        1e-12,
    );
    assert!(rb.solve.entropy.conduction_generation_j_k > 0.0);
}

#[test]
fn exact_events_and_unequal_representable_halves_keep_fine_only_accounts() {
    let law = Law {
        datum: 0.0,
        latent: false,
        phase: false,
    };
    let event = TimedEvent {
        time_s: 0.35,
        event: Event::Boundary(
            CellId(0),
            Boundary {
                ambient_k: 300.0,
                htc_w_m2_k: 0.0,
                source_w: 0.0,
                reservoir: Some(Reservoir {
                    area_m2: 1.0,
                    conductance_w_k: 0.0,
                    temperature_k: 500.0,
                }),
            },
        ),
    };
    let (thermal, observer) = setup(&law, &[400.0], Some(200.0), vec![event]);
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, report) = chain
        .advance_adaptive_to(&thermal, 0.8, policy(0.8, 1e6), &mut workspace, || Ok(()))
        .unwrap();
    let mut exact = 400.0;
    let mut reservoir_entropy = 0.0;
    for (dt, bath) in [
        (0.175, 200.0),
        (0.175, 200.0),
        (0.225, 500.0),
        (0.225, 500.0),
    ] {
        let old = exact;
        exact = (old + dt * bath) / (1.0 + dt);
        reservoir_entropy += (exact - old) / bath;
    }
    near(next.temperature_k(CellId(0)).unwrap(), exact, 1e-12);
    near(report.solve.energy.external_heat_j, exact - 400.0, 1e-10);
    near(
        report.solve.entropy.reservoir_entropy_input_j_k,
        reservoir_entropy,
        1e-12,
    );
    near(
        report.solve.entropy.thermal_storage_change_j_k,
        (exact / 400.0).ln(),
        1e-12,
    );
    assert_eq!(report.accepted_comparisons, 2);
    assert_eq!(report.solve.energy.substeps, 4);
    assert_eq!(report.solve.energy.events, 1);
    assert!(report.solve.residual_evaluations >= 18);

    // At this clock magnitude, [start,start+3] splits into spans2 and1.
    // The fine-error correction is therefore1.25, not the equal-halves value1.
    let (mut thermal, observer) = setup(&law, &[400.0], Some(300.0), vec![]);
    thermal.time_s = 2.0_f64.powi(52);
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let (next, report) = chain
        .advance_adaptive_to(
            &thermal,
            thermal.time_s + 3.0,
            policy(3.0, 11.0),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    let fine = ((400.0 + 2.0 * 300.0) / 3.0 + 300.0) / 2.0;
    let coarse = (400.0 + 3.0 * 300.0) / 4.0;
    near(next.temperature_k(CellId(0)).unwrap(), fine, 1e-12);
    near(
        report.last_error.unwrap().normalized,
        1.25 * (fine - coarse).abs() / 11.0,
        1e-12,
    );
    assert_eq!(report.accepted_comparisons, 1);
}

#[test]
fn missing_observables_work_floors_and_late_cancellation_preserve_input() {
    let law = Law {
        datum: 0.0,
        latent: false,
        phase: false,
    };
    let (thermal, observer) = setup(&law, &[400.0], Some(300.0), vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let mut missing = policy(1.0, 1.0);
    missing.phase_fraction_absolute_tolerance = Some(0.01);
    assert!(matches!(
        chain.advance_adaptive_to(&thermal, 1.0, missing, &mut workspace, || Ok(())),
        Err(AdaptiveImplicitConductionError::Invalid(_))
    ));
    let mut floor = policy(1.0, 1e-8);
    floor.step.minimum_dt_s = 1.0;
    assert!(matches!(
        chain.advance_adaptive_to(&thermal, 1.0, floor, &mut workspace, || Ok(())),
        Err(AdaptiveImplicitConductionError::TimeAccuracy { .. })
    ));
    let mut attempts = policy(1.0, 1e-8);
    attempts.maximum_trials = 1;
    assert!(matches!(
        chain.advance_adaptive_to(&thermal, 1.0, attempts, &mut workspace, || Ok(())),
        Err(AdaptiveImplicitConductionError::TrialBudget)
    ));
    let mut work = policy(1.0, 1000.0);
    work.step.maximum_residual_evaluations = 4;
    assert!(matches!(
        chain.advance_adaptive_to(&thermal, 1.0, work, &mut workspace, || Ok(())),
        Err(AdaptiveImplicitConductionError::Implicit(
            ImplicitConductionError::Thermal(DepositionError::Budget(_))
        ))
    ));
    let mut checks = 0;
    chain
        .advance_adaptive_to(&thermal, 1.0, policy(1.0, 1000.0), &mut workspace, || {
            checks += 1;
            Ok(())
        })
        .unwrap();
    let last = checks;
    checks = 0;
    assert!(matches!(
        chain.advance_adaptive_to(&thermal, 1.0, policy(1.0, 1000.0), &mut workspace, || {
            checks += 1;
            if checks == last {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(AdaptiveImplicitConductionError::Implicit(
            ImplicitConductionError::Thermal(DepositionError::Cancelled)
        ))
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 400.0);
    assert_eq!(
        thermal.cells()[&CellId(0)]
            .boundary
            .reservoir
            .unwrap()
            .conductance_w_k,
        0.0
    );
}
