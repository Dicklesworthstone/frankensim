use super::*;

struct Entropic;
impl ThermalLaw for Entropic {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        SENSIBLE.temperature_k(h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        Ok(Some((self.temperature_k(b)? / self.temperature_k(a)?).ln()))
    }
}
static ENTROPIC: Entropic = Entropic;
fn body(mass: f64, source: f64) -> DepositionSimulation<'static> {
    let mut c = cell(1, mass, 10.0);
    c.boundary = Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 2.0,
        source_w: source,
        reservoir: Some(Reservoir {
            area_m2: 1.0,
            conductance_w_k: 3.0,
            temperature_k: 330.0,
        }),
    };
    simulation(&ENTROPIC, vec![c], None, vec![])
}
fn coefficients(
    state: &DepositionSimulation<'_>,
    extra: f64,
) -> Result<BoundaryStepCoefficients<BoundaryConductances>, DepositionError> {
    let t = state.temperature_k(CellId(1))?;
    let boundaries = BoundaryConductances {
        ambient_w_k: [(CellId(1), t / 320.0)].into(),
        reservoir_w_k: [(CellId(1), t / 80.0)].into(),
    };
    Ok(BoundaryStepCoefficients {
        bulk: StepCoefficients {
            sources_w: if extra == 0.0 {
                BTreeMap::new()
            } else {
                [(CellId(1), extra)].into()
            },
            conductances_w_k: BTreeMap::new(),
            payload: boundaries.clone(),
        },
        boundaries,
    })
}

#[test]
fn g1_current_boundary_conductances_share_energy_entropy_and_native_decay_bound() {
    let old = body(2.0, 0.0);
    let original = old.cells()[&CellId(1)].clone();
    let step = BackwardEulerStep::new(&old, 1.0, || Ok(())).unwrap();
    // Endpoint T=320: ambient Q=-20J, prescribed reservoir Q=+40J.
    // Both are actual physical boundary exchange, with no generated source.
    let trial = step
        .evaluate_with_boundary_conductances_and_entropy(
            &[20.0],
            |state| coefficients(state, 0.0),
            || Ok(()),
        )
        .unwrap();
    assert!(trial.maximum_cell_residual_j() < 1e-9);
    let entropy = *trial.entropy_audit().unwrap();
    near(
        entropy.thermal_storage_change_j_k,
        2.0 * (320.0_f64 / 310.0).ln(),
    );
    near(
        entropy.reservoir_entropy_input_j_k,
        -20.0 / 300.0 + 40.0 / 330.0,
    );
    near(
        entropy.conduction_generation_j_k,
        (-20.0 / 320.0) * (-20.0 / 300.0) + (40.0 / 320.0) * (10.0 / 330.0),
    );
    assert_eq!(entropy.unclassified_source_at_cell_j_k, 0.0);
    assert!(!entropy.has_unclassified_sources);
    let (next, used, audit) = trial.accept(1e-9, 1e-9, || Ok(())).unwrap();
    near(audit.external_heat_j, 20.0);
    near(audit.stored_energy_change_j, 20.0);
    assert_eq!(next.cells()[&CellId(1)].boundary, original.boundary);
    assert_eq!(next.cells()[&CellId(1)].geometry, original.geometry);
    assert_eq!(next.cells()[&CellId(1)].mass_kg, original.mass_kg);
    assert_eq!(old.cells()[&CellId(1)], original);
    near(
        next.external_step_power_with_boundary_conductances(
            &BTreeMap::new(),
            &BTreeMap::new(),
            &used,
            || Ok(()),
        )
        .unwrap(),
        20.0,
    );
    let rates = next
        .heat_rates_with_boundaries(100.0, &used, &mut || Ok(()))
        .unwrap();
    near(rates.explicit_dt_s, 0.8 * 2.0 / (1.0 + 4.0));
    near(rates.external_w, 20.0);
}

#[test]
fn g3_area_coefficients_do_not_rescale_applied_watts_or_change_ordinary_callbacks() {
    let old = body(3.0, 3.0);
    let step = BackwardEulerStep::new(&old, 1.0, || Ok(())).unwrap();
    let trial = step
        .evaluate_with_boundary_conductances_and_entropy(
            &[30.0],
            |state| coefficients(state, 7.0),
            || Ok(()),
        )
        .unwrap();
    near(
        trial
            .entropy_audit()
            .unwrap()
            .unclassified_source_at_cell_j_k,
        10.0 / 320.0,
    );
    let (next, _, audit) = trial.accept(1e-9, 1e-9, || Ok(())).unwrap();
    near(audit.external_heat_j, 30.0);
    assert_eq!(next.cells()[&CellId(1)].boundary.source_w, 3.0);
    let ordinary = step.evaluate(&[0.0], empty, || Ok(())).unwrap();
    let geometric = step
        .evaluate_with_boundary_conductances(
            &[0.0],
            |state| empty(state).map(BoundaryStepCoefficients::from),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(ordinary.residual_j(), geometric.residual_j());
    assert_eq!(ordinary.audit(), geometric.audit());
}

#[test]
fn g4_unknown_or_uncovered_boundaries_and_late_cancellation_return_no_candidate() {
    let old = body(2.0, 0.0);
    let original = old.cells()[&CellId(1)].clone();
    let step = BackwardEulerStep::new(&old, 1.0, || Ok(())).unwrap();
    for (id, g) in [(CellId(9), 1.0), (CellId(1), -1.0), (CellId(1), f64::NAN)] {
        let boundary = BoundaryConductances {
            ambient_w_k: [(id, g)].into(),
            ..Default::default()
        };
        assert!(
            old.external_step_power_with_boundary_conductances(
                &BTreeMap::new(),
                &BTreeMap::new(),
                &boundary,
                || Ok(())
            )
            .is_err()
        );
    }
    let no_reservoir = simulation(&ENTROPIC, vec![cell(1, 2.0, 10.0)], None, vec![]);
    let missing = BoundaryConductances {
        reservoir_w_k: [(CellId(1), 0.0)].into(),
        ..Default::default()
    };
    assert!(
        no_reservoir
            .external_step_power_with_boundary_conductances(
                &BTreeMap::new(),
                &BTreeMap::new(),
                &missing,
                || Ok(())
            )
            .is_err()
    );
    let mut covered = original.clone();
    covered.boundary.reservoir.as_mut().unwrap().area_m2 = 6.0;
    let covered = simulation(&ENTROPIC, vec![covered], None, vec![]);
    let nonzero = BoundaryConductances {
        ambient_w_k: [(CellId(1), 1.0)].into(),
        ..Default::default()
    };
    assert!(
        covered
            .external_step_power_with_boundary_conductances(
                &BTreeMap::new(),
                &BTreeMap::new(),
                &nonzero,
                || Ok(())
            )
            .is_err()
    );
    let mut polls = 0;
    step.evaluate_with_boundary_conductances(
        &[20.0],
        |s| coefficients(s, 0.0),
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    let mut calls = 0;
    assert!(matches!(
        step.evaluate_with_boundary_conductances(
            &[20.0],
            |s| coefficients(s, 0.0),
            || {
                calls += 1;
                if calls == polls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }
        ),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&CellId(1)], original);
    let (_, _, audit) = step
        .evaluate_with_boundary_conductances(&[20.0], |s| coefficients(s, 0.0), || Ok(()))
        .unwrap()
        .accept(1e-9, 1e-9, || Ok(()))
        .unwrap();
    near(audit.external_heat_j, 20.0);
}
