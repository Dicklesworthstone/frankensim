use super::*;

#[path = "tests/boundaries.rs"]
mod boundaries;

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(-100.0..=1.0e6).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(300.0 + h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
static SENSIBLE: Sensible = Sensible;
struct Latent;
impl ThermalLaw for Latent {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=1000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if h <= 20.0 {
            300.0 + h / 2.0
        } else if h <= 120.0 {
            310.0
        } else {
            310.0 + (h - 120.0) / 4.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.5
    }
}
static LATENT: Latent = Latent;
fn boundary() -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: 0.0,
    }
}
fn cell(id: u64, mass: f64, h: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(1),
        mass_kg: mass,
        specific_h_j_kg: h,
        boundary: boundary(),
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64 + 1.0, 0.0, 0.0],
            volume_m3: 1.0,
            surface_area_m2: 6.0,
        },
    }
}
fn simulation(
    law: &'static dyn ThermalLaw,
    cells: Vec<DepositedCell>,
    link: Option<f64>,
    later: Vec<TimedEvent>,
) -> DepositionSimulation<'static> {
    let mut events: Vec<_> = cells
        .into_iter()
        .map(|c| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(c),
        })
        .collect();
    if let Some(g) = link {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 9,
                a: CellId(1),
                b: CellId(2),
                area_m2: 1.0,
                conductance_w_k: g,
            }),
        });
    }
    events.extend(later);
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 20,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn empty(_: &DepositionSimulation<'_>) -> Result<StepCoefficients<()>, DepositionError> {
    Ok(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    })
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a} != {b}");
}

#[test]
fn g1_heterogeneous_two_body_backward_euler_exceeds_explicit_cfl_and_conserves() {
    let old = simulation(
        &SENSIBLE,
        vec![cell(1, 1.0, 10.0), cell(2, 2.0, 0.0)],
        Some(4.0),
        vec![],
    );
    let step = BackwardEulerStep::new(&old, 2.0, || Ok(())).unwrap();
    assert_eq!(step.cell_ids(), &[CellId(1), CellId(2)]);
    let exchange = 80.0 / 13.0;
    let trial = step
        .evaluate(&[-exchange, exchange], empty, || Ok(()))
        .unwrap();
    assert!(trial.maximum_cell_residual_j() < 1e-9);
    let (next, (), audit) = trial.accept(1e-9, 1e-9, || Ok(())).unwrap();
    near(next.temperature_k(CellId(1)).unwrap(), 310.0 - exchange);
    near(
        next.temperature_k(CellId(2)).unwrap(),
        300.0 + exchange / 2.0,
    );
    near(next.energy_j(), old.energy_j());
    assert_eq!(next.mass_kg(), 3.0);
    assert_eq!(audit.substeps, 1);
    assert_eq!(next.time_s(), 2.0);
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn g1_latent_energy_and_source_payload_use_the_same_endpoint_not_kelvin_updates() {
    let old = simulation(&LATENT, vec![cell(1, 1.0, 40.0)], None, vec![]);
    let step = BackwardEulerStep::new(&old, 0.5, || Ok(())).unwrap();
    let trial = step
        .evaluate(
            &[40.0],
            |state| {
                assert_eq!(state.time_s(), 0.5);
                assert_eq!(state.cells()[&CellId(1)].specific_h_j_kg, 80.0);
                assert_eq!(state.temperature_k(CellId(1))?, 310.0);
                Ok::<_, DepositionError>(StepCoefficients {
                    sources_w: [(CellId(1), 80.0)].into(),
                    conductances_w_k: BTreeMap::new(),
                    payload: state.cells()[&CellId(1)].clone(),
                })
            },
            || Ok(()),
        )
        .unwrap();
    let (next, payload, audit) = trial.accept(0.0, 0.0, || Ok(())).unwrap();
    assert_eq!(payload, next.cells()[&CellId(1)]);
    assert_eq!(next.cells()[&CellId(1)].boundary.source_w, 0.0);
    assert_eq!(audit.external_heat_j, 40.0);
    assert_eq!(
        next.temperature_k(CellId(1)).unwrap(),
        old.temperature_k(CellId(1)).unwrap()
    );
    assert_eq!(old.cells()[&CellId(1)].specific_h_j_kg, 40.0);
}

#[test]
fn g0_local_residuals_cannot_hide_inside_a_zero_global_balance() {
    let old = simulation(
        &SENSIBLE,
        vec![cell(1, 1.0, 0.0), cell(2, 1.0, 0.0)],
        None,
        vec![],
    );
    let step = BackwardEulerStep::new(&old, 1.0, || Ok(())).unwrap();
    let trial = step.evaluate(&[1.0, -1.0], empty, || Ok(())).unwrap();
    assert_eq!(trial.audit().residual_j, 0.0);
    assert_eq!(trial.maximum_cell_residual_j(), 1.0);
    assert!(matches!(
        trial.accept(0.5, 0.0, || Ok(())),
        Err(DepositionError::EnergyBalance(_))
    ));
    let trial = step.evaluate(&[0.25, 0.25], empty, || Ok(())).unwrap();
    assert!(matches!(
        trial.accept(0.3, 0.3, || Ok(())),
        Err(DepositionError::EnergyBalance(_))
    ));
}

#[test]
fn g3_nonlinear_endpoint_transport_reuses_flux_and_keeps_authored_contact() {
    let old = simulation(
        &SENSIBLE,
        vec![cell(1, 1.0, 0.0), cell(2, 1.0, 10.0)],
        Some(0.5),
        vec![],
    );
    let step = BackwardEulerStep::new(&old, 1.0 / 16.0, || Ok(())).unwrap();
    let trial = step
        .evaluate(
            &[1.0, -1.0],
            |state| {
                let g = 1.0 + state.temperature_k(CellId(1))? - 300.0;
                Ok::<_, DepositionError>(StepCoefficients {
                    sources_w: BTreeMap::new(),
                    conductances_w_k: [(9, g)].into(),
                    payload: g,
                })
            },
            || Ok(()),
        )
        .unwrap();
    assert_eq!(*trial.payload(), 2.0);
    let (next, g, audit) = trial.accept(0.0, 0.0, || Ok(())).unwrap();
    assert_eq!(g, 2.0);
    assert_eq!(next.links()[&9].conductance_w_k, 0.5);
    assert_eq!(old.contact_patch(9), next.contact_patch(9));
    assert_eq!(audit.external_heat_j, 0.0);
}

#[test]
fn g0_events_cannot_be_skipped_or_silently_applied_to_other_physics() {
    let later = TimedEvent {
        time_s: 0.5,
        event: Event::Boundary(
            CellId(1),
            Boundary {
                source_w: 2.0,
                ..boundary()
            },
        ),
    };
    let old = simulation(&SENSIBLE, vec![cell(1, 1.0, 0.0)], None, vec![later]);
    assert!(BackwardEulerStep::new(&old, 0.6, || Ok(())).is_err());
    let step = BackwardEulerStep::new(&old, 0.5, || Ok(())).unwrap();
    let (left, (), audit) = step
        .evaluate(&[0.0], empty, || Ok(()))
        .unwrap()
        .accept(0.0, 0.0, || Ok(()))
        .unwrap();
    assert_eq!(audit.events, 0);
    assert_eq!(left.cells()[&CellId(1)].boundary.source_w, 0.0);
    assert!(BackwardEulerStep::new(&left, 0.75, || Ok(())).is_err());
    let (right, event_audit) = left
        .advance_to(
            0.5,
            StepConfig {
                max_dt_s: 0.1,
                energy_tolerance_j: 0.0,
            },
            || Ok(()),
        )
        .unwrap();
    assert_eq!(event_audit.events, 1);
    assert_eq!(right.cells()[&CellId(1)].boundary.source_w, 2.0);
    assert!(BackwardEulerStep::new(&right, 0.75, || Ok(())).is_ok());
}

#[test]
fn g0_resolvable_increment_is_zero_exactly_when_the_stored_coordinate_cannot_move() {
    // 1.0 has spacing 2^-52 above and 2^-53 below: half a spacing ties back to 1.0 and the next float past it moves h.
    let up_tie = f64::EPSILON / 2.0;
    let down_tie = f64::EPSILON / 4.0;
    let zero = 0.0_f64.to_bits();
    for dead in [up_tie, up_tie / 2.0, 1e-30, f64::MIN_POSITIVE, -down_tie, -down_tie / 2.0, -1e-30, -0.0] {
        assert_eq!(resolvable_increment_j(1.0, dead, 1.0).to_bits(), zero, "{dead:e}");
    }
    for live in [up_tie.next_up(), (-down_tie).next_down(), 1.0, -0.5, f64::EPSILON] {
        assert_ne!(1.0 + live, 1.0, "{live:e}");
        assert_eq!(resolvable_increment_j(1.0, live, 1.0).to_bits(), live.to_bits(), "{live:e}");
    }
    // The boundary scales with the cell mass: 2 kg halves the specific increment, so the joule threshold doubles.
    assert_eq!(resolvable_increment_j(1.0, 2.0 * up_tie, 2.0).to_bits(), zero);
    let live_2kg = 2.0 * up_tie.next_up();
    assert_eq!(resolvable_increment_j(1.0, live_2kg, 2.0).to_bits(), live_2kg.to_bits());
    // A nonfinite value is not a dead increment: it passes through for the owner's refusal.
    assert!(resolvable_increment_j(1.0, f64::NAN, 1.0).is_nan());
    assert_eq!(resolvable_increment_j(1.0, f64::INFINITY, 1.0), f64::INFINITY);
    // Dry-contact fixture magnitudes (h about 5.9e3 J/kg, 3.1e-6 kg cells): a 1e-9 J difference step moves h, 1e-20 J cannot.
    assert_eq!(resolvable_increment_j(5863.0, 1e-20, 3.12e-6).to_bits(), zero);
    assert_eq!(resolvable_increment_j(5863.0, 1e-9, 3.12e-6).to_bits(), 1e-9_f64.to_bits());
}

#[test]
fn g4_a_projected_iterate_evaluates_while_the_raw_subresolution_increment_is_still_refused() {
    let old = simulation(&LATENT, vec![cell(1, 1.0, 40.0)], None, vec![]);
    let step = BackwardEulerStep::new(&old, 0.5, || Ok(())).unwrap();
    // The guard is unchanged: a nonzero increment that cannot move h is refused when passed to evaluate directly.
    assert!(step.evaluate(&[1e-30], empty, || Ok(())).is_err());
    // The projected iterate is the exactly representable zero and evaluates; a resolvable one is left bit for bit.
    let projected = step.resolvable_increments_j(&[1e-30]);
    assert_eq!(projected[0].to_bits(), 0.0_f64.to_bits());
    assert!(step.evaluate(&projected, empty, || Ok(())).is_ok());
    let live = step.resolvable_increments_j(&[0.25]);
    assert_eq!(live[0].to_bits(), 0.25_f64.to_bits());
    assert!(step.evaluate(&live, empty, || Ok(())).is_ok());
    // Nonfinite and wrong-length vectors are not repaired: evaluate still refuses them.
    assert!(step.evaluate(&step.resolvable_increments_j(&[f64::NAN]), empty, || Ok(())).is_err());
    assert!(step.evaluate(&step.resolvable_increments_j(&[0.0, 1e-30]), empty, || Ok(())).is_err());
}

#[test]
fn g4_invalid_trials_and_every_cancelled_evaluation_leave_the_old_history_untouched() {
    let old = simulation(&LATENT, vec![cell(1, 1.0, 40.0)], None, vec![]);
    let step = BackwardEulerStep::new(&old, 0.5, || Ok(())).unwrap();
    for vector in [vec![], vec![f64::NAN], vec![-50.0], vec![1e-30]] {
        assert!(step.evaluate(&vector, empty, || Ok(())).is_err());
    }
    assert!(
        step.evaluate(
            &[0.0],
            |_| Ok::<_, DepositionError>(StepCoefficients {
                sources_w: [(CellId(2), 0.0)].into(),
                conductances_w_k: BTreeMap::new(),
                payload: (),
            }),
            || Ok(())
        )
        .is_err()
    );
    let mut count = 0;
    step.evaluate(&[0.0], empty, || {
        count += 1;
        Ok(())
    })
    .unwrap();
    for stop in 1..=count {
        let mut seen = 0;
        let refused = step.evaluate(&[0.0], empty, || {
            seen += 1;
            if seen == stop {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(refused, Err(DepositionError::Cancelled)));
    }
    let trial = step.evaluate(&[0.0], empty, || Ok(())).unwrap();
    assert!(matches!(
        trial.accept(0.0, 0.0, || Err(DepositionError::Cancelled)),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&CellId(1)].specific_h_j_kg, 40.0);
}

#[test]
fn external_rate_query_reuses_native_boundaries_without_internal_transfers() {
    let mut hot = cell(1, 1.0, 10.0);
    hot.boundary = Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 2.0,
        reservoir: Some(Reservoir {
            area_m2: 1.0,
            conductance_w_k: 3.0,
            temperature_k: 320.0,
        }),
        source_w: 4.0,
    };
    let old = simulation(&SENSIBLE, vec![hot, cell(2, 1.0, 0.0)], Some(0.5), vec![]);
    let sources = [(CellId(1), 5.0), (CellId(2), -2.0)].into();
    let links = [(9, 2.0)].into();
    // Cell 1 has 6 - 1 contact - 1 reservoir = 4 m2 exposed.
    // Ambient -80 W, reservoir +30 W, authored +4 W, additions +3 W.
    near(
        old.external_step_power_w(&sources, &links, || Ok(()))
            .unwrap(),
        -43.0,
    );
    let step = BackwardEulerStep::new(&old, 0.5, || Ok(())).unwrap();
    let trial = step
        .evaluate(
            &[0.0, 0.0],
            |_| {
                Ok::<_, DepositionError>(StepCoefficients {
                    sources_w: sources.clone(),
                    conductances_w_k: links.clone(),
                    payload: (),
                })
            },
            || Ok(()),
        )
        .unwrap();
    near(trial.audit().external_heat_j, -21.5);
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&CellId(1)].boundary.source_w, 4.0);
    assert_eq!(old.links()[&9].conductance_w_k, 0.5);
}

#[test]
fn external_rate_query_refuses_invalid_maps_and_late_cancellation_without_mutation() {
    let old = simulation(&SENSIBLE, vec![cell(1, 1.0, 10.0)], None, vec![]);
    let sources = [(CellId(1), 5.0)].into();
    let links = BTreeMap::new();
    assert!(
        old.external_step_power_w(&[(CellId(9), 1.0)].into(), &links, || Ok(()))
            .is_err()
    );
    let mut polls = 0;
    old.external_step_power_w(&sources, &links, || {
        polls += 1;
        Ok(())
    })
    .unwrap();
    let mut count = 0;
    assert!(matches!(
        old.external_step_power_w(&sources, &links, || {
            count += 1;
            if count == polls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&CellId(1)].specific_h_j_kg, 10.0);
    assert_eq!(old.cells()[&CellId(1)].boundary.source_w, 0.0);
}

#[test]
fn resolved_endpoint_reservoir_uses_original_heat_and_actual_bath_entropy() {
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
            Ok(Some(fs_math::det::ln(
                self.temperature_k(b)? / self.temperature_k(a)?,
            )))
        }
    }
    static LAW: Entropic = Entropic;
    let mut sample = cell(1, 2.0, 10.0);
    sample.boundary.reservoir = Some(Reservoir {
        area_m2: 1.0,
        conductance_w_k: 0.25,
        temperature_k: 325.0,
    });
    let old = simulation(&LAW, vec![sample], None, vec![]);
    let step = BackwardEulerStep::new(&old, 0.5, || Ok(())).unwrap();
    let trial = step
        .evaluate_with_entropy_and_reservoirs(
            &[10.0],
            |state| {
                let g = state.cells()[&CellId(1)].specific_h_j_kg / 7.5;
                Ok::<_, DepositionError>((
                    StepCoefficients {
                        sources_w: BTreeMap::new(),
                        conductances_w_k: BTreeMap::new(),
                        payload: g,
                    },
                    [(CellId(1), g)].into(),
                ))
            },
            || Ok(()),
        )
        .unwrap();
    assert_eq!(trial.maximum_cell_residual_j(), 0.0);
    assert_eq!(*trial.payload(), 2.0);
    let entropy = *trial.entropy_audit().unwrap();
    near(
        entropy.thermal_storage_change_j_k,
        2.0 * (315.0_f64 / 310.0).ln(),
    );
    near(entropy.reservoir_entropy_input_j_k, 10.0 / 325.0);
    near(
        entropy.conduction_generation_j_k,
        10.0 * (325.0 - 315.0) / (315.0 * 325.0),
    );
    assert!(!entropy.has_unclassified_sources);
    assert_eq!(entropy.unclassified_source_at_cell_j_k, 0.0);
    let (next, _, audit) = trial.accept(0.0, 0.0, || Ok(())).unwrap();
    assert_eq!(audit.external_heat_j, 10.0);
    assert_eq!(next.temperature_k(CellId(1)).unwrap(), 315.0);
    assert_eq!(
        next.cells()[&CellId(1)].boundary,
        old.cells()[&CellId(1)].boundary
    );
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.temperature_k(CellId(1)).unwrap(), 310.0);
}

#[test]
fn endpoint_reservoir_admission_and_late_cancellation_preserve_original_history() {
    let mut sample = cell(1, 1.0, 10.0);
    sample.boundary.reservoir = Some(Reservoir {
        area_m2: 1.0,
        conductance_w_k: 3.0,
        temperature_k: 320.0,
    });
    let old = simulation(&SENSIBLE, vec![sample], None, vec![]);
    let step = BackwardEulerStep::new(&old, 1.0, || Ok(())).unwrap();
    for (id, g) in [(CellId(1), -1.0), (CellId(1), f64::NAN), (CellId(2), 0.0)] {
        assert!(
            step.evaluate_with_reservoirs(
                &[0.0],
                |state| { Ok::<_, DepositionError>((empty(state)?, [(id, g)].into())) },
                || Ok(())
            )
            .is_err()
        );
    }
    let insulating = |state: &DepositionSimulation<'_>| {
        Ok::<_, DepositionError>((empty(state)?, [(CellId(1), 0.0)].into()))
    };
    let mut polls = 0;
    let trial = step
        .evaluate_with_reservoirs(&[0.0], insulating, || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    assert_eq!(trial.audit().external_heat_j, 0.0);
    let mut seen = 0;
    assert!(matches!(
        step.evaluate_with_reservoirs(&[0.0], insulating, || {
            seen += 1;
            if seen == polls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    ));
    let (next, (), _) = trial.accept(0.0, 0.0, || Ok(())).unwrap();
    assert_eq!(
        next.cells()[&CellId(1)]
            .boundary
            .reservoir
            .unwrap()
            .conductance_w_k,
        3.0
    );
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&CellId(1)].specific_h_j_kg, 10.0);
}
