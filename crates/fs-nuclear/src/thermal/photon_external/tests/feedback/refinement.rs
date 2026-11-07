//! G1/G3/G4: complete native photon/thermal time refinement.

use super::*;

fn time_controls() -> PhotonTimeRefinement {
    PhotonTimeRefinement {
        initial_max_dt_s: 1.0,
        minimum_max_dt_s: 1e-6,
        maximum_halvings: 12,
        maximum_total_substeps: 4000,
        maximum_total_transport_preparations: 4000,
        maximum_total_transport_work: 200_000_000,
        specific_energy_absolute_tolerance_j_kg: 2e-3,
        temperature_absolute_tolerance_k: 1e-3,
        photon_energy_absolute_tolerance_j: 5e-3,
        relative_tolerance: 0.0,
        conservation_tolerance_j: 1e-8,
    }
}

fn loose_controls() -> PhotonTimeRefinement {
    PhotonTimeRefinement {
        specific_energy_absolute_tolerance_j_kg: 100.0,
        temperature_absolute_tolerance_k: 100.0,
        photon_energy_absolute_tolerance_j: 100.0,
        maximum_halvings: 2,
        ..time_controls()
    }
}

fn native(
    law: &'static dyn ThermalLaw,
    stored_h: f64,
    htc: f64,
    later: Vec<TimedEvent>,
) -> DepositionSimulation<'static> {
    let mut events = vec![];
    for (id, material, start, mass) in [(11, 1, 0.0, 2.0), (29, 2, 1.0, 3.0)] {
        let mut c = cell(id, material, start, mass, 0.0);
        c.specific_h_j_kg = stored_h;
        if id == 11 {
            c.boundary.htc_w_m2_k = htc;
        }
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(c),
        });
    }
    events.extend(later);
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law), (MaterialId(2), law)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, controls(1.0), || Ok(()))
    .unwrap()
    .0
}

#[test]
fn g1_refined_photon_feedback_converges_and_returns_original_cell_heat() {
    let original = thermal(&LAW, 2.0, 0.0);
    let owner = feedback(
        &original,
        temperature_layers(true, 320.0),
        vec![steady()],
        None,
        PhotonFeedbackBudget::default(),
    )
    .unwrap();
    let c = time_controls();
    let (next, accepted, report) = owner.advance_refined(&original, 4.0, c, || Ok(())).unwrap();
    // Independent analytic solution: exp(0.05*theta) obeys a linear ODE.
    let theta = ((-0.5_f64).exp() + (1.0 - (-0.5_f64).exp()) * 0.4_f64.exp()).ln() / 0.05;
    let exact_rear = 300.0 + (32.0 - 4.0 * theta) * (1.0 - (-1.0_f64).exp()) / 6.0;
    assert!((next.temperature_k(CellId(11)).unwrap() - 300.0 - theta).abs() < 1.1e-3);
    assert!((next.temperature_k(CellId(29)).unwrap() - exact_rear).abs() < 1.1e-3);
    assert!(report.halvings >= 2);
    assert!(report.normalized_difference <= 1.0);
    assert!(
        report
            .observed_order
            .is_some_and(|p| (0.8..1.2).contains(&p))
    );
    assert!(!report.roundoff_limited);
    assert!(report.total_substeps > accepted.audit().substeps);
    assert_eq!(
        report.total_transport_preparations,
        report.total_substeps - (report.halvings + 1),
    );
    assert!(report.total_transport_work > accepted.feedback_work().unwrap().transport_work);
    close(
        report.receiver_thermal_credit_j[&CellId(11)],
        4.0 * (next.temperature_k(CellId(11)).unwrap() - 300.0),
    );
    close(
        report.receiver_thermal_credit_j[&CellId(29)],
        6.0 * (next.temperature_k(CellId(29)).unwrap() - 300.0),
    );
    close(
        report.receiver_thermal_credit_j.values().sum(),
        accepted.audit().thermal_credit_j,
    );
    close(report.thermal.external_heat_j, accepted.audit().absorbed_j);
    close(accepted.audit().incident_j, 32.0);
    close(
        accepted.audit().absorbed_j + accepted.audit().escaped_j,
        32.0,
    );
    // Refinement publishes an actual native trajectory, with no extrapolation.
    let (direct, direct_source, _) = owner
        .advance_to(
            &original,
            4.0,
            controls(report.accepted_max_dt_s),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(direct.cells(), next.cells());
    assert_eq!(direct_source.audit(), accepted.audit());
    assert_eq!(owner.time_s(), 0.0);
    assert_eq!(original.temperature_k(CellId(11)).unwrap(), 300.0);
}

#[test]
fn g3_affine_source_knots_native_events_and_restart_keep_original_accounts() {
    let original = native(
        &LAW,
        600.0,
        0.0,
        vec![TimedEvent {
            time_s: 1.25,
            event: Event::Boundary(
                CellId(11),
                Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 2.0,
                },
            ),
        }],
    );
    let owner = bind(&original, vec![pulse()], false);
    let c = PhotonTimeRefinement {
        initial_max_dt_s: 100.0,
        maximum_total_transport_preparations: 0,
        maximum_total_transport_work: 0,
        ..time_controls()
    };
    let (whole, accepted, report) = owner.advance_refined(&original, 3.0, c, || Ok(())).unwrap();
    assert_eq!(report.halvings, 2);
    assert!(report.roundoff_limited);
    assert_eq!(report.observed_order, None);
    assert!(report.accepted_max_dt_s <= 0.4375);
    assert!(report.largest_accepted_interval_s <= 0.4375);
    assert_eq!(report.thermal.events, 1);
    assert_eq!(report.total_transport_preparations, 0);
    assert_eq!(report.total_transport_work, 0);
    close(accepted.audit().incident_j, 45.0);
    close(
        report.thermal.external_heat_j,
        accepted.audit().absorbed_j + 3.5,
    );
    let (middle, prefix, prefix_report) = owner
        .advance_refined(&original, 1.25, c, || Ok(()))
        .unwrap();
    let (split, suffix, suffix_report) =
        prefix.advance_refined(&middle, 3.0, c, || Ok(())).unwrap();
    for id in [CellId(11), CellId(29)] {
        close(
            whole.temperature_k(id).unwrap(),
            split.temperature_k(id).unwrap(),
        );
        close(
            prefix_report.receiver_thermal_credit_j[&id]
                + suffix_report.receiver_thermal_credit_j[&id],
            report.receiver_thermal_credit_j[&id],
        );
    }
    close(prefix.audit().incident_j + suffix.audit().incident_j, 45.0);
    close(
        prefix.audit().escaped_j + suffix.audit().escaped_j,
        accepted.audit().escaped_j,
    );
}

#[test]
fn g4_discarded_trajectories_spend_work_and_returned_budgets_are_original() {
    let original = thermal(&LAW, 2.0, 0.0);
    let owner = feedback(
        &original,
        temperature_layers(true, 320.0),
        vec![steady()],
        None,
        PhotonFeedbackBudget::default(),
    )
    .unwrap();
    let cost = owner.feedback_work().unwrap().transport_work;
    let c = PhotonTimeRefinement {
        maximum_total_substeps: 7,
        maximum_total_transport_preparations: 4,
        maximum_total_transport_work: 4 * cost,
        ..loose_controls()
    };
    let (next, accepted, report) = owner.advance_refined(&original, 1.0, c, || Ok(())).unwrap();
    assert_eq!(report.total_substeps, 7); // 1 + 2 + 4 actual native steps.
    assert_eq!(report.total_transport_preparations, 4); // 0 + 1 + 3 rebuilds.
    assert_eq!(report.total_transport_work, 4 * cost);
    assert_eq!(accepted.feedback_work().unwrap().preparations, 3);
    // The accepted trial had only three rebuilds left; that temporary ceiling
    // must not replace the original source budget on its next physical advance.
    let (_, continued, _) = accepted
        .advance_to(&next, 2.0, controls(0.25), || Ok(()))
        .unwrap();
    assert_eq!(continued.feedback_work().unwrap().preparations, 4);
    for limited in [
        PhotonTimeRefinement {
            maximum_total_substeps: 6,
            ..c
        },
        PhotonTimeRefinement {
            maximum_total_transport_preparations: 3,
            ..c
        },
        PhotonTimeRefinement {
            maximum_total_transport_work: 3 * cost,
            ..c
        },
    ] {
        assert!(matches!(
            owner.advance_refined(&original, 1.0, limited, || Ok(())),
            Err(FissionThermalError::Radiation(RadiationError::Budget(_)))
                | Err(FissionThermalError::Thermal(DepositionError::Budget(_)))
        ));
    }
    let per_run = feedback(
        &original,
        temperature_layers(true, 320.0),
        vec![steady()],
        None,
        PhotonFeedbackBudget {
            max_transport_preparations: 2,
            ..PhotonFeedbackBudget::default()
        },
    )
    .unwrap();
    assert!(matches!(
        per_run.advance_refined(&original, 1.0, loose_controls(), || Ok(())),
        Err(FissionThermalError::Radiation(RadiationError::Budget(_)))
    ));
    assert_eq!(owner.audit(), ExternalPhotonHeatingAudit::default());
    assert_eq!(owner.feedback_work().unwrap().preparations, 1);
}

#[test]
fn g3_native_cfl_cannot_hide_an_inactive_time_refinement_axis() {
    // G=32 W/K, C=4 J/K, native explicit bound 0.8*C/G = 0.1 s.
    let original = native(&LAW, 600.0, 16.0 / 3.0, vec![]);
    let owner = bind(&original, vec![steady()], false);
    let (next, accepted, report) = owner
        .advance_refined(
            &original,
            0.2,
            PhotonTimeRefinement {
                initial_max_dt_s: 100.0,
                specific_energy_absolute_tolerance_j_kg: 0.02,
                temperature_absolute_tolerance_k: 0.01,
                photon_energy_absolute_tolerance_j: 0.01,
                ..time_controls()
            },
            || Ok(()),
        )
        .unwrap();
    assert!(report.accepted_max_dt_s <= 0.025 * (1.0 + 1e-13));
    assert!(report.largest_accepted_interval_s <= 0.025 * (1.0 + 1e-13));
    assert!(accepted.audit().substeps >= 8);
    assert!(report.total_substeps >= 14);
    assert!(report.maximum_temperature_difference_k > 0.0);
    let exact_rise = 8.0 * (1.0 - (-0.5_f64).exp()) / 32.0 * (1.0 - (-1.6_f64).exp());
    assert!((next.temperature_k(CellId(11)).unwrap() - 300.0 - exact_rise).abs() < 0.004);
}

struct ShiftedCaloric;
impl ThermalLaw for ShiftedCaloric {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([119; 32])
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        LAW.temperature_k(h - 100_000.0)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.5
    }
}
static SHIFTED: ShiftedCaloric = ShiftedCaloric;

#[test]
fn g3_caloric_reference_and_latent_storage_do_not_hide_time_error() {
    let c = PhotonTimeRefinement {
        specific_energy_absolute_tolerance_j_kg: 1e-7,
        temperature_absolute_tolerance_k: 1e6,
        photon_energy_absolute_tolerance_j: 1e6,
        relative_tolerance: 1e-3,
        ..time_controls()
    };
    let mut reports = vec![];
    for (law, h) in [(&LAW as &dyn ThermalLaw, 600.0), (&SHIFTED, 100_600.0)] {
        let original = native(law, h, 0.0, vec![]);
        let owner = feedback(
            &original,
            temperature_layers(true, 320.0),
            vec![steady()],
            None,
            PhotonFeedbackBudget::default(),
        )
        .unwrap();
        let (_, _, report) = owner.advance_refined(&original, 4.0, c, || Ok(())).unwrap();
        assert!(report.halvings > 2);
        assert!(report.normalized_difference <= 1.0);
        reports.push(report);
    }
    assert_eq!(reports[0].halvings, reports[1].halvings);
    assert!((reports[0].normalized_difference - reports[1].normalized_difference).abs() < 1e-5);

    let original = thermal(&LATENT, 2.0, 0.0);
    let owner = feedback(
        &original,
        temperature_layers(true, 320.0),
        vec![steady()],
        None,
        PhotonFeedbackBudget::default(),
    )
    .unwrap();
    let (next, accepted, report) = owner
        .advance_refined(
            &original,
            2.0,
            PhotonTimeRefinement {
                maximum_total_transport_preparations: 0,
                maximum_total_transport_work: 0,
                ..time_controls()
            },
            || Ok(()),
        )
        .unwrap();
    assert!(report.roundoff_limited);
    assert_eq!(report.observed_order, None);
    assert_eq!(report.total_transport_preparations, 0);
    assert!(Arc::ptr_eq(owner.response(), accepted.response()));
    for id in [CellId(11), CellId(29)] {
        assert_eq!(next.temperature_k(id).unwrap(), 300.0);
        assert!(next.cells()[&id].specific_h_j_kg > 600.0);
        assert!(report.receiver_thermal_credit_j[&id] > 0.0);
    }
}

#[test]
fn g4_accuracy_floor_domain_and_cancellation_refuse_without_partial_publication() {
    let original = thermal(&LAW, 2.0, 0.0);
    let owner = feedback(
        &original,
        temperature_layers(true, 320.0),
        vec![steady()],
        None,
        PhotonFeedbackBudget::default(),
    )
    .unwrap();
    for c in [
        PhotonTimeRefinement {
            specific_energy_absolute_tolerance_j_kg: 1e-15,
            temperature_absolute_tolerance_k: 1e-15,
            photon_energy_absolute_tolerance_j: 1e-15,
            ..loose_controls()
        },
        PhotonTimeRefinement {
            minimum_max_dt_s: 0.3,
            ..loose_controls()
        },
    ] {
        assert!(matches!(
            owner.advance_refined(&original, 1.0, c, || Ok(())),
            Err(FissionThermalError::PhotonTimeConvergence { .. })
        ));
    }
    let narrow = feedback(
        &original,
        temperature_layers(true, 300.01),
        vec![steady()],
        None,
        PhotonFeedbackBudget::default(),
    )
    .unwrap();
    assert!(matches!(
        narrow.advance_refined(&original, 1.0, time_controls(), || Ok(())),
        Err(FissionThermalError::PhotonTemperatureDomain { .. })
    ));
    let (_, idle, idle_report) = owner
        .advance_refined(
            &original,
            0.0,
            PhotonTimeRefinement {
                maximum_total_substeps: 0,
                maximum_total_transport_preparations: 0,
                maximum_total_transport_work: 0,
                ..time_controls()
            },
            || Ok(()),
        )
        .unwrap();
    assert_eq!(idle_report.total_substeps, 0);
    assert_eq!(idle_report.halvings, 0);
    assert_eq!(idle.audit(), ExternalPhotonHeatingAudit::default());

    let mut calls = 0;
    owner
        .advance_refined(&original, 1.0, loose_controls(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    let expected = DepositionError::InvalidInput("original refinement callback");
    let response = Arc::clone(owner.response());
    let initial_work = owner.feedback_work();
    for stop in [1, calls / 2, calls] {
        let mut visited = 0;
        assert!(matches!(
            owner.advance_refined(&original, 1.0, loose_controls(), || {
                visited += 1;
                if visited == stop { Err(expected.clone()) } else { Ok(()) }
            }),
            Err(FissionThermalError::Thermal(error)) if error == expected
        ));
        assert_eq!(owner.time_s(), 0.0);
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(owner.audit(), ExternalPhotonHeatingAudit::default());
        assert_eq!(owner.feedback_work(), initial_work);
        assert!(Arc::ptr_eq(&response, owner.response()));
    }
}
