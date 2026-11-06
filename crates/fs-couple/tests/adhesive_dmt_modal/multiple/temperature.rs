//! Original per-contact temperatures retain one complete native mechanical owner.
use super::*;
#[path = "../../../../fs-contact/tests/adhesive_temperature/common.rs"]
mod thermal;
use fs_contact::normal_patch::adhesive::temperature::AdhesiveTemperatures;
use fs_contact::normal_patch::dmt::coupled::DmtCoupledError;
use fs_couple::render::schedule::force::coupled::contact::adhesive_dmt::multiple::temperature::{
    AdhesiveDmtSourceRefresh, MultipleAdhesiveDmtSourceRefreshReport,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_math::dd::Dd;

const COLD: AdhesiveTemperatures = AdhesiveTemperatures {
    surface_a_k: 300.0,
    surface_b_k: 300.0,
    interface_k: 300.0,
};
const TEMPERATURES: [AdhesiveTemperatures; 2] = [
    AdhesiveTemperatures {
        surface_a_k: 310.0,
        surface_b_k: 290.0,
        interface_k: 315.0,
    },
    AdhesiveTemperatures {
        surface_a_k: 295.0,
        surface_b_k: 305.0,
        interface_k: 288.0,
    },
];

fn shear(normal: &SourceDmtSphereContact) -> SourceAdhesiveShear {
    SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

// Dimensional potential from the authored bulk/interface curves, independent
// of both the refresh receipt and native adapter's source-work account.
fn potential(delta: f64, w300: f64, t: AdhesiveTemperatures) -> f64 {
    let ea = 1e11 * (1.0 - 0.01 * (t.surface_a_k - 300.0));
    let eb = 2e11 * (1.0 - 0.01 * (t.surface_b_k - 300.0));
    let e = 1.0 / ((1.0 - 0.25_f64.powi(2)) / ea + (1.0 - 0.25_f64.powi(2)) / eb);
    let w = w300 * (1.0 - 0.03 * (t.interface_k - 300.0));
    let pull = 2.0 * PI * 1e-6 * w;
    if delta >= 0.0 {
        (8.0 / 15.0) * e * 1e-3 * delta.powf(2.5) - pull * (delta + 0.5 * Z)
    } else if delta > -Z {
        -0.5 * pull * (delta + Z).powi(2) / Z
    } else {
        0.0
    }
}

fn mean_area(a: f64, b: f64) -> f64 {
    let approach = if a >= 0.0 && b >= 0.0 {
        0.5 * (a + b)
    } else if a <= 0.0 && b <= 0.0 {
        0.0
    } else {
        0.5 * a.max(b).powi(2) / (a - b).abs()
    };
    PI * 1e-6 * approach
}

#[test]
fn independent_temperatures_preserve_all_native_histories_and_drive_successor_shear() {
    let f0 = thermal::Fixture::with_options(1e5, 0.01, true, 320.0);
    let f1 = thermal::Fixture::with_options(1e5, 0.014, true, 320.0);
    let normal = [f0.dmt(Z), f1.dmt(Z)];
    let shears = normal.each_ref().map(shear);
    let make_owner = || {
        MultipleAdhesiveDmtModalSystem::new(
            rich_network(),
            vec![
                contact(&normal[0], &shears[0], [0.20 * Z, 0.0, 0.0], columns()),
                contact(
                    &normal[1],
                    &shears[1],
                    [0.17 * Z, 0.0, 0.0],
                    second_columns(),
                ),
            ],
            controls(),
            &gate(),
        )
        .unwrap()
    };
    let mut owner = make_owner();
    let mut replay = make_owner();
    let rates = [[0.0, 1e-3, -0.5e-3], [0.0, -0.7e-3, 0.3e-3]];
    owner.step(&[0.0; 6], &rates).unwrap();
    replay.step(&[0.0; 6], &rates).unwrap();
    let native = bits(owner.network());
    let frame = owner.last_frame().unwrap().clone();
    let native_frame = owner.network().last_frame().unwrap().clone();
    let before_states = [owner.state(0).unwrap(), owner.state(1).unwrap()];
    let compliance = owner.compliance_m_per_n().to_vec();
    let maps = [
        owner.work_columns(0).unwrap().clone(),
        owner.work_columns(1).unwrap().clone(),
    ];
    assert!(owner.network().material_memory().iter().any(|x| *x != 0.0));
    assert!(owner.network().radiation_state().iter().any(|x| *x != 0.0));
    assert!(
        before_states
            .iter()
            .all(|s| s.regime() == DmtCoupledRegime::Sliding)
    );
    let refreshed = [
        normal[0]
            .refresh_at_held_indentation(
                before_states[0].position_m()[0],
                f0.cards(),
                TEMPERATURES[0],
            )
            .unwrap(),
        normal[1]
            .refresh_at_held_indentation(
                before_states[1].position_m()[0],
                f1.cards(),
                TEMPERATURES[1],
            )
            .unwrap(),
    ];
    let next = [
        shears[0]
            .rebind_dmt_temperature(&normal[0], &refreshed[0])
            .unwrap(),
        shears[1]
            .rebind_dmt_temperature(&normal[1], &refreshed[1])
            .unwrap(),
    ];
    let requests = [
        AdhesiveDmtSourceRefresh {
            contact: 0,
            refreshed: &refreshed[0],
            next_shear: &next[0],
        },
        AdhesiveDmtSourceRefresh {
            contact: 1,
            refreshed: &refreshed[1],
            next_shear: &next[1],
        },
    ];
    let before_energy = owner.total_energy_j().unwrap();
    let admitted: Result<MultipleAdhesiveDmtSourceRefreshReport, MultipleAdhesiveDmtModalError> =
        owner.apply_source_refreshes_with_admission(
            &[requests[1], requests[0]],
            &gate(),
            |report, candidates| {
                assert_eq!(report.sample, 1);
                assert_eq!(candidates.len(), 2);
                for (expected, (index, line, state)) in candidates.iter().enumerate() {
                    assert_eq!(*index, expected);
                    assert!(std::ptr::eq(line.normal(), refreshed[expected].source()));
                    assert!(std::ptr::eq(line.shear(), &next[expected]));
                    assert_eq!(state.position_m(), before_states[expected].position_m());
                    assert_eq!(state.regime(), before_states[expected].regime());
                }
                Ok(report.clone())
            },
        );
    let report = admitted.unwrap();
    let replay_report = replay.apply_source_refreshes(&requests, &gate()).unwrap();
    assert_eq!(report, replay_report);
    assert_eq!(bits(owner.network()), native);
    assert_eq!(owner.last_frame().unwrap(), &frame);
    assert_eq!(owner.network().last_frame().unwrap(), &native_frame);
    assert_eq!(owner.compliance_m_per_n(), compliance.as_slice());
    assert_eq!(owner.config(), controls());
    assert_eq!(owner.samples_rendered(), 1);
    assert_eq!(report.energy_before_j, before_energy);
    assert_eq!(report.energy_after_j, owner.total_energy_j().unwrap());
    assert_eq!(
        report.network_energy_j,
        owner.network().total_energy_j().unwrap()
    );
    let mut expected_work = Dd::ZERO;
    let mut expected_before = Dd::ZERO;
    let mut expected_after = Dd::ZERO;
    let mut contact_scale = 0.0;
    for i in 0..2 {
        assert!(std::ptr::eq(
            owner.normal(i).unwrap(),
            refreshed[i].source()
        ));
        assert!(std::ptr::eq(owner.shear(i).unwrap(), &next[i]));
        assert_eq!(
            owner.state(i).unwrap().position_m(),
            before_states[i].position_m()
        );
        assert_eq!(owner.state(i).unwrap().regime(), before_states[i].regime());
        assert_eq!(owner.work_columns(i).unwrap(), &maps[i]);
        let delta = before_states[i].position_m()[0];
        let before = potential(delta, [0.01, 0.014][i], COLD);
        let after = potential(delta, [0.01, 0.014][i], TEMPERATURES[i]);
        let change = &report.contacts[i];
        assert_eq!(change.contact, i);
        assert_eq!(change.original_model_identity, normal[i].identity());
        assert_eq!(change.model_identity, refreshed[i].source().identity());
        close(change.surface_energy_before_j, before, 1e-30, 2e-13);
        close(change.surface_energy_after_j, after, 1e-30, 2e-13);
        close(change.parameter_work_j, after - before, 1e-30, 2e-13);
        assert!(change.energy_residual_j.abs() <= change.energy_tolerance_j);
        expected_work = expected_work + Dd::from_f64(after) - Dd::from_f64(before);
        expected_before = expected_before + Dd::from_f64(before);
        expected_after = expected_after + Dd::from_f64(after);
        contact_scale += before.abs().max(after.abs()).max((after - before).abs());
    }
    close(
        report.surface_energy_before_j,
        expected_before.to_f64(),
        1e-30,
        2e-13,
    );
    close(
        report.surface_energy_after_j,
        expected_after.to_f64(),
        1e-30,
        2e-13,
    );
    close(
        report.parameter_work_j,
        expected_work.to_f64(),
        1e-30,
        2e-13,
    );
    assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
    assert!(report.energy_tolerance_j < 130.0 * f64::EPSILON * contact_scale);
    for _ in 0..5 {
        let before = owner.total_energy_j().unwrap();
        let frame = owner.step(&[0.0; 6], &rates).unwrap().clone();
        assert_eq!(frame, *replay.step(&[0.0; 6], &rates).unwrap());
        assert_eq!(bits(owner.network()), bits(replay.network()));
        for (i, contact) in frame.contacts.iter().enumerate() {
            let a = contact.position_before_m;
            let b = contact.position_after_m;
            let area = mean_area(a[0], b[0]);
            let factor = 1.0 - 0.01 * (TEMPERATURES[i].interface_k - 300.0);
            close(contact.path.mean_real_area_m2, area, 1e-30, 2e-13);
            close(
                contact.path.static_capacity_n,
                3000.0 * factor * area,
                1e-30,
                2e-13,
            );
            close(
                contact.path.kinetic_capacity_n,
                2000.0 * factor * area,
                1e-30,
                2e-13,
            );
            close(
                contact.friction_dissipation_j,
                2000.0 * factor * area * (b[1] - a[1]).hypot(b[2] - a[2]),
                1e-30,
                2e-8,
            );
        }
        close(
            owner.total_energy_j().unwrap()
                + frame.network_dissipation_j
                + frame.friction_dissipation_j
                - frame.prescribed_work_j,
            before,
            1e-25,
            2e-8,
        );
    }
    // A later selected transition follows the retained hot source at the now
    // actual geometry. Reusing its initial receipt cannot represent this return
    // to the original temperature after real native motion.
    let position = owner.state(0).unwrap().position_m();
    let cooling = owner
        .normal(0)
        .unwrap()
        .refresh_at_held_indentation(position[0], f0.cards(), COLD)
        .unwrap();
    let cool_shear = owner
        .shear(0)
        .unwrap()
        .rebind_dmt_temperature(owner.normal(0).unwrap(), &cooling)
        .unwrap();
    let cool = [AdhesiveDmtSourceRefresh {
        contact: 0,
        refreshed: &cooling,
        next_shear: &cool_shear,
    }];
    let native = bits(owner.network());
    let frame = owner.last_frame().unwrap().clone();
    let unchanged_neighbour = owner.state(1);
    let report = owner.apply_source_refreshes(&cool, &gate()).unwrap();
    assert_eq!(
        report,
        replay.apply_source_refreshes(&cool, &gate()).unwrap()
    );
    assert_eq!(bits(owner.network()), native);
    assert_eq!(owner.last_frame().unwrap(), &frame);
    assert_eq!(owner.state(0).unwrap().position_m(), position);
    assert_eq!(owner.state(1), unchanged_neighbour);
    close(
        report.parameter_work_j,
        potential(position[0], 0.01, COLD) - potential(position[0], 0.01, TEMPERATURES[0]),
        1e-30,
        2e-13,
    );
    assert_eq!(
        owner.step(&[0.0; 6], &rates).unwrap(),
        replay.step(&[0.0; 6], &rates).unwrap()
    );
    assert_eq!(bits(owner.network()), bits(replay.network()));
}

#[test]
fn a_selected_refresh_retains_untouched_sources_and_contact_gap_and_separated_branches() {
    let fixture = thermal::Fixture::with_options(1e5, 0.01, false, 320.0);
    let normal = fixture.dmt(Z);
    let shear = shear(&normal);
    let indentation = [0.2 * Z, -0.4 * Z, -2.0 * Z];
    for selected in 0..3 {
        let mut owner = MultipleAdhesiveDmtModalSystem::new(
            rich_network(),
            indentation
                .iter()
                .map(|&d| contact(&normal, &shear, [d, 0.0, 0.0], columns()))
                .collect(),
            controls(),
            &gate(),
        )
        .unwrap();
        let native = bits(owner.network());
        let states: Vec<_> = (0..3).map(|i| owner.state(i).unwrap()).collect();
        let refreshed = normal
            .refresh_at_held_indentation(indentation[selected], fixture.cards(), TEMPERATURES[0])
            .unwrap();
        let next = shear.rebind_dmt_temperature(&normal, &refreshed).unwrap();
        let report = owner
            .apply_source_refreshes(
                &[AdhesiveDmtSourceRefresh {
                    contact: selected,
                    refreshed: &refreshed,
                    next_shear: &next,
                }],
                &gate(),
            )
            .unwrap();
        assert_eq!(bits(owner.network()), native);
        assert_eq!(owner.samples_rendered(), 0);
        assert!(owner.last_frame().is_none());
        assert_eq!(report.contacts.len(), 1);
        assert_eq!(report.contacts[0].contact, selected);
        assert_eq!(report.energy_after_j, owner.total_energy_j().unwrap());
        for i in 0..3 {
            assert_eq!(owner.state(i).unwrap().position_m(), states[i].position_m());
            assert_eq!(owner.state(i).unwrap().regime(), states[i].regime());
            if i != selected {
                assert_eq!(owner.state(i).unwrap(), states[i]);
                assert!(std::ptr::eq(owner.normal(i).unwrap(), &normal));
                assert!(std::ptr::eq(owner.shear(i).unwrap(), &shear));
            }
        }
        close(
            report.parameter_work_j,
            potential(indentation[selected], 0.01, TEMPERATURES[0])
                - potential(indentation[selected], 0.01, COLD),
            1e-30,
            2e-13,
        );
        if selected == 2 {
            assert_eq!(report.parameter_work_j, 0.0);
            assert_eq!(report.energy_tolerance_j, 0.0);
        }
    }
}

#[test]
fn locally_admissible_cooling_is_refused_by_the_complete_successor_curvature_screen() {
    let fixture = thermal::Fixture::with_options(1e5, 0.01, false, 320.0);
    let normal = fixture.dmt(Z);
    let shear = shear(&normal);
    let initial = [0.2 * Z, 0.0, 0.0];
    let base = make(&normal, &shear, [0.0; 3], initial, false, RATE);
    let a = -0.5 * normal.evaluate(0.0).unwrap().one_sided_tangents_n_per_m[0];
    let c = 0.25 * PI * normal.effective_radius_m() * shear.law().kinetic_shear_strength_pa();
    let eta = 0.5 * (a + a.hypot(2.0 * c));
    let scale = (0.38 / (eta * base.compliance_m_per_n()[0][0])).sqrt();
    let maps = || columns().map(|v| v.into_iter().map(|x| x * scale).collect());
    let local = base
        .compliance_m_per_n()
        .map(|r| r.map(|s| s * scale * scale));
    let mut owner = MultipleAdhesiveDmtModalSystem::new(
        network([0.0; 3], false, RATE),
        vec![
            contact(&normal, &shear, initial, maps()),
            contact(&normal, &shear, initial, maps()),
        ],
        controls(),
        &gate(),
    )
    .unwrap();
    let cooled = normal
        .refresh_at_held_indentation(
            initial[0],
            fixture.cards(),
            AdhesiveTemperatures {
                surface_a_k: 280.0,
                surface_b_k: 280.0,
                interface_k: 280.0,
            },
        )
        .unwrap();
    let cold_shear = shear.rebind_dmt_temperature(&normal, &cooled).unwrap();
    let line = DmtCoupledLoadLine::new(&normal, &shear, local, config()).unwrap();
    line.with_temperature_refresh(
        line.initialize(initial).unwrap(),
        &cooled,
        &cold_shear,
        &gate(),
    )
    .unwrap();
    let request = |contact| AdhesiveDmtSourceRefresh {
        contact,
        refreshed: &cooled,
        next_shear: &cold_shear,
    };
    let native = bits(owner.network());
    let states = [owner.state(0), owner.state(1)];
    for batch in [vec![request(0)], vec![request(1), request(0)]] {
        let mut called = false;
        let refused: Result<(), MultipleAdhesiveDmtModalError> = owner
            .apply_source_refreshes_with_admission(&batch, &gate(), |_, _| {
                called = true;
                Ok(())
            });
        assert!(matches!(
            refused,
            Err(MultipleAdhesiveDmtModalError::CoupledStability)
        ));
        assert!(!called);
        assert_eq!(bits(owner.network()), native);
        assert_eq!([owner.state(0), owner.state(1)], states);
        for i in 0..2 {
            assert!(std::ptr::eq(owner.normal(i).unwrap(), &normal));
            assert!(std::ptr::eq(owner.shear(i).unwrap(), &shear));
        }
    }
    let warmed = normal
        .refresh_at_held_indentation(initial[0], fixture.cards(), TEMPERATURES[0])
        .unwrap();
    let warm_shear = shear.rebind_dmt_temperature(&normal, &warmed).unwrap();
    owner
        .apply_source_refreshes(
            &[AdhesiveDmtSourceRefresh {
                contact: 0,
                refreshed: &warmed,
                next_shear: &warm_shear,
            }],
            &gate(),
        )
        .unwrap();
    assert_eq!(bits(owner.network()), native);
    assert_eq!(owner.state(1), states[1]);
}

#[test]
fn stale_foreign_duplicate_and_receiver_refusals_preserve_exact_batch_replay() {
    let f0 = thermal::Fixture::with_options(1e5, 0.01, false, 320.0);
    let f1 = thermal::Fixture::with_options(1e5, 0.014, false, 320.0);
    let normal = [f0.dmt(Z), f1.dmt(Z)];
    let shears = normal.each_ref().map(shear);
    let make_owner = || {
        MultipleAdhesiveDmtModalSystem::new(
            rich_network(),
            vec![
                contact(&normal[0], &shears[0], [0.2 * Z, 0.0, 0.0], columns()),
                contact(
                    &normal[1],
                    &shears[1],
                    [-0.4 * Z, 0.0, 0.0],
                    second_columns(),
                ),
            ],
            controls(),
            &gate(),
        )
        .unwrap()
    };
    let mut owner = make_owner();
    let mut replay = make_owner();
    let rates = [[0.0, 1e-3, -0.5e-3], [0.0; 3]];
    owner.step(&[0.0; 6], &rates).unwrap();
    replay.step(&[0.0; 6], &rates).unwrap();
    let original = [owner.state(0).unwrap(), owner.state(1).unwrap()];
    let native = bits(owner.network());
    let frame = owner.last_frame().unwrap().clone();
    let native_frame = owner.network().last_frame().unwrap().clone();
    let energy = owner.total_energy_j().unwrap();
    let unchanged = |system: &MultipleAdhesiveDmtModalSystem<'_>| {
        assert_eq!(bits(system.network()), native);
        assert_eq!(system.last_frame().unwrap(), &frame);
        assert_eq!(system.network().last_frame().unwrap(), &native_frame);
        assert_eq!(system.total_energy_j().unwrap(), energy);
        for i in 0..2 {
            assert_eq!(system.state(i).unwrap(), original[i]);
            assert!(std::ptr::eq(system.normal(i).unwrap(), &normal[i]));
            assert!(std::ptr::eq(system.shear(i).unwrap(), &shears[i]));
        }
    };
    let refreshed = [
        normal[0]
            .refresh_at_held_indentation(original[0].position_m()[0], f0.cards(), TEMPERATURES[0])
            .unwrap(),
        normal[1]
            .refresh_at_held_indentation(original[1].position_m()[0], f1.cards(), TEMPERATURES[1])
            .unwrap(),
    ];
    let next = [
        shears[0]
            .rebind_dmt_temperature(&normal[0], &refreshed[0])
            .unwrap(),
        shears[1]
            .rebind_dmt_temperature(&normal[1], &refreshed[1])
            .unwrap(),
    ];
    let good = [
        AdhesiveDmtSourceRefresh {
            contact: 0,
            refreshed: &refreshed[0],
            next_shear: &next[0],
        },
        AdhesiveDmtSourceRefresh {
            contact: 1,
            refreshed: &refreshed[1],
            next_shear: &next[1],
        },
    ];
    let stale = normal[1]
        .refresh_at_held_indentation(
            original[1].position_m()[0] + 0.01 * Z,
            f1.cards(),
            TEMPERATURES[1],
        )
        .unwrap();
    let foreign = normal[0]
        .refresh_at_held_indentation(original[1].position_m()[0], f0.cards(), TEMPERATURES[0])
        .unwrap();
    let foreign_shear = shears[0]
        .rebind_dmt_temperature(&normal[0], &foreign)
        .unwrap();
    let held_load = normal[0]
        .refresh_at_held_load(
            normal[0]
                .evaluate(original[0].position_m()[0])
                .unwrap()
                .normal_force_n,
            f0.cards(),
            TEMPERATURES[0],
        )
        .unwrap();
    let bad_second = [
        AdhesiveDmtSourceRefresh {
            contact: 1,
            refreshed: &stale,
            next_shear: &next[1],
        },
        AdhesiveDmtSourceRefresh {
            contact: 1,
            refreshed: &foreign,
            next_shear: &foreign_shear,
        },
        AdhesiveDmtSourceRefresh {
            contact: 1,
            refreshed: &refreshed[1],
            next_shear: &shears[1],
        },
    ];
    for bad in bad_second {
        assert!(matches!(
            owner.apply_source_refreshes(&[good[0], bad], &gate()),
            Err(MultipleAdhesiveDmtModalError::Contact {
                contact: 1,
                source: DmtCoupledError::StateMismatch
            })
        ));
        unchanged(&owner);
    }
    assert!(matches!(
        owner.apply_source_refreshes(
            &[AdhesiveDmtSourceRefresh {
                contact: 0,
                refreshed: &held_load,
                next_shear: &next[0],
            }],
            &gate()
        ),
        Err(MultipleAdhesiveDmtModalError::Contact {
            contact: 0,
            source: DmtCoupledError::StateMismatch
        })
    ));
    unchanged(&owner);
    for bad in [
        vec![],
        vec![good[0], good[0]],
        vec![good[0], good[1], good[0]],
        vec![AdhesiveDmtSourceRefresh {
            contact: 2,
            ..good[1]
        }],
    ] {
        assert!(owner.apply_source_refreshes(&bad, &gate()).is_err());
        unchanged(&owner);
    }
    let cancelled = gate();
    cancelled.request();
    assert!(matches!(
        owner.apply_source_refreshes(&good, &cancelled),
        Err(MultipleAdhesiveDmtModalError::Native(
            ModalCouplingError::Cancelled
        ))
    ));
    unchanged(&owner);
    let refused: Result<(), MultipleAdhesiveDmtModalError> = owner
        .apply_source_refreshes_with_admission(&good, &gate(), |report, candidates| {
            assert_eq!(report.contacts.len(), 2);
            assert_eq!(candidates.len(), 2);
            Err(MultipleAdhesiveDmtModalError::EvaluationBudget)
        });
    assert!(refused.is_err());
    unchanged(&owner);
    let cancelled = gate();
    let refused: Result<(), MultipleAdhesiveDmtModalError> = owner
        .apply_source_refreshes_with_admission(&good, &cancelled, |_, _| {
            cancelled.request();
            Ok(())
        });
    assert!(matches!(
        refused,
        Err(MultipleAdhesiveDmtModalError::Native(
            ModalCouplingError::Cancelled
        ))
    ));
    unchanged(&owner);
    let report = owner
        .apply_source_refreshes(&[good[1], good[0]], &gate())
        .unwrap();
    assert_eq!(
        report,
        replay.apply_source_refreshes(&good, &gate()).unwrap()
    );
    let accepted = [owner.state(0), owner.state(1)];
    let accepted_energy = owner.total_energy_j().unwrap();
    assert!(matches!(
        owner.apply_source_refreshes(&good, &gate()),
        Err(MultipleAdhesiveDmtModalError::Contact {
            contact: 0,
            source: DmtCoupledError::StateMismatch
        })
    ));
    assert_eq!([owner.state(0), owner.state(1)], accepted);
    assert_eq!(owner.total_energy_j().unwrap(), accepted_energy);
    assert_eq!(bits(owner.network()), native);
    assert_eq!(owner.last_frame().unwrap(), &frame);
    assert_eq!(
        owner.step(&[0.0; 6], &rates).unwrap(),
        replay.step(&[0.0; 6], &rates).unwrap()
    );
    assert_eq!(bits(owner.network()), bits(replay.network()));
}
