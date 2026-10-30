//! Finite lever arms, actual elastic network and joint all-site transaction.
use super::*;
use fs_couple::cohesive_contact::patch::CohesivePatchConfig;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::cohesive::patch::{
    CohesiveModalSite, CohesivePatchModalConfig, CohesivePatchModalError, CohesivePatchModalSystem,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError,
};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};
const MA: f64 = 2.0;
const MB: f64 = 3.0;
const IA: f64 = 0.003;
const IB: f64 = 0.004;
const LOCATIONS: [f64; 3] = [-0.02, 0.0, 0.02];

fn bar_columns(x: f64) -> [Vec<f64>; 3] {
    let mut columns: [Vec<f64>; 3] = core::array::from_fn(|_| vec![0.0; 8]);
    for i in 0..3 {
        columns[i][i] = -1.0 / det::sqrt(MA);
        columns[i][4 + i] = 1.0 / det::sqrt(MB);
    }
    columns[0][3] = x / det::sqrt(IA);
    columns[0][7] = -x / det::sqrt(IB);
    columns
}

fn bars(relative_velocity: [f64; 4], rate: u32, ceiling: f64, elastic: bool) -> CoupledModalSystem {
    let models = [(MA, IA), (MB, IB)]
        .into_iter()
        .enumerate()
        .map(|(side, (mass, inertia))| {
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                rate,
                vec![
                    ModalAcousticMode {
                        angular_frequency_rad_s: 0.0,
                        damping_ratio: 0.0,
                        pressure_per_modal_velocity: C64::ZERO
                    };
                    4
                ],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap();
            let states: Vec<_> = (0..4)
                .map(|i| {
                    let (own, other) = if i == 3 {
                        (inertia, if side == 0 { IB } else { IA })
                    } else {
                        (mass, if side == 0 { MB } else { MA })
                    };
                    ModalAcousticState {
                        displacement_m_sqrt_kg: 0.0,
                        velocity_m_sqrt_kg_per_s: det::sqrt(own)
                            * (if side == 0 {
                                -other / (own + other)
                            } else {
                                other / (own + other)
                            })
                            * relative_velocity[i],
                    }
                })
                .collect();
            model.restore_states(&states).unwrap();
            model
        })
        .collect();
    CoupledModalSystem::new(
        models,
        if elastic {
            vec![ModalConnection {
                left: ModalAttachment {
                    component: 0,
                    shapes: vec![
                        0.7 / det::sqrt(MA),
                        0.3 / det::sqrt(MA),
                        -0.2 / det::sqrt(MA),
                        0.002 / det::sqrt(IA),
                    ],
                },
                right: ModalAttachment {
                    component: 1,
                    shapes: vec![
                        0.7 / det::sqrt(MB),
                        0.3 / det::sqrt(MB),
                        -0.2 / det::sqrt(MB),
                        0.002 / det::sqrt(IB),
                    ],
                },
                stiffness_n_m: 8e6,
                damping_n_s_m: 0.0,
                rest_extension_m: 0.0,
            }]
        } else {
            vec![]
        },
        ModalCouplingConfig {
            max_modes: 8,
            max_connections: 1,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: ceiling,
            maximum_abs_pressure_pa: 1000.0,
            maximum_abs_connection_force_n: 1e5,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-10,
        },
        &CancelGate::new(),
    )
    .unwrap()
}

fn patch<'a>(
    source: &'a ResolvedBilinearCohesive,
    ceiling: f64,
    rate: u32,
) -> CohesivePatchModalSystem<'a> {
    patch_with_connection(source, ceiling, rate, true)
}
fn patch_with_connection<'a>(
    source: &'a ResolvedBilinearCohesive,
    ceiling: f64,
    rate: u32,
    elastic: bool,
) -> CohesivePatchModalSystem<'a> {
    let areas = [0.25, 0.4, 0.35];
    let sites = LOCATIONS
        .iter()
        .enumerate()
        .map(|(p, &x)| {
            CohesiveModalSite::from_source(
                source,
                receiving(source),
                AREA * areas[p],
                bar_columns(x),
            )
            .unwrap()
        })
        .collect();
    CohesivePatchModalSystem::new(
        bars([0.008, 0.01, -0.005, -2.0], rate, ceiling, elastic),
        sites,
        CohesivePatchModalConfig {
            patch: CohesivePatchConfig {
                max_points: 8,
                contact: local_config().contact,
            },
            max_setup_terms: 100_000,
            force_absolute_tolerance_n: 1e-9,
            force_relative_tolerance: 1e-9,
        },
        &CancelGate::new(),
    )
    .unwrap()
}
fn enable<'a>(
    s: &mut CohesivePatchModalSystem<'a>,
    source: &'a ResolvedBilinearCohesive,
    friction: &'a ResidualCoulombSource,
) {
    let saved = (
        states(s),
        s.points().to_vec(),
        s.samples_rendered(),
        s.last_frame().copied(),
        s.last_point_frames().map(|x| x.to_vec()),
    );
    s.enable_frictional_damage(
        &[friction; 3],
        &[receiving(source); 3],
        &[ReclosureActivation::AtRest; 3],
        config(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(
        (
            states(s),
            s.points().to_vec(),
            s.samples_rendered(),
            s.last_frame().copied(),
            s.last_point_frames().map(|x| x.to_vec())
        ),
        saved
    );
    for p in 0..3 {
        assert!(core::ptr::eq(s.sites()[p].source().unwrap(), source));
        assert!(core::ptr::eq(
            s.damage_friction_source(p).unwrap(),
            friction
        ));
    }
}
fn states(s: &CohesivePatchModalSystem<'_>) -> Vec<ModalAcousticState> {
    s.components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}
fn momentum(s: &CohesivePatchModalSystem<'_>) -> [f64; 4] {
    core::array::from_fn(|i| {
        det::sqrt(if i == 3 { IA } else { MA })
            * s.components()[0].states()[i].velocity_m_sqrt_kg_per_s
            + det::sqrt(if i == 3 { IB } else { MB })
                * s.components()[1].states()[i].velocity_m_sqrt_kg_per_s
    })
}

#[test]
fn g1_original_rank_deficient_elastic_patch_peels_and_slides_with_force_torque_and_work_closure() {
    let (source, friction) = sources();
    let mut s = patch(&source, 1.0, 5000);
    let matrix = s.compliance_m_per_n().to_vec();
    // Nine site coordinates share four relative body coordinates. The complete
    // reciprocal response has actual inter-site normal/tangent coupling.
    assert!(matrix[(3 * 0) * 9 + 3 * 2 + 1].abs() > 1e-12);
    assert!(matrix[(3 * 0 + 1) * 9 + 3 * 2 + 2].abs() > 1e-12);
    enable(&mut s, &source, &friction);
    assert_eq!(s.compliance_m_per_n(), matrix);
    let initial = s.total_energy_j().unwrap();
    let original_momentum = momentum(&s);
    let mut loss = 0.0;
    let mut work = 0.0;
    let mut partial_sliding = false;
    let mut edge_ahead = false;
    for _ in 0..250 {
        let old = s.points().to_vec();
        let frame = *s.step(&[0.0; 8], &[[0.0; 3]; 3]).unwrap();
        let observations = s.last_point_frames().unwrap();
        let mut friction_sum = 0.0;
        for p in 0..3 {
            let a = observations[p];
            let extra = a.damage_friction.unwrap();
            let (material, mean, weighted) = old[p]
                .trial_over_time_with_damage_weighted_compression(
                    a.jump_m,
                    s.sample_period_s(),
                    1e-8,
                )
                .unwrap();
            assert_eq!(s.points()[p], material.equilibrium.point);
            assert!(a.damage >= old[p].damage());
            close(
                extra.damaged_compression_force_n,
                a.area_m2 * weighted,
                1e-11,
            );
            close(
                a.cohesive_energy_j,
                a.area_m2 * s.points()[p].stored_energy_j_m2(),
                1e-12,
            );
            if extra.regime == DamageFrictionRegime::Sliding {
                let slip =
                    (a.jump_m[1] - old[p].jump_m()[1]).hypot(a.jump_m[2] - old[p].jump_m()[2]);
                close(
                    extra.force_n[0].hypot(extra.force_n[1]),
                    0.3 * a.area_m2 * weighted,
                    1e-10,
                );
                close(
                    extra.dissipation_j,
                    0.3 * a.area_m2 * weighted * slip,
                    1e-12,
                );
                for axis in 1..3 {
                    close(
                        a.mean_force_n[axis],
                        a.area_m2 * mean[axis] + extra.force_n[axis - 1],
                        1e-9,
                    );
                }
                partial_sliding |= a.damage > 0.0 && a.damage < 1.0 && extra.dissipation_j > 0.0;
            }
            friction_sum += extra.dissipation_j;
            assert!(extra.work_defect_j.abs() < 1e-12);
        }
        close(frame.friction_dissipation_j, friction_sum, 1e-15);
        edge_ahead |= observations[2].damage > observations[0].damage + 0.1;
        loss += frame.fracture_dissipation_j
            + frame.viscous_dissipation_j
            + frame.friction_dissipation_j
            + frame.network_dissipation_j;
        work += frame.prescribed_work_j + frame.external_work_j;
        close(s.total_energy_j().unwrap() + loss, initial + work, 5e-11);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        for i in 0..4 {
            close(momentum(&s)[i], original_momentum[i], 1e-13);
        }
        assert!(frame.material_evaluations <= config().friction.contact.max_material_evaluations);
    }
    assert!(partial_sliding && edge_ahead);
    assert!(s.points()[2].failed());
    let history = s.points()[2].dissipated_j_m2();
    let close_rate = (-2e-4 - s.points()[2].jump_m()[0]) / s.sample_period_s();
    let frame = *s.step(&[0.0; 8], &[[close_rate, 0.05, -0.02]; 3]).unwrap();
    assert!(s.points()[2].jump_m()[0] < 0.0);
    assert!(s.points()[2].failed());
    assert_eq!(s.points()[2].dissipated_j_m2(), history);
    assert_eq!(
        s.last_point_frames().unwrap()[2].fracture_dissipation_j,
        0.0
    );
    assert!(frame.friction_dissipation_j > 0.0);
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
}

#[test]
fn g4_all_sites_sources_body_and_clock_roll_back_on_late_refusal_and_retry_exactly() {
    let (source, friction) = sources();
    let mut s = patch(&source, 0.02, 5000);
    let saved = (
        states(&s),
        s.points().to_vec(),
        s.samples_rendered(),
        s.last_frame().copied(),
    );
    let mut receiving_sites = [receiving(&source); 3];
    receiving_sites[2].swap(0, 1);
    assert!(
        s.enable_frictional_damage(
            &[&friction; 3],
            &receiving_sites,
            &[ReclosureActivation::AtRest; 3],
            config(),
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(
        (
            states(&s),
            s.points().to_vec(),
            s.samples_rendered(),
            s.last_frame().copied()
        ),
        saved
    );
    assert!((0..3).all(|p| s.damage_friction_source(p).is_none()));
    enable(&mut s, &source, &friction);
    let mut replay = patch(&source, 0.02, 5000);
    enable(&mut replay, &source, &friction);
    for _ in 0..140 {
        s.step(&[0.0; 8], &[[0.0; 3]; 3]).unwrap();
        replay.step(&[0.0; 8], &[[0.0; 3]; 3]).unwrap();
    }
    let saved = (
        states(&s),
        s.points().to_vec(),
        s.samples_rendered(),
        s.last_frame().copied(),
        s.last_point_frames().unwrap().to_vec(),
        (0..3)
            .map(|p| s.damage_friction_occupancy(p))
            .collect::<Vec<_>>(),
    );
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        s.step_under_gate(&[0.0; 8], &[[0.0; 3]; 3], &cancelled)
            .is_err()
    );
    let refused = s.step(&[0.0; 8], &[[-40.0, 0.1, 0.05]; 3]);
    assert!(
        matches!(
            refused,
            Err(CohesivePatchModalError::Network(
                ModalCouplingError::Budget {
                    what: "total energy including finite cohesive patch",
                    ..
                }
            ))
        ),
        "{refused:?}"
    );
    assert_eq!(
        (
            states(&s),
            s.points().to_vec(),
            s.samples_rendered(),
            s.last_frame().copied(),
            s.last_point_frames().unwrap().to_vec(),
            (0..3)
                .map(|p| s.damage_friction_occupancy(p))
                .collect::<Vec<_>>()
        ),
        saved
    );
    assert_eq!(
        *s.step(&[0.0; 8], &[[0.0; 3]; 3]).unwrap(),
        *replay.step(&[0.0; 8], &[[0.0; 3]; 3]).unwrap()
    );
    assert_eq!(states(&s), states(&replay));
    assert_eq!(s.points(), replay.points());
    assert_eq!(s.last_point_frames(), replay.last_point_frames());
}

// On the explicitly unconnected free-coordinate network the original free
// propagation is exactly v*dt. This constructs prescribed apparatus velocities
// for a known force/jump pair using the COMPLETE native patch compliance.
fn drive_for(
    s: &CohesivePatchModalSystem<'_>,
    target: &[[f64; 3]],
    force: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    let free = free_for(s.compliance_m_per_n(), target, force);
    let modes = states(s);
    (0..3)
        .map(|p| {
            core::array::from_fn(|i| {
                let velocity = s.sites()[p].work_columns()[i]
                    .iter()
                    .zip(&modes)
                    .map(|(b, q)| b * q.velocity_m_sqrt_kg_per_s)
                    .sum::<f64>();
                (free[p][i] - s.points()[p].jump_m()[i]) / s.sample_period_s() - velocity
            })
        })
        .collect()
}

#[test]
fn g1_zero_net_patch_slip_keeps_each_moving_endpoints_kinetic_occupancy() {
    let (source, friction) = sources();
    let mut s = patch_with_connection(&source, 1.0, 5000, false);
    for normal in [4e-4, -1e-4] {
        let target = vec![[normal, 0.0, 0.0]; 3];
        let forces: Vec<_> = (0..3)
            .map(|p| {
                let (_, mean) = s.points()[p]
                    .trial_over_time_with_mean_traction(target[p], s.sample_period_s(), 1e-8)
                    .unwrap();
                mean.map(|x| x * s.sites()[p].area_m2())
            })
            .collect();
        let drive = drive_for(&s, &target, &forces);
        s.step(&[0.0; 8], &drive).unwrap();
    }
    assert!(
        s.points()
            .iter()
            .all(|p| p.damage() > 0.0 && p.damage() < 1.0 && p.jump_m()[0] < 0.0)
    );
    enable(&mut s, &source, &friction);
    let old: Vec<_> = s.points().iter().map(|p| p.jump_m()).collect();
    let forces: Vec<_> = (0..3)
        .map(|p| {
            let pressure = -s.sites()[p].area_m2() * COMPRESSION * old[p][0];
            [-pressure, 0.2 * s.points()[p].damage() * pressure, 0.0]
        })
        .collect();
    let drive = drive_for(&s, &old, &forces);
    let frame = *s.step(&[0.0; 8], &drive).unwrap();
    assert_eq!(frame.friction_dissipation_j, 0.0);
    for p in 0..3 {
        let observation = s.last_point_frames().unwrap()[p];
        assert_eq!(
            observation.damage_friction.unwrap().regime,
            DamageFrictionRegime::Sticking
        );
        close(observation.jump_m[1], old[p][1], 1e-14);
        assert_eq!(
            s.damage_friction_occupancy(p),
            Some(DamageFrictionRegime::Sliding)
        );
    }
    // The next zero-slip demand would fit the static cone but exceeds kinetic
    // capacity. Actual moving endpoints must keep the latter at EVERY site.
    let current: Vec<_> = s.points().iter().map(|p| p.jump_m()).collect();
    let forces: Vec<_> = (0..3)
        .map(|p| {
            let pressure = -s.sites()[p].area_m2() * COMPRESSION * current[p][0];
            [-pressure, 0.45 * s.points()[p].damage() * pressure, 0.0]
        })
        .collect();
    let drive = drive_for(&s, &current, &forces);
    let frame = *s.step(&[0.0; 8], &drive).unwrap();
    assert!(frame.friction_dissipation_j > 0.0);
    for p in 0..3 {
        let extra = s.last_point_frames().unwrap()[p].damage_friction.unwrap();
        assert_eq!(extra.regime, DamageFrictionRegime::Sliding);
        close(
            extra.force_n[0],
            0.3 * extra.damaged_compression_force_n,
            1e-10,
        );
    }
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
}
