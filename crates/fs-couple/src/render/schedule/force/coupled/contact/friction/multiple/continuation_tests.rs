//! Branch continuation and root disclosure in the simultaneous friction owner
//! (IRP pexfp.152).
//!
//! Two light attachments, each a single mode driven at constant speed into
//! stick-slip, are bowed jointly. Before every accepted step the test rebuilds
//! each port's load line from the network's own prediction, including the
//! cross-port reaction committed in the same step, and checks that the frame's
//! disclosed root count and committed index match an independent dense scan of
//! that line. Three-root steps must be taken stuck while sticking and on the
//! high branch while sliding, for both Stribeck and hyperbolic laws.

use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
use crate::render::schedule::force::coupled::{CoupledModalSystem, ModalCouplingConfig};
use fs_math::c64::C64;

fn budgets() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 8,
        max_connections: 4,
        max_setup_terms: 4096,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6,
        maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-8,
    }
}
fn config() -> ModalFrictionConfig {
    ModalFrictionConfig {
        max_iterations: 96,
        maximum_force_n: 50.0,
        maximum_slip_speed_m_s: 100.0,
        force_absolute_tolerance_n: 1e-12,
        force_relative_tolerance: 1e-10,
        sticking_tolerance_m: 1e-12,
    }
}
fn model(hz: f64) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(
        48_000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: core::f64::consts::TAU * hz,
            damping_ratio: 0.01,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap()
}
fn attach(component: usize) -> ModalAttachment {
    // A light mode (1/sqrt(1e-4 kg)): a large per-step mobility, so a 1 N load
    // opens a three-root window.
    ModalAttachment {
        component,
        shapes: vec![100.0],
    }
}
fn system(law: FrictionLaw) -> MultiFrictionModalSystem {
    let gate = CancelGate::new_clock_free();
    let network =
        CoupledModalSystem::new(vec![model(180.0), model(260.0)], vec![], budgets(), &gate)
            .unwrap();
    let core = FrictionModalSystem::new(
        network,
        ModalFrictionPort {
            left: attach(0),
            right: None,
            law: law.clone(),
        },
        config(),
        &gate,
    )
    .unwrap();
    core.with_additional_ports(
        vec![(
            ModalFrictionPort {
                left: attach(1),
                right: None,
                law,
            },
            config(),
        )],
        MultiFrictionConfig {
            max_ports: 4,
            max_sweeps: 128,
            max_setup_terms: 1_000_000,
        },
        &gate,
    )
    .unwrap()
}

/// Sliding roots of g on (0, s] by a dense scan, each refined by bisection.
fn scan_roots(s: f64, m: f64, normal: f64, law: &FrictionLaw) -> Vec<f64> {
    let g = |v: f64| {
        let mu = if v == 0.0 {
            law.zero_speed_kinetic_coefficient().unwrap()
        } else {
            law.kinetic_coefficient(v).unwrap()
        };
        v - s + m * normal * mu
    };
    let cells = 50_000_u32;
    let mut roots = Vec::new();
    let mut prev = (0.0, g(0.0));
    for i in 1..=cells {
        let v = s * f64::from(i) / f64::from(cells);
        let cur = g(v);
        if (prev.1 < 0.0) != (cur < 0.0) {
            let (mut lo, mut hi, below) = (prev.0, v, prev.1 < 0.0);
            loop {
                let mid = lo + 0.5 * (hi - lo);
                if mid <= lo || mid >= hi {
                    break;
                }
                if (g(mid) < 0.0) == below {
                    lo = mid
                } else {
                    hi = mid
                }
            }
            roots.push(hi);
        }
        prev = (v, cur);
    }
    roots
}

fn exercise(law: FrictionLaw) -> [usize; 2] {
    let mut system = system(law.clone());
    let drives = [
        FrictionDrive {
            speed_m_s: 0.2,
            normal_force_n: 1.0,
        },
        FrictionDrive {
            speed_m_s: -0.15,
            normal_force_n: 1.0,
        },
    ];
    let static_mu = law.static_coefficient().unwrap();
    let mut three_root = [0_usize; 2];
    for _ in 0..4_000 {
        // Recompute each port's held-force load line from the prediction plus
        // the cross-port reaction committed in this same step.
        system
            .core
            .network
            .prepare_prediction(&[0.0, 0.0], None)
            .unwrap();
        let mut free = [0.0_f64; 2];
        for i in 0..2 {
            free[i] = drives[i].speed_m_s * system.sample_period_s()
                - dot(&system.points[i].column, &system.core.network.free_delta).unwrap();
        }
        let frame = system.step(&[0.0, 0.0], &drives).unwrap().clone();
        let p = 2;
        for i in 0..2 {
            let m = system.compliance[i * p + i] / system.sample_period_s();
            // The same step's cross-port reaction shifts this port's free slip.
            let v_free = (free[i] - system.compliance[i * p + (1 - i)] * system.tractions[1 - i])
                / system.sample_period_s();
            let port = &frame.ports[i];
            let roots = port
                .roots
                .expect("a loaded Stribeck/hyperbolic port discloses roots");
            let cone = v_free.abs() / m <= static_mu * drives[i].normal_force_n;
            let sliding = scan_roots(v_free.abs(), m, drives[i].normal_force_n, &law);
            assert_eq!(
                usize::from(roots.found),
                usize::from(cone) + sliding.len(),
                "port {i} intersections at v_free {v_free}"
            );
            assert!(roots.selected < roots.found);
            if port.regime == ModalFrictionRegime::Sticking {
                assert!(
                    cone && roots.selected == 0,
                    "port {i} stuck outside the cone"
                );
            } else {
                // The committed slide names a real sliding intersection (cone
                // first, then the sliding roots by speed). When three sliding
                // roots exist, sliding[1] is the unstable middle and is never
                // selected.
                let index = usize::from(roots.selected) - usize::from(cone);
                if sliding.len() == 3 {
                    assert_ne!(index, 1, "port {i} slide took the unstable middle root");
                }
                let expected = v_free.signum()
                    * law.kinetic_coefficient(sliding[index]).unwrap()
                    * drives[i].normal_force_n;
                assert!(
                    (port.traction_n - expected).abs() <= 1e-10 * expected.abs().max(1.0),
                    "port {i} traction {} vs selected root {expected}",
                    port.traction_n
                );
            }
            if roots.found == 3 {
                three_root[i] += 1;
            }
        }
    }
    three_root
}

#[test]
fn g1_multi_owner_discloses_roots_for_stribeck_and_hyperbolic() {
    for law in [
        FrictionLaw::Stribeck {
            static_mu: 1.2,
            kinetic_mu: 0.35,
            characteristic_speed: 0.05,
            viscous_per_speed: 0.0,
        },
        FrictionLaw::Hyperbolic {
            static_mu: 1.2,
            kinetic_mu: 0.35,
            characteristic_speed: 0.02,
        },
    ] {
        let three = exercise(law.clone());
        println!(
            "{{\"test\":\"multi_branch_continuation\",\"law\":\"{law:?}\",\"three_root_steps_per_port\":{three:?}}}"
        );
        assert!(
            three[0] > 0 && three[1] > 0,
            "both ports must visit a three-root window: {three:?}"
        );
    }
}

#[test]
fn g0_coulomb_multi_port_discloses_nothing() {
    let mut system = system(FrictionLaw::Coulomb {
        static_mu: 0.8,
        kinetic_mu: 0.4,
    });
    let frame = system
        .step(
            &[0.0, 0.0],
            &[
                FrictionDrive {
                    speed_m_s: 0.2,
                    normal_force_n: 1.0,
                },
                FrictionDrive {
                    speed_m_s: 0.2,
                    normal_force_n: 0.0,
                },
            ],
        )
        .unwrap();
    assert!(frame.ports.iter().all(|p| p.roots.is_none()));
}

fn residual_system(g: Vec<f64>) -> MultiFrictionModalSystem {
    let joint = system(FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 });
    let maps = joint.points.iter().map(|p| surfaces::geometry::FrictionGeometry {
        left: p.port.left.clone(), right: p.port.right.clone(),
    }).collect();
    joint.with_friction_residual_ports(maps, g, &CancelGate::new_clock_free()).unwrap()
}

#[test]
fn friction_residual_cross_terms_match_the_analytic_first_sticking_step() {
    // Independent retained modes; the omitted spectrum supplies the only cross-term.
    let g = vec![2e-5, 5e-6, 5e-6, 4e-5];
    let mut joint = residual_system(g.clone());
    let a: Vec<f64> = joint.compliance.iter().zip(&g).map(|(s,g)| s+g).collect();
    let drive = [FrictionDrive { speed_m_s: 0.0, normal_force_n: 1.0 },
        FrictionDrive { speed_m_s: 0.02, normal_force_n: 1.0 }];
    let distance = drive[1].speed_m_s*joint.sample_period_s();
    let determinant = a[0]*a[3]-a[1]*a[2];
    let expected = [-a[1]*distance/determinant, a[0]*distance/determinant];
    let f = joint.step(&[0.0,0.0], &drive).unwrap();
    for (port, force) in f.ports.iter().zip(expected) {
        assert_eq!(port.regime, ModalFrictionRegime::Sticking);
        assert!((port.traction_n-force).abs() < 1e-10);
    }
    let storage = 0.5*(g[0]*expected[0]*expected[0]
        +2.0*g[1]*expected[0]*expected[1]+g[3]*expected[1]*expected[1]);
    assert!((f.friction_residual_energy_j-storage).abs() < 1e-14);
    assert_eq!(f.friction_residual_energy_j, f.friction_residual_loss_j);
    assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
}

#[test]
fn friction_residual_and_material_drive_share_one_load_line_and_separate_losses() {
    let g = vec![2e-5, 5e-6, 5e-6, 4e-5];
    let gate = CancelGate::new_clock_free();
    let mut joint = residual_system(g.clone())
        .with_series_compliant_drives(vec![(0, 4e4, 0.4), (1, 7e4, 0.8)], &gate)
        .unwrap();
    let dt = joint.sample_period_s();
    // Independent two-contact static solve: retained response + string tail +
    // hair compliance. Cross-port coupling belongs only to the original string.
    let a = joint.compliance[0] + g[0] + 1.0 / (0.5 * 4e4 + 0.4 / dt);
    let b = joint.compliance[1] + g[1];
    let d = joint.compliance[3] + g[3] + 1.0 / (0.5 * 7e4 + 0.8 / dt);
    let distance = 0.02 * dt;
    let determinant = a * d - b * b;
    let expected = [-b * distance / determinant, a * distance / determinant];
    let drive = [FrictionDrive { speed_m_s: 0.0, normal_force_n: 1.0 },
        FrictionDrive { speed_m_s: 0.02, normal_force_n: 1.0 }];
    let frame = joint.step(&[0.0, 0.0], &drive).unwrap().clone();
    let mut hair_storage = 0.0;
    let mut hair_loss = 0.0;
    for (i, (k, c)) in [(4e4, 0.4), (7e4, 0.8)].into_iter().enumerate() {
        let stretch = expected[i] / (0.5 * k + c / dt);
        let stored = 0.5 * k * stretch * stretch;
        let heat = c * stretch * stretch / dt;
        assert_eq!(frame.ports[i].regime, ModalFrictionRegime::Sticking);
        assert!((frame.ports[i].traction_n - expected[i]).abs() < 1e-10);
        assert!((frame.ports[i].series_drive_dissipation_j - heat).abs() < 1e-14);
        assert_eq!(frame.ports[i].dissipation_j, 0.0);
        hair_storage += stored;
        hair_loss += heat;
    }
    let tail_storage = 0.5 * (g[0] * expected[0].powi(2)
        + 2.0 * g[1] * expected[0] * expected[1] + g[3] * expected[1].powi(2));
    assert!((frame.friction_residual_energy_j - tail_storage).abs() < 1e-14);
    assert_eq!(frame.friction_residual_energy_j, frame.friction_residual_loss_j);
    assert!((frame.friction_dissipation_j - hair_loss).abs() < 1e-14);
    let network_storage = joint.core.network.total_energy_j().unwrap();
    assert!((joint.total_energy_j().unwrap() - network_storage - hair_storage - tail_storage).abs() < 1e-14);
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);

    // A released interface must unload both stores rather than discard them.
    // The original series state and residual force carry the history through
    // reversal and release; numerical tail loss must never become hair heat.
    for i in 0..96 {
        let speed = if i < 32 { 0.02 } else { -0.03 };
        let normal = if i < 64 { 0.1 } else { 0.0 };
        let frame = joint.step(&[0.0, 0.0], &[
            FrictionDrive { speed_m_s: speed, normal_force_n: normal },
            FrictionDrive { speed_m_s: -speed, normal_force_n: normal },
        ]).unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert!(frame.ports.iter().all(|p| p.series_drive_dissipation_j >= 0.0));
    }
    assert_eq!(joint.last_frame().unwrap().friction_residual_energy_j, 0.0);
}

#[test]
fn friction_residual_history_closes_work_through_drive_reversal_and_release() {
    let g = vec![2e-5, -5e-6, -5e-6, 4e-5];
    let mut joint = residual_system(g.clone());
    let mut old = [0.0;2];
    let mut numerical = 0.0;
    for step in 0..600 {
        let speed = if step < 200 {0.02} else {-0.03};
        let normal = if step < 400 {0.1} else {0.0};
        let drives = [FrictionDrive {speed_m_s:speed,normal_force_n:normal},
            FrictionDrive {speed_m_s:0.5*speed,normal_force_n:normal}];
        let before: Vec<f64> = joint.components().iter().map(|m|m.states()[0].displacement_m_sqrt_kg).collect();
        let f = joint.step(&[0.0,0.0], &drives).unwrap().clone();
        let force = [f.ports[0].traction_n,f.ports[1].traction_n];
        let delta = [force[0]-old[0],force[1]-old[1]];
        let energy = |x:[f64;2]|0.5*(g[0]*x[0]*x[0]+2.0*g[1]*x[0]*x[1]+g[3]*x[1]*x[1]);
        assert!((f.friction_residual_energy_j-energy(force)).abs()<1e-14);
        assert!((f.friction_residual_loss_j-energy(delta)).abs()<1e-14);
        assert!(f.energy_residual_j.abs()<=f.energy_tolerance_j);
        for i in 0..2 {
            let displacement = 100.0*(joint.components()[i].states()[0].displacement_m_sqrt_kg-before[i]);
            let slip = drives[i].speed_m_s*joint.sample_period_s()-displacement-g[2*i]*delta[0]-g[2*i+1]*delta[1];
            assert!((f.ports[i].slip_distance_m-slip).abs()<1e-13);
        }
        numerical += f.friction_residual_loss_j;
        old = force;
    }
    assert!(numerical>0.0 && joint.last_frame().unwrap().friction_residual_energy_j==0.0);
    let accepted=joint.last_frame().unwrap().clone();
    let energy=joint.total_energy_j().unwrap();
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(joint.step_under_gate(&[0.0,0.0], &[FrictionDrive {speed_m_s:0.1,normal_force_n:1.0};2], &cancelled).is_err());
    assert_eq!(joint.last_frame(),Some(&accepted));
    assert_eq!(joint.total_energy_j().unwrap(),energy);
}

#[test]
fn friction_residual_refuses_nonpassive_matrices_and_moving_maps() {
    let gate=CancelGate::new_clock_free();
    for g in [vec![1e-5],vec![-1e-5,0.0,0.0,1e-5],vec![1e-5,2e-5,2e-5,1e-5],
        vec![1e-5,0.0,1e-6,1e-5],vec![f64::NAN,0.0,0.0,1e-5]] {
        let joint=system(FrictionLaw::Coulomb {static_mu:0.8,kinetic_mu:0.4});
        let maps=joint.points.iter().map(|p|surfaces::geometry::FrictionGeometry {
            left:p.port.left.clone(),right:p.port.right.clone(),
        }).collect();
        assert!(joint.with_friction_residual_ports(maps,g,&gate).is_err());
    }
    let joint=residual_system(vec![1e-5,0.0,0.0,1e-5]);
    assert!(joint.validate_friction_geometry_config(surfaces::geometry::FrictionGeometryConfig {
        max_update_terms:1_000_000,
    }).is_err());
}

// ----- Series compliant drive (bow-hair compliance, M4.4-W2) -----

fn stribeck() -> FrictionLaw {
    FrictionLaw::Stribeck {
        static_mu: 1.2,
        kinetic_mu: 0.35,
        characteristic_speed: 0.05,
        viscous_per_speed: 0.0,
    }
}
fn two_drives() -> [FrictionDrive; 2] {
    [
        FrictionDrive {
            speed_m_s: 0.2,
            normal_force_n: 1.0,
        },
        FrictionDrive {
            speed_m_s: -0.15,
            normal_force_n: 1.0,
        },
    ]
}

#[test]
fn g3_soft_series_drive_shifts_the_motion_and_keeps_the_energy_ledger() {
    // A soft, damped hair changes the motion relative to the rigid drive, yet
    // every joint step still closes the energy balance: the series spring energy
    // and the dashpot loss are both in the census, so the residual stays within
    // the frame's own tolerance (a step that did not balance would refuse).
    let gate = CancelGate::new_clock_free();
    let mut soft = system(stribeck())
        .with_series_compliant_drives(vec![(0, 2.0e4, 8.0), (1, 2.0e4, 8.0)], &gate)
        .unwrap();
    let mut rigid = system(stribeck());
    let drives = two_drives();
    let mut differed = false;
    for _ in 0..2_000 {
        let s = soft.step(&[0.0, 0.0], &drives).unwrap().clone();
        let r = rigid.step(&[0.0, 0.0], &drives).unwrap().clone();
        assert!(
            s.energy_residual_j.abs() <= s.energy_tolerance_j,
            "soft-drive energy residual {} exceeds tolerance {}",
            s.energy_residual_j,
            s.energy_tolerance_j
        );
        for i in 0..2 {
            if (s.ports[i].traction_n - r.ports[i].traction_n).abs()
                > 1e-4 * r.ports[i].traction_n.abs().max(1.0)
            {
                differed = true;
            }
        }
    }
    assert!(
        differed,
        "a soft hair must change the motion vs the rigid drive"
    );
}

#[test]
fn g3b_undamped_soft_series_drive_balances_with_spring_energy_only() {
    // Undamped hair: no dashpot loss, so the balance rests on the series spring
    // energy alone being accounted in the ledger.
    let gate = CancelGate::new_clock_free();
    let mut soft = system(stribeck())
        .with_series_compliant_drives(vec![(0, 1.0e4, 0.0), (1, 1.0e4, 0.0)], &gate)
        .unwrap();
    let drives = two_drives();
    for _ in 0..1_500 {
        let s = soft.step(&[0.0, 0.0], &drives).unwrap().clone();
        assert!(s.energy_residual_j.abs() <= s.energy_tolerance_j);
    }
}

#[test]
fn g4_series_drive_binding_refusals() {
    let gate = CancelGate::new_clock_free();
    // Unknown port.
    assert!(
        system(stribeck())
            .with_series_compliant_drives(vec![(9, 1.0e5, 1.0)], &gate)
            .is_err()
    );
    // Two bindings on one port.
    assert!(
        system(stribeck())
            .with_series_compliant_drives(vec![(0, 1.0e5, 1.0), (0, 1.0e5, 1.0)], &gate)
            .is_err()
    );
    // Non-physical stiffness.
    assert!(
        system(stribeck())
            .with_series_compliant_drives(vec![(0, -1.0, 1.0)], &gate)
            .is_err()
    );
    // After a sample has been rendered.
    let mut rendered = system(stribeck());
    rendered.step(&[0.0, 0.0], &two_drives()).unwrap();
    assert!(
        rendered
            .with_series_compliant_drives(vec![(0, 1.0e5, 1.0)], &gate)
            .is_err()
    );
}
