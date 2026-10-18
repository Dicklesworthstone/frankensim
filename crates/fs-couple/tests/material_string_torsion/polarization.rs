//! G1/G3 physical two-plane response on the SAME material/body/friction path.
//! Sources are the parent's explicit synthetic material fixtures, not measured
//! violin specimens. Rotation and work checks test wiring, not timbre realism.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::force::{
    StringForceCoupling, StringForcePort,
};
use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::TransverseDirections;

fn directions(bow: [f64; 2], bridge: [f64; 2]) -> TransverseDirections {
    TransverseDirections::try_new(bow, bridge).unwrap()
}
fn polarized(
    cx: &Cx<'_>,
    gate: &CancelGate,
    input: MaterialStringSource,
    directions: TransverseDirections,
) -> ReciprocalStringBank {
    ReciprocalStringBank::from_material_plate(
        cx,
        vec![input.with_transverse_polarizations(directions).unwrap()],
        gas(),
        body(1),
        law(),
        options(),
        gate,
    )
    .unwrap()
}
fn near(a: f64, b: f64, floor: f64) {
    assert!(
        (a - b).abs() <= 2e-9 * a.abs().max(b.abs()) + floor,
        "{a:e} != {b:e}"
    );
}

#[test]
fn both_planes_keep_the_same_material_mass_stiffness_loss_and_original_torsion() {
    with_cx(|cx, gate| {
        let mut old = radiating(cx, gate, 80e9);
        let mut bank = polarized(
            cx,
            gate,
            source(80e9, 1000.0).with_torsion("G", "eta_s", 3).unwrap(),
            directions([1.0, 0.0], [1.0, 0.0]),
        );
        let modes = bank.system().components()[0].modes();
        assert_eq!(modes.len(), 11);
        assert_eq!(&modes[..4], &modes[4..8]);
        assert_eq!(&modes[..4], &old.system().components()[0].modes()[..4]);
        assert_eq!(&modes[8..], &old.system().components()[0].modes()[4..]);
        assert_eq!(
            bank.cards()[0].linear_density_kg_m,
            old.cards()[0].linear_density_kg_m
        );
        let row = &bank.system().port().left.shapes;
        assert_eq!(&row[..4], &old.system().port().left.shapes[..4]);
        assert!(row[4..8].iter().all(|v| *v == 0.0));
        assert_eq!(&row[8..], &old.system().port().left.shapes[4..]);
        let MaterialStringSource::WithPolarizations(p) = &bank.source().strings()[0].source else {
            panic!("missing original two-plane source")
        };
        assert!(matches!(p.source(), MaterialStringSource::WithTorsion(_)));
        let old_zero = vec![0.0; old.system().mode_count()];
        let new_zero = vec![0.0; bank.system().mode_count()];
        for _ in 0..128 {
            let a = old.step(&old_zero, drive()).unwrap().observer_pressure_pa;
            let b = bank.step(&new_zero, drive()).unwrap().observer_pressure_pa;
            near(a, b, 1e-13);
        }
        assert!(
            bank.system().components()[0].states()[4..8]
                .iter()
                .all(|s| s.displacement_m_sqrt_kg == 0.0 && s.velocity_m_sqrt_kg_per_s == 0.0)
        );
    });
}

#[test]
fn common_transverse_basis_rotation_preserves_motion_sound_and_friction_work() {
    with_cx(|cx, gate| {
        let src = source(80e9, 1000.0).with_torsion("G", "eta_s", 3).unwrap();
        let mut a = polarized(cx, gate, src.clone(), directions([1.0, 0.0], [0.6, 0.8]));
        // Rotate ALL physical directions by the same orthogonal matrix.
        let mut b = polarized(cx, gate, src, directions([0.6, 0.8], [-0.28, 0.96]));
        let zero = vec![0.0; a.system().mode_count()];
        let mut second_plane = 0.0_f64;
        for n in 0..256 {
            let drive = if n < 128 {
                FrictionDrive {
                    speed_m_s: 0.01,
                    normal_force_n: 1.0,
                }
            } else {
                FrictionDrive {
                    speed_m_s: 0.0,
                    normal_force_n: 0.0,
                }
            };
            let fa = *a.step(&zero, drive).unwrap();
            let fb = *b.step(&zero, drive).unwrap();
            near(fa.observer_pressure_pa, fb.observer_pressure_pa, 1e-12);
            near(fa.drive_work_j, fb.drive_work_j, 1e-17);
            near(fa.friction_dissipation_j, fb.friction_dissipation_j, 1e-17);
            assert_eq!(fa.regime, fb.regime);
            for f in [fa, fb] {
                assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
                assert!(f.friction_dissipation_j >= 0.0);
            }
            let qa = a.system().components()[0].states();
            let qb = b.system().components()[0].states();
            for k in 0..4 {
                let x = qa[k].displacement_m_sqrt_kg;
                let y = qa[4 + k].displacement_m_sqrt_kg;
                near(qb[k].displacement_m_sqrt_kg, 0.6 * x - 0.8 * y, 1e-14);
                near(qb[4 + k].displacement_m_sqrt_kg, 0.8 * x + 0.6 * y, 1e-14);
                second_plane = second_plane.max(y.abs());
            }
            for k in 8..11 {
                near(
                    qa[k].displacement_m_sqrt_kg,
                    qb[k].displacement_m_sqrt_kg,
                    1e-14,
                );
            }
        }
        assert!(
            second_plane > 1e-10,
            "oblique bridge did not return force into the unbowed plane"
        );
        let initial = a.system().total_energy_j().unwrap();
        for _ in 0..128 {
            let f = a
                .step(
                    &zero,
                    FrictionDrive {
                        speed_m_s: 0.0,
                        normal_force_n: 0.0,
                    },
                )
                .unwrap();
            assert_eq!(f.drive_work_j, 0.0);
        }
        assert!(a.system().total_energy_j().unwrap() <= initial + 1e-12);
    });
}

#[test]
fn physical_bow_direction_changes_bridge_sound_and_a_decoupled_plane_stays_silent() {
    with_cx(|cx, gate| {
        let mut parallel = polarized(
            cx,
            gate,
            source(80e9, 1000.0),
            directions([1.0, 0.0], [1.0, 0.0]),
        );
        let mut oblique = polarized(
            cx,
            gate,
            source(80e9, 1000.0),
            directions([0.6, 0.8], [1.0, 0.0]),
        );
        let mut perpendicular = polarized(
            cx,
            gate,
            source(80e9, 1000.0),
            directions([0.0, 1.0], [1.0, 0.0]),
        );
        let zero = vec![0.0; parallel.system().mode_count()];
        let (mut audible, mut changed) = (0.0_f64, 0.0_f64);
        for _ in 0..256 {
            let a = parallel.step(&zero, drive()).unwrap().observer_pressure_pa;
            let b = oblique.step(&zero, drive()).unwrap().observer_pressure_pa;
            let dark = perpendicular
                .step(&zero, drive())
                .unwrap()
                .observer_pressure_pa;
            assert_eq!(
                dark, 0.0,
                "rank-one bridge invented a perpendicular reaction"
            );
            audible = audible.max(a.abs());
            changed = changed.max((a - b).abs());
        }
        assert!(audible > 1e-10 && changed > 1e-10);
        assert!(
            perpendicular.system().components()[0].states()[4..8]
                .iter()
                .any(|s| s.velocity_m_sqrt_kg_per_s.abs() > 1e-8)
        );
        assert!(
            oblique.system().components()[0].states()[..8]
                .iter()
                .all(|s| s.displacement_m_sqrt_kg.is_finite())
        );
    });
}

#[test]
fn centreline_and_surface_loads_share_both_planes_without_duplicating_torque_or_body_force() {
    with_cx(|cx, gate| {
        let source = source(80e9, 1000.0).with_torsion("G", "eta_s", 3).unwrap();
        let bank = polarized(cx, gate, source, directions([0.6, 0.8], [1.0, 0.0]));
        let surface = bank.system().port().left.shapes.clone();
        let bank = bank
            .with_additional_bows(
                vec![],
                MultiFrictionConfig {
                    max_ports: 4,
                    max_sweeps: 128,
                    max_setup_terms: 65_536,
                },
                gate,
            )
            .unwrap();
        let ports = [
            StringForceCoupling::Centreline,
            StringForceCoupling::Surface,
        ]
        .map(|coupling| StringForcePort {
            location: StringBowLocation {
                string: 0,
                station_fraction: 0.2,
            },
            coupling,
        });
        let columns = bank.string_force_columns(&ports, 4096, gate).unwrap();
        assert_eq!(&columns[1][..11], surface);
        assert_eq!(&columns[0][..8], &surface[..8]);
        assert!(columns[0][8..].iter().all(|x| *x == 0.0));
        assert!(columns[1][11..].iter().all(|x| *x == 0.0));
        for k in 0..4 {
            near(columns[0][k] * 0.8, columns[0][4 + k] * 0.6, 1e-14);
        }
        let squared: f64 = columns[0][..8].iter().map(|x| x * x).sum();
        let card = &bank.cards()[0];
        let expected: f64 = (1..=4)
            .map(|n| {
                let value = (n as f64 * core::f64::consts::PI * 0.2).sin();
                value * value / (0.5 * card.linear_density_kg_m * card.length_m)
            })
            .sum();
        near(squared, expected, 1e-12);
    });
}

#[test]
fn torsion_and_polarization_construction_order_is_identical_and_invalid_geometry_refuses() {
    for direction in [
        [0.0, 0.0],
        [2.0, 0.0],
        [f64::NAN, 0.0],
        [f64::INFINITY, 0.0],
    ] {
        assert!(TransverseDirections::try_new(direction, [1.0, 0.0]).is_err());
        assert!(TransverseDirections::try_new([1.0, 0.0], direction).is_err());
    }
    let d = directions([0.6, 0.8], [1.0, 0.0]);
    assert!(
        source(80e9, 1000.0)
            .with_transverse_polarizations(d)
            .unwrap()
            .with_transverse_polarizations(d)
            .is_err()
    );
    with_cx(|cx, gate| {
        let a = source(80e9, 1000.0)
            .with_torsion("G", "eta_s", 3)
            .unwrap()
            .with_transverse_polarizations(d)
            .unwrap();
        let b = source(80e9, 1000.0)
            .with_transverse_polarizations(d)
            .unwrap()
            .with_torsion("G", "eta_s", 3)
            .unwrap();
        let a = raw(cx, gate, a, options()).unwrap();
        let b = raw(cx, gate, b, options()).unwrap();
        assert_eq!(
            a.system().components()[0].modes(),
            b.system().components()[0].modes()
        );
        assert_eq!(
            a.system().port().left.component,
            b.system().port().left.component
        );
        assert_eq!(a.system().port().left.shapes, b.system().port().left.shapes);
        assert_eq!(a.system().port().law, b.system().port().law);
        let mut cap = options();
        cap.network.max_modes = 8;
        assert!(
            raw(
                cx,
                gate,
                source(80e9, 1000.0)
                    .with_transverse_polarizations(d)
                    .unwrap(),
                cap
            )
            .is_err()
        );
        assert!(
            source(80e9, 1000.0)
                .with_transverse_polarizations(d)
                .unwrap()
                .into_nonlinear(
                    cx,
                    gas(),
                    48_000,
                    0.1,
                    vec![fs_couple::modal_acoustic_time::ModalAcousticState::default(); 8],
                    gate
                )
                .is_err()
        );
    });
}

#[test]
fn two_planes_and_rotational_bridge_share_one_covariant_passive_response() {
    with_cx(|cx, gate| {
        let make = |d| {
            let bank = polarized(
                cx,
                gate,
                source(80e9, 1000.0).with_torsion("G", "eta_s", 3).unwrap(),
                d,
            );
            let moving = &bank.source().strings()[0].connections[0].body;
            let fixed = ModalAttachment {
                component: moving.component,
                shapes: vec![0.0; moving.shapes.len()],
            };
            let rotation =
                ModalRotationAttachment::from_displacement_footprints(moving, &fixed, 0.1).unwrap();
            bank.with_rotational_body_connections(vec![rotational_link(rotation)], gate)
                .unwrap()
        };
        let mut a = make(directions([1.0, 0.0], [0.6, 0.8]));
        let mut b = make(directions([0.6, 0.8], [-0.28, 0.96]));
        let mut transverse_only = polarized(
            cx,
            gate,
            source(80e9, 1000.0).with_torsion("G", "eta_s", 3).unwrap(),
            directions([1.0, 0.0], [0.6, 0.8]),
        );
        assert_eq!(a.source().rotational_body_connections().len(), 1);
        assert_eq!(
            a.system().acceleration_pressure_weights(),
            transverse_only.system().acceleration_pressure_weights()
        );
        assert_eq!(a.system().components()[0].modes().len(), 11);
        let mut external = vec![0.0; a.system().mode_count()];
        // Drive only the first material torsion coordinate. The second
        // transverse plane occupies 4..8 and must never be mistaken for twist.
        external[8] = 0.002;
        let released = FrictionDrive {
            speed_m_s: 0.0,
            normal_force_n: 0.0,
        };
        let mut peak = 0.0_f64;
        let mut loss = 0.0;
        let mut release_energy = 0.0;
        for sample in 0..192 {
            if sample == 32 {
                release_energy = a.system().total_energy_j().unwrap();
                external[8] = 0.0;
            }
            let fa = *a.step(&external, released).unwrap();
            let fb = *b.step(&external, released).unwrap();
            let dark = *transverse_only.step(&external, released).unwrap();
            assert_eq!(
                dark.observer_pressure_pa, 0.0,
                "torsion acquired a direct radiation path without the joint"
            );
            near(fa.observer_pressure_pa, fb.observer_pressure_pa, 1e-12);
            near(fa.external_work_j, fb.external_work_j, 1e-17);
            assert_eq!(fa.drive_work_j, 0.0);
            assert!(fa.energy_residual_j.abs() <= fa.energy_tolerance_j);
            assert!(fb.energy_residual_j.abs() <= fb.energy_tolerance_j);
            if sample >= 32 {
                loss += fa.network_dissipation_j;
                assert_eq!(fa.external_work_j, 0.0);
                peak = peak.max(fa.observer_pressure_pa.abs());
            }
            let qa = a.system().components()[0].states();
            let qb = b.system().components()[0].states();
            for k in 0..4 {
                let x = qa[k].displacement_m_sqrt_kg;
                let y = qa[4 + k].displacement_m_sqrt_kg;
                near(qb[k].displacement_m_sqrt_kg, 0.6 * x - 0.8 * y, 1e-14);
                near(qb[4 + k].displacement_m_sqrt_kg, 0.8 * x + 0.6 * y, 1e-14);
            }
            for k in 8..11 {
                near(
                    qa[k].displacement_m_sqrt_kg,
                    qb[k].displacement_m_sqrt_kg,
                    1e-14,
                );
            }
        }
        assert!(
            peak > 1e-12,
            "source-derived torque did not reach body pressure"
        );
        assert!(loss > 0.0);
        near(
            a.system().total_energy_j().unwrap() + loss,
            release_energy,
            1e-12,
        );
    });
}

#[test]
fn material_stops_preserve_two_plane_covariance_and_never_add_radial_torque() {
    use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::stops::{
        MaterialStringStop, StringStopSide,
    };
    use fs_couple::render::schedule::force::coupled::contact::ModalContactConfig;
    use fs_couple::render::schedule::force::coupled::contact::friction::multiple::obstacles::FrictionObstacleConfig;
    use fs_dcontact::Obstacle;
    with_cx(|cx, gate| {
        let port = StringForcePort {
            location: StringBowLocation {
                string: 0,
                station_fraction: 0.43,
            },
            coupling: StringForceCoupling::Centreline,
        };
        let make = |directions| {
            polarized(
                cx,
                gate,
                source(80e9, 1000.0).with_torsion("G", "eta_s", 3).unwrap(),
                directions,
            )
            .with_additional_bows(
                vec![],
                MultiFrictionConfig {
                    max_ports: 4,
                    max_sweeps: 128,
                    max_setup_terms: 65_536,
                },
                gate,
            )
            .unwrap()
            .with_string_stops(
                vec![MaterialStringStop {
                    location: port.location,
                    side: StringStopSide::Positive,
                    body: None,
                    law: Obstacle::new(
                        vec![-1.0],
                        1,
                        1,
                        vec![0.0],
                        vec![1.0],
                        5000.0,
                        1.0,
                        "synthetic two-plane contact fixture".into(),
                    )
                    .unwrap()
                    .with_internal_loss(0.2)
                    .unwrap(),
                    limits: ModalContactConfig {
                        max_iterations: 96,
                        maximum_force_n: 10.0,
                        maximum_penetration_m: 0.001,
                        force_absolute_tolerance_n: 1e-11,
                        force_relative_tolerance: 1e-9,
                    },
                }],
                FrictionObstacleConfig {
                    max_contacts: 4,
                    max_setup_terms: 65_536,
                },
                gate,
            )
            .unwrap()
        };
        let mut a = make(directions([0.6, 0.8], [1.0, 0.0]));
        let mut b = make(directions([-0.28, 0.96], [0.6, 0.8]));
        let force = |bank: &MultiBowStringBank| {
            bank.string_force_columns(&[port], 4096, gate).unwrap()[0]
                .iter()
                .map(|x| 0.001 * x)
                .collect::<Vec<_>>()
        };
        let fa = force(&a);
        let fb = force(&b);
        let normal = &a.system().obstacle(0).unwrap().left.shapes;
        assert_eq!(normal.len(), 11);
        assert!(normal[8..].iter().all(|x| *x == 0.0));
        assert!(normal[..4].iter().any(|x| *x != 0.0));
        assert!(normal[4..8].iter().any(|x| *x != 0.0));
        let zero = vec![0.0; fa.len()];
        let released = [FrictionDrive {
            speed_m_s: 0.0,
            normal_force_n: 0.0,
        }];
        let mut engaged = false;
        for sample in 0..96 {
            let x = a
                .step(if sample < 32 { &fa } else { &zero }, &released)
                .unwrap()
                .clone();
            let y = b
                .step(if sample < 32 { &fb } else { &zero }, &released)
                .unwrap()
                .clone();
            near(x.observer_pressure_pa, y.observer_pressure_pa, 1e-12);
            assert!(x.energy_residual_j.abs() <= x.energy_tolerance_j);
            assert!(y.energy_residual_j.abs() <= y.energy_tolerance_j);
            engaged |= a.system().last_obstacle_frames().unwrap()[0].normal_force_n > 0.0;
        }
        assert!(
            engaged,
            "source-derived two-plane normal contact did not engage"
        );
        assert!(
            a.system().components()[0].states()[8..]
                .iter()
                .all(|s| s.displacement_m_sqrt_kg == 0.0 && s.velocity_m_sqrt_kg_per_s == 0.0)
        );
    });
}
