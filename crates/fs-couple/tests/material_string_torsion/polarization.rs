//! G1/G3 physical two-plane response on the SAME material/body/friction path.
//! Sources are the parent's explicit synthetic material fixtures, not measured
//! violin specimens. Rotation and work checks test wiring, not timbre realism.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::force::{
    StringForceCoupling, StringForcePort,
};
use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::{
    StringBodyTranslationConnection, TransverseDirection, TransverseDirections,
};

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

// Two physical receiver displacement coordinates in one two-mode body basis.
// The mode shapes are an explicit synthetic fixture, not a measured violin.
fn two_axis_bank(cx: &Cx<'_>, gate: &CancelGate, rho: f64, axes: TransverseDirections)
    -> ReciprocalStringBank
{
    let input = MaterialStringSource::Homogeneous(with_uniform_circular_material_state(
        template(), 0.0003, &material(rho, 80e9, 1000.0, 200_000.0),
    ).unwrap().with_kelvin_voigt_bending_loss().unwrap())
        .with_torsion("G", "eta_s", 3).unwrap()
        .with_transverse_polarizations(axes).unwrap();
    let receiver = ModalAcousticTimeModel::try_new(48_000,
        [450.0, 700.0].into_iter().map(|angular_frequency_rad_s| ModalAcousticMode {
            angular_frequency_rad_s, damping_ratio: 0.01,
            pressure_per_modal_velocity: fs_math::c64::C64::new(1.0, 0.0),
        }).collect(), ModalAcousticTimeBudget::audible_reference()).unwrap();
    ReciprocalStringBank::new(cx, vec![MaterialStringMember { source: input,
        connections: vec![InteriorStringConnection { station_fraction: 0.8,
            body: ModalAttachment { component: 0, shapes: vec![1.0, 0.25] },
            stiffness_n_m: 3000.0, damping_n_s_m: 0.1, rest_extension_m: 0.0 }],
    }], gas(), vec![receiver], law(), options(), gate).unwrap()
}

fn second_axis(direction: TransverseDirection) -> StringBodyTranslationConnection {
    StringBodyTranslationConnection { string: 0, station_fraction: 0.8, direction,
        body: ModalAttachment { component: 0, shapes: vec![0.35, 0.8] },
        stiffness_n_m: 1200.0, damping_n_s_m: 0.04, rest_extension_m: 0.0 }
}

#[test]
fn g3_second_bridge_axis_radiates_reciprocally_and_a_material_swap_changes_the_sound() {
    with_cx(|cx, gate| {
        let axes = directions([1.0, 0.0], [1.0, 0.0]);
        let mut dark = two_axis_bank(cx, gate, 7800.0, axes);
        let mut coupled = two_axis_bank(cx, gate, 7800.0, axes)
            .with_translational_body_connections(vec![second_axis(TransverseDirection::Y)], gate).unwrap();
        let mut heavier = two_axis_bank(cx, gate, 15600.0, axes)
            .with_translational_body_connections(vec![second_axis(TransverseDirection::Y)], gate).unwrap();
        assert_eq!(coupled.source().translational_body_connections(), heavier.source().translational_body_connections());
        assert_eq!(coupled.cards()[0].tension_n, heavier.cards()[0].tension_n);
        assert_eq!(coupled.cards()[0].bending_stiffness_n_m2, heavier.cards()[0].bending_stiffness_n_m2);
        near(heavier.cards()[0].fundamental_hz()/coupled.cards()[0].fundamental_hz(),
            1.0/2.0_f64.sqrt(), 1e-14);
        let released = FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 };
        let (mut peak, mut change) = (0.0_f64, 0.0_f64);
        let (mut work, mut loss) = (0.0, 0.0);
        for n in 0..384 {
            // The same physical transverse point force, entirely in plane Y.
            let force = |bank: &ReciprocalStringBank| {
                let mut force = vec![0.0; bank.system().mode_count()];
                if n < 32 {
                    let card = &bank.cards()[0];
                    let normalization = (0.5*card.length_m*card.linear_density_kg_m).sqrt();
                    for k in 0..4 {
                        force[4+k] = 0.01*((k+1) as f64*core::f64::consts::PI*0.3).sin()/normalization;
                    }
                }
                force
            };
            let fd = *dark.step(&force(&dark), released).unwrap();
            let a = *coupled.step(&force(&coupled), released).unwrap();
            let b = *heavier.step(&force(&heavier), released).unwrap();
            assert_eq!(fd.observer_pressure_pa, 0.0, "the absent bridge axis must stay disconnected");
            peak = peak.max(a.observer_pressure_pa.abs());
            change = change.max((a.observer_pressure_pa-b.observer_pressure_pa).abs());
            work += a.external_work_j;
            loss += a.network_dissipation_j;
            for frame in [a, b] {
                assert_eq!(frame.drive_work_j, 0.0);
                assert_eq!(frame.friction_dissipation_j, 0.0);
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            }
            assert!(coupled.system().components()[0].states()[8..].iter().all(|s|
                s.displacement_m_sqrt_kg == 0.0 && s.velocity_m_sqrt_kg_per_s == 0.0),
                "centreline bridge force must not invent torque");
        }
        assert!(peak > 1e-8 && change > 1e-8, "new axis or material did not reach physical sound");
        assert!(coupled.system().components()[0].states()[..4].iter()
            .any(|s| s.velocity_m_sqrt_kg_per_s.abs() > 1e-8), "body did not return force into plane X");
        assert!((coupled.system().total_energy_j().unwrap()+loss-work).abs() < 2e-10);
    });
}

#[test]
fn g1_independent_bridge_axes_are_invariant_under_a_common_transverse_basis_rotation() {
    with_cx(|cx, gate| {
        let mut a = two_axis_bank(cx, gate, 7800.0, directions([1.0, 0.0], [1.0, 0.0]))
            .with_translational_body_connections(vec![second_axis(TransverseDirection::Y)], gate).unwrap();
        let mut b = two_axis_bank(cx, gate, 7800.0, directions([0.6, 0.8], [0.6, 0.8]))
            .with_translational_body_connections(vec![second_axis(
                TransverseDirection::try_new([-0.8, 0.6]).unwrap())], gate).unwrap();
        let zero = vec![0.0; a.system().mode_count()];
        for _ in 0..192 {
            let fa = *a.step(&zero, drive()).unwrap();
            let fb = *b.step(&zero, drive()).unwrap();
            near(fa.observer_pressure_pa, fb.observer_pressure_pa, 1e-12);
            near(fa.drive_work_j, fb.drive_work_j, 1e-17);
            near(fa.network_dissipation_j, fb.network_dissipation_j, 1e-17);
            let qa = a.system().components()[0].states();
            let qb = b.system().components()[0].states();
            for k in 0..4 {
                let (x,y) = (qa[k].displacement_m_sqrt_kg, qa[4+k].displacement_m_sqrt_kg);
                near(qb[k].displacement_m_sqrt_kg, 0.6*x-0.8*y, 1e-14);
                near(qb[4+k].displacement_m_sqrt_kg, 0.8*x+0.6*y, 1e-14);
            }
        }
    });
}

#[test]
fn g0_directed_bridge_refuses_missing_planes_invalid_laws_and_live_rebinding() {
    with_cx(|cx, gate| {
        let make = || two_axis_bank(cx, gate, 7800.0, directions([1.0, 0.0], [1.0, 0.0]));
        let good = second_axis(TransverseDirection::Y);
        let mut one_plane_link = good.clone();
        one_plane_link.body.shapes = vec![1.0];
        let missing_plane = raw(cx, gate, source(80e9,1000.0), options()).unwrap()
            .with_translational_body_connections(vec![one_plane_link], gate).err().unwrap();
        assert!(missing_plane.to_string().contains("second-plane"), "{missing_plane}");
        for case in 0..7 {
            let mut bad = good.clone();
            match case {
                0 => bad.string = 1,
                1 => bad.station_fraction = 1.0,
                2 => bad.body.component = 1,
                3 => bad.body.shapes = vec![1.0],
                4 => bad.stiffness_n_m = -1.0,
                5 => bad.damping_n_s_m = f64::NAN,
                _ => bad.rest_extension_m = f64::INFINITY,
            }
            assert!(make().with_translational_body_connections(vec![bad], gate).is_err(), "case {case}");
        }
        assert!(make().with_translational_body_connections(vec![good.clone()], gate).unwrap()
            .with_translational_body_connections(vec![good.clone()], gate).is_err());
        // Residual flexibility admits an axis orthogonal to the string's bridge, in either order,
        // and refuses one with an omitted-mode cross term to it.
        assert!(make().with_residual_flexibility(gate).unwrap()
            .with_translational_body_connections(vec![good.clone()], gate).is_ok());
        assert!(make().with_translational_body_connections(vec![good.clone()], gate).unwrap()
            .with_residual_flexibility(gate).is_ok());
        let oblique = second_axis(TransverseDirection::try_new([0.6, 0.8]).unwrap());
        let refused = make().with_residual_flexibility(gate).unwrap()
            .with_translational_body_connections(vec![oblique.clone()], gate).err().unwrap();
        assert!(refused.to_string().contains("orthogonal"), "{refused}");
        assert!(make().with_translational_body_connections(vec![oblique], gate).unwrap()
            .with_residual_flexibility(gate).is_err());
        let mut moving = make();
        moving.step(&vec![0.0; moving.system().mode_count()], drive()).unwrap();
        assert!(moving.with_translational_body_connections(vec![good], gate).is_err());
    });
}

// The second plane's directed bridge is the mirror image of the first plane's scalar bridge: same
// station, law and body map. With residual flexibility both see the same omitted-mode compliance, so
// a force pattern in plane Y must sound as the same pattern in plane X does. A directed spring left
// uncondensed is stiffer by 1 + k G and breaks the mirror. The binding order must not matter.
#[test]
fn g4_residual_flexibility_condenses_an_orthogonal_directed_bridge_like_its_scalar_mirror() {
    use fs_couple::bowed_string::runtime::reciprocal::bank::residual::residual_compliance_m_n;
    with_cx(|cx, gate| {
        let axes = directions([1.0, 0.0], [1.0, 0.0]);
        let mirror = StringBodyTranslationConnection { string: 0, station_fraction: 0.8,
            direction: TransverseDirection::Y, body: ModalAttachment { component: 0, shapes: vec![1.0, 0.25] },
            stiffness_n_m: 3000.0, damping_n_s_m: 0.1, rest_extension_m: 0.0 };
        let residual_first = || two_axis_bank(cx, gate, 7800.0, axes).with_residual_flexibility(gate).unwrap()
            .with_translational_body_connections(vec![mirror.clone()], gate).unwrap();
        let (mut x, mut y) = (residual_first(), residual_first());
        let mut late = two_axis_bank(cx, gate, 7800.0, axes)
            .with_translational_body_connections(vec![mirror.clone()], gate).unwrap()
            .with_residual_flexibility(gate).unwrap();
        let kg = 3000.0*residual_compliance_m_n(&x.cards()[0], 0.8, 0.8).unwrap();
        assert!(kg > 1e-5, "k G = {kg:e} is too small for the mirror to tell condensed from not");
        let released = FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 };
        let (mut peak, mut work, mut loss) = (0.0_f64, 0.0, 0.0);
        for n in 0..384 {
            let pattern = |bank: &ReciprocalStringBank, plane: usize| {
                let mut force = vec![0.0; bank.system().mode_count()];
                if n < 32 {
                    let card = &bank.cards()[0];
                    let normalization = (0.5*card.length_m*card.linear_density_kg_m).sqrt();
                    for k in 0..4 {
                        force[4*plane+k] = 0.01*((k+1) as f64*core::f64::consts::PI*0.3).sin()/normalization;
                    }
                }
                force
            };
            let fx = *x.step(&pattern(&x, 0), released).unwrap();
            let fy = *y.step(&pattern(&y, 1), released).unwrap();
            let fl = *late.step(&pattern(&late, 1), released).unwrap();
            near(fx.observer_pressure_pa, fy.observer_pressure_pa, 1e-15);
            assert_eq!(fy.observer_pressure_pa.to_bits(), fl.observer_pressure_pa.to_bits(), "binding order changed the sound");
            for frame in [fx, fy] {
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            }
            peak = peak.max(fy.observer_pressure_pa.abs());
            work += fy.external_work_j;
            loss += fy.network_dissipation_j;
        }
        assert!(peak > 1e-8, "plane Y did not reach the body");
        assert!((y.system().total_energy_j().unwrap()+loss-work).abs() < 2e-10);
    });
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
                .is_ok()
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
