//! Material geometry -> compliant stop -> shared bow/body mechanics and pressure.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::stops::{
    MaterialStringStop, StringStopSide,
};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::force::{
    StringForceCoupling, StringForcePort, StringForceProgram,
};
use fs_couple::render::schedule::force::{ForceInitialization, ForceRenderConfig, ModalForceEvent};
use fs_couple::render::schedule::force::coupled::ModalAttachment;
use fs_couple::render::schedule::force::coupled::contact::ModalContactConfig;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::obstacles::FrictionObstacleConfig;
use fs_dcontact::Obstacle;

fn limits() -> FrictionObstacleConfig { FrictionObstacleConfig { max_contacts: 8, max_setup_terms: 100_000 } }
fn port() -> StringForcePort { StringForcePort {
    location: StringBowLocation { string: 0, station_fraction: 0.43 }, coupling: StringForceCoupling::Centreline,
} }
fn stop() -> MaterialStringStop {
    MaterialStringStop { location: port().location, side: StringStopSide::Positive, body: None,
        law: Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 5000.0, 1.0,
            "synthetic compliant stop, not measured finger tissue".into()).unwrap().with_internal_loss(0.2).unwrap(),
        limits: ModalContactConfig { max_iterations: 96, maximum_force_n: 10.0, maximum_penetration_m: 0.001,
            force_absolute_tolerance_n: 1e-11, force_relative_tolerance: 1e-9 } }
}
fn stopped(cx: &Cx<'_>, gate: &CancelGate) -> MultiBowStringBank {
    multi(cx, gate).with_string_stops(vec![stop()], limits(), gate).unwrap()
}
fn released() -> [FrictionDrive; 2] { [FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 }; 2] }
fn forced(bank: &MultiBowStringBank, gate: &CancelGate) -> Vec<f64> {
    bank.string_force_columns(&[port()], 1000, gate).unwrap()[0].iter().map(|x| 0.001*x).collect()
}

#[test]
fn stops_use_original_material_centreline_maps_without_inventing_radial_torque() {
    with_cx(|cx, gate| {
        let build_density = |rho| {
            let binding = interface(0.4, 10.0);
            let bank = ReciprocalStringBank::from_material_plate(cx, vec![source_at_density(rho)],
                GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap(), body(),
                binding.card().friction_law().clone(), options(), gate).unwrap().with_bow_interface(binding).unwrap();
            bank.with_additional_bows(vec![], MultiFrictionConfig { max_ports: 4, max_sweeps: 128,
                max_setup_terms: 100_000 }, gate).unwrap().with_string_stops(vec![stop()], limits(), gate).unwrap()
        };
        let a = build_density(7800.0); let b = build_density(15600.0);
        let expected = a.string_force_columns(&[port()], 1000, gate).unwrap();
        let normal = &a.system().obstacle(0).unwrap().left.shapes;
        assert_eq!(normal, &expected[0][..7]);
        assert!(normal[4..].iter().all(|x| *x == 0.0));
        assert!(a.system().port(0).unwrap().left.shapes[4..].iter().any(|x| *x != 0.0));
        for (x, y) in normal[..4].iter().zip(&b.system().obstacle(0).unwrap().left.shapes[..4]) {
            assert!((y/x-0.5_f64.sqrt()).abs() < 1e-12);
        }
        assert_eq!(a.system().components()[1].modes(), b.system().components()[1].modes());
        assert_eq!(a.source().stops()[0].location, port().location);
        assert_eq!(a.bow_interface(0), Some(&interface(0.4, 10.0)));
    });
}

#[test]
fn physical_stop_changes_string_motion_and_body_sound_and_release_keeps_ringdown() {
    with_cx(|cx, gate| {
        let mut a = multi(cx, gate); let mut b = stopped(cx, gate);
        let load = forced(&b, gate); let zero = vec![0.0; load.len()];
        let (mut hit, mut difference) = (false, 0.0);
        for _ in 0..128 {
            let fa = a.step(&load, &drivers()).unwrap().clone();
            let fb = b.step(&load, &drivers()).unwrap().clone();
            hit |= b.system().last_obstacle_frames().unwrap()[0].normal_force_n > 0.0;
            difference += (fa.observer_pressure_pa-fb.observer_pressure_pa).abs();
            assert!(fb.energy_residual_j.abs() <= fb.energy_tolerance_j);
        }
        assert!(hit && difference > 1e-15);
        assert_ne!(a.system().components()[0].states(), b.system().components()[0].states());
        let mut energy = b.system().total_energy_j().unwrap(); assert!(energy > 0.0);
        for _ in 0..128 {
            let f = b.step(&zero, &released()).unwrap().clone();
            let next = b.system().total_energy_j().unwrap();
            assert!(next <= energy+f.energy_tolerance_j); energy = next;
            assert_eq!(f.external_work_j, 0.0); assert_eq!(f.drive_work_j, 0.0);
        }
        assert!(energy > 0.0);
    });
}

#[test]
fn opposite_side_and_receiver_use_the_same_signed_displacement_force_pair() {
    with_cx(|cx, gate| {
        let bank = multi(cx, gate);
        let receiver = bank.source().bank().strings()[0].connections[0].body.clone();
        let expected = bank.string_force_columns(&[port()], 1000, gate).unwrap();
        let mut stop = stop(); stop.side = StringStopSide::Negative; stop.body = Some(receiver.clone());
        let bank = bank.with_string_stops(vec![stop], limits(), gate).unwrap();
        let contact = bank.system().obstacle(0).unwrap();
        assert_eq!(contact.left.component, 0); assert_eq!(contact.right.component, 1);
        assert_eq!(contact.left.shapes, expected[0][..7].iter().map(|x| -x).collect::<Vec<_>>());
        assert_eq!(contact.right.shapes, receiver.shapes.iter().map(|x| -x).collect::<Vec<_>>());
        assert_eq!(bank.source().stops()[0].body.as_ref().unwrap(), &receiver);
        assert_eq!(bank.system().components().len(), 2, "no invented receiving mass or duplicate body");
    });
}

fn schedule() -> GestureSchedule {
    GestureSchedule::try_new(1000, ["first", "second"].into_iter().zip(drivers()).zip([0.2, 0.31])
        .map(|((id, d), station)| GestureTrack { id: id.into(), target: GestureTarget::BowStroke { string: 0 },
            initial: GestureValue::Bow { velocity_m_per_s: d.speed_m_s, normal_force_n: d.normal_force_n, station },
            events: vec![] }).collect()).unwrap()
}
fn program(initialization: ForceInitialization) -> StringForceProgram {
    StringForceProgram { ports: vec![port()], initial_force_n: vec![0.001],
        events: vec![ModalForceEvent { sample: 64, voice: 0, port: 0, force_n: 0.0 }], initialization,
        limits: ForceRenderConfig { sample_rate_hz: 48000, max_block: 64, max_events: 16,
            max_controls: 1024, max_projection_terms: 10000 } }
}

#[test]
fn stopped_material_bows_and_physical_force_history_share_direct_and_scheduled_pressure() {
    with_cx(|cx, gate| {
        let mut direct = stopped(cx, gate); let load = forced(&direct, gate); let zero = vec![0.0; load.len()];
        let expected: Vec<_> = (0..96).map(|i| direct.step(if i < 64 { &load } else { &zero }, &drivers())
            .unwrap().observer_pressure_pa.to_bits()).collect();
        for partition in [1, 37, 64] {
            let mut performance = stopped(cx, gate).into_scheduled(schedule(), vec!["first".into(), "second".into()],
                MultiBowPerformanceConfig { samples: 96, max_block: 64, max_work: 10000, max_events: 64 }).unwrap()
                .with_string_forces(program(ForceInitialization::RetainState), gate).unwrap();
            let mut pressure = vec![0.0; 96];
            for block in pressure.chunks_mut(partition) { performance.block(block).unwrap(); }
            assert_eq!(pressure.iter().map(|x| x.to_bits()).collect::<Vec<_>>(), expected);
            assert_eq!(performance.bank().source().stops().len(), 1);
            for i in 0..2 { assert_eq!(performance.bank().bow_interface(i), Some(&interface(0.4, 10.0))); }
            assert_eq!(performance.samples_rendered(), 96);
        }
    });
}

#[test]
fn unsupported_stop_geometry_and_linear_equilibrium_never_become_silent_fallbacks() {
    with_cx(|cx, gate| {
        for which in 0..4 {
            let mut stop = stop();
            match which {
                0 => stop.location.station_fraction = 1.0,
                1 => stop.location.string = 7,
                2 => stop.body = Some(ModalAttachment { component: 99, shapes: vec![1.0, 0.0] }),
                _ => stop.body = Some(ModalAttachment { component: 0, shapes: vec![0.0, 0.0] }),
            }
            assert!(multi(cx, gate).with_string_stops(vec![stop], limits(), gate).is_err());
        }
        assert!(stopped(cx, gate).with_string_stops(vec![stop()], limits(), gate).is_err());
        let performance = stopped(cx, gate).into_scheduled(schedule(), vec!["first".into(), "second".into()],
            MultiBowPerformanceConfig { samples: 96, max_block: 64, max_work: 10000, max_events: 64 }).unwrap();
        assert!(performance.with_string_forces(program(ForceInitialization::StaticPreload), gate).is_err());
    });
}
