//! G0/G3 geometry interventions in the existing source-bound mechanical performer.
//! All material/contact coefficients here are synthetic, not instrument calibration.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::{
    TransverseDirection as Direction, TransverseDirections,
};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::stops::{
    MaterialStringStop, StringStopSide,
};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::force::{
    StringForceCoupling, StringForcePort, StringForceProgram,
};
use fs_couple::render::schedule::force::{ForceInitialization, ForceRenderConfig, ModalForceEvent};
use fs_couple::render::schedule::force::coupled::contact::ModalContactConfig;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::obstacles::FrictionObstacleConfig;
use fs_dcontact::Obstacle;

fn port(direction: Direction) -> StringForcePort {
    StringForcePort { location: StringBowLocation { string: 0, station_fraction: 0.43 },
        coupling: StringForceCoupling::CentrelineAlong(direction) }
}
fn stop(direction: Direction) -> MaterialStringStop {
    MaterialStringStop { location: port(direction).location, side: StringStopSide::Along(direction), body: None,
        law: Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 5000.0, 1.0,
            "synthetic normal, independent of bow tangent".into()).unwrap().with_internal_loss(0.2).unwrap(),
        limits: ModalContactConfig { max_iterations: 96, maximum_force_n: 10.0, maximum_penetration_m: 0.001,
            force_absolute_tolerance_n: 1e-11, force_relative_tolerance: 1e-9 } }
}
fn limits() -> FrictionObstacleConfig {
    FrictionObstacleConfig { max_contacts: 4, max_setup_terms: 100_000 }
}
fn bank(cx: &Cx<'_>, gate: &CancelGate, rho: f64, bow: [f64; 2], bridge: [f64; 2]) -> MultiBowStringBank {
    let source = source_at_density(rho).with_transverse_polarizations(
        TransverseDirections::try_new(bow, bridge).unwrap()).unwrap();
    let binding = interface(0.4, 10.0);
    ReciprocalStringBank::from_material_plate(cx, vec![source],
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap(),
        body(), binding.card().friction_law().clone(), options(), gate).unwrap()
        .with_bow_interface(binding).unwrap()
        .with_additional_bows(vec![], MultiFrictionConfig { max_ports: 4, max_sweeps: 128,
            max_setup_terms: 100_000 }, gate).unwrap()
}
fn unrotated(cx: &Cx<'_>, gate: &CancelGate) -> MultiBowStringBank {
    bank(cx, gate, 7800.0, [1.0, 0.0], [0.6, 0.8])
}
fn load(bank: &MultiBowStringBank, direction: Direction, gate: &CancelGate) -> Vec<f64> {
    bank.string_force_columns(&[port(direction)], 1000, gate).unwrap()[0]
        .iter().map(|x| 0.001*x).collect()
}
fn release() -> [FrictionDrive; 1] { [FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 }] }

#[test]
fn unit_directions_and_missing_transverse_planes_refuse_instead_of_normalizing() {
    for direction in [[0.0, 0.0], [2.0, 0.0], [f64::INFINITY, 0.0], [f64::NAN, 1.0], [f64::MAX, 1.0]] {
        assert!(Direction::try_new(direction).is_err());
    }
    assert_eq!(Direction::try_new([0.6, 0.8]).unwrap().components(), [0.6, 0.8]);
    with_cx(|cx, gate| {
        let one_plane = multi(cx, gate);
        assert!(one_plane.string_force_columns(&[port(Direction::Y)], 1000, gate).is_err());
        let negative_x = Direction::try_new([-1.0, 0.0]).unwrap();
        let rows = one_plane.string_force_columns(&[port(Direction::X), port(negative_x)], 1000, gate).unwrap();
        for (x, y) in rows[0].iter().zip(&rows[1]) { assert_eq!(*x, -*y); }
        assert!(one_plane.with_string_stops(vec![stop(Direction::Y)], limits(), gate).is_err());
    });
}

#[test]
fn independent_normal_and_load_maps_keep_material_mass_and_zero_radial_torque() {
    with_cx(|cx, gate| {
        let a = unrotated(cx, gate);
        let b = bank(cx, gate, 15_600.0, [1.0, 0.0], [0.6, 0.8]);
        let mut legacy = port(Direction::X); legacy.coupling = StringForceCoupling::Centreline;
        let rows = a.string_force_columns(&[legacy, port(Direction::X), port(Direction::Y)], 1000, gate).unwrap();
        assert_eq!(rows[0], rows[1]);
        assert!(rows[2][..4].iter().all(|x| *x == 0.0));
        assert_eq!(&rows[1][..4], &rows[2][4..8]);
        assert!(rows[2][8..].iter().all(|x| *x == 0.0));
        let heavy = b.string_force_columns(&[port(Direction::Y)], 1000, gate).unwrap();
        for (light, heavy) in rows[2][4..8].iter().zip(&heavy[0][4..8]) {
            assert!((heavy/light - 0.5_f64.sqrt()).abs() < 1e-12);
        }
        let a = a.with_string_stops(vec![stop(Direction::Y)], limits(), gate).unwrap();
        let normal = &a.system().obstacle(0).unwrap().left.shapes;
        assert_eq!(normal, &rows[2][..11]);
        let bow = &a.system().port(0).unwrap().left.shapes;
        assert!(bow[4..8].iter().all(|x| *x == 0.0));
        assert!(bow[8..].iter().any(|x| *x != 0.0));
        assert_eq!(a.bow_interface(0), Some(&interface(0.4, 10.0)));
    });
}

#[test]
fn receiver_normal_is_retained_without_a_hidden_sign_or_body_duplication() {
    with_cx(|cx, gate| {
        let a = unrotated(cx, gate);
        let direction = Direction::try_new([0.0, -1.0]).unwrap();
        let mut s = stop(direction);
        let mut receiver = a.source().bank().strings()[0].connections[0].body.clone();
        for x in &mut receiver.shapes { *x = -*x; }
        s.body = Some(receiver.clone());
        let expected = a.string_force_columns(&[port(direction)], 1000, gate).unwrap();
        let a = a.with_string_stops(vec![s], limits(), gate).unwrap();
        let pair = a.system().obstacle(0).unwrap();
        assert_eq!(pair.left.shapes, expected[0][..11]);
        assert_eq!(pair.right.shapes, receiver.shapes);
        assert_eq!(pair.right.component, 1);
        assert_eq!(a.system().components().len(), 2);
    });
}

#[test]
fn rotating_only_the_stop_normal_changes_mechanics_and_physical_body_pressure() {
    with_cx(|cx, gate| {
        let mut a = unrotated(cx, gate).with_string_stops(vec![stop(Direction::X)], limits(), gate).unwrap();
        let mut b = unrotated(cx, gate).with_string_stops(vec![stop(Direction::Y)], limits(), gate).unwrap();
        let force = load(&b, Direction::Y, gate);
        let mut difference = 0.0; let mut hit = false;
        for _ in 0..128 {
            let fa = a.step(&force, &release()).unwrap().clone();
            let fb = b.step(&force, &release()).unwrap().clone();
            difference += (fa.observer_pressure_pa-fb.observer_pressure_pa).abs();
            hit |= b.system().last_obstacle_frames().unwrap()[0].normal_force_n > 0.0;
            assert!(fb.energy_residual_j.abs() <= fb.energy_tolerance_j);
        }
        assert!(hit && difference > 1e-15);
        assert_ne!(a.system().components()[0].states(), b.system().components()[0].states());
        // Centreline loading/normal contact never becomes a torsional drive.
        assert!(b.system().components()[0].states()[8..].iter().all(|s|
            s.displacement_m_sqrt_kg == 0.0 && s.velocity_m_sqrt_kg_per_s == 0.0));
    });
}

#[test]
fn rotating_every_physical_direction_preserves_the_instrument_response() {
    with_cx(|cx, gate| {
        let left = Direction::try_new([-1.0, 0.0]).unwrap();
        let mut a = unrotated(cx, gate).with_string_stops(vec![stop(Direction::Y)], limits(), gate).unwrap();
        let mut b = bank(cx, gate, 7800.0, [0.0, 1.0], [-0.8, 0.6])
            .with_string_stops(vec![stop(left)], limits(), gate).unwrap();
        let fa = load(&a, Direction::Y, gate); let fb = load(&b, left, gate);
        let drive = [drivers()[0]];
        let mut peak = 0.0_f64; let mut error = 0.0_f64;
        for _ in 0..96 {
            let a = a.step(&fa, &drive).unwrap().clone();
            let b = b.step(&fb, &drive).unwrap().clone();
            peak = peak.max(a.observer_pressure_pa.abs()).max(b.observer_pressure_pa.abs());
            error = error.max((a.observer_pressure_pa-b.observer_pressure_pa).abs());
            assert!((a.network_energy_j-b.network_energy_j).abs() <= a.energy_tolerance_j+b.energy_tolerance_j);
        }
        assert!(peak > 1e-12);
        assert!(error <= 1e-7*peak + 1e-12, "basis rotation changed pressure: {error}/{peak}");
    });
}

#[test]
fn directed_stopped_preload_release_and_callback_partitions_use_the_same_solver() {
    with_cx(|cx, gate| {
        let construct = || unrotated(cx, gate).with_string_stops(vec![stop(Direction::Y)], limits(), gate).unwrap();
        let direct = construct(); let force = load(&direct, Direction::Y, gate);
        let zero = vec![0.0; force.len()];
        let (_, mut direct) = direct.into_parts();
        let energy = direct.initialize_released_equilibrium(&force, gate).unwrap(); assert!(energy > 0.0);
        let expected: Vec<_> = (0..96).map(|i| direct.step(if i < 32 { &force } else { &zero }, &release())
            .unwrap().observer_pressure_pa.to_bits()).collect();
        let schedule = GestureSchedule::try_new(700, vec![GestureTrack { id: "bow".into(),
            target: GestureTarget::BowStroke { string: 0 }, initial: GestureValue::Bow {
                velocity_m_per_s: 0.0, normal_force_n: 0.0, station: 0.2 }, events: vec![] }]).unwrap();
        for partition in [1, 37, 64] {
            let program = StringForceProgram { ports: vec![port(Direction::Y)], initial_force_n: vec![0.001],
                events: vec![ModalForceEvent { sample: 32, voice: 0, port: 0, force_n: 0.0 }],
                initialization: ForceInitialization::StaticPreload,
                limits: ForceRenderConfig { sample_rate_hz: 48000, max_block: 64, max_events: 16,
                    max_controls: 1024, max_projection_terms: 10000 } };
            let mut p = construct().into_scheduled(schedule.clone(), vec!["bow".into()], MultiBowPerformanceConfig {
                samples: 96, max_block: 64, max_work: 10000, max_events: 64 }).unwrap()
                .with_string_forces(program, gate).unwrap();
            let mut output = vec![0.0; 96]; for block in output.chunks_mut(partition) { p.block(block).unwrap(); }
            assert_eq!(output.iter().map(|p| p.to_bits()).collect::<Vec<_>>(), expected);
            assert_eq!(p.bank().source().stops()[0].side, StringStopSide::Along(Direction::Y));
        }
    });
}
