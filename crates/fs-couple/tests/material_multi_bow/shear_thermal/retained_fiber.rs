//! G1/G3: source fibers, retained string films and finite heat share one step.
use super::*;
use super::super::fiber_ribbon::{normal_fiber_bank, support_limits};
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::{
    obstacles::geometry::ObstacleGeometry,
    surfaces::geometry::{FrictionGeometry, FrictionGeometryConfig},
};

fn films(card: &MaterialCard) -> Vec<(usize, ResolvedShearLayer)> {
    (0..3).map(|port| (port, layer(card))).collect()
}
fn cells(curve: &EquilibriumEnthalpyPhaseCurve) -> DepositionSimulation<'_> {
    let original = thermal(curve, false).cells()[&FILM].clone();
    let events = (0..3).map(|index| {
        let mut cell = original.clone();
        cell.id = CellId(10 + index);
        cell.segment_id = index;
        cell.body_id = index;
        TimedEvent { time_s: 0.0, event: Event::Birth(cell) }
    }).collect();
    DepositionSimulation::new(events, &[(MaterialId(0), curve as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 3, max_links: 0, max_events: 3, max_substeps: 4 })
        .unwrap().advance_to(0.0, StepConfig { max_dt_s: 1e-5, energy_tolerance_j: 1e-14 },
            || Ok(())).unwrap().0
}
fn bindings<'a>(card: &'a MaterialCard, curve: &'a EquilibriumEnthalpyPhaseCurve)
    -> Vec<LayerThermalCell<'a>> {
    (0..3).map(|port| LayerThermalCell { port, cell: CellId(10 + port as u64),
        material: card, caloric: curve, temperature_axis: "T" }).collect()
}
fn state(bank: &MultiBowStringBank) -> Vec<Vec<ModalAcousticState>> {
    bank.system().components().iter().map(|m| m.states().to_vec()).collect()
}

#[test]
fn retained_string_films_follow_accelerating_counterface_and_keep_thermal_history() {
    with_cx(|cx, gate| {
        let (card, capacity) = coating(true, false);
        let make = || normal_fiber_bank(cx, gate, 40.0, -1e-6)
            .with_bow_shear_layers_attached_to_strings(films(&card), gate).unwrap();
        let mut fixed = make();
        let mut moving = make().with_fiber_support_motion(0.0, support_limits(), gate).unwrap();
        let mut heat = cells(capacity.curve());
        let mut fixed_heat = cells(capacity.curve());
        let bound = bindings(&card, capacity.curve());
        let force = vec![0.0; moving.system().mode_count()];
        let initial_energy = moving.system().total_energy_j().unwrap();
        let mut balance = initial_energy;
        let original_left = moving.system().port(0).unwrap().left.clone();
        let original_right = moving.system().port(0).unwrap().right.clone();
        let mut plastic_heat = 0.0;
        let mut pressure_difference = 0.0;
        // Stay in the retained-compression interval; this law deliberately
        // refuses later separation rather than discarding stored film strain.
        for sample in 1..=4 {
            let velocity = 0.001 * sample as f64 / 4.0;
            let drives = [drive(velocity, 0.0); 3];
            let frame = moving.step_with_layer_thermal_network(&mut heat, &bound,
                controls(3), &[0.0; 3], &force, &drives, gate, other).unwrap();
            let stationary = fixed.step_with_layer_thermal_network(&mut fixed_heat, &bound,
                controls(3), &[0.0; 3], &force, &drives, gate, other).unwrap();
            let mechanical = &frame.mechanics;
            let support = moving.fiber_support_frame().unwrap();
            let normal = moving.system().last_obstacle_frames().unwrap();
            for port in 0..3 {
                assert_eq!(mechanical.ports[port].normal_force_n, normal[port].normal_force_n);
                assert!(normal[port].normal_force_n > 0.0);
                assert!(moving.system().material_shear_layer_attached_to_left(port));
                assert_eq!(moving.bow_shear_layer(port).unwrap().temperature_kelvin(),
                    heat.temperature_k(CellId(10 + port as u64)).unwrap());
            }
            // Replace support-frame inertial and normal geometry work once;
            // film plastic work is retained in finite heat rather than discarded.
            balance += mechanical.external_work_j - support.inertial_coordinate_work_j
                + support.support_work_j + mechanical.obstacle_work_j
                - support.normal_attachment_work_j - mechanical.network_dissipation_j
                - mechanical.obstacle_dissipation_j - mechanical.obstacle_residual_loss_j;
            assert!((support.laboratory_mechanical_energy_j + total_heat(&heat) - balance).abs() < 2e-9);
            assert!(frame.total_energy_residual_j.abs() <= controls(3).total_energy_tolerance_j);
            plastic_heat += frame.routed_plastic_heat_j;
            pressure_difference += (mechanical.observer_pressure_pa
                - stationary.mechanics.observer_pressure_pa).abs();
        }
        assert_eq!(moving.system().port(0).unwrap().left, original_left);
        assert_ne!(moving.system().port(0).unwrap().right, original_right);
        assert!(moving.fiber_support_frame().unwrap().translation_m > 0.0);
        assert!(plastic_heat > 0.0 && total_heat(&heat) > 0.0 && pressure_difference > 0.0);
        assert!(moving.system().material_shear_state(0).unwrap().plastic_slip_m()[0] != 0.0);
    });
}

#[test]
fn retained_film_rejects_film_side_motion_and_rolls_back_late_thermal_failure() {
    with_cx(|cx, gate| {
        let (card, capacity) = coating(true, false);
        assert!(normal_fiber_bank(cx, gate, 40.0, -1e-6)
            .with_bow_shear_layers(films(&card), gate).unwrap()
            .with_fiber_support_motion(0.0, support_limits(), gate).is_err());
        let make = || normal_fiber_bank(cx, gate, 40.0, -1e-6)
            .with_bow_shear_layers_attached_to_strings(films(&card), gate).unwrap()
            .with_fiber_support_motion(0.0, support_limits(), gate).unwrap();
        let (_, mut generic) = make().into_parts();
        let maps = (0..3).map(|i| { let p = generic.port(i).unwrap();
            FrictionGeometry { left: p.left.clone(), right: p.right.clone() } }).collect::<Vec<_>>();
        let normals = (0..3).map(|i| { let p = generic.obstacle(i).unwrap();
            ObstacleGeometry { left: p.left.clone(), right: p.right.clone() } }).collect::<Vec<_>>();
        let zero = vec![0.0; generic.mode_count()];
        let controls_geometry = FrictionGeometryConfig { max_update_terms: 1_048_576 };
        let mut moved = maps.clone();
        moved[0].left.shapes[0] += 1e-6;
        assert!(generic.step_with_friction_geometry(&zero, &[drive(0.001, 0.0); 3],
            &moved, Some(&normals), None, controls_geometry, gate).is_err());
        let mut moved_normals = normals.clone();
        moved_normals[0].left.shapes[0] += 1e-6;
        assert!(generic.step_with_friction_geometry(&zero, &[drive(0.001, 0.0); 3],
            &maps, Some(&moved_normals), None, controls_geometry, gate).is_err());
        assert_eq!(generic.samples_rendered(), 0);
        let mut actual = make();
        let mut expected = make();
        let mut heat = cells(capacity.curve());
        let mut expected_heat = cells(capacity.curve());
        let bound = bindings(&card, capacity.curve());
        let forces = vec![0.0; actual.system().mode_count()];
        let drives = [drive(0.001, 0.0); 3];
        actual.step_with_layer_thermal_network(&mut heat, &bound, controls(3), &[0.0; 3],
            &forces, &drives, gate, other).unwrap();
        expected.step_with_layer_thermal_network(&mut expected_heat, &bound, controls(3), &[0.0; 3],
            &forces, &drives, gate, other).unwrap();
        let before = state(&actual);
        let support = actual.fiber_support_frame().copied();
        let cells = heat.cells().clone();
        let history = actual.system().material_shear_state(0).unwrap().clone();
        let frame = actual.system().last_frame().unwrap().clone();
        let maps = (0..3).map(|i| actual.system().port(i).unwrap().left.clone()).collect::<Vec<_>>();
        let mut stations = actual.bow_station_fractions().unwrap().to_vec();
        for station in &mut stations { *station += 1e-5; }
        assert!(actual.step_with_bow_stations(&forces, &drives, &stations, None).is_err());
        let failure = |_: &DepositionSimulation<'_>| -> Result<StepCoefficients<()>, &'static str> {
            Err("manufactured late source refusal after moving mechanics")
        };
        assert!(actual.step_with_layer_thermal_network(&mut heat, &bound, controls(3), &[0.0; 3],
            &forces, &drives, gate, failure).is_err());
        assert_eq!(state(&actual), before);
        assert_eq!(actual.fiber_support_frame().copied(), support);
        assert_eq!(heat.cells(), &cells);
        assert_eq!(actual.system().last_frame(), Some(&frame));
        assert_eq!(actual.system().material_shear_state(0), Some(&history));
        for (i, map) in maps.iter().enumerate() { assert_eq!(&actual.system().port(i).unwrap().left, map); }
        let a = actual.step_with_layer_thermal_network(&mut heat, &bound, controls(3), &[0.0; 3],
            &forces, &drives, gate, other).unwrap();
        let b = expected.step_with_layer_thermal_network(&mut expected_heat, &bound, controls(3), &[0.0; 3],
            &forces, &drives, gate, other).unwrap();
        assert_eq!(a.mechanics, b.mechanics);
        assert_eq!(state(&actual), state(&expected));
        assert_eq!(heat.cells(), expected_heat.cells());
        assert_eq!(actual.fiber_support_frame(), expected.fiber_support_frame());
    });
}
