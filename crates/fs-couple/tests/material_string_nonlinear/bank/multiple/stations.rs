//! Source mass/EA, moving finger work, fixed normal-loaded bow and gas history.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::stops::stations::StopStationMotionLimits;

fn limits() -> StopStationMotionLimits {
    StopStationMotionLimits { maximum_speed_m_s: 2.0, maximum_step_m: 1e-4, max_update_terms: 1_048_576 }
}
fn bank(cx: &Cx<'_>, gate: &CancelGate, young: f64, medium: GasState, slope: f64) -> MultiBowStringBank {
    let source = physical_bank(cx, gate, young, medium, 1)
        .with_geometric_stretching(slope, nonlinear_config(), gate).unwrap();
    let body = source.source().strings()[0].connections[0].body.clone();
    let modes: Vec<_> = source.system().components().iter().map(|m| m.modes().to_vec()).collect();
    let weights = source.system().acceleration_pressure_weights().unwrap().to_vec();
    let stops = [(0.2, [0.8, -0.6], None), (0.43, [1.0, 0.0], Some(body))].into_iter()
        .map(|(station_fraction, normal, body)| MaterialStringStop {
            location: StringBowLocation { string: 0, station_fraction },
            side: StringStopSide::Along(TransverseDirection::try_new(normal).unwrap()), body,
            law: Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 5000.0, 1.0,
                "synthetic moving finger and fixed bow surface; not measured tissue".into()).unwrap()
                .with_internal_loss(0.2).unwrap(),
            limits: ModalContactConfig { max_iterations: 96, maximum_force_n: 10.0,
                maximum_penetration_m: 0.001, force_absolute_tolerance_n: 1e-11,
                force_relative_tolerance: 1e-10 },
        }).collect();
    let bank = source.with_additional_bows(vec![], joint(), gate).unwrap()
        .with_string_stops(stops, FrictionObstacleConfig { max_contacts: 2, max_setup_terms: 1_048_576 }, gate).unwrap()
        .with_bow_contact_loads(vec![BowContactLoad { bow: 0, stop: 0 }], gate).unwrap()
        .with_stop_motion(ObstacleMotionLimits { maximum_abs_offset_m: 0.001, maximum_speed_m_s: 0.1 }, gate).unwrap()
        .with_stop_station_motion(limits(), gate).unwrap();
    assert_eq!(bank.system().components().iter().map(|m| m.modes().to_vec()).collect::<Vec<_>>(), modes);
    assert_eq!(bank.system().acceleration_pressure_weights(), Some(weights.as_slice()));
    bank
}
fn snapshot(bank: &MultiBowStringBank) -> Vec<Vec<fs_couple::modal_acoustic_time::ModalAcousticState>> {
    bank.system().components().iter().map(|m| m.states().to_vec()).collect()
}

#[test]
fn moving_finger_keeps_shared_strain_actual_bow_load_and_full_source_history() {
    with_cx(|cx, gate| {
        let media = [gas(), GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 60000.0).unwrap()];
        let young = [200e9, 50e9];
        let mut banks = [bank(cx, gate, young[0], media[0], 0.2), bank(cx, gate, young[1], media[1], 0.2)];
        let mut work = [0.0; 2]; let mut loss = [0.0; 2];
        let zero = vec![0.0; banks[0].system().mode_count()];
        let drive = [FrictionDrive { speed_m_s: 0.04, normal_force_n: 0.0 }];
        let mut old_station = 0.43;
        let mut old_offsets = [0.0, 0.0];
        let (mut attachment_work, mut stretch, mut difference) = (0.0, 0.0_f64, 0.0);
        for sample in 0..320 {
            let progress = if sample < 160 { (sample+1) as f64/160.0 } else { (319-sample) as f64/160.0 };
            let station = 0.43+0.004*progress;
            let press = if sample < 64 { (sample+1) as f64/64.0 }
                else if sample < 256 { 1.0 } else { 1.0-2.0*(sample-255) as f64/64.0 };
            let offsets = [-3e-5*press, -1e-5*press];
            let mut pressure = [0.0; 2];
            for (index, bank) in banks.iter_mut().enumerate() {
                let before = snapshot(bank);
                let frame = bank.step_with_stop_stations_under_gate(&zero, &drive, &[0.2, station], Some(&offsets), gate)
                    .unwrap().clone();
                pressure[index] = frame.observer_pressure_pa;
                let contacts = bank.system().last_obstacle_frames().unwrap();
                assert_eq!(frame.ports[0].normal_force_n, contacts[0].normal_force_n);
                assert_eq!(contacts[0].attachment_work_j, 0.0);
                assert_eq!(bank.stop_station_fractions().unwrap(), &[0.2, station]);
                assert_eq!(bank.source().stops()[1].location.station_fraction, 0.43);
                let card = &bank.cards()[0]; let q = bank.system().components()[0].states();
                let mass = card.linear_density_kg_m*card.length_m*0.5;
                // Independently integrate the loaded string's two sine fields.
                let strain = 0.25*(core::f64::consts::PI/card.length_m).powi(2)
                    *(q[0].displacement_m_sqrt_kg.powi(2)+q[1].displacement_m_sqrt_kg.powi(2))/mass;
                let ea = young[index]*core::f64::consts::PI*0.0003_f64.powi(2);
                let geometric = 0.5*ea*card.length_m*strain*strain;
                assert!((frame.nonlinear.unwrap().nonlinear_energy_j-geometric).abs() < 1e-15);
                assert!((bank.stretching_observation(0).unwrap().tension_n-card.tension_n-ea*strain).abs() < 1e-11);
                stretch = stretch.max(geometric);
                let phi = (core::f64::consts::PI*station).sin()/mass.sqrt();
                let old_phi = (core::f64::consts::PI*old_station).sin()/mass.sqrt();
                let actual = bank.system().obstacle_normal_column(1).unwrap();
                assert_eq!(actual[0], phi); assert_eq!(actual[1], 0.0); assert_eq!(actual[2], 0.0);
                let receiver = &bank.system().obstacle(1).unwrap().right.shapes;
                for (i, &b) in receiver.iter().enumerate() { assert_eq!(actual[3+i], -b); }
                assert!(actual[3+receiver.len()..].iter().all(|x| *x == 0.0));
                let migration = contacts[1].normal_force_n*(phi-old_phi)
                    *0.5*(q[0].displacement_m_sqrt_kg+before[0][0].displacement_m_sqrt_kg);
                assert!((contacts[1].attachment_work_j-migration).abs() < 1e-15);
                attachment_work += migration.abs();
                let normal_work = -(0..2).map(|i| contacts[i].normal_force_n*(offsets[i]-old_offsets[i])).sum::<f64>();
                assert!((frame.obstacle_work_j-migration-normal_work).abs() < 1e-15);
                // Contact potential includes the original receiver motion; the
                // normal law has not silently become a fixed-ground penalty.
                let closure = phi*q[0].displacement_m_sqrt_kg-offsets[1]
                    - receiver.iter().zip(bank.system().components()[1].states())
                        .map(|(b,s)| b*s.displacement_m_sqrt_kg).sum::<f64>();
                assert!((contacts[1].stored_energy_j-2500.0*closure.max(0.0).powi(2)).abs() < 1e-14);
                work[index] += frame.external_work_j+frame.drive_work_j+normal_work+migration;
                loss[index] += frame.network_dissipation_j+frame.friction_dissipation_j+frame.obstacle_dissipation_j;
                assert!((bank.system().total_energy_j().unwrap()+loss[index]-work[index]).abs() < 1e-9);
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            }
            difference += (pressure[0]-pressure[1]).abs();
            old_station = station; old_offsets = offsets;
        }
        assert!(attachment_work > 0.0 && stretch > 0.0 && difference > 1e-12);
        for (i, bank) in banks.iter().enumerate() {
            assert_eq!(bank.source().bank().gas(), media[i]);
            assert_eq!(bank.source().bank().cylinder_gas_loading().unwrap().strings().len(), 1);
            assert_ne!(bank.system().components()[0].states()[2].velocity_m_sqrt_kg_per_s, 0.0);
            assert!(bank.system().components()[1].states().iter().any(|s| s.velocity_m_sqrt_kg_per_s != 0.0));
            assert!(bank.system().components()[2..].iter().flat_map(|m| m.states())
                .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0), "migration discarded causal gas history");
        }
    });
}

#[test]
fn station_limits_mapped_normals_and_material_domain_refuse_without_advancing_geometry() {
    with_cx(|cx, gate| {
        let mut actual = bank(cx, gate, 200e9, gas(), 1e-15);
        let mut reference = bank(cx, gate, 200e9, gas(), 1e-15);
        let zero = vec![0.0; actual.system().mode_count()];
        let drives = [FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 }];
        actual.step(&zero, &drives).unwrap(); reference.step(&zero, &drives).unwrap();
        let before = snapshot(&actual); let frame = actual.system().last_frame().unwrap().clone();
        let normal = actual.system().last_obstacle_frames().unwrap().to_vec();
        let columns: Vec<_> = (0..2).map(|i| actual.system().obstacle_normal_column(i).unwrap().to_vec()).collect();
        // This move passes geometry admission and stages real normal reactions,
        // then leaves the retained material stretching domain in the DG trial.
        let error = actual.step_with_stop_stations_under_gate(&zero, &drives, &[0.2, 0.43001], Some(&[-1e-6, -1e-6]), gate).unwrap_err();
        assert!(error.to_string().contains("continuous-span slope envelope"), "{error}");
        for bad in [vec![], vec![0.2], vec![0.2, f64::NAN], vec![0.2, 1.0], vec![0.2, 0.9], vec![0.20001, 0.43]] {
            assert!(actual.step_with_stop_stations(&zero, &drives, &bad, None).is_err());
        }
        let cancel = CancelGate::new_clock_free(); cancel.request();
        assert!(actual.step_with_stop_stations_under_gate(&zero, &drives, &[0.2, 0.43001], None, &cancel).is_err());
        assert_eq!(snapshot(&actual), before); assert_eq!(actual.system().last_frame(), Some(&frame));
        assert_eq!(actual.system().last_obstacle_frames().unwrap(), normal);
        assert_eq!(actual.stop_station_fractions().unwrap(), &[0.2, 0.43]);
        assert_eq!(actual.system().obstacle_offsets_m().unwrap(), &[0.0, 0.0]);
        for (i, column) in columns.iter().enumerate() { assert_eq!(actual.system().obstacle_normal_column(i).unwrap(), column); }
        // Rejected motion is not a queued command; valid retries retain DG and
        // match a realization that never saw the invalid endpoint.
        for station in [0.43001, 0.43002, 0.43002, 0.43] {
            assert_eq!(actual.step_with_stop_stations(&zero, &drives, &[0.2, station], Some(&[0.0, 0.0])).unwrap(),
                reference.step_with_stop_stations(&zero, &drives, &[0.2, station], Some(&[0.0, 0.0])).unwrap());
            assert_eq!(snapshot(&actual), snapshot(&reference));
            assert!(actual.system().last_frame().unwrap().nonlinear.is_some());
        }
    });
}
