//! Retained twist and omitted bending share their original receiving body.
use super::*;

fn rotation(rest: f64) -> StringBodyRotationConnection {
    let mut link = rotational_link(ModalRotationAttachment {
        component: 0,
        angular_shapes: vec![2.0, -1.0],
    });
    link.station_fraction = 0.73;
    link.rest_angle_rad = rest;
    link
}

fn bank(
    cx: &Cx<'_>,
    gate: &CancelGate,
    rotation_first: bool,
    rest: f64,
    torque: bool,
    max_connections: usize,
) -> Result<ReciprocalStringBank, ReciprocalStringError> {
    let original = two_axis_bank(cx, gate, 7800.0, directions([0.6, 0.8], [1.0, 0.0]));
    let mut config = options();
    config.network.max_connections = max_connections;
    let mut bank = ReciprocalStringBank::new(
        cx,
        original.source().strings().to_vec(),
        gas(),
        original.system().components()[1..].to_vec(),
        law(),
        config,
        gate,
    )?;
    let mut second = second_axis(TransverseDirection::try_new([0.6, 0.8]).unwrap());
    second.station_fraction = 0.67;
    if rotation_first {
        if torque {
            bank = bank.with_rotational_body_connections(vec![rotation(rest)], gate)?;
        }
        bank.with_translational_body_connections(vec![second], gate)?
            .with_residual_flexibility(gate)
    } else {
        bank = bank
            .with_residual_flexibility(gate)?
            .with_translational_body_connections(vec![second], gate)?;
        if torque {
            bank = bank.with_rotational_body_connections(vec![rotation(rest)], gate)?;
        }
        Ok(bank)
    }
}

#[test]
fn g0_residual_torque_binding_orders_preserve_rest_energy_and_native_column_budgets() {
    with_cx(|cx, gate| {
        let rest = 0.001;
        let mut before = bank(cx, gate, true, rest, true, 4).unwrap();
        let mut after = bank(cx, gate, false, rest, true, 4).unwrap();
        assert!(before.source().residual_flexibility());
        assert_eq!(before.source().rotational_body_connections(), &[rotation(rest)]);
        assert_eq!(
            before.source().rotational_body_connections(),
            after.source().rotational_body_connections()
        );
        near(
            before.system().total_energy_j().unwrap(),
            0.5 * rotation(rest).stiffness_n_m_per_rad * rest * rest,
            1e-22,
        );
        assert_eq!(
            before.system().total_energy_j().unwrap(),
            after.system().total_energy_j().unwrap()
        );
        let released = FrictionDrive {
            speed_m_s: 0.0,
            normal_force_n: 0.0,
        };
        let zero = vec![0.0; before.system().mode_count()];
        for _ in 0..32 {
            let a = *before.step(&zero, released).unwrap();
            let b = *after.step(&zero, released).unwrap();
            assert_eq!(a, b);
            assert_eq!(a.external_work_j, 0.0);
            assert_eq!(a.drive_work_j, 0.0);
            assert!(a.energy_residual_j.abs() <= a.energy_tolerance_j);
            for (a, b) in before
                .system()
                .components()
                .iter()
                .zip(after.system().components())
            {
                assert_eq!(a.states(), b.states());
            }
        }
        // Three native factors describe the nonorthogonal bending supports;
        // the unchanged torque spring/dashpot consumes its own fourth column.
        for rotation_first in [false, true] {
            assert!(bank(cx, gate, rotation_first, 0.0, true, 3).is_err());
        }
    });
}

#[test]
fn g1_torque_only_excitation_reaches_both_planes_through_the_shared_body_after_release() {
    with_cx(|cx, gate| {
        let mut coupled = multi(bank(cx, gate, false, 0.0, true, 4).unwrap(), gate)
            .with_bow_residual_flexibility(gate)
            .unwrap();
        let mut detached = multi(bank(cx, gate, false, 0.0, false, 4).unwrap(), gate)
            .with_bow_residual_flexibility(gate)
            .unwrap();
        let card = &coupled.cards()[0];
        let radius = 0.0003_f64;
        let rotary_inertia = 7800.0 * core::f64::consts::PI * radius.powi(4) / 2.0;
        let normalization = (0.5 * card.length_m * rotary_inertia).sqrt();
        let mut force = vec![0.0; coupled.system().mode_count()];
        for n in 1..=3 {
            force[7 + n] = 2e-8 * (n as f64 * core::f64::consts::PI * 0.3).sin()
                / normalization;
        }
        let released = [FrictionDrive {
            speed_m_s: 0.0,
            normal_force_n: 0.0,
        }];
        let (mut pressure, mut x, mut y) = (0.0_f64, 0.0_f64, 0.0_f64);
        let (mut work, mut loss) = (0.0, 0.0);
        for sample in 0..192 {
            if sample == 32 {
                force.fill(0.0);
            }
            let a = coupled.step(&force, &released).unwrap().clone();
            let b = detached.step(&force, &released).unwrap().clone();
            assert_eq!(b.observer_pressure_pa, 0.0);
            assert!(detached.system().components()[0].states()[..8]
                .iter()
                .all(|q| q.displacement_m_sqrt_kg == 0.0 && q.velocity_m_sqrt_kg_per_s == 0.0));
            assert_eq!(a.drive_work_j, 0.0);
            assert!(a.energy_residual_j.abs() <= a.energy_tolerance_j);
            work += a.external_work_j;
            loss += a.network_dissipation_j;
            if sample >= 32 {
                assert_eq!(a.external_work_j, 0.0);
                pressure = pressure.max(a.observer_pressure_pa.abs());
                let state = coupled.system().components()[0].states();
                x = x.max(state[..4].iter().map(|q| q.displacement_m_sqrt_kg.abs()).fold(0.0, f64::max));
                y = y.max(state[4..8].iter().map(|q| q.displacement_m_sqrt_kg.abs()).fold(0.0, f64::max));
            }
        }
        assert!(pressure > 0.0 && x > 0.0 && y > 0.0);
        near(coupled.system().total_energy_j().unwrap() + loss, work, 1e-15);
    });
}

#[test]
fn g1_torque_support_preserves_bow_twist_and_normal_bending_residual_columns() {
    with_cx(|cx, gate| {
        let make = |torque| bank(cx, gate, false, 0.0, torque, 4).unwrap();
        let original = make(true).system().port().left.shapes[8..].to_vec();
        let bowed = multi(make(true), gate).with_bow_residual_flexibility(gate).unwrap();
        let unbound = multi(make(false), gate).with_bow_residual_flexibility(gate).unwrap();
        assert_eq!(bowed.system().port(0).unwrap().left.shapes[8..], original);
        assert!(original.iter().any(|x| *x != 0.0));
        assert_eq!(
            bowed.system().friction_port_column(0),
            unbound.system().friction_port_column(0)
        );
        assert_eq!(
            bowed.system().friction_residual_compliance(),
            unbound.system().friction_residual_compliance()
        );
        let mut normal = stop(TransverseDirection::try_new([0.6, 0.8]).unwrap());
        normal.location.station_fraction = 0.2;
        let mut stopped = multi(make(true), gate)
            .with_string_stops(vec![normal.clone()], stop_limits(), gate).unwrap();
        let reference = multi(make(false), gate)
            .with_string_stops(vec![normal], stop_limits(), gate).unwrap();
        assert_eq!(stopped.system().obstacle_normal_column(0), reference.system().obstacle_normal_column(0));
        assert_eq!(stopped.system().obstacle_residual_compliance(), reference.system().obstacle_residual_compliance());
        let normal = stopped.system().obstacle_normal_column(0).unwrap();
        let bow = bowed.system().friction_port_column(0).unwrap();
        assert_eq!(&normal[..8], &bow[..8]);
        assert_eq!(&normal[11..], &bow[11..]);
        assert!(normal[8..11].iter().all(|x| *x == 0.0));
        let mut external = vec![0.0; stopped.system().mode_count()];
        external[8] = 0.002;
        let released = [FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 }];
        for sample in 0..64 {
            if sample == 16 { external.fill(0.0); }
            let frame = stopped.step(&external, &released).unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    });
}
