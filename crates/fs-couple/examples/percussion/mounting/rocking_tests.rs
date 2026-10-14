use super::*;
use super::super::{Spec, row};
use crate::{Experiment, Stroke, mallets, shaft_playing, specimen};
use fs_exec::CancelGate;

const CARD: &str = include_str!("../estimated-rocking-mount.fsmount");
fn experiment(spec: &Spec, audio: bool) -> Experiment {
    let mut shell = specimen::Specimen::reference(); shell.azimuths = 8;
    let mut result = crate::splash_with_mount(512, 2e-6, audio,
        Stroke { speed_m_s: 0.0, position_m: Some([0.068, 0.01]) }, Some(shell), &[],
        Some(Stroke { speed_m_s: 0.0, position_m: Some([-0.055, 0.012]) }), None,
        &mallets::Selection::default(), &shaft_playing::Selection::default(), None,
        Some(0.03), Some(spec)).unwrap();
    result.system = result.system.into_analytic_nonlinear().unwrap(); result
}
#[test]
fn rocking_input_is_optional_and_refuses_invalid_principal_inertia_or_resistance() {
    let spec = Spec::parse(CARD).unwrap(); assert!(spec.rocking.is_some());
    let record = "rocking,0,0,0.00008,0.00009,0.15,0.18,0.001,0.001\n";
    assert!(Spec::parse(&CARD.replace(record, "")).unwrap().rocking.is_none());
    assert!(Spec::parse(&format!("{CARD}{record}")).is_err());
    for bad in [CARD.replace("0.00008", "0"), CARD.replace("0.00009", "NaN"),
        CARD.replace("0.15,0.18,0.001", "-0.15,0.18,0.001"),
        CARD.replace("0.18,0.001,0.001", "0.18,-0.001,0.001")] {
        assert!(Spec::parse(&bad).is_err());
    }
}
#[test]
fn noncentral_contacts_produce_the_right_signed_moment_and_energy_units() {
    let support = Support { centre_xy_m: [0.1, -0.2], inertia_kg_m2: [0.04, 0.09],
        stiffness_nm_rad: [2.0, 3.0], resistance_nm_s_rad: [0.1, 0.2] };
    support.validate().unwrap();
    let row = support.vertical_row([0.13, -0.16]).unwrap();
    assert!((row[0] - 0.2).abs() < 1e-14);
    assert!((row[1] + 0.1).abs() < 1e-14);
    let angular_velocity = [0.7, -0.8]; let force = 2.5;
    let p = [angular_velocity[0] * 0.04_f64.sqrt(), angular_velocity[1] * 0.09_f64.sqrt()];
    let modal_power = force * (row[0]*p[0] + row[1]*p[1]);
    let physical_power = force * (0.04*angular_velocity[0] - 0.03*angular_velocity[1]);
    assert!((modal_power - physical_power).abs() < 1e-14);
    assert!((0.5*(p[0]*p[0]+p[1]*p[1])
        - 0.5*(0.04*angular_velocity[0].powi(2)+0.09*angular_velocity[1].powi(2))).abs() < 1e-14);
}
#[test]
fn rocking_preserves_both_sticks_shell_observers_and_acoustic_source_addresses() {
    let rocking = Spec::parse(CARD).unwrap(); let mut translating = rocking.clone(); translating.rocking = None;
    let a = experiment(&translating, true); let b = experiment(&rocking, true);
    let n = a.force.len(); assert_eq!(b.force.len(), n + 2);
    assert_eq!(a.second_stick.unwrap().coordinate, b.second_stick.unwrap().coordinate);
    assert_eq!(a.acoustics.as_ref().unwrap().state_modes(), b.acoustics.as_ref().unwrap().state_modes());
    assert_eq!(&b.observer_a[..n], &a.observer_a); assert_eq!(&b.observer_b[..n], &a.observer_b);
    assert_eq!(&b.observer_a[n..], &[0.0; 2]); assert_eq!(&b.observer_b[n..], &[0.0; 2]);
    let port = rocking.port(b.force.len()).unwrap(); assert_eq!(port.translation.coordinate, n - 1);
    assert_eq!(port.angular_motion(b.system.state()).unwrap(), Some([[0.0; 2]; 2]));
}
#[test]
fn physical_support_torque_changes_the_shell_and_respects_force_work_and_rollback() {
    let spec = Spec::parse(CARD).unwrap(); let mut a = experiment(&spec, false);
    let port = spec.port(a.force.len()).unwrap();
    let weights = port.inverse_sqrt_inertia.unwrap();
    let mut heavy = spec.clone(); heavy.rocking.as_mut().unwrap().inertia_kg_m2[0] *= 10.0;
    let mut b = experiment(&heavy, false); let other = heavy.port(b.force.len()).unwrap();
    let gate = CancelGate::new_clock_free();
    let mut supplied = 0.0; let mut balance = 0.0_f64;
    for tick in 0..512 {
        // Equal PHYSICAL applied moments, not equal mass-normalized forces.
        let torque = if tick < 64 {0.001} else {0.0};
        a.force[port.translation.coordinate + 1] = torque * weights[0];
        b.force[other.translation.coordinate + 1] = torque * other.inverse_sqrt_inertia.unwrap()[0];
        if tick == 32 {
            let before = a.system.state().to_vec(); let cancel = CancelGate::new_clock_free(); cancel.request();
            assert!(a.system.step(&a.force, &cancel).is_err()); assert_eq!(a.system.state(), before.as_slice());
        }
        let f = a.system.step(&a.force, &gate).unwrap(); b.system.step(&b.force, &gate).unwrap();
        assert_eq!(f.sample, tick + 1); supplied += f.supplied_work_j;
        balance = balance.max(f.balance_residual_j.abs());
    }
    assert!(supplied > 0.0 && balance < 1e-7);
    assert!(port.angular_motion(a.system.state()).unwrap().unwrap()[0][0].abs() > 1e-9);
    let end = 2*port.translation.coordinate;
    assert!(a.system.state()[2..end].iter().zip(&b.system.state()[2..end])
        .any(|(a,b)| (a-b).abs() > 1e-12), "inertia must affect actual shell motion, not just reporting");
    let before = a.system.state().to_vec(); assert!(a.system.step(&a.force, &gate).is_err());
    assert_eq!(a.system.state(), before.as_slice());
}
#[test]
fn angular_observers_do_not_read_felt_or_radiation_memory_as_body_motion() {
    let spec = Spec::parse(CARD).unwrap(); let port = spec.port(8).unwrap();
    let mut state = vec![0.0; 24];
    state[16..].fill(1e6); // private histories follow all eight q/p pairs.
    let scales = port.inverse_sqrt_inertia.unwrap();
    state[12] = 0.001/scales[0]; state[13] = 0.02/scales[0];
    state[14] = -0.002/scales[1]; state[15] = 0.03/scales[1];
    let angles = port.angular_motion(&state).unwrap().unwrap();
    assert!((angles[0][0] - 0.001).abs() < 1e-15);
    assert!((angles[1][0] + 0.002).abs() < 1e-15);
    let mut output = Vec::new(); port.header(&mut output).unwrap(); output.push(b'\n');
    row(port, &state, &mut output).unwrap();
    let text = String::from_utf8(output).unwrap();
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines[0].split(',').count(), lines[1].split(',').count());
    assert!(text.contains("mount_rotation_y_rad"));
    assert!(port.angular_motion(&state[..15]).is_err());
}
