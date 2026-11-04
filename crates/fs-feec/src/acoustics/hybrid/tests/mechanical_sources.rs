//! A real structural source must survive loaded scattering and hybrid closure.

use super::*;
use crate::acoustics::thermoviscous::ThermoviscousThermodynamics;
use crate::acoustics::thermoviscous::moving_wall::MovingWall;
use crate::highorder::simplex::SimplexSpace;
use fs_modal::time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};

#[test]
fn nonzero_wall_force_reaches_helmholtz_and_closes_original_source_work() {
    let size = [0.004, 0.0008, 0.001];
    let (hm, hp) = box_mesh([2, 2, 2], size, 0.0);
    let (lm, lp) = box_mesh([2, 2, 2], size, size[0]);
    let hf = HelmholtzFe::new(&hm, &hp, 2, medium()).unwrap();
    let lf = HelmholtzFe::new(&lm, &lp, 2, medium()).unwrap();
    let hs = shapes(&hf, 0.0, size[0], size[1], false, false);
    let ls = shapes(&lf, size[0], size[0], size[1], false, true);
    let h = InteriorMultiport::new(&hf, &complete(&hf), &hs, false, &|| false).unwrap();
    let lb = boundary(&lf, &ls, false);
    let moving_faces = lb
        .no_slip_faces
        .iter()
        .copied()
        .filter(|face| face.iter().all(|&v| lp[v as usize][1] == size[1]))
        .collect();

    // Same conforming hierarchical trace as the original moving-wall battery:
    // only the top face's interior vertex moves; stationary intersections and
    // all hierarchical edge coefficients remain zero.
    let space = SimplexSpace::new(&lm, 2);
    let mut phi = vec![[0.0; 3]; space.ndof];
    for (i, point) in lp.iter().enumerate() {
        if point[1] == size[1]
            && point[0] > size[0]
            && point[0] < 2.0 * size[0]
            && point[2] > 0.0
            && point[2] < size[2]
        {
            phi[i][1] = 100.0;
        }
    }
    let natural_omega = 3000.0;
    let damping_ratio = 0.02;
    let structure = ModalAcousticTimeModel::try_new_with_free_coordinates(
        48_000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: natural_omega,
            damping_ratio,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let l = ThermoviscousFe::new_with_moving_wall(
        &lm,
        &lp,
        gas(),
        ThermoviscousThermodynamics::ideal_gas(gas()).unwrap(),
        &lb,
        MovingWall {
            faces: moving_faces,
            velocity_shapes: vec![phi],
            structure,
        },
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap();
    let hybrid = HybridInterior::new(&h, &l, &[1], &[0], &|| false).unwrap();
    assert!((hybrid.interface().flow_map()[0] - 1.0).abs() < 1e-14);
    let frequency_hz = 700.0;
    let cp = C64::from_re(gas().cp_j_kg_k);
    let characteristic = medium().density_kg_m3 * medium().sound_speed_m_s / (size[1] * size[2]);
    let hr = [0.7 * characteristic, 1.3 * characteristic];
    let lr = [0.9 * characteristic, 1.1 * characteristic];
    let force = C64::new(1.0, 0.3);
    let prepared = hybrid
        .prepare(frequency_hz, cp, &hr, &lr, &|| false)
        .unwrap();
    let actual = prepared
        .response_with_sources(&[C64::ZERO; 2], &[force], &[], &|| false)
        .unwrap();

    // Independent physical-port oracle. Unknowns are Q_H0, Q_H1, P_L0,
    // P_L1. The structural drive contributes Q_force to Q_L=Y_L P_L+Q_force.
    // External loads give P+RQ=0; the interface gives P_H1=P_L0,
    // Q_H1+Q_L0=0. Neither a hybrid nor a scattering matrix enters this solve.
    let unloaded = l
        .prepare_with_heat_capacity(frequency_hz, cp, &mut || false)
        .unwrap();
    let source_only = unloaded
        .response_with_wall_forces(&[C64::ZERO; 2], &[force], &mut || false)
        .unwrap();
    let z = h.impedance(frequency_hz, &|| false).unwrap();
    let y = l
        .admittance_with_heat_capacity(frequency_hz, cp, &mut || false)
        .unwrap();
    let a = [
        z[0] + C64::from_re(hr[0]),
        z[1],
        C64::ZERO,
        C64::ZERO,
        z[2],
        z[3],
        -C64::ONE,
        C64::ZERO,
        C64::ZERO,
        C64::ONE,
        y[0],
        y[1],
        C64::ZERO,
        C64::ZERO,
        y[2].scale(lr[1]),
        C64::ONE + y[3].scale(lr[1]),
    ];
    let mut oracle = [
        C64::ZERO,
        C64::ZERO,
        -source_only.inward_flow_m3_s[0],
        -source_only.inward_flow_m3_s[1].scale(lr[1]),
    ];
    lu_complex_with_cancel(&a, 4, || false)
        .unwrap()
        .solve_with_cancel(&mut oracle, || false)
        .unwrap();
    let expected = unloaded
        .response_with_wall_forces(&oracle[2..], &[force], &mut || false)
        .unwrap();
    for i in 0..2 {
        assert!(relative(actual.helmholtz.inward_flow_m3_s[i], oracle[i]) < 2e-8);
        assert!(relative(actual.thermoviscous.pressure_drive_pa[i], oracle[2 + i]) < 2e-8);
    }
    for (observed, reference) in [
        (&actual.thermoviscous.pressure_pa, &expected.pressure_pa),
        (&actual.thermoviscous.temperature_k, &expected.temperature_k),
        (
            &actual.thermoviscous.inward_flow_m3_s,
            &expected.inward_flow_m3_s,
        ),
    ] {
        assert_eq!(observed.len(), reference.len());
        let scale = reference.iter().map(|z| z.abs()).fold(1e-30_f64, f64::max);
        assert!(
            observed
                .iter()
                .zip(reference)
                .all(|(&a, &b)| (a - b).abs() < 2e-8 * scale)
        );
    }
    let u = actual
        .thermoviscous
        .moving_wall
        .modal_velocity_m_sqrt_kg_per_s[0];
    let expected_u = expected.moving_wall.modal_velocity_m_sqrt_kg_per_s[0];
    assert!(u.abs() > 0.0 && relative(u, expected_u) < 2e-8);
    let omega = std::f64::consts::TAU * frequency_hz;
    let structural_impedance = C64::new(
        2.0 * damping_ratio * natural_omega,
        omega - natural_omega * natural_omega / omega,
    );
    let reaction = actual.thermoviscous.moving_wall.fluid_force_n_per_sqrt_kg[0];
    assert!(relative(force + reaction, structural_impedance * u) < 2e-8);

    // With no incident wave, the wall alone supplies emission and both owners'
    // physical losses. Its work must be counted exactly once, independently of
    // the negative external-port input and the cancelling interface work.
    let source_work = 0.5 * (force.conj() * u).re;
    let emitted = 0.5 * actual.outgoing.iter().map(|b| b.norm_sq()).sum::<f64>();
    let power = actual.power;
    assert!(source_work > 0.0 && emitted > 0.0);
    assert!(
        actual.outgoing[0].abs() > 0.0,
        "the Helmholtz region must receive the source"
    );
    assert!(power.helmholtz_interface_w > 0.0 && power.thermoviscous_interface_w < 0.0);
    assert!(power.thermoviscous.structural_damping_w > 0.0);
    assert_eq!(power.thermoviscous.heat_source_work_w, 0.0);
    let scale = source_work + emitted + power.dissipated_w();
    assert!((power.source_input_w - source_work).abs() < 2e-9 * scale);
    assert!((actual.thermoviscous.moving_wall.external_input_w - source_work).abs() < 2e-9 * scale);
    assert!((power.external_input_w + emitted).abs() < 2e-9 * scale);
    assert!(power.interface_residual_w().abs() < 2e-9 * scale);
    assert!((source_work - emitted - power.dissipated_w()).abs() < 2e-8 * scale);
    assert!(power.balance_residual_w().abs() < 2e-8 * scale);
}
