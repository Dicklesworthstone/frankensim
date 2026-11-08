//! G0/G1/G3/G4 physical loading-spring integration, using authored coefficients.
use super::*;
use fs_contact::normal_patch::{
    PlasticLayerShearDrive, PlasticLayerShearDriveConfig, PlasticLayerShearDriveError,
};

fn config(length: f64) -> PlasticLayerShearDriveConfig {
    PlasticLayerShearDriveConfig {
        stiffness_n_per_m: [
            [3000.0, 150.0, 80.0],
            [150.0, 2500.0, -100.0],
            [80.0, -100.0, 2000.0],
        ],
        spin_length_m: length,
        maximum_rotation_rad: 0.1,
        generalized_force_tolerance_n: 1e-10,
        energy_tolerance_j: 2e-11,
    }
}

fn fixture(
    cells: &[usize],
    laws: Vec<ThirdBodyLayerContact>,
    matrix: Vec<f64>,
    spring: PlasticLayerShearDriveConfig,
    cx: &Cx<'_>,
) -> PlasticLayerShearDrive {
    let contact = FiniteGapPlasticLayerContact::new(
        normal(cells, cx),
        cells
            .iter()
            .copied()
            .zip(laws)
            .map(|(cell_index, law)| FiniteGapPlasticLayerSite { cell_index, law })
            .collect(),
        response(cells.len(), matrix),
        cx,
    )
    .unwrap();
    PlasticLayerShearDrive::new(contact, spring, cx).unwrap()
}

fn cross_response() -> Vec<f64> {
    vec![
        0.0004, 0.00005, 0.00012, -0.00003, 0.00005, 0.0003, 0.00002, 0.00008, 0.00012, 0.00002,
        0.0005, -0.00004, -0.00003, 0.00008, -0.00004, 0.00035,
    ]
}
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn spring_energy(k: [[f64; 3]; 3], e: [f64; 3]) -> f64 {
    0.5 * dot3(e, k.map(|row| dot3(row, e)))
}

/// Independent seven-unknown force/carriage solve. This keeps the spring as
/// three extra equilibrium equations; it does NOT form K^-1 or condense it.
fn elastic_reference(
    body: &[f64],
    k: [[f64; 3]; 3],
    stiffness: [f64; 2],
    a: [[f64; 3]; 4],
    old_layer: [f64; 4],
    old_spring: [f64; 3],
    stage: [f64; 3],
    free: [f64; 4],
) -> Vec<f64> {
    let mut matrix = vec![0.0; 49];
    let mut rhs = vec![0.0; 7];
    for i in 0..4 {
        rhs[i] = 2.0 * old_layer[i] + free[i];
        for j in 0..4 {
            matrix[7 * i + j] = body[4 * i + j];
        }
        matrix[7 * i + i] += 2.0 / stiffness[i / 2];
        for j in 0..3 {
            matrix[7 * i + 4 + j] = -a[i][j];
        }
    }
    for i in 0..3 {
        rhs[4 + i] = dot3(
            k[i],
            core::array::from_fn(|j| old_spring[j] + 0.5 * stage[j]),
        );
        for j in 0..4 {
            matrix[7 * (4 + i) + j] = a[j][i];
        }
        for j in 0..3 {
            matrix[7 * (4 + i) + 4 + j] = 0.5 * k[i][j];
        }
    }
    solve(matrix, rhs)
}

#[test]
fn coupled_translation_twist_and_full_body_recoil_match_uncondensed_equilibrium() {
    with_cx(|cx, _| {
        let spring = config(1e-3);
        let body = cross_response();
        let stiffness = [1000.0, 1500.0];
        let drive = fixture(
            &[11, 13],
            stiffness.into_iter().map(elastic).collect(),
            body.clone(),
            spring,
            cx,
        );
        let mut before = drive
            .advance_normal(&drive.initial_state(), 1e-5, cx)
            .unwrap()
            .candidate;
        let a = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 1.0],
        ];
        let mut old_layer = [0.0; 4];
        let mut old_spring = [0.0; 3];
        let mut carriage = [0.0; 3];
        for (stage, free) in [
            ([2e-4, 1e-4, 0.03], [[2e-5, -1e-5], [-1e-5, 3e-5]]),
            ([-1e-4, 0.0, -0.02], [[0.0; 2]; 2]),
            ([0.0; 3], [[1e-5, 0.0], [0.0, -2e-5]]),
        ] {
            let saved = before.clone();
            let scaled_stage = [stage[0], stage[1], stage[2] * spring.spin_length_m];
            let flat = [free[0][0], free[0][1], free[1][0], free[1][1]];
            let expected = elastic_reference(
                &body,
                spring.stiffness_n_per_m,
                stiffness,
                a,
                old_layer,
                old_spring,
                scaled_stage,
                flat,
            );
            let step = drive
                .advance(&before, stage, &free, 0.001, controls(), cx)
                .unwrap();
            assert_eq!(before, saved);
            assert_eq!(
                step.candidate.contact_state().normal_state(),
                before.contact_state().normal_state()
            );
            let force: Vec<_> = step
                .shear
                .contacts
                .iter()
                .flat_map(|s| s.traction_n)
                .collect();
            let mut body_work = 0.0;
            let mut layer_change = 0.0;
            for i in 0..4 {
                close(force[i], expected[i], 2e-10);
                let recoil: f64 = (0..4).map(|j| body[4 * i + j] * expected[j]).sum();
                body_work += expected[i] * (recoil - flat[i]);
                let slip = dot3(a[i], [expected[4], expected[5], expected[6]]) + flat[i] - recoil;
                let new = old_layer[i] + slip;
                layer_change += 0.5 * stiffness[i / 2] * (new * new - old_layer[i] * old_layer[i]);
                close(
                    step.candidate.contact_state().layer_states()[i / 2].elastic_slip_m()[i % 2],
                    new,
                    5e-13,
                );
                old_layer[i] = new;
            }
            let new_spring =
                core::array::from_fn(|j| old_spring[j] + scaled_stage[j] - expected[4 + j]);
            let energy = spring_energy(spring.stiffness_n_per_m, new_spring)
                - spring_energy(spring.stiffness_n_per_m, old_spring);
            let wrench = [
                expected[0] + expected[2],
                expected[1] + expected[3],
                1e-3 * (expected[3] - expected[1]),
            ];
            close(step.stage_work_j, dot3(wrench, stage), 2e-13);
            close(step.spring_stored_change_j, energy, 2e-13);
            close(step.shear.stored_change_j, layer_change, 2e-13);
            close(step.body_response_work_j, body_work, 2e-13);
            close(step.shear.dissipated_j, 0.0, 1e-15);
            close(step.balance_residual_j, 0.0, spring.energy_tolerance_j);
            for j in 0..3 {
                carriage[j] += expected[4 + j];
                let length = if j == 2 { spring.spin_length_m } else { 1.0 };
                close(
                    step.carriage_increment_m_rad[j],
                    expected[4 + j] / length,
                    5e-10,
                );
                close(
                    step.candidate.carriage_position_m_rad()[j],
                    carriage[j] / length,
                    5e-10,
                );
                close(
                    step.candidate.spring_extension_m_rad()[j],
                    new_spring[j] / length,
                    5e-10,
                );
                close(step.mean_wrench_n_n_m[j], wrench[j], 2e-10);
                close(step.spring_equilibrium_residual_n[j], 0.0, 1e-10);
            }
            close(
                step.candidate.spring_energy_j(),
                spring_energy(spring.stiffness_n_per_m, new_spring),
                2e-13,
            );
            old_spring = new_spring;
            before = step.candidate;
        }
    });
}

#[test]
fn stationary_stage_relaxes_then_reversal_retains_plastic_maxwell_and_kelvin_memory() {
    with_cx(|cx, _| {
        let (k, km, eta, kelvin, tau, yield_force, ks) =
            (1000.0_f64, 500.0, 10.0, 0.05_f64, 0.1, 0.2, 2000.0);
        let law = ThirdBodyLayerContact::new_tensor(
            ShearStiffness::positive_definite([[k, 0.0], [0.0, 1500.0]]).unwrap(),
            vec![TensorLayerRelaxationArm {
                stiffness: ShearStiffness::positive_definite([[km, 0.0], [0.0, 700.0]]).unwrap(),
                relaxation_time_s: tau,
            }],
            yield_force,
            eta,
            300.0,
            0.02,
        )
        .unwrap()
        .with_kelvin_viscosity(
            ShearViscosity::from_gram_rows(vec![[kelvin.sqrt(), 0.0], [0.0, 0.3]]).unwrap(),
        )
        .unwrap();
        let mut spring = config(1e-3);
        spring.stiffness_n_per_m = [[ks, 0.0, 0.0], [0.0, 3000.0, 0.0], [0.0, 0.0, 4000.0]];
        let drive = fixture(&[12], vec![law], vec![0.0; 4], spring, cx);
        let mut before = drive
            .advance_normal(&drive.initial_state(), 1e-5, cx)
            .unwrap()
            .candidate;
        let (mut e, mut r, mut p, mut es) = (0.0, 0.0, 0.0, 0.0);
        let mut maximum_plastic = 0.0_f64;
        for (index, (stage, dt)) in [(0.001, 0.001), (0.0, 0.002), (0.0, 0.004), (-0.0025, 0.002)]
            .into_iter()
            .enumerate()
        {
            // Independent one-axis midpoint recurrence with explicitly eliminated
            // Maxwell branch. The plastic branch has a closed scalar force solve.
            let weight = 1.0 / (1.0 + 0.5 * dt / tau);
            let h = k * e + km * weight * r;
            let a = 0.5 * (k + km * weight) + kelvin / dt;
            let flow = dt / eta;
            let trial = (h + a * (stage + 2.0 * es)) / (1.0 + 2.0 * a / ks);
            let sign = trial.signum();
            let magnitude = if trial.abs() > yield_force {
                (h + a * (stage + 2.0 * es) + sign * a * flow * yield_force)
                    / (1.0 + 2.0 * a / ks + a * flow)
            } else {
                trial
            };
            let dp = if magnitude.abs() > yield_force {
                flow * (magnitude.abs() - yield_force) * sign
            } else {
                0.0
            };
            let next_es = 2.0 * magnitude / ks - es;
            let du = stage - (next_es - es);
            let de = du - dp;
            let next_e = e + de;
            let next_r = 2.0 * weight * (r + 0.5 * de) - r;
            let next_p = p + dp;
            let plastic = yield_force * dp.abs() + eta * dp * dp / dt;
            let relaxation = dt / tau * km * (0.5 * r + 0.5 * next_r).powi(2);
            let viscous = kelvin * de * de / dt;
            let spring_change = 0.5 * ks * (next_es * next_es - es * es);
            let layer_change =
                0.5 * k * (next_e * next_e - e * e) + 0.5 * km * (next_r * next_r - r * r);
            let saved = before.clone();
            let step = drive
                .advance(&before, [stage, 0.0, 0.0], &[[0.0; 2]], dt, controls(), cx)
                .unwrap();
            assert_eq!(before, saved);
            let state = &step.candidate.contact_state().layer_states()[0];
            close(state.elastic_slip_m()[0], next_e, 2e-12);
            close(state.plastic_slip_m()[0], next_p, 2e-12);
            close(state.relaxation_slip_m()[0][0], next_r, 2e-12);
            close(step.mean_wrench_n_n_m[0], magnitude, 2e-9);
            close(step.candidate.spring_extension_m_rad()[0], next_es, 2e-12);
            close(step.dissipation[0].plastic_j, plastic, 2e-12);
            close(step.dissipation[0].relaxation_j, relaxation, 2e-12);
            close(step.dissipation[0].kelvin_j, viscous, 2e-12);
            close(step.shear.stored_change_j, layer_change, 2e-12);
            close(
                step.shear.dissipated_j,
                plastic + relaxation + viscous,
                2e-12,
            );
            close(step.spring_stored_change_j, spring_change, 2e-12);
            close(step.stage_work_j, magnitude * stage, 2e-12);
            close(step.body_response_work_j, 0.0, 1e-15);
            close(step.balance_residual_j, 0.0, spring.energy_tolerance_j);
            if index == 1 || index == 2 {
                assert!(spring_change < 0.0 && plastic > 0.0 && relaxation > 0.0);
                close(step.stage_work_j, 0.0, 0.0);
            }
            maximum_plastic = maximum_plastic.max(next_p);
            if index == 3 {
                assert!(magnitude < 0.0 && next_p < maximum_plastic);
                assert!(next_p != 0.0); // reversal does not erase the accumulated history
            }
            (e, r, p, es) = (next_e, next_r, next_p, next_es);
            before = step.candidate;
        }
    });
}

#[test]
fn changing_the_authored_spin_length_preserves_physical_forces_motion_and_work() {
    with_cx(|cx, _| {
        let first = config(1e-3);
        let scale = 2.3;
        let mut second = first;
        second.spin_length_m *= scale;
        for i in 0..3 {
            for j in 0..3 {
                if i == 2 {
                    second.stiffness_n_per_m[i][j] /= scale;
                }
                if j == 2 {
                    second.stiffness_n_per_m[i][j] /= scale;
                }
            }
        }
        let a = fixture(
            &[11, 13],
            vec![elastic(1000.0), elastic(1500.0)],
            cross_response(),
            first,
            cx,
        );
        let b = fixture(
            &[11, 13],
            vec![elastic(1000.0), elastic(1500.0)],
            cross_response(),
            second,
            cx,
        );
        let mut sa = a
            .advance_normal(&a.initial_state(), 1e-5, cx)
            .unwrap()
            .candidate;
        let mut sb = b
            .advance_normal(&b.initial_state(), 1e-5, cx)
            .unwrap()
            .candidate;
        for stage in [[1e-4, 2e-4, 0.04], [-5e-5, -1e-4, -0.025]] {
            let aa = a
                .advance(&sa, stage, &[[0.0; 2]; 2], 0.002, controls(), cx)
                .unwrap();
            let bb = b
                .advance(&sb, stage, &[[0.0; 2]; 2], 0.002, controls(), cx)
                .unwrap();
            for j in 0..3 {
                close(aa.mean_wrench_n_n_m[j], bb.mean_wrench_n_n_m[j], 2e-10);
                close(
                    aa.candidate.carriage_position_m_rad()[j],
                    bb.candidate.carriage_position_m_rad()[j],
                    3e-10,
                );
                close(
                    aa.candidate.spring_extension_m_rad()[j],
                    bb.candidate.spring_extension_m_rad()[j],
                    3e-10,
                );
            }
            close(aa.stage_work_j, bb.stage_work_j, 1e-13);
            close(aa.spring_stored_change_j, bb.spring_stored_change_j, 1e-13);
            close(aa.body_response_work_j, bb.body_response_work_j, 1e-13);
            close(aa.shear.stored_change_j, bb.shear.stored_change_j, 1e-13);
            sa = aa.candidate;
            sb = bb.candidate;
        }
    });
}

#[test]
fn native_normal_updates_and_late_refusals_preserve_the_complete_fixture_checkpoint() {
    with_cx(|cx, gate| {
        let spring = config(1e-3);
        let drive = fixture(&[12], vec![elastic(1000.0)], vec![0.0; 4], spring, cx);
        assert!(matches!(
            drive.advance(
                &drive.initial_state(),
                [0.0; 3],
                &[[0.0; 2]],
                0.001,
                controls(),
                cx
            ),
            Err(PlasticLayerShearDriveError::Contact(
                FiniteGapPlasticLayerError::DetachedLayer { .. }
            ))
        ));
        let loaded = drive
            .advance_normal(&drive.initial_state(), 1e-5, cx)
            .unwrap()
            .candidate;
        let sheared = drive
            .advance(
                &loaded,
                [1e-4, 0.0, 0.0],
                &[[0.0; 2]],
                0.001,
                controls(),
                cx,
            )
            .unwrap()
            .candidate;
        let saved = sheared.clone();
        let target = 0.5;
        let normal = drive
            .advance_normal_at_load(&sheared, target, load_controls(), cx)
            .unwrap();
        close(normal.normal.normal.normal_force_n, target, 2e-10);
        assert_eq!(
            normal.candidate.contact_state().layer_states(),
            sheared.contact_state().layer_states()
        );
        assert_eq!(
            normal.candidate.spring_extension_m_rad(),
            sheared.spring_extension_m_rad()
        );
        assert_eq!(
            normal.candidate.carriage_position_m_rad(),
            sheared.carriage_position_m_rad()
        );
        let normal_drive = FiniteGapPlasticDrive {
            position_m: 2e-5,
            stiffness_n_per_m: 1e5,
        };
        let continued = drive
            .advance_normal_with_drive(&normal.candidate, normal_drive, load_controls(), cx)
            .unwrap();
        assert_eq!(
            continued.candidate.spring_extension_m_rad(),
            sheared.spring_extension_m_rad()
        );
        assert_eq!(
            continued.candidate.contact_state().layer_states(),
            sheared.contact_state().layer_states()
        );
        assert!(drive.advance_normal(&sheared, -1e-4, cx).is_err());
        assert!(matches!(
            drive.advance(
                &sheared,
                [0.0, 0.0, 0.2],
                &[[0.0; 2]],
                0.001,
                controls(),
                cx
            ),
            Err(PlasticLayerShearDriveError::Rotation { .. })
        ));
        let mut altered = spring;
        altered.stiffness_n_per_m[0][0] *= 2.0;
        let other = fixture(&[12], vec![elastic(1000.0)], vec![0.0; 4], altered, cx);
        assert!(matches!(
            other.advance(&sheared, [0.0; 3], &[[0.0; 2]], 0.001, controls(), cx),
            Err(PlasticLayerShearDriveError::Invalid(_))
        ));
        for bad in [
            [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[1.0, 0.1, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
        ] {
            let mut changed = spring;
            changed.stiffness_n_per_m = bad;
            assert!(PlasticLayerShearDrive::new(drive.contact().clone(), changed, cx).is_err());
        }
        gate.request();
        assert!(
            drive
                .advance(&sheared, [0.0; 3], &[[0.0; 2]], 0.001, controls(), cx)
                .is_err()
        );
        assert_eq!(sheared, saved);
    });
}
