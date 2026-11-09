//! Pressure-feedback paths on original, heterogeneous accepted bulk-hardness fields.

use super::*;
use fs_contact::normal_patch::plastic_source::history_field::integration::{
    HistoryHardnessIntegrationError, WearIntegrationControls, WearIntegrationReceipt,
};

fn accurate() -> WearIntegrationControls {
    WearIntegrationControls {
        absolute_depth_tolerance_m: 1.0e-10,
        absolute_pressure_tolerance_pa: 10.0,
        relative_tolerance: 1.0e-5,
        maximum_accepted_substeps: 4096,
        maximum_step_evaluations: 12000,
    }
}
fn closes(receipt: &WearIntegrationReceipt) {
    let scale = receipt.elastic_energy_change_j.abs() + receipt.plastic_work_j.abs()
        + receipt.normal_work_j.abs() + receipt.parameter_work_j.abs();
    assert!(receipt.closure_residual_j.abs() <= 128.0 * f64::EPSILON * scale.max(1.0e-30));
    assert!(receipt.maximum_accepted_normalized_error <= 1.0);
}
fn relative_close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * expected.abs(),
        "{actual:e} != {expected:e}");
}

#[test]
fn g1_adaptive_mapped_wear_preserves_finite_contact_and_converges_to_exponential() {
    with_cx(|cx, _| {
        let fixture = Fixture::new(false);
        let fields = [fixture.field(0, &[0.01]), fixture.field(1, &[0.02])];
        let own_beta = std::array::from_fn::<_, 2, _>(|side| {
            K[side] / fields[side].hardness_pa(0).unwrap()
        });
        let beta = own_beta[0] + own_beta[1];
        let pristine = fixture.owner(fields.clone(), vec![[0, 0]; 25], false, cx).unwrap();
        let c = compliance(&pristine);
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: AREA / c,
            position_m: 2.0 * c * 0.004 / AREA,
        };
        for boundary in [
            Control::Approach(5.0e-5),
            Control::Drive { drive, controls: load_controls() },
        ] {
            let before = pristine.advance_normal(boundary, cx).unwrap().candidate;
            let saved = before.clone();
            let p0 = before.contact_state().pressure_pa()[12];
            let c_effective = match boundary {
                Control::Drive { drive, .. } => c + AREA / drive.stiffness_n_per_m,
                _ => c,
            };
            let slip = local_slip(2.0 * c_effective / beta);
            let coarse = before.advance_wear(&slip, boundary, cx).unwrap();
            assert_eq!(coarse.candidate.contact_state().pressure_pa()[12], 0.0);
            let mut loose_controls = accurate();
            loose_controls.absolute_depth_tolerance_m = 1.0e-8;
            loose_controls.absolute_pressure_tolerance_pa = 1000.0;
            loose_controls.relative_tolerance = 1.0e-3;
            let loose = before.integrate_wear(&slip, boundary, loose_controls, cx).unwrap();
            let fine = before.integrate_wear(&slip, boundary, accurate(), cx).unwrap();
            // With one loaded cell, dz/ds = beta*p and dp/dz = -1/C_eff.
            // Native indentation is retained, so p(s)=p0 exp(-beta*s/C_eff).
            let exact_pressure = p0 * (-2.0_f64).exp();
            let actual_pressure = fine.candidate.contact_state().pressure_pa()[12];
            relative_close(actual_pressure, exact_pressure, 0.01);
            assert!((actual_pressure - exact_pressure).abs()
                < (loose.candidate.contact_state().pressure_pa()[12] - exact_pressure).abs());
            assert!(fine.receipt.accepted_substeps > loose.receipt.accepted_substeps);
            assert!(fine.receipt.rejected_intervals > 0);
            let exact_recession = c_effective * p0 * (1.0 - (-2.0_f64).exp());
            for side in 0..2 {
                let expected = exact_recession * own_beta[side] / beta;
                relative_close(fine.candidate.recession_m()[side][12], expected, 0.01);
                relative_close(fine.candidate.volumes_m3()[side], AREA * expected, 0.01);
                relative_close(fine.receipt.removed_volume_m3[side],
                    fine.candidate.volumes_m3()[side], 1.0e-11);
                relative_close(fine.receipt.removed_volume_m3[side],
                    own_beta[side] * fine.receipt.pressure_weighted_sliding_j, 1.0e-11);
                assert_eq!(fine.candidate.recession_m()[side][0], 0.0);
            }
            if let Control::Drive { drive, .. } = boundary {
                let a0 = before.contact_state().approach_m();
                let a1 = fine.candidate.contact_state().approach_m();
                let expected_work = 0.5 * drive.stiffness_n_per_m
                    * ((drive.position_m - a0).powi(2) - (drive.position_m - a1).powi(2));
                relative_close(fine.receipt.normal_work_j, expected_work, 1.0e-10);
            } else {
                assert_eq!(fine.receipt.normal_work_j, 0.0);
            }
            assert_eq!(fine.candidate.contact_state().plastic_depth_m(),
                before.contact_state().plastic_depth_m());
            assert_eq!(fine.receipt.plastic_work_j, 0.0);
            assert_eq!(fine.candidate.hardness_fields(), &fields);
            assert_eq!(fine.candidate.cell_to_material_point(), before.cell_to_material_point());
            assert_eq!(fine.source_identity, before.source_identity());
            assert_eq!(fine.receipt.controls, accurate());
            closes(&fine.receipt);
            assert_eq!(before, saved);
        }
    });
}

#[test]
fn g1_nonuniform_native_hardness_and_slip_follow_coupled_maintained_load_oracle() {
    with_cx(|cx, _| {
        let fixture = Fixture::new(false);
        let fields = [
            fixture.field(0, &[0.0, 0.02]),
            fixture.field(1, &[0.0, 0.01]),
        ];
        let mut map = vec![[0, 0]; 25];
        map[12] = [0, 1];
        map[13] = [1, 0];
        let pristine = fixture.owner(fields.clone(), map.clone(), true, cx).unwrap();
        let width = pristine.model().config().grid.cell_width_m;
        let modulus = pristine.model().config().reduced_modulus_pa;
        // Independent rectangular Boussinesq integral at the adjacent cell centre.
        let neighbor = width * (3.0 * (1.0_f64 / 3.0).asinh()
            + 3.0_f64.asinh() - 2.0 * 1.0_f64.asinh())
            / (std::f64::consts::PI * modulus);
        let difference = compliance(&pristine) - neighbor;
        let force = 0.006;
        let boundary = Control::Load { force_n: force, controls: load_controls() };
        let before = pristine.advance_normal(boundary, cx).unwrap().candidate;
        let saved = before.clone();
        let own_beta = [12, 13].map(|cell| {
            std::array::from_fn::<_, 2, _>(|side| {
                K[side] / fields[side].hardness_pa(map[cell][side]).unwrap()
            })
        });
        let beta = own_beta.map(|sides| sides[0] + sides[1]);
        let distance = [2.0 * difference / beta[0], difference / beta[1]];
        let mut slip = local_slip(distance[0]);
        slip[13] = distance[1];
        let result = before.integrate_wear(&slip, boundary, accurate(), cx).unwrap();

        // At maintained load p1+p2=P. With a=beta1*s1, b=beta2*s2,
        // (z1-z2)'=P(a-b)/2-(a+b)(z1-z2)/(2*(Cself-Cneighbor)).
        let p_sum = force / AREA;
        let a = beta[0] * distance[0];
        let b = beta[1] * distance[1];
        let decay = (a + b) / (2.0 * difference);
        let equilibrium = p_sum * b / (a + b);
        let pressure0 = equilibrium + (0.5 * p_sum - equilibrium) * (-decay).exp();
        let integral0 = equilibrium
            + (0.5 * p_sum - equilibrium) * (1.0 - (-decay).exp()) / decay;
        let integrals = [integral0, p_sum - integral0];
        relative_close(result.candidate.contact_state().pressure_pa()[12], pressure0, 0.01);
        relative_close(result.candidate.contact_state().pressure_pa()[13], p_sum - pressure0, 0.01);
        relative_close(
            AREA * result.candidate.contact_state().pressure_pa().iter().sum::<f64>(),
            force, 1.0e-9,
        );
        for side in 0..2 {
            let mut expected_volume = 0.0;
            for (i, cell) in [12, 13].into_iter().enumerate() {
                let expected = own_beta[i][side] * distance[i] * integrals[i];
                relative_close(result.candidate.recession_m()[side][cell], expected, 0.01);
                expected_volume += AREA * expected;
            }
            relative_close(result.candidate.volumes_m3()[side], expected_volume, 0.01);
            relative_close(result.receipt.removed_volume_m3[side],
                result.candidate.volumes_m3()[side], 1.0e-11);
        }
        let normal_work = force * (result.candidate.contact_state().approach_m()
            - before.contact_state().approach_m());
        relative_close(result.receipt.normal_work_j, normal_work, 1.0e-10);
        assert_eq!(result.candidate.hardness_fields(), &fields);
        assert_eq!(result.candidate.cell_to_material_point(), map);
        closes(&result.receipt);
        assert_eq!(before, saved);
    });
}

#[test]
fn g4_adaptive_history_wear_discards_partial_paths_and_admits_zero_path_boundary() {
    with_cx(|cx, gate| {
        let fixture = Fixture::new(false);
        let before = fixture.owner(
            [fixture.field(0, &[0.01]), fixture.fixed[1].clone()],
            vec![[0, 0]; 25], false, cx,
        ).unwrap().advance_normal(Control::Approach(5.0e-5), cx).unwrap().candidate;
        let saved = before.clone();
        let boundary = Control::Approach(before.contact_state().approach_m());
        let slip = local_slip(1.0e-3);
        let mut limited = accurate();
        limited.maximum_step_evaluations = 2;
        assert!(matches!(
            before.integrate_wear(&slip, boundary, limited, cx),
            Err(HistoryHardnessIntegrationError::BudgetExceeded {
                step_evaluations: 2, ..
            })
        ));
        limited = accurate();
        limited.maximum_accepted_substeps = 2;
        assert!(matches!(
            before.integrate_wear(&slip, boundary, limited, cx),
            Err(HistoryHardnessIntegrationError::BudgetExceeded {
                accepted_substeps: 2, ..
            })
        ));
        let zero = vec![0.0; 25];
        let unchanged = before.integrate_wear(&zero, boundary, accurate(), cx).unwrap();
        assert_eq!(unchanged.receipt.accepted_substeps, 1);
        assert_eq!(unchanged.receipt.step_evaluations, 1);
        assert_eq!(unchanged.receipt.removed_volume_m3, [0.0; 2]);
        assert_eq!(unchanged.candidate.recession_m(), before.recession_m());
        assert_eq!(unchanged.candidate.volumes_m3(), before.volumes_m3());
        assert!(before.integrate_wear(&zero, Control::Approach(4.0e-5), accurate(), cx).is_err());
        let mut bad = slip.clone();
        bad[24] = -1.0;
        assert!(before.integrate_wear(&bad, boundary, accurate(), cx).is_err());
        gate.request();
        assert!(before.integrate_wear(&slip, boundary, accurate(), cx).is_err());
        assert_eq!(before, saved);
    });
}
