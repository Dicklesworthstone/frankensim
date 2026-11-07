//! Actual native stick/slip, depth hardness and retained normal histories.

#[path = "../source_plastic_wear/sliding/field_fixture.rs"]
mod field_fixture;

use super::*;
use field_fixture::{controls, motion, normal_controls};
use fs_contact::normal_patch::plastic_source::wear::sliding::PlasticSlidingInput;
use fs_contact::normal_patch::plastic_source::wear::sliding::hardness::{
    HardnessSlidingControl as Control, HardnessSlidingError, SourceDepthSlidingWear,
};
use fs_contact::normal_patch::{FiniteGapPlasticContinuationKind, FiniteGapPlasticDrive};
use fs_tribo::FrictionRegime;

fn fixture(
    gradient: f64,
    support: [f64; 2],
) -> (
    Fixture,
    fs_contact::interface_binding::thermal_friction::SourceThermalCoulomb,
) {
    let mut fixture = Fixture::new(gradient);
    let (normal, coefficients, friction) =
        field_fixture::bind(&fixture.normal, &fixture.keys, support);
    fixture.normal = normal;
    fixture.coefficients = coefficients;
    (fixture, friction)
}

#[test]
fn g1_native_stick_slip_drives_original_depth_integrals_and_both_normal_accounts() {
    with_cx(|cx, _| {
        let (fixture, friction) = fixture(GRADIENT, [1.0e4, 4.0e6]);
        let pristine = fixture.owner(true, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx);
        let before = SourceDepthSlidingWear::new(pristine, friction)
            .unwrap()
            .advance_normal(Control::Approach(5.0e-5), cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let mut input = motion();
        input[12] = PlasticSlidingInput {
            released_increment_m: [1.0e-5, 0.0],
            compliance_m_per_n: 0.01,
        };
        input[13] = PlasticSlidingInput {
            released_increment_m: [1.0e-4, 0.0],
            compliance_m_per_n: 0.01,
        };
        input[0].released_increment_m = [1.0, 0.0];
        let next = before
            .advance(&input, 0.3, controls(), Control::Approach(5.0e-5), cx)
            .unwrap();
        assert_eq!(
            next.sliding.cells[12].response.regime,
            FrictionRegime::Sticking
        );
        assert_eq!(next.sliding.unsigned_slip_m[12], 0.0);
        assert_eq!(next.recession_increment_m[0][12], 0.0);
        assert_eq!(next.recession_increment_m[1][12], 0.0);
        let actual = 1.0e-4 - 0.01 * 0.4 * H0 * AREA;
        close(next.sliding.unsigned_slip_m[13], actual);
        close(
            next.sliding.cells[13].response.traction_n[0],
            0.4 * H0 * AREA,
        );
        let q = K[0] * H0 * actual;
        let a = 2.0 * q / (H0 + (H0 * H0 + 2.0 * GRADIENT * q).sqrt());
        let b = K[1] * H0 * actual / HB;
        for (side, expected) in [a, b].into_iter().enumerate() {
            close(next.recession_increment_m[side][13], expected);
            close(next.wear[side].increment_m3, AREA * expected);
            close(
                next.candidate.wear_state().volumes_m3()[side],
                AREA * expected,
            );
            assert_eq!(next.candidate.wear_state().recession_m()[side][0], 0.0);
        }
        close(next.sliding.frictional_work_j, 0.4 * H0 * AREA * actual);
        close(next.sliding.pressure_weighted_sliding_j, H0 * AREA * actual);
        assert_eq!(next.sliding.friction_sources.len(), 1);
        assert_eq!(
            next.sliding.cells[12].source_index,
            next.sliding.cells[13].source_index
        );
        assert_eq!(next.sliding.cells[0].source_index, None);
        assert_eq!(next.normal.len(), 2);
        assert_eq!(
            next.normal[0].normal.kind,
            FiniteGapPlasticContinuationKind::SurfaceRecession
        );
        assert_eq!(
            next.normal[1].normal.kind,
            FiniteGapPlasticContinuationKind::MaterialRefresh
        );
        assert_eq!(
            next.normal[0].normal.predecessor_state_identity,
            before.wear_state().contact_state().identity()
        );
        assert_eq!(
            next.normal[1].normal.predecessor_state_identity,
            next.normal[0].normal.successor_state_identity
        );
        assert_eq!(
            next.normal[1].normal.successor_state_identity,
            next.candidate.wear_state().contact_state().identity()
        );
        assert!(next.candidate.wear_state().contact_state().pressure_pa()[13] < H0);
        assert_eq!(before, saved);
    });
}

#[test]
fn g1_depth_sliding_reversal_accumulates_actual_wear_at_all_three_boundaries() {
    with_cx(|cx, _| {
        let (fixture, friction) = fixture(0.0, [1.0e4, 4.0e6]);
        let before = SourceDepthSlidingWear::new(
            fixture.owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx),
            friction,
        )
        .unwrap()
        .advance_normal(
            Control::Load {
                force_n: 0.005,
                controls: normal_controls(),
            },
            cx,
        )
        .unwrap()
        .candidate;
        let approach = before.wear_state().contact_state().approach_m();
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: 1000.0,
            position_m: approach + 0.005 / 1000.0,
        };
        for boundary in [
            Control::Approach(approach),
            Control::Load {
                force_n: 0.005,
                controls: normal_controls(),
            },
            Control::Drive {
                drive,
                controls: normal_controls(),
            },
        ] {
            let distance = 1.0e-5;
            let mut input = motion();
            input[12].released_increment_m = [distance, 0.0];
            let forward = before
                .advance(&input, 0.1, controls(), boundary, cx)
                .unwrap();
            let retained = forward
                .candidate
                .advance_normal(boundary, cx)
                .unwrap()
                .candidate;
            assert_eq!(
                retained.tangential_states(),
                forward.candidate.tangential_states()
            );
            input[12].released_increment_m[0] = -distance;
            let reverse = retained
                .advance(&input, 0.1, controls(), boundary, cx)
                .unwrap();
            assert_eq!(reverse.candidate.tangential_states()[12].slip_m(), [0.0; 2]);
            assert!(reverse.sliding.cells[12].response.traction_n[0] < 0.0);
            let p_sum = before.wear_state().contact_state().pressure_pa()[12]
                + retained.wear_state().contact_state().pressure_pa()[12];
            for (side, h) in [H0, HB].into_iter().enumerate() {
                close(
                    reverse.candidate.wear_state().volumes_m3()[side],
                    K[side] * AREA * p_sum * distance / h,
                );
                assert!(
                    forward.wear[side].increment_m3 > 0.0 && reverse.wear[side].increment_m3 > 0.0
                );
            }
            close(
                forward.sliding.frictional_work_j + reverse.sliding.frictional_work_j,
                0.4 * AREA * p_sum * distance,
            );
            assert_eq!(forward.normal.len(), 2);
            assert_eq!(reverse.normal.len(), 2);
        }
    });
}

#[test]
fn g4_depth_sliding_late_source_wear_and_cancellation_refusals_keep_shear_and_inventory() {
    with_cx(|cx, gate| {
        for support in [[1.0e4, 0.9e6], [1.0e4, 4.0e6]] {
            let (fixture, friction) = fixture(GRADIENT, support);
            let before = SourceDepthSlidingWear::new(
                fixture.owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx),
                friction,
            )
            .unwrap()
            .advance_normal(Control::Approach(5.0e-5), cx)
            .unwrap()
            .candidate;
            let saved = before.clone();
            let mut input = motion();
            input[12].released_increment_m = [1.0, 0.0];
            let error = before
                .advance(&input, 0.1, controls(), Control::Approach(5.0e-5), cx)
                .unwrap_err();
            if support[1] < H0 {
                assert!(matches!(error, HardnessSlidingError::Sliding(_)));
            } else {
                assert!(matches!(error, HardnessSlidingError::Depth(_)));
            }
            assert_eq!(before, saved);
            if support[1] > H0 {
                input[12].released_increment_m = [1.0e-5, 0.0];
                assert!(
                    before
                        .advance(&input, 0.1, controls(), Control::Approach(4.0e-5), cx)
                        .is_err()
                );
                gate.request();
                assert!(
                    before
                        .advance(&input, 0.1, controls(), Control::Approach(5.0e-5), cx)
                        .is_err()
                );
                assert_eq!(before, saved);
                return;
            }
        }
    });
}
