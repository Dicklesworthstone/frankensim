//! Full mechanical recoil changes actual stick/slip and own-hardness wear.

#[path = "../sliding/field_fixture.rs"]
mod field_fixture;

use super::*;
use field_fixture::{controls, motion};
use fs_contact::normal_patch::plastic_source::wear::sliding::hardness::SourceHistorySlidingWear;
use fs_contact::normal_patch::plastic_source::wear::sliding::hardness::coupled::{
    CoupledHardnessSlidingError, HardnessSlidingCoupling,
};
use fs_tribo::FrictionRegime;
use fs_tribo::tangential::matrix::MatrixTangentialConfig;
use fs_tribo::tangential::multiple::{MultipleTangentialConfig, PreparedTangentialContacts};

fn fixture() -> (
    Fixture,
    fs_contact::interface_binding::thermal_friction::SourceThermalCoulomb,
) {
    let mut fixture = Fixture::new(false);
    let (normal, coefficients, friction) =
        field_fixture::bind(&fixture.normal, &fixture.keys, [1.0e4, 4.0e6]);
    fixture.normal = normal;
    fixture.coefficients = coefficients;
    (fixture, friction)
}
fn response(matrix: Vec<f64>, evaluations: usize) -> PreparedTangentialContacts {
    PreparedTangentialContacts::new(
        2,
        matrix,
        MultipleTangentialConfig {
            maximum_contacts: 2,
            maximum_matrix_terms: 16,
            maximum_sweeps: 256,
            maximum_law_evaluations: evaluations,
            local: MatrixTangentialConfig {
                maximum_evaluations: 512,
                displacement_tolerance_m: 2.0e-13,
                energy_tolerance_j: 2.0e-13,
                matrix_relative_tolerance: 1.0e-12,
            },
            total_energy_tolerance_j: 4.0e-13,
        },
        || false,
    )
    .unwrap()
}
fn diagonal() -> Vec<f64> {
    (0..16)
        .map(|i| if i / 4 == i % 4 { 0.01 } else { 0.0 })
        .collect()
}
fn rotated_response(axis: [f64; 2]) -> Vec<f64> {
    let transverse = [-axis[1], axis[0]];
    let mut matrix = vec![0.0; 16];
    for i in 0..2 {
        for j in 0..2 {
            let longitudinal = if i == j { 0.01 } else { 0.008 };
            let lateral = if i == j { 0.02 } else { 0.003 };
            for a in 0..2 {
                for b in 0..2 {
                    matrix[(2 * i + a) * 4 + 2 * j + b] =
                        longitudinal * axis[a] * axis[b] + lateral * transverse[a] * transverse[b];
                }
            }
        }
    }
    matrix
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn g1_cross_cell_recoil_prevents_fictitious_slip_and_wear_in_both_tangent_frames() {
    with_cx(|cx, _| {
        let (fixture, friction) = fixture();
        let fields = [fixture.field(0, &[0.0]), fixture.field(1, &[0.01, 0.02])];
        let mut map = vec![[0, 0]; 25];
        map[13] = [0, 1];
        let before = SourceHistorySlidingWear::new(
            fixture
                .owner(fields.clone(), map.clone(), true, cx)
                .unwrap(),
            friction,
        )
        .unwrap()
        .advance_normal(Control::Approach(5.0e-5), cx)
        .unwrap()
        .candidate;
        let saved = before.clone();
        let normal = before.wear_state().contact_state().pressure_pa()[12] * AREA;
        near(normal, 0.01, 1.0e-14);
        let force1 = 0.4 * normal;
        let force2 = (7.0e-5 - 0.008 * force1) / 0.01;
        let actual1 = 1.0e-4 - 0.01 * force1 - 0.008 * force2;
        assert!(force2.abs() < 0.6 * normal && actual1 > 0.0);
        for axis in [[1.0, 0.0], [0.6, 0.8]] {
            let matrix = rotated_response(axis);
            let coupling =
                HardnessSlidingCoupling::new(vec![12, 13], response(matrix.clone(), 65536), cx)
                    .unwrap();
            let released = [1.0e-4, 7.0e-5].map(|s| [s * axis[0], s * axis[1]]);
            let next = before
                .advance_coupled(
                    &coupling,
                    &released,
                    0.3,
                    controls(),
                    Control::Approach(5.0e-5),
                    cx,
                )
                .unwrap();
            let shear = &next.sliding.shear;
            assert_eq!(shear.contacts[0].regime, FrictionRegime::Sliding);
            assert_eq!(shear.contacts[1].regime, FrictionRegime::Sticking);
            assert_eq!(next.sliding.unsigned_slip_m[13], 0.0);
            near(next.sliding.unsigned_slip_m[12], actual1, 1.0e-12);
            for (site, force) in [force1, force2].into_iter().enumerate() {
                for a in 0..2 {
                    near(shear.contacts[site].traction_n[a], force * axis[a], 1.0e-10);
                    near(
                        next.candidate.tangential_states()[[12, 13][site]].slip_m()[a],
                        if site == 0 { actual1 * axis[a] } else { 0.0 },
                        1.0e-12,
                    );
                }
            }
            for side in 0..2 {
                let h = fields[side].hardness_pa(map[12][side]).unwrap();
                let depth = K[side] * (normal / AREA) * actual1 / h;
                near(next.recession_increment_m[side][12], depth, 1.0e-7 * depth);
                near(
                    next.wear[side].increment_m3,
                    AREA * depth,
                    1.0e-7 * AREA * depth,
                );
                assert_eq!(next.recession_increment_m[side][13], 0.0);
            }
            near(shear.dissipated_j, force1 * actual1, 1.0e-12);
            near(
                next.sliding.pressure_weighted_sliding_j,
                normal * actual1,
                1.0e-12,
            );
            let forces: Vec<_> = shear.contacts.iter().flat_map(|c| c.traction_n).collect();
            let free: Vec<_> = released.into_iter().flatten().collect();
            let mut work = 0.0;
            for row in 0..4 {
                let actual =
                    free[row] - (0..4).map(|j| matrix[row * 4 + j] * forces[j]).sum::<f64>();
                near(
                    actual,
                    shear.contacts[row / 2].slip_increment_m[row % 2],
                    4.0e-13,
                );
                work += forces[row] * actual;
            }
            near(work, shear.dissipated_j, 4.0e-13);
            assert!(shear.body_work_residual_j.abs() <= 4.0e-13);
            assert_eq!(next.sliding.friction_sources.len(), 1);
            assert_eq!(next.sliding.cell_source_indices, vec![Some(0), Some(0)]);
            assert_eq!(
                next.candidate.tangential_states()[0],
                before.tangential_states()[0]
            );
            assert_eq!(next.candidate.wear_state().hardness_fields(), &fields);
            assert_eq!(next.candidate.wear_state().cell_to_material_point(), map);
            assert_eq!(next.source_identity, before.source_identity());
            assert_eq!(before, saved);
        }
        // Dropping the off-diagonal recoil falsely predicts sliding at site 2.
        let mut independent = motion();
        independent[12].released_increment_m = [1.0e-4, 0.0];
        independent[13].released_increment_m = [7.0e-5, 0.0];
        independent[12].compliance_m_per_n = 0.01;
        independent[13].compliance_m_per_n = 0.01;
        let wrong = before
            .advance(&independent, 0.3, controls(), Control::Approach(5.0e-5), cx)
            .unwrap();
        assert!(wrong.recession_increment_m[0][13] > 0.0);
    });
}

#[test]
fn g3_coupled_diagonal_reduction_retains_all_original_normal_boundaries() {
    with_cx(|cx, _| {
        let (fixture, friction) = fixture();
        let fields = [fixture.field(0, &[0.01]), fixture.fixed[1].clone()];
        let before = SourceHistorySlidingWear::new(
            fixture.owner(fields, vec![[0, 0]; 25], true, cx).unwrap(),
            friction,
        )
        .unwrap()
        .advance_normal(
            Control::Load {
                force_n: 0.006,
                controls: load_controls(),
            },
            cx,
        )
        .unwrap()
        .candidate;
        let approach = before.wear_state().contact_state().approach_m();
        let coupling =
            HardnessSlidingCoupling::new(vec![12, 13], response(diagonal(), 65536), cx).unwrap();
        let released = [[1.0e-4, 0.0], [0.0, -8.0e-5]];
        for control in [
            Control::Approach(approach),
            Control::Load {
                force_n: 0.006,
                controls: load_controls(),
            },
            Control::Drive {
                drive: FiniteGapPlasticDrive {
                    stiffness_n_per_m: 1000.0,
                    position_m: approach + 0.006 / 1000.0,
                },
                controls: load_controls(),
            },
        ] {
            let mut scalar = motion();
            for (site, cell) in [12, 13].into_iter().enumerate() {
                scalar[cell].released_increment_m = released[site];
                scalar[cell].compliance_m_per_n = 0.01;
            }
            let reference = before
                .advance(&scalar, 0.2, controls(), control, cx)
                .unwrap();
            let coupled = before
                .advance_coupled(&coupling, &released, 0.2, controls(), control, cx)
                .unwrap();
            let policy = coupling.mechanical_response().config();
            let slip_tolerance = 2.0 * policy.local.displacement_tolerance_m;
            let normal_tolerance = 4.0
                * before
                    .wear_state()
                    .model()
                    .config()
                    .complementarity_tolerance_m;
            let depth_tolerance = slip_tolerance.max(normal_tolerance);
            let pressure_tolerance = 4.0 * depth_tolerance / compliance(before.wear_state());
            let energy_tolerance =
                policy.total_energy_tolerance_j + 2.0 * policy.local.energy_tolerance_j;
            let volume_tolerance = AREA * 25.0 * depth_tolerance;
            let actual = coupled.candidate.wear_state();
            let expected = reference.candidate.wear_state();
            assert_eq!(coupled.source_identity, reference.source_identity);
            assert_eq!(
                actual.source_identity(),
                before.wear_state().source_identity()
            );
            assert_eq!(
                actual.hardness_fields(),
                before.wear_state().hardness_fields()
            );
            assert_eq!(
                actual.cell_to_material_point(),
                before.wear_state().cell_to_material_point()
            );
            assert_eq!(actual.normal_source(), before.wear_state().normal_source());
            assert_eq!(
                actual.coefficient_source(),
                before.wear_state().coefficient_source()
            );
            assert_eq!(
                actual.model().config().hardness_pa,
                before.wear_state().model().config().hardness_pa
            );
            near(
                actual.contact_state().approach_m(),
                expected.contact_state().approach_m(),
                depth_tolerance,
            );
            for cell in 0..25 {
                near(
                    actual.contact_state().pressure_pa()[cell],
                    expected.contact_state().pressure_pa()[cell],
                    pressure_tolerance,
                );
                near(
                    actual.contact_state().plastic_depth_m()[cell],
                    expected.contact_state().plastic_depth_m()[cell],
                    depth_tolerance,
                );
                for side in 0..2 {
                    near(
                        actual.recession_m()[side][cell],
                        expected.recession_m()[side][cell],
                        depth_tolerance,
                    );
                    near(
                        coupled.recession_increment_m[side][cell],
                        reference.recession_increment_m[side][cell],
                        depth_tolerance,
                    );
                }
            }
            for side in 0..2 {
                near(
                    actual.volumes_m3()[side],
                    expected.volumes_m3()[side],
                    volume_tolerance,
                );
                near(
                    coupled.wear[side].increment_m3,
                    reference.wear[side].increment_m3,
                    volume_tolerance,
                );
                near(
                    coupled.wear[side].total_volume_m3,
                    reference.wear[side].total_volume_m3,
                    volume_tolerance,
                );
            }
            for (site, cell) in [12, 13].into_iter().enumerate() {
                let a = &coupled.candidate.tangential_states()[cell];
                let b = &reference.candidate.tangential_states()[cell];
                assert_eq!(a.contact_temperature_k(), b.contact_temperature_k());
                assert_eq!(
                    coupled.sliding.shear.contacts[site].regime,
                    reference.sliding.cells[cell].response.regime
                );
                for axis in 0..2 {
                    near(a.slip_m()[axis], b.slip_m()[axis], slip_tolerance);
                    near(
                        coupled.sliding.shear.contacts[site].traction_n[axis],
                        reference.sliding.cells[cell].response.traction_n[axis],
                        slip_tolerance / 0.01 + 2.0 * controls().force_absolute_tolerance_n,
                    );
                }
            }
            let a = &coupled.normal.normal;
            let b = &reference.normal[0].normal;
            assert_eq!(a.kind, b.kind);
            assert_eq!(
                std::mem::discriminant(&coupled.normal.boundary),
                std::mem::discriminant(&reference.normal[0].boundary)
            );
            for (actual, expected) in [
                (a.normal_work_j, b.normal_work_j),
                (a.elastic_energy_change_j, b.elastic_energy_change_j),
                (a.parameter_work_j, b.parameter_work_j),
                (a.response.plastic_work_j, b.response.plastic_work_j),
            ] {
                near(actual, expected, energy_tolerance);
            }
            near(
                coupled.sliding.shear.dissipated_j,
                reference.sliding.frictional_work_j,
                energy_tolerance,
            );
        }
    });
}

#[test]
fn g4_registration_native_solver_and_late_wear_refuse_without_partial_state() {
    with_cx(|cx, gate| {
        let (fixture, friction) = fixture();
        let before = SourceHistorySlidingWear::new(
            fixture
                .owner(
                    [fixture.field(0, &[0.0]), fixture.fixed[1].clone()],
                    vec![[0, 0]; 25],
                    true,
                    cx,
                )
                .unwrap(),
            friction,
        )
        .unwrap()
        .advance_normal(Control::Approach(5.0e-5), cx)
        .unwrap()
        .candidate;
        let saved = before.clone();
        assert!(
            HardnessSlidingCoupling::new(vec![12, 12], response(diagonal(), 65536), cx).is_err()
        );
        let bad =
            HardnessSlidingCoupling::new(vec![12, 25], response(diagonal(), 65536), cx).unwrap();
        let control = Control::Approach(5.0e-5);
        let released = [[1.0e-4, 0.0], [7.0e-5, 0.0]];
        assert!(
            before
                .advance_coupled(&bad, &released, 0.1, controls(), control, cx)
                .is_err()
        );
        let limited =
            HardnessSlidingCoupling::new(vec![12, 13], response(diagonal(), 1), cx).unwrap();
        assert!(matches!(
            before.advance_coupled(&limited, &released, 0.1, controls(), control, cx),
            Err(CoupledHardnessSlidingError::Shear(_))
        ));
        let coupling =
            HardnessSlidingCoupling::new(vec![12, 13], response(diagonal(), 65536), cx).unwrap();
        assert!(matches!(
            before.advance_coupled(&coupling, &[[1.0, 0.0]; 2], 0.1, controls(), control, cx),
            Err(CoupledHardnessSlidingError::Source(_))
        ));
        assert_eq!(before, saved);
        gate.request();
        assert!(
            before
                .advance_coupled(&coupling, &released, 0.1, controls(), control, cx)
                .is_err()
        );
        assert_eq!(before, saved);
    });
}
