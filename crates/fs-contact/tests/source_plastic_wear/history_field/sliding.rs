//! Actual-slip coupling preserves each original native material-point hardness.

#[path = "../sliding/field_fixture.rs"]
mod field_fixture;

use super::*;
use field_fixture::{controls, motion, normal_controls};
use fs_contact::normal_patch::plastic_source::wear::sliding::hardness::{
    HardnessSlidingError, SourceHistorySlidingWear,
};

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

#[test]
fn g1_actual_sliding_wears_each_mapped_body_with_its_own_native_history_hardness() {
    with_cx(|cx, _| {
        let (fixture, friction) = fixture();
        let fields = [
            fixture.field(0, &[0.0, 0.02]),
            fixture.field(1, &[0.0, 0.01]),
        ];
        let mut map = vec![[0, 0]; 25];
        map[12] = [0, 1];
        map[13] = [1, 0];
        let owner = fixture
            .owner(fields.clone(), map.clone(), true, cx)
            .unwrap();
        let before = SourceHistorySlidingWear::new(owner, friction)
            .unwrap()
            .advance_normal(Control::Approach(5.0e-5), cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let mut input = motion();
        input[12].released_increment_m = [1.0e-5, 0.0];
        input[13].released_increment_m = [0.0, 2.0e-5];
        let next = before
            .advance(&input, 0.3, controls(), Control::Approach(5.0e-5), cx)
            .unwrap();
        let mut measure = 0.0;
        for (cell, distance) in [(12, 1.0e-5), (13, 2.0e-5)] {
            let pressure = before.wear_state().contact_state().pressure_pa()[cell];
            close(next.sliding.unsigned_slip_m[cell], distance);
            measure += pressure * AREA * distance;
            for side in 0..2 {
                let h = fields[side].hardness_pa(map[cell][side]).unwrap();
                close(
                    next.recession_increment_m[side][cell],
                    K[side] * pressure * distance / h,
                );
            }
            let source_index = next.sliding.cells[cell].source_index.unwrap();
            assert_eq!(
                next.sliding.friction_sources[source_index].query().axes()["p"],
                pressure
            );
        }
        for side in 0..2 {
            close(
                next.wear[side].increment_m3,
                AREA * (next.recession_increment_m[side][12]
                    + next.recession_increment_m[side][13]),
            );
            close(
                next.candidate.wear_state().volumes_m3()[side],
                next.wear[side].increment_m3,
            );
        }
        close(next.sliding.pressure_weighted_sliding_j, measure);
        close(next.sliding.frictional_work_j, 0.4 * measure);
        assert_eq!(next.normal.len(), 1);
        assert_eq!(
            next.normal[0].normal.predecessor_state_identity,
            before.wear_state().contact_state().identity()
        );
        assert_eq!(
            next.normal[0].normal.successor_state_identity,
            next.candidate.wear_state().contact_state().identity()
        );
        assert_eq!(next.candidate.wear_state().hardness_fields(), &fields);
        assert_eq!(next.candidate.wear_state().cell_to_material_point(), map);
        assert_eq!(next.sliding.cells[0].source_index, None);
        assert_eq!(before, saved);
    });
}

#[test]
fn g1_field_refresh_preserves_shear_and_reversal_wears_with_the_updated_body_hardness() {
    with_cx(|cx, _| {
        let (fixture, friction) = fixture();
        let initial_fields = [fixture.field(0, &[0.01]), fixture.fixed[1].clone()];
        let owner = fixture
            .owner(initial_fields.clone(), vec![[0, 0]; 25], false, cx)
            .unwrap();
        let force = 0.004;
        let boundary = Control::Load {
            force_n: force,
            controls: normal_controls(),
        };
        let before = SourceHistorySlidingWear::new(owner, friction)
            .unwrap()
            .advance_normal(boundary, cx)
            .unwrap()
            .candidate;
        let mut input = motion();
        let distance = 1.0e-5;
        input[12].released_increment_m = [distance, 0.0];
        let forward = before
            .advance(&input, 0.1, controls(), boundary, cx)
            .unwrap();
        let fields = [fixture.field(0, &[0.02]), fixture.fixed[1].clone()];
        let refreshed = forward
            .candidate
            .refresh_hardness(fields.clone(), boundary, cx)
            .unwrap();
        assert_eq!(
            refreshed.candidate.tangential_states(),
            forward.candidate.tangential_states()
        );
        assert_eq!(
            refreshed.candidate.wear_state().volumes_m3(),
            forward.candidate.wear_state().volumes_m3()
        );
        assert_eq!(
            refreshed.candidate.wear_state().recession_m(),
            forward.candidate.wear_state().recession_m()
        );
        let saved = refreshed.candidate.clone();
        assert!(matches!(
            refreshed.candidate.refresh_hardness(
                [fixture.field(0, &[0.0]), fixture.fixed[1].clone()],
                boundary,
                cx,
            ),
            Err(HardnessSlidingError::History(_))
        ));
        assert_eq!(refreshed.candidate, saved);
        let retained = refreshed
            .candidate
            .advance_normal(boundary, cx)
            .unwrap()
            .candidate;
        assert_eq!(
            retained.tangential_states(),
            refreshed.candidate.tangential_states()
        );
        input[12].released_increment_m[0] = -distance;
        let reverse = retained
            .advance(&input, 0.1, controls(), boundary, cx)
            .unwrap();
        assert_eq!(reverse.candidate.tangential_states()[12].slip_m(), [0.0; 2]);
        assert!(reverse.wear[0].increment_m3 < forward.wear[0].increment_m3);
        for side in 0..2 {
            let h0 = initial_fields[side].hardness_pa(0).unwrap();
            let h1 = fields[side].hardness_pa(0).unwrap();
            close(
                reverse.candidate.wear_state().volumes_m3()[side],
                K[side] * force * distance * (1.0 / h0 + 1.0 / h1),
            );
        }
        close(
            forward.sliding.frictional_work_j + reverse.sliding.frictional_work_j,
            2.0 * 0.4 * force * distance,
        );
        assert_eq!(forward.normal.len(), 1);
        assert_eq!(reverse.normal.len(), 1);
    });
}
