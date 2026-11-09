//! Actual shear survives independently accepted depth/history continuation.

#[path = "../../source_plastic_wear/sliding/field_fixture.rs"]
mod field_fixture;

use super::refresh::{HistoryFixture, compliance, controls as normal_controls};
use super::*;
use fs_contact::normal_patch::FiniteGapPlasticDrive;
use fs_contact::normal_patch::plastic_source::wear::sliding::hardness::history::DepthSlidingHistoryRefresh;
use fs_contact::normal_patch::plastic_source::wear::sliding::hardness::{
    HardnessSlidingBoundaryReceipt as Boundary, HardnessSlidingControl as Control,
    HardnessSlidingError, SourceDepthSlidingWear,
};
use fs_contact::normal_patch::plastic_source::wear::sliding::{
    PlasticSlidingInput, SourcePlasticSlidingError,
};

fn prepared(support: [f64; 2], cx: &Cx<'_>) -> (HistoryFixture, SourceDepthSlidingWear) {
    let mut fixture = HistoryFixture::new();
    let (normal, coefficients, friction) =
        field_fixture::bind(&fixture.fixture.normal, &fixture.fixture.keys, support);
    fixture.fixture.normal = normal;
    fixture.fixture.coefficients = coefficients;
    let boundary = Control::Approach(5e-5);
    let loaded = SourceDepthSlidingWear::new(
        fixture
            .fixture
            .owner(false, FiniteGapPlasticLinearSolver::Dense, 1e-4, cx),
        friction,
    )
    .unwrap()
    .advance_normal(boundary, cx)
    .unwrap()
    .candidate;
    let mut input = field_fixture::motion();
    input[12].released_increment_m = [1e-5, 0.0];
    let worn = loaded
        .advance(&input, 0.2, field_fixture::controls(), boundary, cx)
        .unwrap()
        .candidate;
    assert_eq!(worn.tangential_states()[12].slip_m(), [1e-5, 0.0]);
    assert!(
        worn.wear_state()
            .volumes_m3()
            .iter()
            .all(|volume| *volume > 0.0)
    );
    (fixture, worn)
}

fn retained(before: &SourceDepthSlidingWear, refresh: &DepthSlidingHistoryRefresh) {
    let next = &refresh.candidate;
    assert_eq!(refresh.original_source_identity, before.source_identity());
    assert_eq!(refresh.successor_source_identity, next.source_identity());
    assert_eq!(next.tangential_states(), before.tangential_states());
    assert_eq!(next.friction_source(), before.friction_source());
    assert_eq!(
        next.interface_temperature_k(),
        before.interface_temperature_k()
    );
    assert_eq!(
        next.wear_state().volumes_m3(),
        before.wear_state().volumes_m3()
    );
    assert_eq!(
        next.wear_state().recession_m(),
        before.wear_state().recession_m()
    );
    assert_eq!(
        refresh.normal.normal.predecessor_state_identity,
        before.wear_state().contact_state().identity()
    );
    assert_eq!(
        refresh.normal.normal.successor_state_identity,
        next.wear_state().contact_state().identity()
    );
    for (cell, &pressure) in next
        .wear_state()
        .contact_state()
        .pressure_pa()
        .iter()
        .enumerate()
    {
        match refresh.cell_source_indices[cell] {
            None => assert_eq!(pressure, 0.0),
            Some(slot) => {
                assert!(pressure > 0.0);
                let source = &refresh.friction_sources[slot];
                assert_eq!(source.state_coordinate("p"), Some(pressure));
                assert_eq!(
                    source.state_coordinate("T"),
                    Some(before.interface_temperature_k())
                );
                assert_eq!(
                    source
                        .property("mu_kinetic")
                        .unwrap()
                        .answer()
                        .receipt
                        .selected,
                    before
                        .friction_source()
                        .initial_source()
                        .property("mu_kinetic")
                        .unwrap()
                        .answer()
                        .receipt
                        .selected
                );
            }
        }
    }
}

#[test]
fn g1_profile_softening_retains_genuine_shear_and_wear_and_reverse_step_uses_it() {
    with_cx(|cx, _| {
        let (fixture, before) = prepared([1e4, 4e6], cx);
        let saved = before.clone();
        let control = Control::Approach(before.wear_state().contact_state().approach_m());
        let z = before.wear_state().recession_m()[0][12];
        let p0 = before.wear_state().contact_state().pressure_pa()[12];
        let expected = fixture.stages[1].hardness_pa_at(z).unwrap();
        assert!(expected < p0);
        let soft = before
            .refresh_hardness_profiles(fixture.profiles(1), control, cx)
            .unwrap();
        retained(&before, &soft);
        assert_ne!(
            soft.successor_source_identity,
            soft.original_source_identity
        );
        close(
            soft.candidate.wear_state().model().config().hardness_pa[12],
            expected,
        );
        close(
            soft.candidate.wear_state().contact_state().pressure_pa()[12],
            expected,
        );
        let dent = compliance(before.wear_state()) * (p0 - expected);
        close(
            soft.candidate
                .wear_state()
                .contact_state()
                .plastic_depth_m()[12]
                - before.wear_state().contact_state().plastic_depth_m()[12],
            dent,
        );
        close(
            soft.normal.normal.response.plastic_work_j,
            AREA * expected * dent,
        );
        close(
            soft.normal.normal.parameter_work_j,
            -0.5 * AREA * compliance(before.wear_state()) * (p0 - expected).powi(2),
        );
        assert_eq!(soft.normal.normal.normal_work_j, 0.0);
        assert_eq!(soft.normal.boundary, Boundary::Approach);
        assert_eq!(soft.friction_sources.len(), 1);
        assert_eq!(soft.cell_source_indices[12], Some(0));
        assert_eq!(soft.cell_source_indices[0], None);

        // A real reverse step consumes the retained signed travel. Refresh itself
        // adds neither this new physical interval nor its friction work/removal.
        let mut reverse = field_fixture::motion();
        reverse[12] = PlasticSlidingInput {
            released_increment_m: [-1e-5, 0.0],
            compliance_m_per_n: 0.0,
        };
        let next = soft
            .candidate
            .advance(&reverse, 0.17, field_fixture::controls(), control, cx)
            .unwrap();
        assert_eq!(next.candidate.tangential_states()[12].slip_m(), [0.0; 2]);
        close(next.sliding.unsigned_slip_m[12], 1e-5);
        close(next.sliding.frictional_work_j, 0.4 * expected * AREA * 1e-5);
        assert!(next.wear.iter().all(|wear| wear.increment_m3 > 0.0));
        assert_eq!(before, saved);
    });
}

#[test]
fn g1_profile_rehardening_delegates_all_normal_boundaries_and_preserves_complete_shear() {
    with_cx(|cx, _| {
        let (fixture, before) = prepared([1e4, 4e6], cx);
        let saved = before.clone();
        let approach = before.wear_state().contact_state().approach_m();
        let force_n = before.wear_state().contact_state().pressure_pa()[12] * AREA;
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: 1000.0,
            position_m: approach + force_n / 1000.0,
        };
        for control in [
            Control::Approach(approach),
            Control::Load {
                force_n,
                controls: normal_controls(),
            },
            Control::Drive {
                drive,
                controls: normal_controls(),
            },
        ] {
            let profiles = fixture.profiles(2);
            let (native, boundary) = match control {
                Control::Approach(_) => (
                    before
                        .wear_state()
                        .refresh_hardness_profiles(profiles.clone(), cx)
                        .unwrap(),
                    Boundary::Approach,
                ),
                Control::Load { force_n, controls } => {
                    let step = before
                        .wear_state()
                        .refresh_hardness_profiles_at_load(profiles.clone(), force_n, controls, cx)
                        .unwrap();
                    (step.history, Boundary::Load(step.load))
                }
                Control::Drive { drive, controls } => {
                    let step = before
                        .wear_state()
                        .refresh_hardness_profiles_with_drive(profiles.clone(), drive, controls, cx)
                        .unwrap();
                    (step.history, Boundary::Drive(step.drive))
                }
            };
            let refreshed = before
                .refresh_hardness_profiles(profiles, control, cx)
                .unwrap();
            retained(&before, &refreshed);
            assert_eq!(refreshed.candidate.wear_state(), &native.candidate);
            assert_eq!(refreshed.normal.normal, native.receipt);
            assert_eq!(refreshed.normal.boundary, boundary);
            assert!(
                refreshed
                    .candidate
                    .wear_state()
                    .model()
                    .config()
                    .hardness_pa[12]
                    > before.wear_state().model().config().hardness_pa[12]
            );
            assert_eq!(
                refreshed
                    .candidate
                    .wear_state()
                    .contact_state()
                    .plastic_depth_m(),
                before.wear_state().contact_state().plastic_depth_m()
            );
            assert_eq!(before, saved);
        }
    });
}

#[test]
fn g4_late_actual_pressure_friction_refusal_leaves_full_shear_and_profile_owner_intact() {
    with_cx(|cx, gate| {
        let (fixture, before) = prepared([9.5e5, 4e6], cx);
        let saved = before.clone();
        let control = Control::Approach(before.wear_state().contact_state().approach_m());
        let native = before
            .wear_state()
            .refresh_hardness_profiles(fixture.profiles(1), cx)
            .unwrap();
        assert!(before.wear_state().contact_state().pressure_pa()[12] > 9.5e5);
        assert!(native.candidate.contact_state().pressure_pa()[12] < 9.5e5);
        // The reference sample is supported; the new actual loaded pressure is not.
        assert!(
            before
                .friction_source()
                .coefficients_at_query(&before.friction_source().initial_source().query(),)
                .is_ok()
        );
        let error = before
            .refresh_hardness_profiles(fixture.profiles(1), control, cx)
            .unwrap_err();
        assert!(matches!(
            error,
            HardnessSlidingError::Sliding(SourcePlasticSlidingError::Cell { index: 12, .. })
        ));
        assert_eq!(before, saved);
        assert!(matches!(
            before.refresh_hardness_profiles(fixture.profiles(1), Control::Approach(4e-5), cx,),
            Err(HardnessSlidingError::Depth(
                DepthHardnessWearError::Invalid(_)
            ))
        ));
        gate.request();
        assert!(
            before
                .refresh_hardness_profiles(fixture.profiles(1), control, cx)
                .is_err()
        );
        assert_eq!(before, saved);
    });
}
