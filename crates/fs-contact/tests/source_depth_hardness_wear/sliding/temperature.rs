//! Independent affine mu(T,p), H(T,z), actual slip and native boundary work.

use super::super::temperature::{COLD, HOT, thermal_fixture};
use super::*;
use fs_contact::normal_patch::plastic_source::wear::sliding::SourcePlasticSlidingError;
use fs_contact::normal_patch::plastic_source::wear::sliding::hardness::temperature::DepthSlidingTemperatureError;

fn prepared(support: [f64; 2], cx: &Cx<'_>) -> (Fixture, SourceDepthSlidingWear) {
    let mut fixture = thermal_fixture(HOT, HOT);
    let (normal, coefficients, friction) =
        field_fixture::bind_thermal(&fixture.normal, &fixture.keys, [COLD, HOT], support);
    fixture.normal = normal;
    fixture.coefficients = coefficients;
    let boundary = Control::Load {
        force_n: 0.003,
        controls: normal_controls(),
    };
    let loaded = SourceDepthSlidingWear::new(
        fixture.owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx),
        friction,
    )
    .unwrap()
    .advance_normal(boundary, cx)
    .unwrap()
    .candidate;
    let mut input = motion();
    input[12].released_increment_m = [1.0e-5, 0.0];
    let before = loaded
        .advance(&input, 0.3, controls(), boundary, cx)
        .unwrap()
        .candidate;
    (fixture, before)
}

fn cell_compliance(source: &SourceDepthHardnessWear) -> f64 {
    4.0 * source.model().config().grid.cell_width_m * (1.0 + 2.0_f64.sqrt()).ln()
        / (std::f64::consts::PI * source.model().config().reduced_modulus_pa)
}

#[test]
fn g1_thermal_shear_composes_original_profiles_actual_slip_and_all_boundary_accounts() {
    with_cx(|cx, _| {
        let (fixture, before) = prepared([1.0e4, 1.2e6], cx);
        let saved = before.clone();
        let a0 = before.wear_state().contact_state().approach_m();
        let p0 = before.wear_state().contact_state().pressure_pa()[12];
        let c0 = cell_compliance(before.wear_state());
        let zsum =
            before.wear_state().recession_m()[0][12] + before.wear_state().recession_m()[1][12];
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: 1000.0,
            position_m: a0 + 0.003 / 1000.0,
        };
        let body_temperatures = [HOT, COLD + 10.0];
        for boundary in [
            Control::Approach(a0),
            Control::Load {
                force_n: 0.003,
                controls: normal_controls(),
            },
            Control::Drive {
                drive,
                controls: normal_controls(),
            },
        ] {
            let refreshed = before
                .refresh_at_temperatures(
                    [&fixture.cards[0], &fixture.cards[1]],
                    body_temperatures,
                    HOT,
                    boundary,
                    cx,
                )
                .unwrap();
            let hot = &refreshed.candidate;
            assert_eq!(refreshed.original_source_identity, before.source_identity());
            assert_ne!(hot.source_identity(), before.source_identity());
            assert_eq!(hot.interface_temperature_k(), HOT);
            assert_eq!(hot.tangential_states(), before.tangential_states());
            assert_eq!(
                hot.tangential_states()[12].contact_temperature_k(),
                Some(COLD)
            );
            assert_eq!(
                hot.wear_state().volumes_m3(),
                before.wear_state().volumes_m3()
            );
            assert_eq!(
                hot.wear_state().recession_m(),
                before.wear_state().recession_m()
            );
            assert_eq!(
                hot.wear_state().contact_state().plastic_depth_m(),
                before.wear_state().contact_state().plastic_depth_m()
            );
            let p = hot.wear_state().contact_state().pressure_pa()[12];
            let a = hot.wear_state().contact_state().approach_m();
            let c = 1.5 * c0;
            close(cell_compliance(hot.wear_state()), c);
            let (expected_p, normal_work) = match boundary {
                Control::Approach(_) => ((a0 - zsum) / c, 0.0),
                Control::Load { force_n, .. } => (force_n / AREA, force_n * (a - a0)),
                Control::Drive { drive, .. } => (
                    (drive.position_m - zsum) / (c + AREA / drive.stiffness_n_per_m),
                    0.5 * drive.stiffness_n_per_m
                        * ((drive.position_m - a0).powi(2) - (drive.position_m - a).powi(2)),
                ),
            };
            close(p, expected_p);
            let elastic_change = 0.5 * AREA * (c * p * p - c0 * p0 * p0);
            let normal = &refreshed.normal.normal;
            close(normal.normal_work_j, normal_work);
            close(normal.elastic_energy_change_j, elastic_change);
            close(normal.parameter_work_j, elastic_change - normal_work);
            assert_eq!(normal.response.plastic_work_j, 0.0);
            assert_eq!(
                normal.predecessor_state_identity,
                before.wear_state().contact_state().identity()
            );
            assert_eq!(
                normal.successor_state_identity,
                hot.wear_state().contact_state().identity()
            );
            assert_eq!(refreshed.friction_sources.len(), 1);
            assert_eq!(refreshed.cell_source_indices[12], Some(0));
            assert_eq!(refreshed.cell_source_indices[0], None);
            let selected = &refreshed.friction_sources[0];
            assert_eq!(selected.state_coordinate("T"), Some(HOT));
            assert_eq!(selected.state_coordinate("p"), Some(p));
            let mu = 0.4 - 0.005 * (HOT - COLD) + 1.0e-8 * (p - 5.0e5);
            close(selected.property("mu_kinetic").unwrap().value_si(), mu);
            assert_eq!(
                selected
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

            // A genuine reverse step now uses the new T/p coefficient. External
            // body compliance subtracts elastic release from the wear distance.
            let mut input = motion();
            input[12] = PlasticSlidingInput {
                released_increment_m: [-4.0e-5, 0.0],
                compliance_m_per_n: 0.01,
            };
            let next = hot.advance(&input, 0.17, controls(), boundary, cx).unwrap();
            let force = mu * p * AREA;
            let actual = 4.0e-5 - 0.01 * force;
            close(next.sliding.cells[12].response.traction_n[0], -force);
            close(next.sliding.unsigned_slip_m[12], actual);
            close(next.sliding.frictional_work_j, force * actual);
            close(next.sliding.pressure_weighted_sliding_j, p * AREA * actual);
            close(
                next.candidate.tangential_states()[12].slip_m()[0],
                1.0e-5 - actual,
            );
            assert_eq!(
                next.candidate.tangential_states()[12].contact_temperature_k(),
                Some(HOT)
            );
            let sample =
                &next.sliding.friction_sources[next.sliding.cells[12].source_index.unwrap()];
            assert_eq!(sample.state_coordinate("T"), Some(HOT));
            assert_eq!(sample.state_coordinate("p"), Some(p));
            for (side, slope) in [-1.0e11, -2.0e11].into_iter().enumerate() {
                let old_depth = hot.wear_state().recession_m()[side][12];
                let h =
                    [H0, HB][side] - 2.5e4 * (body_temperatures[side] - COLD) + slope * old_depth;
                close(
                    hot.wear_state().hardness_profiles()[side]
                        .hardness_pa_at(old_depth)
                        .unwrap(),
                    h,
                );
                let q = 2.0 * K[side] * p * actual;
                let depth = 2.0 * q / (h + (h * h + 2.0 * slope * q).sqrt());
                close(next.recession_increment_m[side][12], depth);
                close(next.wear[side].increment_m3, AREA * depth);
            }
            assert_eq!(next.normal.len(), 2);
            assert_eq!(
                next.normal[0].normal.predecessor_state_identity,
                normal.successor_state_identity
            );
            assert_eq!(
                next.normal[1].normal.predecessor_state_identity,
                next.normal[0].normal.successor_state_identity
            );
            assert_eq!(before, saved);
        }
    });
}

#[test]
fn g4_thermal_shear_refuses_actual_successor_pressure_after_native_normal_accepts() {
    with_cx(|cx, gate| {
        let (fixture, before) = prepared([2.5e5, 1.2e6], cx);
        let saved = before.clone();
        let native = before
            .wear_state()
            .refresh_at_temperatures([&fixture.cards[0], &fixture.cards[1]], [HOT; 2], HOT, cx)
            .unwrap();
        assert!(native.candidate.contact_state().pressure_pa()[12] < 2.5e5);
        let reference = before
            .friction_source()
            .initial_source()
            .query()
            .with("T", HOT)
            .unwrap();
        assert!(
            before
                .friction_source()
                .coefficients_at_query(&reference)
                .is_ok()
        );
        let boundary = Control::Approach(before.wear_state().contact_state().approach_m());
        let error = before
            .refresh_at_temperatures(
                [&fixture.cards[0], &fixture.cards[1]],
                [HOT; 2],
                HOT,
                boundary,
                cx,
            )
            .unwrap_err();
        assert!(matches!(
            error,
            DepthSlidingTemperatureError::Sliding(SourcePlasticSlidingError::Cell {
                index: 12,
                ..
            })
        ));
        assert_eq!(before, saved);
        gate.request();
        assert!(
            before
                .refresh_at_temperatures(
                    [&fixture.cards[0], &fixture.cards[1]],
                    [COLD; 2],
                    COLD,
                    boundary,
                    cx,
                )
                .is_err()
        );
        assert_eq!(before, saved);
    });
}
