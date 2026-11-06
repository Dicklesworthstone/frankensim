//! Native stress and growth endpoints share one selected ligament limit.
use super::*;
use fs_material::fatigue::crack_growth::{CentreCrackTensileError, FiniteWidthCentreCrack};
use fs_material::strength::IsotropicTensileLimit;

fn bound_gauge(source: &ParisSourceFixture) -> TetCrackGrowthGauge {
    let mut selected = gauge();
    selected.material_state_identity = source.bundle.identity().0;
    selected.material_source = Some(source.bundle.clone());
    selected.law = source.law.clone();
    selected.geometry = FiniteWidthCentreCrack::new(0.01)
        .unwrap()
        .tabulate(vec![0.0005, 0.001, 0.0025, 0.004], 4)
        .unwrap()
        .into();
    selected
}
fn strength(source: &ParisSourceFixture, name: &str) -> IsotropicTensileLimit {
    IsotropicTensileLimit::resolve(
        &source.card,
        &source.bundle.query(),
        name,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn native_overload_and_failed_source_binding_preserve_samples_time_and_growth() {
    with_cx(false, |cx| {
        let source = paris_source_fixture(0);
        let foreign = paris_source_fixture(1);
        let fixed = fixed();
        let problem = problem(&source.material, &fixed);
        let geometry = geometry(cx);
        let mut path = path(&problem, &source.thermal, vec![bound_gauge(&source)], cx);
        assert!(
            path.set_net_section_tensile_limit(
                0,
                strength(&foreign, "selected_net_section_limit"),
                cx
            )
            .is_err()
        );
        assert!(path.net_section_tensile_limit(0).is_none());
        with_cx(true, |cancelled| {
            assert!(matches!(
                path.set_net_section_tensile_limit(
                    0,
                    strength(&source, "selected_net_section_limit"),
                    cancelled
                ),
                Err(TetCrackGrowthError::Cancelled)
            ));
        });
        path.set_net_section_tensile_limit(0, strength(&source, "selected_net_section_limit"), cx)
            .unwrap();
        let guard = path.net_section_tensile_limit(0).unwrap();
        assert_eq!(
            guard.limit().resolved().card_identity(),
            source.card.content_hash()
        );
        path.advance(0.0, &load(&geometry, 2e6, cx), cx).unwrap();
        let before = path.states();
        assert!(matches!(
            path.advance(1.0, &load(&geometry, 10e6, cx), cx),
            Err(TetCrackGrowthError::NetSectionTensile {
                gauge: 0,
                error: CentreCrackTensileError::LimitReached(_)
            })
        ));
        assert_eq!(path.states(), before);
        assert_eq!(path.observed_samples(), 1);
        assert_eq!(path.last_time_s(), Some(0.0));
        path.advance(1.0, &load(&geometry, 8e6, cx), cx).unwrap();
    });
}

#[test]
fn net_limit_forecast_and_late_gauge_crossing_preserve_the_reusable_cycle() {
    with_cx(false, |cx| {
        let source = paris_source_fixture(0);
        let fixed = fixed();
        let problem = problem(&source.material, &fixed);
        let geometry = geometry(cx);
        let selected = bound_gauge(&source);
        let mut path = path(
            &problem,
            &source.thermal,
            vec![selected.clone(), selected],
            cx,
        );
        path.set_net_section_tensile_limit(1, strength(&source, "selected_net_section_limit"), cx)
            .unwrap();
        for (i, stress) in [2e6, 8e6, 2e6].into_iter().enumerate() {
            path.advance(i as f64, &load(&geometry, stress, cx), cx)
                .unwrap();
        }
        let before = path.states();
        let forecast = path.cycles_to_net_section_tensile_limit(1, cx).unwrap();
        close(forecast.target_length_m, 0.005 * (1.0 - 8e6 / 12e6), 1e-15);
        assert_eq!(path.states(), before);
        let crossing = forecast.additional_cycles.ceil() as u64;
        assert!(crossing > 1 && crossing < 100);
        assert!(matches!(
            path.close_repeated_cycle(crossing, cx),
            Err(TetCrackGrowthError::NetSectionTensile {
                gauge: 1,
                error: CentreCrackTensileError::LimitReached(_)
            })
        ));
        assert_eq!(path.states(), before);
        assert_eq!(path.observed_samples(), 3);
        assert_eq!(path.last_time_s(), Some(2.0));
        assert_eq!(
            path.cycles_to_net_section_tensile_limit(1, cx).unwrap(),
            forecast
        );
        path.close_repeated_cycle(crossing - 1, cx).unwrap();
        assert!(path.states()[1].crack_length_m < forecast.target_length_m);
        assert_eq!(path.states()[0], path.states()[1]);
        // Previously consumed unguarded history cannot acquire a guard later.
        assert!(
            path.set_net_section_tensile_limit(
                0,
                strength(&source, "selected_net_section_limit"),
                cx
            )
            .is_err()
        );
    });
}

#[test]
fn transient_initial_limit_and_peak_overload_refuse_without_advancing_motion() {
    with_cx(false, |cx| {
        let source = paris_source_fixture(0);
        let fixed = fixed();
        let problem = problem(&source.material, &fixed);
        let geometry = geometry(cx);
        let forcing = load(&geometry, 8e6, cx);
        let mut initial = ZERO;
        initial[1][0] = 2e6 / 200e9;
        let native = forcing
            .begin_transient(
                &problem,
                Some(TetThermalStrainField::Uniform(&source.thermal)),
                0.0,
                &initial,
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let mut original = native.clone();
        let mut observed = native
            .with_crack_growth(FRAME, vec![bound_gauge(&source)], config(), cx)
            .unwrap();
        assert!(matches!(
            observed.set_net_section_tensile_limit(
                0,
                strength(&source, "low_net_section_limit"),
                cx
            ),
            Err(TetCrackGrowthError::NetSectionTensile {
                gauge: 0,
                error: CentreCrackTensileError::LimitReached(_)
            })
        ));
        assert!(observed.net_section_tensile_limit(0).is_none());
        assert_eq!(observed.mechanics().state(), original.state());
        observed
            .set_net_section_tensile_limit(0, strength(&source, "selected_net_section_limit"), cx)
            .unwrap();
        // Newmark's pi/2 oscillator rotation: 2 -> 8 -> 14 MPa. The second
        // sampled peak crosses the 9.6 MPa current gross-section capacity.
        let h = 2.0 / (10.0_f64 * 200e9 / 7800.0).sqrt();
        let expected = original.advance(&forcing, h, cx).unwrap();
        let accepted = observed.advance(&forcing, h, cx).unwrap();
        assert_eq!(accepted.mechanics, expected);
        close(accepted.normal_stress_pa[0], 8e6, 1e-5);
        let before = observed.states();
        let response = observed.recover_response(cx).unwrap();
        assert!(matches!(
            observed.advance(&forcing, h, cx),
            Err(TetCrackGrowthError::NetSectionTensile {
                gauge: 0,
                error: CentreCrackTensileError::LimitReached(_)
            })
        ));
        assert_eq!(observed.mechanics().state(), original.state());
        assert_eq!(observed.recover_response(cx).unwrap(), response);
        assert_eq!(observed.states(), before);
        assert_eq!(observed.observed_samples(), 2);
    });
}
