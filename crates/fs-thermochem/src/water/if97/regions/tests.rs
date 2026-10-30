use super::*;

fn pressure_controls() -> DensePressureSolve {
    DensePressureSolve {
        pressure_tolerance_pa: 1e-4,
        density_tolerance_kg_m3: 1e-8,
        maximum_iterations: 128,
    }
}
fn junction_controls() -> If97JunctionSolve {
    If97JunctionSolve {
        pressure: pressure_controls(),
        temperature_tolerance_k: 1e-8,
        maximum_iterations: 64,
    }
}
fn state(p: f64, t: f64) -> If97FluidState {
    If97FluidState::at_pressure_temperature(
        Pressure::new(p),
        Temperature::new(t),
        pressure_controls(),
        || false,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "actual={actual:.15e}, expected={expected:.15e}, tolerance={tolerance:.8e}"
    );
}

#[test]
fn regional_routing_preserves_original_states_and_density_work() {
    // G0: published-table regions and pressure-cycle states use their original
    // evaluators. Source wrapping introduces no caloric reference correction.
    for (p, t, expected) in [
        (3e6, 300.0, If97FluidRegion::Region1),
        (3_500.0, 700.0, If97FluidRegion::Region2),
        (20e6, 600.0, If97FluidRegion::Region1),
        (20e6, 630.0, If97FluidRegion::Region3),
        (20e6, 645.0, If97FluidRegion::Region3),
        (20e6, 823.15, If97FluidRegion::Region2),
        (25e6, 600.0, If97FluidRegion::Region1),
        (25e6, 650.0, If97FluidRegion::Region3),
        (25e6, 823.15, If97FluidRegion::Region2),
        (100e6, 900.0, If97FluidRegion::Region2),
    ] {
        let point = state(p, t);
        assert_eq!(point.region(), expected);
        assert_eq!(point.requested_pressure().value(), p);
        assert_eq!(point.temperature().value(), t);
        near(
            point.pressure_residual_pa(),
            point.pressure().value() - p,
            0.0,
        );
        near(
            point.pressure_residual_pa(),
            0.0,
            pressure_controls().pressure_tolerance_pa,
        );
        let explicit = If97FluidState::at_pressure_temperature_in_region(
            Pressure::new(p),
            Temperature::new(t),
            expected,
            pressure_controls(),
            || false,
        )
        .unwrap();
        assert_eq!(point, explicit);
        near(
            point.specific_enthalpy().value() - point.specific_internal_energy().value(),
            point.pressure().value() * point.specific_volume().value(),
            2e-9,
        );
        assert!(point.cp().value() >= point.cv().value());
        assert!(point.cv().value() > 0.0);
        assert!(point.sound_speed().value() > 0.0);
        match point {
            If97FluidState::Gibbs(source) => {
                assert_eq!(point.specific_enthalpy(), source.specific_enthalpy());
                assert_eq!(point.specific_entropy(), source.specific_entropy());
                assert_eq!(point.iterations(), 0);
                assert_eq!(point.source_evaluations(), 0);
            }
            If97FluidState::Dense(source) => {
                assert_eq!(point.specific_enthalpy(), source.state.specific_enthalpy());
                assert_eq!(point.specific_entropy(), source.state.specific_entropy());
                assert!(point.iterations() > 0);
                assert!(point.source_evaluations() > 0);
            }
        }
    }
    assert_eq!(state(25e6, 623.15).region(), If97FluidRegion::Region1);
    assert_eq!(
        state(25e6, 623.15_f64.next_up()).region(),
        If97FluidRegion::Region3
    );
    let boundary_p = boundary23_pressure_raw(700.0);
    assert_eq!(state(boundary_p, 700.0).region(), If97FluidRegion::Region2);
    assert_eq!(state(25e6, 650.0).phase(), If97FluidPhase::Supercritical);
    assert_eq!(
        state(25e6, 647.096).phase(),
        If97FluidPhase::CriticalIsotherm
    );
}

#[test]
fn saturation_routes_both_original_pairs_and_refuses_both_ambiguous_anchors() {
    // G0: no phase branch is inferred from Region 4 inverse roundoff.
    for t in [273.16, 400.0, 623.15, 630.0, 640.0, 647.0] {
        let pair = If97FluidSaturationPair::at_temperature(
            Temperature::new(t),
            pressure_controls(),
            || false,
        )
        .unwrap();
        assert_eq!(pair.temperature().value(), t);
        assert!(pair.liquid().density().value() > pair.vapor().density().value());
        assert!(pair.latent_enthalpy().value() > 0.0);
        assert!(pair.latent_entropy().value() > 0.0);
        near(
            pair.latent_enthalpy().value(),
            pair.vapor().specific_enthalpy().value() - pair.liquid().specific_enthalpy().value(),
            0.0,
        );
        near(
            pair.gibbs_energy_gap().value(),
            pair.vapor().specific_gibbs_energy().value()
                - pair.liquid().specific_gibbs_energy().value(),
            0.0,
        );
        assert_eq!(
            If97FluidState::at_pressure_temperature(
                pair.pressure(),
                pair.temperature(),
                pressure_controls(),
                || false
            ),
            Err(If97FluidError::AmbiguousSaturation),
        );
        let at_p =
            If97FluidSaturationPair::at_pressure(pair.pressure(), pressure_controls(), || false)
                .unwrap();
        assert_eq!(at_p.pressure(), pair.pressure());
        assert!(at_p.temperature().value() >= MINIMUM_STABLE_TEMPERATURE_K);
        if t == MINIMUM_STABLE_TEMPERATURE_K || t == REGION1_MAXIMUM_TEMPERATURE_K {
            assert_eq!(at_p.temperature(), pair.temperature());
            assert_eq!(at_p, pair);
        }
        assert_eq!(
            If97FluidState::at_pressure_temperature(
                at_p.pressure(),
                at_p.temperature(),
                pressure_controls(),
                || false
            ),
            Err(If97FluidError::AmbiguousSaturation),
        );
        let expected = if t <= 623.15 {
            If97FluidRegion::Region1
        } else {
            If97FluidRegion::Region3
        };
        assert_eq!(pair.liquid().region(), expected);
        assert!(pair.gibbs_energy_gap().value() != 0.0);
    }
    let join = saturation_pressure(Temperature::new(623.15)).unwrap();
    assert!(matches!(
        If97FluidSaturationPair::at_pressure(join, pressure_controls(), || false).unwrap(),
        If97FluidSaturationPair::Gibbs(_),
    ));
}

#[test]
fn regional_junctions_keep_signed_caloric_gaps_and_admitted_one_sided_sources() {
    // G0/G3: a source-fit junction has no latent inventory. The signed original
    // gaps can be positive OR negative; section 10/Table 43 permits both.
    for p in [20e6, 25e6, 50e6, 100e6] {
        for kind in [If97JunctionKind::LiquidDense, If97JunctionKind::DenseVapor] {
            let pair =
                If97RegionalJunction::new(Pressure::new(p), kind, junction_controls(), || false)
                    .unwrap();
            let [lo, hi] = pair.temperature_bounds_k();
            assert!(hi >= lo && hi - lo <= junction_controls().temperature_tolerance_k);
            assert_eq!(pair.pressure().value(), p);
            assert_eq!(pair.lower().requested_pressure().value(), p);
            assert_eq!(pair.upper().requested_pressure().value(), p);
            assert!(
                pair.pressure_residuals_pa()
                    .into_iter()
                    .all(|x| x.abs() <= 1e-4)
            );
            near(
                pair.enthalpy_change().value(),
                pair.upper().specific_enthalpy().value() - pair.lower().specific_enthalpy().value(),
                0.0,
            );
            near(
                pair.internal_energy_change().value() + p * pair.volume_change().value(),
                pair.enthalpy_change().value(),
                2e-6,
            );
            assert!(pair.enthalpy_change().value().abs() < 200.0);
            assert!(pair.entropy_change().value().abs() < 0.2);
            assert!(pair.gibbs_energy_change().value().abs() < 6.0);
            if kind == If97JunctionKind::LiquidDense {
                assert_eq!([lo, hi], [623.15; 2]);
                assert_eq!(pair.lower().region(), If97FluidRegion::Region1);
                assert_eq!(pair.upper().region(), If97FluidRegion::Region3);
                assert_eq!(pair.boundary_iterations(), 0);
                if p == 20e6 {
                    assert!(pair.enthalpy_change().value() > 5.0);
                }
                if p == 25e6 {
                    assert!(pair.enthalpy_change().value() < -6.0);
                }
            } else {
                assert_eq!(pair.lower().region(), If97FluidRegion::Region3);
                assert_eq!(pair.upper().region(), If97FluidRegion::Region2);
                assert!(boundary23_pressure_raw(lo) <= pair.lower().pressure().value());
                assert!(hi > 863.15 || p <= boundary23_pressure_raw(hi));
                if p == 20e6 {
                    assert!(pair.enthalpy_change().value() < -35.0);
                }
                if p == 25e6 {
                    assert!(pair.enthalpy_change().value() > 0.05);
                }
            }
        }
    }
}

#[test]
fn common_pressure_derivatives_follow_original_maxwell_relations() {
    // G1: finite-pressure perturbations are independent of the common derivative
    // accessors and exercise both Gibbs and density-inverted Helmholtz sources.
    for (p, t) in [(25e6, 600.0), (25e6, 650.0), (25e6, 823.15)] {
        let center = state(p, t);
        let step = 100.0;
        let below = state(p - step, t);
        let above = state(p + step, t);
        near(
            center.dh_dp_at_temperature().value(),
            (above.specific_enthalpy().value() - below.specific_enthalpy().value()) / (2.0 * step),
            2e-8,
        );
        near(
            center.ds_dp_at_temperature().value(),
            (above.specific_entropy().value() - below.specific_entropy().value()) / (2.0 * step),
            4e-11,
        );
        near(
            center.dv_dp_at_temperature().value(),
            (above.specific_volume().value() - below.specific_volume().value()) / (2.0 * step),
            2e-14,
        );
        near(
            center.density().value() * center.specific_volume().value(),
            1.0,
            2e-16,
        );
    }
}

#[test]
fn domain_critical_boundary_work_and_late_cancellation_refuse_without_a_state() {
    for (p, t) in [
        (1e5, 273.15),
        (100e6_f64.next_up(), 700.0),
        (1e5, 1073.15_f64.next_up()),
    ] {
        assert_eq!(
            If97FluidState::at_pressure_temperature(
                Pressure::new(p),
                Temperature::new(t),
                pressure_controls(),
                || false
            ),
            Err(If97FluidError::OutsideDomain),
        );
    }
    assert_eq!(
        If97FluidState::at_pressure_temperature(
            Pressure::new(CRITICAL_PRESSURE_PA),
            Temperature::new(CRITICAL_TEMPERATURE_K),
            pressure_controls(),
            || false
        ),
        Err(If97FluidError::CriticalSingularity),
    );
    assert_eq!(
        If97FluidSaturationPair::at_temperature(
            Temperature::new(CRITICAL_TEMPERATURE_K),
            pressure_controls(),
            || false
        ),
        Err(If97FluidError::OutsideDomain),
    );
    assert_eq!(
        If97RegionalJunction::new(
            Pressure::new(1e6),
            If97JunctionKind::DenseVapor,
            junction_controls(),
            || false
        ),
        Err(If97FluidError::OutsideDomain),
    );
    let mut short = junction_controls();
    short.maximum_iterations = 1;
    assert!(matches!(
        If97RegionalJunction::new(
            Pressure::new(25e6),
            If97JunctionKind::DenseVapor,
            short,
            || false
        ),
        Err(If97FluidError::BoundaryNotConverged { iterations: 1, .. }),
    ));
    assert_eq!(
        If97FluidState::at_pressure_temperature(
            Pressure::new(25e6),
            Temperature::new(650.0),
            pressure_controls(),
            || true
        ),
        Err(If97FluidError::Cancelled),
    );
    let mut polls = 0;
    If97RegionalJunction::new(
        Pressure::new(25e6),
        If97JunctionKind::DenseVapor,
        junction_controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut repeated = 0;
    assert_eq!(
        If97RegionalJunction::new(
            Pressure::new(25e6),
            If97JunctionKind::DenseVapor,
            junction_controls(),
            || {
                repeated += 1;
                repeated >= polls
            }
        ),
        Err(If97FluidError::Cancelled),
    );
    assert_eq!(polls, repeated);
}
