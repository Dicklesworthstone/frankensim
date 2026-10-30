use super::*;
use fs_thermochem::water::if97::regions::If97FluidRegion;

fn controls(lo: f64, hi: f64) -> RegionalWaterControls {
    RegionalWaterControls {
        minimum_temperature_k: lo,
        maximum_temperature_k: hi,
        enthalpy_tolerance_j_kg: 1e-6,
        maximum_iterations: 96,
        maximum_saturation_entropy_defect_j_kg_k: 0.1,
        pressure: DensePressureSolve {
            pressure_tolerance_pa: 1e-4,
            density_tolerance_kg_m3: 1e-8,
            maximum_iterations: 128,
        },
        junctions: RegionalJunctionControls {
            maximum_enthalpy_jump_j_kg: 200.0,
            maximum_entropy_jump_j_kg_k: 0.25,
            maximum_volume_jump_m3_kg: 1e-5,
            maximum_temperature_gap_k: 1e-8,
            maximum_iterations: 96,
        },
    }
}
fn chart(p: f64, lo: f64, hi: f64) -> RegionalWater {
    RegionalWater::new(Pressure::new(p), controls(lo, hi), || false).unwrap()
}
fn at_t(chart: &RegionalWater, t: f64) -> RegionalWaterState {
    chart
        .state_at_temperature(Temperature::new(t), || false)
        .unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!(
        a.is_finite() && (a - b).abs() <= tol,
        "{a:.15e} != {b:.15e}; delta={:.8e}; limit={tol:.8e}",
        a - b
    );
}
fn entropy_controls() -> WaterEntropySolve {
    WaterEntropySolve {
        entropy_tolerance_j_kg_k: 1e-8,
        maximum_iterations: 96,
    }
}

#[test]
fn cold_feed_to_steam_and_supercritical_heat_close_original_energy_across_both_junctions() {
    for p in [20e6, 25e6] {
        let chart = chart(p, 290.0, 873.15);
        let cold = at_t(&chart, 300.0);
        let hot = at_t(&chart, 823.15);
        assert_eq!(
            cold.single_phase_source().unwrap().region(),
            If97FluidRegion::Region1
        );
        assert_eq!(
            hot.single_phase_source().unwrap().region(),
            If97FluidRegion::Region2
        );
        let mass = 0.5;
        let heat = mass * (hot.specific_enthalpy_j_kg() - cold.specific_enthalpy_j_kg());
        let step = chart
            .advance_heat(cold, mass, heat, 2e-5, || false)
            .unwrap();
        near(step.state.temperature_k(), 823.15, 2e-8);
        near(step.enthalpy_change_j, heat, 2e-5);
        near(
            step.internal_energy_change_j + step.boundary_work_j,
            heat,
            2e-5,
        );
        near(
            step.sensible_enthalpy_change_j
                + step.latent_enthalpy_change_j
                + step.regional_enthalpy_matching_change_j,
            step.enthalpy_change_j,
            1e-9,
        );
        let matching: f64 = chart.junctions().map(|j| j.enthalpy_change().value()).sum();
        near(
            step.regional_enthalpy_matching_change_j,
            mass * matching,
            1e-10,
        );
        assert_eq!(chart.junctions().count(), 2);
        if p < CRITICAL_PRESSURE_PA {
            near(
                step.latent_enthalpy_change_j,
                mass * chart.saturation().unwrap().latent_enthalpy().value(),
                2e-6,
            );
            assert_eq!(hot.vapor_mass_fraction(), Some(1.0));
        } else {
            assert_eq!(step.latent_enthalpy_change_j, 0.0);
            assert_eq!(hot.phase(), RegionalWaterPhase::Supercritical);
            assert_eq!(hot.vapor_mass_fraction(), None);
        }
        let inverse = chart
            .state_at_specific_entropy(hot.specific_entropy_j_kg_k(), entropy_controls(), || false)
            .unwrap();
        near(inverse.state().temperature_k(), 823.15, 2e-8);
        let back = chart
            .advance_heat(step.state, mass, -heat, 2e-5, || false)
            .unwrap();
        near(back.state.temperature_k(), 300.0, 2e-8);
        near(
            step.regional_enthalpy_matching_change_j + back.regional_enthalpy_matching_change_j,
            0.0,
            1e-12,
        );
        near(step.boundary_work_j + back.boundary_work_j, 0.0, 1e-6);
    }
}

#[test]
fn caloric_gaps_and_overlaps_are_explicit_and_heat_continues_in_its_retained_direction() {
    // The h and s jump signs can differ at the same boundary. Every coordinate
    // must be classified using its actual original one-sided endpoint values.
    for p in [20e6, 25e6] {
        let chart = chart(p, 300.0, 873.15);
        for junction in chart.junctions() {
            let (a, b) = (junction.lower(), junction.upper());
            for coordinate in [WaterCoordinate::Enthalpy, WaterCoordinate::Entropy] {
                let read = |s: If97FluidState| {
                    if coordinate == WaterCoordinate::Enthalpy {
                        s.specific_enthalpy().value()
                    } else {
                        s.specific_entropy().value()
                    }
                };
                let target = read(a) + 0.5 * (read(b) - read(a));
                let result = if coordinate == WaterCoordinate::Enthalpy {
                    chart
                        .state_at_specific_enthalpy(target, || false)
                        .map(|s| s.into_state())
                } else {
                    chart
                        .state_at_specific_entropy(target, entropy_controls(), || false)
                        .map(|s| s.into_state())
                };
                if read(b) > read(a) {
                    assert!(
                        matches!(result, Err(RegionalWaterError::CaloricGap { .. })),
                        "{p} {coordinate:?}: {result:?}"
                    );
                } else {
                    assert!(
                        matches!(
                            result,
                            Err(RegionalWaterError::AmbiguousCaloricState { candidates: 2, .. })
                        ),
                        "{p} {coordinate:?}: {result:?}"
                    );
                }
            }
            if junction.enthalpy_change().value() < 0.0 {
                let target = 0.5 * (a.specific_enthalpy().value() + b.specific_enthalpy().value());
                let lower = at_t(&chart, a.temperature().value() - 1.0);
                let upper = at_t(&chart, b.temperature().value() + 1.0);
                let heating = chart
                    .advance_heat(
                        lower,
                        1.0,
                        target - lower.specific_enthalpy_j_kg(),
                        2e-5,
                        || false,
                    )
                    .unwrap();
                let cooling = chart
                    .advance_heat(
                        upper,
                        1.0,
                        target - upper.specific_enthalpy_j_kg(),
                        2e-5,
                        || false,
                    )
                    .unwrap();
                assert_eq!(heating.inversion_candidates, 2);
                assert_eq!(cooling.inversion_candidates, 2);
                assert!(heating.state.temperature_k() > lower.temperature_k());
                assert!(heating.state.temperature_k() < a.temperature().value());
                assert!(cooling.state.temperature_k() < upper.temperature_k());
                assert!(cooling.state.temperature_k() > b.temperature().value());
                assert_eq!(heating.regional_enthalpy_matching_change_j, 0.0);
                assert_eq!(cooling.regional_enthalpy_matching_change_j, 0.0);
            }
        }
    }
}

#[test]
fn independent_piece_cp_and_cp_over_t_integrals_separate_sensible_latent_and_matching() {
    // G1: numerical quadrature uses source derivatives and never differences
    // source calorics. Matching defects and physical latent changes are separate.
    for p in [20e6, 25e6] {
        let chart = chart(p, 300.0, 823.15);
        let mut h = 0.0;
        let mut s = 0.0;
        let panels = 1024;
        for piece in chart.pieces.into_iter().flatten() {
            let lo = piece.lower.temperature().value();
            let hi = piece.upper.temperature().value();
            if lo == hi {
                continue;
            }
            let width = (hi - lo) / f64::from(panels);
            for i in 0..=panels {
                let source = if i == 0 {
                    piece.lower
                } else if i == panels {
                    piece.upper
                } else {
                    If97FluidState::at_pressure_temperature_in_region(
                        chart.pressure,
                        Temperature::new(lo + f64::from(i) * width),
                        piece.lower.region(),
                        chart.controls.pressure,
                        || false,
                    )
                    .unwrap()
                };
                let weight = if i == 0 || i == panels {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                h += weight * source.cp().value() * width / 3.0;
                s += weight * source.cp().value() / source.temperature().value() * width / 3.0;
            }
        }
        let cold = at_t(&chart, 300.0);
        let hot = at_t(&chart, 823.15);
        let heat = hot.specific_enthalpy_j_kg() - cold.specific_enthalpy_j_kg();
        let step = chart.advance_heat(cold, 1.0, heat, 2e-5, || false).unwrap();
        near(step.sensible_enthalpy_change_j, h, 0.03);
        near(step.sensible_entropy_change_j_k, s, 5e-5);
        near(
            step.source_entropy_change_j_k,
            step.sensible_entropy_change_j_k
                + step.latent_entropy_change_j_k
                + step.regional_entropy_matching_change_j_k,
            1e-10,
        );
        near(
            step.source_entropy_change_j_k,
            hot.specific_entropy_j_kg_k() - cold.specific_entropy_j_kg_k(),
            1e-7,
        );
    }
}

#[test]
fn unused_junctions_do_not_block_cold_liquid_dense_boiling_or_very_low_pressure_vapor() {
    let mut c = controls(300.0, 400.0);
    c.junctions.maximum_enthalpy_jump_j_kg = 0.0;
    c.junctions.maximum_entropy_jump_j_kg_k = 0.0;
    c.junctions.maximum_volume_jump_m3_kg = 0.0;
    c.maximum_saturation_entropy_defect_j_kg_k = 0.0;
    let cold = RegionalWater::new(Pressure::new(25e6), c, || false).unwrap();
    assert_eq!(cold.junctions().count(), 0);
    assert!(cold.saturation().is_none());
    assert_eq!(
        at_t(&cold, 350.0).single_phase_source().unwrap().region(),
        If97FluidRegion::Region1
    );
    c.minimum_temperature_k = 273.16;
    let vapor = RegionalWater::new(Pressure::new(100.0), c, || false).unwrap();
    assert!(vapor.saturation().is_none());
    assert_eq!(
        at_t(&vapor, 300.0).single_phase_source().unwrap().region(),
        If97FluidRegion::Region2
    );
    c.minimum_temperature_k = 625.0;
    c.maximum_temperature_k = 645.0;
    c.maximum_saturation_entropy_defect_j_kg_k = 0.1;
    let dense = RegionalWater::new(Pressure::new(20e6), c, || false).unwrap();
    assert_eq!(dense.junctions().count(), 0);
    assert!(dense.saturation().is_some());
    let mixture = dense.state_at_quality(0.4, || false).unwrap();
    assert!(mixture.vapor_volume_fraction().unwrap() > 0.4);
    let mass = 0.2;
    let liquid = dense.state_at_quality(0.0, || false).unwrap();
    let heat = mass * (mixture.specific_enthalpy_j_kg() - liquid.specific_enthalpy_j_kg());
    let step = dense
        .advance_heat(liquid, mass, heat, 2e-5, || false)
        .unwrap();
    near(step.phase_masses_kg.unwrap()[1], mass * 0.4, 1e-15);
    assert_eq!(step.regional_enthalpy_matching_change_j, 0.0);
    near(step.latent_enthalpy_change_j, heat, 2e-5);

    for temperature in [273.16, 623.15, 630.0] {
        let p = saturation_pressure(Temperature::new(temperature)).unwrap();
        let endpoint =
            RegionalWater::new(p, controls(temperature, temperature + 1.0), || false).unwrap();
        assert_eq!(
            endpoint.saturation().unwrap().temperature().value(),
            temperature
        );
        assert_eq!(endpoint.saturation().unwrap().pressure(), p);
        let plateau = RegionalWater::new(p, controls(temperature, temperature), || false).unwrap();
        let liquid = plateau.state_at_quality(0.0, || false).unwrap();
        let heat = plateau.saturation().unwrap().latent_enthalpy().value();
        let vapor = plateau
            .advance_heat(liquid, 1.0, heat, 2e-5, || false)
            .unwrap();
        assert_eq!(vapor.state.vapor_quality(), Some(1.0));
    }
    let mut invalid = controls(300.0, 400.0);
    invalid.pressure.pressure_tolerance_pa = f64::NAN;
    assert!(matches!(
        RegionalWater::new(Pressure::new(25e6), invalid, || false),
        Err(RegionalWaterError::Invalid(_))
    ));
}

#[test]
fn heating_subdivision_and_tiny_phase_changes_preserve_the_original_accounts() {
    let chart = chart(20e6, 290.0, 873.15);
    let cold = at_t(&chart, 300.0);
    let hot = at_t(&chart, 823.15);
    let heat = hot.specific_enthalpy_j_kg() - cold.specific_enthalpy_j_kg();
    let direct = chart.advance_heat(cold, 1.0, heat, 2e-5, || false).unwrap();
    let mut state = cold;
    let mut matching = 0.0;
    let mut latent = 0.0;
    let mut work = 0.0;
    for _ in 0..20 {
        let step = chart
            .advance_heat(state, 1.0, heat / 20.0, 2e-5, || false)
            .unwrap();
        state = step.state;
        matching += step.regional_enthalpy_matching_change_j;
        latent += step.latent_enthalpy_change_j;
        work += step.boundary_work_j;
    }
    near(state.temperature_k(), hot.temperature_k(), 2e-7);
    near(matching, direct.regional_enthalpy_matching_change_j, 1e-10);
    near(latent, direct.latent_enthalpy_change_j, 2e-6);
    near(work, direct.boundary_work_j, 2e-5);
    let tiny = chart
        .advance_heat(cold, 0.5, 0.001, 2e-6, || false)
        .unwrap();
    assert!(tiny.state.temperature_k() > cold.temperature_k());
    let zero = chart.advance_heat(cold, 0.5, 0.0, 0.0, || false).unwrap();
    assert_eq!(zero.state, cold);
    assert_eq!(zero.enthalpy_change_j, 0.0);
    // A tiny but representable latent increment must not disappear when its
    // quality is compared with an integer path-segment number.
    let low = self::chart(10_000.0, 273.16, 500.0);
    let liquid = low.state_at_quality(0.0, || false).unwrap();
    let q = liquid.specific_enthalpy_j_kg().next_up() - liquid.specific_enthalpy_j_kg();
    let tiny = low.advance_heat(liquid, 1.0, q, q, || false).unwrap();
    assert!(tiny.state.vapor_quality().unwrap() > 0.0);
    assert_eq!(tiny.state.temperature_k(), liquid.temperature_k());
}

#[test]
fn matching_work_domain_energy_and_cancellation_budgets_refuse_atomically() {
    let chart = chart(20e6, 290.0, 873.15);
    let state = at_t(&chart, 300.0);
    let mut c = controls(290.0, 873.15);
    c.junctions.maximum_enthalpy_jump_j_kg = 0.0;
    assert!(matches!(
        RegionalWater::new(Pressure::new(20e6), c, || false),
        Err(RegionalWaterError::JunctionBudget { .. })
    ));
    assert_eq!(
        chart.state_at_quality(0.5, || true),
        Err(RegionalWaterError::Cancelled)
    );
    assert_eq!(
        chart.state_at_temperature(chart.saturation().unwrap().temperature(), || false),
        Err(RegionalWaterError::AmbiguousSaturation)
    );
    assert!(matches!(
        chart.advance_heat(state, 0.0, 1.0, 1.0, || false),
        Err(RegionalWaterError::Invalid(_))
    ));
    assert!(matches!(
        chart.advance_heat(state, 1.0, f64::MIN_POSITIVE, 1.0, || false),
        Err(RegionalWaterError::Unrepresentable(_))
    ));
    let other = self::chart(25e6, 290.0, 873.15);
    assert_eq!(
        other.advance_heat(state, 1.0, 1.0, 1.0, || false),
        Err(RegionalWaterError::ChartMismatch)
    );
    let mut c = controls(290.0, 873.15);
    c.maximum_iterations = 1;
    let short = RegionalWater::new(Pressure::new(20e6), c, || false).unwrap();
    assert!(matches!(
        short.state_at_specific_enthalpy(at_t(&chart, 750.0).specific_enthalpy_j_kg(), || false),
        Err(RegionalWaterError::NotConverged { .. })
    ));
    let mut polls = 0;
    let result = chart.advance_heat(state, 1.0, 10_000.0, 1e-5, || {
        polls += 1;
        polls >= 8
    });
    assert_eq!(result, Err(RegionalWaterError::Cancelled));
    assert_eq!(state, at_t(&chart, 300.0));
    assert_eq!(
        RegionalWater::new(
            Pressure::new(CRITICAL_PRESSURE_PA),
            controls(300.0, 800.0),
            || false
        ),
        Err(RegionalWaterError::OutsideDomain)
    );
}
