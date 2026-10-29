use super::*;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};

fn controls() -> WaterChartControls {
    WaterChartControls {
        minimum_temperature_k: 273.16,
        maximum_temperature_k: 1073.15,
        enthalpy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
        maximum_entropy_seam_j_kg_k: 0.1,
        heat_capacity: If97HeatCapacityBoundSolve {
            maximum_intervals: 65535,
        },
    }
}
fn chart(p: f64) -> IsobaricWater {
    IsobaricWater::new(Pressure::new(p), controls(), || false).unwrap()
}
fn at_t(law: &IsobaricWater, t: f64) -> WaterState {
    law.state_at_temperature(Temperature::new(t), || false)
        .unwrap()
}
fn at_h(law: &IsobaricWater, h: f64) -> WaterState {
    law.state_at_specific_enthalpy(h, || false)
        .unwrap()
        .into_state()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}
fn entropy_solve() -> WaterEntropySolve {
    WaterEntropySolve {
        entropy_tolerance_j_kg_k: 1e-9,
        maximum_iterations: 128,
    }
}

#[test]
fn source_h_and_raw_s_round_trip_across_liquid_coexistence_and_superheat() {
    let maximum_pressure =
        fs_thermochem::water::if97::saturation_pressure(Temperature::new(623.15))
            .unwrap()
            .value();
    for pressure in [101325.0, 3e6, maximum_pressure] {
        let law = chart(pressure);
        let ts = law.saturation().temperature().value();
        let mut states = [273.16, 0.5 * (273.16 + ts), ts - 0.25, ts + 25.0, 1073.15]
            .map(|t| at_t(&law, t))
            .to_vec();
        let l = law.saturation().liquid().specific_enthalpy().value();
        let v = law.saturation().vapor().specific_enthalpy().value();
        for quality in [0.0, 0.125, 0.5, 1.0] {
            let state = at_h(&law, l + quality * (v - l));
            near(state.vapor_quality().unwrap(), quality, 2e-14);
            near(
                state.specific_volume_m3_kg(),
                (1.0 - quality) * law.saturation().liquid().specific_volume().value()
                    + quality * law.saturation().vapor().specific_volume().value(),
                1e-13,
            );
            states.push(state);
        }
        for state in states {
            assert_eq!(state.law(), &law);
            assert_eq!(state.pressure_pa(), pressure);
            let h = law
                .state_at_specific_enthalpy(state.specific_enthalpy_j_kg(), || false)
                .unwrap();
            near(h.state().temperature_k(), state.temperature_k(), 1e-8);
            assert!(h.enthalpy_residual_j_kg().abs() <= controls().enthalpy_tolerance_j_kg);
            let s = law
                .state_at_specific_entropy(state.specific_entropy_j_kg_k(), entropy_solve(), || {
                    false
                })
                .unwrap();
            near(s.state().temperature_k(), state.temperature_k(), 1e-8);
            near(
                s.state().vapor_mass_fraction(),
                state.vapor_mass_fraction(),
                2e-12,
            );
            assert!(s.entropy_residual_j_kg_k().abs() <= entropy_solve().entropy_tolerance_j_kg_k);
            near(
                state.specific_enthalpy_j_kg() - state.specific_internal_energy_j_kg(),
                pressure * state.specific_volume_m3_kg(),
                2e-9,
            );
        }
        assert_eq!(at_t(&law, 300.0).phase(), WaterPhase::CompressedLiquid);
        assert_eq!(at_t(&law, 900.0).phase(), WaterPhase::SuperheatedVapor);
        assert_eq!(at_t(&law, 900.0).vapor_quality(), None);
    }
}

#[test]
fn complete_heating_condensation_and_reentry_keep_heat_work_and_exact_phase_knots() {
    let law = chart(101325.0);
    let mass = 0.25;
    let initial = at_t(&law, 300.0);
    let hot = at_t(&law, 700.0);
    let heat = mass * (hot.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg());
    let direct = law
        .advance_heat(initial, mass, heat, 1e-7, || false)
        .unwrap();
    near(direct.state.temperature_k(), 700.0, 1e-8);
    assert_eq!(direct.state.phase(), WaterPhase::SuperheatedVapor);
    assert!(direct.boundary_work_j > 0.05 * heat);
    near(
        direct.internal_energy_change_j + direct.boundary_work_j,
        heat,
        1e-7,
    );
    near(direct.enthalpy_change_j, heat, 1e-7);
    near(direct.represented_enthalpy_residual_j, 0.0, 1e-7);
    near(direct.represented_energy_residual_j, 0.0, 1e-7);
    let l = law.saturation().liquid().specific_enthalpy().value();
    let v = law.saturation().vapor().specific_enthalpy().value();
    let mut current = initial;
    let mut total = 0.0;
    for (h, phase) in [
        (l, WaterPhase::SaturatedLiquid),
        (v, WaterPhase::SaturatedVapor),
        (hot.specific_enthalpy_j_kg(), WaterPhase::SuperheatedVapor),
        (v, WaterPhase::SaturatedVapor),
        (l, WaterPhase::SaturatedLiquid),
        (
            initial.specific_enthalpy_j_kg(),
            WaterPhase::CompressedLiquid,
        ),
    ] {
        let q = mass * (h - current.specific_enthalpy_j_kg());
        let step = law.advance_heat(current, mass, q, 1e-7, || false).unwrap();
        assert_eq!(step.state.phase(), phase);
        if h == l {
            assert_eq!(step.phase_masses_kg, [mass, 0.0]);
        }
        if h == v {
            assert_eq!(step.phase_masses_kg, [0.0, mass]);
        }
        near(
            step.internal_energy_change_j + step.boundary_work_j,
            q,
            1e-7,
        );
        near(step.enthalpy_change_j, q, 1e-7);
        total += step.enthalpy_change_j;
        current = step.state;
    }
    near(current.temperature_k(), 300.0, 1e-8);
    near(total, 0.0, 2e-7);
    assert_eq!(
        law.advance_heat(initial, mass, 0.0, 0.0, || false)
            .unwrap()
            .state,
        initial
    );
    assert!(
        law.advance_heat(initial, mass, 1e-30, 1.0, || false)
            .is_err()
    );
    let small = law.advance_heat(hot, mass, 0.001, 1e-8, || false).unwrap();
    assert!(small.state.temperature_k() > hot.temperature_k());
    near(
        small.internal_energy_change_j + small.boundary_work_j,
        0.001,
        1e-8,
    );
}

#[test]
fn explicit_entropy_seam_retains_raw_if97_and_a_path_additive_thermal_primitive() {
    let law = chart(101325.0);
    let ts = law.saturation().temperature().value();
    let l = law.saturation().liquid();
    let v = law.saturation().vapor();
    let latent = v.specific_enthalpy().value() - l.specific_enthalpy().value();
    let raw_ds = v.specific_entropy().value() - l.specific_entropy().value();
    let defect = latent / ts - raw_ds;
    assert_ne!(defect, 0.0);
    assert_eq!(law.saturation_entropy_defect_j_kg_k(), defect);
    let h = |x: f64| l.specific_enthalpy().value() + x * latent;
    for x in [0.0, 0.25, 0.75, 1.0] {
        let state = at_h(&law, h(x));
        near(
            state.specific_entropy_j_kg_k(),
            l.specific_entropy().value() + x * raw_ds,
            2e-12,
        );
        near(state.thermal_entropy_correction_j_kg_k(), x * defect, 2e-14);
        near(
            state.thermal_specific_entropy_j_kg_k(),
            l.specific_entropy().value() + x * latent / ts,
            2e-12,
        );
    }
    near(
        ThermalLaw::specific_entropy_change_j_kg_k(&law, h(0.0), h(1.0))
            .unwrap()
            .unwrap(),
        latent / ts,
        1e-10,
    );
    let knots = [
        at_t(&law, 300.0).specific_enthalpy_j_kg(),
        h(0.0),
        h(0.4),
        h(1.0),
        at_t(&law, 700.0).specific_enthalpy_j_kg(),
    ];
    let whole = law
        .specific_entropy_change_j_kg_k(knots[0], knots[4], || false)
        .unwrap();
    let parts: f64 = knots
        .windows(2)
        .map(|p| {
            law.specific_entropy_change_j_kg_k(p[0], p[1], || false)
                .unwrap()
        })
        .sum();
    near(whole, parts, 1e-10);
    assert_eq!(
        law.specific_entropy_change_j_kg_k(knots[4], knots[0], || false)
            .unwrap(),
        -whole
    );
    let raw_vapor = law
        .state_at_specific_entropy(v.specific_entropy().value(), entropy_solve(), || false)
        .unwrap();
    assert_eq!(raw_vapor.state().phase(), WaterPhase::SaturatedVapor);
    assert_eq!(raw_vapor.entropy_residual_j_kg_k(), 0.0);
    let rejected = WaterChartControls {
        maximum_entropy_seam_j_kg_k: defect.abs() * 0.5,
        ..controls()
    };
    assert!(matches!(
        IsobaricWater::new(Pressure::new(101325.0), rejected, || false),
        Err(WaterError::EntropySeam { .. })
    ));
}

#[test]
fn bounds_domains_inverse_work_and_cancellation_are_explicit() {
    let law = chart(3e6);
    let ts = law.saturation().temperature().value();
    let slope = law.maximum_temperature_slope();
    assert!(slope.is_finite() && slope > 0.0);
    for bound in law.heat_capacity_bounds().iter().flatten() {
        assert!(bound.minimum_heat_capacity_j_kg_k() > 0.0);
        assert!(bound.interval_evaluations() <= controls().heat_capacity.maximum_intervals);
    }
    for (region, lo, hi) in [
        (If97Region::Region1, 273.16, ts),
        (If97Region::Region2, ts, 1073.15),
    ] {
        for i in 0..=16 {
            let t = lo + (hi - lo) * (i as f64 / 16.0);
            let source =
                If97State::at_pressure_temperature(Pressure::new(3e6), Temperature::new(t), region)
                    .unwrap();
            assert!(slope >= 1.0 / source.cp().value());
        }
    }
    let limited = WaterChartControls {
        heat_capacity: If97HeatCapacityBoundSolve {
            maximum_intervals: 1,
        },
        ..controls()
    };
    assert!(matches!(
        IsobaricWater::new(Pressure::new(16e6), limited, || false),
        Err(WaterError::HeatCapacity(
            If97HeatCapacityBoundError::NotProven { .. }
        ))
    ));
    for invalid in [
        WaterChartControls {
            minimum_temperature_k: 273.15,
            ..controls()
        },
        WaterChartControls {
            maximum_temperature_k: 1074.0,
            ..controls()
        },
    ] {
        assert!(IsobaricWater::new(Pressure::new(101325.0), invalid, || false).is_err());
    }
    assert!(IsobaricWater::new(Pressure::new(20e6), controls(), || false).is_err());
    assert_eq!(
        law.state_at_temperature(Temperature::new(ts), || false),
        Err(WaterError::AmbiguousSaturation)
    );
    let [lo, hi] = law.enthalpy_bounds_j_kg();
    assert!(law.state_at_specific_enthalpy(lo - 1.0, || false).is_err());
    assert!(law.state_at_specific_enthalpy(hi + 1.0, || false).is_err());
    let limited = IsobaricWater::new(
        Pressure::new(3e6),
        WaterChartControls {
            maximum_iterations: 1,
            enthalpy_tolerance_j_kg: 1e-12,
            ..controls()
        },
        || false,
    )
    .unwrap();
    let target = at_t(&law, 700.0).specific_enthalpy_j_kg();
    assert!(matches!(
        limited.state_at_specific_enthalpy(target, || false),
        Err(WaterError::NotConverged { .. })
    ));
    assert_eq!(
        law.state_at_specific_enthalpy(target, || true),
        Err(WaterError::Cancelled)
    );
    let current = at_t(&law, 300.0);
    let saved = current;
    let mut polls = 0;
    assert!(matches!(
        law.advance_heat(current, 0.5, 1000.0, 1e-6, || {
            polls += 1;
            polls == 3
        }),
        Err(WaterError::Cancelled)
    ));
    assert_eq!(current, saved);
    assert!(matches!(
        chart(101325.0).advance_heat(current, 0.5, 1.0, 1e-6, || false),
        Err(WaterError::ChartMismatch)
    ));
}

fn boundary(power: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: power,
    }
}
fn birth(id: u64, state: WaterState, mass: f64, power: f64) -> TimedEvent {
    TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(id),
            segment_id: id,
            body_id: id,
            material_id: MaterialId(1),
            mass_kg: mass,
            specific_h_j_kg: state.specific_enthalpy_j_kg(),
            geometry: CellGeometry {
                start_m: [id as f64, 0.0, 0.0],
                end_m: [id as f64, 0.0, 1.0],
                volume_m3: mass * state.specific_volume_m3_kg(),
                surface_area_m2: 1.0,
            },
            boundary: boundary(power),
        }),
    }
}
fn budget() -> DepositionBudget {
    DepositionBudget {
        max_cells: 2,
        max_links: 1,
        max_events: 4,
        max_substeps: 2048,
    }
}
fn config() -> StepConfig {
    StepConfig {
        max_dt_s: 1.0,
        energy_tolerance_j: 1e-6,
    }
}

#[test]
fn original_thermal_history_heats_liquid_boils_superheats_and_condenses() {
    let law = chart(101325.0);
    let cold = at_t(&law, 300.0);
    let hot = at_t(&law, 700.0);
    let mass = 0.05;
    let heat = mass * (hot.specific_enthalpy_j_kg() - cold.specific_enthalpy_j_kg());
    let power = heat / 100.0;
    let initial = DepositionSimulation::new(
        vec![
            birth(1, cold, mass, power),
            TimedEvent {
                time_s: 100.0,
                event: Event::Boundary(CellId(1), boundary(-power)),
            },
        ],
        &[(MaterialId(1), &law)],
        budget(),
    )
    .unwrap();
    let (active, _) = initial.advance_to(0.0, config(), || Ok(())).unwrap();
    let latent_middle = 0.5
        * (law.saturation().liquid().specific_enthalpy().value()
            + law.saturation().vapor().specific_enthalpy().value());
    let middle_time = mass * (latent_middle - cold.specific_enthalpy_j_kg()) / power;
    let (boiling, _, _) = active
        .advance_to_with_entropy(middle_time, config(), || Ok(()))
        .unwrap();
    let state = at_h(&law, boiling.cells()[&CellId(1)].specific_h_j_kg);
    assert_eq!(state.phase(), WaterPhase::LiquidVapor);
    near(state.vapor_mass_fraction(), 0.5, 1e-12);
    let (heated, audit, entropy) = active
        .advance_to_with_entropy(100.0, config(), || Ok(()))
        .unwrap();
    let final_state = at_h(&law, heated.cells()[&CellId(1)].specific_h_j_kg);
    assert_eq!(final_state.phase(), WaterPhase::SuperheatedVapor);
    near(final_state.temperature_k(), 700.0, 1e-7);
    near(audit.external_heat_j, heat, 1e-7);
    near(audit.residual_j, 0.0, 1e-6);
    near(
        entropy.thermal_storage_change_j_k,
        mass * (hot.thermal_specific_entropy_j_kg_k() - cold.thermal_specific_entropy_j_kg_k()),
        1e-8,
    );
    assert_eq!(
        heated.cells()[&CellId(1)].geometry.volume_m3,
        active.cells()[&CellId(1)].geometry.volume_m3
    );
    assert!(
        mass * final_state.specific_volume_m3_kg()
            > 100.0 * heated.cells()[&CellId(1)].geometry.volume_m3
    );
    let (cooled, audit, _) = heated
        .advance_to_with_entropy(200.0, config(), || Ok(()))
        .unwrap();
    near(audit.external_heat_j, -heat, 1e-7);
    let final_state = at_h(&law, cooled.cells()[&CellId(1)].specific_h_j_kg);
    assert_eq!(final_state.phase(), WaterPhase::CompressedLiquid);
    near(final_state.temperature_k(), 300.0, 1e-7);
    assert!(matches!(
        heated.advance_to(200.0, config(), || Err(DepositionError::Cancelled)),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(active.time_s(), 0.0);
    assert_eq!(heated.time_s(), 100.0);
}

#[test]
fn finite_steam_contact_drives_real_boiling_with_closed_enthalpy_and_pressure_work() {
    let law = chart(101325.0);
    let cold = at_t(&law, 300.0);
    let hot = at_t(&law, 700.0);
    let (cold_mass, hot_mass) = (0.05, 0.5);
    let initial = DepositionSimulation::new(
        vec![
            birth(1, cold, cold_mass, 0.0),
            birth(2, hot, hot_mass, 0.0),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 1,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.5,
                    conductance_w_k: 1.0,
                }),
            },
        ],
        &[(MaterialId(1), &law)],
        budget(),
    )
    .unwrap();
    let (active, _) = initial.advance_to(0.0, config(), || Ok(())).unwrap();
    let (next, audit, entropy) = active
        .advance_to_with_entropy(1000.0, config(), || Ok(()))
        .unwrap();
    assert_eq!(audit.external_heat_j, 0.0);
    near(audit.residual_j, 0.0, 1e-6);
    let cold_final = at_h(&law, next.cells()[&CellId(1)].specific_h_j_kg);
    let hot_final = at_h(&law, next.cells()[&CellId(2)].specific_h_j_kg);
    assert_eq!(cold_final.phase(), WaterPhase::SuperheatedVapor);
    assert!(cold_final.temperature_k() > law.saturation().temperature().value());
    assert!(
        cold_final.temperature_k() < hot_final.temperature_k() && hot_final.temperature_k() < 700.0
    );
    assert!(entropy.thermal_storage_change_j_k > 0.0 && entropy.conduction_generation_j_k > 0.0);
    let mut closure = 0.0;
    for (mass, a, b) in [(cold_mass, cold, cold_final), (hot_mass, hot, hot_final)] {
        let q = mass * (b.specific_enthalpy_j_kg() - a.specific_enthalpy_j_kg());
        let step = law.advance_heat(a, mass, q, 1e-6, || false).unwrap();
        closure += step.internal_energy_change_j + step.boundary_work_j;
        near(step.state.temperature_k(), b.temperature_k(), 1e-8);
    }
    near(closure, 0.0, 2e-6);
}
