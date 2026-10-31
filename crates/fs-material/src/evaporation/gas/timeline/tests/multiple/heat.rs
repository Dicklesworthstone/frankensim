//! Joint gas heat feeds original wall/phase evolution without renewing budgets.
use super::*;
use crate::evaporation::gas::chamber::heat::GasLiquidChamberHeatSolve;
use crate::evaporation::gas::timeline::multiple::{
    GasLiquidTimelineHeatSurface, advance_deposition_with_gas_and_heat as advance_heat,
};

fn heated<'a>(
    id: u64,
    binding: &'a GasLiquidBinding,
    coefficient: f64,
) -> GasLiquidTimelineHeatSurface<'a> {
    GasLiquidTimelineHeatSurface {
        surface: selected(id, binding),
        heat_transfer_coefficient_w_m2_k: coefficient,
        heat_controls: film_controls(),
    }
}
fn heat_only<'a>(
    id: u64,
    binding: &'a GasLiquidBinding,
    coefficient: f64,
) -> GasLiquidTimelineHeatSurface<'a> {
    let mut result = heated(id, binding, coefficient);
    result.surface.mass_transfer_coefficient_m_s = 0.0;
    result
}
fn solve() -> GasLiquidChamberHeatSolve {
    GasLiquidChamberHeatSolve {
        maximum_iterations: 128,
    }
}

#[test]
fn joint_heat_only_matches_independent_nonlinear_liquid_and_finite_gas_equations() {
    let source = catalog(0.0);
    let initial_gas = gas(&source, 1e-3);
    let binding = GasLiquidBinding {
        species: id(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let original = two_films(
        vec![
            film_birth(1, 0.0, 280.0, 0.0),
            film_birth(2, 0.0, 320.0, 0.0),
        ],
        100,
    )
    .advance_to(0.0, super::heat(0.2), || Ok(()))
    .unwrap()
    .0;
    let surfaces = [heat_only(1, &binding, 120.0), heat_only(2, &binding, 340.0)];
    let result = advance_heat(
        &original,
        &surfaces,
        &initial_gas,
        &source,
        0.2,
        super::heat(0.2),
        joint_controls(),
        solve(),
        || Ok(()),
    )
    .unwrap();
    // Independent algebra: the authored liquid h(T) is quadratic, while the
    // retained dry carrier has cv = (3.5 - 1) R/M. No production inverse is used.
    let capacity =
        initial_gas.mass_kg() * 2.5 * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
    let responses = |tg: f64| {
        let mut liquid = [0.0; 2];
        let mut residual = capacity * (tg - 300.0);
        for (i, (t0, coefficient)) in [(280.0, 120.0), (320.0, 340.0)].into_iter().enumerate() {
            let k = coefficient * 1e-3 * 0.2;
            let a = 0.01 * R * 0.0005;
            let b = 0.01 * R * 4.0 + k;
            let c = 0.01 * R * (4.0 * t0 + 0.0005 * t0 * t0) + k * tg;
            liquid[i] = 2.0 * c / (b + (b * b + 4.0 * a * c).sqrt());
            residual += k * (tg - liquid[i]);
        }
        (residual, liquid)
    };
    let (mut lo, mut hi) = (280.0, 320.0);
    for _ in 0..80 {
        let mid = lo + 0.5 * (hi - lo);
        if responses(mid).0 > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let expected = lo + 0.5 * (hi - lo);
    near(result.gas().state().state().temperature, expected, 2e-8);
    for (index, t) in responses(expected).1.into_iter().enumerate() {
        near(
            result
                .thermal()
                .temperature_k(CellId(index as u64 + 1))
                .unwrap(),
            t,
            2e-8,
        );
    }
    let audit = result.audit();
    assert!(audit.cells[&CellId(1)].sensible_heat_into_liquid_j > 0.0);
    assert!(audit.cells[&CellId(2)].sensible_heat_into_liquid_j < 0.0);
    assert_eq!(audit.transferred_mass_kg, 0.0);
    assert_eq!(audit.energy_into_gas_j, 0.0);
    assert_eq!(audit.phase_steps, 1);
    assert_eq!(result.thermal().mass_kg(), original.mass_kg());
    assert_eq!(result.gas().mass_kg(), initial_gas.mass_kg());
    near(
        audit.thermal_energy_change_j,
        audit.sensible_heat_into_liquid_j,
        1e-8,
    );
    near(
        audit.gas_energy_change_j,
        -audit.sensible_heat_into_liquid_j,
        1e-8,
    );
    assert!(audit.sum_absolute_energy_residual_j <= 1e-6);
    assert!(
        audit
            .last_heat
            .as_ref()
            .unwrap()
            .maximum_temperature_equation_residual_k
            <= 1e-8
    );
}

#[test]
fn wall_heat_sensible_heat_and_opposing_phase_transfer_share_one_storage_account() {
    let mut endpoints = Vec::new();
    for datum in [0.0, 1e8] {
        let source = catalog(datum);
        let initial_gas = ClosedGasParcel::new(
            source
                .with_composition(
                    Composition::new(CompositionBasis::MoleFraction, vec![0.985, 0.015]).unwrap(),
                )
                .unwrap(),
            point(300.0, 100000.0),
            1e-3,
        )
        .unwrap();
        let binding = GasLiquidBinding {
            species: id(),
            liquid_enthalpy_offset_j_kg: datum,
        };
        let original = two_films(
            vec![
                film_birth(1, 0.0, 280.0, 0.5),
                film_birth(2, 0.0, 320.0, 1.5),
            ],
            1000,
        )
        .advance_to(0.0, super::heat(0.025), || Ok(()))
        .unwrap()
        .0;
        let surfaces = [heated(2, &binding, 200.0), heated(1, &binding, 200.0)];
        let run = |inputs: &[GasLiquidTimelineHeatSurface<'_>]| {
            advance_heat(
                &original,
                inputs,
                &initial_gas,
                &source,
                0.5,
                super::heat(0.025),
                joint_controls(),
                solve(),
                || Ok(()),
            )
            .unwrap()
        };
        let result = run(&surfaces);
        let reordered = run(&[surfaces[1], surfaces[0]]);
        assert_eq!(result.gas(), reordered.gas());
        assert_eq!(result.thermal().cells(), reordered.thermal().cells());
        assert_eq!(result.audit(), reordered.audit());
        let audit = result.audit();
        assert!(audit.cells[&CellId(1)].condensed_mass_kg > 0.0);
        assert!(audit.cells[&CellId(2)].evaporated_mass_kg > 0.0);
        assert!(audit.cells[&CellId(1)].sensible_heat_into_liquid_j > 0.0);
        assert!(audit.cells[&CellId(2)].sensible_heat_into_liquid_j < 0.0);
        near(audit.thermal.external_heat_j, 1.0, 1e-12);
        near(
            audit.thermal_energy_change_j + audit.gas_energy_change_j,
            1.0,
            1e-6,
        );
        near(audit.mass_residual_kg, 0.0, 1e-14);
        assert!(audit.sum_absolute_energy_residual_j <= 1e-6);
        endpoints.push((
            result.gas().state().state().temperature,
            result.thermal().temperature_k(CellId(1)).unwrap(),
            result.thermal().temperature_k(CellId(2)).unwrap(),
            audit.transferred_mass_kg,
        ));
    }
    near(endpoints[0].0, endpoints[1].0, 2e-8);
    near(endpoints[0].1, endpoints[1].1, 2e-8);
    near(endpoints[0].2, endpoints[1].2, 2e-8);
    near(endpoints[0].3, endpoints[1].3, 1e-14);
}

#[test]
fn reciprocal_heat_respects_each_birth_and_preceding_exposure_at_covering() {
    let source = catalog(0.0);
    let initial_gas = ClosedGasParcel::new(
        gas(&source, 1e-3).model().clone(),
        point(340.0, 100000.0),
        1e-3,
    )
    .unwrap();
    let binding = GasLiquidBinding {
        species: id(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let mut covered = boundary(0.0);
    covered.reservoir = Some(Reservoir {
        area_m2: 1e-3,
        conductance_w_k: 0.0,
        temperature_k: 300.0,
    });
    let original = two_films(
        vec![
            film_birth(1, 0.1, 300.0, 0.0),
            film_birth(2, 0.2, 300.0, 0.0),
            TimedEvent {
                time_s: 0.3,
                event: Event::Boundary(CellId(1), covered),
            },
        ],
        1000,
    );
    let surfaces = [
        heat_only(1, &binding, 1000.0),
        heat_only(2, &binding, 1000.0),
    ];
    let run = |target| {
        advance_heat(
            &original,
            &surfaces,
            &initial_gas,
            &source,
            target,
            super::heat(1.0),
            joint_controls(),
            solve(),
            || Ok(()),
        )
        .unwrap()
    };
    let born = run(0.1);
    assert_eq!(born.gas(), &initial_gas);
    assert_eq!(born.audit().sensible_heat_into_liquid_j, 0.0);
    assert!(born.audit().last_heat.is_none());
    let at_cover = run(0.3);
    let result = run(0.4);
    assert_eq!(
        result.thermal().cells()[&CellId(1)].specific_h_j_kg,
        at_cover.thermal().cells()[&CellId(1)].specific_h_j_kg
    );
    assert_eq!(
        result.audit().cells[&CellId(1)].sensible_heat_into_liquid_j,
        at_cover.audit().cells[&CellId(1)].sensible_heat_into_liquid_j
    );
    assert!(
        result.audit().cells[&CellId(2)].sensible_heat_into_liquid_j
            > at_cover.audit().cells[&CellId(2)].sensible_heat_into_liquid_j
    );
    let audit = result.audit();
    assert_eq!(audit.thermal.events, 3);
    assert_eq!(audit.transferred_mass_kg, 0.0);
    near(
        audit.thermal_energy_change_j + audit.gas_energy_change_j - audit.incoming_energy_j,
        0.0,
        1e-6,
    );
    assert!(audit.sum_absolute_energy_residual_j <= 1e-6);
}

#[test]
fn complete_heat_retry_cancellation_and_zero_heat_use_the_same_original_clock() {
    let source = catalog(0.0);
    let initial_gas = ClosedGasParcel::new(
        gas(&source, 1e-3).model().clone(),
        point(340.0, 100000.0),
        1e-3,
    )
    .unwrap();
    let binding = GasLiquidBinding {
        species: id(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let events = vec![
        film_birth(1, 0.0, 300.0, 0.0),
        film_birth(2, 0.0, 300.0, 0.0),
    ];
    let original = two_films(events.clone(), 1000);
    let surfaces = [
        heat_only(1, &binding, 1000.0),
        heat_only(2, &binding, 1000.0),
    ];
    let mut polls = 0;
    let result = advance_heat(
        &original,
        &surfaces,
        &initial_gas,
        &source,
        0.5,
        super::heat(0.5),
        joint_controls(),
        solve(),
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    assert!(result.audit().rejected_intervals > 0);
    assert!(
        result.audit().attempted_thermal_substeps
            >= result.audit().thermal.substeps + result.audit().rejected_intervals
    );
    assert_eq!(result.audit().thermal.events, 2);
    assert_eq!(result.thermal().time_s(), 0.5);
    let mut seen = 0;
    assert!(matches!(
        advance_heat(
            &original,
            &surfaces,
            &initial_gas,
            &source,
            0.5,
            super::heat(0.5),
            joint_controls(),
            solve(),
            || {
                seen += 1;
                if seen == polls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }
        ),
        Err(GasEvaporationTimelineError::Thermal(
            DepositionError::Cancelled
        ))
    ));
    let limited = two_films(events, 1);
    assert!(matches!(
        advance_heat(
            &limited,
            &surfaces,
            &initial_gas,
            &source,
            0.5,
            super::heat(0.5),
            joint_controls(),
            solve(),
            || Ok(())
        ),
        Err(GasEvaporationTimelineError::Thermal(
            DepositionError::Budget(_)
        ))
    ));
    let mut disabled = surfaces;
    for input in &mut disabled {
        input.heat_transfer_coefficient_w_m2_k = 0.0;
        input.heat_controls.maximum_iterations = 0;
    }
    let zero = advance_heat(
        &original,
        &disabled,
        &initial_gas,
        &source,
        0.5,
        super::heat(0.5),
        joint_controls(),
        GasLiquidChamberHeatSolve {
            maximum_iterations: 0,
        },
        || Ok(()),
    )
    .unwrap();
    let phase = advance_multiple(
        &original,
        &[surfaces[0].surface, surfaces[1].surface],
        &initial_gas,
        &source,
        0.5,
        super::heat(0.5),
        joint_controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(zero.gas(), phase.gas());
    assert_eq!(zero.thermal().cells(), phase.thermal().cells());
    assert_eq!(zero.audit(), phase.audit());
    let mut mismatched = surfaces;
    mismatched[0].heat_controls.pure_liquid_enthalpy_range_j_kg[0] += 1.0;
    assert!(
        advance_heat(
            &original,
            &mismatched,
            &initial_gas,
            &source,
            0.0,
            super::heat(0.5),
            joint_controls(),
            solve(),
            || Ok(())
        )
        .is_err()
    );
    assert!(original.cells().is_empty());
    assert_eq!(original.time_s(), 0.0);
}

#[test]
fn complete_wall_joint_heat_and_phase_split_converges_under_time_refinement() {
    let source = catalog(0.0);
    let initial_gas = ClosedGasParcel::new(
        gas(&source, 1e-3).model().clone(),
        point(340.0, 100000.0),
        1e-3,
    )
    .unwrap();
    let binding = GasLiquidBinding {
        species: id(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let original = two_films(
        vec![
            film_birth(1, 0.0, 300.0, 1.0),
            film_birth(2, 0.0, 310.0, 2.0),
        ],
        10000,
    );
    let surfaces = [heated(1, &binding, 100.0), heated(2, &binding, 150.0)];
    let run = |dt| {
        let result = advance_heat(
            &original,
            &surfaces,
            &initial_gas,
            &source,
            1.0,
            super::heat(dt),
            joint_controls(),
            solve(),
            || Ok(()),
        )
        .unwrap();
        let audit = result.audit();
        assert!(audit.transferred_mass_kg > 0.0);
        assert!(audit.sensible_heat_into_liquid_j > 0.0);
        near(audit.thermal.external_heat_j, 3.0, 1e-12);
        near(
            audit.thermal_energy_change_j + audit.gas_energy_change_j - audit.incoming_energy_j,
            3.0,
            1e-6,
        );
        assert!(audit.sum_absolute_energy_residual_j <= 1e-6);
        result.gas().state().state().temperature
    };
    let reference = run(0.00625);
    let a = (run(0.2) - reference).abs();
    let b = (run(0.1) - reference).abs();
    let c = (run(0.05) - reference).abs();
    assert!(a > 1e-5);
    assert!(b < 0.75 * a, "{a}, {b}");
    assert!(c < 0.75 * b, "{b}, {c}");
}
