//! Native births and topology events around the original joint gas/heat owner.
use super::*;
use crate::evaporation::gas::implicit::timeline::advance_deposition_with_gas_timeline;
use fs_thermal::enthalpy::Reservoir;
use fs_thermal::evaporation::{EvaporationError, EvaporationLaw, EvaporationProperties};

// Binary-exact areas make the final zero exposure a physical boundary case,
// without asking the native geometry owner to round a negative remainder.
const AREA: f64 = 1.0 / 1024.0;

fn scheduled<'m>(
    events: Vec<TimedEvent>,
    law: &'m dyn ThermalLaw,
    max_substeps: usize,
) -> DepositionSimulation<'m> {
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 4,
            max_links: 2,
            max_events: 16,
            max_substeps,
        },
    )
    .unwrap()
}

fn event_policy() -> AdaptiveImplicitGasControls {
    let mut config = policy(1.0, 1e-5);
    config.initial_dt_s = 0.5;
    config.maximum_dt_s = 0.5;
    config
}

fn at_events<'m>(thermal: &DepositionSimulation<'m>, time_s: f64) -> DepositionSimulation<'m> {
    thermal
        .advance_to(
            time_s,
            StepConfig {
                max_dt_s: 0.5,
                energy_tolerance_j: 1e-6,
            },
            || Ok(()),
        )
        .unwrap()
        .0
}

fn birth(time_s: f64, id: u64, temperature: f64, source_w: f64, shift: f64) -> TimedEvent {
    let mut deposited = cell(id, temperature, source_w, AREA);
    deposited.specific_h_j_kg += shift;
    TimedEvent {
        time_s,
        event: Event::Birth(deposited),
    }
}

fn history(shift: f64) -> Vec<TimedEvent> {
    let mut covered = cell(2, 280.0, 2.0, AREA).boundary;
    covered.reservoir = Some(Reservoir {
        area_m2: 0.75 * AREA,
        conductance_w_k: 0.0,
        temperature_k: 300.0,
    });
    vec![
        birth(0.5, 1, 350.0, 1.0, shift),
        birth(1.0, 2, 280.0, 2.0, shift),
        TimedEvent {
            time_s: 1.0,
            event: Event::Contact(ThermalLink {
                id: 9,
                a: CellId(1),
                b: CellId(2),
                area_m2: 0.25 * AREA,
                conductance_w_k: 0.1,
            }),
        },
        TimedEvent {
            time_s: 1.5,
            event: Event::Boundary(CellId(2), covered),
        },
        birth(2.0, 3, 300.0, 0.0, shift),
    ]
}

fn declarations(binding: &GasLiquidBinding) -> [ImplicitGasSurface<'_>; 3] {
    // Authored descriptor order need not match birth order or canonical cell IDs.
    [3, 2, 1].map(|id| surface(id, binding, 0.02, 20.0))
}

struct PhaseQueries(Cell<usize>);
impl EvaporationLaw for PhaseQueries {
    fn temperature_range_k(&self) -> [f64; 2] {
        self.0.set(self.0.get() + 1);
        Phase.temperature_range_k()
    }
    fn properties(&self, temperature: f64) -> Result<EvaporationProperties, EvaporationError> {
        self.0.set(self.0.get() + 1);
        Phase.properties(temperature)
    }
}

#[test]
fn empty_geometry_and_heat_only_preludes_activate_liquid_only_at_its_actual_birth() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let initial_gas = gas(&source, 300.0, 0.0, 0.001);
    let spy = PhaseQueries(Cell::new(0));
    let mut declared = surface(1, &binding, 0.02, 20.0);
    declared.law = &spy;
    let original = scheduled(vec![birth(1.0, 1, 300.0, 2.0, 0.0)], &ORIGINAL_LIQUID, 64);
    let calls = Cell::new(0usize);
    for target in [0.5, 1.0] {
        let candidate = advance_deposition_with_gas_timeline(
            &original,
            &[declared],
            &initial_gas,
            &source,
            target,
            event_policy(),
            |_, _, _| {
                calls.set(calls.get() + 1);
                Ok(coefficients(()))
            },
            || Ok(()),
        )
        .unwrap();
        assert_eq!(candidate.thermal().time_s(), target);
        assert_eq!(candidate.gas(), &initial_gas);
        assert_eq!(candidate.report().accepted_intervals, 0);
        assert_eq!(candidate.report().comparisons, 0);
        assert_eq!(candidate.report().gas_iterations, 0);
        assert!(candidate.report().surface_totals.is_empty());
        assert!(candidate.report().last_step.is_none());
        assert!(candidate.payload().is_none());
        assert_eq!(calls.get(), 0);
        if target < 1.0 {
            assert!(candidate.thermal().cells().is_empty());
            assert_eq!(spy.0.get(), 0);
            assert_eq!(candidate.report().events, 0);
            assert_eq!(candidate.report().incoming_mass_kg, 0.0);
        } else {
            assert!(
                spy.0.get() > 0,
                "target birth still admits its phase source"
            );
            assert_eq!(candidate.thermal().cells().len(), 1);
            assert_eq!(candidate.report().events, 1);
            assert_eq!(candidate.report().incoming_mass_kg, LIQUID_MASS);
            assert_eq!(
                candidate.report().incoming_energy_j,
                LIQUID_MASS * liquid(300.0).specific_enthalpy_j_kg
            );
        }
    }
    let candidate = advance_deposition_with_gas_timeline(
        &original,
        &[declared],
        &initial_gas,
        &source,
        1.5,
        event_policy(),
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    let born = at_events(&original, 1.0);
    let manual = advance_deposition_with_gas_adaptive(
        &born,
        &[declared],
        &initial_gas,
        &source,
        1.5,
        event_policy(),
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    assert_eq!(candidate.thermal().cells(), manual.thermal().cells());
    assert_eq!(candidate.gas(), manual.gas());
    assert!(candidate.gas().species_mass_kg(&vapor()).unwrap() > 0.0);
    near(candidate.report().external_heat_j, 1.0, 1e-10);
    near(candidate.report().mass_residual_kg, 0.0, 1e-14);

    // An unselected cell has real heat evolution before the selected birth;
    // this interval must use the original heat network and leave gas frozen.
    let heat_first = scheduled(
        vec![
            birth(0.0, 9, 300.0, 3.0, 0.0),
            birth(1.0, 1, 300.0, 2.0, 0.0),
        ],
        &ORIGINAL_LIQUID,
        64,
    );
    let heated = advance_deposition_with_gas_timeline(
        &heat_first,
        &[declared],
        &initial_gas,
        &source,
        1.0,
        event_policy(),
        |thermal, liquids, gas| {
            assert!(liquids.is_empty());
            assert_eq!(gas, &initial_gas);
            assert_eq!(thermal.cells().len(), 1);
            Ok(coefficients(thermal.time_s()))
        },
        || Ok(()),
    )
    .unwrap();
    assert_eq!(heated.gas(), &initial_gas);
    assert_eq!(heated.report().gas_iterations, 0);
    assert_eq!(heated.report().events, 2);
    assert_eq!(heated.payload(), Some(&1.0));
    assert!(heated.report().surface_totals.is_empty());
    near(
        heated.thermal().cells()[&CellId(9)].specific_h_j_kg,
        liquid(300.0).specific_enthalpy_j_kg + 3.0 / LIQUID_MASS,
        1e-8,
    );
    assert_eq!(heated.thermal().cells()[&CellId(1)].mass_kg, LIQUID_MASS);
    near(heated.report().external_heat_j, 3.0, 1e-10);
}

#[test]
fn births_contact_and_covering_match_manual_event_split_composition_and_preserve_final_payload() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let surfaces = declarations(&binding);
    let original = scheduled(history(0.0), &ORIGINAL_LIQUID, 64);
    let initial_gas = gas(&source, 320.0, 1600.0, 0.001);
    let candidate = advance_deposition_with_gas_timeline(
        &original,
        &surfaces,
        &initial_gas,
        &source,
        2.0,
        event_policy(),
        |thermal, liquids, _| {
            let ids: Vec<_> = liquids.keys().copied().collect();
            assert_eq!(
                ids,
                if thermal.time_s() <= 1.0 {
                    vec![CellId(1)]
                } else {
                    vec![CellId(1), CellId(2)]
                }
            );
            Ok(coefficients((thermal.time_s(), ids, thermal.cells().len())))
        },
        || Ok(()),
    )
    .unwrap();

    let mut manual_thermal = at_events(&original, 0.5);
    let mut manual_gas = initial_gas.clone();
    let mut manual_transfer = BTreeMap::<CellId, f64>::new();
    let mut intervals = 0usize;
    for target in [1.0, 1.5, 2.0] {
        let active: Vec<_> = surfaces
            .iter()
            .copied()
            .filter(|surface| manual_thermal.cells().contains_key(&surface.cell_id))
            .collect();
        let part = advance_deposition_with_gas_adaptive(
            &manual_thermal,
            &active,
            &manual_gas,
            &source,
            target,
            event_policy(),
            empty_heat,
            || Ok(()),
        )
        .unwrap();
        for (id, total) in &part.report().surface_totals {
            *manual_transfer.entry(*id).or_default() += total.transferred_mass_kg;
        }
        intervals += part.report().accepted_intervals;
        let (thermal, gas, _, _, _) = part.into_parts();
        manual_thermal = at_events(&thermal, target);
        manual_gas = gas;
    }
    assert_eq!(candidate.thermal().cells(), manual_thermal.cells());
    assert_eq!(candidate.thermal().links(), manual_thermal.links());
    assert_eq!(candidate.gas(), &manual_gas);
    assert!(candidate.thermal().pending_events().is_empty());
    assert_eq!(candidate.thermal().exposed_area_m2(CellId(2)).unwrap(), 0.0);
    assert_eq!(
        candidate.payload(),
        Some(&(2.0, vec![CellId(1), CellId(2)], 2))
    );
    let report = candidate.report();
    assert_eq!(report.events, 5);
    assert_eq!(report.accepted_intervals, intervals);
    near(report.incoming_mass_kg, 3.0 * LIQUID_MASS, 1e-16);
    near(report.external_heat_j, 3.5, 1e-10);
    near(report.mass_residual_kg, 0.0, 1e-14);
    near(report.energy_residual_j, 0.0, 1e-7);
    assert!(report.surface_totals[&CellId(1)].transferred_mass_kg > 0.0);
    assert!(report.surface_totals[&CellId(2)].transferred_mass_kg < 0.0);
    for (id, mass) in manual_transfer {
        near(report.surface_totals[&id].transferred_mass_kg, mass, 1e-15);
        near(
            LIQUID_MASS - candidate.thermal().cells()[&id].mass_kg,
            mass,
            1e-14,
        );
    }
    let final_surfaces = &report.last_step.as_ref().unwrap().thermal.thermal.surfaces;
    assert_eq!(final_surfaces[&CellId(2)].transferred_mass_kg, 0.0);
    assert!(!final_surfaces.contains_key(&CellId(3)));
    assert!(!report.surface_totals.contains_key(&CellId(3)));
    assert_eq!(candidate.thermal().cells()[&CellId(3)].mass_kg, LIQUID_MASS);
    assert_eq!(
        candidate.thermal().cells()[&CellId(3)].specific_h_j_kg,
        liquid(300.0).specific_enthalpy_j_kg
    );
}

struct ShiftedLiquid(f64);
impl ThermalLaw for ShiftedLiquid {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        ORIGINAL_LIQUID.temperature_k(h - self.0)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        ORIGINAL_LIQUID.maximum_temperature_slope()
    }
}

#[test]
fn birth_energy_keeps_native_and_liquid_to_gas_datums_separate_with_invariant_physics() {
    let mut baseline: Option<(Vec<f64>, Vec<f64>, f64, f64)> = None;
    for (gas_shift, liquid_shift) in [(0.0, 0.0), (1e12, 1_048_576.0)] {
        let source = catalog(gas_shift);
        let binding = binding(gas_shift - liquid_shift);
        let law = ShiftedLiquid(liquid_shift);
        let original = scheduled(history(liquid_shift), &law, 64);
        let initial_gas = gas(&source, 320.0, 1600.0, 0.001);
        let mut surfaces = declarations(&binding);
        for surface in &mut surfaces {
            for h in &mut surface.pure_liquid_enthalpy_range_j_kg {
                *h += liquid_shift;
            }
        }
        let mut config = event_policy();
        config.solve.phase.latent_compatibility_tolerance_j_kg = 0.01;
        let candidate = advance_deposition_with_gas_timeline(
            &original,
            &surfaces,
            &initial_gas,
            &source,
            2.0,
            config,
            empty_heat,
            || Ok(()),
        )
        .unwrap();
        let report = candidate.report();
        let birth_energy: f64 = [350.0, 280.0, 300.0]
            .into_iter()
            .map(|t| LIQUID_MASS * (liquid(t).specific_enthalpy_j_kg + liquid_shift))
            .sum();
        near(report.incoming_energy_j, birth_energy, 1e-9);
        near(
            report.liquid_datum_incoming_energy_j,
            3.0 * LIQUID_MASS * (gas_shift - liquid_shift),
            1e-5,
        );
        near(
            report.incoming_energy_j + report.liquid_datum_incoming_energy_j,
            [350.0, 280.0, 300.0]
                .into_iter()
                .map(|t| LIQUID_MASS * (liquid(t).specific_enthalpy_j_kg + gas_shift))
                .sum(),
            1e-5,
        );

        let temperatures: Vec<_> = [CellId(1), CellId(2), CellId(3)]
            .into_iter()
            .map(|id| candidate.thermal().temperature_k(id).unwrap())
            .collect();
        let masses: Vec<_> = candidate
            .thermal()
            .cells()
            .values()
            .map(|cell| cell.mass_kg)
            .collect();
        let gas_t = candidate.gas().state().state().temperature;
        let vapor_mass = candidate.gas().species_mass_kg(&vapor()).unwrap();
        let vapor_initial = initial_gas.species_mass_kg(&vapor()).unwrap();
        let thermal_change: f64 = temperatures
            .iter()
            .zip(&masses)
            .zip([350.0, 280.0, 300.0])
            .map(|((&t, &m), birth_t)| {
                m * liquid(t).specific_enthalpy_j_kg
                    - LIQUID_MASS * liquid(birth_t).specific_enthalpy_j_kg
            })
            .sum();
        // Integrate the original NASA9 Cv without subtracting huge formation
        // totals; actual birth masses and the source-relative vapor inventory
        // independently determine the whole-window physical balance.
        let gas_change =
            initial_gas.species_mass_kg(&carrier()).unwrap() * CARRIER_R * 2.5 * (gas_t - 320.0)
                + vapor_mass * R * (gas_t - 320.0) * (3.0 + 0.0005 * (gas_t + 320.0))
                + (vapor_mass - vapor_initial) * vapor_u(320.0);
        near(thermal_change + gas_change, 3.5, 1e-6);
        near(report.mass_residual_kg, 0.0, 1e-14);
        assert!(report.absolute_energy_residual_j <= config.solve.thermal.total_energy_tolerance_j);
        if let Some((old_temperatures, old_masses, old_gas_t, old_vapor_mass)) = &baseline {
            for (a, b) in temperatures.iter().zip(old_temperatures) {
                near(*a, *b, 2e-8);
            }
            for (a, b) in masses.iter().zip(old_masses) {
                near(*a, *b, 1e-13);
            }
            near(gas_t, *old_gas_t, 2e-8);
            near(vapor_mass, *old_vapor_mass, 1e-13);
        } else {
            baseline = Some((temperatures, masses, gas_t, vapor_mass));
        }
    }
}

#[test]
fn registry_validation_precedes_sources_and_pure_liquid_admission_waits_for_birth() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let initial_gas = gas(&source, 300.0, 1000.0, 0.001);
    let original = scheduled(vec![birth(1.0, 1, 280.0, 0.0, 0.0)], &ORIGINAL_LIQUID, 64);
    let snapshot = original.clone();
    let spy = PhaseQueries(Cell::new(0));
    let mut declared = surface(1, &binding, 0.02, 20.0);
    declared.law = &spy;
    declared.pure_liquid_enthalpy_range_j_kg = [
        liquid(300.0).specific_enthalpy_j_kg,
        liquid(400.0).specific_enthalpy_j_kg,
    ];
    let before = advance_deposition_with_gas_timeline(
        &original,
        &[declared],
        &initial_gas,
        &source,
        0.5,
        event_policy(),
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    assert_eq!(before.gas(), &initial_gas);
    assert_eq!(spy.0.get(), 0);
    let refused = advance_deposition_with_gas_timeline(
        &original,
        &[declared],
        &initial_gas,
        &source,
        1.0,
        event_policy(),
        empty_heat,
        || Ok(()),
    );
    assert!(matches!(
        refused,
        Err(AdaptiveImplicitGasError::Implicit(
            ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
                DepositionError::MaterialDomain
            ))
        ))
    ));
    assert_eq!(
        spy.0.get(),
        0,
        "the original h branch refuses before querying a phase source"
    );

    // Even an already-active first descriptor must not query its law before
    // a later unknown or duplicate declaration has passed structural checks.
    let active = thermal(vec![cell(1, 300.0, 0.0, AREA)], None);
    let mut valid = surface(1, &binding, 0.02, 20.0);
    valid.law = &spy;
    let mut unknown = valid;
    unknown.cell_id = CellId(99);
    for surfaces in [[valid, unknown], [valid, valid]] {
        let calls = Cell::new(0usize);
        let refused = advance_deposition_with_gas_timeline(
            &active,
            &surfaces,
            &initial_gas,
            &source,
            0.5,
            event_policy(),
            |_, _, _| {
                calls.set(calls.get() + 1);
                Ok(coefficients(()))
            },
            || Ok(()),
        );
        assert!(matches!(refused, Err(AdaptiveImplicitGasError::Invalid(_))));
        assert_eq!(calls.get(), 0);
        assert_eq!(spy.0.get(), 0);
    }
    unchanged(&original, &snapshot, &initial_gas, &initial_gas.clone());
    assert_eq!(original.pending_events(), snapshot.pending_events());
}

#[test]
fn event_cuts_do_not_renew_work_and_late_source_or_final_cancel_refuses_atomically() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let surfaces = declarations(&binding);
    let original = scheduled(history(0.0), &ORIGINAL_LIQUID, 64);
    let snapshot = original.clone();
    let initial_gas = gas(&source, 320.0, 1600.0, 0.001);
    let gas_snapshot = initial_gas.clone();
    let config = event_policy();
    let born = at_events(&original, 0.5);
    let first = advance_deposition_with_gas_adaptive(
        &born,
        &[surfaces[2]],
        &initial_gas,
        &source,
        1.0,
        config,
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    assert_eq!(first.report().accepted_intervals, 2);
    for budget in 0..4 {
        let mut limited = config;
        let thermal_limit = if budget == 0 { 2 } else { 64 };
        match budget {
            1 => limited.solve.thermal.maximum_evaluations = first.report().residual_evaluations,
            2 => limited.maximum_gas_iterations = first.report().gas_iterations,
            3 => limited.maximum_comparisons = 1,
            _ => (),
        }
        let input = scheduled(history(0.0), &ORIGINAL_LIQUID, thermal_limit);
        let input_snapshot = input.clone();
        let mut queries = Vec::new();
        let refused = advance_deposition_with_gas_timeline(
            &input,
            &surfaces,
            &initial_gas,
            &source,
            2.0,
            limited,
            |thermal, _, _| {
                queries.push(thermal.time_s());
                Ok(coefficients(()))
            },
            || Ok(()),
        );
        if budget == 3 {
            assert!(matches!(
                refused,
                Err(AdaptiveImplicitGasError::ComparisonBudget)
            ));
        } else {
            let expected = match budget {
                0 => "thermal substeps",
                1 => "adaptive finite-gas whole-window solver work",
                _ => "adaptive finite-gas whole-window gas iterations",
            };
            match refused {
                Err(AdaptiveImplicitGasError::Implicit(ImplicitGasError::Thermal(
                    ThermalEvaporationError::Thermal(DepositionError::Budget(actual)),
                ))) => assert_eq!(actual, expected),
                _ => panic!("event boundaries must retain the spent whole-window budget"),
            }
        }
        assert!(!queries.is_empty());
        assert!(queries.iter().all(|time| *time <= 1.0));
        unchanged(&input, &input_snapshot, &initial_gas, &gas_snapshot);
        assert_eq!(input.pending_events(), input_snapshot.pending_events());
    }
    let mut late = false;
    let refused = advance_deposition_with_gas_timeline(
        &original,
        &surfaces,
        &initial_gas,
        &source,
        2.0,
        config,
        |thermal, liquids, _| {
            if thermal.time_s() > 1.5 {
                assert!(liquids.contains_key(&CellId(2)));
                assert_eq!(thermal.exposed_area_m2(CellId(2))?, 0.0);
                late = true;
                Err(DepositionError::InvalidInput("heat source refused after covering").into())
            } else {
                Ok(coefficients(()))
            }
        },
        || Ok(()),
    );
    assert!(late);
    assert!(matches!(
        refused,
        Err(AdaptiveImplicitGasError::Implicit(
            ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
                DepositionError::InvalidInput("heat source refused after covering")
            ))
        ))
    ));
    let mut polls = 0usize;
    advance_deposition_with_gas_timeline(
        &original,
        &surfaces,
        &initial_gas,
        &source,
        2.0,
        config,
        empty_heat,
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    let mut count = 0usize;
    let cancelled = advance_deposition_with_gas_timeline(
        &original,
        &surfaces,
        &initial_gas,
        &source,
        2.0,
        config,
        empty_heat,
        || {
            count += 1;
            if count == polls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        cancelled,
        Err(AdaptiveImplicitGasError::Implicit(
            ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(DepositionError::Cancelled))
        ))
    ));
    assert_eq!(count, polls);
    unchanged(&original, &snapshot, &initial_gas, &gas_snapshot);
    assert_eq!(original.pending_events(), snapshot.pending_events());
}
