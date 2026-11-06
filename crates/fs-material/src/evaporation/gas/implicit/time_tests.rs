//! Dry-species admission, independent continuous limits, and event boundaries.
#[path = "adaptive_tests.rs"]
mod adaptive;

use super::*;

#[test]
fn initially_absent_vapor_reenters_from_the_original_catalog_without_trace_inventory() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let dry = gas(&source, 300.0, 0.0, 1e-4);
    // Retain only the actual dry carrier in the parcel. The separate original
    // catalog must supply the absent volatile's NASA9 and transport sources.
    let carrier_mass = dry.species_mass_kg(&carrier()).unwrap();
    let initial = dry
        .at_species_inventory_temperature(
            &source,
            &BTreeMap::from([(carrier(), carrier_mass)]),
            300.0,
        )
        .unwrap();
    assert_eq!(initial.model().thermodynamics().components().len(), 1);
    assert!(
        !initial
            .species_inventory_kg()
            .unwrap()
            .contains_key(&vapor())
    );
    assert_eq!(initial.species_partial_pressure_pa(&vapor()).unwrap(), 0.0);
    let old = thermal(vec![cell(1, 300.0, 2.0, 0.001)], None);
    let snapshot = old.clone();
    let gas_snapshot = initial.clone();
    let mut trial_vapor = Vec::new();
    let candidate = solve_deposition_with_gas(
        &old,
        &[surface(1, &binding, 0.02, 20.0)],
        &initial,
        &source,
        2.0,
        controls(),
        |_, _, gas| {
            trial_vapor.push(gas.species_mass_kg(&vapor())?);
            Ok(coefficients(()))
        },
        || Ok(()),
    )
    .unwrap();
    assert_eq!(trial_vapor[0], 0.0);
    assert!(trial_vapor.iter().skip(1).any(|&mass| mass > 0.0));
    let report = candidate.report();
    let receipt = &report.thermal.thermal.surfaces[&CellId(1)];
    let gas = candidate.gas();
    assert!(receipt.transferred_mass_kg > 1e-8);
    assert!(receipt.liquid.mass_kg < LIQUID_MASS);
    assert!(gas.species_partial_pressure_pa(&vapor()).unwrap() > 0.0);
    assert_eq!(gas.model().thermodynamics().components().len(), 2);
    assert_eq!(gas.model().water_species(), Some(&vapor()));
    assert_eq!(gas.volume_m3(), initial.volume_m3());
    assert_eq!(gas.point().axes().get("local-lot"), Some(&11.0));
    assert!(std::ptr::eq(candidate.catalog(), &source));
    near(
        gas.species_mass_kg(&carrier()).unwrap(),
        carrier_mass,
        1e-14,
    );
    near(
        gas.species_mass_kg(&vapor()).unwrap(),
        receipt.transferred_mass_kg,
        1e-14,
    );
    near(
        candidate.thermal().mass_kg() + gas.mass_kg(),
        old.mass_kg() + initial.mass_kg(),
        1e-14,
    );
    let gas_energy = gas.mass_kg() * gas.state().specific_internal_energy_j_kg()
        - initial.mass_kg() * initial.state().specific_internal_energy_j_kg();
    near(
        candidate.thermal().energy_j() - old.energy_j() + gas_energy,
        4.0,
        1e-7,
    );
    near(report.external_heat_j, 4.0, 1e-10);
    near(report.energy_residual_j, 0.0, 1e-7);
    near(report.gas_mass_residual_kg, 0.0, 1e-14);
    unchanged(&old, &snapshot, &initial, &gas_snapshot);
}

struct Continuous {
    carrier_mass: f64,
    initial_vapor_mass: f64,
    volume: f64,
    source_w: f64,
}
impl Continuous {
    const AREA: f64 = 0.001;
    const MASS_COEFFICIENT: f64 = 0.02;
    const HEAT_COEFFICIENT: f64 = 60.0;
    const DURATION: f64 = 10.0;

    // Three physical unknowns: liquid temperature, volatile gas mass, and gas
    // temperature. Actual liquid mass follows the exact shared mass invariant.
    // Independent energy differentiation gives m_l*c_l*T_l' = Q-r*(h_up-h_l)
    // and C_g*T_g' = r*(h_up-u_v)-Q_film. The upstream h is selected by the
    // signed stream, so condensation retains hot-gas sensible vapor enthalpy.
    fn derivative(&self, state: [f64; 3]) -> [f64; 3] {
        let [liquid_t, vapor_mass, gas_t] = state;
        let liquid_mass = LIQUID_MASS - (vapor_mass - self.initial_vapor_mass);
        let pressure = (self.carrier_mass * CARRIER_R + vapor_mass * R) * gas_t / self.volume;
        let vapor_pressure = vapor_mass * R * gas_t / self.volume;
        let rate = Self::AREA
            * Self::MASS_COEFFICIENT
            * unit_film_rate(liquid_t, gas_t, pressure, vapor_pressure);
        let heat = Self::AREA * Self::HEAT_COEFFICIENT * (gas_t - liquid_t);
        let upstream = vapor_h(if rate >= 0.0 { liquid_t } else { gas_t });
        let liquid_capacity = liquid_mass * R * (4.0 + 0.001 * liquid_t);
        let gas_capacity =
            self.carrier_mass * CARRIER_R * 2.5 + vapor_mass * R * (3.0 + 0.001 * gas_t);
        assert!(liquid_mass > 0.0 && vapor_mass >= 0.0);
        assert!(liquid_capacity > 0.0 && gas_capacity > 0.0);
        [
            (self.source_w + heat - rate * (upstream - vapor_h(liquid_t) + LATENT))
                / liquid_capacity,
            rate,
            (rate * (upstream - vapor_u(gas_t)) - heat) / gas_capacity,
        ]
    }

    fn energy(&self, [liquid_t, vapor_mass, gas_t]: [f64; 3]) -> f64 {
        let liquid_mass = LIQUID_MASS - (vapor_mass - self.initial_vapor_mass);
        liquid_mass * (vapor_h(liquid_t) - LATENT)
            + self.carrier_mass * CARRIER_R * 2.5 * gas_t
            + vapor_mass * vapor_u(gas_t)
    }

    // A small test-only classical RK4 oracle. No production owner or implicit
    // residual is called: only the explicit original fixture equations above.
    fn integrate(&self, mut state: [f64; 3], steps: usize) -> [f64; 3] {
        let dt = Self::DURATION / steps as f64;
        let increment = |state: [f64; 3], slope: [f64; 3], span: f64| {
            std::array::from_fn(|i| state[i] + span * slope[i])
        };
        for _ in 0..steps {
            let a = self.derivative(state);
            let b = self.derivative(increment(state, a, 0.5 * dt));
            let c = self.derivative(increment(state, b, 0.5 * dt));
            let d = self.derivative(increment(state, c, dt));
            state = std::array::from_fn(|i| {
                state[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0
            });
        }
        state
    }
}

#[test]
fn evaporation_and_condensation_refine_at_first_order_toward_the_independent_continuous_balance() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    // Both directions are active throughout their window. The second case has
    // a large gas/liquid temperature difference, making the condensation
    // upstream sensible enthalpy and variable gas capacity consequential.
    for (gas_t, vapor_pressure, power) in [(300.0, 1000.0, 2.0), (340.0, 5000.0, 0.0)] {
        let initial_gas = gas(&source, gas_t, vapor_pressure, 0.001);
        let initial_thermal = thermal(vec![cell(1, 300.0, power, Continuous::AREA)], None);
        let oracle = Continuous {
            carrier_mass: initial_gas.species_mass_kg(&carrier()).unwrap(),
            initial_vapor_mass: initial_gas.species_mass_kg(&vapor()).unwrap(),
            volume: initial_gas.volume_m3(),
            source_w: power,
        };
        let start = [300.0, oracle.initial_vapor_mass, gas_t];
        let reference = oracle.integrate(start, 2048);
        let coarser_reference = oracle.integrate(start, 1024);
        near(reference[0], coarser_reference[0], 1e-9);
        near(reference[1], coarser_reference[1], 1e-15);
        near(reference[2], coarser_reference[2], 1e-9);
        near(
            oracle.energy(reference) - oracle.energy(start),
            power * Continuous::DURATION,
            1e-7,
        );
        let direction = if power > 0.0 { 1.0 } else { -1.0 };
        assert!(direction * oracle.derivative(start)[1] > 0.0);
        assert!(direction * oracle.derivative(reference)[1] > 0.0);
        let mut errors = Vec::new();
        for steps in [4, 8, 16] {
            let mut thermal = initial_thermal.clone();
            let mut gas = initial_gas.clone();
            let mut external_heat = 0.0;
            for step in 0..steps {
                let target = Continuous::DURATION * (step + 1) as f64 / steps as f64;
                let candidate = solve_deposition_with_gas(
                    &thermal,
                    &[surface(
                        1,
                        &binding,
                        Continuous::MASS_COEFFICIENT,
                        Continuous::HEAT_COEFFICIENT,
                    )],
                    &gas,
                    &source,
                    target,
                    controls(),
                    empty_heat,
                    || Ok(()),
                )
                .unwrap();
                let report = candidate.report();
                assert!(
                    direction * report.thermal.thermal.surfaces[&CellId(1)].transferred_mass_kg
                        > 0.0
                );
                near(report.energy_residual_j, 0.0, 1e-7);
                external_heat += report.external_heat_j;
                let (next_thermal, next_gas, _, (), _) = candidate.into_parts();
                thermal = next_thermal;
                gas = next_gas;
            }
            let state = [
                thermal.temperature_k(CellId(1)).unwrap(),
                gas.species_mass_kg(&vapor()).unwrap(),
                gas.state().state().temperature,
            ];
            let error: [f64; 3] = std::array::from_fn(|i| (state[i] - reference[i]).abs());
            errors.push(error);
            near(external_heat, power * Continuous::DURATION, 1e-9);
            near(
                thermal.mass_kg() + gas.mass_kg(),
                initial_thermal.mass_kg() + initial_gas.mass_kg(),
                1e-13,
            );
            let gas_energy = gas.mass_kg() * gas.state().specific_internal_energy_j_kg()
                - initial_gas.mass_kg() * initial_gas.state().specific_internal_energy_j_kg();
            near(
                thermal.energy_j() - initial_thermal.energy_j() + gas_energy,
                power * Continuous::DURATION,
                1e-7,
            );
            assert_eq!(thermal.time_s(), Continuous::DURATION);
        }
        assert!(errors[0][0] > 1e-4 && errors[0][1] > 1e-9 && errors[0][2] > 1e-2);
        // Check every physical coordinate separately. A temperature-only norm
        // must not conceal a wrong finite-humidity/mass limit.
        for pair in errors.windows(2) {
            for coordinate in 0..3 {
                let ratio = pair[1][coordinate] / pair[0][coordinate];
                assert!(
                    (0.4..0.65).contains(&ratio),
                    "non-first-order refinement: {errors:?}"
                );
            }
        }
    }
}

#[test]
fn interior_events_refuse_and_exact_endpoint_events_stay_pending_for_both_owners() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let initial_gas = gas(&source, 300.0, 1000.0, 0.001);
    let initial_cell = cell(1, 300.0, 1.0, 0.001);
    let mut changed = initial_cell.boundary;
    changed.source_w = 3.0;
    let original = DepositionSimulation::new(
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(initial_cell),
            },
            TimedEvent {
                time_s: 1.0,
                event: Event::Boundary(CellId(1), changed),
            },
            TimedEvent {
                time_s: 1.0,
                event: Event::Birth(cell(2, 310.0, 0.0, 0.001)),
            },
        ],
        &[(MaterialId(1), &ORIGINAL_LIQUID)],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 3,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-6,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let snapshot = original.clone();
    let gas_snapshot = initial_gas.clone();
    let surfaces = [surface(1, &binding, 0.02, 20.0)];
    assert!(matches!(
        solve_deposition_with_gas(
            &original,
            &surfaces,
            &initial_gas,
            &source,
            2.0,
            controls(),
            empty_heat,
            || Ok(()),
        ),
        Err(ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput(
                "process pending events or split implicit step at the next event"
            )
        )))
    ));
    let candidate = solve_deposition_with_gas(
        &original,
        &surfaces,
        &initial_gas,
        &source,
        1.0,
        controls(),
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    assert_eq!(candidate.thermal().time_s(), 1.0);
    assert_eq!(
        candidate.thermal().pending_events(),
        original.pending_events()
    );
    assert_eq!(candidate.thermal().cells().len(), 1);
    assert_eq!(
        candidate.thermal().cells()[&CellId(1)].boundary.source_w,
        1.0
    );
    near(candidate.report().external_heat_j, 1.0, 1e-10);
    assert!(matches!(
        solve_deposition_with_gas(
            candidate.thermal(),
            &surfaces,
            candidate.gas(),
            &source,
            2.0,
            controls(),
            empty_heat,
            || Ok(()),
        ),
        Err(ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput(
                "process pending events or split implicit step at the next event"
            )
        )))
    ));
    let retained_liquid = candidate.thermal().cells()[&CellId(1)].mass_kg;
    let (after_events, events) = candidate
        .thermal()
        .advance_to(
            1.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-6,
            },
            || Ok(()),
        )
        .unwrap();
    assert!(after_events.pending_events().is_empty());
    assert_eq!(after_events.cells().len(), 2);
    assert_eq!(after_events.cells()[&CellId(1)].boundary.source_w, 3.0);
    assert_eq!(after_events.cells()[&CellId(1)].mass_kg, retained_liquid);
    assert_eq!(events.incoming_mass_kg, LIQUID_MASS);
    assert_eq!(events.external_heat_j, 0.0);
    assert_eq!(original.pending_events(), snapshot.pending_events());
    unchanged(&original, &snapshot, &initial_gas, &gas_snapshot);
}
