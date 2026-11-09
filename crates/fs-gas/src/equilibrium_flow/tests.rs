use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::MolarMass;
use fs_thermochem::equilibrium::composition::EquilibriumSolveV1;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, NASA9_LAW_ID_V1,
    NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, StandardStateConventionV1, StandardStatePhaseV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K, required_nasa9_parameter_dims,
};
use std::collections::BTreeMap;

mod nozzle;
mod nozzle_prepared;
#[cfg(feature = "material")]
mod nozzle_material;

const R: f64 = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
const T0: f64 = 1000.0;
const P0: f64 = 100000.0;

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}

// Genuine NASA-9 source cards with constant Cp and declared reference datums.
// These are input models, not substitutes for the original equilibrium solver.
fn model(
    name: &str,
    mass: f64,
    cp_r: f64,
    formation: f64,
    entropy_offset: f64,
) -> Nasa9StandardStateModelV1 {
    let mut coefficients = [0.0; 9];
    coefficients[2] = cp_r;
    coefficients[7] = formation / R;
    coefficients[8] = cp_r * (1.0 - T0.ln()) + entropy_offset / R;
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(name.to_string(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: P0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.to_string()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 6000.0),
        sources: Vec::new(),
        provenance: Provenance {
            source: "manufactured reacting machinery oracle".into(),
            license: "test fixture".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        id(name),
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(P0),
            ElementalReferenceIdV1::new("reacting-flow-reference").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn owner(models: Vec<Nasa9StandardStateModelV1>, row: Vec<i128>) -> IdealGasEquilibriumV1 {
    let species: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let matrix = ElementalMatrix::new(
        vec![ElementId::new("H").unwrap()],
        species.clone(),
        vec![row],
    )
    .unwrap();
    let charge = ChargeVector::new(species.clone(), vec![0; species.len()]).unwrap();
    IdealGasEquilibriumV1::new(matrix, charge, models).unwrap()
}
fn single() -> IdealGasEquilibriumV1 {
    owner(vec![model("A", 0.01, 4.0, 0.0, 0.0)], vec![1])
}
fn dissociation() -> IdealGasEquilibriumV1 {
    owner(
        vec![
            model("A", 0.001, 2.5, 0.0, 0.0),
            model("A2", 0.002, 3.5, 0.0, 0.0),
        ],
        vec![1, 2],
    )
}
fn controls() -> EquilibriumFlowSolve {
    let equilibrium = EquilibriumSolveV1 {
        relative_balance_tolerance: 1e-13,
        dimensionless_potential_tolerance: 1e-12,
        ..EquilibriumSolveV1::default()
    };
    EquilibriumFlowSolve {
        enthalpy: EquilibriumEnthalpySolveV1 {
            equilibrium,
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 3000.0,
            enthalpy_tolerance_j: 1e-7,
            maximum_temperature_iterations: 80,
        },
        isentropic: EquilibriumEntropySolveV1 {
            equilibrium,
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 3000.0,
            entropy_tolerance_j_k: 1e-10,
            maximum_temperature_iterations: 80,
        },
        energy_tolerance_j: 2e-7,
        entropy_production_tolerance_j_k: 1e-9,
    }
}
fn solve(
    chemistry: &IdealGasEquilibriumV1,
    feed: &[(SpeciesId, Amount)],
    pressure: f64,
    process: EquilibriumPressureProcess,
) -> EquilibriumPressureChange<EquilibriumThermodynamicStateV1> {
    equilibrium_pressure_change(
        chemistry,
        Temperature::new(T0),
        Pressure::new(P0),
        feed,
        Pressure::new(pressure),
        process,
        controls(),
        || false,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:e}"
    );
}
fn bisect(target: f64, value: impl Fn(f64) -> f64) -> f64 {
    let (mut low, mut high) = (200.0, 3000.0);
    assert!(value(low) <= target && value(high) >= target);
    for _ in 0..85 {
        let middle = 0.5 * (low + high);
        if value(middle) < target {
            low = middle;
        } else {
            high = middle;
        }
    }
    0.5 * (low + high)
}

#[test]
fn frozen_limit_matches_analytic_compressors_expanders_and_stream_rates() {
    let chemistry = single();
    let feed = [(id("A"), Amount::new(2.0))];
    for (pressure, compressor) in [(4.0 * P0, true), (P0 / 4.0, false)] {
        for efficiency in [1.0, 0.72] {
            let process = if compressor {
                EquilibriumPressureProcess::Compressor {
                    isentropic_efficiency: efficiency,
                }
            } else {
                EquilibriumPressureProcess::Expander {
                    isentropic_efficiency: efficiency,
                }
            };
            let result = solve(&chemistry, &feed, pressure, process);
            let ideal_temperature = T0 * (pressure / P0).powf(0.25);
            let ideal_work = 2.0 * 4.0 * R * (ideal_temperature - T0);
            let work = if compressor {
                ideal_work / efficiency
            } else {
                ideal_work * efficiency
            };
            let temperature = T0 + work / (2.0 * 4.0 * R);
            let entropy = 2.0 * R * (4.0 * (temperature / T0).ln() - (pressure / P0).ln());
            close(
                result.outlet.equilibrium().temperature().value(),
                temperature,
                2e-9,
            );
            close(result.balance.shaft_work_j, work, 2e-7);
            close(result.balance.entropy_production_j_k, entropy, 2e-9);
            close(
                result
                    .isentropic_outlet
                    .as_ref()
                    .unwrap()
                    .equilibrium()
                    .temperature()
                    .value(),
                ideal_temperature,
                2e-9,
            );
            assert!(result.balance.energy_residual_j.abs() <= controls().energy_tolerance_j);
            if efficiency < 1.0 {
                assert!(result.balance.entropy_production_j_k > 0.0);
            }
            let rates = result.balance.rates(MassFlowRate::new(0.6)).unwrap();
            close(rates.shaft_power.value(), work * 30.0, 1e-5);
            close(rates.entropy_production_w_k, entropy * 30.0, 1e-7);
            assert_eq!(
                rates.energy_residual.value(),
                result.balance.energy_residual_j * 30.0
            );
            assert!(result.inlet.equilibrium().is_bound_to(&chemistry));
            assert!(result.outlet.equilibrium().is_bound_to(&chemistry));
        }
    }
}

fn alpha(temperature: f64, pressure: f64) -> f64 {
    let kp = (temperature / T0).powf(1.5);
    (kp / (kp + 4.0 * pressure / P0)).sqrt()
}
fn dissociation_h(temperature: f64, pressure: f64) -> f64 {
    R * temperature * (3.5 + 1.5 * alpha(temperature, pressure))
}
fn dissociation_s(temperature: f64, pressure: f64) -> f64 {
    let a = alpha(temperature, pressure);
    let n_a = 2.0 * a;
    let n_a2 = 1.0 - a;
    let n = n_a + n_a2;
    R * (n_a * (2.5 * (1.0 + (temperature / T0).ln()) - (n_a / n * pressure / P0).ln())
        + n_a2 * (3.5 * (1.0 + (temperature / T0).ln()) - (n_a2 / n * pressure / P0).ln()))
}
fn check_products(state: &EquilibriumThermodynamicStateV1, chemistry: &IdealGasEquilibriumV1) {
    let t = state.equilibrium().temperature().value();
    let p = state.equilibrium().pressure().value();
    let species = state.equilibrium().species();
    let amount = |name: &str| {
        species
            .iter()
            .find(|s| s.species.as_str() == name)
            .unwrap()
            .amount
            .value()
    };
    close(amount("A"), 2.0 * alpha(t, p), 5e-11);
    close(amount("A2"), 1.0 - alpha(t, p), 5e-11);
    close(amount("A") + 2.0 * amount("A2"), 2.0, 3e-12);
    assert!(state.equilibrium().is_bound_to(chemistry));
    assert!(species.iter().all(|s| s.source.is_some()));
}

#[test]
fn reacting_throttle_cools_with_zero_work_and_conserved_elements() {
    let chemistry = dissociation();
    let feed = [(id("A2"), Amount::new(1.0))];
    let p = 0.1 * P0;
    let result = solve(&chemistry, &feed, p, EquilibriumPressureProcess::Throttle);
    let t = bisect(dissociation_h(T0, P0), |t| dissociation_h(t, p));
    close(result.outlet.equilibrium().temperature().value(), t, 2e-8);
    close(
        result.balance.entropy_production_j_k,
        dissociation_s(t, p) - dissociation_s(T0, P0),
        2e-8,
    );
    assert!(t < T0 && result.balance.entropy_production_j_k > 0.0);
    assert_eq!(result.balance.shaft_work_j, 0.0);
    assert_eq!(result.balance.isentropic_enthalpy_change_j, None);
    assert_eq!(result.isentropic_outlet, None);
    assert!(result.balance.enthalpy_change_j.abs() <= controls().energy_tolerance_j);
    assert_eq!(result.balance.mass_kg, 0.002);
    check_products(&result.inlet, &chemistry);
    check_products(&result.outlet, &chemistry);
    let rates = result.balance.rates(MassFlowRate::new(2.0)).unwrap();
    assert_eq!(rates.shaft_power.value(), 0.0);
    assert!(rates.entropy_production_w_k > 0.0);
}

#[test]
fn reacting_machines_match_independent_mass_action_enthalpy_and_entropy_roots() {
    let chemistry = dissociation();
    let feed = [(id("A2"), Amount::new(1.0))];
    for (p, compressor) in [(5.0 * P0, true), (0.2 * P0, false)] {
        let ts = bisect(dissociation_s(T0, P0), |t| dissociation_s(t, p));
        let ideal = dissociation_h(ts, p) - dissociation_h(T0, P0);
        // Frozen composition predicts a distinctly different reversible T.
        let a0 = alpha(T0, P0);
        let frozen_exponent = (1.0 + a0) / (3.5 + 1.5 * a0);
        assert!((ts - T0 * (p / P0).powf(frozen_exponent)).abs() > 1.0);
        for efficiency in [1.0, 0.75] {
            let process = if compressor {
                EquilibriumPressureProcess::Compressor {
                    isentropic_efficiency: efficiency,
                }
            } else {
                EquilibriumPressureProcess::Expander {
                    isentropic_efficiency: efficiency,
                }
            };
            let work = if compressor {
                ideal / efficiency
            } else {
                ideal * efficiency
            };
            let t = bisect(dissociation_h(T0, P0) + work, |t| dissociation_h(t, p));
            let result = solve(&chemistry, &feed, p, process);
            close(
                result.balance.isentropic_enthalpy_change_j.unwrap(),
                ideal,
                3e-7,
            );
            close(result.balance.shaft_work_j, work, 4e-7);
            close(result.outlet.equilibrium().temperature().value(), t, 2e-8);
            close(
                result.balance.entropy_production_j_k,
                dissociation_s(t, p) - dissociation_s(T0, P0),
                2e-8,
            );
            if efficiency == 1.0 {
                assert!(
                    result.balance.entropy_production_j_k.abs()
                        <= controls().entropy_production_tolerance_j_k
                );
            } else {
                assert!(result.balance.entropy_production_j_k > 0.0);
            }
            assert!(
                result.balance.relative_mass_residual.abs()
                    <= controls().enthalpy.equilibrium.relative_balance_tolerance
            );
            check_products(&result.outlet, &chemistry);
            check_products(result.isentropic_outlet.as_ref().unwrap(), &chemistry);
        }
    }
}

#[test]
fn small_machine_work_survives_large_formation_and_entropy_references() {
    let feed = [(id("A"), Amount::new(1.0))];
    let p = P0 * (1.0 + 1e-8);
    let expected = 4.0 * R * T0 * (0.25 * (p / P0).ln()).exp_m1();
    for datum in [-1e11, 0.0, 1e11] {
        let chemistry = owner(vec![model("A", 0.01, 4.0, datum, datum / T0)], vec![1]);
        let mut c = controls();
        c.enthalpy.enthalpy_tolerance_j = 1e-9;
        c.energy_tolerance_j = 2e-9;
        let result = equilibrium_pressure_change(
            &chemistry,
            Temperature::new(T0),
            Pressure::new(P0),
            &feed,
            Pressure::new(p),
            EquilibriumPressureProcess::Compressor {
                isentropic_efficiency: 1.0,
            },
            c,
            || false,
        )
        .unwrap();
        close(result.balance.shaft_work_j, expected, 3e-9);
        assert!(result.balance.energy_residual_j.abs() <= 2e-9);
        assert!(result.balance.entropy_production_j_k.abs() <= 1e-9);
    }
}

#[test]
fn equal_pressure_is_exact_identity_and_avoids_sp_or_hp_work() {
    let chemistry = dissociation();
    let feed = [(id("A2"), Amount::new(1.0))];
    for process in [
        EquilibriumPressureProcess::Throttle,
        EquilibriumPressureProcess::Compressor {
            isentropic_efficiency: 0.4,
        },
        EquilibriumPressureProcess::Expander {
            isentropic_efficiency: 0.4,
        },
    ] {
        let result = solve(&chemistry, &feed, P0, process);
        assert_eq!(result.inlet, result.outlet);
        assert_eq!(result.balance.shaft_work_j, 0.0);
        assert_eq!(result.balance.enthalpy_change_j, 0.0);
        assert_eq!(result.balance.energy_residual_j, 0.0);
        assert_eq!(result.balance.entropy_production_j_k, 0.0);
        assert_eq!(result.balance.temperature_iterations, 0);
        assert_eq!(
            result.balance.composition_iterations,
            result.inlet.equilibrium().composition_iterations()
        );
        assert_eq!(result.isentropic_outlet, None);
    }
}

#[test]
fn refusals_preserve_original_source_budgets_domains_and_process_direction() {
    let chemistry = dissociation();
    let feed = [(id("A2"), Amount::new(1.0))];
    let call = |p, process, c| {
        equilibrium_pressure_change(
            &chemistry,
            Temperature::new(T0),
            Pressure::new(P0),
            &feed,
            Pressure::new(p),
            process,
            c,
            || false,
        )
    };
    assert!(matches!(
        call(2.0 * P0, EquilibriumPressureProcess::Throttle, controls()),
        Err(EquilibriumFlowError::Invalid(_))
    ));
    for efficiency in [0.0, -0.5, 1.1, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            call(
                2.0 * P0,
                EquilibriumPressureProcess::Compressor {
                    isentropic_efficiency: efficiency
                },
                controls()
            ),
            Err(EquilibriumFlowError::Invalid(_))
        ));
    }
    let compressor = EquilibriumPressureProcess::Compressor {
        isentropic_efficiency: 0.75,
    };
    assert!(matches!(
        call(0.5 * P0, compressor, controls()),
        Err(EquilibriumFlowError::Invalid(_))
    ));
    let mut c = controls();
    c.isentropic.equilibrium.maximum_iterations -= 1;
    assert!(matches!(
        call(2.0 * P0, compressor, c),
        Err(EquilibriumFlowError::Invalid(_))
    ));
    c = controls();
    c.enthalpy.minimum_temperature_k = 150.0;
    assert!(matches!(
        call(0.1 * P0, EquilibriumPressureProcess::Throttle, c),
        Err(EquilibriumFlowError::Chemistry {
            stage: EquilibriumFlowStage::Inlet,
            ..
        })
    ));
    c = controls();
    c.isentropic.maximum_temperature_iterations = 1;
    assert!(matches!(
        call(5.0 * P0, compressor, c),
        Err(EquilibriumFlowError::Chemistry {
            stage: EquilibriumFlowStage::Isentropic,
            ..
        })
    ));
    let result = solve(
        &chemistry,
        &feed,
        0.5 * P0,
        EquilibriumPressureProcess::Throttle,
    );
    assert!(result.balance.rates(MassFlowRate::new(0.0)).is_err());
    assert!(result.balance.rates(MassFlowRate::new(-1.0)).is_err());
    assert!(result.balance.rates(MassFlowRate::new(f64::MAX)).is_err());
}

#[test]
fn independent_physical_energy_gate_does_not_replace_native_residuals() {
    let chemistry = dissociation();
    let feed = [(id("A2"), Amount::new(1.0))];
    let mut c = controls();
    // Deliberately admit a gross native HP residual. The process still owes
    // its independent tight energy closure; this is never a heat source.
    c.enthalpy.enthalpy_tolerance_j = 1e5;
    c.energy_tolerance_j = 1e-20;
    let result = equilibrium_pressure_change(
        &chemistry,
        Temperature::new(T0),
        Pressure::new(P0),
        &feed,
        Pressure::new(12345.0),
        EquilibriumPressureProcess::Throttle,
        c,
        || false,
    );
    assert!(matches!(result, Err(EquilibriumFlowError::Energy { .. })));
}

#[test]
fn cancellation_through_final_publication_is_atomic_and_retry_is_deterministic() {
    let chemistry = dissociation();
    let feed = [(id("A2"), Amount::new(1.0))];
    let process = EquilibriumPressureProcess::Expander {
        isentropic_efficiency: 0.8,
    };
    let call = |cancelled: &mut dyn FnMut() -> bool| {
        equilibrium_pressure_change(
            &chemistry,
            Temperature::new(T0),
            Pressure::new(P0),
            &feed,
            Pressure::new(P0 / 2.0),
            process,
            controls(),
            cancelled,
        )
    };
    let mut count = 0;
    let expected = call(&mut || {
        count += 1;
        false
    })
    .unwrap();
    assert!(count > 20);
    for stop in [1, count / 2, count] {
        let mut polls = 0;
        assert!(matches!(
            call(&mut || {
                polls += 1;
                polls == stop
            }),
            Err(EquilibriumFlowError::Cancelled)
        ));
        assert_eq!(polls, stop);
    }
    let actual = call(&mut || false);
    assert_eq!(actual.unwrap(), expected);
    assert_eq!(feed[0].1.value(), 1.0);
}
