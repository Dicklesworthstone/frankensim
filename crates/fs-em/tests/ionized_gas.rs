//! G1/G3: equilibrium ionization through the original FEEC current/heat owner.
mod common;
#[path = "../../fs-material/tests/common/ionized_gas.rs"]
mod fixture;
use common::{config, ends, layered, limits, near};
use fs_em::conduction::{ConductionControls, ConductionError, PlanarConduction};
use fs_material::ionized_gas::IonizedGasConductivity;

fn controls() -> ConductionControls {
    ConductionControls {
        linear: config(),
        current_tolerance_a: 1e-6,
        power_tolerance_w: 1e-6,
    }
}

#[test]
fn g1_ionized_material_uses_original_field_current_joule_and_circuit_laws() {
    for t in [6000.0, 8000.0] {
        let gas = fixture::material(t, 100_000.0);
        let sigma = gas.conductivity_s_m();
        let problem =
            PlanarConduction::new(layered(0.5, 1.0, 0.2), vec![gas; 4], limits(), &mut || {
                false
            })
            .unwrap();
        let field = problem
            .solve(&ends(2.0), controls(), &mut || false)
            .unwrap();
        // Unit length and width, depth=0.2: E=(-2,0), P=sigma E² volume.
        for (e, j) in field
            .electric_field_v_m()
            .iter()
            .zip(field.current_density_a_m2())
        {
            near(e[0], -2.0);
            near(e[1], 0.0);
            near(j[0], -2.0 * sigma);
            near(j[1], 0.0);
        }
        near(field.power_w(), 0.8 * sigma);
        assert!(field.power_balance_residual_w().abs() < 1e-6);
        let resistor = problem
            .resistance(&[2, 5], &[0, 3], controls(), &mut || false)
            .unwrap();
        near(resistor.resistance_ohm(), 1.0 / (0.2 * sigma));
        let reversed = resistor
            .field_at_voltage(-2.0, 1e-6, 1e-6, &mut || false)
            .unwrap();
        near(reversed.power_w(), field.power_w());
    }
}

#[test]
fn g0_actual_pde_and_circuit_fields_cannot_cross_the_weak_drift_domain() {
    let gas = fixture::material(8000.0, 100_000.0);
    let excessive = 2.0 * gas.maximum_electric_field_v_m();
    let problem =
        PlanarConduction::new(layered(0.5, 1.0, 0.2), vec![gas; 4], limits(), &mut || {
            false
        })
        .unwrap();
    assert!(matches!(
        problem.solve(&ends(excessive), controls(), &mut || false),
        Err(ConductionError::Invalid(_))
    ));
    let resistor = problem
        .resistance(&[2, 5], &[0, 3], controls(), &mut || false)
        .unwrap();
    assert!(matches!(
        resistor.field_at_voltage(excessive, 1e-6, 1e-6, &mut || false),
        Err(ConductionError::Invalid(_))
    ));
}

#[test]
fn g1_small_physical_drive_does_not_require_an_out_of_domain_unit_voltage_basis() {
    let chemistry = fixture::chemistry();
    let equilibrium = fixture::equilibrium(&chemistry, 8000.0, 100_000.0);
    let mut domain = fixture::limits();
    domain.maximum_drift_thermal_ratio = 1e-8;
    let gas = IonizedGasConductivity::new(
        &chemistry,
        equilibrium,
        fixture::rates(8000.0, 100_000.0, 1.0),
        domain,
        || false,
    )
    .unwrap();
    let sigma = gas.conductivity_s_m();
    let limit = gas.maximum_electric_field_v_m();
    assert!(limit < 0.001);
    let problem =
        PlanarConduction::new(layered(0.5, 1.0, 0.2), vec![gas; 4], limits(), &mut || {
            false
        })
        .unwrap();
    assert!(
        problem
            .resistance(&[2, 5], &[0, 3], controls(), &mut || false)
            .is_err()
    );
    let reference_v = 0.1 * limit;
    let resistor = problem
        .resistance_at_reference_voltage(&[2, 5], &[0, 3], reference_v, controls(), &mut || false)
        .unwrap();
    fixture::relative(resistor.reference_voltage_v(), reference_v, 1e-14);
    fixture::relative(
        resistor.reference_field().electric_field_v_m()[0][0],
        -reference_v,
        1e-10,
    );
    fixture::relative(resistor.resistance_ohm(), 1.0 / (0.2 * sigma), 1e-10);
    let actual_v = 0.5 * limit;
    let field = resistor
        .field_at_voltage(actual_v, 1e-8, 1e-10, &mut || false)
        .unwrap();
    fixture::relative(field.power_w(), 0.2 * sigma * actual_v * actual_v, 1e-10);
    fixture::relative(field.current_density_a_m2()[0][0], -sigma * actual_v, 1e-10);
    let current = 0.2 * sigma * actual_v;
    let driven = problem
        .solve_current_with_reference_voltage(
            &[2, 5],
            &[0, 3],
            current,
            reference_v,
            controls(),
            &mut || false,
        )
        .unwrap();
    fixture::relative(driven.power_w(), field.power_w(), 1e-10);
    assert!(
        resistor
            .field_at_voltage(2.0 * limit, 1e-8, 1e-10, &mut || false)
            .is_err()
    );
    assert!(
        problem
            .solve_current_with_reference_voltage(
                &[2, 5],
                &[0, 3],
                4.0 * current,
                reference_v,
                controls(),
                &mut || false,
            )
            .is_err()
    );
}

#[test]
fn g1_finite_rate_endpoint_drives_native_current_without_an_lte_projection() {
    use fs_gas::chemistry::{IsochoricKineticsSolve, advance_adiabatic};
    use fs_qty::{Amount, Temperature, Time, Volume};
    use fs_thermochem::kinetics::{
        ArrheniusRateV1, ElementaryReactionV1, IdealGasKineticsV1, ReactionProgressRateQuantityV1,
    };
    use fs_thermochem::{MolarEnergyQuantityV1, ReactionId};
    let chemistry = fixture::chemistry();
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("manufactured-electron-impact").unwrap(),
        vec![(fixture::id("H"), 1), (fixture::id("e-"), 1)],
        vec![(fixture::id("H+"), 1), (fixture::id("e-"), 2)],
        ArrheniusRateV1::new(
            Temperature::new(16000.0),
            ReactionProgressRateQuantityV1::new(20.0),
            0.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap(),
    )
    .unwrap();
    let mechanism = IdealGasKineticsV1::new(
        chemistry,
        "manufactured-ionized-transport",
        1,
        vec![reaction],
        1e-12,
    )
    .unwrap();
    let endpoint = advance_adiabatic(
        &mechanism,
        &[
            (fixture::id("H"), Amount::new(0.999)),
            (fixture::id("H+"), Amount::new(0.001)),
            (fixture::id("e-"), Amount::new(0.001)),
        ],
        Volume::new(1.0),
        Temperature::new(16000.0),
        Time::new(0.01),
        IsochoricKineticsSolve {
            minimum_temperature_k: 2000.0,
            maximum_temperature_k: 20000.0,
            energy_tolerance_j: 1e-6,
            relative_tolerance: 1e-5,
            absolute_amount_tolerance_mol: 1e-12,
            entropy_tolerance_j_k: 1e-8,
            initial_step_s: 1e-4,
            maximum_step_s: 0.002,
            maximum_steps: 10000,
            ..IsochoricKineticsSolve::default()
        },
        || false,
    )
    .unwrap();
    let state = endpoint.frozen_state(1e-12, || false).unwrap();
    let electron_density = state.number_density_per_m3(2).unwrap();
    assert!(endpoint.amounts()[2].value() > 0.001);
    assert!(endpoint.temperature().value() < 16000.0);
    let equilibrium = fixture::equilibrium(
        mechanism.thermodynamics(),
        endpoint.temperature().value(),
        endpoint.pressure().value(),
    );
    let kinetic_fraction =
        endpoint.amounts()[2].value() / endpoint.amounts().iter().map(|n| n.value()).sum::<f64>();
    assert!(equilibrium.species()[2].mole_fraction > 10.0 * kinetic_fraction);
    let electron = fixture::electron_relaxation(&state, 1e-20);
    let expected_sigma = electron_density * fixture::E * fixture::E
        / (fixture::ME * electron.rate_per_s())
        + electron_density * fixture::E * fixture::E / ((fixture::MH - fixture::ME) * 1e9);
    let ion = fixture::source_two_temperature(
        "H+",
        endpoint.temperature().value(),
        endpoint.temperature().value(),
        state.pressure().value(),
        1e9,
    );
    let gas =
        IonizedGasConductivity::from_state(state, vec![electron, ion], fixture::limits(), || false)
            .unwrap();
    assert!(gas.equilibrium().is_none());
    assert_eq!(gas.state().amounts(), endpoint.amounts());
    fixture::relative(gas.conductivity_s_m(), expected_sigma, 1e-13);
    let problem =
        PlanarConduction::new(layered(0.5, 1.0, 0.2), vec![gas; 4], limits(), &mut || {
            false
        })
        .unwrap();
    let field = problem
        .solve(&ends(2.0), controls(), &mut || false)
        .unwrap();
    fixture::relative(
        field.current_density_a_m2()[0][0],
        -2.0 * expected_sigma,
        1e-12,
    );
    fixture::relative(field.power_w(), 0.8 * expected_sigma, 1e-12);
    let circuit = problem
        .resistance(&[2, 5], &[0, 3], controls(), &mut || false)
        .unwrap();
    fixture::relative(
        circuit.resistance_ohm(),
        1.0 / (0.2 * expected_sigma),
        1e-12,
    );
    let reversed = circuit
        .field_at_voltage(-2.0, 1e-6, 1e-6, &mut || false)
        .unwrap();
    fixture::relative(reversed.power_w(), field.power_w(), 1e-12);
}

#[test]
fn g3_electron_temperature_changes_cross_section_drag_at_fixed_populations() {
    let mut electron_conductivities = Vec::new();
    for te in [6000.0, 12000.0] {
        let state = fixture::frozen_state(0.001, 6000.0, te, 1.0);
        let electron = fixture::electron_relaxation(&state, 1e-20);
        let ion = fixture::source_two_temperature("H+", 6000.0, te, state.pressure().value(), 1e9);
        let gas = IonizedGasConductivity::from_state(
            state,
            vec![electron, ion],
            fixture::limits(),
            || false,
        )
        .unwrap();
        electron_conductivities.push(
            gas.carriers()
                .iter()
                .find(|c| c.species() == &fixture::id("e-"))
                .unwrap()
                .conductivity_s_m(),
        );
        let sigma = gas.conductivity_s_m();
        let problem =
            PlanarConduction::new(layered(0.5, 1.0, 0.2), vec![gas; 4], limits(), &mut || {
                false
            })
            .unwrap();
        let field = problem
            .solve(&ends(1.0), controls(), &mut || false)
            .unwrap();
        fixture::relative(field.power_w(), 0.2 * sigma, 1e-12);
    }
    // Constant sigma_mt gives nu_e proportional to sqrt(Te). Actual populations
    // stay fixed; there is no fictitious LTE increase in carrier abundance.
    fixture::relative(
        electron_conductivities[0],
        2.0_f64.sqrt() * electron_conductivities[1],
        2e-9,
    );
}

#[test]
fn g1_native_joule_work_heats_both_pools_and_refreshes_the_same_field_law() {
    use fs_gas::chemistry::energy::{ElectronHeavyEnergySolve, advance_frozen_energy};
    use fs_qty::{Power, Temperature, Time};
    use fs_thermochem::electron_collision::MaxwellianRateControlsV1;
    use fs_thermochem::equilibrium::composition::state::FrozenIonizedGasStateV1;
    let initial = fixture::frozen_state(0.001, 6000.0, 6000.0, 0.2);
    let material = |state: &FrozenIonizedGasStateV1| {
        let electron = fixture::electron_relaxation(state, 1e-20);
        let ion = fixture::source_two_temperature(
            "H+",
            state.heavy_temperature().value(),
            state.electron_temperature().value(),
            state.pressure().value(),
            1e9,
        );
        IonizedGasConductivity::from_state(
            state.clone(),
            vec![electron, ion],
            fixture::limits(),
            || false,
        )
        .unwrap()
    };
    let gas = material(&initial);
    let initial_sigma = gas.conductivity_s_m();
    let problem =
        PlanarConduction::new(layered(0.5, 1.0, 0.2), vec![gas; 4], limits(), &mut || {
            false
        })
        .unwrap();
    let field = problem
        .solve(&ends(20.0), controls(), &mut || false)
        .unwrap();
    let mut electron_power = 0.0;
    let mut heavy_power = 0.0;
    let mut total_volume = 0.0;
    for (cell, e) in field.electric_field_v_m().iter().enumerate() {
        let heat = problem.conductors()[cell].joule_heating_w_m3(*e).unwrap();
        let volume = problem.mesh().face_measure(cell).unwrap();
        electron_power += volume * heat.electron_power_density_w_m3();
        heavy_power += volume * heat.heavy_power_density_w_m3();
        total_volume += volume;
        fixture::relative(
            volume * heat.total_power_density_w_m3(),
            field.cell_power_w()[cell],
            1e-13,
        );
    }
    fixture::relative(total_volume, initial.volume().value(), 1e-14);
    fixture::relative(electron_power + heavy_power, field.power_w(), 1e-13);
    assert!(electron_power > heavy_power && heavy_power > 0.0);
    let mut energy_controls = ElectronHeavyEnergySolve::new(
        Temperature::new(2000.0),
        Temperature::new(20000.0),
        0.2,
        0.1,
        MaxwellianRateControlsV1 {
            relative_tolerance: 1e-10,
            maximum_subdivisions: 1024,
        },
    );
    energy_controls.energy_tolerance_j = 1e-7;
    energy_controls.temperature_tolerance_k = 1e-7;
    energy_controls.initial_step_s = 1e-8;
    energy_controls.maximum_step_s = 1e-7;
    let duration = 1e-6;
    let advanced = advance_frozen_energy(
        &initial,
        &[fixture::momentum_cross_section("H", 1e-20)],
        Power::new(electron_power),
        Power::new(heavy_power),
        Time::new(duration),
        energy_controls,
        || false,
    )
    .unwrap();
    let next = advanced.state();
    assert_eq!(next.amounts(), initial.amounts());
    assert_eq!(next.volume(), initial.volume());
    assert!(next.electron_temperature().value() > 6000.0);
    assert!(next.heavy_temperature().value() > 6000.0);
    assert!(
        (next.internal_energy().value()
            - initial.internal_energy().value()
            - (electron_power + heavy_power) * duration)
            .abs()
            < 1e-6
    );
    let updated = material(next);
    assert!(updated.conductivity_s_m() < initial_sigma);
    let refreshed = PlanarConduction::new(
        layered(0.5, 1.0, 0.2),
        vec![updated; 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let next_field = refreshed
        .solve(&ends(20.0), controls(), &mut || false)
        .unwrap();
    assert!(next_field.power_w() < field.power_w());
}
