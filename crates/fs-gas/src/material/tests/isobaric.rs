use super::*;
use crate::chemistry::IsochoricKineticsMethod;
use crate::chemistry::heat::{HeatBath, ThermalConductance, material::FiniteHeatBody};
use crate::material::isobaric::MaterialIsobaricReactor;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_qty::Mass;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn g1_reaction_uses_enthalpy_cp_and_actual_expansion_work() {
    let (bank, mechanism, mut parcel) = fixture(-24_000.0, 2000.0);
    let initial = parcel.clone();
    let mut reactor =
        MaterialIsobaricReactor::new(bank.clone(), mechanism.clone(), controls(), || false)
            .unwrap();
    let bytes = reactor.workspace_bytes();
    let report = reactor
        .advance_adiabatic(&mut parcel, Time::new(0.05), || false)
        .unwrap();
    let fraction = parcel
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions()[1];
    let temperature = parcel.state().state().temperature;
    close(temperature, 400.0 + 24_000.0 * fraction / (4.0 * R), 1e-6);
    close(
        parcel.volume_m3() / initial.volume_m3(),
        temperature / 400.0,
        1e-12,
    );
    assert_eq!(parcel.mass_kg().to_bits(), initial.mass_kg().to_bits());
    assert_eq!(
        parcel.state().state().pressure.to_bits(),
        initial.state().state().pressure.to_bits()
    );
    assert_eq!(parcel.point().axes()["sample_id"], 7.0);
    assert_eq!(
        parcel.model().thermodynamics().components(),
        bank.thermodynamics().components()
    );
    assert_eq!(
        parcel.model().transport_materials(),
        bank.transport_materials()
    );
    assert!(parcel.state().state().sound_speed > initial.state().state().sound_speed);
    assert!(
        parcel.state().state().thermal_conductivity > initial.state().state().thermal_conductivity
    );
    parcel.state().state().validate_coherence().unwrap();
    close(report.material.enthalpy_change_j, 0.0, 2e-8);
    let actual_work = initial.state().state().pressure * (parcel.volume_m3() - initial.volume_m3());
    let actual_storage = parcel.mass_kg()
        * (parcel.state().specific_internal_energy_j_kg()
            - initial.state().specific_internal_energy_j_kg());
    assert!(actual_work > 0.0);
    close(report.material.energy.work_by_gas_j, actual_work, 1e-12);
    close(actual_storage + actual_work, 0.0, 2e-8);
    close(report.balance.freeze_boundary_work_residual_j, 0.0, 2e-8);
    assert_eq!(report.chemistry.balance().heat_input_j, 0.0);
    assert_eq!(report.chemistry.elapsed_time(), Time::new(0.05));
    // A rigid vessel uses Cv=3R. Comparing both independent balances prevents
    // a constant-pressure wrapper from silently reusing the rigid U constraint.
    let mut rigid = initial.clone();
    advance_closed_parcel_adiabatic(
        &mut rigid,
        &bank,
        &mechanism,
        Time::new(0.05),
        controls(),
        || false,
    )
    .unwrap();
    let rigid_fraction = rigid
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions()[1];
    close(
        rigid.state().state().temperature,
        400.0 + 24_000.0 * rigid_fraction / (3.0 * R),
        1e-6,
    );
    close(
        (temperature - 400.0)
            / fraction
            / ((rigid.state().state().temperature - 400.0) / rigid_fraction),
        0.75,
        1e-9,
    );
    // Reusing the prepared owner on an independent parcel gives the same state.
    assert_eq!(reactor.workspace_bytes(), bytes);
    let accepted = parcel;
    let mut replay = initial;
    reactor
        .advance_adiabatic(&mut replay, Time::new(0.05), || false)
        .unwrap();
    assert_eq!(replay, accepted);
    assert_eq!(reactor.workspace_bytes(), bytes);
}

#[test]
fn g1_bath_heat_and_temperature_dependent_reaction_share_the_same_time() {
    let (bank, original, mut parcel) = fixture(0.0, 2000.0);
    let rate = ArrheniusRateV1::new(
        Temperature::new(400.0),
        ReactionProgressRateQuantityV1::new(4.0),
        1.0,
        MolarEnergyQuantityV1::new(0.0),
    )
    .unwrap();
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("temperature-driven-isomerization").unwrap(),
        vec![(id("A"), 1)],
        vec![(id("B"), 1)],
        rate,
    )
    .unwrap();
    let mechanism = IdealGasKineticsV1::new(
        original.thermodynamics().clone(),
        "manufactured-thermal-rate",
        1,
        vec![reaction],
        0.0,
    )
    .unwrap();
    let initial = parcel.clone();
    let mut control = controls();
    control.chemistry.method = IsochoricKineticsMethod::Sdirk2;
    let mut reactor = MaterialIsobaricReactor::new(bank, mechanism, control, || false).unwrap();
    let (duration, conductance, bath_temperature) = (0.1, 0.8, 700.0);
    let bath = HeatBath::new(
        ThermalConductance::new(conductance),
        Temperature::new(bath_temperature),
    )
    .unwrap();
    let report = reactor
        .advance_with_heat_bath(&mut parcel, Time::new(duration), bath, || false)
        .unwrap();
    // Equal source enthalpies make the temperature ODE independent of extent.
    // k(T)=4*T/400 has an independently integrable history on that exact ODE.
    let amount = initial.mass_kg() / 0.01;
    let capacity = amount * 4.0 * R;
    let decay = conductance / capacity;
    let exponential = (-decay * duration).exp();
    let expected_t = bath_temperature + (400.0 - bath_temperature) * exponential;
    let rate_integral = 4.0 / 400.0
        * (bath_temperature * duration + (400.0 - bath_temperature) * (1.0 - exponential) / decay);
    let expected_b = 0.5 * (1.0 - (-2.0 * rate_integral).exp());
    close(parcel.state().state().temperature, expected_t, 3e-4);
    close(
        parcel
            .model()
            .thermodynamics()
            .mole_composition()
            .fractions()[1],
        expected_b,
        3e-6,
    );
    close(
        report.chemistry.balance().heat_input_j,
        capacity * (parcel.state().state().temperature - 400.0),
        2e-8,
    );
    close(report.material.enthalpy_residual_j, 0.0, 2e-8);
    close(report.material.energy.energy_residual_j, 0.0, 2e-8);
    assert!(report.material.energy.work_by_gas_j > 0.0);
    assert!(report.chemistry.balance().entropy_production_j_k >= -1e-9);
    assert_eq!(
        parcel.state().state().pressure,
        initial.state().state().pressure
    );
}

fn calorimeter(bank: &ThermochemicalGasModel) -> EquilibriumEnthalpyPhaseCurve {
    let identity = bank.thermodynamics().components()[0].regions()[0]
        .card()
        .content_hash();
    let knot = |h, temperature_k, liquid_mass_fraction, bulk_density_kg_m3| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k,
        liquid_mass_fraction,
        bulk_density_kg_m3,
    };
    EquilibriumEnthalpyPhaseCurve::try_new(
        identity,
        vec![
            knot(0.0, 200.0, 0.0, 1000.0),
            knot(10_000.0, 300.0, 0.0, 900.0),
            knot(30_000.0, 300.0, 1.0, 700.0),
            knot(50_000.0, 500.0, 1.0, 600.0),
        ],
    )
    .unwrap()
}

#[test]
fn g1_finite_body_latent_heat_and_both_volume_works_close_together() {
    let (bank, mechanism, mut parcel) = fixture(0.0, 2000.0);
    let chart = calorimeter(&bank);
    let body = FiniteHeatBody::new(
        &chart,
        Mass::new(0.001),
        Pressure::new(2e5),
        15_000.0,
        ThermalConductance::new(1.0),
    )
    .unwrap();
    let initial = parcel.clone();
    let mut control = controls();
    control.chemistry.method = IsochoricKineticsMethod::Sdirk2;
    let mut reactor = MaterialIsobaricReactor::new(bank, mechanism, control, || false).unwrap();
    let step = reactor
        .advance_with_finite_body(&mut parcel, Time::new(0.05), body, || false)
        .unwrap();
    let capacity = initial.mass_kg() / 0.01 * 4.0 * R;
    let expected_t = 300.0 + 100.0 * (-0.05 / capacity).exp();
    close(parcel.state().state().temperature, expected_t, 3e-4);
    assert_eq!(step.body.state().temperature_k(), 300.0);
    assert!(step.body.state().liquid_mass_fraction() > body.state().liquid_mass_fraction());
    assert!(step.body.state().liquid_mass_fraction() < 1.0);
    assert!(std::ptr::eq(step.body.curve(), &chart));
    assert!(step.gas.material.energy.work_by_gas_j < 0.0);
    assert!(step.body_report.boundary_work_j > 0.0);
    close(
        step.body_report.enthalpy_change_j,
        -step.gas.chemistry.balance().heat_input_j,
        2e-8,
    );
    let actual_gas_u = initial.mass_kg()
        * (parcel.state().specific_internal_energy_j_kg()
            - initial.state().specific_internal_energy_j_kg());
    let actual_gas_work =
        initial.state().state().pressure * (parcel.volume_m3() - initial.volume_m3());
    close(
        actual_gas_u
            + step.body_report.internal_energy_change_j
            + actual_gas_work
            + step.body_report.boundary_work_j,
        0.0,
        3e-8,
    );
    close(step.combined_enthalpy_residual_j, 0.0, 3e-8);
    close(step.combined_energy_residual_j, 0.0, 3e-8);
    close(
        step.gas.chemistry.balance().entropy_transfer_j_k,
        -step.body_entropy_change_j_k,
        1e-14,
    );
}

#[test]
fn g0_zero_time_source_identity_and_endpoint_refusals_preserve_material_state() {
    let (bank, mechanism, mut parcel) = fixture(0.0, 2000.0);
    let initial = parcel.clone();
    let mut reactor =
        MaterialIsobaricReactor::new(bank.clone(), mechanism.clone(), controls(), || false)
            .unwrap();
    let step = reactor
        .advance_adiabatic(&mut parcel, Time::new(0.0), || false)
        .unwrap();
    assert_eq!(step.material.enthalpy_change_j, 0.0);
    assert_eq!(parcel, initial);
    let wrong = ThermochemicalGasModel::new(
        bank.thermodynamics().clone(),
        vec![transport("A", 2.0, 2000.0), transport("B", 5.0, 2000.0)],
        None,
    )
    .unwrap();
    let mut wrong = MaterialIsobaricReactor::new(wrong, mechanism, controls(), || false).unwrap();
    assert!(matches!(
        wrong.advance_adiabatic(&mut parcel, Time::new(0.1), || false),
        Err(MaterialKineticsError::SourceMismatch(_))
    ));
    assert_eq!(parcel, initial);
    let (bank, mechanism, mut hot) = fixture(-24_000.0, 450.0);
    let initial = hot.clone();
    let mut reactor = MaterialIsobaricReactor::new(bank, mechanism, controls(), || false).unwrap();
    assert!(
        reactor
            .advance_adiabatic(&mut hot, Time::new(0.05), || false)
            .is_err()
    );
    assert_eq!(hot, initial);
}

#[test]
fn g4_cancelled_joint_heat_step_preserves_both_endpoints_and_prepared_reuse() {
    let (bank, mechanism, original) = fixture(0.0, 2000.0);
    let chart = calorimeter(&bank);
    let body = FiniteHeatBody::new(
        &chart,
        Mass::new(0.001),
        Pressure::new(2e5),
        15_000.0,
        ThermalConductance::new(1.0),
    )
    .unwrap();
    let mut reactor = MaterialIsobaricReactor::new(bank, mechanism, controls(), || false).unwrap();
    let mut accepted = original.clone();
    let mut polls = 0;
    let expected_body = reactor
        .advance_with_finite_body(&mut accepted, Time::new(0.01), body, || {
            polls += 1;
            false
        })
        .unwrap()
        .body;
    assert!(polls > 10);
    for stop in [1, polls / 2, polls] {
        let mut parcel = original.clone();
        let mut count = 0;
        assert_eq!(
            reactor
                .advance_with_finite_body(&mut parcel, Time::new(0.01), body, || {
                    count += 1;
                    count >= stop
                },)
                .unwrap_err(),
            MaterialKineticsError::Cancelled
        );
        assert_eq!(parcel, original);
        assert_eq!(body.state().specific_enthalpy_j_kg(), 15_000.0);
        let replay = reactor
            .advance_with_finite_body(&mut parcel, Time::new(0.01), body, || false)
            .unwrap();
        assert_eq!(replay.body, expected_body);
        assert_eq!(parcel, accepted);
    }
}
