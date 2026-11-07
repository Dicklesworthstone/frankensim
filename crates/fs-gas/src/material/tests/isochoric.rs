use super::*;
use crate::chemistry::IsochoricKineticsMethod;
use crate::chemistry::heat::{HeatBath, ThermalConductance, material::FiniteHeatBody};
use crate::material::isochoric::{MaterialIsochoricReactor, MaterialIsochoricStep};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase};
use fs_qty::Mass;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

// Both isomers have Cp/R=4+0.004T, hence Cv/R=3+0.004T.
// Their formation difference releases binding joules per mole of A -> B.
fn variable_fixture(
    binding: f64,
    maximum_transport_temperature: f64,
) -> (ThermochemicalGasModel, IdealGasKineticsV1, ClosedGasParcel) {
    let (base, original, initial) = fixture(-binding, maximum_transport_temperature);
    let models: Vec<_> = original
        .thermodynamics()
        .models()
        .iter()
        .map(|model| {
            let mut card = model.regions()[0].card().clone();
            card.parameters.get_mut("a3").unwrap().value = 0.004;
            Nasa9StandardStateModelV1::new(
                model.species().clone(),
                model.molar_mass(),
                model.convention().clone(),
                vec![Nasa9RegionV1::from_card(card).unwrap()],
            )
            .unwrap()
        })
        .collect();
    let source = IdealGasEquilibriumV1::new(
        original.thermodynamics().elemental_matrix().clone(),
        original.thermodynamics().charge_vector().clone(),
        models.clone(),
    )
    .unwrap();
    let mechanism = IdealGasKineticsV1::new(
        source,
        "nonlinear-Cv-rigid-fixture",
        1,
        vec![original.reactions()[0].clone().irreversible()],
        0.0,
    )
    .unwrap();
    let bank = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            models,
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
        )
        .unwrap(),
        base.transport_materials().to_vec(),
        None,
    )
    .unwrap();
    let parcel = ClosedGasParcel::new(
        bank.with_mole_amounts(&[Amount::new(1.0), Amount::new(0.0)])
            .unwrap(),
        initial.point().clone(),
        initial.volume_m3(),
    )
    .unwrap();
    (bank, mechanism, parcel)
}

fn nonlinear_body(bank: &ThermochemicalGasModel) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity_exact(
        bank.thermodynamics().components()[0].regions()[0]
            .card()
            .content_hash(),
        SolidLiquidPhase::Solid,
        0.0,
        &[
            HeatCapacityKnot {
                temperature_k: 200.0,
                specific_heat_capacity_j_kg_k: 100.0,
                bulk_density_kg_m3: 1000.0,
            },
            HeatCapacityKnot {
                temperature_k: 1200.0,
                specific_heat_capacity_j_kg_k: 300.0,
                bulk_density_kg_m3: 800.0,
            },
        ],
    )
    .unwrap()
}

fn body(chart: &EquilibriumEnthalpyPhaseCurve, temperature: f64) -> FiniteHeatBody<'_> {
    let delta = temperature - 200.0;
    FiniteHeatBody::new(
        chart,
        Mass::new(0.001),
        Pressure::new(2e5),
        100.0 * delta + 0.1 * delta * delta,
        ThermalConductance::new(2.0),
    )
    .unwrap()
}

// Independent original-capacity energy equations, with exact first-order
// chemical release and scalar RK4. No production caloric inverse is used.
fn reference(amount: f64, binding: f64, boundary: f64, finite: bool) -> [f64; 2] {
    let rate = |time: f64, t: [f64; 2]| {
        let heat = 2.0 * (t[1] - t[0]);
        [
            (amount * binding * 4.0 * (-4.0 * time).exp() + heat)
                / (amount * R * (3.0 + 0.004 * t[0])),
            if finite {
                -heat / (0.001 * (60.0 + 0.2 * t[1]))
            } else {
                0.0
            },
        ]
    };
    let mut temperature = [400.0, boundary];
    let dt = 0.1 / 4096.0;
    for step in 0..4096 {
        let time = step as f64 * dt;
        let shift = |rate: [f64; 2], fraction: f64| {
            std::array::from_fn(|i| temperature[i] + fraction * dt * rate[i])
        };
        let a = rate(time, temperature);
        let b = rate(time + 0.5 * dt, shift(a, 0.5));
        let c = rate(time + 0.5 * dt, shift(b, 0.5));
        let d = rate(time + dt, shift(c, 1.0));
        for i in 0..2 {
            temperature[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    temperature
}

fn check_material(
    initial: &ClosedGasParcel,
    parcel: &ClosedGasParcel,
    report: &MaterialIsochoricStep<'_>,
    binding: f64,
) {
    let state = parcel.state().state();
    let fraction = parcel
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions()[1];
    let amount = initial.mass_kg() / 0.01;
    let delta = state.temperature - 400.0;
    let energy = amount * R * delta * (3.0 + 0.002 * (state.temperature + 400.0))
        - amount * binding * fraction;
    close(report.material.stored_energy_change_j, energy, 3e-8);
    close(report.material.supplied_heat_j, energy, 3e-8);
    assert_eq!(
        report.material.supplied_heat_j,
        report.chemistry.balance().heat_input_j
    );
    assert_eq!(report.material.work_by_gas_j, 0.0);
    assert_eq!(parcel.mass_kg().to_bits(), initial.mass_kg().to_bits());
    assert_eq!(parcel.volume_m3().to_bits(), initial.volume_m3().to_bits());
    close(
        state.pressure / initial.state().state().pressure,
        state.temperature / 400.0,
        1e-10,
    );
    assert_eq!(parcel.point().axes()["sample_id"], 7.0);
    assert!(state.thermal_conductivity > initial.state().state().thermal_conductivity);
    state.validate_coherence().unwrap();
    close(report.balance.freeze_internal_energy_residual_j, 0.0, 3e-8);
    assert!(report.chemistry.balance().entropy_production_j_k >= -1e-9);
}

#[test]
fn g1_g3_bath_and_original_variable_cv_reaction_publish_one_material_endpoint() {
    for method in [
        IsochoricKineticsMethod::BackwardEuler,
        IsochoricKineticsMethod::Sdirk2,
    ] {
        for (binding, boundary) in [(0.0, 700.0), (0.0, 300.0), (24_000.0, 300.0)] {
            let (bank, mechanism, mut parcel) = variable_fixture(binding, 2000.0);
            let initial = parcel.clone();
            let mut control = controls();
            control.chemistry.method = method;
            control.chemistry.relative_tolerance = 2e-7;
            let mut reactor =
                MaterialIsochoricReactor::new(bank.clone(), mechanism, control, || false).unwrap();
            let bytes = reactor.workspace_bytes();
            let bath =
                HeatBath::new(ThermalConductance::new(2.0), Temperature::new(boundary)).unwrap();
            let report = reactor
                .advance_with_heat_bath(&mut parcel, Time::new(0.1), bath, || false)
                .unwrap();
            let expected = reference(initial.mass_kg() / 0.01, binding, boundary, false);
            let tolerance = if method == IsochoricKineticsMethod::BackwardEuler {
                0.12
            } else {
                0.002
            };
            close(parcel.state().state().temperature, expected[0], tolerance);
            check_material(&initial, &parcel, &report, binding);
            assert_eq!(
                report.material.supplied_heat_j.is_sign_positive(),
                boundary > 400.0
            );
            assert_eq!(
                parcel.model().transport_materials(),
                bank.transport_materials()
            );
            assert_eq!(
                parcel.model().thermodynamics().components(),
                bank.thermodynamics().components()
            );
            assert_eq!(reactor.workspace_bytes(), bytes);
            let mut replay = initial;
            reactor
                .advance_with_heat_bath(&mut replay, Time::new(0.1), bath, || false)
                .unwrap();
            assert_eq!(replay, parcel);
        }
    }
}

#[test]
fn g1_finite_source_cp_and_rigid_cv_close_heat_body_work_and_material_response() {
    for method in [
        IsochoricKineticsMethod::BackwardEuler,
        IsochoricKineticsMethod::Sdirk2,
    ] {
        for (binding, boundary) in [(0.0, 800.0), (0.0, 300.0), (24_000.0, 300.0)] {
            let (bank, mechanism, mut parcel) = variable_fixture(binding, 2000.0);
            let initial = parcel.clone();
            let chart = nonlinear_body(&bank);
            assert_eq!(chart.knots().len(), 2);
            let original_body = body(&chart, boundary);
            let mut control = controls();
            control.chemistry.method = method;
            control.chemistry.relative_tolerance = 2e-7;
            let mut reactor =
                MaterialIsochoricReactor::new(bank, mechanism, control, || false).unwrap();
            let report = reactor
                .advance_with_finite_body(&mut parcel, Time::new(0.1), original_body, || false)
                .unwrap();
            let expected = reference(initial.mass_kg() / 0.01, binding, boundary, true);
            let tolerance = if method == IsochoricKineticsMethod::BackwardEuler {
                0.12
            } else {
                0.002
            };
            close(parcel.state().state().temperature, expected[0], tolerance);
            close(report.body.state().temperature_k(), expected[1], tolerance);
            check_material(&initial, &parcel, &report.gas, binding);
            assert!(std::ptr::eq(report.body.curve(), &chart));
            assert!(
                report.body_report.boundary_work_j * report.body_report.enthalpy_change_j > 0.0
            );
            close(report.combined_heat_residual_j, 0.0, 3e-8);
            close(report.combined_energy_residual_j, 0.0, 3e-8);
            close(
                report.gas.material.stored_energy_change_j
                    + report.body_report.internal_energy_change_j
                    + report.body_report.boundary_work_j,
                0.0,
                3e-8,
            );
            assert_eq!(
                report.gas.chemistry.balance().entropy_transfer_j_k,
                -report.body_entropy_change_j_k
            );
        }
    }
}

#[test]
fn g0_g4_zero_time_late_transport_refusal_and_cancellation_preserve_joint_state() {
    let (bank, mechanism, original) = variable_fixture(0.0, 2000.0);
    let chart = nonlinear_body(&bank);
    let original_body = body(&chart, 800.0);
    let mut reactor =
        MaterialIsochoricReactor::new(bank.clone(), mechanism.clone(), controls(), || false)
            .unwrap();
    let mut parcel = original.clone();
    reactor
        .advance_adiabatic(&mut parcel, Time::new(0.0), || false)
        .unwrap();
    let zero = reactor
        .advance_with_finite_body(&mut parcel, Time::new(0.0), original_body, || false)
        .unwrap();
    assert_eq!(zero.body, original_body);
    assert_eq!(zero.gas.material.stored_energy_change_j, 0.0);
    assert_eq!(parcel, original);
    let mut polls = 0;
    let expected_body = reactor
        .advance_with_finite_body(&mut parcel, Time::new(0.1), original_body, || {
            polls += 1;
            false
        })
        .unwrap()
        .body;
    let accepted = parcel.clone();
    for stop in [1, polls / 2, polls] {
        parcel = original.clone();
        let mut count = 0;
        assert_eq!(
            reactor
                .advance_with_finite_body(&mut parcel, Time::new(0.1), original_body, || {
                    count += 1;
                    count >= stop
                })
                .unwrap_err(),
            MaterialKineticsError::Cancelled
        );
        assert_eq!(parcel, original);
        let replay = reactor
            .advance_with_finite_body(&mut parcel, Time::new(0.1), original_body, || false)
            .unwrap();
        assert_eq!(replay.body, expected_body);
        assert_eq!(parcel, accepted);
    }
    let wrong = ThermochemicalGasModel::new(
        bank.thermodynamics().clone(),
        vec![transport("A", 2.0, 2000.0), transport("B", 5.0, 2000.0)],
        None,
    )
    .unwrap();
    let mut wrong = MaterialIsochoricReactor::new(wrong, mechanism, controls(), || false).unwrap();
    parcel = original.clone();
    assert!(matches!(
        wrong.advance_adiabatic(&mut parcel, Time::new(0.0), || false),
        Err(MaterialKineticsError::SourceMismatch(_))
    ));
    assert_eq!(parcel, original);
    let (bank, mechanism, mut limited) = variable_fixture(0.0, 450.0);
    let before = limited.clone();
    let mut limited_reactor =
        MaterialIsochoricReactor::new(bank, mechanism, controls(), || false).unwrap();
    assert!(matches!(
        limited_reactor.advance_with_finite_body(
            &mut limited,
            Time::new(0.1),
            original_body,
            || false
        ),
        Err(MaterialKineticsError::Energy(_)) | Err(MaterialKineticsError::Material(_))
    ));
    assert_eq!(limited, before);
    limited_reactor
        .advance_with_finite_body(&mut limited, Time::new(0.001), original_body, || false)
        .unwrap();
    assert!(limited.state().state().temperature > 400.0);
}
