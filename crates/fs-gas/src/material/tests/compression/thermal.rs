//! Independent work/heat trajectories on the common source-bearing gas parcel.
mod adaptive;

use super::*;
use crate::chemistry::IsochoricKineticsMethod;
use crate::chemistry::heat::{HeatBath, ThermalConductance, material::FiniteHeatBody};
use crate::material::compression::thermal::{
    MaterialThermalCompressionReactor, MaterialThermalCompressionSolve,
};
use crate::material::isochoric::MaterialIsochoricReactor;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase,
};
use fs_qty::Mass;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && expected.is_finite() && (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn kinetics() -> MaterialKineticsSolve {
    let mut control = controls();
    control.chemistry.method = IsochoricKineticsMethod::Sdirk2;
    control.chemistry.relative_tolerance = 2e-9;
    control
}

fn path(substeps: usize) -> MaterialThermalCompressionSolve {
    MaterialThermalCompressionSolve {
        volume: compression_controls(substeps).volume,
        substeps,
    }
}

// Constant Cv and irreversible A -> B, with prescribed log-volume rate and
// finite conductance: T' = (a-b)T + b*Tb + c*exp(-4t).
// Return T, compression work, and bath heat from an independent closed form.
fn bath_reference(final_volume: f64, duration: f64, bath: f64, g: f64, release: f64) -> [f64; 3] {
    let capacity = 0.01 * 3.0 * R;
    let a = (0.01 / final_volume).ln() / (3.0 * duration);
    let b = g / capacity;
    let lambda = a - b;
    let c = 4.0 * release / (3.0 * R);
    let exponential = (lambda * duration).exp();
    let chemical = (-4.0 * duration).exp();
    let integral_exponential = (lambda * duration).exp_m1() / lambda;
    let temperature = 400.0 * exponential
        + b * bath * integral_exponential
        + c * (exponential - chemical) / (lambda + 4.0);
    let integral_temperature = 400.0 * integral_exponential
        + b * bath * (integral_exponential - duration) / lambda
        + c * (integral_exponential - (1.0 - chemical) / 4.0) / (lambda + 4.0);
    [
        temperature,
        capacity * a * integral_temperature,
        g * (bath * duration - integral_temperature),
    ]
}

#[test]
fn reacting_bath_compression_and_expansion_refine_temperature_work_and_signed_heat() {
    for (final_volume, bath_temperature) in [(0.005, 300.0), (0.02, 700.0)] {
        let duration = 0.15;
        let release = 12_000.0;
        let expected = bath_reference(final_volume, duration, bath_temperature, 0.2, release);
        let (bank, mechanism, initial) = fixture(-release, 2000.0);
        let mut reactor =
            MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
        let workspace = reactor.workspace_bytes();
        let bath = HeatBath::new(
            ThermalConductance::new(0.2),
            Temperature::new(bath_temperature),
        )
        .unwrap();
        let mut temperature_errors = Vec::new();
        let mut work_errors = Vec::new();
        for substeps in [4, 8, 16] {
            let mut parcel = initial.clone();
            let report = reactor
                .advance_with_heat_bath(
                    &mut parcel,
                    Volume::new(final_volume),
                    Time::new(duration),
                    bath,
                    path(substeps),
                    || false,
                )
                .unwrap();
            temperature_errors.push((parcel.state().state().temperature - expected[0]).abs());
            work_errors.push((report.boundary_work_on_gas_j - expected[1]).abs());
            let actual_u = parcel.mass_kg()
                * (parcel.state().specific_internal_energy_j_kg()
                    - initial.state().specific_internal_energy_j_kg());
            close(report.internal_energy_change_j, actual_u, 2e-6);
            close(
                actual_u,
                report.boundary_work_on_gas_j + report.heat_input_j + report.energy_residual_j,
                2e-6,
            );
            close(report.energy_residual_j, 0.0, 2e-7);
            close(
                report.thermal_partner_entropy_change_j_k,
                -report.heat_input_j / bath_temperature,
                2e-10,
            );
            assert!(report.entropy_production_j_k >= -1e-7);
            assert_eq!(report.substeps, substeps);
            assert!(report.chemical_accepted_steps >= substeps);
            assert_eq!(parcel.volume_m3().to_bits(), final_volume.to_bits());
            assert_eq!(parcel.mass_kg().to_bits(), initial.mass_kg().to_bits());
            assert_eq!(report.elapsed_time.value(), duration);
            assert_eq!(
                report.heat_input_j.is_sign_positive(),
                bath_temperature > 400.0
            );
            assert_eq!(
                report.boundary_work_on_gas_j.is_sign_positive(),
                final_volume < 0.01
            );
            close(
                parcel
                    .model()
                    .thermodynamics()
                    .mole_composition()
                    .fractions()[1],
                1.0 - (-4.0 * duration).exp(),
                5e-7,
            );
            assert_eq!(
                parcel.model().transport_materials(),
                reactor.bank().transport_materials()
            );
            parcel.state().state().validate_coherence().unwrap();
            if substeps == 16 {
                close(parcel.state().state().temperature, expected[0], 0.013);
                close(report.boundary_work_on_gas_j, expected[1], 0.005);
                close(report.heat_input_j, expected[2], 0.002);
            }
        }
        for errors in [&temperature_errors, &work_errors] {
            for pair in errors.windows(2) {
                assert!(
                    pair[0] / pair[1] > 3.8,
                    "second-order refinement: {errors:?}"
                );
            }
        }
        assert_eq!(reactor.workspace_bytes(), workspace);
    }
}

// Both species retain Cp/R = 4 + 0.004*T; formation(B) = -release.
fn variable_fixture(release: f64) -> (ThermochemicalGasModel, IdealGasKineticsV1, ClosedGasParcel) {
    let (base, original, initial) = fixture(-release, 2000.0);
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
    let mechanism = IdealGasKineticsV1::new(
        IdealGasEquilibriumV1::new(
            original.thermodynamics().elemental_matrix().clone(),
            original.thermodynamics().charge_vector().clone(),
            models.clone(),
        )
        .unwrap(),
        "variable-Cv-compression-fixture",
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

fn body_curve(bank: &ThermochemicalGasModel) -> EquilibriumEnthalpyPhaseCurve {
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

fn finite_body(
    curve: &EquilibriumEnthalpyPhaseCurve,
    temperature: f64,
    g: f64,
) -> FiniteHeatBody<'_> {
    let delta = temperature - 200.0;
    FiniteHeatBody::new(
        curve,
        Mass::new(0.001),
        Pressure::new(2e5),
        100.0 * delta + 0.1 * delta * delta,
        ThermalConductance::new(g),
    )
    .unwrap()
}

// Independent variable-Cv/Cp first laws. No production thermodynamic inverse,
// chemistry step or volume map participates in this RK4 reference.
fn finite_body_reference(final_volume: f64, boundary: f64) -> [f64; 4] {
    let duration = 0.1;
    let compression_rate = (0.01 / final_volume).ln() / duration;
    let rate = |time: f64, state: [f64; 4]| {
        let heat = 2.0 * (state[1] - state[0]);
        let work = 0.01 * R * state[0] * compression_rate;
        [
            (work + 0.01 * 24_000.0 * 4.0 * (-4.0 * time).exp() + heat)
                / (0.01 * R * (3.0 + 0.004 * state[0])),
            -heat / (0.001 * (60.0 + 0.2 * state[1])),
            work,
            heat,
        ]
    };
    let mut state = [400.0, boundary, 0.0, 0.0];
    let dt = duration / 4096.0;
    for index in 0..4096 {
        let time = index as f64 * dt;
        let shift = |rate: [f64; 4], fraction: f64| {
            std::array::from_fn(|i| state[i] + fraction * dt * rate[i])
        };
        let a = rate(time, state);
        let b = rate(time + 0.5 * dt, shift(a, 0.5));
        let c = rate(time + 0.5 * dt, shift(b, 0.5));
        let d = rate(time + dt, shift(c, 1.0));
        for i in 0..4 {
            state[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}

#[test]
fn finite_body_retains_variable_source_capacity_density_and_both_work_accounts() {
    for (final_volume, boundary) in [(0.005, 800.0), (0.02, 300.0)] {
        let (bank, mechanism, mut parcel) = variable_fixture(24_000.0);
        let initial = parcel.clone();
        let curve = body_curve(&bank);
        let body = finite_body(&curve, boundary, 2.0);
        let expected = finite_body_reference(final_volume, boundary);
        let mut reactor =
            MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
        let report = reactor
            .advance_with_finite_body(
                &mut parcel,
                Volume::new(final_volume),
                Time::new(0.1),
                body,
                path(64),
                || false,
            )
            .unwrap();
        close(parcel.state().state().temperature, expected[0], 0.006);
        close(report.body.state().temperature_k(), expected[1], 0.006);
        close(report.gas.boundary_work_on_gas_j, expected[2], 0.003);
        close(report.gas.heat_input_j, expected[3], 0.003);
        let gas_u = parcel.mass_kg()
            * (parcel.state().specific_internal_energy_j_kg()
                - initial.state().specific_internal_energy_j_kg());
        close(gas_u, report.gas.internal_energy_change_j, 3e-6);
        close(
            gas_u + report.body_internal_energy_change_j + report.body_boundary_work_j,
            report.gas.boundary_work_on_gas_j,
            3e-6,
        );
        close(
            report.body_enthalpy_change_j,
            -report.gas.heat_input_j,
            3e-6,
        );
        close(report.combined_heat_residual_j, 0.0, 3e-6);
        close(report.combined_energy_residual_j, 0.0, 3e-6);
        close(
            report.body_entropy_change_j_k,
            body.mass().value()
                * curve
                    .specific_entropy_change_j_kg_k(
                        body.state().specific_enthalpy_j_kg(),
                        report.body.state().specific_enthalpy_j_kg(),
                    )
                    .unwrap(),
            3e-9,
        );
        assert!(report.body_boundary_work_j * report.body_enthalpy_change_j > 0.0);
        assert_eq!(report.body.mass(), body.mass());
        assert_eq!(report.body.pressure(), body.pressure());
        assert!(std::ptr::eq(report.body.curve(), &curve));
        assert_eq!(parcel.volume_m3(), final_volume);
        parcel.state().state().validate_coherence().unwrap();
    }
}

#[test]
fn compression_can_melt_a_finite_body_without_relabelling_latent_heat_as_gas_work() {
    let (bank, mechanism, mut parcel) = fixture(0.0, 2000.0);
    let initial = parcel.clone();
    let knot = |h, t, fraction, rho| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: fraction,
        bulk_density_kg_m3: rho,
    };
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        bank.thermodynamics().components()[0].regions()[0]
            .card()
            .content_hash(),
        vec![
            knot(0.0, 200.0, 0.0, 1100.0),
            knot(20_000.0, 400.0, 0.0, 1000.0),
            knot(120_000.0, 400.0, 1.0, 800.0),
            knot(200_000.0, 1200.0, 1.0, 700.0),
        ],
    )
    .unwrap();
    let body = FiniteHeatBody::new(
        &curve,
        Mass::new(0.001),
        Pressure::new(2e5),
        20_000.0,
        ThermalConductance::new(1.0),
    )
    .unwrap();
    let expected = bath_reference(0.005, 0.1, 400.0, 1.0, 0.0);
    let mut reactor =
        MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
    let report = reactor
        .advance_with_finite_body(
            &mut parcel,
            Volume::new(0.005),
            Time::new(0.1),
            body,
            path(64),
            || false,
        )
        .unwrap();
    close(parcel.state().state().temperature, expected[0], 0.004);
    close(report.gas.boundary_work_on_gas_j, expected[1], 0.001);
    close(report.gas.heat_input_j, expected[2], 0.001);
    assert_eq!(report.body.state().temperature_k(), 400.0);
    let fraction = report.body.state().liquid_mass_fraction();
    assert!(fraction > 0.0 && fraction < 1.0);
    close(
        fraction,
        -expected[2] / (body.mass().value() * 100_000.0),
        1e-5,
    );
    assert!(report.body_boundary_work_j > 0.0);
    close(
        report.body_internal_energy_change_j + report.body_boundary_work_j,
        -report.gas.heat_input_j,
        3e-6,
    );
    let gas_u = parcel.mass_kg()
        * (parcel.state().specific_internal_energy_j_kg()
            - initial.state().specific_internal_energy_j_kg());
    close(
        gas_u + report.body_internal_energy_change_j + report.body_boundary_work_j,
        report.gas.boundary_work_on_gas_j,
        3e-6,
    );
}

#[test]
fn fixed_volume_matches_the_native_reactor_and_zero_duration_preserves_both_states() {
    let (bank, mechanism, initial) = fixture(-12_000.0, 2000.0);
    let curve = body_curve(&bank);
    let body = finite_body(&curve, 800.0, 2.0);
    let bath = HeatBath::new(ThermalConductance::new(0.2), Temperature::new(300.0)).unwrap();
    let mut reactor =
        MaterialThermalCompressionReactor::new(bank.clone(), mechanism.clone(), kinetics(), || {
            false
        })
        .unwrap();
    let mut native = MaterialIsochoricReactor::new(bank, mechanism, kinetics(), || false).unwrap();
    let mut fixed = initial.clone();
    let report = reactor
        .advance_with_heat_bath(
            &mut fixed,
            Volume::new(initial.volume_m3()),
            Time::new(0.05),
            bath,
            path(4096),
            || false,
        )
        .unwrap();
    let mut reference = initial.clone();
    let original = native
        .advance_with_heat_bath(&mut reference, Time::new(0.05), bath, || false)
        .unwrap();
    assert_eq!(fixed, reference);
    assert_eq!(report.substeps, 1);
    assert_eq!(
        report.chemical_accepted_steps,
        original.chemistry.balance().accepted_steps
    );
    assert_eq!(report.heat_input_j, original.material.supplied_heat_j);
    assert_eq!(report.boundary_work_on_gas_j, 0.0);
    let mut fixed = initial.clone();
    let report = reactor
        .advance_with_finite_body(
            &mut fixed,
            Volume::new(initial.volume_m3()),
            Time::new(0.05),
            body,
            path(4096),
            || false,
        )
        .unwrap();
    let mut reference = initial.clone();
    let original = native
        .advance_with_finite_body(&mut reference, Time::new(0.05), body, || false)
        .unwrap();
    assert_eq!(fixed, reference);
    assert_eq!(report.body, original.body);
    assert_eq!(report.gas.substeps, 1);
    assert_eq!(
        report.body_enthalpy_change_j,
        original.body_report.enthalpy_change_j
    );
    // Insulation recovers the original adiabatic compression owner with the
    // same implicit chemistry controls, without adding a second heat source.
    let mut insulated = initial.clone();
    let zero_bath = HeatBath::new(ThermalConductance::new(0.0), Temperature::new(800.0)).unwrap();
    let insulated_report = reactor
        .advance_with_heat_bath(
            &mut insulated,
            Volume::new(0.007),
            Time::new(0.05),
            zero_bath,
            path(4),
            || false,
        )
        .unwrap();
    let mut adiabatic = initial.clone();
    let original = advance_closed_parcel_compression(
        &mut adiabatic,
        reactor.bank(),
        reactor.mechanism(),
        Volume::new(0.007),
        Time::new(0.05),
        MaterialCompressionSolve {
            kinetics: kinetics(),
            volume: path(4).volume,
            substeps: 4,
        },
        || false,
    )
    .unwrap();
    close(
        insulated.state().state().temperature,
        adiabatic.state().state().temperature,
        1e-6,
    );
    close(
        insulated_report.boundary_work_on_gas_j,
        original.boundary_work_on_gas_j,
        1e-6,
    );
    assert_eq!(insulated_report.heat_input_j, 0.0);

    let mut unchanged = initial.clone();
    let zero = reactor
        .advance_with_finite_body(
            &mut unchanged,
            Volume::new(initial.volume_m3()),
            Time::new(0.0),
            body,
            path(3),
            || false,
        )
        .unwrap();
    assert_eq!(unchanged, initial);
    assert_eq!(zero.body, body);
    assert_eq!(zero.gas.heat_input_j, 0.0);
    assert_eq!(zero.gas.internal_energy_change_j, 0.0);
    assert_eq!(zero.body_enthalpy_change_j, 0.0);
    for (volume, time, steps) in [(0.005, 0.0, 3), (0.005, 0.1, 0), (0.005, 0.1, 4097)] {
        assert!(matches!(
            reactor.advance_with_heat_bath(
                &mut unchanged,
                Volume::new(volume),
                Time::new(time),
                bath,
                path(steps),
                || false,
            ),
            Err(MaterialCompressionError::Invalid(_))
        ));
        assert_eq!(unchanged, initial);
    }
}

#[test]
fn finite_body_strokes_roll_back_late_cancellation_and_source_failure_then_replay() {
    let (bank, mechanism, initial) = variable_fixture(24_000.0);
    let curve = body_curve(&bank);
    let body = finite_body(&curve, 800.0, 2.0);
    let mut reactor =
        MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
    let mut accepted = initial.clone();
    let mut polls = 0;
    let expected = reactor
        .advance_with_finite_body(
            &mut accepted,
            Volume::new(0.005),
            Time::new(0.1),
            body,
            path(8),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut parcel = initial.clone();
        let mut calls = 0;
        assert_eq!(
            reactor
                .advance_with_finite_body(
                    &mut parcel,
                    Volume::new(0.005),
                    Time::new(0.1),
                    body,
                    path(8),
                    || {
                        calls += 1;
                        calls >= stop
                    },
                )
                .unwrap_err(),
            MaterialCompressionError::Cancelled,
        );
        assert_eq!(parcel, initial);
        let replay = reactor
            .advance_with_finite_body(
                &mut parcel,
                Volume::new(0.005),
                Time::new(0.1),
                body,
                path(8),
                || false,
            )
            .unwrap();
        assert_eq!(parcel, accepted);
        assert_eq!(replay, expected);
    }

    let (bank, mechanism, initial) = fixture(0.0, 450.0);
    let curve = body_curve(&bank);
    let body = finite_body(&curve, 300.0, 0.05);
    let mut reactor =
        MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
    let first_volume = 0.01 * (-(2.0_f64).ln() / 8.0).exp();
    let mut prefix = initial.clone();
    reactor
        .advance_with_finite_body(
            &mut prefix,
            Volume::new(first_volume),
            Time::new(0.1 / 8.0),
            body,
            path(1),
            || false,
        )
        .unwrap();
    assert!(prefix.state().state().temperature > 400.0);
    let mut failed = initial.clone();
    assert!(matches!(
        reactor.advance_with_finite_body(
            &mut failed,
            Volume::new(0.005),
            Time::new(0.1),
            body,
            path(8),
            || false,
        ),
        Err(MaterialCompressionError::Volume(_)) | Err(MaterialCompressionError::Kinetics(_))
    ));
    assert_eq!(failed, initial);
    let retry = reactor
        .advance_with_finite_body(
            &mut failed,
            Volume::new(first_volume),
            Time::new(0.1 / 8.0),
            body,
            path(1),
            || false,
        )
        .unwrap();
    assert_eq!(failed, prefix);
    assert!(retry.body.state().temperature_k() > body.state().temperature_k());
}
