//! Conservative real-gas midpoint behavior and analytic constitutive tangents.
use super::*;
use crate::thermochemical_gas::tests::model_in_domains;

fn solve_controls() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 210.0,
        maximum_temperature_k: 450.0,
        energy_tolerance_j_kg: 1.0e-7,
        maximum_iterations: 128,
    }
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual / expected - 1.0).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, relative tolerance={tolerance:e}"
    );
}

#[test]
fn real_gas_midpoint_matches_independent_dense_heat_work_entropy_and_derivatives() {
    let gas = parcel();
    let saved = gas.clone();
    // Independent 75-digit Decimal PR p(T,v), analytic departure and NASA Cv
    // integration. Entropy integrates in T,v without the production G/fugacity
    // identity. Values include PR energy's dependence on BOTH densities.
    for (volume, heat, tm, pm, t1, p1, work, stored, entropy, dpdv, dpdq) in [
        (
            0.008,
            2500.0,
            318.350_370_608_637_3,
            5_937_732.629_253_424,
            337.151_945_135_981_43,
            7_126_636.252_254_585,
            -11_875.465_258_506_847,
            14_375.465_258_506_847,
            7.727_422_975_232_885,
            -518_856_750.773_195_8,
            26.192_838_930_502_605,
        ),
        (
            0.012,
            -1500.0,
            286.590_190_866_755_34,
            4_326_347.291_170_364,
            273.457_011_230_247,
            3_769_661.020_966_815,
            8_652.694_582_340_73,
            -10_152.694_582_340_73,
            -5.167_228_839_700_856,
            -278_289_279.776_018_3,
            19.476_461_365_723_275,
        ),
    ] {
        let candidate = gas
            .midpoint_heat_volume(volume, heat, solve_controls(), || false)
            .unwrap();
        let step = candidate.step();
        near(candidate.midpoint().temperature_k(), tm, 2.0e-8);
        near(candidate.midpoint().pressure_pa(), pm, 1.0e-3);
        near(candidate.endpoint().state().temperature_k(), t1, 2.0e-8);
        near(candidate.endpoint().state().pressure_pa(), p1, 1.0e-3);
        near(step.work_by_gas_j, work, 2.0e-5);
        near(step.stored_energy_change_j, stored, 2.0e-5);
        near(step.entropy_change_j_k, entropy, 2.0e-8);
        relative(step.pressure_final_volume_derivative_pa_m3, dpdv, 2.0e-10);
        relative(step.pressure_heat_derivative_pa_j, dpdq, 2.0e-10);
        assert!(
            step.energy_residual_j.abs() <= gas.mass_kg() * solve_controls().energy_tolerance_j_kg
        );
        assert!(step.midpoint_energy_residual_j_kg.abs() <= solve_controls().energy_tolerance_j_kg);
        assert_eq!(candidate.endpoint().mass_kg(), gas.mass_kg());
        assert_eq!(candidate.endpoint().volume_m3(), volume);
        assert_eq!(candidate.endpoint().density_kg_m3(), gas.mass_kg() / volume);
        assert_eq!(candidate.endpoint().model(), gas.model());
        assert_eq!(
            candidate.endpoint().point().axes().get("material-lot"),
            Some(&12.0)
        );
        for state in [candidate.midpoint(), candidate.endpoint().state()] {
            assert_eq!(
                state.transport_state().unwrap().density_kg_m3(),
                state.density_kg_m3()
            );
            assert_eq!(
                state.acoustic_state().unwrap().density_kg_m3(),
                state.density_kg_m3()
            );
        }
        // A temperature average is not the conservative midpoint for this law.
        assert!((tm - 0.5 * (300.0 + t1)).abs() > 0.05);
        assert_eq!(gas, saved);
    }
}

#[test]
fn real_gas_midpoint_discrete_tangents_match_fresh_source_trials() {
    let gas = parcel();
    for volume in [0.008, 0.01, 0.012] {
        let heat = 1000.0;
        let candidate = gas
            .midpoint_heat_volume(volume, heat, solve_controls(), || false)
            .unwrap();
        let pressure = |v, q| {
            gas.midpoint_heat_volume(v, q, solve_controls(), || false)
                .unwrap()
                .midpoint()
                .pressure_pa()
        };
        let dv = 1.0e-7;
        let dq = 0.5;
        relative(
            candidate.step().pressure_final_volume_derivative_pa_m3,
            (pressure(volume + dv, heat) - pressure(volume - dv, heat)) / (2.0 * dv),
            2.0e-6,
        );
        relative(
            candidate.step().pressure_heat_derivative_pa_j,
            (pressure(volume, heat + dq) - pressure(volume, heat - dq)) / (2.0 * dq),
            2.0e-6,
        );
    }
}

#[test]
fn real_gas_midpoint_dilute_limit_matches_discrete_gamma_law() {
    let source = model_in_domains(
        vec![0.75, 0.25],
        0.0,
        [200.0, 2000.0],
        [200.0, 2000.0],
        [0.01, 1.0e7],
    );
    let gas = IsochoricRealGasParcel::new(
        bound(source),
        point(300.0, 0.1),
        0.01,
        PengRobinsonPhaseV1::Vapor,
    )
    .unwrap();
    let r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.022;
    let cv = 2.25 * r;
    for volume in [0.008, 0.012] {
        for q in [-2000.0, 0.0, 5000.0] {
            let candidate = gas
                .midpoint_heat_volume(volume, gas.mass_kg() * q, solve_controls(), || false)
                .unwrap();
            let vm = 0.5 * (gas.volume_m3() + volume);
            let ratio = (volume - gas.volume_m3()) / vm;
            let denominator = cv + 0.5 * ratio * r;
            let tm = (cv * 300.0 + 0.5 * q) / denominator;
            let pm = gas.mass_kg() * r * tm / vm;
            relative(candidate.midpoint().temperature_k(), tm, 2.0e-8);
            relative(candidate.midpoint().pressure_pa(), pm, 2.0e-8);
            relative(
                candidate.endpoint().state().temperature_k(),
                2.0 * tm - 300.0,
                2.0e-8,
            );
            relative(
                candidate.step().work_by_gas_j,
                pm * (volume - gas.volume_m3()),
                2.0e-8,
            );
            relative(
                candidate.step().pressure_final_volume_derivative_pa_m3,
                -pm / (2.0 * vm) * (cv + r) / denominator,
                2.0e-8,
            );
            relative(
                candidate.step().pressure_heat_derivative_pa_j,
                0.5 * r / (vm * denominator),
                2.0e-8,
            );
        }
    }
}

#[test]
fn real_gas_midpoint_preserves_reference_invariance_and_fixed_volume_heating() {
    let gas = parcel();
    let identity = gas
        .midpoint_heat_volume(gas.volume_m3(), 0.0, solve_controls(), || false)
        .unwrap();
    assert_eq!(identity.endpoint(), &gas);
    assert_eq!(identity.midpoint(), gas.state());
    assert_eq!(identity.step().iterations, 0);
    for heat in [-2500.0, 2500.0] {
        let mut direct = gas.clone();
        direct.add_heat(heat, solve_controls(), || false).unwrap();
        let candidate = gas
            .midpoint_heat_volume(gas.volume_m3(), heat, solve_controls(), || false)
            .unwrap();
        near(
            candidate.endpoint().state().temperature_k(),
            direct.state().temperature_k(),
            2.0e-8,
        );
        near(
            candidate.endpoint().state().pressure_pa(),
            direct.state().pressure_pa(),
            1.0e-3,
        );
        assert_eq!(candidate.endpoint().density_kg_m3(), gas.density_kg_m3());
        assert_eq!(candidate.step().work_by_gas_j, 0.0);
    }
    let plain = super::volume::with_references(0.0, 0.0);
    let shifted = super::volume::with_references(1.0e20, 1.0e24);
    for volume in [0.008, 0.012] {
        let a = plain
            .midpoint_heat_volume(volume, 2500.0, solve_controls(), || false)
            .unwrap();
        let b = shifted
            .midpoint_heat_volume(volume, 2500.0, solve_controls(), || false)
            .unwrap();
        assert_eq!(a.step(), b.step());
        assert_eq!(a.midpoint().temperature_k(), b.midpoint().temperature_k());
        assert_eq!(a.midpoint().pressure_pa(), b.midpoint().pressure_pa());
        assert_eq!(
            a.endpoint().state().temperature_k(),
            b.endpoint().state().temperature_k()
        );
        assert_eq!(
            a.endpoint().state().pressure_pa(),
            b.endpoint().state().pressure_pa()
        );
        assert_eq!(
            a.endpoint().state().sound_speed_m_s(),
            b.endpoint().state().sound_speed_m_s()
        );
    }
}

#[test]
fn real_gas_midpoint_adiabatic_refinement_approaches_original_isentrope() {
    let initial = parcel();
    let exact = initial
        .at_volume_isentropic(
            0.008,
            crate::thermochemical_gas::energy::entropy::GasEntropySolve {
                minimum_temperature_k: 210.0,
                maximum_temperature_k: 450.0,
                entropy_tolerance_j_kg_k: 1.0e-9,
                maximum_iterations: 128,
            },
            || false,
        )
        .unwrap()
        .0;
    let mut errors = Vec::new();
    for count in [4, 8, 16] {
        let mut gas = initial.clone();
        let mut work = 0.0;
        for i in 1..=count {
            let volume = 0.01 - 0.002 * f64::from(i) / f64::from(count);
            let candidate = gas
                .midpoint_heat_volume(volume, 0.0, solve_controls(), || false)
                .unwrap();
            work += candidate.step().work_by_gas_j;
            gas = candidate.into_parts().0;
        }
        let mass = initial.model().mechanics().molar_mass().value();
        let increment = initial
            .model()
            .source()
            .thermodynamics()
            .specific_internal_energy_change(
                Temperature::new(300.0),
                Temperature::new(gas.state().temperature_k()),
            )
            .unwrap()
            .value()
            + (gas.state().root().residual_internal_energy_j_mol()
                - initial.state().root().residual_internal_energy_j_mol())
                / mass;
        near(
            initial.mass_kg() * increment + work,
            0.0,
            2.0 * f64::from(count) * initial.mass_kg() * solve_controls().energy_tolerance_j_kg,
        );
        errors.push((gas.state().temperature_k() - exact.state().temperature_k()).abs());
        assert_eq!(gas.mass_kg(), initial.mass_kg());
    }
    assert!(errors[0] / errors[1] > 3.7, "midpoint errors={errors:?}");
    assert!(errors[1] / errors[2] > 3.7, "midpoint errors={errors:?}");
}

#[test]
fn real_gas_midpoint_refuses_both_stage_domains_and_cancellation_before_publication() {
    let gas = parcel();
    let saved = gas.clone();
    for (volume, heat) in [
        (0.0, 0.0),
        (f64::NAN, 0.0),
        (0.01, f64::INFINITY),
        (0.01, 1.0e12),
        (0.01, -1.0e12),
        (1.0e-12, 0.0),
    ] {
        assert!(
            gas.midpoint_heat_volume(volume, heat, solve_controls(), || false)
                .is_err()
        );
        assert_eq!(gas, saved);
    }
    assert!(
        gas.midpoint_heat_volume(
            0.008,
            2500.0,
            GasEnergySolve {
                maximum_iterations: 1,
                ..solve_controls()
            },
            || false
        )
        .is_err()
    );
    assert!(
        gas.midpoint_heat_volume(
            gas.volume_m3(),
            0.0,
            GasEnergySolve {
                maximum_iterations: 0,
                ..solve_controls()
            },
            || false
        )
        .is_err()
    );
    // T_mid is about 371 K and T_final=440 K. A 340 K transport ceiling
    // refuses the midpoint; a 400 K ceiling admits it but refuses the endpoint.
    for ceiling in [340.0, 400.0] {
        let source = model_in_domains(
            vec![0.75, 0.25],
            0.001,
            [200.0, 2000.0],
            [200.0, ceiling],
            [1000.0, 1.0e7],
        );
        let bounded = IsochoricRealGasParcel::new(
            bound(source),
            point(300.0, 5.0e6),
            0.01,
            PengRobinsonPhaseV1::Vapor,
        )
        .unwrap();
        let heat = bounded.mass_kg()
            * bounded
                .model()
                .specific_internal_energy_change_at_density(bounded.density_kg_m3(), 300.0, 440.0)
                .unwrap();
        assert!(matches!(
            bounded.midpoint_heat_volume(bounded.volume_m3(), heat, solve_controls(), || false),
            Err(GasEnergyError::Material(_))
        ));
    }
    let mut total_polls = 0;
    let expected = gas
        .midpoint_heat_volume(0.008, 2500.0, solve_controls(), || {
            total_polls += 1;
            false
        })
        .unwrap();
    for stop in [1, total_polls / 2, total_polls] {
        let mut polls = 0;
        assert!(matches!(
            gas.midpoint_heat_volume(0.008, 2500.0, solve_controls(), || {
                polls += 1;
                polls >= stop
            }),
            Err(GasEnergyError::Cancelled)
        ));
        assert_eq!(gas, saved);
    }
    assert_eq!(
        gas.midpoint_heat_volume(0.008, 2500.0, solve_controls(), || false)
            .unwrap(),
        expected
    );
}
