//! G1/G3: finite original PR/NASA calorics, entropy, work and implicit tangents.
use super::*;
use crate::thermochemical_gas::energy::entropy::GasEntropySolve;
use crate::thermochemical_gas::energy::parcel::ClosedGasParcel;
use crate::thermochemical_gas::nonideal::energy::parcel::discrete_gradient::RealGasDiscreteGradientCandidate;
use crate::thermochemical_gas::tests::model_in_domains;

fn solve_controls() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 210.0,
        maximum_temperature_k: 450.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}
fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual / expected - 1.0).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, relative tolerance {tolerance:e}"
    );
}
fn balances(old: &IsochoricRealGasParcel, candidate: &RealGasDiscreteGradientCandidate) {
    let step = candidate.step();
    let tolerance = old.mass_kg() * solve_controls().energy_tolerance_j_kg;
    assert!(step.energy_residual_j.abs() <= tolerance);
    assert!(step.gibbs_residual_j.abs() <= tolerance);
    assert!(
        step.entropy_residual_j_k.abs()
            <= tolerance / candidate.thermodynamic_mean().temperature_k()
    );
    assert_eq!(
        step.work_by_gas_j,
        step.work_pressure_pa * (step.final_volume_m3 - step.initial_volume_m3)
    );
    assert_eq!(candidate.endpoint().mass_kg(), old.mass_kg());
    assert_eq!(candidate.endpoint().model(), old.model());
    assert!(step.pressure_final_volume_derivative_pa_m3 < 0.0);
    assert!(step.pressure_heat_derivative_pa_j.is_finite());
}

#[test]
fn dense_gibbs_steps_match_independent_coordinate_integrals_and_keep_real_sampling_states() {
    let gas = parcel();
    let saved = gas.clone();
    // Independent closed-form PR attraction/integral and NASA Cv/Cv-over-T
    // evaluation, with a separate bisection root and no production state calls.
    // These synthetic double-precision oracles are not measured fluid evidence.
    for (volume, heat, t1, theta, pressure, p1, work, entropy, dpdv, dpdq) in [
        (
            0.008,
            2500.0,
            337.2438649035121,
            318.28786114882445,
            5959161.917185798,
            7129046.013348099,
            -11918.323834371597,
            7.854525117535189,
            -547189260.9171518,
            25.25343398601428,
        ),
        (
            0.012,
            -1500.0,
            273.4160571467964,
            286.5175863363745,
            4335651.228358884,
            3768970.7565968344,
            8671.302456717769,
            -5.235280735050517,
            -270226093.1508604,
            20.147722998580537,
        ),
        (
            0.01,
            30000.0,
            364.77082151482114,
            331.4210007396757,
            5645146.278260178,
            6323168.417925255,
            0.0,
            90.51930907530019,
            -408581208.73420644,
            20.734433761530276,
        ),
    ] {
        let candidate = gas
            .discrete_gradient_heat_volume(volume, heat, solve_controls(), || false)
            .unwrap();
        balances(&gas, &candidate);
        let step = candidate.step();
        near(candidate.endpoint().state().temperature_k(), t1, 2e-8);
        near(candidate.endpoint().state().pressure_pa(), p1, 2e-3);
        near(candidate.thermodynamic_mean().temperature_k(), theta, 2e-8);
        near(step.work_pressure_pa, pressure, 2e-3);
        near(step.work_by_gas_j, work, 2e-5);
        near(step.stored_energy_change_j, heat - work, 2e-5);
        near(step.entropy_change_j_k, entropy, 2e-8);
        relative(step.pressure_final_volume_derivative_pa_m3, dpdv, 2e-9);
        relative(step.pressure_heat_derivative_pa_j, dpdq, 2e-9);
        // The original material pressure is not overwritten by the conjugate
        // pressure, including at fixed volume with appreciable heating.
        assert!(
            (step.work_pressure_pa - candidate.thermodynamic_mean().pressure_pa()).abs() > 1000.0
        );
        let mean_density = if volume == gas.volume_m3() {
            gas.density_kg_m3()
        } else {
            gas.mass_kg() / f64::midpoint(gas.volume_m3(), volume)
        };
        assert_eq!(
            candidate.thermodynamic_mean().density_kg_m3().to_bits(),
            mean_density.to_bits()
        );
        for state in [candidate.thermodynamic_mean(), candidate.endpoint().state()] {
            assert_eq!(
                state.transport_state().unwrap().density_kg_m3(),
                state.density_kg_m3()
            );
            assert_eq!(
                state.acoustic_state().unwrap().density_kg_m3(),
                state.density_kg_m3()
            );
        }
        assert_eq!(
            candidate.endpoint().point().axes().get("material-lot"),
            Some(&12.0)
        );
        assert_eq!(gas, saved);
    }
}

#[test]
fn adiabatic_steps_agree_with_original_entropy_inversion_and_reverse_the_same_path() {
    let gas = parcel();
    let entropy_controls = GasEntropySolve {
        minimum_temperature_k: 210.0,
        maximum_temperature_k: 450.0,
        entropy_tolerance_j_kg_k: 1e-10,
        maximum_iterations: 128,
    };
    for volume in [0.008, 0.012] {
        let candidate = gas
            .discrete_gradient_heat_volume(volume, 0.0, solve_controls(), || false)
            .unwrap();
        balances(&gas, &candidate);
        let (reference, reference_step) = gas
            .at_volume_isentropic(volume, entropy_controls, || false)
            .unwrap();
        near(
            candidate.endpoint().state().temperature_k(),
            reference.state().temperature_k(),
            2e-8,
        );
        near(
            candidate.step().work_by_gas_j,
            reference_step.work_by_gas_j,
            2e-5,
        );
        let reversed = candidate
            .endpoint()
            .discrete_gradient_heat_volume(gas.volume_m3(), 0.0, solve_controls(), || false)
            .unwrap();
        balances(candidate.endpoint(), &reversed);
        near(
            reversed.endpoint().state().temperature_k(),
            gas.state().temperature_k(),
            2e-8,
        );
        relative(
            reversed.step().work_pressure_pa,
            candidate.step().work_pressure_pa,
            2e-10,
        );
        near(
            reversed.step().work_by_gas_j + candidate.step().work_by_gas_j,
            0.0,
            2.0 * gas.mass_kg() * solve_controls().energy_tolerance_j_kg,
        );
    }
}

#[test]
fn isothermal_real_gas_work_includes_density_energy_and_entropy_changes() {
    let gas = parcel();
    for volume in [0.008, 0.012] {
        let rho1 = gas.mass_kg() / volume;
        let ds = gas
            .model()
            .specific_entropy_change_between_densities(gas.density_kg_m3(), 300.0, rho1, 300.0)
            .unwrap();
        let heat = gas.mass_kg() * 300.0 * ds;
        let candidate = gas
            .discrete_gradient_heat_volume(volume, heat, solve_controls(), || false)
            .unwrap();
        balances(&gas, &candidate);
        assert_eq!(candidate.endpoint().state().temperature_k(), 300.0);
        assert_eq!(candidate.thermodynamic_mean().temperature_k(), 300.0);
        assert!(candidate.step().stored_energy_change_j.abs() > 100.0);
        near(
            candidate.step().entropy_change_j_k,
            gas.mass_kg() * ds,
            1e-12,
        );
        near(
            candidate.step().stored_energy_change_j + candidate.step().work_by_gas_j,
            heat,
            gas.mass_kg() * solve_controls().energy_tolerance_j_kg,
        );
    }
}

#[test]
fn implicit_tangents_match_fresh_original_source_solves_through_both_coincident_limits() {
    let gas = parcel();
    for (volume, heat) in [
        (0.008, 1000.0),
        (0.01, 1000.0),
        (0.012, 1000.0),
        (0.010000001, 0.001),
        (0.009999999, 0.001),
        (0.01, 0.001),
    ] {
        let candidate = gas
            .discrete_gradient_heat_volume(volume, heat, solve_controls(), || false)
            .unwrap();
        balances(&gas, &candidate);
        let pressure = |v, q| {
            gas.discrete_gradient_heat_volume(v, q, solve_controls(), || false)
                .unwrap()
                .step()
                .work_pressure_pa
        };
        let dv = 1e-7;
        let dq = 0.5;
        relative(
            candidate.step().pressure_final_volume_derivative_pa_m3,
            (pressure(volume + dv, heat) - pressure(volume - dv, heat)) / (2.0 * dv),
            2e-6,
        );
        relative(
            candidate.step().pressure_heat_derivative_pa_j,
            (pressure(volume, heat + dq) - pressure(volume, heat - dq)) / (2.0 * dq),
            2e-6,
        );
    }
}

#[test]
fn fixed_volume_matches_original_heating_and_zero_step_retains_the_entire_inventory() {
    let gas = parcel();
    let quiet = gas
        .discrete_gradient_heat_volume(gas.volume_m3(), 0.0, solve_controls(), || false)
        .unwrap();
    assert_eq!(quiet.endpoint(), &gas);
    assert_eq!(quiet.thermodynamic_mean(), gas.state());
    assert_eq!(quiet.step().work_pressure_pa, gas.state().pressure_pa());
    assert_eq!(quiet.step().energy_residual_j, 0.0);
    assert_eq!(quiet.step().entropy_change_j_k, 0.0);
    assert_eq!(quiet.step().iterations, 0);
    for heat in [-5000.0, 30000.0, gas.mass_kg() * 1e-3] {
        let mut reference = gas.clone();
        reference
            .add_heat(heat, solve_controls(), || false)
            .unwrap();
        let candidate = gas
            .discrete_gradient_heat_volume(gas.volume_m3(), heat, solve_controls(), || false)
            .unwrap();
        balances(&gas, &candidate);
        near(
            candidate.endpoint().state().temperature_k(),
            reference.state().temperature_k(),
            2e-8,
        );
        near(
            candidate.endpoint().state().pressure_pa(),
            reference.state().pressure_pa(),
            2e-3,
        );
        assert_eq!(
            candidate.endpoint().density_kg_m3().to_bits(),
            gas.density_kg_m3().to_bits()
        );
        assert_eq!(
            candidate.thermodynamic_mean().density_kg_m3().to_bits(),
            gas.density_kg_m3().to_bits()
        );
        assert_eq!(candidate.step().work_by_gas_j, 0.0);
    }
}

#[test]
fn small_coupled_changes_are_invariant_to_large_formation_energy_offsets() {
    let initialize = |reference| {
        let source = ThermochemicalGasModel::new(
            FrozenIdealGasMixtureModelV1::new(
                vec![
                    species_with_energy_reference("A", 0.028, 3.5, 0.001, reference),
                    species_with_energy_reference("B", 0.004, 2.5, 0.0, reference),
                ],
                Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
            )
            .unwrap(),
            vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)],
            None,
        )
        .unwrap();
        IsochoricRealGasParcel::new(
            bound(source),
            point(300.0, 5e6),
            0.01,
            PengRobinsonPhaseV1::Vapor,
        )
        .unwrap()
    };
    let plain = initialize(0.0);
    let shifted = initialize(1e20);
    for volume in [0.01, 0.0100000001, 0.0099999999] {
        let heat = plain.mass_kg() * 1e-3;
        let first = plain
            .discrete_gradient_heat_volume(volume, heat, solve_controls(), || false)
            .unwrap();
        let second = shifted
            .discrete_gradient_heat_volume(volume, heat, solve_controls(), || false)
            .unwrap();
        assert_eq!(first.step(), second.step());
        assert_eq!(
            first.endpoint().state().temperature_k(),
            second.endpoint().state().temperature_k()
        );
        assert_eq!(
            first.endpoint().state().pressure_pa(),
            second.endpoint().state().pressure_pa()
        );
        assert_eq!(
            first.thermodynamic_mean().temperature_k(),
            second.thermodynamic_mean().temperature_k()
        );
    }
}

#[test]
fn dilute_original_pr_limit_recovers_the_ideal_discrete_gibbs_rule_and_its_tangents() {
    let source = model_in_domains(
        vec![0.75, 0.25],
        0.0,
        [200.0, 2000.0],
        [200.0, 2000.0],
        [0.01, 1e7],
    );
    let ideal = ClosedGasParcel::new(source.clone(), point(300.0, 0.1), 0.01).unwrap();
    let real = IsochoricRealGasParcel::new(
        bound(source),
        point(300.0, 0.1),
        0.01,
        PengRobinsonPhaseV1::Vapor,
    )
    .unwrap();
    for volume in [0.008, 0.012] {
        for specific_heat in [-2000.0, 0.0, 5000.0] {
            let a = real
                .discrete_gradient_heat_volume(
                    volume,
                    real.mass_kg() * specific_heat,
                    solve_controls(),
                    || false,
                )
                .unwrap();
            let b = ideal
                .discrete_gradient_heat_volume(
                    volume,
                    ideal.mass_kg() * specific_heat,
                    solve_controls(),
                    || false,
                )
                .unwrap();
            relative(
                a.endpoint().state().temperature_k(),
                b.endpoint().state().state().temperature,
                2e-7,
            );
            relative(
                a.thermodynamic_mean().temperature_k(),
                b.thermodynamic_mean().state().temperature,
                2e-7,
            );
            relative(a.step().work_pressure_pa, b.step().work_pressure_pa, 2e-7);
            relative(
                a.step().pressure_final_volume_derivative_pa_m3,
                b.step().pressure_final_volume_derivative_pa_m3,
                2e-7,
            );
            relative(
                a.step().pressure_heat_derivative_pa_j,
                b.step().pressure_heat_derivative_pa_j,
                2e-7,
            );
            assert_eq!(
                b.step().work_pressure_pa,
                b.thermodynamic_mean().state().pressure
            );
        }
    }
}

#[test]
fn four_stable_corners_cannot_silently_cross_a_peng_robinson_spinodal_rectangle() {
    // Pure A. fs-thermochem refuses a zero declared fraction (mixture.rs:357-364, d571a4691d;
    // CONTRACT.md:2214-2215: absent species are omitted, not represented by zero), so the
    // former `vec![1.0, 0.0]` source refused before the test began. `with_composition` omits
    // the explicit zero B fraction and keeps A's own NASA-9 and transport cards. The
    // Peng-Robinson binding is `bound()`'s A record alone, with no pair interaction. A zero B
    // fraction adds nothing to the PR attraction or covolume, so this is the same pure-A
    // fluid, and the test keeps its corner claim: each of the four (T, rho) corners, T in
    // {150, 151} K and rho in {liquid, vapor} at 150 K and 1 MPa, is admissible as a
    // mechanically stable state on its own, yet the rectangle they span crosses the PR
    // spinodal, so the discrete gradient must refuse it instead of accepting its corners.
    let pure_a_fractions =
        Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap();
    let source = model_in_domains(
        vec![0.5, 0.5],
        0.001,
        [140.0, 220.0],
        [140.0, 220.0],
        [1e4, 1e8],
    )
    .with_composition(pure_a_fractions)
    .unwrap();
    let pure_a = ThermochemicalPengRobinsonModel::new(
        source,
        vec![PengRobinsonSpeciesBinding {
            species: SpeciesId::new("A").unwrap(),
            critical_temperature: Temperature::new(190.564),
            critical_pressure: Pressure::new(4.5992e6),
            acentric_factor: 0.01142,
        }],
        vec![],
    )
    .unwrap();
    let gas =
        IsochoricRealGasParcel::new(pure_a, point(150.0, 1e6), 0.01, PengRobinsonPhaseV1::Liquid)
            .unwrap();
    let vapor = gas
        .state()
        .mechanics()
        .select(PengRobinsonPhaseV1::Vapor)
        .mass_density_kg_per_m3();
    assert!(gas.density_kg_m3() > 5.0 * vapor);
    for t in [150.0, 151.0] {
        for rho in [gas.density_kg_m3(), vapor] {
            assert!(
                gas.model()
                    .evaluate_density_temperature(&point(t, 1e6), rho)
                    .is_ok()
            );
        }
    }
    let limited = GasEnergySolve {
        minimum_temperature_k: 150.0,
        maximum_temperature_k: 151.0,
        ..solve_controls()
    };
    let error = gas
        .discrete_gradient_heat_volume(gas.mass_kg() / vapor, 0.0, limited, || false)
        .unwrap_err();
    assert!(matches!(&error, GasEnergyError::Material(_)));
    assert!(
        error
            .to_string()
            .contains("real-gas Gibbs source rectangle")
    );
}

#[test]
fn source_domain_exhaustion_and_late_cancellation_return_no_partial_candidate() {
    let gas = parcel();
    let saved = gas.clone();
    assert!(matches!(
        gas.discrete_gradient_heat_volume(
            0.008,
            2500.0,
            GasEnergySolve {
                maximum_iterations: 1,
                ..solve_controls()
            },
            || false
        ),
        Err(GasEnergyError::NotConverged { .. })
    ));
    for (volume, heat) in [(0.0, 0.0), (0.01, f64::NAN), (0.01, 1e12)] {
        assert!(
            gas.discrete_gradient_heat_volume(volume, heat, solve_controls(), || false)
                .is_err()
        );
        assert_eq!(gas, saved);
    }
    assert!(
        gas.discrete_gradient_heat_volume(
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
    let dense = IsochoricRealGasParcel::new(
        bound(model(vec![0.75, 0.25], 0.001)),
        point(300.0, 9e6),
        0.01,
        PengRobinsonPhaseV1::Vapor,
    )
    .unwrap();
    let heat = dense.mass_kg()
        * dense
            .model()
            .specific_internal_energy_change_at_density(dense.density_kg_m3(), 300.0, 440.0)
            .unwrap();
    assert!(matches!(
        dense.discrete_gradient_heat_volume(dense.volume_m3(), heat, solve_controls(), || false),
        Err(GasEnergyError::Material(_))
    ));
    let mut calls = 0;
    gas.discrete_gradient_heat_volume(0.008, 2500.0, solve_controls(), || {
        calls += 1;
        false
    })
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        gas.discrete_gradient_heat_volume(0.008, 2500.0, solve_controls(), || {
            seen += 1;
            seen == calls
        }),
        Err(GasEnergyError::Cancelled)
    ));
    assert!(matches!(
        gas.discrete_gradient_heat_volume(0.008, 2500.0, solve_controls(), || true),
        Err(GasEnergyError::Cancelled)
    ));
    assert_eq!(gas, saved);
}
