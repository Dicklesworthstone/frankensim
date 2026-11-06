use super::*;
use fs_thermochem::fusion::{KELVIN_PER_KEV, NuclearSpeciesV1, SPEED_OF_LIGHT_M_S};
pub(super) mod history;
mod particles;
#[path = "../../examples/dt_fusion/source.rs"]
mod source;
mod transport;

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}
fn controls() -> DtFusionSolve {
    let mut c = DtFusionSolve::new(
        Temperature::new(0.2 * KELVIN_PER_KEV),
        Temperature::new(100.0 * KELVIN_PER_KEV),
        Temperature::new(1e4),
        Temperature::new(20.0 * KELVIN_PER_KEV),
        DtFusionPlasmaLimitsV1 {
            maximum_coulomb_coupling: 0.01,
            maximum_thermal_phase_space_density: 0.01,
            maximum_thermal_energy_rest_energy_ratio: 0.05,
        },
    );
    c.relative_tolerance = 1e-9;
    c.absolute_amount_tolerance_mol = 1e-20;
    c.energy_tolerance_j = 1e-8;
    c.initial_step_s = 1e-3;
    c.maximum_step_s = 0.1;
    c
}
fn heating(fraction: f64) -> DtFusionHeating {
    DtFusionHeating {
        alpha_electron_fraction: fraction,
        electron_power: Power::new(0.0),
        heavy_power: Power::new(0.0),
    }
}
fn run(
    density: [f64; 3],
    volume: f64,
    temperatures: [f64; 2],
    fraction: f64,
    duration: f64,
) -> DtFusionResult {
    let (model, state) = source::plasma(density, volume, temperatures);
    advance_dt_fusion(
        &state,
        &model,
        heating(fraction),
        Time::new(duration),
        controls(),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_frozen_depletion_matches_balanced_and_unequal_analytic_solutions() {
    let rate = 1.1e-22;
    let a = 5e19 / AVOGADRO;
    let (remaining, extent) = solve::frozen_depletion(a, 0.0, rate, 1.0, 1000.0).unwrap();
    let expected = a / (1.0 + rate * 5e19 * 1000.0);
    close(remaining, expected, 8e-15, 0.0);
    close(remaining + extent, a, 8e-15, 0.0);
    let (a, delta, volume, dt) = (2e-4, 7e-4, 0.01, 2.0);
    let ratio = a / (a + delta) * (-delta * AVOGADRO * rate * dt / volume).exp();
    let expected = delta * ratio / (1.0 - ratio);
    let (whole, whole_extent) = solve::frozen_depletion(a, delta, rate, volume, dt).unwrap();
    let (half, x1) = solve::frozen_depletion(a, delta, rate, volume, dt / 2.0).unwrap();
    let (twice, x2) = solve::frozen_depletion(half, delta, rate, volume, dt / 2.0).unwrap();
    close(whole, expected, 2e-14, 0.0);
    close(twice, whole, 2e-14, 0.0);
    close(x1 + x2, whole_extent, 2e-14, 0.0);
    // At an exponent of 40 the limiting fuel remains a representable trace,
    // even though its consumed amount rounds almost to the original stock.
    let trace_dt = 40.0 * volume / (delta * AVOGADRO * rate);
    let (trace, consumed) = solve::frozen_depletion(a, delta, rate, volume, trace_dt).unwrap();
    assert!(trace > 0.0 && trace < a * 1e-16);
    assert!(consumed <= a);
    // Regression: independent fractions previously rounded the denominator
    // below one and falsely consumed over 0.1 mol from a valid 0.1/0.4 mol pair.
    let large_dt = 100.0 / (0.3 * AVOGADRO * rate);
    let (remaining, extent) = solve::frozen_depletion(0.1, 0.3, rate, 1.0, large_dt).unwrap();
    assert!(remaining > 0.0 && extent <= 0.1);
    assert_eq!(
        solve::frozen_depletion(0.0, delta, rate, volume, dt).unwrap(),
        (0.0, 0.0)
    );
}

#[test]
fn g1_balanced_burn_matches_independent_ode_oracle() {
    // Independent fixed-step RK4, platform cbrt/exp, with compensated physical
    // increments. Step-halving 16/32/64 shows fourth-order convergence; 1024,
    // 2048,4096,8192 agree at reported precision. A five-point numerical slope
    // of ln(rate) also reproduces these endpoints. This tests the declared
    // nonrelativistic model, not an observed plasma or source uncertainty.
    let r = run([5e19, 5e19, 0.0], 1.0, [10.0, 5.0], 0.8, 0.1);
    close(
        r.balance().reaction_extent_mol,
        4.790395236840328e-8,
        3e-7,
        1e-18,
    );
    close(
        r.state().heavy_temperature().value(),
        117545084.91871865,
        3e-8,
        0.01,
    );
    close(
        r.state().electron_temperature().value(),
        64381117.94537633,
        3e-8,
        0.01,
    );
    close(
        r.balance().neutron_energy_escaped_j,
        65093.633606311225,
        3e-7,
        1e-7,
    );
    close(
        r.balance().alpha_energy_deposited_j,
        16460.42701906106,
        3e-7,
        1e-7,
    );
    for i in r.plasma().indices()[..2].iter().copied() {
        close(
            r.state().amounts()[i].value(),
            8.297904940632394e-5,
            3e-10,
            1e-18,
        );
    }
}

#[test]
fn g1_unequal_burn_matches_independent_ode_oracle() {
    let r = run([3e20, 7e20, 0.0], 0.01, [10.0, 3.0], 0.65, 0.02);
    close(
        r.balance().reaction_extent_mol,
        8.313099089216964e-9,
        3e-7,
        1e-19,
    );
    close(
        r.state().heavy_temperature().value(),
        120716893.98918392,
        3e-8,
        0.01,
    );
    close(
        r.state().electron_temperature().value(),
        43779564.615858406,
        3e-8,
        0.01,
    );
    close(
        r.balance().neutron_energy_escaped_j,
        11296.376682967757,
        3e-7,
        1e-7,
    );
    close(
        r.balance().alpha_energy_deposited_j,
        2856.6722509200677,
        3e-7,
        1e-7,
    );
    let [d, t, _] = r.plasma().indices();
    close(
        r.state().amounts()[d].value(),
        4.973304102432323e-6,
        3e-10,
        1e-19,
    );
    close(
        r.state().amounts()[t].value(),
        1.161546037112771e-5,
        3e-10,
        1e-19,
    );
}

#[test]
fn g0_burn_closes_source_energy_charge_baryons_and_mass_defect() {
    let r = run([5e19, 5e19, 2e18], 1.0, [10.0, 5.0], 0.8, 0.1);
    let initial = r.initial_state();
    let state = r.state();
    let b = r.balance();
    let emission = r.neutron_emission().unwrap();
    assert_eq!(emission.duration(), r.elapsed_time());
    assert_eq!(emission.total_amount(), r.neutron_amount());
    assert_eq!(emission.total_kinetic_energy(), r.escaped_neutron_energy());
    assert_eq!(
        emission.channels()[0].reaction(),
        FusionReactionV1::DeuteriumTritium,
    );
    assert_eq!(emission.channels()[1].amount().value(), 0.0);
    assert_eq!(emission.channels()[1].mean_particle_energy(), None);
    close(
        emission.channels()[0]
            .mean_particle_energy()
            .unwrap()
            .value()
            * emission.total_particle_count(),
        b.neutron_energy_escaped_j,
        4e-15,
        0.0,
    );
    let [d, t, he] = r.plasma().indices();
    let e = r.plasma().electron_index();
    assert_eq!(initial.amounts()[e], state.amounts()[e]);
    assert_eq!(initial.thermodynamics(), state.thermodynamics());
    for i in [d, t] {
        close(
            initial.amounts()[i].value() - state.amounts()[i].value(),
            b.reaction_extent_mol,
            1e-10,
            1e-18,
        );
    }
    close(
        state.amounts()[he].value() - initial.amounts()[he].value(),
        b.reaction_extent_mol,
        1e-10,
        1e-18,
    );
    assert!(state.pressure().value() > initial.pressure().value());
    assert!(b.relative_baryon_residual.abs() < 1e-12);
    assert!(b.relative_charge_residual.abs() < 1e-12);
    assert!(b.total_energy_residual_j.abs() < 1e-6);
    assert!(b.heavy_energy_residual_j.abs() < 1e-6);
    assert!(b.electron_energy_residual_j.abs() < 1e-6);
    close(
        b.alpha_energy_deposited_j + b.neutron_energy_escaped_j - b.reacting_ion_energy_removed_j,
        b.nuclear_energy_released_j,
        2e-13,
        1e-8,
    );
    assert!(
        b.neutron_energy_escaped_j
            > r.neutron_amount().value()
                * AVOGADRO
                * FusionReactionV1::DeuteriumTritium.neutron_q_j()
    );
    let rest_loss = initial.mass().value()
        - state.mass().value()
        - b.reaction_extent_mol * AVOGADRO * NuclearSpeciesV1::Neutron.mass_kg();
    close(
        rest_loss * SPEED_OF_LIGHT_M_S * SPEED_OF_LIGHT_M_S,
        b.nuclear_energy_released_j,
        2e-9,
        1e-4,
    );
    assert!(b.maximum_normalized_local_error <= 1.0);
}

#[test]
fn g1_no_fuel_preserves_composition_and_uses_original_thermal_capacities() {
    let (model, state) = source::plasma([0.0, 5e19, 2e19], 1.0, [10.0, 5.0]);
    let h = DtFusionHeating {
        alpha_electron_fraction: 0.5,
        electron_power: Power::new(1000.0),
        heavy_power: Power::new(-100.0),
    };
    let r = advance_dt_fusion(&state, &model, h, Time::new(0.02), controls(), || false).unwrap();
    assert_eq!(r.state().amounts(), state.amounts());
    assert_eq!(r.neutron_amount().value(), 0.0);
    assert_eq!(r.escaped_neutron_energy().value(), 0.0);
    assert!(
        r.neutron_emission()
            .unwrap()
            .channels()
            .iter()
            .all(|channel| channel.mean_particle_energy().is_none())
    );
    close(
        r.state().electron_temperature().value() - state.electron_temperature().value(),
        20.0 / state.electron_heat_capacity_j_k(),
        1e-10,
        1e-7,
    );
    close(
        r.state().heavy_temperature().value() - state.heavy_temperature().value(),
        -2.0 / state.heavy_heat_capacity_j_k(),
        1e-10,
        1e-7,
    );
    let zero = advance_dt_fusion(&state, &model, h, Time::new(0.0), controls(), || false).unwrap();
    assert_eq!(zero.state(), &state);
    assert_eq!(zero.balance().accepted_steps, 0);
    assert_eq!(zero.neutron_emission().unwrap().duration(), Time::new(0.0));
    // A zero-time request at either inclusive rate endpoint must retain the
    // exact original state despite independently rounded NASA-9 calorics.
    for ion_kev in [0.2, 100.0] {
        let (model, endpoint) = source::plasma([5e19, 5e19, 0.0], 1.0, [ion_kev, 5.0]);
        let zero =
            advance_dt_fusion(&endpoint, &model, h, Time::new(0.0), controls(), || false).unwrap();
        assert_eq!(zero.state(), &endpoint);
        assert_eq!(zero.balance().accepted_steps, 0);
    }
}

#[test]
fn g3_volume_scaling_preserves_temperature_and_scales_extensive_accounts() {
    let a = run([5e19, 5e19, 0.0], 1.0, [10.0, 5.0], 0.8, 0.1);
    let b = run([5e19, 5e19, 0.0], 0.01, [10.0, 5.0], 0.8, 0.1);
    close(
        a.state().heavy_temperature().value(),
        b.state().heavy_temperature().value(),
        1e-8,
        0.01,
    );
    close(
        a.state().electron_temperature().value(),
        b.state().electron_temperature().value(),
        1e-8,
        0.01,
    );
    close(
        a.balance().reaction_extent_mol,
        100.0 * b.balance().reaction_extent_mol,
        1e-7,
        1e-18,
    );
    close(
        a.balance().neutron_energy_escaped_j,
        100.0 * b.balance().neutron_energy_escaped_j,
        1e-7,
        1e-7,
    );
}

#[test]
fn g3_rate_uses_ion_temperature_and_thermalization_split_is_explicit() {
    let cold_e = run([5e19, 5e19, 0.0], 1.0, [10.0, 2.0], 0.0, 0.001);
    let warm_e = run([5e19, 5e19, 0.0], 1.0, [10.0, 5.0], 0.0, 0.001);
    close(
        cold_e.balance().reaction_extent_mol,
        warm_e.balance().reaction_extent_mol,
        1e-9,
        1e-20,
    );
    close(
        cold_e.state().electron_temperature().value(),
        cold_e.initial_state().electron_temperature().value(),
        1e-14,
        1e-7,
    );
    assert!(
        cold_e.state().heavy_temperature().value()
            > cold_e.initial_state().heavy_temperature().value()
    );
    let electron_heating = run([5e19, 5e19, 0.0], 1.0, [10.0, 5.0], 1.0, 0.001);
    assert!(
        electron_heating.state().electron_temperature().value()
            > electron_heating
                .initial_state()
                .electron_temperature()
                .value()
    );
    assert!(
        electron_heating.state().heavy_temperature().value()
            < electron_heating.initial_state().heavy_temperature().value()
    );
}

#[test]
fn g4_cancellation_budget_and_source_domain_refusals_are_atomic() {
    let (model, initial) = source::plasma([5e19, 5e19, 0.0], 1.0, [10.0, 5.0]);
    let snapshot = initial.clone();
    let mut polls = 0;
    let result = advance_dt_fusion(
        &initial,
        &model,
        heating(0.8),
        Time::new(0.1),
        controls(),
        || {
            polls += 1;
            polls > 50
        },
    );
    assert!(matches!(result, Err(DtFusionError::Cancelled)));
    let mut budget = controls();
    budget.maximum_steps = 1;
    budget.initial_step_s = 1e-6;
    budget.maximum_step_s = 1e-6;
    assert!(matches!(
        advance_dt_fusion(
            &initial,
            &model,
            heating(0.8),
            Time::new(0.1),
            budget,
            || false
        ),
        Err(DtFusionError::StepBudget)
    ));
    let mut outside = controls();
    outside.minimum_heavy_temperature_k = 0.1 * KELVIN_PER_KEV;
    assert!(matches!(
        advance_dt_fusion(
            &initial,
            &model,
            heating(0.8),
            Time::new(0.1),
            outside,
            || false
        ),
        Err(DtFusionError::Reaction(
            FusionErrorV1::OutsideTemperatureDomain { .. }
        ))
    ));
    let mut invalid = heating(0.8);
    invalid.alpha_electron_fraction = 1.1;
    assert!(matches!(
        advance_dt_fusion(
            &initial,
            &model,
            invalid,
            Time::new(0.1),
            controls(),
            || false
        ),
        Err(DtFusionError::Invalid(_))
    ));
    assert_eq!(initial, snapshot);
}

#[test]
fn g5_repeat_and_restart_retain_physical_state_and_accounts() {
    let (model, initial) = source::plasma([5e19, 5e19, 0.0], 1.0, [10.0, 5.0]);
    let advance = |state: &FrozenIonizedGasStateV1, duration| {
        advance_dt_fusion(
            state,
            &model,
            heating(0.8),
            Time::new(duration),
            controls(),
            || false,
        )
        .unwrap()
    };
    let whole = advance(&initial, 0.1);
    assert_eq!(whole, advance(&initial, 0.1));
    let first = advance(&initial, 0.05);
    let second = advance(first.state(), 0.05);
    close(
        second.state().heavy_temperature().value(),
        whole.state().heavy_temperature().value(),
        2e-8,
        0.01,
    );
    close(
        second.state().electron_temperature().value(),
        whole.state().electron_temperature().value(),
        2e-8,
        0.01,
    );
    close(
        first.balance().reaction_extent_mol + second.balance().reaction_extent_mol,
        whole.balance().reaction_extent_mol,
        2e-7,
        1e-18,
    );
    close(
        first.balance().neutron_energy_escaped_j + second.balance().neutron_energy_escaped_j,
        whole.balance().neutron_energy_escaped_j,
        2e-7,
        1e-7,
    );
}

#[test]
fn g1_burn_state_reuses_existing_coulomb_thermal_owner_in_its_domain() {
    use crate::chemistry::energy::{ElectronHeavyEnergySolve, advance_frozen_energy_with_coulomb};
    use fs_thermochem::coulomb::{CoulombLimitsV1, CoulombModelV1, CoulombScreeningV1};
    use fs_thermochem::electron_collision::MaxwellianRateControlsV1;
    // A deliberately cold electron pool keeps the original CLASSICAL Coulomb
    // cutoff admissible. This tests interface reuse, not a typical operating
    // fusion plasma; hot-electron transport needs a quantum-cutoff owner.
    let burn = run([5e19, 5e19, 0.0], 1.0, [10.0, 0.005], 0.0, 1e-5);
    let state = burn.state();
    let collision = CoulombModelV1::new(
        state,
        CoulombScreeningV1::StaticAllChargedDebye,
        CoulombLimitsV1 {
            minimum_coulomb_logarithm: 4.0,
            minimum_particles_per_debye_sphere: 10.0,
            maximum_coulomb_coupling: 0.1,
            maximum_thermal_phase_space_density: 0.1,
            maximum_thermal_energy_rest_energy_ratio: 0.01,
            maximum_quantum_to_classical_ratio: 0.8,
            maximum_relative_drift_thermal_ratio: 0.1,
        },
        || false,
    )
    .unwrap();
    let mut c = ElectronHeavyEnergySolve::new(
        Temperature::new(1e4),
        Temperature::new(2e9),
        1.0,
        0.1,
        MaxwellianRateControlsV1 {
            relative_tolerance: 1e-8,
            maximum_subdivisions: 1024,
        },
    );
    c.energy_tolerance_j = 1e-6;
    c.initial_step_s = 1e-11;
    c.maximum_step_s = 1e-10;
    let result = advance_frozen_energy_with_coulomb(
        state,
        &[],
        &collision,
        Power::new(0.0),
        Power::new(0.0),
        Time::new(1e-10),
        c,
        || false,
    )
    .unwrap();
    assert_eq!(result.state().amounts(), state.amounts());
    assert_eq!(result.state().thermodynamics(), state.thermodynamics());
    assert!(result.state().electron_temperature().value() > state.electron_temperature().value());
    assert!(result.state().heavy_temperature().value() < state.heavy_temperature().value());
    close(
        result.state().internal_energy().value(),
        state.internal_energy().value(),
        1e-10,
        1e-6,
    );
}
