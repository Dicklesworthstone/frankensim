use super::super::WaterCellControls;
use super::*;
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::phase::water::isochoric::IsochoricWaterControls;
use fs_material::phase::water::{WaterPhase, acoustics::HomogeneousEquilibriumWater};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::{If97Region, If97SaturationPair, If97State};

fn cell_controls(lo: f64, hi: f64) -> WaterCellControls {
    WaterCellControls {
        equilibrium: IsochoricWaterControls {
            minimum_temperature_k: lo,
            maximum_temperature_k: hi,
            minimum_pressure_pa: 100.0,
            maximum_pressure_pa: 100e6,
            volume_tolerance_m3_kg: 1e-13,
            internal_energy_tolerance_j_kg: 1e-7,
            maximum_pressure_iterations: 96,
            maximum_temperature_iterations: 96,
        },
        maximum_relative_mass_exchange: 0.1,
        maximum_inflow_parcels: 4,
        mass_tolerance_kg: 1e-13,
        energy_tolerance_j: 1e-7,
        volume_tolerance_m3: 1e-13,
    }
}
fn controls() -> WaterVolumeControls {
    WaterVolumeControls {
        maximum_step_s: 1.0,
        maximum_relative_volume_change: 0.1,
        maximum_temperature_increment_k: 20.0,
        energy_tolerance_j: 1e-7,
        maximum_iterations: 96,
    }
}
fn mixed(temperature: f64, quality: f64) -> EquilibriumWaterCell {
    let pair = If97SaturationPair::at_temperature(Temperature::new(temperature)).unwrap();
    let v = (1.0 - quality) * pair.liquid().specific_volume().value()
        + quality * pair.vapor().specific_volume().value();
    EquilibriumWaterCell::at_temperature(
        0.05 * v,
        0.05,
        Temperature::new(temperature),
        0.0,
        cell_controls(temperature - 20.0, temperature + 20.0),
        || false,
    )
    .unwrap()
}
fn steam() -> EquilibriumWaterCell {
    let source = If97State::at_pressure_temperature(
        Pressure::new(100_000.0),
        Temperature::new(500.0),
        If97Region::Region2,
    )
    .unwrap();
    EquilibriumWaterCell::new(
        0.01 * source.specific_volume().value(),
        0.01,
        0.01 * source.specific_internal_energy().value(),
        0.0,
        cell_controls(480.0, 520.0),
        || false,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn rigid_limit_is_the_original_heat_update_and_volume_proposals_are_atomic() {
    let original = mixed(400.0, 0.02);
    let mut rigid = original;
    let mut driven = original;
    let v = original.state().volume_m3();
    let expected = rigid.advance(&[], 0.0, 10.0, 0.5, || false).unwrap();
    let actual = driven
        .advance_volume(v, 10.0, 0.5, controls(), || false)
        .unwrap();
    assert_eq!(expected, actual.water);
    assert_eq!(driven, rigid);
    assert_eq!(actual.water.boundary_work_j, 0.0);

    let mut candidate_owner = original;
    let mut calls = 0;
    {
        let candidate = candidate_owner
            .propose_volume(v * 1.02, 0.0, 0.5, controls(), || {
                calls += 1;
                false
            })
            .unwrap();
        assert!(candidate.state().volume_m3() > v);
        assert_eq!(candidate.state().mass_kg(), original.state().mass_kg());
    }
    assert_eq!(candidate_owner, original);
    let mut polls = 0;
    let refused = candidate_owner.propose_volume(v * 1.02, 0.0, 0.5, controls(), || {
        polls += 1;
        polls == calls
    });
    assert!(matches!(
        refused,
        Err(WaterVolumeError::Water(WaterCellError::Cancelled))
    ));
    drop(refused);
    assert_eq!(candidate_owner, original);
    assert!(
        candidate_owner
            .advance_volume(v * 1.2, 0.0, 0.5, controls(), || false)
            .is_err()
    );
    assert_eq!(candidate_owner, original);
    let mut short = controls();
    short.maximum_temperature_increment_k = 1e-6;
    assert!(
        candidate_owner
            .advance_volume(v * 1.02, 100.0, 0.5, short, || false)
            .is_err()
    );
    assert_eq!(candidate_owner, original);

    let mut a = original;
    let mut b = original;
    let sa = a
        .advance_volume(v * 1.02, 0.0, 0.5, controls(), || false)
        .unwrap();
    let sb = b
        .advance_volume(v * 1.02, 0.0, 0.5, controls(), || false)
        .unwrap();
    assert_eq!(sa, sb);
    assert_eq!(a, b);
}

#[test]
fn isothermal_phase_expansion_and_condensation_use_original_latent_enthalpy_once() {
    let mut cell = mixed(400.0, 0.01);
    let initial = cell.state();
    let pair = initial.source_state().saturation_pair().unwrap();
    let mass = initial.mass_kg();
    let dx = 0.0005;
    let dv = mass
        * dx
        * (pair.vapor().specific_volume().value() - pair.liquid().specific_volume().value());
    // Independent coexistence identity: at fixed T,p, Q=Delta m_v*h_fg.
    let heat = mass * dx * pair.latent_enthalpy().value();
    let expected_work = pair.pressure().value() * dv;
    let step = cell
        .advance_volume(initial.volume_m3() + dv, heat, 0.5, controls(), || false)
        .unwrap();
    near(step.water.state.source_state().temperature_k(), 400.0, 2e-8);
    near(
        step.water.state.phase_masses_kg()[1],
        initial.phase_masses_kg()[1] + mass * dx,
        2e-13,
    );
    near(step.water.boundary_work_j, expected_work, 1e-8);
    near(
        step.water.internal_energy_change_j,
        heat - expected_work,
        1e-7,
    );
    assert_eq!(step.water.mass_change_kg, 0.0);
    near(step.water.energy_residual_j, 0.0, 1e-7);
    near(step.source_energy_root_residual_j, 0.0, 1e-7);
    // Original entropy is retained: the finite Gibbs seam is not patched.
    let native_entropy = mass * dx * pair.latent_entropy().value();
    near(step.water.source_entropy_change_j_k, native_entropy, 1e-9);

    let back = cell
        .advance_volume(initial.volume_m3(), -heat, 1.0, controls(), || false)
        .unwrap();
    near(back.water.boundary_work_j, -expected_work, 1e-8);
    near(cell.state().source_state().temperature_k(), 400.0, 2e-8);
    near(
        cell.state().internal_energy_j(),
        initial.internal_energy_j(),
        1e-7,
    );
    near(
        cell.state().phase_masses_kg()[1],
        initial.phase_masses_kg()[1],
        2e-13,
    );
}

#[test]
fn adiabatic_single_phase_expansion_converges_to_independent_entropy_invariant() {
    let original = steam();
    let initial = original.state();
    let entropy = initial.source_state().specific_entropy_j_kg_k();
    let mut errors = Vec::new();
    for count in [2, 4, 8] {
        let mut cell = original;
        let mut work = 0.0;
        for i in 1..=count {
            let v = initial.volume_m3() * (1.0 + 0.08 * i as f64 / count as f64);
            let step = cell
                .advance_volume(v, 0.0, i as f64 / count as f64, controls(), || false)
                .unwrap();
            assert!(step.water.boundary_work_j > 0.0);
            work += step.water.boundary_work_j;
        }
        assert!(
            cell.state().source_state().temperature_k() < initial.source_state().temperature_k()
        );
        assert!(cell.state().source_state().pressure_pa() < initial.source_state().pressure_pa());
        near(
            cell.state().internal_energy_j() - initial.internal_energy_j() + work,
            0.0,
            3e-7,
        );
        errors.push((cell.state().source_state().specific_entropy_j_kg_k() - entropy).abs());
    }
    // This checks an independent thermodynamic invariant of the original
    // single-phase Gibbs source, not merely the discrete work equation.
    assert!(errors[1] < 0.35 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.35 * errors[1], "{errors:?}");
    let mut refused = original;
    let mut budget = controls();
    budget.maximum_iterations = 1;
    assert!(
        refused
            .advance_volume(initial.volume_m3() * 1.08, 0.0, 1.0, budget, || false)
            .is_err()
    );
    assert_eq!(refused, original);
}

#[test]
fn sealed_liquid_expansion_cavitates_and_recompression_recovers_liquid() {
    let liquid = If97State::at_pressure_temperature(
        Pressure::new(100_000.0),
        Temperature::new(300.0),
        If97Region::Region1,
    )
    .unwrap();
    let mut cell = EquilibriumWaterCell::new(
        liquid.specific_volume().value(),
        1.0,
        liquid.specific_internal_energy().value(),
        0.0,
        cell_controls(299.5, 300.5),
        || false,
    )
    .unwrap();
    let initial = cell.state();
    let expanded = cell
        .advance_volume(initial.volume_m3() * 1.002, 0.0, 0.5, controls(), || false)
        .unwrap();
    assert_eq!(
        expanded.water.state.source_state().phase(),
        WaterPhase::LiquidVapor
    );
    assert!(expanded.water.state.phase_masses_kg()[1] > 0.0);
    assert_eq!(expanded.water.state.mass_kg(), initial.mass_kg());
    assert!(expanded.water.state.source_state().pressure_pa() < 10_000.0);
    assert!(expanded.water.internal_energy_change_j < 0.0);
    near(expanded.water.energy_residual_j, 0.0, 1e-7);
    let compressed = cell
        .advance_volume(initial.volume_m3() * 0.9998, 0.0, 1.0, controls(), || false)
        .unwrap();
    assert_eq!(
        compressed.water.state.source_state().phase(),
        WaterPhase::CompressedLiquid
    );
    assert_eq!(compressed.water.state.phase_masses_kg()[1], 0.0);
    assert_eq!(compressed.water.state.mass_kg(), initial.mass_kg());
    assert!(compressed.water.boundary_work_j < 0.0);
    assert!(compressed.water.internal_energy_change_j > 0.0);
}

#[test]
fn finite_adiabatic_volume_step_matches_native_equilibrium_acoustic_tangent() {
    let mut cell = mixed(400.0, 0.02);
    let initial = cell.state();
    let source = initial.source_state();
    let acoustic = HomogeneousEquilibriumWater::try_new(
        source.saturation_pair().unwrap(),
        source.vapor_mass_fraction(),
    )
    .unwrap()
    .acoustic_state()
    .unwrap();
    let dv = initial.volume_m3() * 1e-5;
    let step = cell
        .advance_volume(initial.volume_m3() + dv, 0.0, 0.5, controls(), || false)
        .unwrap();
    let observed = -(step.water.state.source_state().pressure_pa() - source.pressure_pa())
        * initial.volume_m3()
        / dv;
    near(
        observed,
        acoustic.adiabatic_bulk_modulus_pa(),
        3e-4 * acoustic.adiabatic_bulk_modulus_pa(),
    );
    assert!(step.water.state.phase_masses_kg()[1] > initial.phase_masses_kg()[1]);
}
