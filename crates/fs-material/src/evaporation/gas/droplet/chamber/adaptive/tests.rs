//! G0/G1/G3/G4: evolving physical geometry, accepted-only accounts and retry.
use super::{
    AdaptiveGasDropletControls, AdaptiveGasDropletError, advance_droplets_with_gas_adaptive,
    advance_droplets_with_gas_and_heat_adaptive,
};
use crate::evaporation::gas::chamber::{GasLiquidChamberSolve, heat::GasLiquidChamberHeatSolve};
use crate::evaporation::gas::droplet::SphericalDropletState;
use crate::evaporation::gas::droplet::chamber::{
    GasDropletChamberCandidate, GasDropletHeatSurface, GasDropletSurface, advance_droplets_with_gas,
};
use crate::evaporation::gas::tests::{Phase, binding, catalog, controls, gas, liquid, temperature};
use crate::evaporation::gas::{GasEvaporationError, GasFilmHeatSolve, GasLiquidBinding};
use crate::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_thermochem::SpeciesId;

fn row<'a>(id: u64, binding: &'a GasLiquidBinding, mass: f64, t: f64) -> GasDropletSurface<'a> {
    GasDropletSurface {
        id,
        law: &Phase,
        binding,
        droplet: SphericalDropletState {
            mass_kg: mass,
            specific_enthalpy_j_kg: liquid(t).specific_enthalpy_j_kg,
        },
        temperature_from_enthalpy: &temperature,
        binary_diffusivity_m2_s: 2e-5,
        sherwood_number: 2.0,
    }
}
fn rows(binding: &GasLiquidBinding) -> [GasDropletSurface<'_>; 2] {
    [
        row(9, binding, 1e-9, 350.0),
        GasDropletSurface {
            binary_diffusivity_m2_s: 3e-5,
            sherwood_number: 4.0,
            ..row(2, binding, 8e-9, 310.0)
        },
    ]
}
fn solve() -> GasLiquidChamberSolve {
    GasLiquidChamberSolve {
        surface: controls(),
        maximum_surfaces: 8,
    }
}
fn time_controls(span: f64) -> AdaptiveGasDropletControls {
    let mut result = AdaptiveGasDropletControls::new(solve());
    result.initial_dt_s = span;
    result.minimum_dt_s = 1e-6;
    result.maximum_dt_s = span;
    result.maximum_comparisons = 4096;
    result.maximum_reductions = 32;
    result.enthalpy_absolute_tolerance_j_kg = 0.2;
    result.temperature_absolute_tolerance_k = 1e-4;
    result.mass_absolute_tolerance_kg = 5e-15;
    result.radius_absolute_tolerance_m = 1e-9;
    result.relative_tolerance = 0.0;
    result.mass_conservation_tolerance_kg = 1e-12;
    result.energy_conservation_tolerance_j = 1e-7;
    result
}
fn loose_controls(span: f64) -> AdaptiveGasDropletControls {
    AdaptiveGasDropletControls {
        enthalpy_absolute_tolerance_j_kg: 2000.0,
        temperature_absolute_tolerance_k: 1.0,
        mass_absolute_tolerance_kg: 1e-10,
        radius_absolute_tolerance_m: 1e-4,
        ..time_controls(span)
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

// Independent time refinement uses only the existing original chamber owner.
// Every call receives the actual previous physical inventories, so the owner
// recomputes area and Sh*D/diameter instead of retaining the first sphere.
fn fixed_reference<'catalog>(
    inputs: &[GasDropletSurface<'_>],
    initial_gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    count: usize,
) -> GasDropletChamberCandidate<'catalog> {
    let mut rows = inputs.to_vec();
    rows.sort_unstable_by_key(|row| row.id);
    let mut gas = initial_gas.clone();
    let mut last = None;
    for _ in 0..count {
        let next = advance_droplets_with_gas(
            &rows,
            &gas,
            catalog,
            duration_s / count as f64,
            solve(),
            || false,
        )
        .unwrap();
        for (row, step) in rows.iter_mut().zip(next.droplets()) {
            assert_eq!(row.id, step.id);
            row.droplet = step.droplet;
        }
        gas = next.gas().clone();
        last = Some(next);
    }
    last.unwrap()
}
fn temperature_error(
    a: &GasDropletChamberCandidate<'_>,
    b: &GasDropletChamberCandidate<'_>,
) -> f64 {
    a.droplets()
        .iter()
        .zip(b.droplets())
        .map(|(a, b)| {
            assert_eq!(a.id, b.id);
            (a.exchange.surface_temperature_k[1] - b.exchange.surface_temperature_k[1]).abs()
        })
        .fold(
            (a.gas().state().state().temperature - b.gas().state().state().temperature).abs(),
            f64::max,
        )
}

#[test]
fn time_refinement_resolves_shrinking_spheres_and_only_accepted_inventory_and_energy() {
    let source = catalog(-1e7);
    let binding = binding(-1e7);
    let original = gas(&source, 320.0, 200.0, 1e-5);
    let saved = original.clone();
    let inputs = rows(&binding);
    let duration_s = 0.4;
    let loose = advance_droplets_with_gas_adaptive(
        &inputs,
        &original,
        &source,
        duration_s,
        loose_controls(duration_s),
        || false,
    )
    .unwrap();
    let controls = time_controls(duration_s);
    let tight = advance_droplets_with_gas_adaptive(
        &inputs,
        &original,
        &source,
        duration_s,
        controls,
        || false,
    )
    .unwrap();
    let reference = fixed_reference(&inputs, &original, &source, duration_s, 512);
    let loose_error = temperature_error(loose.endpoint(), &reference);
    let tight_error = temperature_error(tight.endpoint(), &reference);
    assert!(
        tight_error < loose_error * 0.25,
        "tight={tight_error:e}, loose={loose_error:e}"
    );
    assert_eq!(loose.report().accepted_comparisons, 1);
    assert!(tight.report().accuracy_rejections > 0);
    assert!(tight.report().accepted_intervals > loose.report().accepted_intervals);
    assert_eq!(
        tight.report().accepted_intervals,
        2 * tight.report().accepted_comparisons
    );
    assert!(tight.report().attempted_steps > 1 + 3 * tight.report().accepted_comparisons);
    assert!(tight.report().maximum_accepted_normalized_error <= 1.0);
    assert_eq!(tight.report().duration_s, duration_s);
    assert_eq!(tight.gas().volume_m3(), original.volume_m3());
    assert!(std::ptr::eq(tight.catalog(), &source));

    let audit = &tight.report().exchange;
    let mut liquid_mass_change = 0.0;
    let mut liquid_energy_change = 0.0;
    for (step, account) in tight.droplets().iter().zip(&audit.droplets) {
        let input = inputs.iter().find(|input| input.id == step.id).unwrap();
        assert_eq!(step.id, account.id);
        let initial_geometry = input
            .droplet
            .geometry(
                input.law,
                temperature(input.droplet.specific_enthalpy_j_kg).unwrap(),
            )
            .unwrap();
        assert!(step.geometry[1].radius_m() < initial_geometry.radius_m());
        // The LAST accepted fine interval already sees an evolved sphere.
        assert!(step.geometry[0].radius_m() < initial_geometry.radius_m());
        assert!(
            step.mass_transfer_coefficient_m_s
                > initial_geometry
                    .mass_transfer_coefficient_at_sherwood_m_s(
                        input.binary_diffusivity_m2_s,
                        input.sherwood_number
                    )
                    .unwrap()
        );
        assert_eq!(
            account.mass_kg,
            [input.droplet.mass_kg, step.droplet.mass_kg]
        );
        let dm = step.droplet.mass_kg - input.droplet.mass_kg;
        near(
            dm + account.transferred_mass_kg,
            0.0,
            controls.mass_conservation_tolerance_kg,
        );
        liquid_mass_change += dm;
        liquid_energy_change += step.droplet.mass_kg
            * (step.droplet.specific_enthalpy_j_kg + binding.liquid_enthalpy_offset_j_kg)
            - input.droplet.mass_kg
                * (input.droplet.specific_enthalpy_j_kg + binding.liquid_enthalpy_offset_j_kg);
    }
    let vapor = &binding.species;
    let carrier = SpeciesId::new("A").unwrap();
    near(
        liquid_mass_change + tight.gas().mass_kg() - original.mass_kg(),
        0.0,
        controls.mass_conservation_tolerance_kg,
    );
    near(
        liquid_mass_change + tight.gas().species_mass_kg(vapor).unwrap()
            - original.species_mass_kg(vapor).unwrap(),
        0.0,
        controls.mass_conservation_tolerance_kg,
    );
    near(
        tight.gas().species_mass_kg(&carrier).unwrap(),
        original.species_mass_kg(&carrier).unwrap(),
        controls.mass_conservation_tolerance_kg,
    );
    near(
        liquid_mass_change + audit.transferred_mass_kg,
        0.0,
        controls.mass_conservation_tolerance_kg,
    );
    assert!(audit.evaporated_mass_kg > 0.0);
    assert_eq!(audit.condensed_mass_kg, 0.0);
    near(audit.evaporated_mass_kg, audit.transferred_mass_kg, 1e-24);
    let species = audit
        .species
        .iter()
        .find(|row| &row.species == vapor)
        .unwrap();
    near(
        species.transferred_mass_kg,
        audit.transferred_mass_kg,
        1e-24,
    );
    near(
        species.closed_mass_residual_kg,
        0.0,
        controls.mass_conservation_tolerance_kg,
    );
    let gas_energy_change = tight.gas().mass_kg()
        * tight.gas().state().specific_internal_energy_j_kg()
        - original.mass_kg() * original.state().specific_internal_energy_j_kg();
    near(
        liquid_energy_change + gas_energy_change,
        0.0,
        controls.energy_conservation_tolerance_j,
    );
    near(
        audit.liquid_enthalpy_change_j + audit.gas_internal_energy_change_j,
        0.0,
        controls.energy_conservation_tolerance_j,
    );
    assert!(audit.sum_absolute_mass_residual_kg <= controls.mass_conservation_tolerance_kg);
    assert!(audit.sum_absolute_energy_residual_j <= controls.energy_conservation_tolerance_j);
    assert_eq!(original, saved);
}

#[test]
fn reciprocal_heat_is_time_resolved_when_mass_is_sealed_and_enthalpy_datum_changes() {
    let reference = -1e7;
    let source = catalog(reference);
    let original = gas(&source, 320.0, 1000.0, 1e-4);
    let saved = original.clone();
    let mut controls = time_controls(0.8);
    // Let enthalpy control the comparison. A large arbitrary h datum must not
    // enlarge a relative tolerance and hide the reciprocal thermal transient.
    controls.temperature_absolute_tolerance_k = 1000.0;
    controls.radius_absolute_tolerance_m = 1.0;
    controls.mass_absolute_tolerance_kg = 1e-4;
    controls.relative_tolerance = 1e-4;
    let run = |datum: f64| {
        let binding = GasLiquidBinding {
            liquid_enthalpy_offset_j_kg: reference - datum,
            ..binding(reference)
        };
        let inverse = |h| temperature(h - datum);
        let heat_controls = GasFilmHeatSolve {
            pure_liquid_enthalpy_range_j_kg: [
                liquid(200.0).specific_enthalpy_j_kg + datum,
                liquid(900.0).specific_enthalpy_j_kg + datum,
            ],
            temperature_tolerance_k: 1e-8,
            maximum_iterations: 128,
        };
        let inputs =
            [row(9, &binding, 1e-6, 280.0), row(2, &binding, 8e-6, 340.0)].map(|surface| {
                GasDropletHeatSurface {
                    surface: GasDropletSurface {
                        droplet: SphericalDropletState {
                            specific_enthalpy_j_kg: surface.droplet.specific_enthalpy_j_kg + datum,
                            ..surface.droplet
                        },
                        temperature_from_enthalpy: &inverse,
                        binary_diffusivity_m2_s: 0.0,
                        ..surface
                    },
                    heat_transfer_coefficient_w_m2_k: 80.0,
                    heat_controls,
                }
            });
        advance_droplets_with_gas_and_heat_adaptive(
            &inputs,
            &original,
            &source,
            0.8,
            controls,
            GasLiquidChamberHeatSolve {
                maximum_iterations: 128,
            },
            || false,
        )
        .unwrap()
    };
    let normal = run(0.0);
    let shifted = run(1e8);
    assert!(normal.report().accuracy_rejections > 0);
    assert!(normal.report().accepted_comparisons > 1);
    assert!(shifted.report().accepted_comparisons > 1);
    for (normal, shifted) in normal.droplets().iter().zip(shifted.droplets()) {
        assert_eq!(normal.id, shifted.id);
        assert_eq!(normal.droplet.mass_kg, shifted.droplet.mass_kg);
        assert_eq!(normal.geometry[1], shifted.geometry[1]);
        near(
            normal.exchange.surface_temperature_k[1],
            shifted.exchange.surface_temperature_k[1],
            2e-6,
        );
        near(
            normal.droplet.specific_enthalpy_j_kg,
            shifted.droplet.specific_enthalpy_j_kg - 1e8,
            0.005,
        );
    }
    near(
        normal.gas().state().state().temperature,
        shifted.gas().state().state().temperature,
        2e-6,
    );
    for candidate in [&normal, &shifted] {
        let audit = &candidate.report().exchange;
        assert_eq!(audit.transferred_mass_kg, 0.0);
        assert_eq!(candidate.gas().mass_kg(), original.mass_kg());
        assert_eq!(audit.energy_into_gas_j, 0.0);
        assert_ne!(audit.sensible_heat_into_liquid_j, 0.0);
        near(
            audit.liquid_enthalpy_change_j + audit.gas_internal_energy_change_j,
            0.0,
            controls.energy_conservation_tolerance_j,
        );
        assert!(audit.sum_absolute_energy_residual_j <= controls.energy_conservation_tolerance_j);
        assert!(candidate.report().maximum_accepted_normalized_error <= 1.0);
    }
    assert_eq!(original, saved);
}

#[test]
fn canonical_droplets_and_resumed_windows_retain_the_same_physical_trajectory() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let original = gas(&source, 320.0, 200.0, 1e-5);
    let inputs = rows(&binding);
    let controls = AdaptiveGasDropletControls {
        minimum_dt_s: 0.125,
        ..loose_controls(0.125)
    };
    let whole =
        advance_droplets_with_gas_adaptive(&inputs, &original, &source, 0.5, controls, || false)
            .unwrap();
    let permuted = advance_droplets_with_gas_adaptive(
        &[inputs[1], inputs[0]],
        &original,
        &source,
        0.5,
        controls,
        || false,
    )
    .unwrap();
    assert_eq!(whole.droplets(), permuted.droplets());
    assert_eq!(whole.gas(), permuted.gas());
    assert_eq!(
        whole.report().accepted_intervals,
        permuted.report().accepted_intervals
    );
    let first =
        advance_droplets_with_gas_adaptive(&inputs, &original, &source, 0.25, controls, || false)
            .unwrap();
    let resumed: Vec<_> = first
        .droplets()
        .iter()
        .map(|step| GasDropletSurface {
            droplet: step.droplet,
            ..*inputs.iter().find(|row| row.id == step.id).unwrap()
        })
        .collect();
    let second =
        advance_droplets_with_gas_adaptive(&resumed, first.gas(), &source, 0.25, controls, || {
            false
        })
        .unwrap();
    assert_eq!(whole.droplets(), second.droplets());
    assert_eq!(whole.gas(), second.gas());
    assert_eq!(
        whole.report().accepted_intervals,
        first.report().accepted_intervals + second.report().accepted_intervals
    );
    near(
        whole.report().exchange.transferred_mass_kg,
        first.report().exchange.transferred_mass_kg + second.report().exchange.transferred_mass_kg,
        1e-24,
    );
    near(
        whole.report().exchange.energy_into_gas_j,
        first.report().exchange.energy_into_gas_j + second.report().exchange.energy_into_gas_j,
        1e-16,
    );
    near(
        whole.report().exchange.latent_energy_j,
        first.report().exchange.latent_energy_j + second.report().exchange.latent_energy_j,
        1e-16,
    );
}

#[test]
fn physical_refinement_work_limits_and_late_cancellation_leave_original_states_retryable() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let original = gas(&source, 320.0, 200.0, 1e-5);
    let saved = original.clone();
    let inputs = rows(&binding);
    let initial_droplets = inputs.map(|input| input.droplet);
    let mut physical = loose_controls(0.4);
    physical.solve.surface.maximum_inventory_fraction = 0.02;
    let refined =
        advance_droplets_with_gas_adaptive(&inputs, &original, &source, 0.4, physical, || false)
            .unwrap();
    assert!(refined.report().physical_rejections > 0);
    assert!(refined.report().accepted_intervals > 2);
    assert_eq!(refined.report().duration_s, 0.4);

    let fixed = AdaptiveGasDropletControls {
        minimum_dt_s: 0.125,
        ..loose_controls(0.125)
    };
    assert!(matches!(
        advance_droplets_with_gas_adaptive(
            &inputs,
            &original,
            &source,
            0.5,
            AdaptiveGasDropletControls {
                maximum_comparisons: 1,
                ..fixed
            },
            || false
        ),
        Err(AdaptiveGasDropletError::ComparisonBudget)
    ));
    let strict = AdaptiveGasDropletControls {
        minimum_dt_s: 0.125,
        enthalpy_absolute_tolerance_j_kg: 1e-30,
        temperature_absolute_tolerance_k: 1e-30,
        mass_absolute_tolerance_kg: 1e-30,
        radius_absolute_tolerance_m: 1e-30,
        ..time_controls(0.125)
    };
    assert!(matches!(
        advance_droplets_with_gas_adaptive(&inputs, &original, &source, 0.125, strict, || false),
        Err(AdaptiveGasDropletError::TimeAccuracy { .. })
    ));
    let mut calls = 0;
    let accepted =
        advance_droplets_with_gas_adaptive(&inputs, &original, &source, 0.25, fixed, || {
            calls += 1;
            false
        })
        .unwrap();
    let last = calls;
    let mut calls = 0;
    let cancelled =
        advance_droplets_with_gas_adaptive(&inputs, &original, &source, 0.25, fixed, || {
            calls += 1;
            calls == last
        });
    assert!(matches!(
        cancelled,
        Err(AdaptiveGasDropletError::Material(
            GasEvaporationError::Cancelled
        )) | Err(AdaptiveGasDropletError::Step {
            error: GasEvaporationError::Cancelled,
            ..
        })
    ));
    assert_eq!(original, saved);
    assert_eq!(inputs.map(|input| input.droplet), initial_droplets);
    let retry =
        advance_droplets_with_gas_adaptive(&inputs, &original, &source, 0.25, fixed, || false)
            .unwrap();
    assert_eq!(retry.droplets(), accepted.droplets());
    assert_eq!(retry.gas(), accepted.gas());
    let zero =
        advance_droplets_with_gas_adaptive(&inputs, &original, &source, 0.0, fixed, || false)
            .unwrap();
    assert_eq!(zero.gas(), &original);
    assert_eq!(zero.report().attempted_steps, 1);
    assert_eq!(zero.report().accepted_intervals, 0);
    assert_eq!(zero.report().exchange.transferred_mass_kg, 0.0);
    assert_eq!(original, saved);
}
