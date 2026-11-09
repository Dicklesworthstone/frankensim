//! G0/G3/G4 original sourced trajectories, accepted fine accounts and refusal.
use super::*;
use crate::graph::thermal::medium::mass_transfer::evaporation::droplet::chamber::adaptive::{
    AdaptiveGasDropletControls, SourceAdaptiveGasDropletError, advance_droplets_with_gas_adaptive,
    advance_droplets_with_gas_and_heat_adaptive,
};

fn adaptive_controls() -> AdaptiveGasDropletControls {
    AdaptiveGasDropletControls {
        solve: chamber_controls(),
        initial_dt_s: 0.2,
        minimum_dt_s: 1e-8,
        maximum_dt_s: 0.2,
        maximum_comparisons: 4096,
        maximum_reductions: 32,
        enthalpy_absolute_tolerance_j_kg: 1.0,
        temperature_absolute_tolerance_k: 1e-3,
        mass_absolute_tolerance_kg: 1e-13,
        radius_absolute_tolerance_m: 1e-8,
        relative_tolerance: 1e-5,
        mass_conservation_tolerance_kg: 1e-11,
        energy_conservation_tolerance_j: 1e-5,
    }
}
fn small_gas(catalog: &ThermochemicalGasModel, t: f64) -> ClosedGasParcel {
    let original = gas(catalog, t, 100_000.0, 0.003);
    ClosedGasParcel::new(original.model().clone(), original.point().clone(), 1e-5).unwrap()
}

#[test]
fn accepted_fine_pair_refreshes_both_original_sources_and_owns_only_its_physical_accounts() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let inverse = |h| inverse(&caloric, h);
    let binding = liquid_binding();
    let original = small_gas(&catalog, 430.0);
    let saved = original.clone();
    let card = diffusion();
    let context = point(400.0, 100_000.0);
    let pinned = SourceSphericalDropletFilm::new(
        &source(&card, &context).source,
        &binding.species,
        &original,
        || false,
    )
    .unwrap();
    let selected = pinned.selected_diffusivity().clone();
    let forced = SourceForcedSphericalDropletFilm::try_new(&pinned, MODEL, 0.5, 0.04).unwrap();
    let rows = [
        SourceGasDropletSurface {
            id: 9,
            law: &LiquidPhase,
            binding: &binding,
            droplet: sphere(&caloric, 390.0, 1e-6),
            temperature_from_enthalpy: &inverse,
            film: SourceDropletFilm::Stagnant(&pinned),
        },
        SourceGasDropletSurface {
            id: 2,
            law: &LiquidPhase,
            binding: &binding,
            droplet: sphere(&caloric, 400.0, 8e-6),
            temperature_from_enthalpy: &inverse,
            film: SourceDropletFilm::Forced(forced),
        },
    ];
    let duration = 0.2;
    // A deliberately loose local estimate admits this one comparison. Its
    // exact original fine trajectory is an independent composition oracle.
    let controls = AdaptiveGasDropletControls {
        enthalpy_absolute_tolerance_j_kg: 1e4,
        temperature_absolute_tolerance_k: 10.0,
        mass_absolute_tolerance_kg: 1e-7,
        radius_absolute_tolerance_m: 1e-4,
        ..adaptive_controls()
    };
    for with_heat in [false, true] {
        let warm = rows.map(|surface| SourceGasDropletHeatSurface {
            surface,
            heat_transfer_coefficient_w_m2_k: 80.0,
            heat_controls: heat_controls(&caloric),
        });
        let original_step = |rows: &[SourceGasDropletSurface<'_, '_>], gas, dt| {
            if with_heat {
                let warm: Vec<_> = rows
                    .iter()
                    .map(|&surface| SourceGasDropletHeatSurface {
                        surface,
                        heat_transfer_coefficient_w_m2_k: 80.0,
                        heat_controls: heat_controls(&caloric),
                    })
                    .collect();
                advance_droplets_with_gas_and_heat(
                    &warm,
                    gas,
                    &catalog,
                    dt,
                    controls.solve,
                    chamber_heat_controls(),
                    || false,
                )
            } else {
                advance_droplets_with_gas(rows, gas, &catalog, dt, controls.solve, || false)
            }
            .unwrap()
        };
        let half = original_step(&rows, &original, duration / 2.0);
        let next_rows: Vec<_> = rows
            .iter()
            .map(|&surface| SourceGasDropletSurface {
                droplet: half
                    .original()
                    .droplets()
                    .iter()
                    .find(|step| step.id == surface.id)
                    .unwrap()
                    .droplet,
                ..surface
            })
            .collect();
        let fine = original_step(&next_rows, half.original().gas(), duration / 2.0);
        let next = if with_heat {
            advance_droplets_with_gas_and_heat_adaptive(
                &warm,
                &original,
                &catalog,
                duration,
                controls,
                chamber_heat_controls(),
                || false,
            )
        } else {
            advance_droplets_with_gas_adaptive(
                &rows,
                &original,
                &catalog,
                duration,
                controls,
                || false,
            )
        }
        .unwrap();
        let accepted = next.original();
        assert_eq!(accepted.endpoint(), fine.original());
        assert_eq!(next.films(), fine.films());
        assert_eq!(accepted.report().duration_s, duration);
        assert_eq!(accepted.report().comparisons, 1);
        assert_eq!(accepted.report().attempted_steps, 4);
        assert_eq!(accepted.report().accepted_intervals, 2);
        assert_eq!(accepted.report().accuracy_rejections, 0);
        assert_eq!(accepted.report().physical_rejections, 0);
        assert!(accepted.report().maximum_accepted_normalized_error <= 1.0);
        assert!(std::ptr::eq(accepted.catalog(), &catalog));
        assert_eq!(accepted.gas().volume_m3(), original.volume_m3());
        near(
            accepted.report().exchange.transferred_mass_kg,
            half.original().step().transferred_mass_kg + fine.original().step().transferred_mass_kg,
            1e-22,
            1e-13,
        );
        near(
            accepted.report().exchange.energy_into_gas_j,
            half.original().step().energy_into_gas_j + fine.original().step().energy_into_gas_j,
            1e-13,
            1e-13,
        );
        for ((endpoint, last), first) in accepted
            .droplets()
            .iter()
            .zip(next.films())
            .zip(half.films())
        {
            assert_eq!(endpoint.id, last.id);
            assert_ne!(last.films[0].geometry(), first.films[0].geometry());
            assert_ne!(
                last.films[0].binary_diffusivity_m2_s(),
                first.films[0].binary_diffusivity_m2_s()
            );
            let actual_final = match &last.films[1] {
                SourceDropletFilmSample::Stagnant(_) => SourceDropletFilmSample::Stagnant(
                    pinned
                        .resolve(endpoint.geometry[1], accepted.gas(), || false)
                        .unwrap(),
                ),
                SourceDropletFilmSample::Forced(_) => SourceDropletFilmSample::Forced(
                    forced
                        .resolve(
                            endpoint.geometry[1],
                            endpoint.exchange.surface_temperature_k[1],
                            accepted.gas(),
                            || false,
                        )
                        .unwrap(),
                ),
            };
            assert_eq!(last.films[1], actual_final);
        }
        let mut dm = 0.0;
        let mut dh = 0.0;
        for result in accepted.droplets() {
            let input = rows.iter().find(|row| row.id == result.id).unwrap();
            dm += result.droplet.mass_kg - input.droplet.mass_kg;
            dh += result.droplet.mass_kg
                * (result.droplet.specific_enthalpy_j_kg + binding.liquid_enthalpy_offset_j_kg)
                - input.droplet.mass_kg
                    * (input.droplet.specific_enthalpy_j_kg + binding.liquid_enthalpy_offset_j_kg);
        }
        near(
            dm + accepted.gas().mass_kg() - original.mass_kg(),
            0.0,
            controls.mass_conservation_tolerance_kg,
            0.0,
        );
        near(
            dm + accepted.gas().species_mass_kg(&binding.species).unwrap()
                - original.species_mass_kg(&binding.species).unwrap(),
            0.0,
            controls.mass_conservation_tolerance_kg,
            0.0,
        );
        near(
            dh + accepted.gas().mass_kg() * accepted.gas().state().specific_internal_energy_j_kg()
                - original.mass_kg() * original.state().specific_internal_energy_j_kg(),
            0.0,
            controls.energy_conservation_tolerance_j,
            0.0,
        );
        assert_eq!(pinned.selected_diffusivity(), &selected);
        assert_eq!(original, saved);
    }
}

#[test]
fn heated_source_trajectory_refines_preserves_order_and_cancels_without_advancing_inputs() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let inverse = |h| inverse(&caloric, h);
    let binding = liquid_binding();
    let original = small_gas(&catalog, 430.0);
    let saved = original.clone();
    let card = diffusion();
    let context = point(400.0, 100_000.0);
    let pinned = SourceSphericalDropletFilm::new(
        &source(&card, &context).source,
        &binding.species,
        &original,
        || false,
    )
    .unwrap();
    let forced = SourceForcedSphericalDropletFilm::try_new(&pinned, MODEL, 0.5, 0.04).unwrap();
    let rows =
        [(9, 390.0, 1e-6), (2, 400.0, 8e-6)].map(|(id, t, mass)| SourceGasDropletHeatSurface {
            surface: SourceGasDropletSurface {
                id,
                law: &LiquidPhase,
                binding: &binding,
                droplet: sphere(&caloric, t, mass),
                temperature_from_enthalpy: &inverse,
                film: if id == 9 {
                    SourceDropletFilm::Stagnant(&pinned)
                } else {
                    SourceDropletFilm::Forced(forced)
                },
            },
            heat_transfer_coefficient_w_m2_k: 80.0,
            heat_controls: heat_controls(&caloric),
        });
    let original_states = rows.map(|row| row.surface.droplet);
    let controls = AdaptiveGasDropletControls {
        initial_dt_s: 0.4,
        maximum_dt_s: 0.4,
        ..adaptive_controls()
    };
    let mut calls = 0usize;
    let next = advance_droplets_with_gas_and_heat_adaptive(
        &rows,
        &original,
        &catalog,
        0.4,
        controls,
        chamber_heat_controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    assert!(next.original().report().accuracy_rejections > 0);
    assert!(next.original().report().accepted_intervals > 2);
    assert!(next.original().report().maximum_accepted_normalized_error <= 1.0);
    let permuted = advance_droplets_with_gas_and_heat_adaptive(
        &[rows[1], rows[0]],
        &original,
        &catalog,
        0.4,
        controls,
        chamber_heat_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(next.original().endpoint(), permuted.original().endpoint());
    assert_eq!(next.films(), permuted.films());
    assert_eq!(next.original().report(), permuted.original().report());
    let last_gate = calls;
    let mut calls = 0usize;
    assert!(
        advance_droplets_with_gas_and_heat_adaptive(
            &rows,
            &original,
            &catalog,
            0.4,
            controls,
            chamber_heat_controls(),
            || {
                calls += 1;
                calls == last_gate
            },
        )
        .is_err()
    );
    let retry = advance_droplets_with_gas_and_heat_adaptive(
        &rows,
        &original,
        &catalog,
        0.4,
        controls,
        chamber_heat_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(next.original().endpoint(), retry.original().endpoint());
    assert_eq!(next.films(), retry.films());
    assert_eq!(rows.map(|row| row.surface.droplet), original_states);
    assert_eq!(original, saved);
}

#[test]
fn zero_time_admits_original_sources_and_late_source_failure_is_not_an_accuracy_retry() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let inverse = |h| inverse(&caloric, h);
    let binding = liquid_binding();
    let original = gas(&catalog, 430.0, 100_000.0, 0.003);
    let saved = original.clone();
    let seed = gas(&catalog, 400.0, 100_000.0, 0.003);
    let context = point(400.0, 100_000.0);
    let card = material(vec![scalar(
        BINARY_DIFFUSIVITY_PROPERTY,
        BINARY_DIFFUSIVITY_DIMS,
        "temperature",
        vec![(200.0, 1.6e-5), (400.0001, 1.6e-5)],
        [200.0, 400.0001],
    )]);
    let pinned = SourceSphericalDropletFilm::new(
        &source(&card, &context).source,
        &binding.species,
        &seed,
        || false,
    )
    .unwrap();
    let selected = pinned.selected_diffusivity().clone();
    let forced = SourceForcedSphericalDropletFilm::try_new(&pinned, MODEL, 0.5, 0.04).unwrap();
    let sealed = KineticEvaporation::try_new(
        &LiquidPhase,
        InterfaceKinetics::HertzKnudsen {
            accommodation_coefficient: 0.0,
        },
    )
    .unwrap();
    let row = SourceGasDropletSurface {
        id: 2,
        law: &sealed,
        binding: &binding,
        droplet: sphere(&caloric, 370.0, 1e-6),
        temperature_from_enthalpy: &inverse,
        film: SourceDropletFilm::Forced(forced),
    };
    let direct =
        advance_droplets_with_gas(&[row], &original, &catalog, 0.0, chamber_controls(), || {
            false
        })
        .unwrap();
    let zero = advance_droplets_with_gas_adaptive(
        &[row],
        &original,
        &catalog,
        0.0,
        adaptive_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(zero.original().endpoint(), direct.original());
    assert_eq!(zero.films(), direct.films());
    assert_eq!(zero.original().report().attempted_steps, 1);
    assert_eq!(zero.original().report().accepted_intervals, 0);
    assert_eq!(zero.original().gas(), &original);
    let no_query = |_| -> Result<f64, EvaporationError> {
        panic!("the pinned volatile identity must precede the liquid query")
    };
    let wrong = GasLiquidBinding {
        species: SpeciesId::new("A").unwrap(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let invalid = SourceGasDropletSurface {
        binding: &wrong,
        temperature_from_enthalpy: &no_query,
        ..row
    };
    assert!(matches!(
        advance_droplets_with_gas_adaptive(
            &[invalid],
            &original,
            &catalog,
            0.0,
            adaptive_controls(),
            || false,
        ),
        Err(SourceAdaptiveGasDropletError::Step {
            time_s: 0.0,
            duration_s: 0.0,
            error: SphericalDropletFilmError::Film(GasFilmError::Invalid(
                "spherical film cannot change its pinned volatile species"
            )),
        })
    ));
    // The stagnant source would sample the unsupported actual bulk T=430 K,
    // so even a request with no physical transfer must retain that refusal.
    let unsupported = SourceGasDropletSurface {
        film: SourceDropletFilm::Stagnant(&pinned),
        ..row
    };
    assert!(matches!(
        advance_droplets_with_gas_adaptive(
            &[unsupported],
            &original,
            &catalog,
            0.0,
            adaptive_controls(),
            || false,
        ),
        Err(SourceAdaptiveGasDropletError::Step {
            error: SphericalDropletFilmError::Film(GasFilmError::Source(_)),
            ..
        })
    ));
    let warm = SourceGasDropletHeatSurface {
        surface: row,
        heat_transfer_coefficient_w_m2_k: 80.0,
        heat_controls: heat_controls(&caloric),
    };
    assert!(matches!(
        advance_droplets_with_gas_and_heat_adaptive(
            &[warm],
            &original,
            &catalog,
            0.1,
            adaptive_controls(),
            chamber_heat_controls(),
            || false,
        ),
        Err(SourceAdaptiveGasDropletError::Step {
            time_s: 0.0,
            duration_s: 0.1,
            error: SphericalDropletFilmError::Film(GasFilmError::Source(_)),
        })
    ));
    assert_eq!(pinned.selected_diffusivity(), &selected);
    assert_eq!(original, saved);
}
