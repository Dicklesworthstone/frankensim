//! Shared original gas, heterogeneous sourced films and all-or-nothing support.
use super::*;
use crate::graph::thermal::medium::mass_transfer::evaporation::droplet::chamber::{
    SourceDropletFilm, SourceDropletFilmSample, SourceGasDropletHeatSurface,
    SourceGasDropletSurface, advance_droplets_with_gas, advance_droplets_with_gas_and_heat,
};
use crate::graph::thermal::medium::mass_transfer::evaporation::droplet::forced::SourceForcedSphericalDropletFilm;
use crate::graph::thermal::medium::mass_transfer::sphere::ForcedSphereMassTransfer;
use fs_material::evaporation::gas::GasEvaporationError;
use fs_material::evaporation::gas::chamber::{
    GasLiquidChamberSolve, heat::GasLiquidChamberHeatSolve,
};
use fs_material::evaporation::gas::droplet::chamber::{
    GasDropletHeatSurface, GasDropletSurface,
    advance_droplets_with_gas_and_heat as advance_material_with_heat,
};
use fs_thermal::evaporation::surface_flux;

const MODEL: ForcedSphereMassTransfer = ForcedSphereMassTransfer::Frossling1938;

fn chamber_controls() -> GasLiquidChamberSolve {
    GasLiquidChamberSolve {
        surface: phase_controls().phase,
        maximum_surfaces: 8,
    }
}
fn chamber_heat_controls() -> GasLiquidChamberHeatSolve {
    GasLiquidChamberHeatSolve {
        maximum_iterations: 128,
    }
}

#[test]
fn mixed_sources_share_the_same_gas_and_are_invariant_to_input_order_in_both_directions() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let inverse = |h| inverse(&caloric, h);
    let binding = liquid_binding();
    let card = diffusion();
    let context = point(250.0, 50_000.0);
    let pair = source(&card, &context).source;
    let kinetic = KineticEvaporation::try_new(
        &LiquidPhase,
        InterfaceKinetics::Schrage {
            accommodation_coefficient: 1e-5,
        },
    )
    .unwrap();
    for (gas_t, vapor) in [(420.0, 0.0), (430.0, 0.02)] {
        let original = gas(&catalog, gas_t, 100_000.0, vapor);
        let saved = original.clone();
        let pinned =
            SourceSphericalDropletFilm::new(&pair, &binding.species, &original, || false).unwrap();
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
                law: &kinetic,
                binding: &binding,
                droplet: sphere(&caloric, 400.0, 8e-6),
                temperature_from_enthalpy: &inverse,
                film: SourceDropletFilm::Forced(forced),
            },
        ];
        let dt = 0.1;
        let next =
            advance_droplets_with_gas(&rows, &original, &catalog, dt, chamber_controls(), || false)
                .unwrap();
        let permuted = advance_droplets_with_gas(
            &[rows[1], rows[0]],
            &original,
            &catalog,
            dt,
            chamber_controls(),
            || false,
        )
        .unwrap();
        assert_eq!(next, permuted);
        assert_eq!(
            next.films().iter().map(|film| film.id).collect::<Vec<_>>(),
            [2, 9]
        );
        let actual = next.original();
        assert_eq!(actual.droplets().len(), rows.len());
        assert_eq!(actual.step().transferred_mass_kg > 0.0, vapor == 0.0);
        for (result, films) in actual.droplets().iter().zip(next.films()) {
            let input = rows.iter().find(|row| row.id == result.id).unwrap();
            let first = &films.films[0];
            let initial_t = inverse(input.droplet.specific_enthalpy_j_kg).unwrap();
            let film_t = match first {
                SourceDropletFilmSample::Stagnant(sample) => {
                    assert_eq!(sample.temperature_k, gas_t);
                    assert_eq!(sample.total_pressure_pa, original.state().state().pressure);
                    assert_eq!(first.sherwood_number(), 2.0);
                    sample.temperature_k
                }
                SourceDropletFilmSample::Forced(sample) => {
                    near(
                        sample.film_temperature_k,
                        0.5 * initial_t + 0.5 * gas_t,
                        1e-12,
                        0.0,
                    );
                    assert_eq!(sample.bulk_temperature_k, gas_t);
                    assert_eq!(sample.total_pressure_pa, original.state().state().pressure);
                    assert!(first.sherwood_number() > 2.0);
                    sample.film_temperature_k
                }
            };
            let expected = surface_flux(
                input.law,
                initial_t,
                GasFilm {
                    temperature_k: film_t,
                    total_pressure_pa: original.state().state().pressure,
                    vapor_partial_pressure_pa: original
                        .species_partial_pressure_pa(&binding.species)
                        .unwrap(),
                    mass_transfer_coefficient_m_s: first.mass_transfer_coefficient_m_s(),
                },
            )
            .unwrap();
            near(
                result.phase.requested_transferred_mass_kg,
                dt * first.geometry().area_m2() * expected.mass_flux_kg_m2_s,
                1e-22,
                2e-13,
            );
            assert_eq!(
                result.binary_diffusivity_m2_s,
                first.binary_diffusivity_m2_s()
            );
            assert_eq!(result.sherwood_number, first.sherwood_number());
            assert_eq!(
                result.mass_transfer_coefficient_m_s,
                first.mass_transfer_coefficient_m_s()
            );
            assert_eq!(films.films[0].geometry(), result.geometry[0]);
            assert_eq!(films.films[1].geometry(), result.geometry[1]);
            assert_eq!(result.droplet.mass_kg < input.droplet.mass_kg, vapor == 0.0);
            let refreshed = match &films.films[1] {
                SourceDropletFilmSample::Stagnant(sample) => {
                    let direct = pinned
                        .resolve(result.geometry[1], actual.gas(), || false)
                        .unwrap();
                    assert_eq!(sample, &direct);
                    &sample.diffusion
                }
                SourceDropletFilmSample::Forced(sample) => {
                    let direct = forced
                        .resolve(
                            result.geometry[1],
                            result.exchange.surface_temperature_k[1],
                            actual.gas(),
                            || false,
                        )
                        .unwrap();
                    assert_eq!(sample, &direct);
                    &sample.diffusion
                }
            };
            assert_eq!(refreshed.species, pair.species);
            assert!(
                !refreshed
                    .material
                    .property(BINARY_DIFFUSIVITY_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .source_hashes
                    .is_empty()
            );
        }
        let final_liquid_mass: f64 = actual
            .droplets()
            .iter()
            .map(|row| row.droplet.mass_kg)
            .sum();
        near(
            final_liquid_mass + actual.gas().species_mass_kg(&binding.species).unwrap(),
            rows.iter().map(|row| row.droplet.mass_kg).sum::<f64>()
                + original.species_mass_kg(&binding.species).unwrap(),
            phase_controls().phase.gas.mass_tolerance_kg,
            0.0,
        );
        near(
            actual.step().liquid_enthalpy_change_j + actual.step().gas_internal_energy_change_j,
            0.0,
            phase_controls().phase.gas.energy_tolerance_j,
            0.0,
        );
        assert_eq!(actual.gas().volume_m3(), original.volume_m3());
        assert_eq!(actual.gas().point().axes().get("local-lot"), Some(&9.0));
        assert!(std::ptr::eq(actual.catalog(), &catalog));
        assert!(std::ptr::eq(pinned.diffusion_source().material, &card));
        assert_eq!(pinned.selected_diffusivity(), &selected);
        assert_eq!(original, saved);
    }
}

#[test]
fn singleton_keeps_the_original_stagnant_and_forced_candidates_and_receipts() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let inverse = |h| inverse(&caloric, h);
    let binding = liquid_binding();
    let original = gas(&catalog, 430.0, 100_000.0, 0.02);
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
    let initial = sphere(&caloric, 400.0, 1e-6);
    for film in [
        SourceDropletFilm::Stagnant(&pinned),
        SourceDropletFilm::Forced(forced),
    ] {
        let row = SourceGasDropletSurface {
            id: 17,
            law: &LiquidPhase,
            binding: &binding,
            droplet: initial,
            temperature_from_enthalpy: &inverse,
            film,
        };
        let next =
            advance_droplets_with_gas(&[row], &original, &catalog, 0.1, chamber_controls(), || {
                false
            })
            .unwrap();
        let (direct, samples) = match film {
            SourceDropletFilm::Stagnant(source) => {
                let (direct, samples) = source
                    .advance_with_gas(
                        &LiquidPhase,
                        &binding,
                        initial,
                        &original,
                        &catalog,
                        inverse,
                        0.1,
                        phase_controls().phase,
                        || false,
                    )
                    .unwrap()
                    .into_parts();
                (direct, samples.map(SourceDropletFilmSample::Stagnant))
            }
            SourceDropletFilm::Forced(source) => {
                let (direct, samples) = source
                    .advance_with_gas(
                        &LiquidPhase,
                        &binding,
                        initial,
                        &original,
                        &catalog,
                        inverse,
                        0.1,
                        phase_controls().phase,
                        || false,
                    )
                    .unwrap()
                    .into_parts();
                (direct, samples.map(SourceDropletFilmSample::Forced))
            }
        };
        let result = &next.original().droplets()[0];
        assert_eq!(result.droplet, direct.droplet());
        assert_eq!(&result.geometry, direct.geometry());
        assert_eq!(next.original().gas(), direct.gas());
        assert_eq!(
            result.exchange.transferred_mass_kg,
            direct.step().transferred_mass_kg
        );
        assert_eq!(next.films()[0].films, samples);
    }
}

#[test]
fn reciprocal_heat_uses_one_finite_gas_and_refreshes_both_selected_films() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let inverse = |h| inverse(&caloric, h);
    let binding = liquid_binding();
    let original = gas(&catalog, 430.0, 100_000.0, 0.003);
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
    let sealed = KineticEvaporation::try_new(
        &LiquidPhase,
        InterfaceKinetics::HertzKnudsen {
            accommodation_coefficient: 0.0,
        },
    )
    .unwrap();
    let rows = [
        SourceGasDropletHeatSurface {
            surface: SourceGasDropletSurface {
                id: 8,
                law: &sealed,
                binding: &binding,
                droplet: sphere(&caloric, 350.0, 1e-6),
                temperature_from_enthalpy: &inverse,
                film: SourceDropletFilm::Stagnant(&pinned),
            },
            heat_transfer_coefficient_w_m2_k: 80.0,
            heat_controls: heat_controls(&caloric),
        },
        SourceGasDropletHeatSurface {
            surface: SourceGasDropletSurface {
                id: 3,
                law: &sealed,
                binding: &binding,
                droplet: sphere(&caloric, 370.0, 8e-6),
                temperature_from_enthalpy: &inverse,
                film: SourceDropletFilm::Forced(forced),
            },
            heat_transfer_coefficient_w_m2_k: 120.0,
            heat_controls: heat_controls(&caloric),
        },
    ];
    let next = advance_droplets_with_gas_and_heat(
        &rows,
        &original,
        &catalog,
        0.1,
        chamber_controls(),
        chamber_heat_controls(),
        || false,
    )
    .unwrap();
    let direct_rows: Vec<_> = rows
        .iter()
        .map(|row| {
            let used = &next
                .films()
                .iter()
                .find(|film| film.id == row.surface.id)
                .unwrap()
                .films[0];
            GasDropletHeatSurface {
                surface: GasDropletSurface {
                    id: row.surface.id,
                    law: row.surface.law,
                    binding: row.surface.binding,
                    droplet: row.surface.droplet,
                    temperature_from_enthalpy: row.surface.temperature_from_enthalpy,
                    binary_diffusivity_m2_s: used.binary_diffusivity_m2_s(),
                    sherwood_number: used.sherwood_number(),
                },
                heat_transfer_coefficient_w_m2_k: row.heat_transfer_coefficient_w_m2_k,
                heat_controls: row.heat_controls,
            }
        })
        .collect();
    let direct = advance_material_with_heat(
        &direct_rows,
        &original,
        &catalog,
        0.1,
        chamber_controls(),
        chamber_heat_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(next.original(), &direct);
    assert_eq!(direct.step().transferred_mass_kg, 0.0);
    assert!(direct.gas().state().state().temperature < original.state().state().temperature);
    for (result, film) in direct.droplets().iter().zip(next.films()) {
        assert!(result.heat.as_ref().unwrap().heat_into_liquid_j > 0.0);
        assert_eq!(
            result.droplet.mass_kg,
            rows.iter()
                .find(|row| row.surface.id == result.id)
                .unwrap()
                .surface
                .droplet
                .mass_kg
        );
        assert_ne!(
            film.films[0].binary_diffusivity_m2_s(),
            film.films[1].binary_diffusivity_m2_s()
        );
    }
    near(
        direct.step().liquid_enthalpy_change_j + direct.step().gas_internal_energy_change_j,
        0.0,
        phase_controls().phase.gas.energy_tolerance_j,
        0.0,
    );
    let permuted = advance_droplets_with_gas_and_heat(
        &[rows[1], rows[0]],
        &original,
        &catalog,
        0.1,
        chamber_controls(),
        chamber_heat_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(next, permuted);
    assert_eq!(original, saved);
}

#[test]
fn changed_identity_duplicate_ids_work_ceiling_and_cancellation_refuse_before_liquid_queries() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let binding = liquid_binding();
    let original = gas(&catalog, 400.0, 100_000.0, 0.003);
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
    let no_query =
        |_| -> Result<f64, EvaporationError> { panic!("admission must precede liquid query") };
    let wrong = GasLiquidBinding {
        species: SpeciesId::new("A").unwrap(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let row = SourceGasDropletSurface {
        id: 1,
        law: &LiquidPhase,
        binding: &wrong,
        droplet: sphere(&caloric, 400.0, 1e-6),
        temperature_from_enthalpy: &no_query,
        film: SourceDropletFilm::Stagnant(&pinned),
    };
    assert!(matches!(
        advance_droplets_with_gas(&[row], &original, &catalog, 0.1, chamber_controls(), || {
            false
        }),
        Err(SphericalDropletFilmError::Film(GasFilmError::Invalid(
            "spherical film cannot change its pinned volatile species"
        )))
    ));
    let row = SourceGasDropletSurface {
        binding: &binding,
        ..row
    };
    assert!(matches!(
        advance_droplets_with_gas(
            &[row, row],
            &original,
            &catalog,
            0.1,
            chamber_controls(),
            || false
        ),
        Err(SphericalDropletFilmError::Film(GasFilmError::Invalid(
            "sourced droplets require unique stable IDs"
        )))
    ));
    assert!(matches!(
        advance_droplets_with_gas(
            &[row],
            &original,
            &catalog,
            0.1,
            GasLiquidChamberSolve {
                maximum_surfaces: 0,
                ..chamber_controls()
            },
            || false
        ),
        Err(SphericalDropletFilmError::Film(GasFilmError::Invalid(_)))
    ));
    assert!(matches!(
        advance_droplets_with_gas(&[row], &original, &catalog, 0.1, chamber_controls(), || {
            true
        }),
        Err(SphericalDropletFilmError::Film(GasFilmError::Cancelled))
    ));
    assert_eq!(original, saved);
}

#[test]
fn one_late_forced_source_refusal_discards_the_complete_shared_heat_candidate() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let inverse = |h| inverse(&caloric, h);
    let binding = liquid_binding();
    let seed = gas(&catalog, 400.0, 100_000.0, 0.003);
    let original = gas(&catalog, 430.0, 100_000.0, 0.003);
    let saved = original.clone();
    let context = point(400.0, 100_000.0);
    let broad_card = diffusion();
    let broad = SourceSphericalDropletFilm::new(
        &source(&broad_card, &context).source,
        &binding.species,
        &original,
        || false,
    )
    .unwrap();
    let restricted_card = material(vec![
        scalar(
            BINARY_DIFFUSIVITY_PROPERTY,
            BINARY_DIFFUSIVITY_DIMS,
            "temperature",
            vec![(200.0, 1.6e-5), (400.0001, 1.6e-5)],
            [200.0, 400.0001],
        ),
        scalar(
            BINARY_DIFFUSIVITY_PROPERTY,
            BINARY_DIFFUSIVITY_DIMS,
            "temperature",
            vec![(400.0002, 2.2e-5), (900.0, 2.2e-5)],
            [400.0002, 900.0],
        ),
    ]);
    let restricted = SourceSphericalDropletFilm::new(
        &source(&restricted_card, &context).source,
        &binding.species,
        &seed,
        || false,
    )
    .unwrap();
    let selected = restricted.selected_diffusivity().clone();
    let forced = SourceForcedSphericalDropletFilm::try_new(&restricted, MODEL, 0.5, 0.04).unwrap();
    let sealed = KineticEvaporation::try_new(
        &LiquidPhase,
        InterfaceKinetics::HertzKnudsen {
            accommodation_coefficient: 0.0,
        },
    )
    .unwrap();
    let rows = [
        SourceGasDropletHeatSurface {
            surface: SourceGasDropletSurface {
                id: 1,
                law: &sealed,
                binding: &binding,
                droplet: sphere(&caloric, 400.0, 8e-6),
                temperature_from_enthalpy: &inverse,
                film: SourceDropletFilm::Stagnant(&broad),
            },
            heat_transfer_coefficient_w_m2_k: 40.0,
            heat_controls: heat_controls(&caloric),
        },
        SourceGasDropletHeatSurface {
            surface: SourceGasDropletSurface {
                id: 2,
                law: &sealed,
                binding: &binding,
                droplet: sphere(&caloric, 370.0, 1e-6),
                temperature_from_enthalpy: &inverse,
                film: SourceDropletFilm::Forced(forced),
            },
            heat_transfer_coefficient_w_m2_k: 80.0,
            heat_controls: heat_controls(&caloric),
        },
    ];
    let original_rows: Vec<_> = rows
        .iter()
        .map(|row| {
            let t = inverse(row.surface.droplet.specific_enthalpy_j_kg).unwrap();
            let geometry = row.surface.droplet.geometry(&sealed, t).unwrap();
            let (diffusivity, sherwood) = match row.surface.film {
                SourceDropletFilm::Stagnant(source) => {
                    let sample = source.resolve(geometry, &original, || false).unwrap();
                    (sample.diffusion.diffusivity_m2_s, 2.0)
                }
                SourceDropletFilm::Forced(source) => {
                    let sample = source.resolve(geometry, t, &original, || false).unwrap();
                    near(sample.film_temperature_k, 400.0, 1e-12, 0.0);
                    (
                        sample.diffusion.diffusivity_m2_s,
                        sample.correlation.sherwood_number(),
                    )
                }
            };
            GasDropletHeatSurface {
                surface: GasDropletSurface {
                    id: row.surface.id,
                    law: &sealed,
                    binding: &binding,
                    droplet: row.surface.droplet,
                    temperature_from_enthalpy: &inverse,
                    binary_diffusivity_m2_s: diffusivity,
                    sherwood_number: sherwood,
                },
                heat_transfer_coefficient_w_m2_k: row.heat_transfer_coefficient_w_m2_k,
                heat_controls: row.heat_controls,
            }
        })
        .collect();
    let physical = advance_material_with_heat(
        &original_rows,
        &original,
        &catalog,
        0.1,
        chamber_controls(),
        chamber_heat_controls(),
        || false,
    )
    .unwrap();
    let endpoint = &physical.droplets()[1];
    let final_mean = 0.5 * endpoint.exchange.surface_temperature_k[1]
        + 0.5 * physical.gas().state().state().temperature;
    assert!(final_mean > 400.0002);
    assert!(matches!(
        advance_droplets_with_gas_and_heat(
            &rows,
            &original,
            &catalog,
            0.1,
            chamber_controls(),
            chamber_heat_controls(),
            || false,
        ),
        Err(SphericalDropletFilmError::Film(GasFilmError::Source(_)))
    ));
    // The same original guard must run on the actual post-heat phase state,
    // even for a sealed interface. Initial/final source checks alone miss it.
    let initial_fraction = LiquidPhase
        .properties(370.0)
        .unwrap()
        .saturation_pressure_pa
        / original.state().state().pressure;
    let narrow = SourceForcedSphericalDropletFilm::try_new(
        &broad,
        MODEL,
        0.5,
        initial_fraction * (1.0 + 1e-5),
    )
    .unwrap();
    let mut narrowed = rows;
    narrowed[1].surface.film = SourceDropletFilm::Forced(narrow);
    assert!(matches!(
        advance_droplets_with_gas_and_heat(
            &narrowed,
            &original,
            &catalog,
            0.1,
            chamber_controls(),
            chamber_heat_controls(),
            || false,
        ),
        Err(SphericalDropletFilmError::Phase(
            GasEvaporationError::Surface(EvaporationError::Invalid(
                "forced sphere phase evaluation leaves the authored dilute vapor domain"
            ))
        ))
    ));
    assert_eq!(restricted.selected_diffusivity(), &selected);
    assert_eq!(original, saved);
}
