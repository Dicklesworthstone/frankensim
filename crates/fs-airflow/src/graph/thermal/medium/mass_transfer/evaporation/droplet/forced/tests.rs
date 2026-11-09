//! Original sources, bounded forced transport and unchanged gas caloric owners.
use super::*;
use crate::graph::thermal::medium::mass_transfer::evaporation::droplet::forced::{
    SourceForcedSphericalDropletFilm, MAX_DILUTE_SPHERE_VAPOR_MOLE_FRACTION,
};
use crate::graph::thermal::medium::mass_transfer::sphere::ForcedSphereMassTransfer;
use fs_material::evaporation::gas::droplet::advance_droplet_with_gas_and_heat_at_sherwood;
use fs_material::fluid_transport::FluidTransport;
use fs_thermal::evaporation::{EvaporationProperties, GasFilm, SurfaceFlux, surface_flux};
use std::cell::RefCell;

const MODEL: ForcedSphereMassTransfer = ForcedSphereMassTransfer::Frossling1938;

struct Observed<'a, L: EvaporationLaw + ?Sized> {
    source: &'a L,
    temperatures: RefCell<Vec<(f64, f64)>>,
}
impl<L: EvaporationLaw + ?Sized> EvaporationLaw for Observed<'_, L> {
    fn temperature_range_k(&self) -> [f64; 2] {
        self.source.temperature_range_k()
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        self.source.properties(t)
    }
    fn interface_kinetics(&self) -> InterfaceKinetics {
        self.source.interface_kinetics()
    }
    fn evaluate_film(
        &self,
        p: EvaporationProperties,
        t: f64,
        film: GasFilm,
        kinetics: InterfaceKinetics,
    ) -> Result<SurfaceFlux, EvaporationError> {
        self.temperatures.borrow_mut().push((t, film.temperature_k));
        self.source.evaluate_film(p, t, film, kinetics)
    }
}

#[test]
fn mean_transport_uses_original_composition_pressure_and_the_same_molecular_claim() {
    let catalog = catalog();
    let original = gas(&catalog, 400.0, 100_000.0, 0.02);
    let saved = original.clone();
    let card = diffusion();
    let context = point(250.0, 50_000.0);
    let pair = source(&card, &context).source;
    let vapor = SpeciesId::new("B").unwrap();
    let pinned = SourceSphericalDropletFilm::new(&pair, &vapor, &original, || false).unwrap();
    let selected = pinned.selected_diffusivity().clone();
    let forced = SourceForcedSphericalDropletFilm::try_new(&pinned, MODEL, 0.5, 0.04).unwrap();
    let sample = forced.resolve(geometry(), 380.0, &original, || false).unwrap();
    let mean = 390.0;
    let query = original.point().clone()
        .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), mean).unwrap();
    let expected_state = original.model().evaluate(&query).unwrap();
    assert_eq!(sample.film_state.as_ref(), Some(&expected_state));
    assert_eq!(sample.transport, expected_state.transport_state().unwrap());
    assert_ne!(sample.transport, original.state().transport_state().unwrap());
    assert_eq!(sample.surface_temperature_k, 380.0);
    assert_eq!(sample.bulk_temperature_k, 400.0);
    assert_eq!(sample.film_temperature_k, mean);
    assert_eq!(sample.total_pressure_pa, original.state().state().pressure);
    assert_eq!(sample.diffusion.material.state_coordinate("temperature"), Some(mean));
    assert_eq!(
        sample.diffusion.material.state_coordinate("pressure"),
        Some(original.state().state().pressure),
    );
    let diffusivity = 1e-5 + (mean - 200.0) / 700.0 * 2e-5;
    near(sample.diffusion.diffusivity_m2_s, diffusivity, 1e-18, 1e-13);
    let nu = sample.transport.kinematic_viscosity_m2_s();
    let re = 0.5 * geometry().diameter_m() / nu;
    let sh = 2.0 + 0.552 * re.sqrt() * (nu / diffusivity).cbrt();
    near(sample.correlation.reynolds(), re, 1e-13, 1e-13);
    near(sample.correlation.sherwood_number(), sh, 1e-13, 1e-13);
    near(
        sample.mass_transfer_coefficient_m_s(),
        sh * diffusivity / geometry().diameter_m(),
        1e-15,
        1e-13,
    );
    assert!(std::ptr::eq(forced.source(), &pinned));
    assert!(std::ptr::eq(pinned.diffusion_source().material, &card));
    assert!(std::ptr::eq(pinned.diffusion_source().context, &context));
    assert_eq!(pinned.selected_diffusivity(), &selected);
    assert_eq!(original, saved);

    let isothermal = forced.resolve(geometry(), 400.0, &original, || false).unwrap();
    assert_eq!(isothermal.film_state, None);
    assert_eq!(isothermal.transport, original.state().transport_state().unwrap());
    let changed = gas(&catalog, 450.0, 130_000.0, 0.01);
    let larger = SphericalDropletGeometry::from_mass(8e-6, 1000.0).unwrap();
    let refreshed = forced.resolve(larger, 440.0, &changed, || false).unwrap();
    assert_eq!(refreshed.film_temperature_k, 445.0);
    assert_eq!(
        refreshed.diffusion.material.state_coordinate("pressure"),
        Some(changed.state().state().pressure),
    );
    assert_ne!(refreshed.diffusion.diffusivity_m2_s, sample.diffusion.diffusivity_m2_s);
    assert_ne!(refreshed.correlation.reynolds(), sample.correlation.reynolds());
    assert_eq!(pinned.selected_diffusivity(), &selected);
}

#[test]
fn both_phase_directions_keep_original_kinetics_surface_temperature_and_conservation() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let binding = liquid_binding();
    let card = diffusion();
    let context = point(400.0, 100_000.0);
    let pair = source(&card, &context).source;
    let kinetic = KineticEvaporation::try_new(
        &LiquidPhase,
        InterfaceKinetics::Schrage { accommodation_coefficient: 1e-5 },
    ).unwrap();
    let carrier = SpeciesId::new("A").unwrap();
    for (gas_t, fraction) in [(400.0, 0.0), (430.0, 0.02)] {
        let original = gas(&catalog, gas_t, 100_000.0, fraction);
        let saved = original.clone();
        let pinned = SourceSphericalDropletFilm::new(
            &pair, &binding.species, &original, || false,
        ).unwrap();
        let forced = SourceForcedSphericalDropletFilm::try_new(
            &pinned, MODEL, 0.5, 0.04,
        ).unwrap();
        let observed = Observed { source: &kinetic, temperatures: RefCell::new(vec![]) };
        let state = sphere(&caloric, 400.0, 1e-6);
        let dt = 0.1;
        let next = forced.advance_with_gas(
            &observed, &binding, state, &original, &catalog,
            |h| inverse(&caloric, h), dt, phase_controls().phase, || false,
        ).unwrap();
        let first = &next.films()[0];
        let direct_flux = surface_flux(
            &kinetic,
            400.0,
            GasFilm {
                temperature_k: first.film_temperature_k,
                total_pressure_pa: original.state().state().pressure,
                vapor_partial_pressure_pa: original.species_partial_pressure_pa(&binding.species)
                    .unwrap(),
                mass_transfer_coefficient_m_s: first.mass_transfer_coefficient_m_s(),
            },
        ).unwrap();
        let candidate = next.original();
        near(
            candidate.step().requested_transferred_mass_kg,
            dt * first.geometry.area_m2() * direct_flux.mass_flux_kg_m2_s,
            1e-22,
            2e-13,
        );
        assert_eq!(candidate.step().transferred_mass_kg > 0.0, fraction == 0.0);
        assert_eq!(
            candidate.geometry()[1].radius_m() < candidate.geometry()[0].radius_m(),
            fraction == 0.0,
        );
        assert_eq!(
            candidate.report().binary_diffusivity_m2_s,
            first.diffusion.diffusivity_m2_s,
        );
        assert_eq!(
            candidate.report().sherwood_number,
            first.correlation.sherwood_number(),
        );
        assert_eq!(
            candidate.report().mass_transfer_coefficient_m_s,
            first.mass_transfer_coefficient_m_s(),
        );
        {
            let seen = observed.temperatures.borrow();
            assert!(!seen.is_empty());
            assert_eq!(seen[0].0, 400.0);
            assert!(seen.iter().all(|(_, film_t)| *film_t == first.film_temperature_k));
        }
        near(
            candidate.droplet().mass_kg + candidate.gas().species_mass_kg(&binding.species).unwrap(),
            state.mass_kg + original.species_mass_kg(&binding.species).unwrap(),
            phase_controls().phase.gas.mass_tolerance_kg,
            0.0,
        );
        near(
            candidate.gas().species_mass_kg(&carrier).unwrap(),
            original.species_mass_kg(&carrier).unwrap(),
            phase_controls().phase.gas.mass_tolerance_kg,
            0.0,
        );
        near(
            candidate.step().liquid_enthalpy_change_j + candidate.step().gas_internal_energy_change_j,
            0.0, phase_controls().phase.gas.energy_tolerance_j, 0.0,
        );
        assert_eq!(candidate.gas().volume_m3(), original.volume_m3());
        assert_eq!(candidate.gas().point().axes().get("local-lot"), Some(&9.0));
        assert!(std::ptr::eq(candidate.catalog(), &catalog));
        assert_eq!(original, saved);
        let again = forced.advance_with_gas(
            &kinetic, &binding, candidate.droplet(), candidate.gas(), &catalog,
            |h| inverse(&caloric, h), dt, phase_controls().phase, || false,
        ).unwrap();
        assert_eq!(&again.films()[0], &next.films()[1]);
    }
}

#[test]
fn independently_authored_heat_preserves_the_original_reciprocal_owner() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let binding = liquid_binding();
    let original = gas(&catalog, 430.0, 100_000.0, 0.003);
    let saved = original.clone();
    let card = diffusion();
    let context = point(400.0, 100_000.0);
    let pair = source(&card, &context).source;
    let pinned = SourceSphericalDropletFilm::new(
        &pair, &binding.species, &original, || false,
    ).unwrap();
    let forced = SourceForcedSphericalDropletFilm::try_new(&pinned, MODEL, 0.5, 0.04).unwrap();
    let sealed = KineticEvaporation::try_new(
        &LiquidPhase,
        InterfaceKinetics::HertzKnudsen { accommodation_coefficient: 0.0 },
    ).unwrap();
    let state = sphere(&caloric, 400.0, 1e-6);
    let next = forced.advance_with_gas_and_heat(
        &sealed, &binding, state, &original, &catalog,
        |h| inverse(&caloric, h), 80.0, 0.1,
        phase_controls().phase, heat_controls(&caloric), || false,
    ).unwrap();
    let first = &next.films()[0];
    let direct = advance_droplet_with_gas_and_heat_at_sherwood(
        &sealed, &binding, state, &original, &catalog,
        |h| inverse(&caloric, h),
        first.diffusion.diffusivity_m2_s, first.correlation.sherwood_number(),
        80.0, 0.1, phase_controls().phase, heat_controls(&caloric), || false,
    ).unwrap();
    assert_eq!(next.original(), &direct);
    assert_eq!(direct.droplet().mass_kg, state.mass_kg);
    assert!(direct.heat().unwrap().heat_into_liquid_j > 0.0);
    assert_eq!(direct.step().transferred_mass_kg, 0.0);
    near(
        direct.step().liquid_enthalpy_change_j + direct.step().gas_internal_energy_change_j,
        0.0, phase_controls().phase.gas.energy_tolerance_j, 0.0,
    );
    assert_ne!(first.film_temperature_k, next.films()[1].film_temperature_k);
    assert_ne!(first.diffusion.diffusivity_m2_s, next.films()[1].diffusion.diffusivity_m2_s);
    assert_eq!(original, saved);
}

struct HighLiquid;
impl EvaporationLaw for HighLiquid {
    fn temperature_range_k(&self) -> [f64; 2] { LiquidPhase.temperature_range_k() }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        let mut result = LiquidPhase.properties(t)?;
        result.saturation_pressure_pa *= 10.0;
        result.saturation_pressure_derivative_pa_k *= 10.0;
        Ok(result)
    }
}

#[test]
fn forced_selection_refuses_invalid_speed_strong_loading_and_changed_species_atomically() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let binding = liquid_binding();
    let original = gas(&catalog, 400.0, 100_000.0, 0.003);
    let saved = original.clone();
    let card = diffusion();
    let context = point(400.0, 100_000.0);
    let pair = source(&card, &context).source;
    let pinned = SourceSphericalDropletFilm::new(
        &pair, &binding.species, &original, || false,
    ).unwrap();
    for speed in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(SourceForcedSphericalDropletFilm::try_new(&pinned, MODEL, speed, 0.04).is_err());
    }
    for bound in [0.0, -1.0, 0.051, f64::NAN, f64::INFINITY] {
        assert!(SourceForcedSphericalDropletFilm::try_new(&pinned, MODEL, 0.5, bound).is_err());
    }
    assert_eq!(MAX_DILUTE_SPHERE_VAPOR_MOLE_FRACTION, 0.05);
    let forced = SourceForcedSphericalDropletFilm::try_new(&pinned, MODEL, 0.5, 0.04).unwrap();
    assert!(matches!(
        forced.resolve(geometry(), 400.0, &original, || true),
        Err(SphericalDropletFilmError::Film(GasFilmError::Cancelled)),
    ));
    let concentrated = gas(&catalog, 400.0, 100_000.0, 0.1);
    assert!(matches!(
        forced.resolve(geometry(), 400.0, &concentrated, || false),
        Err(SphericalDropletFilmError::Film(GasFilmError::OutsideDilute { .. })),
    ));
    assert!(matches!(
        forced.advance_with_gas(
            &HighLiquid, &binding, sphere(&caloric, 400.0, 1e-6), &original, &catalog,
            |h| inverse(&caloric, h), 0.1, phase_controls().phase, || false,
        ),
        Err(SphericalDropletFilmError::Film(GasFilmError::OutsideDilute { .. })),
    ));
    let wrong = GasLiquidBinding {
        species: SpeciesId::new("A").unwrap(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    assert!(forced.advance_with_gas(
        &LiquidPhase, &wrong, sphere(&caloric, 400.0, 1e-6), &original, &catalog,
        |_| panic!("changed volatile must refuse before inversion"),
        0.1, phase_controls().phase, || false,
    ).is_err());
    assert!(forced.advance_with_gas_and_heat(
        &LiquidPhase, &wrong, sphere(&caloric, 400.0, 1e-6), &original, &catalog,
        |_| panic!("changed volatile must refuse before heat"),
        80.0, 0.1, phase_controls().phase, heat_controls(&caloric), || false,
    ).is_err());
    let unsupported = SourceForcedSphericalDropletFilm::try_new(
        &pinned, MODEL, 1e-6, 0.04,
    ).unwrap();
    assert!(matches!(
        unsupported.resolve(geometry(), 400.0, &original, || false),
        Err(SphericalDropletFilmError::Film(GasFilmError::OutsideCorrelation { .. })),
    ));
    assert_eq!(original, saved);
}

#[test]
fn late_mean_source_and_reynolds_refusals_discard_the_complete_heat_candidate() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let binding = liquid_binding();
    let seed = gas(&catalog, 400.0, 100_000.0, 0.003);
    let original = gas(&catalog, 430.0, 100_000.0, 0.003);
    let saved = original.clone();
    let state = sphere(&caloric, 370.0, 1e-6); // Initial mean film T is 400 K.
    let sealed = KineticEvaporation::try_new(
        &LiquidPhase,
        InterfaceKinetics::HertzKnudsen { accommodation_coefficient: 0.0 },
    ).unwrap();
    let context = point(400.0, 100_000.0);
    let card = material(vec![
        scalar(
            BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS, "temperature",
            vec![(200.0, 1.6e-5), (400.0001, 1.6e-5)], [200.0, 400.0001],
        ),
        scalar(
            BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS, "temperature",
            vec![(400.0002, 2.2e-5), (900.0, 2.2e-5)], [400.0002, 900.0],
        ),
    ]);
    let pair = source(&card, &context).source;
    let pinned = SourceSphericalDropletFilm::new(
        &pair, &binding.species, &seed, || false,
    ).unwrap();
    let selected = pinned.selected_diffusivity().clone();
    let forced = SourceForcedSphericalDropletFilm::try_new(&pinned, MODEL, 0.5, 0.04).unwrap();
    let first = forced.resolve(geometry(), 370.0, &original, || false).unwrap();
    let direct = advance_droplet_with_gas_and_heat_at_sherwood(
        &sealed, &binding, state, &original, &catalog,
        |h| inverse(&caloric, h),
        first.diffusion.diffusivity_m2_s, first.correlation.sherwood_number(),
        80.0, 0.1, phase_controls().phase, heat_controls(&caloric), || false,
    ).unwrap();
    let final_mean = 0.5 * direct.step().surface_temperature_k[1]
        + 0.5 * direct.gas().state().state().temperature;
    assert!(final_mean > 400.0002);
    assert!(matches!(
        forced.advance_with_gas_and_heat(
            &sealed, &binding, state, &original, &catalog,
            |h| inverse(&caloric, h), 80.0, 0.1,
            phase_controls().phase, heat_controls(&caloric), || false,
        ),
        Err(SphericalDropletFilmError::Film(GasFilmError::Source(_))),
    ));
    assert_eq!(pinned.selected_diffusivity(), &selected);
    let fresh = SourceSphericalDropletFilm::new(
        &pair, &binding.species, direct.gas(), || false,
    ).unwrap();
    let reselection = SourceForcedSphericalDropletFilm::try_new(
        &fresh, MODEL, 0.5, 0.04,
    ).unwrap();
    assert!(reselection.resolve(
        direct.geometry()[1], direct.step().surface_temperature_k[1],
        direct.gas(), || false,
    ).is_ok());

    let supported = diffusion();
    let pair = source(&supported, &context).source;
    let pinned = SourceSphericalDropletFilm::new(
        &pair, &binding.species, &original, || false,
    ).unwrap();
    let probe = SourceForcedSphericalDropletFilm::try_new(&pinned, MODEL, 0.5, 0.04).unwrap();
    let first = probe.resolve(geometry(), 370.0, &original, || false).unwrap();
    let speed = 2.000001 * first.transport.kinematic_viscosity_m2_s()
        / first.geometry.diameter_m();
    let boundary = SourceForcedSphericalDropletFilm::try_new(
        &pinned, MODEL, speed, 0.04,
    ).unwrap();
    assert!(boundary.resolve(geometry(), 370.0, &original, || false).is_ok());
    assert!(matches!(
        boundary.advance_with_gas_and_heat(
            &sealed, &binding, state, &original, &catalog,
            |h| inverse(&caloric, h), 80.0, 0.1,
            phase_controls().phase, heat_controls(&caloric), || false,
        ),
        Err(SphericalDropletFilmError::Film(GasFilmError::OutsideCorrelation { .. })),
    ));

    // The phase hook also admits the actual post-heat state, not only the
    // interval's original and accepted endpoints.
    let initial_fraction = LiquidPhase.properties(370.0).unwrap().saturation_pressure_pa
        / original.state().state().pressure;
    let narrow = SourceForcedSphericalDropletFilm::try_new(
        &pinned, MODEL, 0.5, initial_fraction * (1.0 + 1e-5),
    ).unwrap();
    assert!(matches!(
        narrow.advance_with_gas_and_heat(
            &sealed, &binding, state, &original, &catalog,
            |h| inverse(&caloric, h), 80.0, 0.1,
            phase_controls().phase, heat_controls(&caloric), || false,
        ),
        Err(SphericalDropletFilmError::Phase(
            fs_material::evaporation::gas::GasEvaporationError::Surface(
                EvaporationError::Invalid(_),
            ),
        )),
    ));
    assert_eq!(original, saved);
}

#[test]
fn forced_original_if97_and_nasa_gas_keep_original_species_energy_and_sources() {
    use fs_material::evaporation::WaterEvaporation;
    use fs_material::phase::water::{IsobaricWater, WaterChartControls};
    use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;

    let nitrogen = nasa_species(
        "N2",
        0.028_013_4,
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/matdb/seed-v1/nitrogen-nasa9-mcbride2002/nasa9.tsv"
        )),
    );
    let vapor = nasa_species(
        "H2O",
        0.018_015_28,
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/matdb/seed-v1/water-vapor-nasa9-mcbride2002/nasa9.tsv"
        )),
    );
    let source = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![nitrogen, vapor],
            Composition::new(CompositionBasis::MoleFraction, vec![0.99, 0.01]).unwrap(),
        )
        .unwrap(),
        vec![
            transport("N2", 1.8e-5, 0.026),
            transport("H2O", 1.0e-5, 0.02),
        ],
        Some(&SpeciesId::new("H2O").unwrap()),
    )
    .unwrap();
    let chart = IsobaricWater::new(
        Pressure::new(100_000.0),
        WaterChartControls {
            minimum_temperature_k: 273.16,
            maximum_temperature_k: 450.0,
            enthalpy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
            maximum_entropy_seam_j_kg_k: 0.1,
            heat_capacity: If97HeatCapacityBoundSolve {
                maximum_intervals: 65535,
            },
        },
        || false,
    )
    .unwrap();
    let water = WaterEvaporation::try_new(chart, [280.0, 340.0]).unwrap();
    let vapor = SpeciesId::new("H2O").unwrap();
    let carrier = SpeciesId::new("N2").unwrap();
    let binding = GasLiquidBinding {
        species: vapor.clone(),
        liquid_enthalpy_offset_j_kg: source
            .species_specific_enthalpy_j_kg(&vapor, 300.0)
            .unwrap()
            - water.vapor_specific_enthalpy_j_kg(300.0).unwrap(),
    };
    let state = SphericalDropletState {
        mass_kg: 1e-6,
        specific_enthalpy_j_kg: water.liquid_state(300.0).unwrap().specific_enthalpy_j_kg(),
    };
    // Explicit manufactured D(T) and accommodation exercise source plumbing;
    // these are not experimental water-evaporation validation data.
    let card = material(vec![scalar(
        BINARY_DIFFUSIVITY_PROPERTY,
        BINARY_DIFFUSIVITY_DIMS,
        "temperature",
        vec![(280.0, 1.8e-5), (340.0, 2.5e-5)],
        [280.0, 340.0],
    )]);
    let context = point(300.0, 100_000.0);
    let pair = BinaryDiffusionMaterial {
        species: [carrier.clone(), vapor.clone()],
        material: &card,
        context: &context,
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    let selected = KineticEvaporation::try_new(
        &water,
        InterfaceKinetics::HertzKnudsen {
            accommodation_coefficient: 1e-5,
        },
    )
    .unwrap();
    let mut solve = phase_controls().phase;
    solve.gas.gas.minimum_temperature_k = 280.0;
    solve.gas.gas.maximum_temperature_k = 340.0;
    solve.gas.gas.energy_tolerance_j_kg = 1e-6;
    solve.gas.energy_tolerance_j = 1e-7;
    solve.maximum_temperature_increment_k = 1.0;
    solve.latent_compatibility_tolerance_j_kg = 5.0;
    for (t, fraction) in [(300.0, 0.01), (310.0, 0.045)] {
        // Composition coordinates follow the original catalog's canonical
        // species order (H2O sorts before N2), never a guessed vapor index.
        let fractions: Vec<_> = source.thermodynamics().components().iter()
            .map(|component| if component.species() == &vapor { fraction } else { 1.0 - fraction })
            .collect();
        let model = source.with_composition(
            Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
        ).unwrap();
        let original = ClosedGasParcel::new(
            model, point(t, 100_000.0).with("local-lot", 9.0).unwrap(), 0.01,
        ).unwrap();
        let saved = original.clone();
        let prepared = SourceSphericalDropletFilm::new(&pair, &vapor, &original, || false).unwrap();
        let forced = SourceForcedSphericalDropletFilm::try_new(
            &prepared, MODEL, 0.5, 0.05,
        ).unwrap();
        let next = forced
            .advance_with_gas(
                &selected,
                &binding,
                state,
                &original,
                &source,
                |h| water.temperature_at_specific_enthalpy_k(h),
                0.01,
                solve,
                || false,
            )
            .unwrap();
        let candidate = next.original();
        assert!(candidate.report().sherwood_number > 2.0);
        assert_eq!(
            candidate.report().binary_diffusivity_m2_s,
            next.films()[0].diffusion.diffusivity_m2_s,
        );
        near(
            next.films()[0].film_temperature_k,
            0.5 * 300.0 + 0.5 * original.state().state().temperature,
            1e-8,
            1e-12,
        );
        assert_eq!(candidate.step().transferred_mass_kg > 0.0, fraction < 0.02);
        assert_eq!(
            candidate.geometry()[1].radius_m() < candidate.geometry()[0].radius_m(),
            fraction < 0.02
        );
        near(
            candidate.gas().species_mass_kg(&vapor).unwrap() + candidate.droplet().mass_kg,
            original.species_mass_kg(&vapor).unwrap() + state.mass_kg,
            solve.gas.mass_tolerance_kg,
            0.0,
        );
        near(
            candidate.gas().species_mass_kg(&carrier).unwrap(),
            original.species_mass_kg(&carrier).unwrap(),
            solve.gas.mass_tolerance_kg,
            0.0,
        );
        near(
            candidate.step().liquid_enthalpy_change_j
                + candidate.step().gas_internal_energy_change_j,
            0.0,
            solve.gas.energy_tolerance_j,
            0.0,
        );
        assert!(candidate.step().energy_residual_j.abs() <= solve.gas.energy_tolerance_j);
        assert_eq!(candidate.gas().volume_m3(), original.volume_m3());
        assert!(std::ptr::eq(candidate.catalog(), &source));
        assert!(std::ptr::eq(selected.source(), &water));
        assert_eq!(
            next.films()[1]
                .diffusion
                .material
                .state_coordinate("temperature"),
            Some(next.films()[1].film_temperature_k)
        );
        assert_eq!(original, saved);
    }
}
