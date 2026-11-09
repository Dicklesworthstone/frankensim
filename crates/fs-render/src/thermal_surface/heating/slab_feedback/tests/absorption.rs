use super::*;
use crate::material_state::thermo_optic::binding::{
    MaterialThermoOpticBinding, MaterialThermoOpticBindingError,
};
use crate::thermal_surface::heating::{
    CoatingSurfaceSource, StatefulMaterialSlabIlluminationSource, SurfaceHeatingBudget,
    SurfaceHeatingTolerance, SurfaceIncidence, SurfaceOpticalHeating,
};
use crate::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, ResolvedMaterialStatePoint, SpectralDielectricAbsorptionStatePoint,
    resolve_spectral_dielectric_absorption_state_point,
};

const PROPERTIES: usize = 47;
const WAVELENGTH: f64 = 11_000.0;
const POWER: f64 = 30.0;

fn sources(
    edit: impl FnOnce(&mut Vec<PropertyClaim>),
) -> (
    Arc<MaterialCard>,
    SpectralOpticalStatePoint,
    ResolvedMaterialStatePoint,
    SpectralDielectricAbsorptionStatePoint,
) {
    let (card, optical) = fixture(|claims| {
        for claim in claims.iter_mut() {
            if SPECTRAL_DIELECTRIC_IOR_PROPERTIES.contains(&claim.key.name()) {
                claim.validity = claim.validity.clone().with_quantity(
                    "thermal",
                    quantity(QuantityKind::AbsoluteTemperature),
                    300.0,
                    300.0,
                );
            }
        }
        edit(claims);
    });
    let mechanical = resolve_material_state_point(
        &card,
        &point(300.0),
        &[ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let loss = resolve_spectral_dielectric_absorption_state_point(
        &card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (card, optical, mechanical, loss)
}

fn pair() -> (Arc<MaterialSlabLaw>, Arc<MaterialSlabLaw>) {
    let (card, optical, mechanical, loss) = sources(|_| {});
    let current = MaterialSlabLaw::try_schott_with_absorption_source(
        Arc::clone(&card),
        optical.clone(),
        &mechanical,
        loss,
        schott::thermal_law(600.0),
        "thermal",
        Some(751.0),
        PROPERTIES,
    )
    .unwrap();
    let frozen = MaterialSlabLaw::try_schott_with_reference_absorption(
        card,
        optical,
        &mechanical,
        schott::thermal_law(600.0),
        "thermal",
        Some(751.0),
        29,
    )
    .unwrap();
    (Arc::new(current), Arc::new(frozen))
}

fn light(
    thermal: &DepositionSimulation<'_>,
    optics: Arc<MaterialSlabLaw>,
) -> SurfaceOpticalHeating<1> {
    let source = StatefulMaterialSlabIlluminationSource {
        source_id: 1,
        receiver: receiver(),
        optics,
        thickness_m: DEPTH,
        incidence: SurfaceIncidence {
            wavelengths_nm: [WAVELENGTH],
            incident_eta: [1.0],
            incidence_cosine: 1.0,
            s_power_fraction: [0.5],
        },
        transmitted_eta: [1.0],
        source_identity: ContentHash([216; 32]),
        power: SpectralPowerSchedule::try_new(
            &[SpectralPowerSegment {
                start_s: 0.0,
                end_s: 0.2,
                start_power_w: [POWER],
                end_power_w: [POWER],
            }],
            1,
        )
        .unwrap(),
    };
    SurfaceOpticalHeating::try_new_coating_stacks(
        &[CoatingSurfaceSource::MaterialIncoherentSlab(&source)],
        thermal,
        SurfaceHeatingBudget {
            max_sources: 1,
            max_total_power_segments: 1,
            max_source_evaluations: 32,
        },
        PROPERTIES,
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
}

fn normal_power(temperature: f64) -> [f64; 3] {
    // Independent sum of the two-interface incoherent geometric series.
    let n = schott::phase(WAVELENGTH, temperature);
    let r = ((n - 1.0) / (n + 1.0)).powi(2);
    let a = (-absorption(WAVELENGTH, temperature) * DEPTH).exp();
    let denominator = 1.0 - r * r * a * a;
    [
        r + (1.0 - r).powi(2) * r * a * a / denominator,
        (1.0 - r).powi(2) * a / denominator,
        (1.0 - r) * (1.0 - a) / (1.0 - r * a),
    ]
}

fn directional_absorptance(
    wavelength: f64,
    cosine: f64,
    temperature: f64,
    current_loss: bool,
) -> [f64; 2] {
    let n = schott::phase(wavelength, temperature);
    let ct = (1.0 - (1.0 - cosine * cosine) / (n * n)).sqrt();
    let rs = ((cosine - n * ct) / (cosine + n * ct)).powi(2);
    let rp = ((n * cosine - ct) / (n * cosine + ct)).powi(2);
    let loss_temperature = if current_loss { temperature } else { 300.0 };
    let attenuation = (-absorption(wavelength, loss_temperature) * DEPTH / ct).exp();
    [rs, rp].map(|r| (1.0 - r) * (1.0 - attenuation) / (1.0 - r * attenuation))
}

fn cooling_power(temperature: f64, current_loss: bool) -> f64 {
    let width = 6000.0 / Q.wavelength_bins as f64;
    let mut watts = 0.0;
    for i in 0..Q.wavelength_bins {
        let wavelength = 8000.0 + (i as f64 + 0.5) * width;
        for j in 0..Q.cosine_bins {
            let cosine = (j as f64 + 0.5) / Q.cosine_bins as f64;
            let [s, p] = directional_absorptance(wavelength, cosine, temperature, current_loss);
            watts += AREA * core::f64::consts::PI * width / Q.cosine_bins as f64
                * cosine
                * (s + p)
                * (planck(wavelength, temperature) - planck(wavelength, 300.0));
        }
    }
    watts
}

#[test]
fn current_absorption_and_reference_phase_drive_native_beer_and_illumination() {
    let (card, original, mechanical, loss) = sources(|_| {});
    assert!(original.at_query_point(&card, &point(400.0)).is_err());
    assert_eq!(loss.resolved().properties().len(), 18);
    let (current, frozen) = pair();
    assert_eq!(current.budget_size(), PROPERTIES);
    assert_eq!(frozen.budget_size(), 29);
    assert_ne!(current.identity(), frozen.identity());
    let hot = current.at_temperature(400.0).unwrap();
    let accepted = mechanical.at_query_point(&card, &point(400.0)).unwrap();
    assert_eq!(hot.accepted_state_identity(), accepted.identity());
    let hot_glass = glass(hot);
    let frozen_glass = glass(frozen.at_temperature(400.0).unwrap());
    assert_eq!(hot_glass.ior(), frozen_glass.ior());
    near(
        hot_glass.ior().eval(WAVELENGTH).unwrap(),
        schott::phase(WAVELENGTH, 400.0),
    );
    near(
        hot_glass.absorption_ref().extinction_per_m(WAVELENGTH).unwrap(),
        absorption(WAVELENGTH, 400.0),
    );
    near(
        frozen_glass.absorption_ref().extinction_per_m(WAVELENGTH).unwrap(),
        absorption(WAVELENGTH, 300.0),
    );
    near(
        hot_glass.absorption_ref().transmittance(WAVELENGTH, DEPTH).unwrap(),
        (-absorption(WAVELENGTH, 400.0) * DEPTH).exp(),
    );

    let thermal = native(400.0);
    let input = light(&thermal, current);
    let mut expected_temperature = 400.0;
    let mut expected = [0.0; 3];
    for _ in 0..2 {
        let fractions = normal_power(expected_temperature);
        near(fractions.iter().sum(), 1.0);
        for channel in 0..3 {
            expected[channel] += 0.1 * POWER * fractions[channel];
        }
        expected_temperature += 0.1 * POWER * fractions[2];
    }
    let (advanced, observed, audit) = thermal
        .advance_observed(&input, 0.2, config(), || Ok(()))
        .unwrap();
    near(advanced.temperature_k(CELL).unwrap(), expected_temperature);
    near(observed.energy().incident_j[0], 0.2 * POWER);
    near(observed.energy().reflected_j[0], expected[0]);
    near(observed.energy().transmitted_j[0], expected[1]);
    near(observed.energy().absorbed_j[0], expected[2]);
    near(observed.energy().thermal_input_j, expected[2]);
    near(audit.external_heat_j, expected[2]);
    near(advanced.energy_j() - thermal.energy_j(), expected[2]);
    let fixed_loss = light(&thermal, frozen);
    let (fixed_state, _, _) = thermal
        .advance_observed(&fixed_loss, 0.2, config(), || Ok(()))
        .unwrap();
    assert!(expected_temperature - fixed_state.temperature_k(CELL).unwrap() > 0.1);
    assert_eq!(thermal.temperature_k(CELL).unwrap(), 400.0);
    assert_eq!(input.energy().thermal_input_j, 0.0);
}

#[test]
fn source_absorption_updates_directional_emission_and_native_radiative_cooling() {
    let (current, frozen) = pair();
    let mut endpoints = [0.0; 2];
    for (lane, optics) in [current, frozen].into_iter().enumerate() {
        let current_loss = lane == 0;
        let thermal = native(400.0);
        let initial = NativeCoatingExchange::try_new(
            &thermal,
            &[source(optics, 300.0)],
            NativeCoatingExchangeBudget {
                max_total_optical_knots: PROPERTIES,
                ..budget()
            },
            || Ok(()),
        )
        .unwrap();
        let sample = initial
            .emitter(10)
            .unwrap()
            .evaluate(10_625.0, 0.6)
            .unwrap()
            .polarized_radiance_w_m2_sr_nm();
        let absorptance = directional_absorptance(10_625.0, 0.6, 400.0, current_loss);
        for pol in 0..2 {
            near(sample[pol], absorptance[pol] * planck(10_625.0, 400.0) * 0.5);
        }
        near(initial.sources_w()[&CELL], -cooling_power(400.0, current_loss));

        let (mut time, mut expected, mut steps) = (0.0, 400.0, 0);
        while time < 0.2 {
            let power = cooling_power(expected, current_loss);
            let conductance = power / (expected - 300.0);
            // Native caloric storage is 1 J/K. Retain the exact native
            // monotonicity clock while independently evaluating the optics.
            let maximum = config().max_dt_s;
            let dt = (maximum / (1.0 + maximum * conductance / 0.8)).min(0.2 - time);
            time = if dt == 0.2 - time { 0.2 } else { time + dt };
            expected -= dt * power;
            steps += 1;
            assert!(steps <= 16);
        }
        let (advanced, observed, audit) = thermal
            .advance_observed(&initial, 0.2, config(), || Ok(()))
            .unwrap();
        endpoints[lane] = advanced.temperature_k(CELL).unwrap();
        near(endpoints[lane], expected);
        near(observed.energy().thermal_external_j, expected - 400.0);
        near(audit.external_heat_j, -observed.energy().outward_j);
        near(advanced.energy_j() - thermal.energy_j(), expected - 400.0);
        near(observed.sources_w()[&CELL], -cooling_power(expected, current_loss));
        assert_eq!(audit.substeps, steps);
        assert_eq!(observed.emitter(10).unwrap().film_temperature_k(0), Some(endpoints[lane]));
        assert_eq!(thermal.temperature_k(CELL).unwrap(), 400.0);
        assert_eq!(initial.energy().outward_j, 0.0);
    }
    assert!(endpoints[1] - endpoints[0] > 0.1);
}

#[test]
fn stale_loss_state_transparent_support_and_source_limits_refuse_transactionally() {
    let (card, optical, mechanical, loss) = sources(|_| {});
    let SpectralOpticalStatePoint::TabulatedDielectric(table) = &optical else {
        panic!("tabulated reference phase")
    };
    let reference = MaterialOpticalBinding::try_spectral_tabulated_dielectric(
        &mechanical,
        table,
        751.0,
        None,
        ID,
    )
    .unwrap();
    let binding = MaterialThermoOpticBinding::try_new(
        reference,
        &mechanical,
        schott::thermal_law(600.0),
        "thermal",
    )
    .unwrap();
    let accepted = mechanical.at_query_point(&card, &point(400.0)).unwrap();
    assert!(matches!(
        binding.at_material_state_with_absorption(&accepted, &loss),
        Err(MaterialThermoOpticBindingError::AbsorptionStateMismatch)
    ));
    let hot_loss = loss.at_query_point(&card, &point(400.0)).unwrap();
    let hot = binding
        .at_material_state_with_absorption(&accepted, &hot_loss)
        .unwrap();
    assert_eq!(hot.accepted_state_identity(), accepted.identity());
    near(
        glass(hot).absorption_ref().extinction_per_m(WAVELENGTH).unwrap(),
        absorption(WAVELENGTH, 400.0),
    );
    assert_eq!(binding.reference_binding(), &reference);
    assert!(matches!(
        MaterialSlabLaw::try_schott_with_absorption_source(
            card,
            optical,
            &mechanical,
            loss,
            schott::thermal_law(600.0),
            "thermal",
            Some(751.0),
            PROPERTIES - 1,
        ),
        Err(SurfaceHeatingError::Budget)
    ));

    // A zero-alpha source is still finite-band. Its hot support becomes
    // disjoint from the original phase band and cannot turn into unbounded clear.
    let (card, optical, mechanical, loss) = sources(|claims| {
        for claim in claims {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                if SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES.contains(&claim.key.name()) {
                    for (_, value) in knots {
                        *value = 0.0;
                    }
                } else if SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES.contains(&claim.key.name()) {
                    knots[1].1 += 12e-6;
                    knots[2].1 += 12e-6;
                }
            }
        }
    });
    let moving = MaterialSlabLaw::try_schott_with_absorption_source(
        card,
        optical,
        &mechanical,
        loss,
        schott::thermal_law(600.0),
        "thermal",
        Some(751.0),
        PROPERTIES,
    )
    .unwrap();
    assert!(moving.at_temperature(300.0).is_ok());
    assert!(moving.at_temperature(450.0).is_err());

    // Mechanics and phase remain admitted beyond 400 K; only the independent
    // original absorption source refuses the attempted accepted endpoint.
    let (card, optical, mechanical, loss) = sources(|claims| {
        for claim in claims {
            if SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES.contains(&claim.key.name()) {
                claim.validity = claim.validity.clone().with_quantity(
                    "thermal",
                    quantity(QuantityKind::AbsoluteTemperature),
                    300.0,
                    400.0,
                );
            }
        }
    });
    let limited = Arc::new(
        MaterialSlabLaw::try_schott_with_absorption_source(
            card,
            optical,
            &mechanical,
            loss,
            schott::thermal_law(600.0),
            "thermal",
            Some(751.0),
            PROPERTIES,
        )
        .unwrap(),
    );
    let thermal = native(400.0);
    let input = light(&thermal, limited);
    assert!(thermal.advance_observed(&input, 0.1, config(), || Ok(())).is_err());
    assert!(matches!(
        thermal.advance_observed(&input, 0.1, config(), || Err(DepositionError::Cancelled)),
        Err(SurfaceHeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(thermal.temperature_k(CELL).unwrap(), 400.0);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(input.time_s(), 0.0);
    assert_eq!(input.energy().thermal_input_j, 0.0);
    assert_eq!(input.response_evaluations(), 0);
}
