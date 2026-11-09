use super::*;
use fs_material::state_point::{DENSITY_PROPERTY, ResolvedMaterialStatePoint};
use fs_material::thermo_optic::{SchottThermoOpticLaw, ThermoOpticCoefficients};

const SOURCE_PROPERTIES: usize = 29;

fn reference_fixture(
    mechanical_maximum_k: f64,
) -> (
    Arc<MaterialCard>,
    SpectralOpticalStatePoint,
    ResolvedMaterialStatePoint,
) {
    let (card, original) = fixture(|claims| {
        for claim in claims {
            let maximum = if claim.key.name() == DENSITY_PROPERTY {
                mechanical_maximum_k
            } else {
                300.0
            };
            // Optics are measured/admitted only at T0. Their underlying curves
            // deliberately contain different hot values that must not be used.
            claim.validity = claim.validity.clone().with_quantity(
                "thermal",
                quantity(QuantityKind::AbsoluteTemperature),
                300.0,
                maximum,
            );
        }
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
    (card, original, mechanical)
}

pub(super) fn thermal_law(maximum_k: f64) -> SchottThermoOpticLaw {
    // Manufactured coefficients make index feedback observable; these are
    // equation-verification data, not a claim about an actual infrared glass.
    SchottThermoOpticLaw::try_new(
        ThermoOpticCoefficients {
            d: [2e-3, 1e-6, 0.0],
            e: [4e-4, 0.0],
            resonance_um: 2.0,
        },
        300.0,
        [300.0, maximum_k],
        [8000.0, 14_000.0],
        ContentHash([214; 32]),
    )
    .unwrap()
}

pub(super) fn schott_law(
    mechanical_maximum_k: f64,
    optical_maximum_k: f64,
) -> Arc<MaterialSlabLaw> {
    let (card, original, mechanical) = reference_fixture(mechanical_maximum_k);
    Arc::new(
        MaterialSlabLaw::try_schott_with_reference_absorption(
            card,
            original,
            &mechanical,
            thermal_law(optical_maximum_k),
            "thermal",
            Some(751.0),
            SOURCE_PROPERTIES,
        )
        .unwrap(),
    )
}

fn schott_budget() -> NativeCoatingExchangeBudget {
    NativeCoatingExchangeBudget {
        max_total_optical_knots: SOURCE_PROPERTIES,
        ..budget()
    }
}

pub(super) fn phase(wavelength_nm: f64, temperature_k: f64) -> f64 {
    let n0 = index(wavelength_nm, 300.0);
    let dt = temperature_k - 300.0;
    let l = wavelength_nm * 0.001;
    n0 + (n0 * n0 - 1.0) / (2.0 * n0) * (2e-3 * dt + 1e-6 * dt * dt + 4e-4 * dt / (l * l - 4.0))
}

fn thermal_absorptance(wavelength: f64, cosine: f64, temperature: f64) -> [f64; 2] {
    // Independent real-interface Fresnel powers and the complete incoherent
    // geometric series, with refraction changing the optical path through alpha.
    let n = phase(wavelength, temperature);
    let ct = (1.0 - (1.0 - cosine * cosine) / (n * n)).sqrt();
    let rs = ((cosine - n * ct) / (cosine + n * ct)).powi(2);
    let rp = ((n * cosine - ct) / (n * cosine + ct)).powi(2);
    let attenuation = (-absorption(wavelength, 300.0) * DEPTH / ct).exp();
    [rs, rp].map(|r| (1.0 - r) * (1.0 - attenuation) / (1.0 - r * attenuation))
}

fn thermal_outward(body_t: f64, optical_t: f64, bath_t: f64) -> f64 {
    let width = 6000.0 / Q.wavelength_bins as f64;
    let mut watts = 0.0;
    for i in 0..Q.wavelength_bins {
        let wavelength = 8000.0 + (i as f64 + 0.5) * width;
        for j in 0..Q.cosine_bins {
            let cosine = (j as f64 + 0.5) / Q.cosine_bins as f64;
            let [s, p] = thermal_absorptance(wavelength, cosine, optical_t);
            watts += AREA * core::f64::consts::PI * width / Q.cosine_bins as f64
                * cosine
                * (s + p)
                * (planck(wavelength, body_t) - planck(wavelength, bath_t));
        }
    }
    watts
}

fn expected_endpoint(bath: f64, feedback: bool) -> (f64, usize) {
    let (mut time, mut temperature, mut steps) = (0.0, 400.0, 0);
    while time < 0.2 {
        let optics = if feedback { temperature } else { 400.0 };
        let power = thermal_outward(temperature, optics, bath);
        let conductance = power / (temperature - bath);
        // The actual caloric owner has 1 J/K storage, and combines its fixed
        // maximum step with the positive radiative monotonicity limit.
        let maximum = config().max_dt_s;
        let dt = (maximum / (1.0 + maximum * conductance / 0.8)).min(0.2 - time);
        time = if dt == 0.2 - time { 0.2 } else { time + dt };
        temperature -= dt * power;
        steps += 1;
        assert!(steps <= 16);
    }
    (temperature, steps)
}

#[test]
fn reference_only_optics_follow_hot_mechanics_without_requerying_phase_or_absorption() {
    let (card, original, mechanical) = reference_fixture(600.0);
    assert!(original.at_query_point(&card, &point(450.0)).is_err());
    let law = MaterialSlabLaw::try_schott_with_reference_absorption(
        Arc::clone(&card),
        original.clone(),
        &mechanical,
        thermal_law(600.0),
        "thermal",
        Some(751.0),
        SOURCE_PROPERTIES,
    )
    .unwrap();
    assert_eq!(law.material_source(), &original);
    assert_eq!(law.budget_size(), SOURCE_PROPERTIES);
    assert_eq!(law.temperature_bounds_k(), (300.0, 600.0));
    let cold = law.at_temperature(300.0).unwrap();
    for t in [450.0, 600.0, 350.0, 300.0] {
        let accepted = mechanical.at_query_point(&card, &point(t)).unwrap();
        assert_eq!(accepted.properties().len(), 1);
        let binding = law.at_temperature(t).unwrap();
        assert_eq!(binding.accepted_state_identity(), accepted.identity());
        let glass = glass(binding);
        for wavelength in [8001.0, 10_625.0, 13_999.0] {
            near(glass.ior().eval(wavelength).unwrap(), phase(wavelength, t));
            near(
                glass.absorption_ref().extinction_per_m(wavelength).unwrap(),
                absorption(wavelength, 300.0),
            );
        }
        assert!(glass.ior().eval(550.0).is_err());
    }
    assert_eq!(law.at_temperature(300.0).unwrap(), cold);
    assert_ne!(
        law.at_temperature(450.0).unwrap().identity(),
        cold.identity()
    );
    let current_claims = MaterialSlabLaw::try_new(
        card,
        original,
        "thermal",
        (300.0, 600.0),
        Some(751.0),
        SOURCE_PROPERTIES,
    )
    .unwrap();
    assert!(current_claims.at_temperature(450.0).is_err());
    assert_ne!(current_claims.identity(), law.identity());
}

#[test]
fn native_radiative_heating_and_cooling_recompute_schott_refraction_on_accepted_clock() {
    let law = schott_law(600.0, 600.0);
    for bath in [600.0, 300.0] {
        let thermal = native(400.0);
        let patch = source(Arc::clone(&law), bath);
        let initial =
            NativeCoatingExchange::try_new(&thermal, &[patch.clone()], schott_budget(), || Ok(()))
                .unwrap();
        let old_optics = initial.emitter(10).unwrap().optical_identity();
        let sample = initial
            .emitter(10)
            .unwrap()
            .evaluate(10_625.0, 0.6)
            .unwrap();
        for (actual, absorptance) in sample
            .polarized_radiance_w_m2_sr_nm()
            .into_iter()
            .zip(thermal_absorptance(10_625.0, 0.6, 400.0))
        {
            near(actual, absorptance * planck(10_625.0, 400.0) * 0.5);
        }
        near(
            initial.sources_w()[&CELL],
            -thermal_outward(400.0, 400.0, bath),
        );
        let (expected, steps) = expected_endpoint(bath, true);
        let (advanced, observed, audit) = thermal
            .advance_observed(&initial, 0.2, config(), || Ok(()))
            .unwrap();
        let actual = advanced.temperature_k(CELL).unwrap();
        near(actual, expected);
        near(advanced.energy_j() - thermal.energy_j(), expected - 400.0);
        near(observed.energy().thermal_external_j, expected - 400.0);
        near(audit.external_heat_j, -observed.energy().outward_j);
        near(observed.energy().balance_residual_j, 0.0);
        near(
            observed.sources_w()[&CELL],
            -thermal_outward(actual, actual, bath),
        );
        assert_eq!(audit.substeps, steps);
        assert_ne!(observed.emitter(10).unwrap().optical_identity(), old_optics);
        assert_eq!(
            observed.emitter(10).unwrap().film_temperature_k(0),
            Some(actual)
        );

        let mut fixed = patch;
        fixed.stack = IsothermalCoatingStack::IncoherentSlab {
            glass: glass(law.at_temperature(400.0).unwrap()),
            thickness_m: DEPTH,
            opposite_ior: vacuum(),
        };
        let frozen =
            NativeCoatingExchange::try_new(&thermal, &[fixed], schott_budget(), || Ok(())).unwrap();
        let (frozen_state, _, _) = thermal
            .advance_observed(&frozen, 0.2, config(), || Ok(()))
            .unwrap();
        let frozen_t = frozen_state.temperature_k(CELL).unwrap();
        near(frozen_t, expected_endpoint(bath, false).0);
        assert!((actual - frozen_t).abs() > 1e-7);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(initial.emitter(10).unwrap().optical_identity(), old_optics);
        assert_eq!(initial.energy().outward_j, 0.0);
    }
}

#[test]
fn source_domains_budgets_and_cancellation_preserve_accepted_native_state() {
    for (mechanical_maximum, optical_maximum) in [(400.0, 600.0), (600.0, 400.0)] {
        let law = schott_law(mechanical_maximum, optical_maximum);
        assert!(law.at_temperature(400.1).is_err());
        let thermal = native(400.0);
        let patch = source(law, 600.0);
        let initial =
            NativeCoatingExchange::try_new(&thermal, &[patch], schott_budget(), || Ok(())).unwrap();
        let identity = initial.emitter(10).unwrap().identity();
        let optics = initial.emitter(10).unwrap().optical_identity();
        let watts = initial.sources_w().clone();
        let work = initial.work();
        assert!(
            thermal
                .advance_observed(&initial, 0.1, config(), || Ok(()))
                .is_err()
        );
        assert!(matches!(
            thermal.advance_observed(&initial, 0.1, config(), || Err(DepositionError::Cancelled)),
            Err(NativeCoatingExchangeError::Thermal(
                DepositionError::Cancelled
            ))
        ));
        assert_eq!(thermal.temperature_k(CELL).unwrap(), 400.0);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.emitter(10).unwrap().identity(), identity);
        assert_eq!(initial.emitter(10).unwrap().optical_identity(), optics);
        assert_eq!(initial.sources_w(), &watts);
        assert_eq!(initial.work(), work);
        assert_eq!(initial.energy().outward_j, 0.0);
    }
    let (card, original, mechanical) = reference_fixture(600.0);
    assert!(matches!(
        MaterialSlabLaw::try_schott_with_reference_absorption(
            card,
            original,
            &mechanical,
            thermal_law(600.0),
            "thermal",
            Some(751.0),
            SOURCE_PROPERTIES - 1,
        ),
        Err(SurfaceHeatingError::Budget)
    ));
    let law = schott_law(600.0, 600.0);
    let thermal = native(400.0);
    let front = source(Arc::clone(&law), 400.0);
    let mut rear = source(law, 400.0);
    rear.source_id = 20;
    // Two physical patches share one original optical/mechanical requirement
    // budget; the separate accepted mechanical bundle is charged exactly once.
    let patches = [front, rear];
    let observed =
        NativeCoatingExchange::try_new(&thermal, &patches, schott_budget(), || Ok(())).unwrap();
    assert_eq!(observed.sources_w()[&CELL], 0.0);
    let mut insufficient = schott_budget();
    insufficient.max_total_optical_knots -= 1;
    assert!(matches!(
        NativeCoatingExchange::try_new(&thermal, &patches, insufficient, || Ok(())),
        Err(NativeCoatingExchangeError::Budget)
    ));
}
