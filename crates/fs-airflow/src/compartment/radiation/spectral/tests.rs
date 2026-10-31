use super::super::soot::SootOpticalDomain;
use super::super::{AbsorptionCoefficient, OpticalPath, OpticalSegment};
use super::*;
use fs_conduction::radiation::{EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, ViewFactorMatrix};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Area, Density};

fn emissivity(value: f64) -> SurfaceEmissivity {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "analytic gray-boundary fixture".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture".into(),
            phase: "solid".into(),
            process: "gray surface".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    SurfaceEmissivity::from_card("fixture", &card, 300.0, SelectionPolicy::SingleClaimOnly).unwrap()
}
fn plates() -> ParticipatingEnclosure {
    ParticipatingEnclosure::new(
        ViewFactorMatrix::infinite_parallel_plates(1.3).unwrap(),
        vec![OpticalPath {
            surfaces: [0, 1],
            exchange_area: Area::new(1.3),
            segments: vec![
                OpticalSegment {
                    layer: Layer::Lower,
                    length: Length::new(0.4),
                },
                OpticalSegment {
                    layer: Layer::Upper,
                    length: Length::new(0.8),
                },
            ],
        }],
        "authored reciprocal normal-chord fixture",
        1e-12,
        || false,
    )
    .unwrap()
}
fn band(lo: f64, hi: f64, k: [f64; 2]) -> SpectralBand {
    SpectralBand {
        wavelengths: [Length::new(lo), Length::new(hi)],
        absorption: k.map(|v| {
            LayerAbsorption::new(
                AbsorptionCoefficient::new(v),
                "authored frozen band fixture",
            )
            .unwrap()
        }),
    }
}
fn transparent() -> SpectralOptics<'static> {
    SpectralOptics::new(
        vec![band(1e-7, 1e-3, [0.0; 2])],
        [None; 2],
        "analytic boundary values are wavelength independent",
    )
    .unwrap()
}
fn band_emissivity(value: f64) -> BandSurfaceEmissivity {
    BandSurfaceEmissivity::new(
        value,
        [Temperature::new(200.0), Temperature::new(2000.0)],
        "authored constant spectral-band fixture; no measured-material claim",
    )
    .unwrap()
}
fn soot() -> RayleighSoot {
    RayleighSoot::new(
        1.75,
        1.03,
        Density::new(1800.0),
        SootOpticalDomain {
            spectral_temperature: [Temperature::new(200.0), Temperature::new(2500.0)],
            maximum_volume_fraction: 0.001,
        },
        "authored optical-index fixture; no physical-material validation",
    )
    .unwrap()
}
fn controls() -> SpectralSolve {
    SpectralSolve {
        heat_tolerance_w: 1e-7,
        radiosity_tolerance_w_m2: 1e-7,
        ..SpectralSolve::default()
    }
}
fn close(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "{actual:.16e} != {expected:.16e}; delta={:.3e}",
        actual - expected
    );
}
fn simpson(mut f: impl FnMut(f64) -> f64, a: f64, b: f64, intervals: usize) -> f64 {
    let h = (b - a) / intervals as f64;
    let mut sum = f(a) + f(b);
    for i in 1..intervals {
        sum += (if i % 2 == 0 { 2.0 } else { 4.0 }) * f(a + i as f64 * h);
    }
    h * sum / 3.0
}
// Independent wavelength form with exact defining SI constants. Normalize its
// bolometric constant to the original conduction owner's rounded sigma.
fn independent_planck_log(wavelength: f64, temperature: f64) -> f64 {
    let h: f64 = 6.626_070_15e-34;
    let c = 299_792_458.0;
    let k: f64 = 1.380_649e-23;
    let sigma = 2.0 * core::f64::consts::PI.powi(5) * k.powi(4) / (15.0 * h.powi(3) * c * c);
    (2.0 * core::f64::consts::PI * h * c * c / wavelength.powi(4))
        / (h * c / (k * wavelength * temperature)).exp_m1()
        * (STEFAN_BOLTZMANN_W_M2_K4 / sigma)
}
fn heat(report: &SpectralRadiationReport<'_>) -> [f64; 4] {
    [
        report.surface_heat_power[0].value(),
        report.surface_heat_power[1].value(),
        report.gas_heat_power[0].value(),
        report.gas_heat_power[1].value(),
    ]
}
fn conserved(report: &SpectralRadiationReport<'_>) {
    close(heat(report).iter().sum(), report.energy_residual_w, 1e-10);
    assert!(report.energy_residual_w.abs() <= 1e-7);
    assert!(
        report
            .surface_heat_error_estimate_w
            .iter()
            .chain(report.gas_heat_error_estimate_w.iter())
            .all(|e| *e <= controls().heat_tolerance_w)
    );
}

#[test]
fn g1_planck_shape_and_omitted_tails_match_independent_integrals() {
    for t in [300.0, 900.0, 1800.0] {
        for wavelength in [0.2e-6, 2e-6, 10e-6, 100e-6] {
            let independent = independent_planck_log(wavelength, t);
            close(
                planck_log_emission(wavelength, t).unwrap(),
                independent,
                4e-13 * independent.max(1.0),
            );
        }
        let integrated = simpson(
            |s| independent_planck_log(s.exp(), t),
            1e-7_f64.ln(),
            1e-3_f64.ln(),
            16_384,
        );
        let coverage = coverage_at(
            [Length::new(1e-7), Length::new(1e-3)],
            Temperature::new(t),
            true,
        )
        .unwrap();
        let total = STEFAN_BOLTZMANN_W_M2_K4 * t.powi(4);
        close(
            integrated + coverage.omitted_blackbody_w_m2,
            total,
            total * 3e-13,
        );
    }
}

#[test]
fn g1_selective_surfaces_match_independent_planck_exchange_with_unlike_spectra() {
    let enclosure = plates();
    let edges = [1e-7, 5e-6, 1e-3];
    let values = [[0.9, 0.15], [0.2, 0.85]];
    let optics = SpectralOptics::new_with_surface_bands(
        (0..2)
            .map(|i| band(edges[i], edges[i + 1], [0.0; 2]))
            .collect(),
        [None; 2],
        values.map(|row| row.map(band_emissivity).to_vec()).to_vec(),
    )
    .unwrap();
    let mut work = enclosure.spectral_workspace().unwrap();
    let expected: f64 = (0..2)
        .map(|i| {
            let resistance = 1.0 / values[i][0] + 1.0 / values[i][1] - 1.0;
            1.3 / resistance
                * simpson(
                    |s| {
                        independent_planck_log(s.exp(), 1200.0)
                            - independent_planck_log(s.exp(), 300.0)
                    },
                    edges[i].ln(),
                    edges[i + 1].ln(),
                    16_384,
                )
        })
        .sum();
    for (temperatures, sign) in [([1200.0, 300.0], 1.0), ([300.0, 1200.0], -1.0)] {
        let report = enclosure
            .evaluate_spectral(
                &temperatures.map(Temperature::new),
                &[],
                [Temperature::new(500.0); 2],
                &optics,
                controls(),
                &mut work,
                || false,
            )
            .unwrap();
        close(report.surface_heat_power[1].value(), sign * expected, 2e-7);
        close(report.surface_heat_power[0].value(), -sign * expected, 2e-7);
        assert_eq!(report.gas_heat_power, [Power::new(0.0); 2]);
        assert!(report.surface_emissivities.is_empty());
        assert_eq!(report.optics, optics);
        assert_eq!(report.factorizations, 2);
        conserved(&report);
    }
    // A single gray value cannot substitute for selective absorption of the hot
    // source and emission of the cold surface. This difference is physical, not
    // a quadrature refinement effect.
    let gray = enclosure
        .evaluate_spectral(
            &[Temperature::new(1200.0), Temperature::new(300.0)],
            &[emissivity(0.55), emissivity(0.5)],
            [Temperature::new(500.0); 2],
            &transparent(),
            controls(),
            &mut work,
            || false,
        )
        .unwrap();
    assert!((gray.surface_heat_power[1].value() - expected).abs() > 1000.0);
    let equilibrium = enclosure
        .evaluate_spectral(
            &[Temperature::new(600.0); 2],
            &[],
            [Temperature::new(600.0); 2],
            &optics,
            controls(),
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(heat(&equilibrium), [0.0; 4]);
}

#[test]
fn g3_constant_surface_bands_recover_gray_surfaces_with_participating_gas_and_soot() {
    let enclosure = plates();
    let soot = soot();
    for population in [
        None,
        Some(SootPopulation {
            model: &soot,
            volume_fraction: 2e-6,
        }),
    ] {
        let gas_bands = vec![band(1e-7, 5e-6, [0.2, 0.9]), band(5e-6, 1e-3, [0.8, 0.1])];
        let gray =
            SpectralOptics::new(gas_bands.clone(), [population, None], "gray fixture").unwrap();
        let spectral = SpectralOptics::new_with_surface_bands(
            gas_bands,
            [population, None],
            vec![vec![band_emissivity(0.55), band_emissivity(0.8)]; 2],
        )
        .unwrap();
        let mut work = enclosure.spectral_workspace().unwrap();
        let temperatures = [Temperature::new(900.0), Temperature::new(300.0)];
        let gas = [Temperature::new(500.0), Temperature::new(700.0)];
        let expected = enclosure
            .evaluate_spectral(
                &temperatures,
                &[emissivity(0.55), emissivity(0.8)],
                gas,
                &gray,
                controls(),
                &mut work,
                || false,
            )
            .unwrap();
        let actual = enclosure
            .evaluate_spectral(
                &temperatures,
                &[],
                gas,
                &spectral,
                controls(),
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(heat(&actual), heat(&expected));
        assert_eq!(
            actual.surface_radiosity_w_m2,
            expected.surface_radiosity_w_m2
        );
        assert_eq!(actual.factorizations, expected.factorizations);
        conserved(&actual);
    }
}

#[test]
fn spectral_surface_sources_refuse_conflicting_axes_and_out_of_domain_temperatures() {
    let range = [Temperature::new(200.0), Temperature::new(600.0)];
    for value in [-0.1, 1.1, f64::NAN] {
        assert!(BandSurfaceEmissivity::new(value, range, "source").is_err());
    }
    assert!(BandSurfaceEmissivity::new(0.5, range, " ").is_err());
    assert!(SpectralOptics::new_with_surface_bands(
        vec![band(1e-7, 1e-3, [0.0; 2])], [None; 2], vec![],
    ).is_err());
    let bounded = BandSurfaceEmissivity::new(0.5, range, "bounded spectral source").unwrap();
    assert_eq!(bounded.value(), 0.5);
    assert_eq!(bounded.temperature_range(), range);
    assert_eq!(bounded.source(), "bounded spectral source");
    let optics = SpectralOptics::new_with_surface_bands(
        vec![band(1e-7, 1e-3, [0.0; 2])],
        [None; 2],
        vec![vec![bounded.clone(); 2]],
    )
    .unwrap();
    let enclosure = plates();
    let mut work = enclosure.spectral_workspace().unwrap();
    assert!(matches!(
        enclosure.evaluate_spectral(
            &[Temperature::new(900.0), Temperature::new(300.0)],
            &[],
            [Temperature::new(500.0); 2],
            &optics,
            controls(),
            &mut work,
            || false,
        ),
        Err(RadiationError::Invalid(
            "surface temperature is outside spectral emissivity validity"
        ))
    ));
    assert!(matches!(
        enclosure.evaluate_spectral(
            &[Temperature::new(500.0); 2],
            &[emissivity(0.5), emissivity(0.5)],
            [Temperature::new(500.0); 2],
            &optics,
            controls(),
            &mut work,
            || false,
        ),
        Err(RadiationError::Invalid(
            "spectral surface bands and gray emissivity cards cannot both own a boundary"
        ))
    ));
    let wrong_axis = SpectralOptics::new_with_surface_bands(
        optics.bands().to_vec(),
        [None; 2],
        vec![vec![bounded; 3]],
    )
    .unwrap();
    assert!(matches!(
        enclosure.evaluate_spectral(
            &[Temperature::new(500.0); 2],
            &[],
            [Temperature::new(500.0); 2],
            &wrong_axis,
            controls(),
            &mut work,
            || false,
        ),
        Err(RadiationError::Invalid(
            "spectral surface axis differs from enclosure"
        ))
    ));
    let replay = enclosure
        .evaluate_spectral(
            &[Temperature::new(500.0); 2],
            &[],
            [Temperature::new(500.0); 2],
            &optics,
            controls(),
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(heat(&replay), [0.0; 4]);
}

#[test]
fn g1_fixed_optical_bands_recover_original_gray_operator_and_reuse_factors() {
    let enclosure = plates();
    let temperatures = [900.0, 300.0, 500.0, 700.0];
    let emissivities = [emissivity(0.55), emissivity(0.8)];
    let optics = SpectralOptics::new(
        vec![band(1e-7, 1e-3, [0.2, 0.9])],
        [None; 2],
        "analytic gray boundary fixture",
    )
    .unwrap();
    let report = enclosure
        .evaluate_spectral(
            &temperatures[..2]
                .iter()
                .copied()
                .map(Temperature::new)
                .collect::<Vec<_>>(),
            &emissivities,
            [Temperature::new(500.0), Temperature::new(700.0)],
            &optics,
            controls(),
            &mut enclosure.spectral_workspace().unwrap(),
            || false,
        )
        .unwrap();
    // A fixed linear transfer operator commutes with integration. The gray
    // oracle receives the independent finite-window blackbody integrals.
    let equivalent = temperatures.map(|t| {
        let integrated = simpson(
            |s| independent_planck_log(s.exp(), t),
            1e-7_f64.ln(),
            1e-3_f64.ln(),
            16_384,
        );
        Temperature::new((integrated / STEFAN_BOLTZMANN_W_M2_K4).sqrt().sqrt())
    });
    let gray = enclosure
        .evaluate(
            &equivalent[..2],
            &emissivities,
            [equivalent[2], equivalent[3]],
            &optics.bands()[0].absorption,
            RadiationSolve::default(),
            &mut enclosure.workspace().unwrap(),
            || false,
        )
        .unwrap();
    for (a, b) in report
        .surface_heat_power
        .iter()
        .chain(report.gas_heat_power.iter())
        .zip(
            gray.surface_heat_power
                .iter()
                .chain(gray.gas_heat_power.iter()),
        )
    {
        close(a.value(), b.value(), 1e-7);
    }
    for (a, b) in report
        .surface_radiosity_w_m2
        .iter()
        .zip(&gray.surface_radiosity_w_m2)
    {
        close(*a, *b, 1e-7);
    }
    assert_eq!(report.factorizations, 1);
    assert_eq!(report.wavelength_evaluations, 0);
    assert_eq!(report.accepted_panels, 0);
    assert!(report.surface_heat_error_estimate_w.iter().all(|&e| e == 0.0));
    assert_eq!(report.gas_heat_error_estimate_w, [0.0; 2]);
    assert!(report.radiosity_error_estimate_w_m2.iter().all(|&e| e == 0.0));
    conserved(&report);

    let split = SpectralOptics::new(
        vec![band(1e-7, 8e-6, [0.2, 0.9]), band(8e-6, 1e-3, [0.2, 0.9])],
        [None; 2],
        "same gray fixture split into identical source bands",
    )
    .unwrap();
    let split_report = enclosure
        .evaluate_spectral(
            &[Temperature::new(900.0), Temperature::new(300.0)],
            &emissivities,
            [Temperature::new(500.0), Temperature::new(700.0)],
            &split,
            controls(),
            &mut enclosure.spectral_workspace().unwrap(),
            || false,
        )
        .unwrap();
    assert_eq!(split_report.factorizations, 2);
    assert_eq!(split_report.wavelength_evaluations, 0);
    assert_eq!(split_report.accepted_panels, 0);
    for (a, b) in heat(&report).into_iter().zip(heat(&split_report)) {
        close(a, b, 1e-7);
    }
}

#[test]
fn g1_absorbing_molecular_band_preserves_spectral_windows_and_kirchhoff_emission() {
    let enclosure = plates();
    let optics = SpectralOptics::new(
        vec![
            band(1e-7, 4e-6, [0.0; 2]),
            band(4e-6, 5e-6, [0.0, 1.7]),
            band(5e-6, 1e-3, [0.0; 2]),
        ],
        [None; 2],
        "authored one-band absorber and black boundaries, not a molecular dataset",
    )
    .unwrap();
    let report = enclosure
        .evaluate_spectral(
            &[Temperature::new(1000.0), Temperature::new(300.0)],
            &[emissivity(1.0), emissivity(1.0)],
            [Temperature::new(400.0), Temperature::new(700.0)],
            &optics,
            controls(),
            &mut enclosure.spectral_workspace().unwrap(),
            || false,
        )
        .unwrap();
    let emitted = [1000.0, 300.0, 700.0].map(|t| {
        simpson(
            |s| independent_planck_log(s.exp(), t),
            4e-6_f64.ln(),
            5e-6_f64.ln(),
            8192,
        )
    });
    let expected = 1.3 * -(-1.7_f64 * 0.8).exp_m1() * (emitted[0] + emitted[1] - 2.0 * emitted[2]);
    close(report.gas_heat_power[1].value(), expected, 1e-7);
    assert_eq!(report.gas_heat_power[0].value(), 0.0);
    assert_eq!(report.factorizations, 3);
    assert_eq!(report.wavelength_evaluations, 0);
    assert_eq!(report.accepted_panels, 0);
    conserved(&report);
}

#[test]
fn g1_rayleigh_soot_matches_independent_two_direction_spectral_transfer() {
    let enclosure = plates();
    let soot = soot();
    let populations = [1e-7, 1e-6];
    let optics = SpectralOptics::new(
        vec![band(1e-7, 1e-3, [0.0; 2])],
        populations.map(|volume_fraction| {
            Some(SootPopulation {
                model: &soot,
                volume_fraction,
            })
        }),
        "authored black plates; molecular absorption explicitly zero",
    )
    .unwrap();
    let temperatures = [400.0, 300.0, 600.0, 1000.0];
    let report = enclosure
        .evaluate_spectral(
            &[Temperature::new(400.0), Temperature::new(300.0)],
            &[emissivity(1.0), emissivity(1.0)],
            [Temperature::new(600.0), Temperature::new(1000.0)],
            &optics,
            controls(),
            &mut enclosure.spectral_workspace().unwrap(),
            || false,
        )
        .unwrap();
    let n: f64 = 1.75;
    let k: f64 = 1.03;
    let c0 = 36.0 * core::f64::consts::PI * n * k
        / ((n * n - k * k + 2.0).powi(2) + 4.0 * n * n * k * k);
    let mut expected = [0.0; 4];
    for channel in 0..4 {
        expected[channel] = simpson(
            |s| {
                let wavelength = s.exp();
                let blackbody = temperatures.map(|t| independent_planck_log(wavelength, t));
                let lengths = [0.4, 0.8];
                let a = [0, 1].map(|i| -(-c0 * populations[i] * lengths[i] / wavelength).exp_m1());
                let mut q = [0.0; 4];
                // Independent forward and reversed formal Beer--Lambert traversal.
                for (from, to, order) in [(0, 1, [0, 1]), (1, 0, [1, 0])] {
                    let mut incoming = blackbody[from];
                    for layer in order {
                        let absorbed = a[layer] * (incoming - blackbody[2 + layer]);
                        q[2 + layer] += absorbed;
                        incoming -= absorbed;
                    }
                    q[to] += incoming - blackbody[to];
                }
                1.3 * q[channel]
            },
            1e-7_f64.ln(),
            1e-3_f64.ln(),
            16_384,
        );
    }
    for (a, b) in heat(&report).into_iter().zip(expected) {
        close(a, b, 1e-7);
    }
    assert_eq!(report.factorizations, report.wavelength_evaluations);
    conserved(&report);

    let gray_absorption = [0, 1].map(|i| {
        soot.planck_from_volume_fraction(populations[i], Temperature::new(temperatures[2 + i]))
            .unwrap()
            .layer_absorption()
            .unwrap()
    });
    let gray = enclosure
        .evaluate(
            &[Temperature::new(400.0), Temperature::new(300.0)],
            &[emissivity(1.0), emissivity(1.0)],
            [Temperature::new(600.0), Temperature::new(1000.0)],
            &gray_absorption,
            RadiationSolve::default(),
            &mut enclosure.workspace().unwrap(),
            || false,
        )
        .unwrap();
    assert!(
        (gray.gas_heat_power[1].value() - report.gas_heat_power[1].value()).abs() > 10.0,
        "a local Planck mean must not be mistaken for nonequilibrium spectral filtering"
    );
}

#[test]
fn g0_equilibrium_is_exact_and_adjacent_temperature_heat_survives() {
    let enclosure = plates();
    let soot = soot();
    let optics = SpectralOptics::new(
        vec![band(1e-7, 1e-3, [0.0, 0.7])],
        [
            Some(SootPopulation {
                model: &soot,
                volume_fraction: 1e-6,
            }),
            None,
        ],
        "gray test boundary",
    )
    .unwrap();
    let report = enclosure
        .evaluate_spectral(
            &[Temperature::new(500.0); 2],
            &[emissivity(0.6), emissivity(0.8)],
            [Temperature::new(500.0); 2],
            &optics,
            controls(),
            &mut enclosure.spectral_workspace().unwrap(),
            || false,
        )
        .unwrap();
    assert_eq!(heat(&report), [0.0; 4]);
    let next = f64::from_bits(500.0_f64.to_bits() + 1);
    let transparent = transparent();
    let close_report = enclosure
        .evaluate_spectral(
            &[Temperature::new(next), Temperature::new(500.0)],
            &[emissivity(1.0), emissivity(1.0)],
            [Temperature::new(250.0), Temperature::new(1800.0)],
            &transparent,
            controls(),
            &mut enclosure.spectral_workspace().unwrap(),
            || false,
        )
        .unwrap();
    assert!(close_report.surface_heat_power[0].value() < 0.0);
    assert!(close_report.surface_heat_power[1].value() > 0.0);
    let expected = 1.3
        * STEFAN_BOLTZMANN_W_M2_K4
        * (next - 500.0)
        * (next + 500.0)
        * (next * next + 500.0 * 500.0);
    close(
        close_report.surface_heat_power[1].value() / expected,
        1.0,
        1e-5,
    );
    assert!(close_report.coverage[2..].iter().all(|c| !c.participating));
}

#[test]
fn g4_spectral_domain_work_cancellation_and_scratch_replay_are_explicit() {
    assert!(SpectralOptics::new(Vec::new(), [None; 2], "fixture").is_err());
    assert!(
        SpectralOptics::new(
            vec![band(1e-7, 1e-6, [0.0; 2]), band(2e-6, 1e-3, [0.0; 2])],
            [None; 2],
            "no hidden transparent gap"
        )
        .is_err()
    );
    assert!(SpectralOptics::new(vec![band(1e-7, 1e-3, [0.0; 2])], [None; 2], "").is_err());
    let too_many = (0..MAX_SPECTRAL_BANDS + 1)
        .map(|i| band(1e-7 * (i + 1) as f64, 1e-7 * (i + 2) as f64, [0.0; 2]))
        .collect();
    assert!(SpectralOptics::new(too_many, [None; 2], "fixture").is_err());
    let enclosure = plates();
    let temperatures = [Temperature::new(600.0), Temperature::new(300.0)];
    let gas = [Temperature::new(700.0), Temperature::new(900.0)];
    let e = [emissivity(0.6), emissivity(0.8)];
    let optics = transparent();
    let mut work = enclosure.spectral_workspace().unwrap();
    let baseline = enclosure
        .evaluate_spectral(
            &temperatures,
            &e,
            gas,
            &optics,
            controls(),
            &mut work,
            || false,
        )
        .unwrap();
    let small =
        SpectralOptics::new(vec![band(1e-6, 2e-6, [0.1; 2])], [None; 2], "fixture").unwrap();
    assert!(matches!(
        enclosure.evaluate_spectral(
            &temperatures,
            &e,
            gas,
            &small,
            controls(),
            &mut work,
            || false
        ),
        Err(RadiationError::SpectralTruncation { .. })
    ));
    // Fixed band optics need no point-wavelength budget at all: their Planck
    // sources are integrated exactly before one linear solve per band.
    let exact_under_tiny_budget = enclosure
        .evaluate_spectral(
            &temperatures,
            &e,
            gas,
            &optics,
            SpectralSolve {
                maximum_evaluations: 24,
                ..controls()
            },
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(exact_under_tiny_budget.wavelength_evaluations, 0);
    assert_eq!(exact_under_tiny_budget.accepted_panels, 0);

    // Wavelength-varying soot still owns the adaptive point-sampling budget.
    let soot_model = soot();
    let varying = SpectralOptics::new(
        vec![band(1e-7, 1e-3, [0.0; 2])],
        [
            Some(SootPopulation {
                model: &soot_model,
                volume_fraction: 1e-6,
            }),
            None,
        ],
        "varying soot fixture",
    )
    .unwrap();
    assert!(matches!(
        enclosure.evaluate_spectral(
            &temperatures,
            &e,
            gas,
            &varying,
            SpectralSolve {
                maximum_evaluations: 24,
                ..controls()
            },
            &mut work,
            || false
        ),
        Err(RadiationError::SpectralBudget {
            evaluations: 24,
            ..
        })
    ));
    for stop in [1, 10, 50] {
        let mut polls = 0;
        assert!(matches!(
            enclosure.evaluate_spectral(
                &temperatures,
                &e,
                gas,
                &varying,
                controls(),
                &mut work,
                || {
                    polls += 1;
                    polls >= stop
                }
            ),
            Err(RadiationError::Cancelled)
        ));
        assert_eq!(polls, stop);
    }
    let replay = enclosure
        .evaluate_spectral(
            &temperatures,
            &e,
            gas,
            &optics,
            controls(),
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(replay, baseline);
    let mut rates = [LayerRates::default(); 2];
    baseline.add_gas_heat(&mut rates).unwrap();
    assert_eq!(rates, [LayerRates::default(); 2]);
}

#[test]
fn perfect_reflectors_have_no_thermal_spectrum_and_gas_can_anchor_their_radiosity() {
    let enclosure = plates();
    let optics = SpectralOptics::new_with_surface_bands(
        vec![band(1e-7, 1e-3, [0.2, 0.7])],
        [None; 2],
        vec![vec![band_emissivity(0.0), band_emissivity(0.6)]],
    )
    .unwrap();
    let controls = SpectralSolve {
        maximum_omitted_blackbody_fraction: 3e-6,
        ..controls()
    };
    let mut work = enclosure.spectral_workspace().unwrap();
    let run = |reflector_temperature, work: &mut SpectralWorkspace| {
        enclosure
            .evaluate_spectral(
                &[
                    Temperature::new(reflector_temperature),
                    Temperature::new(450.0),
                ],
                &[],
                [Temperature::new(1200.0), Temperature::new(600.0)],
                &optics,
                controls,
                work,
                || false,
            )
            .unwrap()
    };
    let cold = run(200.0, &mut work);
    let hot = run(2000.0, &mut work);
    close(cold.surface_heat_power[0].value(), 0.0, 1e-8);
    assert_eq!(heat(&cold), heat(&hot));
    assert_eq!(cold.surface_radiosity_w_m2, hot.surface_radiosity_w_m2);
    assert!(!cold.coverage[0].participating);
    assert!(cold.coverage[0].longwave_fraction > controls.maximum_omitted_blackbody_fraction);
    assert!(cold.surface_heat_power[1].value() > 0.0);
    conserved(&cold);

    // Both boundaries may reflect perfectly when absorbing gas fixes the
    // radiation field. Only the two gas accounts exchange thermal energy.
    let gas_anchored = SpectralOptics::new_with_surface_bands(
        optics.bands().to_vec(),
        [None; 2],
        vec![vec![band_emissivity(0.0); 2]],
    )
    .unwrap();
    let report = enclosure
        .evaluate_spectral(
            &[Temperature::new(200.0), Temperature::new(2000.0)],
            &[],
            [Temperature::new(1200.0), Temperature::new(600.0)],
            &gas_anchored,
            controls,
            &mut work,
            || false,
        )
        .unwrap();
    for surface in &report.surface_heat_power {
        close(surface.value(), 0.0, 1e-8);
    }
    assert!(report.coverage[..2].iter().all(|c| !c.participating));
    assert!(report.gas_heat_power[0].value() < 0.0);
    assert!(report.gas_heat_power[1].value() > 0.0);
    conserved(&report);
}

#[test]
fn lossless_surface_components_are_refused_even_beside_an_emitting_component() {
    use fs_conduction::radiation::{ViewFactorEvidence, ViewFactorTolerance};

    // Reciprocal nonbinary factors form a closed three-mirror cycle, separate
    // from two ordinary emitting plates. A global Planck reference exists, but
    // cannot give the disconnected mirror component a unique radiosity.
    let chords = [(0, 1, 0.4), (0, 2, 0.7), (1, 2, 0.9), (3, 4, 2.1)];
    let mut areas = vec![0.0; 5];
    for &(a, b, area) in &chords {
        areas[a] += area;
        areas[b] += area;
    }
    let mut factors = vec![vec![0.0; 5]; 5];
    for &(a, b, area) in &chords {
        factors[a][b] = area / areas[a];
        factors[b][a] = area / areas[b];
    }
    let geometry = ViewFactorMatrix::admit(
        areas,
        factors,
        ViewFactorEvidence::Analytic {
            geometry: "authored disconnected reciprocal mirror-cycle fixture".into(),
        },
        ViewFactorTolerance::default(),
    )
    .unwrap();
    let enclosure = ParticipatingEnclosure::new(
        geometry,
        chords
            .into_iter()
            .map(|(a, b, area)| OpticalPath {
                surfaces: [a, b],
                exchange_area: Area::new(area),
                segments: Vec::new(),
            })
            .collect(),
        "authored disconnected reciprocal optical paths",
        1e-12,
        || false,
    )
    .unwrap();
    let optics = SpectralOptics::new_with_surface_bands(
        vec![band(1e-7, 1e-3, [0.0; 2])],
        [None; 2],
        vec![
            [0.0, 0.0, 0.0, 0.6, 0.7]
                .into_iter()
                .map(band_emissivity)
                .collect(),
        ],
    )
    .unwrap();
    assert!(matches!(
        enclosure.evaluate_spectral(
            &[500.0, 600.0, 700.0, 900.0, 300.0].map(Temperature::new),
            &[],
            [Temperature::new(500.0); 2],
            &optics,
            controls(),
            &mut enclosure.spectral_workspace().unwrap(),
            || false,
        ),
        Err(RadiationError::Invalid(
            "closed lossless surface component has indeterminate radiosity"
        ))
    ));

    let mirrors = SpectralOptics::new_with_surface_bands(
        optics.bands().to_vec(),
        [None; 2],
        vec![vec![band_emissivity(0.0); 2]],
    )
    .unwrap();
    let plates = plates();
    assert!(matches!(
        plates.evaluate_spectral(
            &[Temperature::new(500.0); 2],
            &[],
            [Temperature::new(500.0); 2],
            &mirrors,
            controls(),
            &mut plates.spectral_workspace().unwrap(),
            || false,
        ),
        Err(RadiationError::Invalid(
            "closed lossless enclosure has no thermal radiation anchor"
        ))
    ));
}
