use super::*;
use crate::compartment::radiation::{
    OpticalPath, OpticalSegment, ParticipatingEnclosure, RadiationSolve,
};
use fs_conduction::radiation::{
    EMISSIVITY_DIMS, STEFAN_BOLTZMANN_W_M2_K4, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity,
    ViewFactorMatrix,
};
use fs_evidence::ValidityDomain;
use fs_gas::multicomponent::compartment::Layer;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_qty::Area;

fn domain() -> SootOpticalDomain {
    SootOpticalDomain {
        spectral_temperature: [Temperature::new(200.0), Temperature::new(2500.0)],
        maximum_volume_fraction: 1e-3,
    }
}
fn soot() -> RayleighSoot {
    RayleighSoot::new(
        1.75,
        1.03,
        Density::new(1800.0),
        domain(),
        "authored constant-index test population; dilute Rayleigh spheres",
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_spectral_absorption_matches_pinned_radcal_rayleigh_law_in_si_units() {
    let model = soot();
    let [n, k] = model.refractive_index();
    // Independently evaluate POD's formula with wavelength in micrometres,
    // including its explicit 1e6 conversion back to inverse SI metres.
    for micrometres in [0.5_f64, 2.0, 10.0, 100.0] {
        let fv = 2.1e-6;
        let expected = 36.0 * core::f64::consts::PI * n * k
            / micrometres
            / ((n * n - k * k + 2.0).powi(2) + (2.0 * n * k).powi(2))
            * fv
            * 1e6;
        let actual = model
            .spectral_coefficient(fv, Length::new(micrometres * 1e-6))
            .unwrap();
        close(actual.value(), expected, 8e-16);
    }
}

#[test]
fn g1_analytic_planck_mean_agrees_with_independent_wavelength_quadrature() {
    let model = soot();
    let fv = 1.9e-6;
    // Direct SI Planck spectrum in wavelength; Simpson integration in log λ.
    // This uses neither the zeta identity nor the production averaging factor.
    // λ in [1 nm, 0.1 m] retains all but <6e-12 of the 300 K denominator.
    let intervals = 8192;
    let lo = 1e-9_f64.ln();
    let hi = 0.1_f64.ln();
    let step = (hi - lo) / intervals as f64;
    for temperature in [300.0, 900.0, 1800.0, 2500.0] {
        let (mut numerator, mut denominator) = (0.0, 0.0);
        for i in 0..=intervals {
            let lambda = (lo + i as f64 * step).exp();
            let exponent =
                6.626_070_15e-34 * 299_792_458.0 / (lambda * 1.380_649e-23 * temperature);
            let intensity = if exponent > 700.0 {
                0.0
            } else {
                2.0 * 6.626_070_15e-34 * 299_792_458.0_f64.powi(2)
                    / lambda.powi(5)
                    / exponent.exp_m1()
            };
            let weight = if i == 0 || i == intervals {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let integrand = intensity * lambda * weight;
            denominator += integrand;
            numerator += integrand
                * model
                    .spectral_coefficient(fv, Length::new(lambda))
                    .unwrap()
                    .value();
        }
        let actual = model
            .planck_from_volume_fraction(fv, Temperature::new(temperature))
            .unwrap();
        close(actual.coefficient.value(), numerator / denominator, 8e-12);
    }
}

#[test]
fn g3_population_and_spectral_temperature_scale_without_gas_density_assumptions() {
    let model = soot();
    let base = model
        .planck_from_mass_concentration(Density::new(0.0036), Temperature::new(700.0))
        .unwrap();
    close(base.volume_fraction, 2e-6, 2e-16);
    assert_eq!(base.mass_concentration.value(), 0.0036);
    let double_population = model
        .planck_from_volume_fraction(4e-6, Temperature::new(700.0))
        .unwrap();
    let double_spectral_t = model
        .planck_from_volume_fraction(2e-6, Temperature::new(1400.0))
        .unwrap();
    close(
        double_population.coefficient.value(),
        2.0 * base.coefficient.value(),
        4e-16,
    );
    close(
        double_spectral_t.coefficient.value(),
        2.0 * base.coefficient.value(),
        4e-16,
    );
    assert_eq!(base.model(), &model);
    assert_eq!(model.particle_density(), Density::new(1800.0));
    assert_eq!(model.domain(), domain());
    let original = base.layer_absorption().unwrap();
    assert_eq!(original.coefficient(), base.coefficient);
    assert!(original.source().contains(model.source()));
    assert!(original.source().contains("spectral_temperature_K=700"));
    assert!(
        original
            .source()
            .contains("soot_mass_concentration_kg_m3=0.0036")
    );
}

#[test]
fn g0_zero_population_or_nonabsorbing_particles_are_exactly_transparent() {
    for model in [
        soot(),
        RayleighSoot::new(
            1.75,
            0.0,
            Density::new(1800.0),
            domain(),
            "nonabsorbing fixture",
        )
        .unwrap(),
    ] {
        let absent = model
            .planck_from_mass_concentration(Density::new(0.0), Temperature::new(200.0))
            .unwrap();
        assert_eq!(absent.coefficient.value(), 0.0);
        assert_eq!(absent.volume_fraction, 0.0);
        assert_eq!(
            absent.layer_absorption().unwrap().coefficient().value(),
            0.0
        );
        if model.refractive_index()[1] == 0.0 {
            let report = model
                .planck_from_volume_fraction(1e-6, Temperature::new(2500.0))
                .unwrap();
            assert_eq!(report.coefficient.value(), 0.0);
            assert_eq!(
                model
                    .spectral_coefficient(1e-6, Length::new(1e-6))
                    .unwrap()
                    .value(),
                0.0
            );
        }
    }
}

#[test]
fn g0_optical_domain_and_unrepresentable_inputs_refuse_without_clipping() {
    let model = soot();
    for fraction in [-1.0, 0.0010001, f64::INFINITY, f64::NAN] {
        assert!(
            model
                .planck_from_volume_fraction(fraction, Temperature::new(1000.0))
                .is_err()
        );
        assert!(
            model
                .spectral_coefficient(fraction, Length::new(1e-6))
                .is_err()
        );
    }
    for t in [0.0, 199.999, 2500.001, f64::INFINITY, f64::NAN] {
        assert!(
            model
                .planck_from_volume_fraction(0.0, Temperature::new(t))
                .is_err()
        );
    }
    for x in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(
            model
                .planck_from_mass_concentration(Density::new(x), Temperature::new(1000.0))
                .is_err()
        );
    }
    for wavelength in [0.0, -1.0, f64::INFINITY, f64::NAN, f64::from_bits(1)] {
        assert!(
            model
                .spectral_coefficient(1e-6, Length::new(wavelength))
                .is_err()
        );
    }
    for (n, k, rho) in [
        (0.0, 1.0, 1800.0),
        (1.0, -1.0, 1800.0),
        (1.0, 1.0, 0.0),
        (f64::NAN, 1.0, 1800.0),
        (1.0, f64::INFINITY, 1800.0),
        (1e200, 1.0, 1800.0),
    ] {
        assert!(RayleighSoot::new(n, k, Density::new(rho), domain(), "fixture").is_err());
    }
    assert!(RayleighSoot::new(1.0, 1.0, Density::new(1800.0), domain(), " ").is_err());
    let invalid = SootOpticalDomain {
        spectral_temperature: [Temperature::new(500.0), Temperature::new(400.0)],
        ..domain()
    };
    assert!(RayleighSoot::new(1.0, 1.0, Density::new(1800.0), invalid, "fixture").is_err());
    assert!(
        model
            .planck_from_mass_concentration(
                Density::new(f64::from_bits(1)),
                Temperature::new(1000.0)
            )
            .is_err()
    );
    assert!(
        model
            .planck_from_volume_fraction(domain().maximum_volume_fraction, Temperature::new(2500.0))
            .is_ok()
    );
}

fn emissivity(value: f64) -> SurfaceEmissivity {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 2500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "analytic gray plate".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture".into(),
            phase: "solid".into(),
            process: "gray plate".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    SurfaceEmissivity::from_card("fixture", &card, 300.0, SelectionPolicy::SingleClaimOnly).unwrap()
}

#[test]
fn g1_population_drives_original_enclosure_exchange_with_exact_pair_accounting() {
    let soot = soot();
    let enclosure = ParticipatingEnclosure::new(
        ViewFactorMatrix::infinite_parallel_plates(1.0).unwrap(),
        vec![OpticalPath {
            surfaces: [0, 1],
            exchange_area: Area::new(1.0),
            segments: vec![OpticalSegment {
                layer: Layer::Upper,
                length: Length::new(0.5),
            }],
        }],
        "authored one-chord fixture",
        1e-12,
        || false,
    )
    .unwrap();
    let temperatures = [Temperature::new(300.0); 2];
    let surfaces = [emissivity(1.0), emissivity(1.0)];
    let gas = [Temperature::new(300.0), Temperature::new(1000.0)];
    let mut work = enclosure.workspace().unwrap();
    let mut preceding_heat = -1.0;
    for concentration in [0.0, 0.00018, 0.0018, 0.018] {
        let upper = soot
            .planck_from_mass_concentration(Density::new(concentration), gas[1])
            .unwrap();
        let absorption = [
            soot.planck_from_volume_fraction(0.0, gas[0])
                .unwrap()
                .layer_absorption()
                .unwrap(),
            upper.layer_absorption().unwrap(),
        ];
        let report = enclosure
            .evaluate(
                &temperatures,
                &surfaces,
                gas,
                &absorption,
                RadiationSolve::default(),
                &mut work,
                || false,
            )
            .unwrap();
        let heat = report.surface_heat_power[0].value();
        let absorbed = -(-upper.coefficient.value() * 0.5).exp_m1();
        let expected =
            absorbed * STEFAN_BOLTZMANN_W_M2_K4 * (1000.0_f64.powi(4) - 300.0_f64.powi(4));
        if concentration == 0.0 {
            assert_eq!(heat, 0.0);
        } else {
            close(heat, expected, 1e-14);
        }
        assert!(heat > preceding_heat);
        preceding_heat = heat;
        assert_eq!(report.surface_heat_power[0], report.surface_heat_power[1]);
        assert_eq!(report.gas_heat_power[0].value(), 0.0);
        assert_eq!(report.gas_heat_power[1].value(), -2.0 * heat);
        assert_eq!(report.energy_residual_w, 0.0);
        assert_eq!(report.absorption, absorption);
    }
}


#[test]
fn g1_particle_opacity_matches_number_density_times_rayleigh_cross_sections() {
    let model = soot();
    let (fraction, radius, wavelength) = (2e-6, 2e-9_f64, 500e-9_f64);
    let optics = model
        .particle_opacity(fraction, Length::new(radius), Length::new(wavelength))
        .unwrap();
    let [n, k] = model.refractive_index();
    let pi = core::f64::consts::PI;
    let number_density = fraction / ((4.0 / 3.0) * pi * radius.powi(3));
    let contrast_squared = ((n * n - k * k - 1.0).powi(2) + (2.0 * n * k).powi(2))
        / ((n * n - k * k + 2.0).powi(2) + (2.0 * n * k).powi(2));
    let cross_section = (8.0 * pi / 3.0)
        * (2.0 * pi / wavelength).powi(4)
        * radius.powi(6)
        * contrast_squared;
    close(optics.scattering_per_m(), number_density * cross_section, 3e-15);
    assert_eq!(
        optics.absorption_per_m(),
        model.spectral_coefficient(fraction, Length::new(wavelength)).unwrap().value()
    );
    assert_eq!(
        optics.extinction_per_m(),
        optics.absorption_per_m() + optics.scattering_per_m()
    );
    assert!(core::ptr::eq(optics.model(), &model));
    assert_eq!(optics.radius().value(), radius);
    assert_eq!(optics.wavelength().value(), wavelength);
    assert_eq!(optics.volume_fraction(), fraction);
}

#[test]
fn g3_particle_opacity_respects_population_radius_and_wavelength_scaling() {
    let model = soot();
    let sample = |fraction, radius, wavelength| {
        model.particle_opacity(fraction, Length::new(radius), Length::new(wavelength)).unwrap()
    };
    let base = sample(2e-6, 2e-9, 1e-6);
    let doubled_mass = sample(4e-6, 2e-9, 1e-6);
    let doubled_radius = sample(2e-6, 4e-9, 1e-6);
    let doubled_wavelength = sample(2e-6, 2e-9, 2e-6);
    close(doubled_mass.absorption_per_m(), 2.0 * base.absorption_per_m(), 2e-15);
    close(doubled_mass.scattering_per_m(), 2.0 * base.scattering_per_m(), 2e-15);
    close(doubled_mass.extinction_per_m(), 2.0 * base.extinction_per_m(), 2e-15);
    close(doubled_radius.absorption_per_m(), base.absorption_per_m(), 2e-15);
    close(doubled_radius.scattering_per_m(), 8.0 * base.scattering_per_m(), 2e-15);
    close(doubled_wavelength.absorption_per_m(), base.absorption_per_m() / 2.0, 2e-15);
    close(doubled_wavelength.scattering_per_m(), base.scattering_per_m() / 16.0, 2e-15);
    let from_mass = model.particle_opacity_from_mass_concentration(
        Density::new(2e-6 * model.particle_density().value()),
        Length::new(2e-9),
        Length::new(1e-6),
    ).unwrap();
    close(from_mass.extinction_per_m(), base.extinction_per_m(), 2e-15);
}

#[test]
fn g3_primary_beam_transfer_partitions_energy_and_composes_paths() {
    let model = soot();
    let optics = model.particle_opacity(2e-6, Length::new(2e-9), Length::new(500e-9)).unwrap();
    let first = optics.direct_beam(Length::new(0.02)).unwrap();
    let second = optics.direct_beam(Length::new(0.03)).unwrap();
    let whole = optics.direct_beam(Length::new(0.05)).unwrap();
    for beam in [first, second, whole] {
        close(beam.transmittance + beam.absorbed_fraction + beam.scattered_fraction, 1.0, 3e-15);
        close(beam.opacity, beam.absorbed_fraction + beam.scattered_fraction, 3e-15);
        assert!(beam.scattered_fraction > 0.0 && beam.absorbed_fraction > 0.0);
    }
    close(whole.transmittance, first.transmittance * second.transmittance, 3e-15);
    close(
        whole.absorbed_fraction,
        first.absorbed_fraction + first.transmittance * second.absorbed_fraction,
        3e-15,
    );
    close(
        whole.scattered_fraction,
        first.scattered_fraction + first.transmittance * second.scattered_fraction,
        3e-15,
    );
}

#[test]
fn g0_primary_beam_resolves_thin_zero_and_opaque_limits() {
    let model = soot();
    let optics = model.particle_opacity(2e-6, Length::new(2e-9), Length::new(500e-9)).unwrap();
    let path = 1e-18 / optics.extinction_per_m();
    let thin = optics.direct_beam(Length::new(path)).unwrap();
    assert_eq!(thin.transmittance, 1.0);
    close(thin.opacity, 1e-18, 2e-15);
    close(thin.absorbed_fraction, optics.absorption_per_m() * path, 2e-15);
    close(thin.scattered_fraction, optics.scattering_per_m() * path, 2e-15);
    let zero = optics.direct_beam(Length::new(0.0)).unwrap();
    assert_eq!((zero.transmittance, zero.opacity, zero.absorbed_fraction, zero.scattered_fraction),
        (1.0, 0.0, 0.0, 0.0));
    assert!(optics.extinction_per_m() > 1.0);
    let opaque = optics.direct_beam(Length::new(f64::MAX)).unwrap();
    assert_eq!((opaque.transmittance, opaque.opacity), (0.0, 1.0));
    close(opaque.absorbed_fraction + opaque.scattered_fraction, 1.0, 2e-15);
}

#[test]
fn g0_nonabsorbing_particles_remove_beam_power_without_depositing_heat() {
    let model = RayleighSoot::new(
        1.5, 0.0, Density::new(1800.0), domain(), "lossless spheres",
    )
    .unwrap();
    let optics = model.particle_opacity(1e-4, Length::new(2e-9), Length::new(500e-9)).unwrap();
    assert_eq!(optics.absorption_per_m(), 0.0);
    assert!(optics.scattering_per_m() > 0.0);
    let beam = optics.direct_beam(Length::new(1.0 / optics.extinction_per_m())).unwrap();
    close(beam.transmittance, (-1.0_f64).exp(), 2e-15);
    assert_eq!(beam.absorbed_fraction, 0.0);
    assert_eq!(beam.scattered_fraction, beam.opacity);

    let matched = RayleighSoot::new(
        1.0, 0.0, Density::new(1800.0), domain(), "index-matched spheres",
    )
    .unwrap();
    for optics in [
        model.particle_opacity(0.0, Length::new(2e-9), Length::new(500e-9)).unwrap(),
        matched.particle_opacity(1e-4, Length::new(2e-9), Length::new(500e-9)).unwrap(),
    ] {
        assert_eq!(optics.extinction_per_m(), 0.0);
        let beam = optics.direct_beam(Length::new(f64::MAX)).unwrap();
        assert_eq!(
            (beam.transmittance, beam.opacity, beam.absorbed_fraction, beam.scattered_fraction),
            (1.0, 0.0, 0.0, 0.0),
        );
    }
}

#[test]
fn g0_particle_opacity_enforces_loading_size_resonance_and_path_domains() {
    let model = soot();
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(model.particle_opacity(1e-6, Length::new(invalid), Length::new(1e-6)).is_err());
        assert!(model.particle_opacity(1e-6, Length::new(1e-9), Length::new(invalid)).is_err());
    }
    for invalid in [-1.0, 2e-3, f64::NAN, f64::INFINITY] {
        assert!(model.particle_opacity(invalid, Length::new(1e-9), Length::new(1e-6)).is_err());
    }
    assert!(model.particle_opacity(1e-6, Length::new(1e-6), Length::new(1e-6)).is_err());
    let high_index = RayleighSoot::new(
        10.0, 0.0, Density::new(1800.0), domain(), "high-index limit",
    )
    .unwrap();
    assert!(high_index.particle_opacity(
        1e-6, Length::new(0.02e-6 / core::f64::consts::TAU), Length::new(1e-6),
    ).is_err());
    let resonant = RayleighSoot::new(
        1e-10, 2.0_f64.sqrt(), Density::new(1800.0), domain(), "resonance limit",
    )
    .unwrap();
    assert!(resonant.particle_opacity(
        1e-6, Length::new(0.01e-6 / core::f64::consts::TAU), Length::new(1e-6),
    ).is_err());
    let optics = model.particle_opacity(1e-6, Length::new(1e-9), Length::new(1e-6)).unwrap();
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(optics.direct_beam(Length::new(invalid)).is_err());
        assert!(model.particle_opacity_from_mass_concentration(
            Density::new(invalid), Length::new(1e-9), Length::new(1e-6),
        ).is_err());
    }
}
