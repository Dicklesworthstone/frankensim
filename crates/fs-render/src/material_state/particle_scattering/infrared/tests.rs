//! Original broadband material claims reach native transport and emission.

use super::*;
use crate::material_state::optical_mixture::deposition::ConstituentDepositionBudget;
use crate::material_state::optical_mixture::{
    MaterialMixtureEmission, MaterialOpticalComponent, MaterialOpticalMixture,
    MaterialOpticalMixtureGrid,
};
use crate::volumes::transport::{CellPalette, CellPaletteBudget, OpticalMixtureBudget};
use fs_qty::semantic::{Composition, CompositionBasis};
use std::sync::Arc;

const WAVELENGTHS: [f64; 3] = [8000.0, 10_000.0, 12_000.0];

fn knots() -> Vec<SpectralComplexIndexKnot> {
    WAVELENGTHS
        .into_iter()
        .enumerate()
        .map(|(index, wavelength_nm)| SpectralComplexIndexKnot {
            wavelength_nm,
            real_index_property: PropertyKey::new(format!("infrared_n_{index}"), Dims::NONE),
            extinction_index_property: PropertyKey::with_quantity(
                format!("infrared_k_{index}"),
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::OpticalExtinctionCoefficient,
                    ValueForm::Static,
                )),
            ),
        })
        .collect()
}

fn infrared_card(lossy: bool) -> MaterialCard {
    let mut values = vec![(
        PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
        1000.0,
        900.0,
    )];
    for knot in knots() {
        let fraction = (knot.wavelength_nm - 8000.0) / 4000.0;
        let n = 1.5 - 0.05 * fraction;
        let k = if lossy { 0.02 + 0.01 * fraction } else { 0.0 };
        values.push((knot.real_index_property, n, n + 0.03));
        values.push((knot.extinction_index_property, k, 1.5 * k));
    }
    let mut claims = ClaimSet::new();
    for (key, cold, hot) in values {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "thermal".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, cold), (400.0, hot)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("thermal", temperature_quantity(), 300.0, 400.0)
                    .with_quantity("pressure", pressure_quantity(), 100_000.0, 100_000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: Vec::new(),
                provenance: Provenance {
                    source: "synthetic infrared particle dispersion; not measured stock".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("synthetic broadband particles; lossy={lossy}"),
            phase: "solid".into(),
            process: "explicit complex-index wavelength bindings".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn infrared_source(card: &MaterialCard, temperature: f64) -> SpectralComplexIndexStatePoint {
    let knots = knots();
    let mut requirements = vec![
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    ];
    for knot in &knots {
        requirements.push(
            ScalarPropertyRequirement::try_with_key(
                &knot.real_index_property,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        );
        requirements.push(
            ScalarPropertyRequirement::try_with_key(
                &knot.extinction_index_property,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
        );
    }
    let resolved = resolve_material_state_point(
        card,
        &point(temperature),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    SpectralComplexIndexStatePoint::try_from_resolved(resolved, &knots).unwrap()
}

fn infrared_binding(source: &SpectralComplexIndexStatePoint) -> MaterialParticleScattering {
    MaterialParticleScattering::try_spectral(
        source.resolved(),
        source,
        &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
        ParticleSizeDistribution::try_new(&[(2e-9, 0.7), (4e-9, 0.3)]).unwrap(),
        1.0,
    )
    .unwrap()
}

fn infrared_spectrum(
    cx: &Cx<'_>,
    source: &SpectralComplexIndexStatePoint,
) -> MaterialParticleSpectrum {
    infrared_binding(source)
        .source_spectrum(
            cx,
            &WAVELENGTHS,
            MiePopulationBudget::default(),
            0.003,
            2048,
        )
        .unwrap()
}

#[test]
fn g1_original_infrared_material_and_temperature_reach_rayleigh_absorption() {
    with_cx(|cx, _| {
        let card = infrared_card(true);
        let cold = infrared_source(&card, 300.0);
        let hot = cold.at_query_point(&card, &point(400.0)).unwrap();
        for source in [&cold, &hot] {
            let binding = infrared_binding(source);
            let prepared = infrared_spectrum(cx, source);
            let warm = (prepared.accepted_temperature_k().unwrap() - 300.0) / 100.0;
            let density = 1000.0 - 100.0 * warm;
            assert_eq!(
                binding.optical_state_identity(),
                source.resolved().identity()
            );
            assert_eq!(binding.particle_density_kg_m3(), density);
            assert_eq!(
                prepared.spectrum().absorption.wavelength_bounds_nm(),
                (8000.0, 12_000.0)
            );
            assert_eq!(
                prepared.spectrum().scattering.wavelength_bounds_nm(),
                (8000.0, 12_000.0)
            );
            for wavelength in WAVELENGTHS {
                let fraction = (wavelength - 8000.0) / 4000.0;
                let n = 1.5 - 0.05 * fraction + 0.03 * warm;
                let k = (0.02 + 0.01 * fraction) * (1.0 + 0.5 * warm);
                let index = binding.complex_index(wavelength).unwrap();
                close(index.re, n, 2e-15);
                close(index.im, k, 2e-15);
                // Independent quasistatic dipole absorption per particle mass.
                // Im[(m²-1)/(m²+2)] = 6*n*k / |m²+2|², vacuum host.
                let denominator = (n * n - k * k + 2.0).powi(2) + (2.0 * n * k).powi(2);
                let expected = 36.0 * core::f64::consts::PI * n * k
                    / (density * wavelength * 1e-9 * denominator);
                close(
                    prepared.spectrum().absorption.at(wavelength).unwrap(),
                    expected,
                    2e-5,
                );
            }
            let midpoint = binding.complex_index(9000.0).unwrap();
            close(midpoint.re, 1.4875 + 0.03 * warm, 2e-15);
            close(midpoint.im, 0.0225 * (1.0 + 0.5 * warm), 2e-15);
        }
        assert_ne!(
            infrared_binding(&cold).content_identity(),
            infrared_binding(&hot).content_identity()
        );
        // The map is part of the constitutive declaration even with identical
        // source receipts. Changing its wavelength assignment cannot alias it.
        let mut remapped = knots();
        remapped[1].wavelength_nm = 11_000.0;
        let remapped =
            SpectralComplexIndexStatePoint::try_from_resolved(cold.resolved().clone(), &remapped)
                .unwrap();
        assert_ne!(
            infrared_binding(&cold).content_identity(),
            infrared_binding(&remapped).content_identity()
        );
    });
}

#[test]
fn g1_infrared_material_cells_keep_original_temperatures_in_kirchhoff_transport() {
    with_cx(|cx, _| {
        let card = infrared_card(true);
        let cold = infrared_source(&card, 300.0);
        let hot = cold.at_query_point(&card, &point(400.0)).unwrap();
        let spectra = vec![infrared_spectrum(cx, &cold), infrared_spectrum(cx, &hot)];
        let originals = spectra.clone();
        let materials = CellPalette::try_new(
            cx,
            spectra,
            Arc::from([0_u32, 1]),
            CellPaletteBudget {
                max_models: 2,
                max_cells: 2,
            },
        )
        .unwrap();
        let density = [2e-5, 7e-5];
        let temperatures = [300.0, 400.0];
        let field = MaterialParticleCellGrid::try_new(
            cx,
            grid(&density),
            materials.clone(),
            Some(&temperatures),
            1e-4,
            budget(),
        )
        .unwrap();
        let transport = field
            .medium()
            .absorption_grid()
            .transport(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                22.0,
                &WAVELENGTHS,
            )
            .unwrap();
        for (lane, wavelength) in WAVELENGTHS.into_iter().enumerate() {
            let tau = [0, 1].map(|i| {
                10.0 * density[i] * originals[i].spectrum().absorption.at(wavelength).unwrap()
            });
            let planck = temperatures.map(|temperature| {
                crate::thermal_surface::planck_radiance_w_m2_sr_nm(temperature, wavelength).unwrap()
            });
            let expected = planck[0] * (-(-tau[0]).exp_m1())
                + (-tau[0]).exp() * planck[1] * (-(-tau[1]).exp_m1());
            assert!(expected > 0.0);
            close(transport.emission_w_m2_sr_nm[lane], expected, 5e-13);
        }
        assert!(matches!(
            MaterialParticleCellGrid::try_new(
                cx,
                grid(&density),
                materials,
                Some(&[300.0, 300.0]),
                1e-4,
                budget()
            ),
            Err(ParticleScatteringError::TemperatureMismatch { index: 1 })
        ));
    });
}

#[test]
fn g1_infrared_mixture_keeps_scattering_particles_out_of_constituent_heat() {
    with_cx(|cx, _| {
        let scattering = infrared_spectrum(cx, &infrared_source(&infrared_card(false), 300.0));
        let absorbing = infrared_spectrum(cx, &infrared_source(&infrared_card(true), 400.0));
        let components = [
            MaterialOpticalComponent::Particles(&scattering),
            MaterialOpticalComponent::Particles(&absorbing),
        ];
        let composition =
            Composition::new(CompositionBasis::MassFraction, vec![0.25, 0.75]).unwrap();
        let mixture = MaterialOpticalMixture::try_new(
            cx,
            &components,
            &composition,
            OpticalMixtureBudget { max_components: 2 },
        )
        .unwrap();
        let mixtures = CellPalette::try_new(
            cx,
            vec![mixture],
            Arc::from([0_u32]),
            CellPaletteBudget {
                max_models: 1,
                max_cells: 1,
            },
        )
        .unwrap();
        let density = [1e-4];
        let field = MaterialOpticalMixtureGrid::try_new_with_emission(
            cx,
            grid(&density),
            mixtures,
            None,
            MaterialMixtureEmission::ConstituentTemperatures,
            0.01,
            budget(),
        )
        .unwrap();
        let powers = [2.0, 3.0, 5.0];
        let report = field
            .deposit_beam_power_by_component(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                12.0,
                &WAVELENGTHS,
                &powers,
                1,
                ConstituentDepositionBudget {
                    max_rows: 2,
                    max_component_evaluations: 6,
                },
            )
            .unwrap();
        assert_eq!(report.cells().len(), 1);
        assert_eq!(report.cells()[0].component_index(), 1);
        assert_eq!(
            report.cells()[0].source_identity(),
            absorbing.binding_identity()
        );
        for (lane, wavelength) in WAVELENGTHS.into_iter().enumerate() {
            assert_eq!(
                scattering.spectrum().absorption.at(wavelength).unwrap(),
                0.0
            );
            let ka = 0.75 * absorbing.spectrum().absorption.at(wavelength).unwrap();
            let ks = 0.25 * scattering.spectrum().scattering.at(wavelength).unwrap()
                + 0.75 * absorbing.spectrum().scattering.at(wavelength).unwrap();
            let expected =
                powers[lane] * (-(-10.0 * density[0] * (ka + ks)).exp_m1()) * ka / (ka + ks);
            close(report.cells()[0].absorbed_power_w()[lane], expected, 1e-12);
            close(report.beam().absorbed_power_w[lane], expected, 1e-12);
            assert!(report.beam().scattered_power_w[lane] > 0.0);
        }
    });
}

#[test]
fn g0_broadband_particle_binding_retains_source_state_support_and_work_refusals() {
    let card = infrared_card(true);
    let cold = infrared_source(&card, 300.0);
    let hot = cold.at_query_point(&card, &point(400.0)).unwrap();
    assert!(matches!(
        MaterialParticleScattering::try_spectral(
            cold.resolved(),
            &hot,
            &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            sizes(),
            1.0
        ),
        Err(ParticleScatteringError::Material(
            MaterialOpticalError::StateMismatch
        ))
    ));
    assert!(cold.at_query_point(&card, &point(401.0)).is_err());
    with_cx(|cx, gate| {
        let binding = infrared_binding(&cold);
        for wavelength in [780.0, 7999.0, 12_001.0] {
            assert!(matches!(
                binding.complex_index(wavelength),
                Err(ParticleScatteringError::OpticalSource(_))
            ));
        }
        assert!(
            binding
                .visible_spectrum(
                    cx,
                    &[380.0, 780.0],
                    MiePopulationBudget::default(),
                    0.003,
                    2048
                )
                .is_err()
        );
        let visible = super::binding(&super::card(Family::Conductor), 300.0, 1.0);
        assert!(
            visible
                .source_spectrum(
                    cx,
                    &WAVELENGTHS,
                    MiePopulationBudget::default(),
                    0.003,
                    2048
                )
                .is_err()
        );
        let limited = MiePopulationBudget {
            max_bins: 1,
            ..MiePopulationBudget::default()
        };
        assert!(matches!(
            binding.source_spectrum(cx, &WAVELENGTHS, limited, 0.003, 2048),
            Err(ParticleScatteringError::Mie(MieError::Budget))
        ));
        gate.request();
        assert!(matches!(
            binding.source_spectrum(
                cx,
                &WAVELENGTHS,
                MiePopulationBudget::default(),
                0.003,
                2048
            ),
            Err(ParticleScatteringError::Mie(MieError::Execution))
        ));
    });
}
