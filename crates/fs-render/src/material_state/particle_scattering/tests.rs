//! Material/source -> dispersive population -> original native grid transport.

#[path = "indexed/tests.rs"]
mod indexed_cells;

use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyValue,
    Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::*;
use fs_qty::{Dims, semantic::SemanticType};

#[derive(Clone, Copy, Debug)]
enum Family {
    Conductor,
    Cauchy,
    Tabulated,
    Sellmeier,
    Rgb,
}

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("thermal", temperature_quantity(), temperature)
        .unwrap()
        .with("pressure", 100_000.0)
        .unwrap()
}

fn card(family: Family) -> MaterialCard {
    let dimless = QuantitySpec::dimensional(Dims::NONE);
    let mut rows = vec![(
        DENSITY_PROPERTY,
        QuantitySpec::dimensional(Density::DIMS),
        1000.0,
        900.0,
    )];
    match family {
        Family::Conductor => {
            for (i, wavelength) in VISIBLE_COMPLEX_IOR_WAVELENGTHS_NM.into_iter().enumerate() {
                let n = 1.6 - (wavelength - 380.0) * 0.0003;
                let k = 0.01 + (wavelength - 380.0) * 0.0001;
                rows.push((VISIBLE_COMPLEX_IOR_ETA_PROPERTIES[i], dimless, n, n + 0.03));
                rows.push((VISIBLE_COMPLEX_IOR_K_PROPERTIES[i], dimless, k, 1.5 * k));
            }
        }
        Family::Tabulated => {
            for (i, wavelength) in VISIBLE_COMPLEX_IOR_WAVELENGTHS_NM.into_iter().enumerate() {
                let n = 1.6 - (wavelength - 380.0) * 0.0003;
                rows.push((VISIBLE_DIELECTRIC_IOR_PROPERTIES[i], dimless, n, n + 0.03));
            }
        }
        Family::Sellmeier => {
            for i in 0..3 {
                let strength = if i == 0 { 1.0 } else { 0.0 };
                let resonance = if i == 0 { 1.0e-14 } else { 0.0 };
                rows.push((
                    VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[i],
                    dimless,
                    strength,
                    strength,
                ));
                rows.push((
                    VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES[i],
                    QuantitySpec::dimensional(Dims([2, 0, 0, 0, 0, 0])),
                    resonance,
                    resonance,
                ));
            }
        }
        Family::Cauchy | Family::Rgb => {
            rows.extend([
                (VISIBLE_DIELECTRIC_CAUCHY_A_PROPERTY, dimless, 1.3, 1.33),
                (
                    VISIBLE_DIELECTRIC_CAUCHY_B_M2_PROPERTY,
                    QuantitySpec::dimensional(Dims([2, 0, 0, 0, 0, 0])),
                    1.0e-14,
                    1.0e-14,
                ),
                (
                    VISIBLE_DIELECTRIC_CAUCHY_C_M4_PROPERTY,
                    QuantitySpec::dimensional(Dims([4, 0, 0, 0, 0, 0])),
                    0.0,
                    0.0,
                ),
            ]);
        }
    }
    if matches!(family, Family::Rgb) {
        for name in VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES {
            rows.push((name, dimless, 0.9, 0.9));
        }
        rows.push((
            VISIBLE_DIELECTRIC_REFERENCE_DISTANCE_M_PROPERTY,
            QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
            1.0,
            1.0,
        ));
    } else if !matches!(family, Family::Conductor) {
        for (i, wavelength) in VISIBLE_COMPLEX_IOR_WAVELENGTHS_NM.into_iter().enumerate() {
            let alpha = 300_000.0 + 1000.0 * (wavelength - 380.0);
            rows.push((
                VISIBLE_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i],
                VISIBLE_DIELECTRIC_ABSORPTION_QUANTITY,
                alpha,
                1.5 * alpha,
            ));
        }
    }
    let mut claims = ClaimSet::new();
    for (name, quantity, cold, hot) in rows {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve {
                    abscissa: "thermal".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, cold), (400.0, hot)],
                    dims: quantity.dims(),
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("thermal", temperature_quantity(), 300.0, 400.0)
                    .with("pressure", 100_000.0, 100_000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: Vec::new(),
                provenance: Provenance {
                    source: "synthetic dispersive-particle fixture, not measured stock".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "declared synthetic particle fixture".into(),
            phase: "solid".into(),
            process: format!("{family:?}"),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn accepted(card: &MaterialCard, t: f64) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(
        card,
        &point(t),
        &[ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn sizes() -> ParticleSizeDistribution {
    ParticleSizeDistribution::try_new(&[(60e-9, 0.8), (120e-9, 0.2)]).unwrap()
}

fn source(card: &MaterialCard, t: f64) -> VisibleOpticalStatePoint {
    resolve_visible_optical_state_point(card, &point(t), MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
}

fn binding(card: &MaterialCard, t: f64, host: f64) -> MaterialParticleScattering {
    MaterialParticleScattering::try_new(
        &accepted(card, t),
        &source(card, t),
        &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
        sizes(),
        host,
    )
    .unwrap()
}

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 72,
                    kernel_id: 9,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    });
}

fn close(got: f64, want: f64, relative: f64) {
    assert!(
        (got - want).abs() <= relative * want.abs().max(1.0e-300),
        "{got} != {want}"
    );
}

#[test]
fn g1_source_dispersion_and_temperature_drive_every_particle_wavelength() {
    with_cx(|cx, _| {
        for family in [
            Family::Conductor,
            Family::Cauchy,
            Family::Tabulated,
            Family::Sellmeier,
        ] {
            let card = card(family);
            for t in [300.0, 400.0] {
                let binding = binding(&card, t, 1.0);
                assert_eq!(
                    binding.particle_density_kg_m3(),
                    if t == 300.0 { 1000.0 } else { 900.0 }
                );
                for wavelength in [380.0, 455.0, 550.0, 780.0] {
                    let warm = (t - 300.0) / 100.0;
                    let n = match family {
                        Family::Conductor | Family::Tabulated => {
                            1.6 - (wavelength - 380.0) * 0.0003 + 0.03 * warm
                        }
                        Family::Cauchy => {
                            1.3 + 0.03 * warm + 1.0e-14 / (wavelength * 1.0e-9).powi(2)
                        }
                        Family::Sellmeier => (1.0
                            + (wavelength * 1.0e-9).powi(2)
                                / ((wavelength * 1.0e-9).powi(2) - 1.0e-14))
                            .sqrt(),
                        Family::Rgb => unreachable!(),
                    };
                    let k = if matches!(family, Family::Conductor) {
                        (0.01 + (wavelength - 380.0) * 0.0001) * (1.0 + 0.5 * warm)
                    } else {
                        (300_000.0 + 1000.0 * (wavelength - 380.0))
                            * (1.0 + 0.5 * warm)
                            * wavelength
                            * 1.0e-9
                            / (4.0 * core::f64::consts::PI)
                    };
                    let index = binding.complex_index(wavelength).unwrap();
                    close(index.re, n, 3.0e-14);
                    close(index.im, k, 3.0e-14);
                    let population: Vec<_> = sizes()
                        .bins()
                        .iter()
                        .map(|&(radius, weight)| {
                            (
                                MieSphere::try_new(
                                    radius,
                                    n,
                                    k,
                                    1.0,
                                    binding.particle_density_kg_m3(),
                                )
                                .unwrap(),
                                weight,
                            )
                        })
                        .collect();
                    let expected = MiePopulation::try_new(&population)
                        .unwrap()
                        .at_wavelength(cx, wavelength, MiePopulationBudget::default())
                        .unwrap();
                    let actual = binding
                        .at_wavelength(cx, wavelength, MiePopulationBudget::default())
                        .unwrap();
                    close(
                        actual.cross_sections().mass_scattering_m2_kg,
                        expected.cross_sections().mass_scattering_m2_kg,
                        2.0e-12,
                    );
                    close(
                        actual.cross_sections().mass_absorption_m2_kg,
                        expected.cross_sections().mass_absorption_m2_kg,
                        2.0e-12,
                    );
                }
            }
            assert_ne!(
                binding(&card, 300.0, 1.0).content_identity(),
                binding(&card, 400.0, 1.0).content_identity()
            );
        }
    });
}

fn budget() -> AbsorbingGridBudget {
    AbsorbingGridBudget {
        max_grid_cells: 2,
        max_steps_per_segment: 8,
    }
}

fn grid(data: &[f64]) -> VolumeGrid<'_> {
    VolumeGrid {
        dims: [data.len(), 1, 1],
        origin: [0.0, 0.0, 0.0],
        cell: [10.0, 1.0, 1.0],
        data,
    }
}

#[test]
fn g1_dispersive_material_spectrum_enters_original_two_cell_transport() {
    with_cx(|cx, _| {
        let binding = binding(&card(Family::Cauchy), 300.0, 1.0);
        let spectrum = binding
            .visible_spectrum(
                cx,
                &VISIBLE_COMPLEX_IOR_WAVELENGTHS_NM,
                MiePopulationBudget::default(),
                0.003,
                2048,
            )
            .unwrap();
        for wavelength in VISIBLE_COMPLEX_IOR_WAVELENGTHS_NM {
            let direct = binding
                .at_wavelength(cx, wavelength, MiePopulationBudget::default())
                .unwrap();
            assert_eq!(
                spectrum.spectrum().scattering.at(wavelength).unwrap(),
                direct.cross_sections().mass_scattering_m2_kg
            );
        }
        let density = [2e-5, 7e-5];
        let medium = spectrum
            .bind_ambient_grid(cx, grid(&density), None, 1e-4, budget())
            .unwrap();
        close(medium.maximum_volume_fraction(), 7e-8, 1e-15);
        let wavelengths = [400.0, 550.0, 760.0];
        let actual = medium
            .medium()
            .transmittance(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                22.0,
                &wavelengths,
            )
            .unwrap();
        for (i, wavelength) in wavelengths.into_iter().enumerate() {
            let extinction = spectrum.spectrum().absorption.at(wavelength).unwrap()
                + spectrum.spectrum().scattering.at(wavelength).unwrap();
            close(
                actual[i],
                (-10.0 * density.iter().sum::<f64>() * extinction).exp(),
                3e-14,
            );
        }
        assert!(!medium.medium().absorption_grid().has_thermal_emission());
        assert!(
            spectrum
                .bind_ambient_grid(cx, grid(&density), Some(&[300.0, 300.0]), 1e-4, budget())
                .is_ok()
        );
        assert!(matches!(
            spectrum.bind_ambient_grid(cx, grid(&density), Some(&[300.0, 301.0]), 1e-4, budget()),
            Err(ParticleScatteringError::TemperatureMismatch { index: 1 })
        ));
        assert!(matches!(
            spectrum.bind_ambient_grid(cx, grid(&density), None, 1e-9, budget()),
            Err(ParticleScatteringError::VolumeFraction { index: 0 })
        ));
    });
}

#[test]
fn g0_particle_binding_refuses_wrong_source_state_domain_and_cancelled_work() {
    let card = card(Family::Cauchy);
    let cold = accepted(&card, 300.0);
    let hot_optics = source(&card, 400.0);
    let key = PropertyKey::new(DENSITY_PROPERTY, Density::DIMS);
    assert!(matches!(
        MaterialParticleScattering::try_new(&cold, &hot_optics, &key, sizes(), 1.0),
        Err(ParticleScatteringError::Material(
            MaterialOpticalError::StateMismatch
        ))
    ));
    assert!(matches!(
        MaterialParticleScattering::try_new(
            &cold,
            &source(&card, 300.0),
            &PropertyKey::new("absent density", Density::DIMS),
            sizes(),
            1.0
        ),
        Err(ParticleScatteringError::Density)
    ));
    let rgb_card = self::card(Family::Rgb);
    assert!(matches!(
        MaterialParticleScattering::try_new(
            &accepted(&rgb_card, 300.0),
            &source(&rgb_card, 300.0),
            &key,
            sizes(),
            1.0
        ),
        Err(ParticleScatteringError::SpectralAbsorptionRequired)
    ));
    with_cx(|cx, gate| {
        let binding = binding(&card, 300.0, 1.33);
        let spectrum = binding
            .visible_spectrum(
                cx,
                &[380.0, 780.0],
                MiePopulationBudget::default(),
                0.003,
                2048,
            )
            .unwrap();
        assert!(matches!(
            spectrum.bind_ambient_grid(cx, grid(&[1e-5]), None, 0.01, budget()),
            Err(ParticleScatteringError::HostIndex)
        ));
        assert!(
            binding
                .at_wavelength(cx, 379.0, MiePopulationBudget::default())
                .is_err()
        );
        gate.request();
        assert!(matches!(
            binding.at_wavelength(cx, 550.0, MiePopulationBudget::default()),
            Err(ParticleScatteringError::Mie(MieError::Execution))
        ));
    });
}
