use super::*;
use crate::conductor::{
    CauchyCoatingLayer, CoatingIorSample, CoherentAbsorbingCoating,
    CoherentDielectricCoating, CoherentPeriodicCoating,
};
use crate::dielectric::{DielectricCoating, DielectricEvent};
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::{PolarizationFrame, PolarizedBeam};
use fs_evidence::ValidityDomain;
use fs_geom::{Point3, Vec3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECTRAL_CONDUCTOR_ETA_PROPERTIES,
    SPECTRAL_CONDUCTOR_K_PROPERTIES, SPECTRAL_CONDUCTOR_K_QUANTITY,
    SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
    SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY, SPECTRAL_DIELECTRIC_IOR_PROPERTIES,
    SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES, ScalarAdmissibility,
    ScalarPropertyRequirement, SpectralOpticalStatePoint,
    VISIBLE_COMPLEX_IOR_ETA_PROPERTIES, VISIBLE_COMPLEX_IOR_K_PROPERTIES,
    VISIBLE_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES, VISIBLE_DIELECTRIC_ABSORPTION_QUANTITY,
    VISIBLE_DIELECTRIC_IOR_PROPERTIES, resolve_material_state_point,
    resolve_spectral_conductor_state_point, resolve_spectral_tabulated_dielectric_state_point,
    resolve_visible_conductor_state_point, resolve_visible_tabulated_dielectric_state_point,
};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};

const FINISH: ContentHash = ContentHash([191; 32]);

fn point(temperature_k: f64) -> QueryPoint {
    QueryPoint::new().with("T", temperature_k).unwrap()
}

fn claim(name: &str, quantity: QuantitySpec, cold: f64, hot: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::with_quantity(name, quantity),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, cold), (600.0, hot)],
            dims: quantity.dims(),
        },
        validity: ValidityDomain::unconstrained().with("T", 300.0, 600.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(),
        provenance: Provenance {
            source: format!("synthetic optical-finish integration fixture: {name}"),
            license: "internal-test-use".into(),
            artifact: None,
        },
    }
}

fn card(spectral_scale: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let index = QuantitySpec::dimensional(Dims::NONE);
    for (i, wavelength_nm) in VISIBLE_COMPLEX_IOR_WAVELENGTHS_NM.into_iter().enumerate() {
        let eta = 0.5 + 0.03 * i as f64;
        let k = 3.0 + 0.1 * i as f64;
        let n = 1.5 + 0.005 * i as f64;
        let absorption = 0.2 + 0.02 * i as f64;
        for source in [
            claim(VISIBLE_COMPLEX_IOR_ETA_PROPERTIES[i], index, eta, eta + 0.5),
            claim(VISIBLE_COMPLEX_IOR_K_PROPERTIES[i], index, k, k + 1.0),
            claim(VISIBLE_DIELECTRIC_IOR_PROPERTIES[i], index, n, n + 0.3),
            claim(
                VISIBLE_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i],
                VISIBLE_DIELECTRIC_ABSORPTION_QUANTITY, absorption, absorption + 0.6,
            ),
            claim(
                SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[i],
                QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                wavelength_nm * spectral_scale * 1e-9,
                wavelength_nm * spectral_scale * 1e-9,
            ),
            claim(SPECTRAL_CONDUCTOR_ETA_PROPERTIES[i], index, eta, eta + 0.5),
            claim(
                SPECTRAL_CONDUCTOR_K_PROPERTIES[i], SPECTRAL_CONDUCTOR_K_QUANTITY,
                k, k + 1.0,
            ),
            claim(SPECTRAL_DIELECTRIC_IOR_PROPERTIES[i], index, n, n + 0.3),
            claim(
                SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i],
                SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY, absorption, absorption + 0.6,
            ),
        ] {
            claims.insert_claim(source).unwrap();
        }
    }
    claims.insert_claim(claim(
        DENSITY_PROPERTY, QuantitySpec::dimensional(Density::DIMS), 2500.0, 2400.0,
    )).unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic shared optical material".into(),
            phase: "solid".into(),
            process: "native finish refresh verification".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    ).unwrap()
}

fn state(card: &MaterialCard, temperature_k: f64) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(
        card,
        &point(temperature_k),
        &[ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY, Density::DIMS, ScalarAdmissibility::StrictlyPositive,
        ).unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()
}

// Only dispatch between the two original source owners; no synthetic optical law.
enum Original {
    Visible(VisibleOpticalStatePoint),
    Spectral(SpectralOpticalStatePoint),
}

impl Original {
    fn refresh(
        &self,
        binding: MaterialOpticalBinding,
        card: &MaterialCard,
        accepted: &ResolvedMaterialStatePoint,
    ) -> Result<(MaterialOpticalBinding, Self), MaterialOpticalRefreshError> {
        match self {
            Self::Visible(source) => binding
                .at_material_state(source, card, accepted)
                .map(|(binding, source)| (binding, Self::Visible(source))),
            Self::Spectral(source) => binding
                .at_spectral_material_state(source, card, accepted)
                .map(|(binding, source)| (binding, Self::Spectral(source))),
        }
    }
}

fn binding(
    card: &MaterialCard,
    spectral: bool,
    conductor: bool,
    max_gap_nm: f64,
) -> (MaterialOpticalBinding, Original) {
    let accepted = state(card, 300.0);
    let selection = MaterialPropertySelection::SingleClaimOnly;
    if spectral {
        if conductor {
            let source = resolve_spectral_conductor_state_point(
                card, &accepted.query(), selection,
            ).unwrap();
            let binding = MaterialOpticalBinding::try_spectral_conductor(
                &accepted, &source, max_gap_nm, 0.15, FINISH,
            ).unwrap();
            (binding, Original::Spectral(SpectralOpticalStatePoint::Conductor(source)))
        } else {
            let source = resolve_spectral_tabulated_dielectric_state_point(
                card, &accepted.query(), selection,
            ).unwrap();
            let binding = MaterialOpticalBinding::try_spectral_tabulated_dielectric(
                &accepted, &source, max_gap_nm, Some(0.15), FINISH,
            ).unwrap();
            (binding, Original::Spectral(SpectralOpticalStatePoint::TabulatedDielectric(source)))
        }
    } else {
        let source = if conductor {
            VisibleOpticalStatePoint::Conductor(resolve_visible_conductor_state_point(
                card, &accepted.query(), selection,
            ).unwrap())
        } else {
            VisibleOpticalStatePoint::TabulatedDielectric(
                resolve_visible_tabulated_dielectric_state_point(
                    card, &accepted.query(), selection,
                ).unwrap(),
            )
        };
        let binding = MaterialOpticalBinding::try_visible(
            &accepted, &source, 0.15, FINISH,
        ).unwrap();
        (binding, Original::Visible(source))
    }
}

fn film_source() -> ConductorSource {
    ConductorSource::try_new(
        ContentHash([192; 32]), ConductorDataStatus::Representative,
    ).unwrap()
}

fn periodic(thickness_m: f64) -> CoherentPeriodicCoating {
    CoherentPeriodicCoating::try_new(
        &[
            CauchyCoatingLayer::try_new(
                CauchyIor::try_new(1.35, 0.001, 0.0).unwrap(), 0.0, thickness_m,
            ).unwrap(),
            CauchyCoatingLayer::try_new(
                CauchyIor::try_new(1.9, 0.002, 0.0).unwrap(), 0.03, thickness_m * 0.7,
            ).unwrap(),
        ],
        2,
        film_source(),
    ).unwrap()
}

fn air() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_new(1.0, 0.0, 0.0).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    )
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

// A refresh re-queries with PinnedByProperty (state_point.rs:1085-1093, 5c6b323c6e) and the optical identity hashes the
// selection policy, so a binding resolved with SingleClaimOnly (`cold`) and the same state reached by a refresh differ in
// that identity by design (fs-material's refresh tests, a044e818f3, measured the same on another owner). The round trip
// is therefore compared with the same pinned path at the same state, and with the physical values of `cold`.
fn assert_pinned_round_trip(
    source: &Original,
    cold: MaterialOpticalBinding,
    back: MaterialOpticalBinding,
    card: &MaterialCard,
) {
    let pinned = source.refresh(cold, card, &state(card, 300.0)).unwrap().0;
    assert_eq!(
        back, pinned,
        "the returned original source supports the next accepted endpoint"
    );
    match (back.appearance(), cold.appearance()) {
        (
            OpticalAppearance::Conductor {
                optics: back,
                surface: back_surface,
            },
            OpticalAppearance::Conductor {
                optics: cold,
                surface: cold_surface,
            },
        ) => {
            // The bulk source identity hashes the resolved bundle with its selection policy, so it is the pinned
            // path's (compared above through the whole binding); every physical value is the cold one.
            assert_eq!(back.samples(), cold.samples());
            assert_eq!(back.max_knot_span_nm(), cold.max_knot_span_nm());
            assert_eq!(back.coating(), cold.coating());
            assert_eq!(back.absorbing_coating(), cold.absorbing_coating());
            assert_eq!(back.periodic_coating(), cold.periodic_coating());
            assert_eq!(back.source().status(), cold.source().status());
            assert_eq!(back_surface, cold_surface);
        }
        (back, cold) => assert_eq!(back, cold),
    }
    assert_eq!(
        back.accepted_state_identity(),
        cold.accepted_state_identity()
    );
}

#[test]
fn accepted_dielectric_refresh_retains_coated_anisotropic_scalar_and_polarized_transport() {
    let card = card(1.0);
    let hot = state(&card, 600.0);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let wo = Vec3::new(0.6, 0.0, 0.8);
    let wi = Vec3::new(0.0, 0.0, 1.0);
    let finish = DielectricSurface::try_anisotropic_rough(
        0.12, 0.45, Vec3::new(1.0, 0.0, 0.0),
    ).unwrap().with_coating(DielectricCoating::try_new(periodic(90e-9)).unwrap());
    let beam = PolarizedBeam::try_new(
        Point3::new(0.0, 0.0, 0.0),
        PolarizationFrame::try_new(wo.scale(-1.0), Vec3::new(0.0, 1.0, 0.0)).unwrap(),
        550.0,
        air(),
        Stokes::unpolarized(1.0),
    ).unwrap();

    for spectral in [false, true] {
        let (bare, source) = binding(&card, spectral, false, 100.0);
        let OpticalAppearance::Dielectric { glass, .. } = bare.appearance() else {
            panic!("native dielectric");
        };
        let cold = bare.with_surface_finish(
            OpticalAppearance::Dielectric { glass, surface: finish }, FINISH,
        ).unwrap();
        assert_ne!(cold.identity(), bare.identity(), "actual finish bytes enter identity");
        let (heated, refreshed) = source.refresh(cold, &card, &hot).unwrap();
        assert_eq!(heated.accepted_state_identity(), hot.identity());
        let back = refreshed
            .refresh(heated, &card, &state(&card, 300.0))
            .unwrap()
            .0;
        assert_pinned_round_trip(&source, cold, back, &card);
        let mut observations = Vec::new();
        for bound in [cold, heated] {
            let Material::Dielectric { glass, surface } = bound.material() else {
                panic!("native material() transport");
            };
            assert_eq!(surface, finish);
            let n = glass.ior().eval(550.0).unwrap();
            let interface = surface.interface(550.0, false).unwrap();
            let scalar = interface.fresnel(0.8, 1.0, n).unwrap();
            let polarized = interface.fresnel_polarized(0.8, 1.0, n).unwrap();
            near(polarized.reflection.apply(Stokes::unpolarized(1.0)).i,
                scalar.reflectance, 4e-14);
            let rough = interface.evaluate_rough(normal, wo, wi, 1.0, n, 0.12).unwrap();
            let rotated = surface.with_tangent_axis(Vec3::new(0.0, 1.0, 0.0)).unwrap();
            let rotated_rough = rotated.interface(550.0, false).unwrap()
                .evaluate_rough(normal, wo, wi, 1.0, n, 0.12).unwrap();
            assert!(rough.value > 0.0);
            assert!((rough.value - rotated_rough.value).abs() > 1e-6,
                "the authored tangent axis reaches actual native reflection");
            let sampled = interface.sample_rough(
                normal, wo, 1.0, n, 0.12, 0.0, 0.0, 0.0,
            ).unwrap().unwrap();
            let carried = beam.sample_interface(normal, glass, surface, [0.0; 3])
                .unwrap().unwrap();
            assert_eq!(carried.event, DielectricEvent::Reflection);
            assert_eq!(carried.event, sampled.event);
            near(carried.pdf, sampled.pdf, 3e-12);
            near(carried.beam.power().i, sampled.radiance_weight, 3e-12);
            assert!(carried.beam.power().degree() > 1e-5);
            assert!(surface.interface(900.0, false).is_err());
            observations.push((
                scalar.reflectance,
                rough.value,
                carried.beam.power().q / carried.beam.power().i,
                glass.absorption().extinction_per_m(550.0).unwrap(),
            ));
        }
        assert!((observations[0].0 - observations[1].0).abs() > 1e-5);
        assert!((observations[0].1 - observations[1].1).abs() > 1e-5);
        assert!((observations[0].2 - observations[1].2).abs() > 1e-5);
        near(observations[1].3 - observations[0].3, 0.6, 2e-14);
    }
}

#[test]
fn accepted_conductor_refresh_retains_each_native_coating_and_polarization() {
    let card = card(1.0);
    let hot = state(&card, 600.0);
    let transparent = CoherentDielectricCoating::try_new(
        CauchyIor::try_new(1.45, 0.001, 0.0).unwrap(), 100e-9, film_source(),
    ).unwrap();
    let absorbing = CoherentAbsorbingCoating::try_new(
        VISIBLE_COMPLEX_IOR_WAVELENGTHS_NM.map(|wavelength| {
            CoatingIorSample::try_new(wavelength, 1.7, 0.08).unwrap()
        }),
        80e-9,
        film_source(),
    ).unwrap();
    let periodic = periodic(70e-9);
    let finish = ConductorSurface::try_rough(0.27).unwrap();

    for spectral in [false, true] {
        let (bare, source) = binding(&card, spectral, true, 100.0);
        let OpticalAppearance::Conductor { optics: original, .. } = bare.appearance() else {
            panic!("native conductor");
        };
        for coated in [
            original.with_coating(transparent),
            original.with_absorbing_coating(absorbing),
            original.with_periodic_coating(periodic),
        ] {
            let cold = bare.with_surface_finish(
                OpticalAppearance::Conductor { optics: coated, surface: finish }, FINISH,
            ).unwrap();
            let (heated, refreshed) = source.refresh(cold, &card, &hot).unwrap();
            let back = refreshed
                .refresh(heated, &card, &state(&card, 300.0))
                .unwrap()
                .0;
            assert_pinned_round_trip(&source, cold, back, &card);
            let mut reflectances = Vec::new();
            for bound in [cold, heated] {
                let Material::Conductor { optics, surface } = bound.material() else {
                    panic!("native material() conductor");
                };
                assert_eq!(surface, finish);
                assert_eq!(optics.coating(), coated.coating());
                assert_eq!(optics.absorbing_coating(), coated.absorbing_coating());
                assert_eq!(optics.periodic_coating(), coated.periodic_coating());
                assert_eq!(optics.max_knot_span_nm(), original.max_knot_span_nm());
                assert_eq!(optics.source().status(), ConductorDataStatus::MaterialStateClaims);
                let scalar = optics.fresnel(550.0, 1.0, 0.65).unwrap();
                let mueller = optics.fresnel_mueller(550.0, 1.0, 0.65).unwrap();
                near(mueller.apply(Stokes::unpolarized(1.0)).i, scalar, 5e-14);
                let polarized = mueller.apply(Stokes::try_new(1.0, 0.2, 0.3, 0.1).unwrap());
                assert!(polarized.i > 0.0 && polarized.i <= 1.0);
                assert!(polarized.degree() <= 1.0 + 2e-14);
                assert!(optics.fresnel(900.0, 1.0, 0.65).is_err());
                reflectances.push(scalar);
            }
            assert!((reflectances[0] - reflectances[1]).abs() > 1e-5);
            let OpticalAppearance::Conductor { optics: next, .. } = heated.appearance() else {
                panic!("native conductor");
            };
            assert_ne!(next.samples(), original.samples());
            assert_ne!(next.source().identity(), original.source().identity());
            assert_ne!(cold.identity(), heated.identity());
        }
    }
}

#[test]
fn finish_replacement_refuses_bulk_family_schema_and_source_substitution_atomically() {
    let card = card(1.0);
    for conductor in [false, true] {
        let (visible, source) = binding(&card, false, conductor, 100.0);
        let (_, spectral_source) = binding(&card, true, conductor, 100.0);
        let (other_family, _) = binding(&card, false, !conductor, 100.0);
        let (hot, hot_source) = source.refresh(visible, &card, &state(&card, 600.0)).unwrap();
        let before = visible;
        for wrong in [other_family.appearance(), hot.appearance()] {
            assert_eq!(
                visible.with_surface_finish(wrong, FINISH),
                Err(MaterialOpticalError::BulkOpticsMismatch),
            );
        }
        assert_eq!(
            visible.with_surface_finish(visible.appearance(), ContentHash([0; 32])),
            Err(MaterialOpticalError::MissingSurfaceIdentity),
        );
        assert!(matches!(
            spectral_source.refresh(visible, &card, &state(&card, 600.0)),
            Err(MaterialOpticalRefreshError::SourceMismatch),
        ));
        assert!(matches!(
            hot_source.refresh(visible, &card, &state(&card, 600.0)),
            Err(MaterialOpticalRefreshError::SourceMismatch),
        ));
        let other_card = self::card(2.0);
        assert!(matches!(
            source.refresh(visible, &other_card, &state(&other_card, 600.0)),
            Err(MaterialOpticalRefreshError::SourceMismatch),
        ));
        if let OpticalAppearance::Conductor { optics, surface } = visible.appearance() {
            let replaced_source = ConductorOptics::try_new_in_band(
                optics.samples(), optics.max_knot_span_nm(), film_source(),
            ).unwrap();
            assert_eq!(
                visible.with_surface_finish(
                    OpticalAppearance::Conductor { optics: replaced_source, surface }, FINISH,
                ),
                Err(MaterialOpticalError::BulkOpticsMismatch),
            );
            let widened_interpolation = ConductorOptics::try_new_in_band(
                optics.samples(), optics.max_knot_span_nm() * 2.0, optics.source(),
            ).unwrap();
            assert_eq!(
                visible.with_surface_finish(
                    OpticalAppearance::Conductor { optics: widened_interpolation, surface }, FINISH,
                ),
                Err(MaterialOpticalError::BulkOpticsMismatch),
            );
        }
        let reidentified = visible.with_surface_finish(
            visible.appearance(), ContentHash([193; 32]),
        ).unwrap();
        assert_eq!(reidentified.material(), visible.material());
        assert_ne!(reidentified.identity(), visible.identity());
        assert_eq!(visible, before);
    }
}

#[test]
fn spectral_finish_refresh_keeps_independent_film_support_even_at_zero_thickness() {
    let card = card(20.0);
    let hot = state(&card, 600.0);
    let (bare, source) = binding(&card, true, true, 2000.0);
    let OpticalAppearance::Conductor { optics, surface } = bare.appearance() else {
        panic!("native infrared conductor");
    };
    for thickness in [0.0, 0.5e-6] {
        let film_knots = core::array::from_fn(|i| {
            CoatingIorSample::try_new_in_band(9000.0 + 500.0 * i as f64, 1.4, 0.08).unwrap()
        });
        let absorbing = CoherentAbsorbingCoating::try_new_in_band(
            film_knots, thickness, 500.0, film_source(),
        ).unwrap();
        let real_knots = core::array::from_fn(|i| {
            DielectricIorSample::try_new_in_band(9000.0 + 500.0 * i as f64, 1.4).unwrap()
        });
        let transparent = CoherentDielectricCoating::try_new(
            TabulatedIor::try_new_in_band(real_knots, 500.0).unwrap(),
            thickness,
            film_source(),
        ).unwrap();
        for coated in [
            optics.with_coating(transparent),
            optics.with_absorbing_coating(absorbing),
        ] {
            let cold = bare.with_surface_finish(
                OpticalAppearance::Conductor { optics: coated, surface }, FINISH,
            ).unwrap();
            let heated = source.refresh(cold, &card, &hot).unwrap().0;
            for bound in [cold, heated] {
                let Material::Conductor { optics, .. } = bound.material() else {
                    panic!("native conductor");
                };
                assert!(optics.fresnel(10_000.0, 1.0, 0.8).is_ok());
                for unsupported in [8000.0, 14_000.0] {
                    assert!(optics.complex_ior(unsupported).is_ok(), "bulk still supports query");
                    assert!(matches!(optics.fresnel(unsupported, 1.0, 0.8),
                        Err(ConductorError::InvalidWavelength)));
                    assert!(matches!(optics.fresnel_mueller(unsupported, 1.0, 0.0),
                        Err(ConductorError::InvalidWavelength)));
                }
            }
        }
    }

    // A compact dielectric coating retains its native visible support, even
    // when the source-backed bulk is infrared and the authored films vanish.
    let (bare, source) = binding(&card, true, false, 2000.0);
    let OpticalAppearance::Dielectric { glass, .. } = bare.appearance() else {
        panic!("native infrared dielectric");
    };
    let finish = DielectricSurface::try_anisotropic_rough(
        0.12, 0.45, Vec3::new(1.0, 0.0, 0.0),
    ).unwrap().with_coating(DielectricCoating::try_new(periodic(0.0)).unwrap());
    let cold = bare.with_surface_finish(
        OpticalAppearance::Dielectric { glass, surface: finish }, FINISH,
    ).unwrap();
    let heated = source.refresh(cold, &card, &hot).unwrap().0;
    for bound in [cold, heated] {
        let Material::Dielectric { glass, surface } = bound.material() else {
            panic!("native dielectric");
        };
        assert_eq!(surface, finish);
        assert!(glass.ior().eval(10_000.0).is_ok());
        assert!(matches!(surface.interface(10_000.0, false),
            Err(DielectricError::InvalidWavelength)));
    }
}
