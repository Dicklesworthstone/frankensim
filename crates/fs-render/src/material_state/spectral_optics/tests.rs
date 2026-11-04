use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECTRAL_CONDUCTOR_ETA_PROPERTIES,
    SPECTRAL_CONDUCTOR_K_PROPERTIES, SPECTRAL_CONDUCTOR_K_QUANTITY,
    SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES, SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
    SPECTRAL_DIELECTRIC_IOR_PROPERTIES, SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
    resolve_spectral_conductor_state_point, resolve_spectral_tabulated_dielectric_state_point,
};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};

const WAVELENGTHS_UM: [f64; 9] = [8.0, 8.25, 9.0, 9.5, 10.0, 10.25, 11.0, 12.0, 14.0];
const MAX_GAP_NM: f64 = 2001.0;
const SURFACE: ContentHash = ContentHash([173; 32]);

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new().with("T", temperature).unwrap()
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
            source: format!("synthetic infrared binding source: {name}"),
            license: "internal-test-use".into(),
            artifact: None,
        },
    }
}

fn card(wavelength_scale: f64, hot_last_wavelength_m: Option<f64>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (index, wavelength_um) in WAVELENGTHS_UM.into_iter().enumerate() {
        let wavelength = wavelength_um * 1.0e-6 * wavelength_scale;
        let wavelength_hot = if index == 8 {
            hot_last_wavelength_m.unwrap_or(wavelength)
        } else {
            wavelength
        };
        let eta = 1.25 + (wavelength_um - 8.0) / 8.0;
        let k = 2.0 + (wavelength_um - 8.0) / 4.0;
        let n = 1.75 - wavelength_um / 64.0;
        let alpha = if index == 4 { 6.25 } else { 0.25 };
        for source in [
            claim(
                SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[index],
                QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                wavelength,
                wavelength_hot,
            ),
            claim(
                SPECTRAL_CONDUCTOR_ETA_PROPERTIES[index],
                QuantitySpec::dimensional(Dims::NONE),
                eta,
                eta + 0.25,
            ),
            claim(
                SPECTRAL_CONDUCTOR_K_PROPERTIES[index],
                SPECTRAL_CONDUCTOR_K_QUANTITY,
                k,
                k + 0.5,
            ),
            claim(
                SPECTRAL_DIELECTRIC_IOR_PROPERTIES[index],
                QuantitySpec::dimensional(Dims::NONE),
                n,
                n + 0.0625,
            ),
            claim(
                SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[index],
                SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
                alpha,
                alpha + 0.5,
            ),
        ] {
            claims.insert_claim(source).unwrap();
        }
    }
    claims
        .insert_claim(claim(
            DENSITY_PROPERTY,
            QuantitySpec::dimensional(Density::DIMS),
            2500.0,
            2400.0,
        ))
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic IR material".into(),
            phase: "solid".into(),
            process: "explicit source spectrum".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn state(card: &MaterialCard, temperature: f64) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(
        card,
        &point(temperature),
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

#[test]
fn g1_ir_conductor_binding_and_refresh_reach_native_complex_fresnel() {
    let card = card(1.0, None);
    let cold = state(&card, 300.0);
    let optical = resolve_spectral_conductor_state_point(
        &card,
        &cold.query(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let binding =
        MaterialOpticalBinding::try_spectral_conductor(&cold, &optical, MAX_GAP_NM, 0.12, SURFACE)
            .unwrap();
    let original = SpectralOpticalStatePoint::Conductor(optical.clone());
    let hot = cold.at_query_point(&card, &point(600.0)).unwrap();
    let (heated, refreshed) = binding
        .at_spectral_material_state(&original, &card, &hot)
        .unwrap();
    let bounds = (
        optical.samples()[0].wavelength_nm,
        optical.samples()[8].wavelength_nm,
    );
    for (bound, source, delta) in [(binding, &original, 0.0), (heated, &refreshed, 1.0)] {
        let Material::Conductor { optics, surface } = bound.material() else {
            panic!("native conductor");
        };
        assert_eq!(optics.wavelength_bounds_nm(), bounds);
        assert_eq!(optics.max_knot_span_nm(), MAX_GAP_NM);
        assert_eq!(optics.source().identity(), source.identity());
        assert_eq!(
            optics.source().status(),
            ConductorDataStatus::MaterialStateClaims
        );
        assert_eq!(surface.roughness_alpha(), 0.12);
        for wavelength_nm in [bounds.0, 10_125.0, bounds.1] {
            let offset_um = wavelength_nm * 1e-3 - 8.0;
            let eta = 1.25 + offset_um / 8.0 + 0.25 * delta;
            let k = 2.0 + offset_um / 4.0 + 0.5 * delta;
            let actual = optics.complex_ior(wavelength_nm).unwrap();
            assert!((actual.eta() - eta).abs() < 2e-14);
            assert!((actual.k() - k).abs() < 2e-14);
            let expected_r = ((eta - 1.0).powi(2) + k * k) / ((eta + 1.0).powi(2) + k * k);
            assert!((optics.fresnel(wavelength_nm, 1.0, 1.0).unwrap() - expected_r).abs() < 2e-14);
        }
        for wavelength_nm in [550.0, bounds.0 - 1.0, bounds.1 + 1.0] {
            assert!(optics.complex_ior(wavelength_nm).is_err());
        }
    }
    assert_eq!(heated.accepted_state_identity(), hot.identity());
    assert_ne!(binding.identity(), heated.identity());
    assert!(matches!(
        heated.at_spectral_material_state(&original, &card, &hot),
        Err(MaterialOpticalRefreshError::SourceMismatch)
    ));
    let mirror =
        MaterialOpticalBinding::try_spectral_conductor(&cold, &optical, MAX_GAP_NM, 0.0, SURFACE)
            .unwrap();
    let Material::Conductor { surface, .. } = mirror
        .at_spectral_material_state(&original, &card, &hot)
        .unwrap()
        .0
        .material()
    else {
        panic!("native mirror");
    };
    assert!(surface.is_delta());
}

#[test]
fn g1_ir_dielectric_keeps_phase_group_snell_and_original_absorption_feature() {
    let card = card(1.0, None);
    let cold = state(&card, 300.0);
    let optical = resolve_spectral_tabulated_dielectric_state_point(
        &card,
        &cold.query(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let original = SpectralOpticalStatePoint::TabulatedDielectric(optical.clone());
    let hot = cold.at_query_point(&card, &point(600.0)).unwrap();
    for roughness in [None, Some(0.15)] {
        let binding = MaterialOpticalBinding::try_spectral_tabulated_dielectric(
            &cold, &optical, MAX_GAP_NM, roughness, SURFACE,
        )
        .unwrap();
        let (heated, _) = binding
            .at_spectral_material_state(&original, &card, &hot)
            .unwrap();
        for (bound, delta) in [(binding, 0.0), (heated, 1.0)] {
            let Material::Dielectric { glass, surface } = bound.material() else {
                panic!("native dielectric");
            };
            assert_eq!(surface.roughness_alpha(), roughness);
            let wavelength_nm = 10_125.0;
            let n = 1.75 - wavelength_nm / 64_000.0 + delta * 0.0625;
            let alpha = 3.25 + delta * 0.5;
            assert!((glass.ior().eval(wavelength_nm).unwrap() - n).abs() < 2e-14);
            assert!(
                (glass.ior().group_index(wavelength_nm).unwrap() - (1.75 + delta * 0.0625)).abs()
                    < 2e-12
            );
            assert!(
                (glass.absorption().extinction_per_m(wavelength_nm).unwrap() - alpha).abs() < 2e-12
            );
            let actual_t = glass
                .absorption()
                .transmittance(wavelength_nm, 0.4)
                .unwrap();
            assert!((actual_t - (-alpha * 0.4_f64).exp()).abs() < 3e-13);
            for (index, expected) in [(4, 6.25), (5, 0.25)] {
                let wavelength_nm = optical.samples()[index].wavelength_m * 1e9;
                assert_eq!(
                    glass.absorption().extinction_per_m(wavelength_nm).unwrap(),
                    expected + delta * 0.5
                );
            }
            let interface = surface.interface(wavelength_nm, false).unwrap();
            let normal = interface.fresnel(1.0, 1.0, n).unwrap();
            assert!((normal.reflectance - ((n - 1.0) / (n + 1.0)).powi(2)).abs() < 2e-14);
            let oblique = interface.fresnel(0.8, 1.0, n).unwrap();
            assert!((oblique.transmitted_cosine - (1.0 - (0.6 / n).powi(2)).sqrt()).abs() < 2e-14);
            let bounds = glass.wavelength_bounds_nm().unwrap();
            for wavelength_nm in [550.0, bounds.0 - 1.0, bounds.1 + 1.0] {
                assert!(glass.ior().eval(wavelength_nm).is_err());
                assert!(glass.absorption().extinction_per_m(wavelength_nm).is_err());
            }
        }
    }
}

#[test]
fn g0_binding_refuses_state_surface_interpolation_and_si_conversion_mismatches() {
    let card = card(1.0, None);
    let cold = state(&card, 300.0);
    let optical = resolve_spectral_conductor_state_point(
        &card,
        &cold.query(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(
        MaterialOpticalBinding::try_spectral_conductor(
            &state(&card, 600.0),
            &optical,
            MAX_GAP_NM,
            0.0,
            SURFACE
        ),
        Err(MaterialOpticalError::StateMismatch)
    );
    assert_eq!(
        MaterialOpticalBinding::try_spectral_conductor(
            &cold,
            &optical,
            MAX_GAP_NM,
            0.0,
            ContentHash([0; 32])
        ),
        Err(MaterialOpticalError::MissingSurfaceIdentity)
    );
    for gap in [0.0, -1.0, f64::INFINITY, f64::NAN, 1000.0] {
        assert!(
            MaterialOpticalBinding::try_spectral_conductor(&cold, &optical, gap, 0.0, SURFACE)
                .is_err()
        );
    }
    let foreign = self::card(2.0, None);
    assert_eq!(
        MaterialOpticalBinding::try_spectral_conductor(
            &state(&foreign, 300.0),
            &optical,
            MAX_GAP_NM,
            0.0,
            SURFACE
        ),
        Err(MaterialOpticalError::StateMismatch)
    );
    let enormous = self::card(1e306, None);
    let enormous_state = state(&enormous, 300.0);
    assert!(
        resolve_spectral_conductor_state_point(
            &enormous,
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .is_err(),
        "SI-to-nanometre overflow refuses at the shared complex-index owner"
    );
    let source = resolve_spectral_tabulated_dielectric_state_point(
        &enormous,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(
        MaterialOpticalBinding::try_spectral_tabulated_dielectric(
            &enormous_state,
            &source,
            f64::MAX,
            None,
            SURFACE
        )
        .is_err()
    );
}

#[test]
fn g0_refresh_keeps_original_interpolation_limit_and_is_atomic_on_a_wider_gap() {
    let card = card(1.0, Some(16e-6));
    let cold = state(&card, 300.0);
    let optical = resolve_spectral_conductor_state_point(
        &card,
        &cold.query(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let original = SpectralOpticalStatePoint::Conductor(optical.clone());
    let binding =
        MaterialOpticalBinding::try_spectral_conductor(&cold, &optical, MAX_GAP_NM, 0.0, SURFACE)
            .unwrap();
    let before = binding;
    let hot = cold.at_query_point(&card, &point(600.0)).unwrap();
    assert!(
        original.at_query_point(&card, &hot.query()).is_ok(),
        "the source is supported at the endpoint"
    );
    assert!(matches!(
        binding.at_spectral_material_state(&original, &card, &hot),
        Err(MaterialOpticalRefreshError::Binding(
            MaterialOpticalError::Conductor(ConductorError::InvalidTable)
        ))
    ));
    assert_eq!(binding, before);
    assert_eq!(binding.accepted_state_identity(), cold.identity());
}

#[test]
fn g0_two_families_adopted_from_one_complete_snapshot_cannot_substitute_on_refresh() {
    let card = card(1.0, None);
    let requirements: Vec<_> = card
        .claims()
        .claims_ordered()
        .map(|(_, source)| {
            ScalarPropertyRequirement::try_with_key(&source.key, ScalarAdmissibility::Finite)
                .unwrap()
        })
        .collect();
    let accepted = resolve_material_state_point(
        &card,
        &point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let conductor = SpectralConductorStatePoint::try_from_material(accepted.clone()).unwrap();
    let dielectric =
        SpectralTabulatedDielectricStatePoint::try_from_material(accepted.clone()).unwrap();
    assert_eq!(
        conductor.resolved().identity(),
        dielectric.resolved().identity()
    );
    let binding = MaterialOpticalBinding::try_spectral_conductor(
        &accepted, &conductor, MAX_GAP_NM, 0.0, SURFACE,
    )
    .unwrap();
    assert!(matches!(
        binding.at_spectral_material_state(
            &SpectralOpticalStatePoint::TabulatedDielectric(dielectric),
            &card,
            &accepted
        ),
        Err(MaterialOpticalRefreshError::SourceMismatch)
    ));
}
