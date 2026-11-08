//! Analytic source-law and original-density integration fixtures.

use super::*;
use crate::volumes::{
    AbsorptionImageSettings, SpectralAbsorptionImageSettings, VolumeGrid,
    render_density_spectral_box,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyValue,
    Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::*;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};

mod spectral;

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

pub(super) fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("thermal", temperature_quantity(), temperature)
        .unwrap()
}

fn curve(key: PropertyKey, cold: f64, hot: f64) -> PropertyClaim {
    let dims = key.quantity().dims();
    PropertyClaim {
        key,
        value: PropertyValue::Curve {
            abscissa: "thermal".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, cold), (400.0, hot)],
            dims,
        },
        validity: ValidityDomain::unconstrained().with_quantity(
            "thermal",
            temperature_quantity(),
            300.0,
            400.0,
        ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(),
        provenance: Provenance {
            source: "synthetic bulk absorption equation fixture, not measured stock".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn card(claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic bulk absorption material".into(),
            phase: "solid".into(),
            process: "explicit source law".into(),
            revision: 0,
        },
        set,
        Vec::new(),
    )
    .unwrap()
}

pub(super) fn density_key() -> PropertyKey {
    PropertyKey::new(DENSITY_PROPERTY, Density::DIMS)
}

pub(super) fn complex_fixture(
    k: [f64; 2],
) -> (
    MaterialCard,
    ResolvedMaterialStatePoint,
    SpectralComplexIndexStatePoint,
) {
    let mut claims = vec![curve(density_key(), 1000.0, 800.0)];
    let mut bindings = Vec::new();
    for (index, (wavelength_nm, loss)) in [8000.0, 12000.0].into_iter().zip(k).enumerate() {
        let n_key = PropertyKey::new(format!("original-n-{index}"), Dims::NONE);
        let k_key = PropertyKey::new(format!("original-k-{index}"), Dims::NONE);
        claims.push(curve(n_key.clone(), 1.5, 1.6));
        claims.push(curve(k_key.clone(), loss, 2.0 * loss));
        bindings.push(SpectralComplexIndexKnot {
            wavelength_nm,
            real_index_property: n_key,
            extinction_index_property: k_key,
        });
    }
    let requirements = claims
        .iter()
        .map(|claim| {
            ScalarPropertyRequirement::try_with_key(&claim.key, ScalarAdmissibility::Finite)
                .unwrap()
        })
        .collect::<Vec<_>>();
    let card = card(claims);
    let state = resolve_material_state_point(
        &card,
        &point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let optical =
        SpectralComplexIndexStatePoint::try_from_resolved(state.clone(), &bindings).unwrap();
    (card, state, optical)
}

fn visible_fixture(rgb: bool) -> (ResolvedMaterialStatePoint, VisibleOpticalStatePoint) {
    let dimless = QuantitySpec::dimensional(Dims::NONE);
    let mut rows = vec![
        (
            DENSITY_PROPERTY,
            QuantitySpec::dimensional(Density::DIMS),
            1000.0,
        ),
        (VISIBLE_DIELECTRIC_CAUCHY_A_PROPERTY, dimless, 1.5),
        (
            VISIBLE_DIELECTRIC_CAUCHY_B_M2_PROPERTY,
            QuantitySpec::dimensional(Dims([2, 0, 0, 0, 0, 0])),
            0.0,
        ),
        (
            VISIBLE_DIELECTRIC_CAUCHY_C_M4_PROPERTY,
            QuantitySpec::dimensional(Dims([4, 0, 0, 0, 0, 0])),
            0.0,
        ),
    ];
    if rgb {
        for name in VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES {
            rows.push((name, dimless, 0.9));
        }
        rows.push((
            VISIBLE_DIELECTRIC_REFERENCE_DISTANCE_M_PROPERTY,
            QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
            1.0,
        ));
    } else {
        for (index, wavelength) in VISIBLE_COMPLEX_IOR_WAVELENGTHS_NM.into_iter().enumerate() {
            rows.push((
                VISIBLE_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[index],
                VISIBLE_DIELECTRIC_ABSORPTION_QUANTITY,
                10.0 + 0.1 * (wavelength - 380.0),
            ));
        }
    }
    let card = card(
        rows.into_iter()
            .map(|(name, quantity, value)| {
                curve(PropertyKey::with_quantity(name, quantity), value, value)
            })
            .collect(),
    );
    let state = resolve_material_state_point(
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
    let optical = resolve_visible_optical_state_point(
        &card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (state, optical)
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8e-13 * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

fn settings() -> SpectralAbsorptionImageSettings {
    SpectralAbsorptionImageSettings {
        image: AbsorptionImageSettings {
            resolution: [1, 1],
            max_pixels: 1,
            max_grid_cells: 1,
            max_cell_visits: 1,
        },
        max_spectral_evaluations: 3,
    }
}

#[test]
fn g1_original_complex_index_law_drives_exact_samples_and_density_opacity() {
    let (_, accepted, optical) = complex_fixture([0.01, 0.03]);
    let binding =
        MaterialAbsorptionBinding::try_spectral(&accepted, &optical, &density_key()).unwrap();
    let wavelengths = [8000.0, 10000.0, 12000.0];
    let samples = binding.sample(&wavelengths).unwrap();
    assert_eq!(binding.source_identity(), optical.identity());
    assert_eq!(samples.source_identity(), optical.identity());
    assert_eq!(samples.accepted_state_identity(), accepted.identity());
    assert_eq!(samples.binding_identity(), binding.identity());
    assert_eq!(samples.source_wavelength_bounds_nm(), [8000.0, 12000.0]);
    assert!(core::ptr::eq(
        binding.density_property(),
        accepted.property_by_key(&density_key()).unwrap(),
    ));
    assert!(core::ptr::eq(
        samples.wavelengths_nm().as_ptr(),
        wavelengths.as_ptr()
    ));

    let data = [1000.0];
    let grid = VolumeGrid::new([1; 3], &data, [0.0; 3], [1.0, 1.0, 1e-5]);
    let image =
        render_density_spectral_box(&grid, settings(), samples.as_density_spectrum(1.0, 1.0))
            .unwrap();
    for (index, (wavelength, k)) in wavelengths.into_iter().zip([0.01, 0.02, 0.03]).enumerate() {
        let alpha = 4.0 * core::f64::consts::PI * k / (wavelength * 1e-9);
        close(samples.mass_absorption_m2_kg()[index], alpha / 1000.0);
        close(image.samples()[index].optical_depth(), alpha * 1e-5);
        close(
            image.samples()[index].transmittance(),
            (-alpha * 1e-5).exp(),
        );
    }
    let interpolated_coefficients =
        0.5 * (samples.mass_absorption_m2_kg()[0] + samples.mass_absorption_m2_kg()[2]);
    assert!(
        (samples.mass_absorption_m2_kg()[1] / interpolated_coefficients - 1.0).abs() > 0.06,
        "interpolating kappa would silently replace the original linear-k law"
    );

    let mut remapped = optical.bindings().to_vec();
    remapped[1].wavelength_nm = 14000.0;
    let other =
        SpectralComplexIndexStatePoint::try_from_resolved(accepted.clone(), &remapped).unwrap();
    assert_eq!(other.resolved().identity(), optical.resolved().identity());
    let other_binding =
        MaterialAbsorptionBinding::try_spectral(&accepted, &other, &density_key()).unwrap();
    assert_ne!(other_binding.source_identity(), binding.source_identity());
    assert_ne!(other_binding.identity(), binding.identity());
}

#[test]
fn g3_refresh_rebinds_original_temperature_density_and_absorption_together() {
    let (card, cold, cold_optics) = complex_fixture([0.01, 0.03]);
    let before =
        MaterialAbsorptionBinding::try_spectral(&cold, &cold_optics, &density_key()).unwrap();
    let hot = cold.at_query_point(&card, &point(400.0)).unwrap();
    let hot_optics = cold_optics.at_query_point(&card, &point(400.0)).unwrap();
    assert_eq!(
        MaterialAbsorptionBinding::try_spectral(&hot, &cold_optics, &density_key()),
        Err(MaterialAbsorptionError::Optical(
            MaterialOpticalError::StateMismatch
        ))
    );
    let after = MaterialAbsorptionBinding::try_spectral(&hot, &hot_optics, &density_key()).unwrap();
    close(
        after.mass_absorption_at(10000.0).unwrap() / before.mass_absorption_at(10000.0).unwrap(),
        2.5,
    );
    assert_eq!(before.material_density_kg_m3(), 1000.0);
    assert_eq!(after.material_density_kg_m3(), 800.0);
    assert_ne!(before.identity(), after.identity());
    for (binding, k) in [(before, 0.02), (after, 0.04)] {
        let samples = binding.sample(&[10000.0]).unwrap();
        let density = [binding.material_density_kg_m3()];
        let grid = VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0, 1.0, 1e-5]);
        let image =
            render_density_spectral_box(&grid, settings(), samples.as_density_spectrum(1.0, 1.0))
                .unwrap();
        let tau = 4.0 * core::f64::consts::PI * k;
        close(image.samples()[0].optical_depth(), tau);
        close(image.samples()[0].transmittance(), (-tau).exp());
    }
    for (original, refreshed) in cold.properties().iter().zip(hot.properties()) {
        assert_eq!(
            original.answer().receipt.selected,
            refreshed.answer().receipt.selected
        );
    }
    assert_eq!(before.accepted_state().identity(), cold.identity());
    assert_eq!(
        before.optical_state().identity(),
        cold_optics.resolved().identity()
    );
}

#[test]
fn g1_visible_dielectric_uses_native_alpha_and_refuses_rgb_or_missing_density() {
    let (accepted, optical) = visible_fixture(false);
    let binding =
        MaterialAbsorptionBinding::try_visible(&accepted, &optical, &density_key()).unwrap();
    let BulkOpticalLaw::Dielectric(glass) = BulkOpticalLaw::from_state(&optical).unwrap() else {
        panic!("native dielectric");
    };
    let native = glass.absorption().extinction_per_m(555.0).unwrap();
    assert_eq!(binding.mass_absorption_at(555.0).unwrap(), native / 1000.0);
    close(native, 27.5);
    assert_eq!(binding.wavelength_bounds_nm(), [380.0, 780.0]);
    assert_eq!(
        MaterialAbsorptionBinding::try_visible(
            &accepted,
            &optical,
            &PropertyKey::new("absent density", Density::DIMS),
        ),
        Err(MaterialAbsorptionError::Density)
    );
    let (rgb_state, rgb_optics) = visible_fixture(true);
    assert_eq!(
        MaterialAbsorptionBinding::try_visible(&rgb_state, &rgb_optics, &density_key()),
        Err(MaterialAbsorptionError::SpectralAbsorptionRequired)
    );
}

#[test]
fn g0_closed_source_support_and_bounded_samples_survive_lossless_materials() {
    let (_, accepted, optical) = complex_fixture([0.0, 0.0]);
    let binding =
        MaterialAbsorptionBinding::try_spectral(&accepted, &optical, &density_key()).unwrap();
    assert_eq!(
        binding
            .sample(&[8000.0, 12000.0])
            .unwrap()
            .mass_absorption_m2_kg(),
        &[0.0, 0.0]
    );
    for wavelength in [550.0, 7999.0, 12001.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            binding.mass_absorption_at(wavelength),
            Err(MaterialAbsorptionError::WavelengthOutsideSource { .. })
        ));
    }
    for wavelengths in [&[][..], &[9000.0, 9000.0], &[10000.0, 9000.0], &[f64::NAN]] {
        assert_eq!(
            binding.sample(wavelengths),
            Err(MaterialAbsorptionError::InvalidWavelengthSamples)
        );
    }
    let excessive = core::array::from_fn::<_, 17, _>(|i| 8000.0 + i as f64);
    assert_eq!(
        binding.sample(&excessive),
        Err(MaterialAbsorptionError::InvalidWavelengthSamples)
    );
}

#[test]
fn g1_complex_conversion_avoids_intermediate_overflow_and_refuses_lost_absorption() {
    let scale = 4.0 * core::f64::consts::PI * 1e9;
    close(
        complex_mass_absorption(1e308, 500.0, 1e308).unwrap(),
        scale / 500.0,
    );
    close(
        complex_mass_absorption(1e-300, 1e-300, 1e300).unwrap(),
        scale * 1e-300,
    );
    assert_eq!(
        complex_mass_absorption(f64::from_bits(1), 10000.0, f64::MAX),
        Err(MaterialAbsorptionError::UnrepresentableCoefficient)
    );
    assert_eq!(
        checked_coefficient(f64::from_bits(1) / f64::MAX, false),
        Err(MaterialAbsorptionError::UnrepresentableCoefficient)
    );
}
