use super::{MaterialThermoOpticBinding, MaterialThermoOpticBindingError};
use crate::dielectric::{
    DielectricEvent, DielectricGlass, DielectricIor, DielectricSurface, SellmeierIor,
    sample_smooth_dielectric,
};
use crate::material_state::{MaterialOpticalBinding, OpticalAppearance};
use crate::tracer::Material;
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_geom::Vec3;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ResolvedMaterialStatePoint,
    SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES, SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
    SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES,
    SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES, ScalarAdmissibility, ScalarPropertyRequirement,
    VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES, VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES,
    resolve_material_state_point, resolve_spectral_sellmeier_state_point,
};
use fs_material::thermo_optic::{SchottThermoOpticLaw, ThermoOpticCoefficients};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec};

const SURFACE: ContentHash = ContentHash([91; 32]);

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

// Finite-band optics are admitted only at fixed external conditions, and an axis is one only if it carries its semantic
// type (state_point.rs:416-444, spectral_optics.rs:482): a pressure typed by dimensions alone is refused. As in
// particle_scattering/tests.rs:38.
fn pressure() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}

fn point(temperature: f64, pressure_pa: f64, typed: bool) -> QueryPoint {
    let query = QueryPoint::new()
        .with_quantity("P", pressure(), pressure_pa)
        .unwrap();
    if typed {
        query
            .with_quantity("T", absolute_temperature(), temperature)
            .unwrap()
    } else {
        query.with("T", temperature).unwrap()
    }
}

fn card(typed: bool, foreign: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut add = |name: &str, quantity: QuantitySpec, value: f64| {
        let domain =
            ValidityDomain::unconstrained().with_quantity("P", pressure(), 90_000.0, 120_000.0);
        let validity = if typed {
            domain.with_quantity("T", absolute_temperature(), 250.0, 450.0)
        } else {
            domain.with("T", 250.0, 450.0)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Scalar {
                    value,
                    dims: quantity.dims(),
                },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: Vec::new(),
                provenance: Provenance {
                    source: format!("manufactured thermal binding verification: {name}"),
                    license: "internal-test-use".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    add(
        DENSITY_PROPERTY,
        QuantitySpec::dimensional(Density::DIMS),
        2500.0,
    );
    // An exactly constant n0 = 1.5 makes the thermal wavelength derivative
    // independent of the implementation's source-dispersion differentiation.
    for i in 0..3 {
        add(
            VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[i],
            QuantitySpec::dimensional(Dims::NONE),
            [1.25, 0.0, 0.0][i],
        );
        add(
            VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES[i],
            QuantitySpec::dimensional(Dims([2, 0, 0, 0, 0, 0])),
            0.0,
        );
    }
    let length = QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0]));
    for (name, value) in SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES
        .into_iter()
        .zip([350e-9, 1700e-9])
    {
        add(name, length, value);
    }
    for i in 0..9 {
        add(
            SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[i],
            length,
            (400.0 + 150.0 * i as f64) * 1e-9,
        );
        add(
            SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i],
            SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
            0.15,
        );
    }
    let genesis = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured thermo-optic glass".into(),
            phase: "solid".into(),
            process: "analytic verification, not a measured specimen".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    if foreign {
        // A later revision of the same state with the same claims. assemble starts a lineage at revision 0
        // (fs-matdb cards.rs:599), so the successor comes from supersede, as in fs-material's refresh tests
        // (state_point.rs:3817).
        MaterialCard::supersede(&genesis, genesis.claims().clone(), Vec::new()).unwrap()
    } else {
        genesis
    }
}

fn state(
    card: &MaterialCard,
    temperature: f64,
    pressure: f64,
    typed: bool,
) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(
        card,
        &point(temperature, pressure, typed),
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

fn reference(
    card: &MaterialCard,
    accepted: &ResolvedMaterialStatePoint,
    roughness: Option<f64>,
    surface: ContentHash,
) -> MaterialOpticalBinding {
    let optical = resolve_spectral_sellmeier_state_point(
        card,
        &accepted.query(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    MaterialOpticalBinding::try_spectral_sellmeier(accepted, &optical, roughness, surface).unwrap()
}

fn law() -> SchottThermoOpticLaw {
    SchottThermoOpticLaw::try_new(
        ThermoOpticCoefficients {
            d: [2e-5, 3e-8, 0.0],
            e: [4e-7, 0.0],
            resonance_um: 0.2,
        },
        300.0,
        [270.0, 420.0],
        [400.0, 1600.0],
        ContentHash([92; 32]),
    )
    .unwrap()
}

fn glass(binding: &MaterialOpticalBinding) -> (DielectricGlass, DielectricSurface) {
    let OpticalAppearance::Dielectric { glass, surface } = binding.appearance() else {
        panic!("dielectric fixture required");
    };
    (glass, surface)
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_accepted_temperature_drives_native_phase_group_and_snell() {
    let card = card(true, false);
    let cold = state(&card, 300.0, 100_000.0, true);
    let original = reference(&card, &cold, None, SURFACE);
    let source = MaterialThermoOpticBinding::try_new(original, &cold, law(), "T").unwrap();
    let hot = cold
        .at_query_point(&card, &point(400.0, 100_000.0, true))
        .unwrap();
    let heated = source.at_material_state(&hot).unwrap();
    assert_eq!(heated.binding().accepted_state_identity(), hot.identity());
    assert_ne!(heated.binding().identity(), source.binding().identity());
    assert_eq!(heated.reference_binding(), &original);
    assert_eq!(heated.reference_state(), &cold);
    assert_eq!(heated.law(), &law());
    assert_eq!(heated.temperature_axis(), "T");
    let Material::Dielectric { glass, .. } = heated.binding().material() else {
        panic!("native tracer dielectric required");
    };
    for wavelength in [550.0, 1000.0, 1500.0] {
        // Independent TIE-19 phase increment and analytic -lambda dn/dlambda.
        let l: f64 = wavelength * 0.001;
        let q = l * l - 0.04;
        let a = (1.5 * 1.5 - 1.0) / (2.0 * 1.5);
        let phase = 1.5 + a * 100.0 * (2e-5 + 3e-8 * 100.0 + 4e-7 / q);
        let group = phase + 2.0 * a * 4e-7 * 100.0 * l * l / (q * q);
        let native = glass.ior().eval(wavelength).unwrap();
        near(native, phase, 4e-15);
        near(glass.ior().group_index(wavelength).unwrap(), group, 4e-15);
        let transmitted = sample_smooth_dielectric(
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.6, 0.0, 0.8),
            1.0,
            native,
            0.99,
        )
        .unwrap();
        assert_eq!(transmitted.event, DielectricEvent::Transmission);
        near(transmitted.direction.x, -0.6 / phase, 4e-15);
        near(
            transmitted.direction.z,
            -(1.0 - (0.6 / phase).powi(2)).sqrt(),
            4e-15,
        );
    }
}

#[test]
fn g3_temperature_cycles_keep_original_absorption_surface_and_reference() {
    let card = card(true, false);
    let cold = state(&card, 300.0, 100_000.0, true);
    for roughness in [None, Some(0.12)] {
        let original = reference(&card, &cold, roughness, SURFACE);
        let (original_glass, original_surface) = glass(&original);
        let source = MaterialThermoOpticBinding::try_new(original, &cold, law(), "T").unwrap();
        let mut observation = source.clone();
        for temperature in [420.0, 270.0, 350.0, 300.0] {
            let accepted = cold
                .at_query_point(&card, &point(temperature, 100_000.0, true))
                .unwrap();
            observation = observation.at_material_state(&accepted).unwrap();
            let direct = source.at_material_state(&accepted).unwrap();
            assert_eq!(observation, direct);
            assert_eq!(observation.reference_binding(), &original);
            let (current, surface) = glass(observation.binding());
            assert_eq!(surface, original_surface);
            assert_eq!(current.absorption_ref(), original_glass.absorption_ref());
            near(
                current
                    .absorption_ref()
                    .transmittance(1000.0, 0.37)
                    .unwrap(),
                (-0.15_f64 * 0.37).exp(),
                3e-15,
            );
        }
        // A refresh re-queries with PinnedByProperty (state_point.rs:1085-1093, 5c6b323c6e) and the accepted
        // state's identity hashes the selection policy, so the state reached by returning to T0 is the pinned
        // twin of `cold` (SingleClaimOnly): its identity differs from cold's by design and no physical value
        // does. Returning to T0 restores the source observed at that same pinned state, with the original
        // appearance, surface and reference.
        let restored = cold
            .at_query_point(&card, &point(300.0, 100_000.0, true))
            .unwrap();
        assert_eq!(
            observation,
            source.at_material_state(&restored).unwrap(),
            "returning to T0 restores the binding of the same pinned state exactly"
        );
        assert_eq!(
            observation.binding().appearance(),
            source.binding().appearance()
        );
        assert_eq!(
            observation.binding().accepted_state_identity(),
            restored.identity()
        );
        let other_surface = reference(&card, &cold, roughness, ContentHash([93; 32]));
        let other = MaterialThermoOpticBinding::try_new(other_surface, &cold, law(), "T").unwrap();
        assert_eq!(other.binding().appearance(), source.binding().appearance());
        assert_ne!(other.binding().identity(), source.binding().identity());
    }
}

#[test]
fn g0_binding_refuses_foreign_conditions_untyped_temperature_and_outside_law() {
    use MaterialThermoOpticBindingError as Error;

    let card = card(true, false);
    let cold = state(&card, 300.0, 100_000.0, true);
    let original = reference(&card, &cold, None, SURFACE);
    let source = MaterialThermoOpticBinding::try_new(original, &cold, law(), "T").unwrap();
    let warmer = state(&card, 310.0, 100_000.0, true);
    assert_eq!(
        MaterialThermoOpticBinding::try_new(original, &warmer, law(), "T"),
        Err(Error::ReferenceStateMismatch)
    );
    assert_eq!(
        MaterialThermoOpticBinding::try_new(
            reference(&card, &warmer, None, SURFACE),
            &warmer,
            law(),
            "T"
        ),
        Err(Error::ReferenceTemperature)
    );
    assert_eq!(
        MaterialThermoOpticBinding::try_new(original, &cold, law(), "reservoir_temperature"),
        Err(Error::TemperatureAxis)
    );
    let compressed = state(&card, 350.0, 110_000.0, true);
    assert_eq!(
        source.at_material_state(&compressed),
        Err(Error::StateMismatch)
    );
    let foreign = self::card(true, true);
    let foreign_state = state(&foreign, 350.0, 100_000.0, true);
    assert_eq!(
        source.at_material_state(&foreign_state),
        Err(Error::StateMismatch)
    );
    // The physical card accepts this state; the narrower optical law does not.
    let outside = state(&card, 430.0, 100_000.0, true);
    assert!(matches!(
        source.at_material_state(&outside),
        Err(Error::Optical(_))
    ));
    assert!(matches!(
        MaterialThermoOpticBinding::try_new(*source.binding(), &cold, law(), "T"),
        Err(Error::Optical(_))
    ));
    let untyped = self::card(false, false);
    let untyped_state = state(&untyped, 300.0, 100_000.0, false);
    assert_eq!(
        MaterialThermoOpticBinding::try_new(
            reference(&untyped, &untyped_state, None, SURFACE),
            &untyped_state,
            law(),
            "T"
        ),
        Err(Error::TemperatureAxis)
    );
    assert_eq!(source.reference_binding(), &original);
    assert_eq!(source.binding().accepted_state_identity(), cold.identity());
    near(
        glass(source.binding()).0.ior().eval(1000.0).unwrap(),
        1.5,
        3e-15,
    );
}

#[test]
fn g1_thermal_binding_feeds_existing_slab_emission_at_accepted_temperature() {
    let card = card(true, false);
    let cold = state(&card, 300.0, 100_000.0, true);
    let source = MaterialThermoOpticBinding::try_new(
        reference(&card, &cold, None, SURFACE),
        &cold,
        law(),
        "T",
    )
    .unwrap();
    let hot = state(&card, 400.0, 100_000.0, true);
    let heated = source.at_material_state(&hot).unwrap();
    let vacuum: DielectricIor = SellmeierIor::try_new_in_band([0.0; 3], [0.0; 3], (350.0, 1700.0))
        .unwrap()
        .into();
    let before = source
        .binding()
        .incoherent_slab_emission(&cold, "T", 0.04, vacuum)
        .unwrap();
    let after = heated
        .binding()
        .incoherent_slab_emission(&hot, "T", 0.04, vacuum)
        .unwrap();
    assert_eq!(before.film_temperature_k(0), Some(300.0));
    assert_eq!(after.film_temperature_k(0), Some(400.0));
    assert_ne!(before.optical_identity(), after.optical_identity());
    assert!(
        heated
            .binding()
            .incoherent_slab_emission(&cold, "T", 0.04, vacuum)
            .is_err()
    );
    let native = glass(heated.binding())
        .0
        .incoherent_slab_thermal_response_between(1500.0, 0.04, 0.7, 1.0, 1.0, 400.0)
        .unwrap();
    let emitted = after.evaluate(1500.0, 0.7).unwrap();
    for pol in 0..2 {
        let expected = native.polarized_radiance_w_m2_sr_nm[pol];
        assert!(expected.is_finite() && expected > 0.0);
        near(
            emitted.polarized_radiance_w_m2_sr_nm()[pol] / expected,
            1.0,
            2e-13,
        );
    }
    assert!(
        after.radiance_w_m2_sr_nm(1500.0, 0.7).unwrap()
            > before.radiance_w_m2_sr_nm(1500.0, 0.7).unwrap()
    );
}
