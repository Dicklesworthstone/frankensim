//! G0/G1 source-to-film coupling with synthetic constitutive inputs.
//! The independent dimensional dispersion check establishes numerical wiring,
//! not empirical film accuracy or source calibration.

use fs_evidence::ValidityDomain;
use fs_lbm::film::material::MaterialFilmLiquid;
use fs_lbm::film::{FilmError, LongWaveLimits};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::surface::{
    LiquidCapillarity, ResolvedLiquidSurfaceState, SURFACE_TENSION_DIMS, SURFACE_TENSION_PROPERTY,
    resolve_liquid_surface_state,
};
use fs_material::liquid::{
    DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS,
    KINEMATIC_VISCOSITY_PROPERTY,
};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Frequency, QuantitySpec, Temperature};

const LIMITS: LongWaveLimits = LongWaveLimits {
    max_kh: 0.1,
    max_re_kh: 0.1,
};

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap()
        .with("surface_environment", 7.0)
        .unwrap()
}

fn card(basis: LiquidViscosityBasis, with_surface: bool, shear_dependent: bool) -> MaterialCard {
    let (viscosity, dims, values) = match basis {
        LiquidViscosityBasis::Dynamic => (
            DYNAMIC_VISCOSITY_PROPERTY,
            DYNAMIC_VISCOSITY_DIMS,
            [0.002, 0.004],
        ),
        LiquidViscosityBasis::Kinematic => (
            KINEMATIC_VISCOSITY_PROPERTY,
            KINEMATIC_VISCOSITY_DIMS,
            [2e-6, 5e-6],
        ),
    };
    let mut rows = vec![
        (DENSITY_PROPERTY, Density::DIMS, [1000.0, 800.0]),
        (viscosity, dims, values),
    ];
    if with_surface {
        rows.push((SURFACE_TENSION_PROPERTY, SURFACE_TENSION_DIMS, [0.08, 0.04]));
    }
    let mut claims = ClaimSet::new();
    for (name, dims, values) in rows {
        let shear = shear_dependent && name == viscosity;
        let mut validity = ValidityDomain::unconstrained()
            .with_quantity(
                "T",
                *point(300.0).axis_quantities().get("T").unwrap(),
                300.0,
                400.0,
            )
            .with_quantity(
                "pressure",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::Pressure,
                    ValueForm::Static,
                )),
                100_000.0,
                100_000.0,
            )
            .with("surface_environment", 7.0, 7.0);
        if shear {
            validity = validity.with_quantity(
                "rate",
                QuantitySpec::dimensional(Frequency::DIMS),
                0.0,
                10.0,
            );
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: if shear { "rate" } else { "T" }.into(),
                    abscissa_dims: if shear {
                        Frequency::DIMS
                    } else {
                        Temperature::DIMS
                    },
                    knots: if shear {
                        vec![(0.0, values[0]), (10.0, values[1])]
                    } else {
                        vec![(300.0, values[0]), (400.0, values[1])]
                    },
                    dims,
                },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic material-film coupling fixture, not measurements".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic film liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn resolve(
    card: &MaterialCard,
    basis: LiquidViscosityBasis,
    capillarity: LiquidCapillarity,
) -> ResolvedLiquidSurfaceState {
    resolve_liquid_surface_state(
        card,
        &point(300.0),
        basis,
        capillarity,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4e-12 * expected.abs().max(1e-30),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn both_original_viscosity_bases_drive_dimensional_film_and_keep_receipts() {
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let card = card(basis, true, false);
        let source = resolve(&card, basis, LiquidCapillarity::Sourced);
        let bound = MaterialFilmLiquid::try_new(&source).unwrap();
        assert!(core::ptr::eq(bound.material(), &source));
        assert_eq!(
            bound.material().material().card_identity(),
            card.content_hash()
        );
        assert_eq!(bound.material().viscosity_basis(), basis);
        assert_eq!(bound.material().material().properties().len(), 3);
        let film = bound.film(2e-6, 4.0, 3.0).unwrap();
        let h = film.thickness_m();
        near(h.powi(3), 3.0 * 2e-6 * 2e-6 / 4.0);
        near(film.reynolds(), 1.0);
        near(film.wall_shear_pa(), 1000.0 * 4.0 * h);
        // Independent h-power coefficients rather than the native factored Re form.
        let k: f64 = 10.0;
        let nu: f64 = 2e-6;
        let expected = (2.0 * 4.0_f64.powi(2) * h.powi(6) / (15.0 * nu.powi(3))
            - 3.0 * h.powi(3) / (3.0 * nu))
            * k.powi(2)
            - 0.08 * h.powi(3) / (3.0 * 1000.0 * nu) * k.powi(4);
        near(
            film.long_wave_mode(k, LIMITS).unwrap().growth_per_s,
            expected,
        );
    }
}

#[test]
fn pinned_temperature_update_changes_all_coefficients_without_losing_context() {
    let card = card(LiquidViscosityBasis::Dynamic, true, false);
    let cold = resolve(
        &card,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
    );
    let hot = ResolvedLiquidSurfaceState::try_from_material(
        cold.material()
            .at_query_point(&card, &point(400.0))
            .unwrap(),
        cold.viscosity_basis(),
        cold.capillarity(),
    )
    .unwrap();
    let first = MaterialFilmLiquid::try_new(&cold).unwrap();
    let second = MaterialFilmLiquid::try_new(&hot).unwrap();
    assert_eq!(second.coefficients().density_kg_m3, 800.0);
    near(second.coefficients().kinematic_viscosity_m2_s, 5e-6);
    assert_eq!(second.coefficients().surface_tension_n_m, 0.04);
    let a = first.film(2e-6, 4.0, 3.0).unwrap();
    let b = second.film(2e-6, 4.0, 3.0).unwrap();
    near(b.thickness_m().powi(3) / a.thickness_m().powi(3), 2.5);
    near(b.reynolds() / a.reynolds(), 0.4);
    assert!(a.long_wave_mode(10.0, LIMITS).unwrap().growth_per_s > 0.0);
    assert!(b.long_wave_mode(10.0, LIMITS).unwrap().growth_per_s < 0.0);
    for old in cold.material().properties() {
        let refreshed = hot.material().property(old.requirement().name()).unwrap();
        assert_eq!(
            refreshed.answer().receipt.selected,
            old.answer().receipt.selected
        );
        assert_eq!(
            refreshed.answer().receipt.axis_quantities,
            *point(400.0).axis_quantities()
        );
        assert!(
            refreshed
                .answer()
                .receipt
                .query_point
                .contains(&("surface_environment".into(), 7.0))
        );
    }
    assert!(
        cold.material()
            .at_query_point(&card, &point(401.0))
            .is_err()
    );
    assert_eq!(first.coefficients().density_kg_m3, 1000.0);
}

#[test]
fn capillarity_is_omitted_only_by_the_original_explicit_declaration() {
    let card = card(LiquidViscosityBasis::Dynamic, false, false);
    let source = resolve(
        &card,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
    );
    let bound = MaterialFilmLiquid::try_new(&source).unwrap();
    assert_eq!(bound.material().capillarity(), LiquidCapillarity::Neglected);
    assert!(
        bound
            .material()
            .material()
            .property(SURFACE_TENSION_PROPERTY)
            .is_none()
    );
    assert_eq!(bound.coefficients().surface_tension_n_m, 0.0);
    near(
        bound
            .film(2e-6, 4.0, 3.0)
            .unwrap()
            .long_wave_mode(10.0, LIMITS)
            .unwrap()
            .growth_per_s,
        2e-6 * 100.0 * (1.2 - 0.75),
    );
    assert!(
        resolve_liquid_surface_state(
            &card,
            &point(300.0),
            LiquidViscosityBasis::Dynamic,
            LiquidCapillarity::Sourced,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}

#[test]
fn apparent_viscosity_requiring_shear_or_pressure_evolution_is_never_frozen() {
    let card = card(LiquidViscosityBasis::Dynamic, true, true);
    let point = point(300.0)
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 1.0)
        .unwrap();
    let requirements = [
        (DENSITY_PROPERTY, Density::DIMS),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS),
        (SURFACE_TENSION_PROPERTY, SURFACE_TENSION_DIMS),
    ]
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    });
    let material = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let shear = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution_with_capillarity(
        material.clone(),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
        "T",
        "rate",
    )
    .unwrap();
    let pressure = ResolvedLiquidSurfaceState::try_for_pressure_shear_rate_evolution(
        material,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
        "T",
        "rate",
        "pressure",
    )
    .unwrap();
    for source in [&shear, &pressure] {
        assert!(matches!(
            MaterialFilmLiquid::try_new(source),
            Err(FilmError::InvalidInput(_))
        ));
    }
}
