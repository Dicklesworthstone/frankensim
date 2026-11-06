//! G0/G1 source-to-film and source-to-jet coupling with synthetic constitutive inputs.
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
        assert!(matches!(
            fs_lbm::jet::material::MaterialCapillaryJet::inviscid(source, 0.001, 0.1),
            Err(fs_lbm::jet::JetError::InvalidInput(_))
        ));
        assert!(matches!(
            fs_lbm::jet::material::MaterialViscousSlenderJet::try_new(source, 0.001, 0.3),
            Err(fs_lbm::jet::JetError::InvalidInput(_))
        ));
        assert!(matches!(
            fs_lbm::jet::material::nonlinear::MaterialPeriodicSlenderJet::try_new(
                source,
                0.01,
                &[0.001; 64],
                &[0.0; 64],
                fs_lbm::jet::nonlinear::SlenderJetLimits {
                    minimum_radius_m: 1e-4,
                    maximum_abs_slope: 0.5,
                    minimum_radius_cells: 2.0,
                },
                &fs_exec::CancelGate::new_clock_free(),
            ),
            Err(fs_lbm::jet::nonlinear::NonlinearJetError::InvalidInput(_))
        ));
    }
}

#[test]
fn cylindrical_jet_reuses_both_source_viscosity_bases_and_original_receipts() {
    use fs_lbm::jet::JetDispersion;
    use fs_lbm::jet::material::MaterialCapillaryJet;

    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let card = card(basis, true, false);
        let source = resolve(&card, basis, LiquidCapillarity::Sourced);
        let bound = MaterialCapillaryJet::inviscid(&source, 0.001, 0.01).unwrap();
        assert!(core::ptr::eq(bound.material(), &source));
        assert_eq!(
            bound.material().material().card_identity(),
            card.content_hash()
        );
        assert_eq!(bound.material().viscosity_basis(), basis);
        assert_eq!(bound.material().material().properties().len(), 3);
        near(
            bound.ohnesorge(),
            0.002 / (1000.0_f64 * 0.08 * 0.001).sqrt(),
        );
        assert_eq!(bound.maximum_ohnesorge(), 0.01);
        let maximum = bound.fastest_mode([0.0, 1000.0]).unwrap();
        near(maximum.mode.kr, 0.697_018_898_306_3);
        // Classical dimensionless inviscid maximum, independent of the material adapter.
        let rate_scale = (0.08_f64 / (1000.0 * 0.001_f64.powi(3))).sqrt();
        assert!((maximum.mode.growth_per_s / rate_scale - 0.343_338_701_861_3).abs() < 2e-12);
    }
}

#[test]
fn cylindrical_jet_source_temperature_changes_growth_and_viscosity_admission() {
    use fs_lbm::jet::material::MaterialCapillaryJet;

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
    let first = MaterialCapillaryJet::inviscid(&cold, 0.001, 0.01).unwrap();
    // The SAME original source at 400 K has more viscosity and less tension:
    // the old inviscid allowance no longer admits it.
    assert!(MaterialCapillaryJet::inviscid(&hot, 0.001, 0.01).is_err());
    let second = MaterialCapillaryJet::inviscid(&hot, 0.001, 0.03).unwrap();
    let a = first.jet().mode(500.0).unwrap();
    let b = second.jet().mode(500.0).unwrap();
    near(
        b.growth_per_s / a.growth_per_s,
        (0.04_f64 / 800.0 / (0.08 / 1000.0)).sqrt(),
    );
    near(
        second.ohnesorge(),
        0.004 / (800.0_f64 * 0.04 * 0.001).sqrt(),
    );
    assert!(core::ptr::eq(second.material(), &hot));
    for previous in cold.material().properties() {
        let refreshed = second
            .material()
            .material()
            .property(previous.requirement().name())
            .unwrap();
        assert_eq!(
            refreshed.answer().receipt.selected,
            previous.answer().receipt.selected
        );
        assert_eq!(
            refreshed.answer().receipt.axis_quantities,
            *point(400.0).axis_quantities()
        );
    }
}

#[test]
fn cylindrical_material_jet_requires_active_tension_and_explicit_valid_budget() {
    use fs_lbm::jet::material::MaterialCapillaryJet;

    let without_surface = card(LiquidViscosityBasis::Dynamic, false, false);
    let neglected = resolve(
        &without_surface,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
    );
    assert!(MaterialCapillaryJet::inviscid(&neglected, 0.001, 0.1).is_err());
    let with_surface = card(LiquidViscosityBasis::Dynamic, true, false);
    let source = resolve(
        &with_surface,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
    );
    for budget in [-1.0, f64::NAN, f64::INFINITY, 0.0, 0.001] {
        assert!(MaterialCapillaryJet::inviscid(&source, 0.001, budget).is_err());
    }
    let bound = MaterialCapillaryJet::inviscid(&source, 0.001, 0.01).unwrap();
    assert!(MaterialCapillaryJet::inviscid(&source, 0.001, bound.ohnesorge()).is_ok());
    // Shrinking the same sourced jet raises Oh; units are physical radius.
    assert!(MaterialCapillaryJet::inviscid(&source, 0.0001, 0.01).is_err());
}

#[test]
fn viscous_jet_uses_both_source_bases_in_the_dimensional_dispersion() {
    use fs_lbm::jet::JetDispersion;
    use fs_lbm::jet::material::MaterialViscousSlenderJet;

    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let card = card(basis, true, false);
        let source = resolve(&card, basis, LiquidCapillarity::Sourced);
        let bound = MaterialViscousSlenderJet::try_new(&source, 0.001, 0.3).unwrap();
        assert!(core::ptr::eq(bound.material(), &source));
        assert_eq!(
            bound.material().material().card_identity(),
            card.content_hash()
        );
        assert_eq!(bound.material().viscosity_basis(), basis);
        assert_eq!(bound.material().material().properties().len(), 3);
        assert_eq!(bound.jet().maximum_kr(), 0.3);
        near(bound.jet().kinematic_viscosity_m2_s(), 2e-6);
        let maximum = bound.fastest_mode([0.0, 250.0]).unwrap();
        let k = maximum.wavenumber_rad_m;
        let s = maximum.mode.growth_per_s;
        let radius = bound.radius_m();
        // Independent dimensional momentum/mass-balance residual.
        near(
            s * s + 3.0 * 2e-6 * k * k * s,
            0.08 * k * k * (1.0 - radius * radius * k * k) / (2.0 * 1000.0 * radius),
        );
        assert!(bound.fastest_mode([0.0, 400.0]).is_err());
    }
}

#[test]
fn viscous_source_temperature_is_retained_when_the_inviscid_budget_refuses() {
    use fs_lbm::jet::JetDispersion;
    use fs_lbm::jet::material::{MaterialCapillaryJet, MaterialViscousSlenderJet};

    let source_card = card(LiquidViscosityBasis::Dynamic, true, false);
    let cold = resolve(
        &source_card,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
    );
    let hot = ResolvedLiquidSurfaceState::try_from_material(
        cold.material()
            .at_query_point(&source_card, &point(400.0))
            .unwrap(),
        cold.viscosity_basis(),
        cold.capillarity(),
    )
    .unwrap();
    let radius = 1e-7;
    for source in [&cold, &hot] {
        assert!(MaterialCapillaryJet::inviscid(source, radius, 0.1).is_err());
    }
    let first = MaterialViscousSlenderJet::try_new(&cold, radius, 0.3).unwrap();
    let second = MaterialViscousSlenderJet::try_new(&hot, radius, 0.3).unwrap();
    near(second.jet().density_kg_m3(), 800.0);
    near(second.jet().surface_tension_n_m(), 0.04);
    near(second.jet().kinematic_viscosity_m2_s(), 5e-6);
    assert!(core::ptr::eq(second.material(), &hot));
    let before = first.fastest_mode([0.0, 2.5e6]).unwrap();
    let after = second.fastest_mode([0.0, 2.5e6]).unwrap();
    assert!(after.mode.growth_per_s < before.mode.growth_per_s);
    for previous in cold.material().properties() {
        let refreshed = second
            .material()
            .material()
            .property(previous.requirement().name())
            .unwrap();
        assert_eq!(
            refreshed.answer().receipt.selected,
            previous.answer().receipt.selected
        );
        assert_eq!(
            refreshed.answer().receipt.axis_quantities,
            *point(400.0).axis_quantities()
        );
    }
    let without_surface = card(LiquidViscosityBasis::Dynamic, false, false);
    let neglected = resolve(
        &without_surface,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
    );
    assert!(MaterialViscousSlenderJet::try_new(&neglected, radius, 0.3).is_err());
    for budget in [0.0, -0.1, 1.01, f64::NAN, f64::INFINITY] {
        assert!(MaterialViscousSlenderJet::try_new(&cold, radius, budget).is_err());
    }
}

fn nonlinear_jet_fixture() -> (
    f64,
    Vec<f64>,
    Vec<f64>,
    fs_lbm::jet::nonlinear::SlenderJetLimits,
) {
    let length = 4.0 * core::f64::consts::PI * 0.001;
    let radii = (0..64)
        .map(|i| 0.001 * (1.0 + 0.05 * fs_math::det::cos(core::f64::consts::TAU * i as f64 / 64.0)))
        .collect();
    let limits = fs_lbm::jet::nonlinear::SlenderJetLimits {
        minimum_radius_m: 1e-4,
        maximum_abs_slope: 0.5,
        minimum_radius_cells: 2.0,
    };
    (length, radii, vec![0.0; 64], limits)
}

fn nonlinear_jet_budget() -> fs_lbm::jet::nonlinear::SlenderJetStepBudget {
    fs_lbm::jet::nonlinear::SlenderJetStepBudget {
        maximum_step_s: 1e-5,
        cfl: 0.1,
        maximum_substeps: 20_000,
    }
}

#[test]
fn nonlinear_jet_evolves_the_original_material_and_conserves_volume() {
    use fs_exec::CancelGate;
    use fs_lbm::jet::material::nonlinear::MaterialPeriodicSlenderJet;
    use fs_lbm::jet::nonlinear::{SlenderJetEvolution, SlenderJetOutcome};

    let (length, radii, velocities, limits) = nonlinear_jet_fixture();
    let gate = CancelGate::new_clock_free();
    let mut reference_state: Option<Vec<f64>> = None;
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let source_card = card(basis, true, false);
        let source = resolve(&source_card, basis, LiquidCapillarity::Sourced);
        let mut bound = MaterialPeriodicSlenderJet::try_new(
            &source,
            length,
            &radii,
            &velocities,
            limits,
            &gate,
        )
        .unwrap();
        assert!(core::ptr::eq(bound.material(), &source));
        assert_eq!(
            bound.material().material().card_identity(),
            source_card.content_hash()
        );
        assert_eq!(bound.material().viscosity_basis(), basis);
        let parameters = bound.state().parameters();
        near(parameters.density_kg_m3, 1000.0);
        near(parameters.kinematic_viscosity_m2_s, 2e-6);
        near(parameters.surface_tension_n_m, 0.08);
        let expected_volume =
            core::f64::consts::PI * length / 64.0 * radii.iter().map(|r| r * r).sum::<f64>();
        near(bound.state().diagnostics().volume_m3, expected_volume);
        // Exercise the same trait consumed by vessel flight evolution.
        let report =
            SlenderJetEvolution::advance(&mut bound, 0.004, nonlinear_jet_budget(), &gate).unwrap();
        assert_eq!(report.outcome, SlenderJetOutcome::Completed);
        near(bound.state().time_s(), 0.004);
        near(report.after.volume_m3, expected_volume);
        assert!(report.after.minimum_radius_m < report.before.minimum_radius_m);
        assert!(bound.state().velocities_m_s().any(|u| u.abs() > 1e-8));
        let evolved = bound.state().radius_squared_m2().to_vec();
        if let Some(previous) = &reference_state {
            for (&actual, &expected) in evolved.iter().zip(previous) {
                near(actual, expected);
            }
        }
        reference_state = Some(evolved);
        assert!(core::ptr::eq(bound.material(), &source));
    }
}

#[test]
fn nonlinear_jet_temperature_changes_the_evolved_shape_from_original_sources() {
    use fs_exec::CancelGate;
    use fs_lbm::jet::material::nonlinear::MaterialPeriodicSlenderJet;

    let source_card = card(LiquidViscosityBasis::Dynamic, true, false);
    let cold = resolve(
        &source_card,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
    );
    let hot = ResolvedLiquidSurfaceState::try_from_material(
        cold.material()
            .at_query_point(&source_card, &point(400.0))
            .unwrap(),
        cold.viscosity_basis(),
        cold.capillarity(),
    )
    .unwrap();
    let (length, radii, velocities, limits) = nonlinear_jet_fixture();
    let gate = CancelGate::new_clock_free();
    let mut first =
        MaterialPeriodicSlenderJet::try_new(&cold, length, &radii, &velocities, limits, &gate)
            .unwrap();
    let mut second =
        MaterialPeriodicSlenderJet::try_new(&hot, length, &radii, &velocities, limits, &gate)
            .unwrap();
    first.advance(0.004, nonlinear_jet_budget(), &gate).unwrap();
    second
        .advance(0.004, nonlinear_jet_budget(), &gate)
        .unwrap();
    let p = second.state().parameters();
    near(p.density_kg_m3, 800.0);
    near(p.kinematic_viscosity_m2_s, 5e-6);
    near(p.surface_tension_n_m, 0.04);
    assert!(
        second.state().diagnostics().minimum_radius_m
            > first.state().diagnostics().minimum_radius_m
    );
    near(
        second.state().diagnostics().volume_m3,
        first.state().diagnostics().volume_m3,
    );
    assert!(core::ptr::eq(second.material(), &hot));
    for previous in cold.material().properties() {
        let refreshed = second
            .material()
            .material()
            .property(previous.requirement().name())
            .unwrap();
        assert_eq!(
            refreshed.answer().receipt.selected,
            previous.answer().receipt.selected
        );
        assert_eq!(
            refreshed.answer().receipt.axis_quantities,
            *point(400.0).axis_quantities()
        );
    }
}

#[test]
fn nonlinear_material_jet_refuses_missing_tension_and_preserves_state_on_cancel() {
    use fs_exec::CancelGate;
    use fs_lbm::jet::material::nonlinear::MaterialPeriodicSlenderJet;
    use fs_lbm::jet::nonlinear::NonlinearJetError;

    let (length, radii, velocities, limits) = nonlinear_jet_fixture();
    let gate = CancelGate::new_clock_free();
    let missing_card = card(LiquidViscosityBasis::Dynamic, false, false);
    let missing = resolve(
        &missing_card,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
    );
    assert!(matches!(
        MaterialPeriodicSlenderJet::try_new(&missing, length, &radii, &velocities, limits, &gate),
        Err(NonlinearJetError::InvalidInput(_))
    ));
    let source_card = card(LiquidViscosityBasis::Dynamic, true, false);
    let source = resolve(
        &source_card,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
    );
    let mut bound =
        MaterialPeriodicSlenderJet::try_new(&source, length, &radii, &velocities, limits, &gate)
            .unwrap();
    let before_radius = bound.state().radius_squared_m2().to_vec();
    let before_flux = bound.state().radius_squared_flux_m3_s().to_vec();
    let before_time = bound.state().time_s();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        bound.advance(0.004, nonlinear_jet_budget(), &cancelled),
        Err(NonlinearJetError::Cancelled)
    ));
    assert_eq!(bound.state().radius_squared_m2(), before_radius);
    assert_eq!(bound.state().radius_squared_flux_m3_s(), before_flux);
    assert_eq!(bound.state().time_s(), before_time);
    assert!(core::ptr::eq(bound.material(), &source));
}
