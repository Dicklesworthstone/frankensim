//! G0/G1/G3/G4 source-to-pour regressions with synthetic material inputs.
//! These exercise dimensional wiring and native coupling, not measured fluids.

#[path = "material_pour/rheology.rs"]
mod rheology;

use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_lbm::freesurface::ContactModel;
use fs_lbm::rheology::Rheology;
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
use fs_material::state_point::{DENSITY_PROPERTY, MaterialPropertySelection};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Frequency, QuantitySpec, Temperature};
use fs_vessel::pour::run_pour_in_frame;
use fs_vessel::{MaterialPour, MaterialPourError, MaterialPourSpec};

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

fn source(
    basis: LiquidViscosityBasis,
    surface: Option<f64>,
) -> (MaterialCard, ResolvedLiquidSurfaceState) {
    let (viscosity, dims, values) = match basis {
        LiquidViscosityBasis::Dynamic => (
            DYNAMIC_VISCOSITY_PROPERTY,
            DYNAMIC_VISCOSITY_DIMS,
            [0.02, 0.024],
        ),
        LiquidViscosityBasis::Kinematic => (
            KINEMATIC_VISCOSITY_PROPERTY,
            KINEMATIC_VISCOSITY_DIMS,
            [2e-5, 3e-5],
        ),
    };
    let mut rows = vec![
        (DENSITY_PROPERTY, Density::DIMS, [1000.0, 800.0]),
        (viscosity, dims, values),
    ];
    if let Some(sigma) = surface {
        rows.push((
            SURFACE_TENSION_PROPERTY,
            SURFACE_TENSION_DIMS,
            [sigma, 0.5 * sigma],
        ));
    }
    let mut claims = ClaimSet::new();
    for (name, dims, values) in rows {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (400.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        *point(300.0).axis_quantities().get("T").unwrap(),
                        300.0,
                        400.0,
                    )
                    .with("surface_environment", 7.0, 7.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic SI pour fixture, not measurements".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic pour liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let state = resolve_liquid_surface_state(
        &card,
        &point(300.0),
        basis,
        if surface.is_some() {
            LiquidCapillarity::Sourced
        } else {
            LiquidCapillarity::Neglected
        },
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (card, state)
}

fn spec() -> MaterialPourSpec {
    MaterialPourSpec {
        cells: [14, 10],
        lip: [7, 4],
        cell_size_m: 0.001,
        timestep_s: 0.001,
        gravity_m_s2: 0.02,
        tilt_final_radians: 0.0004,
        pivot_m: [0.007, 0.001],
        steps: 32,
        maximum_cells: 140,
        maximum_speed_m_s: 0.05,
        maximum_density_excursion: 0.1,
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4e-12 * expected.abs().max(1e-30),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn g1_source_pour_retains_receipts_and_the_exact_native_field_in_si_units() {
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let (card, liquid) = source(basis, Some(0.002));
        let input = spec();
        let prepared = MaterialPour::try_new(&liquid, input).unwrap();
        assert!(core::ptr::eq(prepared.scaling().material(), &liquid));
        assert_eq!(
            prepared.scaling().material().material().card_identity(),
            card.content_hash()
        );
        assert_eq!(prepared.pivot_lattice(), [6.5, 0.5]);
        let dx = input.cell_size_m;
        let dt = input.timestep_s;
        near(prepared.native_rig().g0, input.gravity_m_s2 * dt * dt / dx);
        near(
            prepared.scaling().flow().viscosity,
            liquid.kinematic_viscosity_m2_s() * dt / (dx * dx),
        );
        near(
            prepared.native_rig().surface_tension,
            liquid.surface_tension_n_m().unwrap() * dt * dt
                / (liquid.density_kg_m3() * dx * dx * dx),
        );
        for contact in [ContactModel::Neutral, ContactModel::Wetting] {
            let actual = prepared
                .run(contact, &CancelGate::new_clock_free())
                .unwrap();
            let direct = run_pour_in_frame(
                prepared.native_rig(),
                contact,
                Rheology::Newtonian {
                    nu: prepared.scaling().flow().viscosity,
                },
                prepared.pivot_lattice(),
            )
            .unwrap();
            assert_eq!(actual.native().mass_field, direct.mass_field);
            assert_eq!(actual.native().mass_drift, direct.mass_drift);
            assert_eq!(actual.native().poured_mass, direct.poured_mass);
            assert_eq!(actual.native().fragments, direct.fragments);
            assert_eq!(actual.native().dribble_cells, direct.dribble_cells);
            assert!(core::ptr::eq(actual.scaling().material(), &liquid));
            let report = actual.report();
            let unit = liquid.density_kg_m3() * dx * dx;
            near(report.duration_s, f64::from(input.steps) * dt);
            near(report.poured_mass_kg_per_m, direct.poured_mass * unit);
            near(
                report.final_mass_kg_per_m,
                direct.mass_field.iter().sum::<f64>() * unit,
            );
            assert!(
                (report.final_mass_kg_per_m / report.initial_mass_kg_per_m - 1.0).abs() < 1e-10
            );
            assert!(
                report.peak_speed_m_s > 0.0 && report.peak_speed_m_s <= input.maximum_speed_m_s
            );
            assert!(report.maximum_density_excursion <= input.maximum_density_excursion);
            for (physical, native) in actual.mass_field_kg_per_m().zip(&direct.mass_field) {
                near(physical, native * unit);
            }
        }
    }
}

#[test]
fn g3_the_original_temperature_source_changes_viscosity_capillarity_and_motion() {
    let (card, cold) = source(LiquidViscosityBasis::Dynamic, Some(0.002));
    let hot = ResolvedLiquidSurfaceState::try_from_material(
        cold.material()
            .at_query_point(&card, &point(400.0))
            .unwrap(),
        cold.viscosity_basis(),
        cold.capillarity(),
    )
    .unwrap();
    let first = MaterialPour::try_new(&cold, spec()).unwrap();
    let second = MaterialPour::try_new(&hot, spec()).unwrap();
    assert_eq!(
        cold.material().card_identity(),
        hot.material().card_identity()
    );
    near(
        second.scaling().flow().viscosity / first.scaling().flow().viscosity,
        1.5,
    );
    near(
        second.native_rig().surface_tension / first.native_rig().surface_tension,
        0.625,
    );
    let gate = CancelGate::new_clock_free();
    let a = first.run(ContactModel::Neutral, &gate).unwrap();
    let b = second.run(ContactModel::Neutral, &gate).unwrap();
    assert_ne!(a.native().mass_field, b.native().mass_field);
    near(
        b.report().initial_mass_kg_per_m / a.report().initial_mass_kg_per_m,
        0.8,
    );
    assert!(core::ptr::eq(b.scaling().material(), &hot));
}

#[test]
fn g3_two_dimensional_mass_scales_with_area_while_time_and_velocity_keep_their_units() {
    let (_, liquid) = source(LiquidViscosityBasis::Kinematic, None);
    let input = spec();
    let rescaled = MaterialPourSpec {
        cell_size_m: 2.0 * input.cell_size_m,
        timestep_s: 4.0 * input.timestep_s,
        gravity_m_s2: input.gravity_m_s2 / 8.0,
        pivot_m: input.pivot_m.map(|value| 2.0 * value),
        maximum_speed_m_s: 0.5 * input.maximum_speed_m_s,
        ..input
    };
    let a = MaterialPour::try_new(&liquid, input).unwrap();
    let b = MaterialPour::try_new(&liquid, rescaled).unwrap();
    assert_eq!(a.native_rig().surface_tension, 0.0);
    assert_eq!(b.native_rig().surface_tension, 0.0);
    let gate = CancelGate::new_clock_free();
    let a = a.run(ContactModel::Neutral, &gate).unwrap();
    let b = b.run(ContactModel::Neutral, &gate).unwrap();
    for (&x, &y) in a.native().mass_field.iter().zip(&b.native().mass_field) {
        near(x, y);
    }
    near(b.report().duration_s / a.report().duration_s, 4.0);
    near(
        b.report().final_mass_kg_per_m / a.report().final_mass_kg_per_m,
        4.0,
    );
    near(b.report().peak_speed_m_s / a.report().peak_speed_m_s, 0.5);
}

#[test]
fn g0_source_rheology_capillary_and_geometric_refusals_are_not_silently_substituted() {
    let (card, liquid) = source(LiquidViscosityBasis::Dynamic, Some(0.002));
    let input = spec();
    for invalid in [
        MaterialPourSpec {
            cells: [4, 10],
            ..input
        },
        MaterialPourSpec {
            lip: [1, 4],
            ..input
        },
        MaterialPourSpec {
            lip: [7, 9],
            ..input
        },
        MaterialPourSpec {
            maximum_cells: 139,
            ..input
        },
        MaterialPourSpec { steps: 0, ..input },
        MaterialPourSpec {
            timestep_s: f64::NAN,
            ..input
        },
        MaterialPourSpec {
            timestep_s: 1e-8,
            ..input
        },
        MaterialPourSpec {
            timestep_s: 0.01,
            maximum_speed_m_s: 0.005,
            ..input
        },
        MaterialPourSpec {
            gravity_m_s2: -9.81,
            ..input
        },
        MaterialPourSpec {
            pivot_m: [f64::INFINITY, 0.0],
            ..input
        },
        MaterialPourSpec {
            maximum_speed_m_s: 1.0,
            ..input
        },
        MaterialPourSpec {
            maximum_density_excursion: 1.0,
            ..input
        },
    ] {
        assert!(
            MaterialPour::try_new(&liquid, invalid).is_err(),
            "admitted {invalid:?}"
        );
    }
    let query = point(300.0)
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0)
        .unwrap();
    let material = liquid.material().at_query_point(&card, &query).unwrap();
    let shear = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution_with_capillarity(
        material.clone(),
        liquid.viscosity_basis(),
        liquid.capillarity(),
        "T",
        "rate",
    )
    .unwrap();
    let pressure = ResolvedLiquidSurfaceState::try_for_pressure_shear_rate_evolution(
        material,
        liquid.viscosity_basis(),
        liquid.capillarity(),
        "T",
        "rate",
        "pressure",
    )
    .unwrap();
    for evolving in [&shear, &pressure] {
        assert!(matches!(
            MaterialPour::try_new(evolving, input),
            Err(MaterialPourError::Liquid(_))
        ));
    }
}

#[test]
fn g4_actual_endpoint_limits_and_cancellation_return_no_partial_result() {
    let (_, liquid) = source(LiquidViscosityBasis::Dynamic, None);
    let original = liquid.clone();
    let slow = MaterialPour::try_new(
        &liquid,
        MaterialPourSpec {
            maximum_speed_m_s: 1e-12,
            ..spec()
        },
    )
    .unwrap();
    assert!(
        matches!(slow.run(ContactModel::Neutral, &CancelGate::new_clock_free()),
        Err(MaterialPourError::State { step, .. }) if step > 0)
    );
    let incompressible = MaterialPour::try_new(
        &liquid,
        MaterialPourSpec {
            maximum_density_excursion: 1e-14,
            ..spec()
        },
    )
    .unwrap();
    assert!(
        matches!(incompressible.run(ContactModel::Neutral, &CancelGate::new_clock_free()),
        Err(MaterialPourError::State { step, .. }) if step > 0)
    );
    let prepared = MaterialPour::try_new(&liquid, spec()).unwrap();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        prepared.run(ContactModel::Neutral, &cancelled),
        Err(MaterialPourError::Cancelled)
    ));
    let gate = CancelGate::new_clock_free();
    let a = prepared.run(ContactModel::Neutral, &gate).unwrap();
    let b = prepared.run(ContactModel::Neutral, &gate).unwrap();
    assert_eq!(a.native().mass_field, b.native().mass_field);
    assert_eq!(a.report(), b.report());
    assert_eq!(liquid, original);
}
