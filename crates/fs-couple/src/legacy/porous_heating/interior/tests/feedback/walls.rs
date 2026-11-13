//! Original temperature-dependent wall sources in the native FE/heat timeline.
//! These authored material curves test composition, not measured wall accuracy.

use super::*;
use crate::bernoulli_aperture::wall::material::{
    AxialConstraint, LocalMaterialWall, WallMaterialLoss, WallRelaxation, YOUNG_VISCOSITY_PROPERTY,
};
use crate::cavity_modal::interior::heterogeneous::AxialGasEmbedding;
use crate::cavity_modal::interior::heterogeneous::mechanical::{
    CylindricalMaterialBoundary, MaterialMechanicalError, MaterialMechanicalHarmonic,
    MaterialMechanicalHarmonicError,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, DynViscosity, Pressure, QuantitySpec, Temperature};

const RADIUS: f64 = 0.1;
const THICKNESS: f64 = 0.0006;
const BAND: [f64; 2] = [100.0, 10_000.0];
const DT: f64 = 0.03125;

#[derive(Clone, Copy)]
enum Loss {
    Kelvin,
    Maxwell,
}

fn card(density_slope: f64, elastic_upper_k: f64) -> MaterialCard {
    let temperatures = [280.0, 300.0, 400.0];
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        (
            "density",
            Density::DIMS,
            temperatures.map(|t| 1000.0 * (1.0 + density_slope * (t - 300.0))),
        ),
        (
            "young_modulus",
            Pressure::DIMS,
            temperatures.map(|t| 2e8 * (1.0 - 0.002 * (t - 300.0))),
        ),
        ("poisson_ratio", Dims::NONE, [0.3; 3]),
        (
            YOUNG_VISCOSITY_PROPERTY,
            DynViscosity::DIMS,
            temperatures.map(|t| 5000.0 * (1.0 + 0.003 * (t - 300.0))),
        ),
        (
            "arm_modulus",
            Pressure::DIMS,
            temperatures.map(|t| 4e7 * (1.0 + 0.004 * (t - 300.0))),
        ),
        (
            "arm_time",
            Dims([0, 0, 1, 0, 0, 0]),
            temperatures.map(|t| 0.0002 * (1.0 + 0.002 * (t - 300.0))),
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: temperatures.into_iter().zip(values).collect(),
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::AbsoluteTemperature,
                            ValueForm::Static,
                        )),
                        280.0,
                        if name == "young_modulus" {
                            elastic_upper_k
                        } else {
                            400.0
                        },
                    )
                    .with("lot", 7.0, 7.0)
                    .with("omega", 0.0, core::f64::consts::TAU * BAND[1]),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored thermal membrane curves; not calibrated material data".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal membrane".into(),
            phase: "solid".into(),
            process: "fixed inventory; original stiffness, viscosity and reservoir laws".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

// Independent source resolution, without the feedback continuation helper.
fn source(card: &MaterialCard, temperature: f64, loss: Loss) -> LocalMaterialWall {
    let requirements = [
        ("density", Density::DIMS),
        ("young_modulus", Pressure::DIMS),
        ("poisson_ratio", Dims::NONE),
        (YOUNG_VISCOSITY_PROPERTY, DynViscosity::DIMS),
        ("arm_modulus", Pressure::DIMS),
        ("arm_time", Dims([0, 0, 1, 0, 0, 0])),
    ]
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap()
    });
    let state = resolve_material_state_point(
        card,
        &QueryPoint::new()
            .with_quantity(
                "T",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                temperature,
            )
            .unwrap()
            .with("lot", 7.0)
            .unwrap()
            .with("omega", core::f64::consts::TAU * FREQUENCY)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let loss = match loss {
        Loss::Kelvin => WallMaterialLoss::KelvinVoigt,
        Loss::Maxwell => WallMaterialLoss::Maxwell(
            WallRelaxation::from_state_point(
                &state,
                &[("arm_modulus".into(), "arm_time".into())],
                BAND,
                0.25,
                0.4,
            )
            .unwrap(),
        ),
    };
    LocalMaterialWall::new(state, THICKNESS, AxialConstraint::StressFree, loss, BAND).unwrap()
}

fn cylindrical(source: &LocalMaterialWall) -> CylindricalMaterialBoundary<'_> {
    // The existing triangular side is a declared cylindrical facet: relative
    // radial error <0.012 and outward-normal error <0.15 at every vertex.
    CylindricalMaterialBoundary {
        faces: vec![[0, 1, 3]],
        source,
        axis: AxialGasEmbedding::new([0.5 * SCALE, RADIUS, 0.0], [0.0, 0.0, 1.0]).unwrap(),
        radius_m: RADIUS,
        relative_geometry_tolerance: 0.16,
    }
}

fn bindings<'t, 'm>(
    caloric: &'t dyn ThermalLaw,
    card: &'m MaterialCard,
    source: &'m LocalMaterialWall,
) -> [TemperatureInteriorWall<'t, 'm>; 1] {
    [TemperatureInteriorWall {
        destination: wall_routes(caloric)[0],
        boundary: cylindrical(source),
        card,
        temperature_axis: "T",
    }]
}

fn wall_events(h: [f64; 2], watts: [f64; 2], wall_h: f64) -> Vec<TimedEvent> {
    let mut events = events(h, watts, true);
    let Event::Birth(wall) = &mut events[2].event else {
        panic!("wall birth")
    };
    wall.specific_h_j_kg = wall_h;
    wall.geometry.end_m = [0.0, THICKNESS, 0.0];
    wall.geometry.volume_m3 = WALL_AREA * THICKNESS;
    events
}

fn gases(thermal: &DepositionSimulation<'_>) -> [GasState; 2] {
    IDS.map(|id| GasState::try_new(&SPEC, thermal.temperature_k(id).unwrap(), 101325.0).unwrap())
}

fn original<'f, 'c, 's>(
    geometry: &'f HelmholtzFe<'c>,
    states: &'s [GasState; 2],
    porous: Option<&'s [PorousFluidLaw; 2]>,
    wall: &'s LocalMaterialWall,
) -> MaterialMechanicalHarmonic<'f, 'c, 's> {
    let sources = core::array::from_fn::<_, 2, _>(|i| match porous {
        Some(laws) => AcousticVolumeRegion::Porous {
            law: &laws[i],
            fluid: &states[i],
        },
        None => AcousticVolumeRegion::Fluid(&states[i]),
    });
    MaterialMechanicalHarmonic::new(
        geometry,
        FREQUENCY,
        &sources,
        &[0, 1],
        &[],
        vec![cylindrical(wall)],
        &[],
        EquivalentFluidControls::default(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn wall_only_temperature_changes_refresh_original_kelvin_and_maxwell_without_resampling_gas() {
    let (mesh, points) = fixture();
    let geometry = geometry(&mesh, &points);
    let caloric = caloric();
    let card = card(0.0, 400.0);
    for (loss, initial_wall_h) in [(Loss::Kelvin, 0.0), (Loss::Maxwell, 30.0 * CP)] {
        let wall = source(&card, 300.0, loss);
        let thermal = initial(wall_events([0.0; 2], [0.0; 2], initial_wall_h), &caloric);
        let gas = fluid([280.0, 500.0]);
        let calls = gas.calls.clone();
        let heating = TemperatureInteriorHeating::try_new_with_material_walls(
            &geometry,
            &thermal,
            &mapping(&caloric, None),
            &bindings(&caloric, &card, &wall),
            gas,
            problem(false),
            3,
            || Ok(()),
        )
        .unwrap();
        let result = heating
            .advance_to(&thermal, 2.0 * DT, config(DT), || Ok(()))
            .unwrap();
        assert_eq!(result.work.fluid_evaluations, 0);
        assert_eq!(result.work.porous_evaluations, 0);
        assert_eq!(result.work.wall_evaluations, 2);
        assert_eq!(result.work.field_solves, 2);
        assert_eq!(calls.get(), 2);
        assert_eq!(result.heating.fluid_states(), heating.fluid_states());
        assert_ne!(
            result.heating.field().pressure_pa,
            heating.field().pressure_pa
        );
        assert!(
            result.step.thermal.temperature_k(WALL).unwrap() > thermal.temperature_k(WALL).unwrap()
        );
        assert!(core::ptr::eq(
            heating.wall_bindings()[0].boundary.source,
            &wall
        ));
        assert!(core::ptr::eq(heating.wall_bindings()[0].card, &card));
        for (owner, state) in [
            (&heating, &thermal),
            (&result.heating, &result.step.thermal),
        ] {
            let expected_wall = source(&card, state.temperature_k(WALL).unwrap(), loss);
            let actual_wall = owner.wall_materials().next().unwrap();
            assert_eq!(actual_wall.material().card_identity(), card.content_hash());
            assert_eq!(
                actual_wall.material().query_point(),
                expected_wall.material().query_point()
            );
            // Refresh receipts record explicit pinning; the independent oracle
            // selects the sole original claims. Compare their actual claims
            // and values, without equating these different selection receipts.
            for name in [
                "density",
                "young_modulus",
                "poisson_ratio",
                YOUNG_VISCOSITY_PROPERTY,
                "arm_modulus",
                "arm_time",
            ] {
                let actual = actual_wall.material().property(name).unwrap();
                let expected = expected_wall.material().property(name).unwrap();
                assert_eq!(actual.value_si(), expected.value_si());
                assert_eq!(
                    actual.answer().receipt.selected,
                    expected.answer().receipt.selected
                );
            }
            let states = gases(state);
            let reference = original(&geometry, &states, None, &expected_wall);
            let field = reference
                .prepared()
                .solve(&problem(false).inward_velocity, &mut || false)
                .unwrap();
            assert_eq!(owner.field().pressure_pa, field.pressure_pa);
            assert_eq!(
                owner.field().walls[0].normal_velocity_m_s,
                field.walls[0].normal_velocity_m_s
            );
            assert_eq!(
                owner.field().walls[0].relaxation_extension,
                field.walls[0].relaxation_extension
            );
            let numerical = reference.prepared().wall(0).unwrap();
            assert_eq!(owner.problem().mechanical_walls[0].wall, numerical.wall);
            assert_eq!(
                owner.problem().mechanical_walls[0].relaxation,
                numerical.relaxation
            );
            let heat = owner.wall_heating().next().unwrap();
            close(heat.area_m2, WALL_AREA, 1e-16);
            close(heat.mass_kg, thermal.cells()[&WALL].mass_kg, 1e-16);
            assert!(heat.absorbed_power_w() > 0.0);
            match loss {
                Loss::Kelvin => assert_eq!(heat.relaxation_dissipation_w, 0.0),
                Loss::Maxwell => {
                    assert_eq!(heat.resistance_dissipation_w, 0.0);
                    let WallMaterialLoss::Maxwell(memory) = actual_wall.loss() else {
                        panic!()
                    };
                    assert_eq!(
                        memory.source_identity(),
                        Some(actual_wall.material().identity())
                    );
                }
            }
        }
        assert_eq!(result.step.acoustic.density_dissipation_j, 0.0);
        assert_eq!(result.step.acoustic.compressional_dissipation_j, 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.cells()[&WALL].specific_h_j_kg, initial_wall_h);
        close(result.step.heat_report.residual_j, 0.0, 2e-11);
    }
}

#[test]
fn simultaneous_material_and_volume_refresh_uses_one_solve_and_accepted_start_maxwell_heat() {
    let (mesh, points) = fixture();
    let geometry = geometry(&mesh, &points);
    let caloric = caloric();
    let card = card(0.0, 400.0);
    let wall = source(&card, 300.0, Loss::Maxwell);
    let porous = laws();
    let thermal = initial(wall_events([0.0, 20.0 * CP], [0.12, 0.1], 0.0), &caloric);
    let heating = TemperatureInteriorHeating::try_new_with_material_walls(
        &geometry,
        &thermal,
        &mapping(&caloric, Some(&porous)),
        &bindings(&caloric, &card, &wall),
        fluid([280.0, 500.0]),
        problem(false),
        3,
        || Ok(()),
    )
    .unwrap();
    let result = heating
        .advance_to(&thermal, 2.0 * DT, config(DT), || Ok(()))
        .unwrap();
    assert_eq!(result.work.fluid_evaluations, 4);
    assert_eq!(result.work.porous_evaluations, 4);
    assert_eq!(result.work.wall_evaluations, 2);
    assert_eq!(result.work.field_solves, 2);
    let mut reference = thermal.clone();
    let mut joules = [0.0; 4];
    for step in 1..=2 {
        let wall = source(&card, reference.temperature_k(WALL).unwrap(), Loss::Maxwell);
        let states = gases(&reference);
        let model = original(&geometry, &states, Some(&porous), &wall);
        let field = model
            .prepared()
            .solve(&problem(false).inward_velocity, &mut || false)
            .unwrap();
        for (sum, power) in joules.iter_mut().zip([
            field.density_dissipation_w,
            field.compressional_dissipation_w,
            field.wall_accounting.resistance_dissipation_w,
            field.wall_accounting.relaxation_dissipation_w,
        ]) {
            *sum += power * DT;
        }
        let frozen = InteriorAcousticHeating::try_new_with_mechanical_walls(
            model.prepared(),
            &field.pressure_pa,
            &field.walls,
            &reference,
            &routes(&caloric),
            &wall_routes(&caloric),
            3,
            || Ok(()),
        )
        .unwrap();
        reference = frozen
            .advance_to(&reference, step as f64 * DT, config(DT), || Ok(()))
            .unwrap()
            .thermal;
    }
    for id in IDS.into_iter().chain([WALL]) {
        close(
            result.step.thermal.cells()[&id].specific_h_j_kg,
            reference.cells()[&id].specific_h_j_kg,
            2e-10,
        );
        assert_eq!(
            result.step.thermal.cells()[&id].mass_kg,
            thermal.cells()[&id].mass_kg
        );
        assert_eq!(
            result.step.thermal.cells()[&id].geometry,
            thermal.cells()[&id].geometry
        );
    }
    let heat = result.step.acoustic;
    for (actual, expected) in [
        heat.density_dissipation_j,
        heat.compressional_dissipation_j,
        heat.wall_resistance_dissipation_j,
        heat.wall_relaxation_dissipation_j,
    ]
    .into_iter()
    .zip(joules)
    {
        close(actual, expected, 2e-12);
    }
    assert!(heat.wall_relaxation_dissipation_j > 0.0);
    close(heat.deposited_j, joules.iter().sum(), 2e-12);
    close(
        result.step.heat_report.external_heat_j,
        heat.deposited_j + 0.27 * 2.0 * DT,
        2e-11,
    );
    close(result.step.heat_report.residual_j, 0.0, 2e-11);
}

#[test]
fn latent_temperature_preserves_material_snapshot_and_factor_while_wall_heat_accumulates() {
    let (mesh, points) = fixture();
    let geometry = geometry(&mesh, &points);
    let caloric = Plateau;
    let card = card(0.0, 400.0);
    let wall = source(&card, 300.0, Loss::Kelvin);
    let porous = laws();
    let thermal = initial(wall_events([0.0; 2], [0.12, 0.1], 0.0), &caloric);
    let gas = fluid([280.0, 500.0]);
    let calls = gas.calls.clone();
    let heating = TemperatureInteriorHeating::try_new_with_material_walls(
        &geometry,
        &thermal,
        &mapping(&caloric, Some(&porous)),
        &bindings(&caloric, &card, &wall),
        gas,
        problem(false),
        3,
        || Ok(()),
    )
    .unwrap();
    let result = heating
        .advance_to(&thermal, 2.0 * DT, config(DT), || Ok(()))
        .unwrap();
    assert_eq!(result.work, InteriorFeedbackWork::default());
    assert_eq!(calls.get(), 2);
    assert!(core::ptr::eq(heating.prepared(), result.heating.prepared()));
    assert!(core::ptr::eq(
        heating.wall_materials().next().unwrap(),
        result.heating.wall_materials().next().unwrap()
    ));
    assert!(result.step.thermal.cells()[&WALL].specific_h_j_kg > 0.0);
    assert!(result.step.acoustic.wall_resistance_dissipation_j > 0.0);
    assert_eq!(result.step.thermal.temperature_k(WALL).unwrap(), 300.0);
}

#[test]
fn original_wall_source_density_domain_and_cancellation_refusals_preserve_accepted_state() {
    let (mesh, points) = fixture();
    let geometry = geometry(&mesh, &points);
    let caloric = caloric();
    let thermal = initial(wall_events([0.0; 2], [0.0; 2], 0.0), &caloric);
    let initial_cells = thermal.cells().clone();
    let original_card = card(0.0, 400.0);
    let wall = source(&original_card, 300.0, Loss::Maxwell);
    let wrong_card = card(0.001, 400.0);
    assert!(matches!(
        TemperatureInteriorHeating::try_new_with_material_walls(
            &geometry,
            &thermal,
            &mapping(&caloric, None),
            &bindings(&caloric, &wrong_card, &wall),
            fluid([280.0, 500.0]),
            problem(false),
            3,
            || Ok(()),
        ),
        Err(InteriorFeedbackError::Wall(
            MaterialMechanicalHarmonicError::Wall(MaterialMechanicalError::Source {
                boundary: 0,
                ..
            })
        ))
    ));
    assert!(matches!(
        TemperatureInteriorHeating::try_new_with_material_walls(
            &geometry,
            &thermal,
            &mapping(&caloric, None),
            &bindings(&caloric, &original_card, &wall),
            fluid([280.0, 500.0]),
            problem(true),
            3,
            || Ok(()),
        ),
        Err(InteriorFeedbackError::Heating(
            InteriorHeatingError::Invalid(_)
        ))
    ));
    for invalid in [card(0.001, 400.0), card(0.0, 300.000001)] {
        let source = source(&invalid, 300.0, Loss::Maxwell);
        let heating = TemperatureInteriorHeating::try_new_with_material_walls(
            &geometry,
            &thermal,
            &mapping(&caloric, None),
            &bindings(&caloric, &invalid, &source),
            fluid([280.0, 500.0]),
            problem(false),
            3,
            || Ok(()),
        )
        .unwrap();
        let pressure = heating.field().pressure_pa.clone();
        let material = heating.wall_materials().next().unwrap().clone();
        assert!(matches!(
            heating.advance_to(&thermal, DT, config(DT), || Ok(())),
            Err(InteriorFeedbackError::Wall(
                MaterialMechanicalHarmonicError::Wall(MaterialMechanicalError::Source {
                    boundary: 0,
                    ..
                })
            ))
        ));
        assert_eq!(heating.field().pressure_pa, pressure);
        assert_eq!(heating.wall_materials().next().unwrap(), &material);
        assert_eq!(thermal.cells(), &initial_cells);
    }
    let heating = TemperatureInteriorHeating::try_new_with_material_walls(
        &geometry,
        &thermal,
        &mapping(&caloric, None),
        &bindings(&caloric, &original_card, &wall),
        fluid([280.0, 500.0]),
        problem(false),
        3,
        || Ok(()),
    )
    .unwrap();
    let material = heating.wall_materials().next().unwrap().clone();
    let mut total_polls = 0;
    let accepted = heating
        .advance_to(&thermal, 2.0 * DT, config(DT), || {
            total_polls += 1;
            Ok(())
        })
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        heating.advance_to(&thermal, 2.0 * DT, config(DT), || {
            polls += 1;
            if polls >= total_polls * 3 / 4 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(InteriorFeedbackError::Heating(
            InteriorHeatingError::Thermal(DepositionError::Cancelled)
        ))
    ));
    assert_eq!(heating.wall_materials().next().unwrap(), &material);
    assert_eq!(thermal.cells(), &initial_cells);
    assert_eq!(thermal.time_s(), 0.0);
    let retry = heating
        .advance_to(&thermal, 2.0 * DT, config(DT), || Ok(()))
        .unwrap();
    assert_eq!(retry.step.thermal.cells(), accepted.step.thermal.cells());
    assert_eq!(
        retry.heating.field().pressure_pa,
        accepted.heating.field().pressure_pa
    );
    assert_eq!(retry.step.acoustic, accepted.step.acoustic);
    assert_eq!(retry.work, accepted.work);
}
