//! Source films act on the original moving cavity and explicit thermal history.
use super::*;
use fs_couple::air_path::input::convection::{
    ForcedGasBoundary, ForcedGasGeometry, GasConvectionBoundary, NaturalConvectionGeometry,
    NaturalGasBoundary,
};
use fs_couple::solid_enthalpy::gas::{SolidGasConvection, SolidGasError};
use fs_couple::solid_enthalpy::thermoelastic::dynamic::cavity::CavityError;
use fs_material::gas_species::input::temperature::IsobaricGasLaw;

fn sphere_gas(t: f64) -> ClosedGasParcel {
    // Pr=0.7185 satisfies the original sphere card's 0.7 floor. Cv=2.5 R
    // preserves the existing independent constant-capacity cavity oracle.
    let (species, transport) = gas_fixture::component("moving-sphere-fixture", 0.027, 0.0);
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![transport],
        None,
    )
    .unwrap();
    ClosedGasParcel::new(model, gas(t, 0.0).point().clone(), GAS_VOLUME).unwrap()
}
fn patches(natural: bool, distance: f64, tolerance: f64) -> [SolidGasContact; 1] {
    contacts(PATCH_AREA).map(|mut contact| {
        contact.geometry.gas_distance_m = 0.0;
        contact.geometry.cell_distance_m = distance;
        contact.convection = Some(SolidGasConvection {
            boundary: if natural {
                GasConvectionBoundary::Natural(NaturalGasBoundary {
                    geometry: NaturalConvectionGeometry::Sphere,
                    characteristic_length_m: 0.4,
                    exposed_area_m2: PATCH_AREA,
                    gravity_m_s2: 0.0,
                })
            } else {
                GasConvectionBoundary::Forced(ForcedGasBoundary {
                    geometry: ForcedGasGeometry::CylinderCrossflow,
                    characteristic_length_m: 0.03,
                    exposed_area_m2: PATCH_AREA,
                    speed_m_s: 3.0,
                })
            },
            maximum_surface_iterations: 96,
            surface_power_tolerance_w: tolerance,
        });
        contact
    })
}
fn launch_film<'s, 'm>(
    source: &'s SolidEnthalpy<'m>,
    thermal: &DepositionSimulation<'_>,
    gas: ClosedGasParcel,
    speed: f64,
    patches: &[SolidGasContact],
    cx: &fs_exec::Cx<'_>,
) -> Result<ThermoelasticCavity<'s, 'm>, CavityError> {
    let face = source
        .geometry()
        .faces()
        .iter()
        .position(|face| {
            face.nodes.contains(&1) && face.nodes.contains(&2) && face.nodes.contains(&3)
        })
        .unwrap();
    let faces = [face];
    let reference = [(3, 0.3)];
    let loads = source
        .geometry()
        .assemble(
            &[SurfaceLoad3::Pressure {
                faces: &faces,
                pascals: gas.state().state().pressure - PRESSURE,
            }],
            VolumeLoad3::None,
            &reference,
            cx,
        )
        .unwrap();
    let body = ThermoelasticSolid::try_new(
        source,
        thermal,
        &loads,
        &fixed(),
        &[[280.0, 360.0]],
        mechanical(),
        solve(),
        cx,
    )
    .unwrap();
    body.begin_cavity(
        thermal,
        gas,
        &faces,
        PRESSURE,
        &reference,
        patches,
        &velocity(speed),
        dynamic_controls(),
        cx,
    )
}

#[test]
fn source_film_series_limit_preserves_native_wall_motion_and_pressure_work() {
    with_cx(|cx| {
        let (card, curve) = wall(10000.0, 0.0, 0.0, 6.0, 1.0);
        let thermal = heat(&curve, 6.0, false);
        let context = context(300.0);
        let source = SolidEnthalpy::try_new(
            &thermal,
            &NODES,
            &TETS,
            vec![source(30, &card, &curve, &context)],
            TetAssemblyBudget::standard(),
            0.05,
            mechanical().mesh_frame,
            cx,
        )
        .unwrap();
        for (speed, gas_t) in [(0.02, 340.0), (-0.02, 290.0)] {
            for distance in [0.0, 0.02] {
                let patches = patches(true, distance, 1e-10);
                let cavity =
                    launch_film(&source, &thermal, sphere_gas(gas_t), speed, &patches, cx).unwrap();
                let before = cavity.solid().nodal_motion(cx).unwrap();
                let original_gas = cavity.gas().clone();
                // Native Churchill sphere at Ra=0 gives h=2*k/D=0.15.
                let g = PATCH_AREA / (distance / 1.0 + 1.0 / 0.15);
                let expected =
                    reference_for_conductance(speed, 0.3, 10000.0, 6.0, g, gas_t, 0.005, 0.01);
                let next = cavity
                    .advance(&thermal, &(), 0.005, 0.01, controls(), [0.0; 2], cx)
                    .unwrap();
                let motion = next.cavity.solid().nodal_motion(cx).unwrap();
                near(
                    next.cavity.gas().state().state().temperature,
                    expected[0],
                    1e-7,
                );
                near(motion[0][1][0], expected[1], 2e-11);
                near(motion[1][1][0], expected[2], 1e-8);
                near(
                    next.thermal.temperature_k(IDS[0]).unwrap(),
                    expected[3],
                    1e-8,
                );
                near(
                    next.cavity.gas().state().state().pressure,
                    expected[4],
                    1e-6,
                );
                let row = &next.balance.contacts[0];
                let film = row.convection.as_ref().unwrap();
                near(row.last_cell_temperature_k, 300.0, 1e-12);
                near(row.last_power_into_solid_w, g * (expected[0] - 300.0), 1e-9);
                near(
                    row.heat_into_solid_j,
                    0.005 * g * (expected[0] - 300.0),
                    1e-10,
                );
                near(row.last_explicit_decay_conductance_w_k, g, 1e-12);
                near(
                    film.film.nusselt.as_ref().unwrap().evidence().value,
                    2.0,
                    0.0,
                );
                near(
                    film.film.wall_temperature_k,
                    expected[0] - row.last_power_into_solid_w / (0.15 * PATCH_AREA),
                    1e-7,
                );
                assert_eq!(
                    row.last_power_into_solid_w.is_sign_positive(),
                    gas_t > 300.0
                );
                assert!(film.surface_balance_residual_w.abs() <= 1e-10);
                near(
                    next.balance.gas_to_solid_heat_j,
                    row.heat_into_solid_j,
                    1e-11,
                );
                let swept = -(expected[1] - before[0][1][0]) / 6.0;
                near(
                    next.balance.gas_pressure_work_j,
                    0.5 * (PRESSURE + expected[4]) * swept,
                    2e-9,
                );
                assert!(
                    next.balance.pressure_work_residual_j.abs() <= controls().energy_tolerance_j
                );
                assert!(next.balance.combined_residual_j.abs() <= controls().energy_tolerance_j);
                assert_eq!(next.cavity.gas().mass_kg(), original_gas.mass_kg());
                assert_eq!(*cavity.gas(), original_gas);
                assert_eq!(cavity.solid().nodal_motion(cx).unwrap(), before);
                assert_eq!(thermal.time_s(), 0.0);
            }
        }
    });
}

#[test]
fn nonlinear_film_tracks_applied_substep_surface_and_uses_actual_flux_decay() {
    with_cx(|cx| {
        let (card, curve) = wall(10000.0, 0.0, 0.0, 6.0, 1.0);
        let mut thermal = heat(&curve, 6.0, true);
        let context = context(300.0);
        let source = SolidEnthalpy::try_new(
            &thermal,
            &NODES,
            &TETS,
            vec![source(30, &card, &curve, &context)],
            TetAssemblyBudget::standard(),
            0.05,
            mechanical().mesh_frame,
            cx,
        )
        .unwrap();
        // Deliberately finite surface allowance makes diagnostic G differ from
        // actual Q/DeltaT; the explicit decay must use the transferred flux.
        let patches = patches(false, 0.02, 0.01);
        let mut cavity =
            launch_film(&source, &thermal, gas(340.0, 0.0), 0.02, &patches, cx).unwrap();
        let mut observer = ReservoirObserver {
            steps: 0,
            reservoir_calls: 0,
            fail_at: None,
        };
        let mut config = controls();
        config.thermal.max_dt_s = 0.001;
        let exhausted = CavitySolve {
            maximum_trials: 1,
            ..config
        };
        assert!(matches!(
            cavity.advance(&thermal, &observer, 0.005, 0.002, exhausted, [0.0; 2], cx),
            Err(CavityError::TrialBudget)
        ));
        assert_eq!(cavity.solid().step(), 0);
        assert_eq!(observer.steps, 0);
        for step in 1..=2 {
            let next = cavity
                .advance(
                    &thermal,
                    &observer,
                    0.005 * step as f64,
                    0.002,
                    config,
                    [0.0; 2],
                    cx,
                )
                .unwrap();
            assert!(next.observer.steps >= observer.steps + 5);
            let row = &next.balance.contacts[0];
            let film = row.convection.as_ref().unwrap();
            let bulk = next.cavity.gas().state().state().temperature;
            let cell = row.last_cell_temperature_k;
            let surface = film.film.wall_temperature_k;
            assert!(cell < surface && surface < bulk);
            // The native harmonic decay budget can leave a short final tail;
            // its applied source still belongs to the pre-tail cell state.
            assert!(next.thermal.temperature_k(IDS[0]).unwrap() > cell + 1e-6);
            near(row.source.query().axes()["T"], cell, 0.0);
            assert_eq!(row.source.card_identity(), card.content_hash());
            let law = IsobaricGasLaw::try_from_model(
                next.cavity.gas().model().clone(),
                next.cavity.gas().point().clone(),
            )
            .unwrap();
            let fresh = patches[0]
                .convection
                .unwrap()
                .boundary
                .evaluate(&law, surface)
                .unwrap();
            assert_eq!(film.film.property_gas, fresh.property_gas);
            assert_eq!(film.film.nusselt, fresh.nusselt);
            near(row.last_power_into_solid_w, fresh.heat_into_wall_w, 0.0);
            near(
                row.last_explicit_decay_conductance_w_k,
                row.last_power_into_solid_w / (bulk - cell),
                0.0,
            );
            assert!(
                (row.last_explicit_decay_conductance_w_k - row.last_conductance_w_k).abs() > 1e-8
            );
            assert!(film.surface_balance_residual_w.abs() > 1e-5);
            assert!(film.surface_balance_residual_w.abs() <= 0.01);
            near(
                film.film.property_gas.state().temperature,
                0.5 * (surface + bulk),
                1e-12,
            );
            near(
                film.film.property_gas.state().pressure,
                next.cavity.gas().state().state().pressure,
                1e-9,
            );
            assert!(
                (next.cavity.gas().state().state().pressure
                    - cavity.gas().state().state().pressure)
                    .abs()
                    > 1e-3
            );
            near(
                next.balance.gas_to_solid_heat_j,
                row.heat_into_solid_j,
                1e-10,
            );
            assert!(next.balance.external_solid_heat_j > 0.09);
            assert!(next.balance.gas_energy_residual_j.abs() <= config.energy_tolerance_j);
            assert!(next.balance.combined_residual_j.abs() <= config.energy_tolerance_j);
            assert_eq!(next.cavity.gas().point().axes().get("lot"), Some(&42.0));
            cavity = next.cavity;
            thermal = next.thermal;
            observer = next.observer;
        }
    });
}

fn limited_wall() -> (MaterialCard, EquilibriumEnthalpyPhaseCurve) {
    let (original, _) = wall(10000.0, 0.0, 0.0, 6.0, 1.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "thermal_conductivity" {
            claim.validity =
                ValidityDomain::unconstrained().with_quantity("T", quantity(), 280.0, 300.2);
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        [280.0, 500.0]
            .map(|temperature_k| EnthalpyPhaseKnot {
                temperature_k,
                specific_enthalpy_j_kg: CP * (temperature_k - 280.0),
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 6.0,
            })
            .to_vec(),
    )
    .unwrap();
    (card, curve)
}

#[test]
fn surface_budget_and_later_source_refusal_publish_no_dynamic_or_thermal_history() {
    with_cx(|cx| {
        let (card, curve) = limited_wall();
        let thermal = heat(&curve, 6.0, true);
        let original_cells = thermal.cells().clone();
        let context = context(300.0);
        let source = SolidEnthalpy::try_new(
            &thermal,
            &NODES,
            &TETS,
            vec![source(30, &card, &curve, &context)],
            TetAssemblyBudget::standard(),
            0.05,
            mechanical().mesh_frame,
            cx,
        )
        .unwrap();
        let patches = patches(false, 0.02, 1e-10);
        let mut bad = patches.clone();
        bad[0]
            .convection
            .as_mut()
            .unwrap()
            .maximum_surface_iterations = 1;
        bad[0]
            .convection
            .as_mut()
            .unwrap()
            .surface_power_tolerance_w = 0.0;
        assert!(
            matches!(launch_film(&source, &thermal, gas(340.0, 0.0), 0.02, &bad, cx),
            Err(CavityError::Contact(SolidGasError::SurfaceBalance { residual_w, iterations: 1 }))
                if residual_w.is_finite() && residual_w != 0.0)
        );
        let cavity = launch_film(&source, &thermal, gas(340.0, 0.0), 0.02, &patches, cx).unwrap();
        let gas_before = cavity.gas().clone();
        let motion_before = cavity.solid().nodal_motion(cx).unwrap();
        let observer = ReservoirObserver {
            steps: 0,
            reservoir_calls: 0,
            fail_at: None,
        };
        let mut config = controls();
        config.thermal.max_dt_s = 0.001;
        // The film is valid initially, then the native reservoir raises the
        // actual wall substep above the conductivity source's 300.2 K ceiling.
        assert!(matches!(
            cavity.advance(&thermal, &observer, 0.01, 0.0, config, [0.0; 2], cx),
            Err(CavityError::Layer(_))
        ));
        assert_eq!(*cavity.gas(), gas_before);
        assert_eq!(cavity.solid().nodal_motion(cx).unwrap(), motion_before);
        assert_eq!(cavity.solid().step(), 0);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.cells(), &original_cells);
        assert_eq!(observer.steps, 0);
        assert_eq!(observer.reservoir_calls, 0);
    });
}
