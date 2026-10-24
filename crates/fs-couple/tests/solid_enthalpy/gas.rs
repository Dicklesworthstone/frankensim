//! G1/G3/G4: one finite gas, two spatial thermal cells and actual solid motion.
use super::*;
use fs_couple::solid_enthalpy::gas::*;
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{GasEnergySolve, IsochoricGasParcel},
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_solver::NewtonKrylovConfig;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis};

// Reuse the existing manufactured NASA9/transport fixture, not a new gas law.
#[allow(dead_code)]
#[path = "../material_exchange/gas.rs"]
mod gas_fixture;

const PAIR_NODES: [[f64; 3]; 5] = [
    [0., 0., 0.],
    [1., 0., 0.],
    [0., 1., 0.],
    [0., 0., 1.],
    [0., 0., -1.],
];
const PAIR_TETS: [[usize; 4]; 2] = [[0, 1, 2, 3], [0, 2, 1, 4]];

fn gas(t: f64, volume: f64) -> IsochoricGasParcel {
    let (species, transport) = gas_fixture::component("spatial-fixture", 0.028, 0.0);
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
    let typed = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let point = QueryPoint::new()
        .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), 100000.0)
        .unwrap()
        .with("material-lot", 7.0)
        .unwrap();
    IsochoricGasParcel::new(model, point, volume).unwrap()
}

fn transport_material(k: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let mut claims = ClaimSet::new();
    for (key, value) in [
        (heat.clone(), 10.0),
        (density.clone(), 6.0),
        (PropertyKey::new("young_modulus", Pressure::DIMS), 1000.0),
        (PropertyKey::new("poisson_ratio", Dims([0; 6])), 0.0),
        (
            PropertyKey::new(
                "linear_thermal_expansion_coefficient",
                INVERSE_TEMPERATURE_DIMS,
            ),
            1e-5,
        ),
        (
            PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            k,
        ),
    ] {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        280.,
                        600.,
                    )
                    .with("lot", 7., 7.),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured spatial gas/solid test".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored specimen".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = caloric_for(&card);
    (card, caloric)
}

fn caloric_for(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(280.),
            upper: &point(600.),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-10,
    )
    .unwrap()
}

fn initial_pair(
    curve: &EquilibriumEnthalpyPhaseCurve,
    temperature: [f64; 2],
    power: [f64; 2],
) -> DepositionSimulation<'_> {
    let mut events = Vec::new();
    for (i, id) in [30, 10].into_iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: 77,
                material_id: MaterialId(0),
                mass_kg: 1.0,
                specific_h_j_kg: 10.0 * (temperature[i] - 280.0),
                geometry: CellGeometry {
                    start_m: [0., 0., i as f64],
                    end_m: [1., 0., i as f64],
                    volume_m3: 1.0 / 6.0,
                    surface_area_m2: 3.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: power[i],
                },
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: CellId(30),
            b: CellId(10),
            area_m2: 0.5,
            conductance_w_k: 0.5,
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), curve as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}

fn contacts(areas: [f64; 2]) -> [SolidGasContact; 2] {
    [(99, 30, areas[0]), (8, 10, areas[1])].map(|(id, cell, area)| SolidGasContact {
        id,
        cell: CellId(cell),
        geometry: GasConductionGeometry {
            area_m2: area,
            gas_distance_m: 0.001,
            cell_distance_m: 0.02,
        },
        conductivity: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        selection: MaterialPropertySelection::SingleClaimOnly,
    })
}
fn solve() -> SolidGasSolve {
    SolidGasSolve {
        maximum_cells: 2,
        maximum_contacts: 4,
        maximum_trials: 500,
        maximum_newton_iterations: 20,
        difference_step_j: 1e-3,
        cell_energy_tolerance_j: 1e-7,
        thermal_energy_tolerance_j: 1e-7,
        gas_energy_tolerance_j: 1e-7,
        total_energy_tolerance_j: 2e-7,
        maximum_gas_temperature_increment_k: 200.0,
        gas: GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1800.0,
            energy_tolerance_j_kg: 1e-9,
            maximum_iterations: 128,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-8,
            relative_tolerance: 1e-14,
            forcing_minimum: 1e-4,
            linear_restart: 4,
            ..NewtonKrylovConfig::default()
        },
    }
}

// Independent elimination of three constant-capacity energy equations. The
// production path has no Cp or hand-written three-reservoir formula.
fn reference(
    t: [f64; 3],
    capacities: [f64; 3],
    g: [f64; 2],
    dt: f64,
    external: [f64; 3],
) -> [f64; 3] {
    let h = 0.5;
    let [c1, c2, cg] = capacities.map(|c| c / dt);
    let [g1, g2] = g;
    let a = c1 + h + g1;
    let b = c2 + h + g2;
    let determinant = a * b - h * h;
    let r1 = c1 * t[0] + external[0] / dt;
    let r2 = c2 * t[1] + external[1] / dt;
    let a1 = (b * r1 + h * r2) / determinant;
    let a2 = (h * r1 + a * r2) / determinant;
    let b1 = (b * g1 + h * g2) / determinant;
    let b2 = (h * g1 + a * g2) / determinant;
    let tg =
        (cg * t[2] + external[2] / dt + g1 * a1 + g2 * a2) / (cg + g1 + g2 - g1 * b1 - g2 * b2);
    [a1 + b1 * tg, a2 + b2 * tg, tg]
}

#[test]
fn finite_gas_and_material_conductivity_drive_nonuniform_motion_with_one_energy_account() {
    with_cx(|cx| {
        let mut observed = Vec::new();
        for (k, volume) in [(1.0, 0.05), (2.0, 0.05), (1.0, 0.2)] {
            let (card, caloric) = transport_material(k);
            let thermal = initial_pair(caloric.curve(), [340., 360.], [8., 0.]);
            let gas = gas(450., volume);
            let saved = gas.clone();
            let context = point(340.);
            let source = SolidEnthalpy::try_new(
                &thermal,
                &PAIR_NODES,
                &PAIR_TETS,
                vec![
                    super::source(30, &card, caloric.curve(), &context),
                    super::source(10, &card, caloric.curve(), &context),
                ],
                TetAssemblyBudget::standard(),
                0.1,
                config().mesh_frame,
                cx,
            )
            .unwrap();
            let next = source
                .step_with_isochoric_gas(
                    &thermal,
                    &gas,
                    2.,
                    4.,
                    &contacts([0.1, 0.2]),
                    solve(),
                    &[0.; 3],
                    cx,
                )
                .unwrap();
            let capacity_gas = 2.5 * 100000.0 * volume / 450.0; // PV=nRT, cv=2.5R in fixture.
            let g = [0.1, 0.2].map(|a| a / (0.001 / 0.03 + 0.02 / k));
            let expected = reference(
                [340., 360., 450.],
                [10., 10., capacity_gas],
                g,
                2.,
                [16., 0., 4.],
            );
            near(
                next.thermal().temperature_k(CellId(30)).unwrap(),
                expected[0],
                2e-8,
            );
            near(
                next.thermal().temperature_k(CellId(10)).unwrap(),
                expected[1],
                2e-8,
            );
            near(next.gas().state().state().temperature, expected[2], 2e-8);
            near(
                next.gas().state().state().pressure,
                100000.0 * expected[2] / 450.0,
                1e-5,
            );
            near(next.report().external_thermal_heat_j, 16., 1e-12);
            near(next.report().external_gas_heat_j, 4., 1e-12);
            assert!(next.report().combined_energy_residual_j.abs() < 2e-7);
            assert!(next.report().maximum_cell_residual_j < 1e-7);
            assert_eq!(
                next.report()
                    .contacts
                    .iter()
                    .map(|p| p.id)
                    .collect::<Vec<_>>(),
                vec![8, 99]
            );
            for sample in &next.report().contacts {
                assert_eq!(sample.source.card_identity(), card.content_hash());
                near(
                    sample
                        .source
                        .property("thermal_conductivity")
                        .unwrap()
                        .value_si(),
                    k,
                    1e-12,
                );
                let i = if sample.cell == CellId(30) { 0 } else { 1 };
                near(
                    sample.solid_to_gas_heat_j,
                    2. * g[i] * (expected[i] - expected[2]),
                    2e-7,
                );
            }
            let loads = source
                .geometry()
                .assemble(&[], VolumeLoad3::None, &[], cx)
                .unwrap();
            let equilibrium = next.state().solve(&loads, &FREE, config(), cx).unwrap();
            let strains = [1e-5 * (expected[0] - 340.), 1e-5 * (expected[1] - 360.)];
            let displacement = equilibrium.solution.displacement.displacement_m();
            near(displacement[1][0], 0.5 * (strains[0] + strains[1]), 1e-10);
            near(displacement[2][1], 0.5 * (strains[0] + strains[1]), 1e-10);
            near(displacement[3][2], strains[0], 1e-10);
            near(displacement[4][2], -strains[1], 1e-10);
            assert_eq!(next.thermal().mass_kg(), 2.0);
            assert_eq!(next.gas().mass_kg(), gas.mass_kg());
            assert_eq!(next.gas().volume_m3(), gas.volume_m3());
            assert_eq!(next.gas().point().axes().get("material-lot"), Some(&7.0));
            assert_eq!(gas, saved);
            assert_eq!(thermal.time_s(), 0.0);
            observed.push(expected);
        }
        assert!((observed[0][0] - observed[1][0]).abs() > 1.0);
        assert!((observed[0][0] - observed[2][0]).abs() > 1.0);
    });
}

#[test]
fn balanced_hot_and_cold_patches_evolve_even_when_net_gas_heat_is_zero() {
    with_cx(|cx| {
        let (card, caloric) = transport_material(1.0);
        let thermal = initial_pair(caloric.curve(), [320., 480.], [0., 0.]);
        let gas = gas(400., 0.1);
        let context = point(320.);
        let source = SolidEnthalpy::try_new(
            &thermal,
            &PAIR_NODES,
            &PAIR_TETS,
            vec![
                super::source(30, &card, caloric.curve(), &context),
                super::source(10, &card, caloric.curve(), &context),
            ],
            TetAssemblyBudget::standard(),
            0.1,
            config().mesh_frame,
            cx,
        )
        .unwrap();
        let mut patches = contacts([0.1, 0.1]);
        let next = source
            .step_with_isochoric_gas(&thermal, &gas, 1., 0., &patches, solve(), &[0.; 3], cx)
            .unwrap();
        let expected_gap = 80.0 / (1.0 + (1.875 + 2.0 * 0.5) / 10.0);
        near(
            next.thermal().temperature_k(CellId(30)).unwrap(),
            400. - expected_gap,
            1e-8,
        );
        near(
            next.thermal().temperature_k(CellId(10)).unwrap(),
            400. + expected_gap,
            1e-8,
        );
        near(next.gas().state().state().temperature, 400., 1e-8);
        near(next.report().solid_to_gas_heat_j, 0., 1e-7);
        near(next.report().gas.stored_energy_change_j, 0., 1e-7);
        assert!(next.report().contacts[0].solid_to_gas_heat_j > 100.0);
        assert!(next.report().contacts[1].solid_to_gas_heat_j < -100.0);
        patches.reverse();
        let reordered = source
            .step_with_isochoric_gas(&thermal, &gas, 1., 0., &patches, solve(), &[0.; 3], cx)
            .unwrap();
        assert_eq!(next.gas(), reordered.gas());
        assert_eq!(next.thermal().cells(), reordered.thermal().cells());
        assert_eq!(
            next.report().solid_to_gas_heat_j,
            reordered.report().solid_to_gas_heat_j
        );
    });
}

#[test]
fn later_intervals_cannot_reselect_a_different_conductivity_claim() {
    with_cx(|cx| {
        let (original, _) = transport_material(1.0);
        let mut claims = ClaimSet::new();
        for (_, claim) in original.claims().claims_ordered() {
            if claim.key.name() != "thermal_conductivity" {
                claims.insert_claim(claim.clone()).unwrap();
            }
        }
        let original_conductivity = original.claims_for("thermal_conductivity")[0].1;
        for (lower, upper, conductivity) in [(280.0, 399.0, 1.0), (400.0, 600.0, 2.0)] {
            let mut claim = original_conductivity.clone();
            claim.validity = ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    lower,
                    upper,
                )
                .with("lot", 7.0, 7.0);
            claim.value = PropertyValue::Scalar {
                value: conductivity,
                dims: THERMAL_CONDUCTIVITY_DIMS,
            };
            claims.insert_claim(claim).unwrap();
        }
        let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
        let caloric = caloric_for(&card);
        let thermal = initial_pair(caloric.curve(), [340.0, 360.0], [600.0, 600.0]);
        let context = point(340.0);
        let source = SolidEnthalpy::try_new(
            &thermal,
            &PAIR_NODES,
            &PAIR_TETS,
            vec![
                super::source(30, &card, caloric.curve(), &context),
                super::source(10, &card, caloric.curve(), &context),
            ],
            TetAssemblyBudget::standard(),
            0.1,
            config().mesh_frame,
            cx,
        )
        .unwrap();
        let warm = thermal
            .advance_to(1.0, thermal_config(), || Ok(()))
            .unwrap()
            .0;
        assert!(warm.temperature_k(CellId(30)).unwrap() >= 400.0);
        let gas = gas(450.0, 0.1);
        // An independent current-state query CAN select the other claim. The
        // original specimen may not silently acquire it for its next interval.
        let independent =
            fs_couple::cavity_phs::piston::thermal::layers::MaterialConductionLayers {
                material: &card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
                conductivity: &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
                selection: MaterialPropertySelection::SingleClaimOnly,
                area_m2: 0.1,
                gas_distance_m: 0.001,
                body_distance_m: 0.02,
            }
            .resolve(
                gas.state(),
                caloric
                    .curve()
                    .state_at_specific_enthalpy(warm.cells()[&CellId(30)].specific_h_j_kg)
                    .unwrap(),
            )
            .unwrap();
        near(
            independent
                .source
                .property("thermal_conductivity")
                .unwrap()
                .value_si(),
            2.0,
            1e-12,
        );
        assert!(matches!(
            source.step_with_isochoric_gas(
                &warm,
                &gas,
                2.0,
                0.0,
                &contacts([0.1, 0.2]),
                solve(),
                &[0.0; 3],
                cx
            ),
            Err(SolidGasError::Source(_))
        ));
        assert_eq!(warm.time_s(), 1.0);
        near(gas.state().state().temperature, 450.0, 1e-12);
    });
}

#[test]
fn failed_source_surface_budget_nonlinear_work_and_cancellation_keep_both_inventories() {
    with_cx(|cx| {
        let (card, caloric) = transport_material(1.0);
        let thermal = initial_pair(caloric.curve(), [340., 360.], [0., 0.]);
        let gas = gas(450., 0.1);
        let saved = gas.clone();
        let original_cells = thermal.cells().clone();
        let context = point(340.);
        let source = SolidEnthalpy::try_new(
            &thermal,
            &PAIR_NODES,
            &PAIR_TETS,
            vec![
                super::source(30, &card, caloric.curve(), &context),
                super::source(10, &card, caloric.curve(), &context),
            ],
            TetAssemblyBudget::standard(),
            0.1,
            config().mesh_frame,
            cx,
        )
        .unwrap();
        let patches = contacts([0.1, 0.2]);
        let limited = SolidGasSolve {
            maximum_trials: 1,
            ..solve()
        };
        assert!(matches!(
            source.step_with_isochoric_gas(&thermal, &gas, 1., 0., &patches, limited, &[0.; 3], cx),
            Err(SolidGasError::TrialBudget)
        ));
        let mut missing = patches.clone();
        missing[0].conductivity =
            PropertyKey::new("absent-conductivity", THERMAL_CONDUCTIVITY_DIMS);
        assert!(
            source
                .step_with_isochoric_gas(&thermal, &gas, 1., 0., &missing, solve(), &[0.; 3], cx)
                .is_err()
        );
        let overlapping = [
            patches[0].clone(),
            SolidGasContact {
                id: 100,
                geometry: GasConductionGeometry {
                    area_m2: 2.5,
                    ..patches[0].geometry
                },
                ..patches[0].clone()
            },
        ];
        assert!(
            source
                .step_with_isochoric_gas(
                    &thermal,
                    &gas,
                    1.,
                    0.,
                    &overlapping,
                    solve(),
                    &[0.; 3],
                    cx
                )
                .is_err()
        );
        assert!(
            source
                .step_with_isochoric_gas(&thermal, &gas, 1., 1e9, &patches, solve(), &[0.; 3], cx)
                .is_err()
        );
        let gate = CancelGate::new();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 7,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(matches!(
                source.step_with_isochoric_gas(
                    &thermal,
                    &gas,
                    1.,
                    0.,
                    &patches,
                    solve(),
                    &[0.; 3],
                    &cancelled
                ),
                Err(SolidGasError::Cancelled)
            ));
        });
        assert_eq!(gas, saved);
        assert_eq!(thermal.cells(), &original_cells);
        assert_eq!(thermal.time_s(), 0.0);
    });
}
