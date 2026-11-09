//! G0/G1/G3/G4: actual finite bodies, one heated molecular snapshot and atomic closure.
//! Synthetic calorics and declared h are analytic fixtures, not measured source laws.

mod trajectory;

use super::*;
use crate::{
    cavity_phs::piston::thermal::HeatExchangeSolve,
    finite_gas_thermal::{
        ConvectionCoefficient, ConvectionError,
        bodies::{FiniteGasBodiesError, FiniteGasBody},
    },
    molecular_exchange::network::thermal::{
        self, MolecularGasThermalBoundary, MolecularNetworkThermalCandidate,
        MolecularNetworkThermalControls, MolecularNetworkThermalError,
    },
};
use fs_blake3::ContentHash;
use fs_conduction::lumped::{LumpedThermalEnvironment, LumpedThermalTransport};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase, UniformEnthalpyStepInput,
};

const WALL_MASS: f64 = 0.005;
const WALL_CP: f64 = 1000.0;
const WALL_AREA: f64 = 0.01;

fn curve(densities: [f64; 2]) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([61; 32]),
        SolidLiquidPhase::Solid,
        [200.0, 1800.0]
            .into_iter()
            .zip(densities)
            .map(|(temperature_k, density)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: WALL_CP * (temperature_k - 200.0),
                temperature_k,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: density,
            })
            .collect(),
    )
    .unwrap()
}
fn wall(curve: &EquilibriumEnthalpyPhaseCurve, temperature: f64) -> UniformEnthalpyStepInput<'_> {
    let initial = curve
        .state_at_specific_enthalpy(WALL_CP * (temperature - 200.0))
        .unwrap();
    UniformEnthalpyStepInput {
        curve,
        initial,
        mass_kg: WALL_MASS,
        volume_m3: WALL_MASS / initial.bulk_density_kg_m3(),
        surface_area_m2: WALL_AREA,
        internal_heat_j: 0.0,
        duration_s: DT,
    }
}
fn body(h: f64) -> LumpedThermalEnvironment {
    LumpedThermalEnvironment {
        temperature_k: T,
        radiation_temperature_k: T,
        convection_w_per_m2_k: h,
        transport: LumpedThermalTransport::try_declared(100.0, 0.0).unwrap(),
        enthalpy_tolerance_j_kg: 1e-10,
        maximum_thermal_residual_j: 1e-9,
    }
}
fn exchange() -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-50.0, 50.0],
        source_tolerance_j: 1e-9,
        energy_tolerance_j: 1e-7,
        maximum_iterations: 128,
    }
}
fn thermal_controls(frame: [f64; 3]) -> MolecularNetworkThermalControls {
    MolecularNetworkThermalControls {
        exchange: budgets(frame),
        maximum_bodies: 16,
        energy_tolerance_j: 1e-6,
        momentum_tolerance_n_s: 1e-11,
    }
}
fn boundary<'a>(gas: u64, bodies: &'a [FiniteGasBody<'a>]) -> MolecularGasThermalBoundary<'a> {
    MolecularGasThermalBoundary {
        gas,
        bodies,
        gas_heat_j: 0.0,
        exchange: exchange(),
        maximum_gas_temperature_increment_k: 200.0,
    }
}
fn scope(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 0x4d4e54,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

// Independent linear-system solution: all walls share one gas endpoint,
// including each distinct wall's internal heat. Valid for PURE and these charts.
fn analytic(
    gas: &MovingGasParcel,
    group: MolecularGasThermalBoundary<'_>,
) -> (f64, BTreeMap<u64, f64>) {
    let cg = gas.parcel().mass_kg() * 1.5 * R / M;
    let mut numerator = cg * temperature(gas) + group.gas_heat_j;
    let mut denominator = cg;
    for wall in group.bodies {
        let c = wall.input.mass_kg * WALL_CP;
        let k = DT * wall.environment.convection_w_per_m2_k * wall.input.surface_area_m2;
        numerator +=
            k * (c * wall.input.initial.temperature_k() + wall.input.internal_heat_j) / (c + k);
        denominator += k * c / (c + k);
    }
    let gas_t = numerator / denominator;
    let walls = group
        .bodies
        .iter()
        .map(|wall| {
            let c = wall.input.mass_kg * WALL_CP;
            let k = DT * wall.environment.convection_w_per_m2_k * wall.input.surface_area_m2;
            (
                wall.id,
                (c * wall.input.initial.temperature_k() + wall.input.internal_heat_j + k * gas_t)
                    / (c + k),
            )
        })
        .collect();
    (gas_t, walls)
}

fn check_total<S>(
    original: &BTreeMap<u64, MovingGasParcel>,
    groups: &[MolecularGasThermalBoundary<'_>],
    candidate: &MolecularNetworkThermalCandidate<S>,
) {
    // These checks intentionally use modest, zero-datum absolute calorics.
    let gas_energy: f64 = original
        .iter()
        .map(|(id, gas)| {
            let next = &candidate.gases()[id];
            internal(next) - internal(gas) + next.kinetic_energy_j() - gas.kinetic_energy_j()
        })
        .sum();
    let wall_energy: f64 = groups
        .iter()
        .flat_map(|g| g.bodies)
        .map(|body| {
            body.input.mass_kg
                * (candidate.bodies()[&body.id].specific_enthalpy_j_kg()
                    - body.input.initial.specific_enthalpy_j_kg())
        })
        .sum();
    let authored: f64 = groups
        .iter()
        .map(|group| {
            group.gas_heat_j
                + group
                    .bodies
                    .iter()
                    .map(|body| body.input.internal_heat_j)
                    .sum::<f64>()
        })
        .sum();
    let radiation: f64 = candidate
        .step()
        .thermal
        .values()
        .map(|x| x.report().radiation_heat_j)
        .sum();
    near(
        candidate.step().external_heat_j,
        authored + radiation,
        2e-12,
    );
    near(gas_energy + wall_energy - authored - radiation, 0.0, 3e-7);
    for d in 0..3 {
        let dp: f64 = original
            .iter()
            .map(|(id, gas)| candidate.gases()[id].momentum_kg_m_s()[d] - gas.momentum_kg_m_s()[d])
            .sum();
        near(dp, 0.0, 3e-13);
    }
    near(
        candidate
            .gases()
            .values()
            .map(|g| g.parcel().mass_kg())
            .sum(),
        original.values().map(|g| g.parcel().mass_kg()).sum(),
        2e-15,
    );
    assert!(
        candidate.step().absolute_energy_defect_j <= thermal_controls([0.0; 3]).energy_tolerance_j
    );
}

#[test]
fn simultaneous_finite_walls_heat_the_original_network_and_close_actual_total_energy() {
    scope(|cx, _| {
        let model = catalog(PURE, 0.0);
        let law = interface(&model, 2000.0);
        let gases = three(&model, [0.0; 3]);
        let saved = gases.clone();
        let chart = curve([1000.0; 2]);
        let hot = body(100.0);
        let cool = body(50.0);
        let other = body(60.0);
        let mut first = [
            FiniteGasBody {
                id: 9,
                environment: &hot,
                input: UniformEnthalpyStepInput {
                    internal_heat_j: 1.5,
                    ..wall(&chart, 600.0)
                },
            },
            FiniteGasBody {
                id: 3,
                environment: &cool,
                input: wall(&chart, 300.0),
            },
        ];
        let second = [FiniteGasBody {
            id: 20,
            environment: &other,
            input: wall(&chart, 900.0),
        }];
        let groups = [
            MolecularGasThermalBoundary {
                gas_heat_j: 2.0,
                ..boundary(2, &first)
            },
            boundary(4, &second),
        ];
        let mut links = [link(4, [2, 4], &law), link(1, [4, 8], &law)];
        links[0].surface.velocity_world_m_s = [2.0, -4.0, 1.0];
        links[1].surface.normal_left_to_right = [0.0, 1.0, 0.0];
        let controls = thermal_controls([7.0, -2.0, 4.0]);
        let candidate =
            thermal::step_candidate(cx, &gases, &links, DT, &solves(&gases), &groups, controls)
                .unwrap();
        for group in groups {
            let (expected_t, expected_walls) = analytic(&gases[&group.gas], group);
            let heated = &candidate.molecular_initial_gases()[&group.gas];
            near(temperature(heated), expected_t, 2e-7);
            assert_eq!(
                heated.parcel().mass_kg(),
                gases[&group.gas].parcel().mass_kg()
            );
            assert_eq!(
                heated.parcel().volume_m3(),
                gases[&group.gas].parcel().volume_m3()
            );
            assert_eq!(
                heated.momentum_kg_m_s(),
                gases[&group.gas].momentum_kg_m_s()
            );
            for (id, expected) in expected_walls {
                near(candidate.bodies()[&id].temperature_k(), expected, 2e-8);
            }
        }
        let native = network::step_candidate(
            candidate.molecular_initial_gases(),
            &links,
            DT,
            &solves(&gases),
            controls.exchange,
            || false,
        )
        .unwrap();
        assert_eq!(candidate.gases(), native.gases());
        assert_eq!(&candidate.step().molecular, native.step());
        assert_eq!(&candidate.molecular_initial_gases()[&8], &gases[&8]);
        assert!(candidate.step().thermal[&2].report().heat_into_gas_j > 0.0);
        assert!(candidate.step().thermal[&4].report().heat_into_gas_j > 0.0);
        check_total(&gases, &groups, &candidate);
        assert_eq!(candidate.step().duration_s, DT);
        assert_eq!(candidate.step().molecular.duration_s, DT);
        assert_eq!(gases, saved);

        first.reverse();
        links.reverse();
        let reordered = [
            boundary(4, &second),
            MolecularGasThermalBoundary {
                gas_heat_j: 2.0,
                ..boundary(2, &first)
            },
        ];
        let again = thermal::step_candidate(
            cx,
            &gases,
            &links,
            DT,
            &solves(&gases),
            &reordered,
            controls,
        )
        .unwrap();
        assert_eq!(candidate, again);
    });
}

#[test]
fn opposite_body_heat_remains_physical_when_net_gas_heat_is_exactly_zero() {
    scope(|cx, _| {
        let gases = BTreeMap::from([(2, gas(catalog(PURE, 0.0), T, P, [7.0, 12.0, -3.0]))]);
        let chart = curve([1000.0; 2]);
        let environment = body(100.0);
        let bodies = [
            FiniteGasBody {
                id: 1,
                environment: &environment,
                input: wall(&chart, 600.0),
            },
            FiniteGasBody {
                id: 2,
                environment: &environment,
                input: wall(&chart, 200.0),
            },
        ];
        let groups = [boundary(2, &bodies)];
        let candidate = thermal::step_candidate(
            cx,
            &gases,
            &[],
            DT,
            &solves(&gases),
            &groups,
            thermal_controls([0.0; 3]),
        )
        .unwrap();
        assert_eq!(candidate.gases(), &gases);
        let report = &candidate.step().thermal[&2];
        assert_eq!(report.report().heat_into_gas_j, 0.0);
        let expected =
            DT * 100.0 * WALL_AREA * 200.0 / (1.0 + DT * 100.0 * WALL_AREA / (WALL_MASS * WALL_CP));
        near(
            report.bodies()[0].report.solid_to_gas_heat_j,
            expected,
            2e-9,
        );
        near(
            report.bodies()[1].report.solid_to_gas_heat_j,
            -expected,
            2e-9,
        );
        assert!(candidate.bodies()[&1].temperature_k() < 600.0);
        assert!(candidate.bodies()[&2].temperature_k() > 200.0);
        assert_eq!(candidate.step().thermal_kinetic_change_j, 0.0);
        check_total(&gases, &groups, &candidate);
    });
}

#[test]
fn new_species_large_formation_datums_and_common_boost_keep_the_same_thermal_physics() {
    scope(|cx, _| {
        let chart = curve([1000.0; 2]);
        let environment = body(100.0);
        let bodies = [FiniteGasBody {
            id: 7,
            environment: &environment,
            input: wall(&chart, 650.0),
        }];
        let groups = [boundary(2, &bodies)];
        let run = |datum: f64, boost: [f64; 3]| {
            let model = catalog(MIXTURE, datum);
            let law = interface(&model, 2000.0);
            let only_a = model
                .with_composition(
                    Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap(),
                )
                .unwrap();
            let only_b = model
                .with_composition(
                    Composition::new(CompositionBasis::MoleFraction, vec![0.0, 1.0]).unwrap(),
                )
                .unwrap();
            let gases = BTreeMap::from([
                (
                    2,
                    gas(
                        only_a,
                        T,
                        180_000.0,
                        core::array::from_fn(|d| [0.0, 30.0, -5.0][d] + boost[d]),
                    ),
                ),
                (
                    4,
                    gas(
                        only_b,
                        700.0,
                        120_000.0,
                        core::array::from_fn(|d| [0.0, -20.0, 4.0][d] + boost[d]),
                    ),
                ),
            ]);
            let mut aperture = link(1, [2, 4], &law);
            aperture.surface.velocity_world_m_s =
                core::array::from_fn(|d| [2.0, -4.0, 1.0][d] + boost[d]);
            thermal::step_candidate(
                cx,
                &gases,
                &[aperture],
                DT,
                &solves(&gases),
                &groups,
                thermal_controls(core::array::from_fn(|d| [7.0, -2.0, 4.0][d] + boost[d])),
            )
            .unwrap()
        };
        let rest = run(0.0, [0.0; 3]);
        let boost = [500.0, -700.0, 300.0];
        let moved = run(1e12, boost);
        let b = SpeciesId::new("B").unwrap();
        assert!(
            !rest.molecular_initial_gases()[&2]
                .parcel()
                .species_inventory_kg()
                .unwrap()
                .contains_key(&b)
        );
        assert!(rest.gases()[&2].parcel().species_inventory_kg().unwrap()[&b] > 0.0);
        near(
            rest.bodies()[&7].temperature_k(),
            moved.bodies()[&7].temperature_k(),
            2e-8,
        );
        for (&id, original) in rest.gases() {
            near(
                temperature(&moved.gases()[&id]),
                temperature(original),
                3e-7,
            );
            near(
                moved.gases()[&id].parcel().mass_kg(),
                original.parcel().mass_kg(),
                2e-15,
            );
            vector_near(
                moved.gases()[&id].velocity_m_s(),
                core::array::from_fn(|d| original.velocity_m_s()[d] + boost[d]),
                3e-10,
            );
        }
    });
}

#[test]
fn final_convection_uses_the_combined_gas_and_drift_and_refusal_is_atomic() {
    scope(|cx, gate| {
        let model = catalog(MIXTURE, 0.0);
        let law = interface(&model, 2000.0);
        let gases = three(&model, [0.0; 3]);
        let saved = gases.clone();
        let chart = curve([1000.0; 2]);
        let environment = body(0.0);
        let bodies = [FiniteGasBody {
            id: 7,
            environment: &environment,
            input: wall(&chart, 600.0),
        }];
        let groups = [boundary(2, &bodies)];
        let links = [link(1, [2, 4], &law)];
        let candidate = thermal::step_candidate_with_convection(
            cx,
            &gases,
            &links,
            DT,
            &solves(&gases),
            &groups,
            thermal_controls([0.0; 3]),
            |id, body, state, velocity| {
                Ok::<_, &'static str>(ConvectionCoefficient {
                    coefficient_w_per_m2_k: state.state().temperature / 4.0,
                    source: (id, state.clone(), velocity, body.input.initial),
                })
            },
        )
        .unwrap();
        let receipt = &candidate.step().final_convection[&7];
        assert_eq!(receipt.source.0, 2);
        assert_eq!(&receipt.source.1, candidate.gases()[&2].parcel().state());
        assert_eq!(receipt.source.2, candidate.gases()[&2].velocity_m_s());
        assert_eq!(receipt.source.3, candidate.bodies()[&7]);
        assert_ne!(receipt.source.3, bodies[0].input.initial);
        near(
            receipt.coefficient_w_per_m2_k,
            temperature(&candidate.gases()[&2]) / 4.0,
            0.0,
        );
        assert_ne!(receipt.source.2, gases[&2].velocity_m_s());
        let heat_receipt = &candidate.step().thermal[&2].bodies()[0].report.convection;
        assert_eq!(
            &heat_receipt.source.1,
            candidate.molecular_initial_gases()[&2].parcel().state()
        );
        assert_eq!(heat_receipt.source.2, gases[&2].velocity_m_s());
        assert_eq!(heat_receipt.source.3, bodies[0].input.initial);
        check_total(&gases, &groups, &candidate);

        let mut calls = 0;
        let refused = thermal::step_candidate_with_convection(
            cx,
            &gases,
            &links,
            DT,
            &solves(&gases),
            &groups,
            thermal_controls([0.0; 3]),
            |id, _, _, velocity| {
                calls += 1;
                if velocity != gases[&id].velocity_m_s() {
                    Err("post-exchange drift outside law")
                } else {
                    Ok(ConvectionCoefficient {
                        coefficient_w_per_m2_k: 100.0,
                        source: (),
                    })
                }
            },
        );
        assert!(calls > 1);
        assert!(matches!(
            refused,
            Err(MolecularNetworkThermalError::Thermal {
                gas: 2,
                error: FiniteGasBodiesError::Body {
                    id: 7,
                    error: ConvectionError::Convection("post-exchange drift outside law")
                }
            })
        ));
        let cancelled = thermal::step_candidate_with_convection(
            cx,
            &gases,
            &links,
            DT,
            &solves(&gases),
            &groups,
            thermal_controls([0.0; 3]),
            |id, _, _, velocity| {
                if velocity != gases[&id].velocity_m_s() {
                    gate.request();
                }
                Ok::<_, &'static str>(ConvectionCoefficient {
                    coefficient_w_per_m2_k: 100.0,
                    source: (),
                })
            },
        );
        assert!(matches!(
            cancelled,
            Err(MolecularNetworkThermalError::Cancelled)
        ));
        assert_eq!(gases, saved);
        assert_eq!(bodies[0].input.initial.temperature_k(), 600.0);
    });
}

#[test]
fn original_source_admission_precedes_heat_and_late_gross_withdrawal_refuses_all_state() {
    scope(|cx, _| {
        let model = catalog(PURE, 0.0);
        let law = interface(&model, 2000.0);
        let mut gases = three(&model, [0.0; 3]);
        gases.insert(2, gas(model.clone(), T, P, [3000.0, 0.0, 0.0]));
        let chart = curve([1000.0; 2]);
        let environment = body(100.0);
        let bodies = [FiniteGasBody {
            id: 7,
            environment: &environment,
            input: wall(&chart, 600.0),
        }];
        let groups = [boundary(2, &bodies)];
        let links = [link(1, [2, 4], &law)];
        let mut calls = 0;
        let refused = thermal::step_candidate_with_convection(
            cx,
            &gases,
            &links,
            DT,
            &solves(&gases),
            &groups,
            thermal_controls([0.0; 3]),
            |_, _, _, _| {
                calls += 1;
                Ok::<_, ()>(ConvectionCoefficient {
                    coefficient_w_per_m2_k: 100.0,
                    source: (),
                })
            },
        );
        assert!(matches!(
            refused,
            Err(MolecularNetworkThermalError::Molecular(
                MolecularNetworkError::Link { id: 1, .. }
            ))
        ));
        assert_eq!(calls, 0);

        gases.insert(2, gas(model, T, P, [0.0; 3]));
        let saved = gases.clone();
        let mut limits = thermal_controls([0.0; 3]);
        limits.exchange.maximum_species_fraction = 1e-12;
        let late = thermal::step_candidate_with_convection(
            cx,
            &gases,
            &links,
            DT,
            &solves(&gases),
            &groups,
            limits,
            |_, _, _, _| {
                calls += 1;
                Ok::<_, ()>(ConvectionCoefficient {
                    coefficient_w_per_m2_k: 100.0,
                    source: (),
                })
            },
        );
        assert!(calls > 1);
        assert!(matches!(
            late,
            Err(MolecularNetworkThermalError::Molecular(
                MolecularNetworkError::Withdrawal { .. }
            ))
        ));
        assert_eq!(gases, saved);
        assert_eq!(bodies[0].input.initial.temperature_k(), 600.0);
    });
}

#[test]
fn duplicate_body_identity_expansion_and_inconsistent_clocks_or_geometry_are_refused() {
    scope(|cx, _| {
        let gases = three(&catalog(PURE, 0.0), [0.0; 3]);
        let chart = curve([1000.0; 2]);
        let environment = body(100.0);
        let input = wall(&chart, 600.0);
        let bodies = [FiniteGasBody {
            id: 7,
            environment: &environment,
            input,
        }];
        let duplicates = [boundary(2, &bodies), boundary(4, &bodies)];
        assert!(matches!(
            thermal::step_candidate(
                cx,
                &gases,
                &[],
                DT,
                &solves(&gases),
                &duplicates,
                thermal_controls([0.0; 3])
            ),
            Err(MolecularNetworkThermalError::Invalid(_))
        ));
        let expanding = curve([1000.0, 950.0]);
        for invalid in [
            UniformEnthalpyStepInput {
                volume_m3: 2.0 * input.volume_m3,
                ..input
            },
            UniformEnthalpyStepInput {
                duration_s: DT * 0.5,
                ..input
            },
            wall(&expanding, 600.0),
        ] {
            let bodies = [FiniteGasBody {
                id: 7,
                environment: &environment,
                input: invalid,
            }];
            assert!(matches!(
                thermal::step_candidate(
                    cx,
                    &gases,
                    &[],
                    DT,
                    &solves(&gases),
                    &[boundary(2, &bodies)],
                    thermal_controls([0.0; 3])
                ),
                Err(MolecularNetworkThermalError::Invalid(_))
            ));
        }
        let mut limits = thermal_controls([0.0; 3]);
        limits.maximum_bodies = 0;
        assert!(matches!(
            thermal::step_candidate(
                cx,
                &gases,
                &[],
                DT,
                &solves(&gases),
                &[boundary(2, &bodies)],
                limits
            ),
            Err(MolecularNetworkThermalError::Invalid(_))
        ));
    });
}

#[test]
fn empty_thermal_set_and_insulated_body_reproduce_native_network_without_extra_heat() {
    scope(|cx, _| {
        let model = catalog(MIXTURE, 0.0);
        let law = interface(&model, 2000.0);
        let gases = three(&model, [0.0; 3]);
        let links = [link(1, [2, 4], &law)];
        let limits = thermal_controls([7.0, -2.0, 4.0]);
        let native =
            network::step_candidate(&gases, &links, DT, &solves(&gases), limits.exchange, || {
                false
            })
            .unwrap();
        let empty = thermal::step_candidate(
            cx,
            &gases,
            &links,
            DT,
            &solves(&gases),
            &[],
            MolecularNetworkThermalControls {
                maximum_bodies: 0,
                ..limits
            },
        )
        .unwrap();
        assert_eq!(empty.gases(), native.gases());
        assert_eq!(&empty.step().molecular, native.step());
        assert!(empty.bodies().is_empty());
        assert!(empty.step().final_convection.is_empty());
        let chart = curve([1000.0; 2]);
        let environment = body(0.0);
        let bodies = [FiniteGasBody {
            id: 7,
            environment: &environment,
            input: wall(&chart, 600.0),
        }];
        let insulated = thermal::step_candidate(
            cx,
            &gases,
            &links,
            DT,
            &solves(&gases),
            &[boundary(2, &bodies)],
            limits,
        )
        .unwrap();
        assert_eq!(insulated.gases(), native.gases());
        assert_eq!(insulated.bodies()[&7], bodies[0].input.initial);
        assert_eq!(insulated.step().external_heat_j, 0.0);
    });
}

#[test]
fn one_fixed_non_cancelling_budget_bounds_all_native_stage_defects() {
    scope(|cx, _| {
        let model = catalog(PURE, 0.0);
        let law = interface(&model, 2000.0);
        let gases = three(&model, [0.0; 3]);
        let chart = curve([1000.0; 2]);
        let mut environment = body(83.0);
        environment.transport = LumpedThermalTransport::try_declared(100.0, 0.2).unwrap();
        environment.radiation_temperature_k = 900.0;
        let bodies = [FiniteGasBody {
            id: 7,
            environment: &environment,
            input: wall(&chart, 613.0),
        }];
        let groups = [boundary(2, &bodies)];
        let links = [link(1, [2, 4], &law)];
        let limits = thermal_controls([7.0, -2.0, 4.0]);
        let accepted =
            thermal::step_candidate(cx, &gases, &links, DT, &solves(&gases), &groups, limits)
                .unwrap();
        assert!(accepted.step().thermal[&2].report().radiation_heat_j > 0.0);
        check_total(&gases, &groups, &accepted);
        let defect = accepted.step().absolute_energy_defect_j;
        assert!(defect > 0.0);
        assert!(defect >= accepted.step().thermal_absolute_energy_defect_j);
        let tight = MolecularNetworkThermalControls {
            energy_tolerance_j: defect * 0.5,
            ..limits
        };
        assert!(matches!(
            thermal::step_candidate(cx, &gases, &links, DT, &solves(&gases), &groups, tight),
            Err(MolecularNetworkThermalError::EnergyBalance { .. })
        ));
    });
}
