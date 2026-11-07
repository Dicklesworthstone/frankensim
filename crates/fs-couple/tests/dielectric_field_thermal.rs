#![cfg(feature = "resolved-dielectric-thermal")]
//! Native distributed polarization loss reaches real finite thermal bodies.
//! Coefficients are synthetic SI fixtures, not measured-material validation.

#[path = "../../fs-conduction/tests/radiation_material/common.rs"]
mod radiative;

#[path = "dielectric_field_thermal/lorentz.rs"]
mod lorentz;

use fs_couple::{
    circuit_enthalpy::{
        CircuitEnthalpyError, CircuitEnthalpySolve,
        conduction::{ConductionRegion, ConductionThermalNetwork},
        contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
        material::{ConductionLink, MaterialCell, MaterialNetwork},
    },
    material_circuit::MaterialCircuitSolve,
};
use fs_em::electrostatic::{
    ElectrostaticControls,
    dispersion::{DielectricFieldState, PlanarDebyeDielectric},
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    conductor::ELECTRICAL_CONDUCTIVITY_DIMS,
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS,
        dielectric::{
            DebyePoleProperties, bulk::ResolvedDebyeMaterial, leaky::ResolvedLeakyDielectric,
        },
    },
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{
        MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
    },
    thermal::source::{EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, ThermalMaterialPoint},
};
use fs_qty::{
    Density, QuantitySpec, Temperature, Time,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::{CheckedCgConfig, NewtonKrylovConfig};
use fs_thermal::enthalpy::*;

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            t,
        )
        .unwrap()
}
fn card(
    epsilon: f64,
    delta: f64,
    tau: f64,
    cp: f64,
    k: f64,
    upper: f64,
    variable: bool,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value, caloric) in [
        ("epsilon_inf", ELECTRIC_PERMITTIVITY_DIMS, epsilon, false),
        ("delta", ELECTRIC_PERMITTIVITY_DIMS, delta, false),
        ("tau", Time::DIMS, tau, false),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            cp,
            true,
        ),
        ("density", Density::DIMS, 1.0, true),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, k, false),
        (SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, 0.7, false),
    ] {
        let varying = variable && name == "epsilon_inf";
        let (value, interpolation) = if caloric || varying {
            (
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (250.0, value),
                        (500.0, if varying { value * 1.2 } else { value }),
                    ],
                    dims,
                },
                InterpolationPolicy::LinearInside,
            )
        } else {
            (
                PropertyValue::Scalar { value, dims },
                InterpolationPolicy::ConstantWithinValidity,
            )
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    250.0,
                    if caloric { 500.0 } else { upper },
                ),
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic spatial dielectric heat fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "finite thermal receiver".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(250.0),
            upper: &point(500.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-9,
    )
    .unwrap()
}

fn with_leakage(base: MaterialCard, sigma: f64, variable: bool) -> MaterialCard {
    let mut claims = base.claims().clone();
    let dims = ELECTRICAL_CONDUCTIVITY_DIMS;
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("sigma", dims),
            value: if variable {
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(250.0, sigma), (500.0, 1.2 * sigma)],
                    dims,
                }
            } else {
                PropertyValue::Scalar { value: sigma, dims }
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                250.0,
                500.0,
            ),
            interpolation: if variable {
                InterpolationPolicy::LinearInside
            } else {
                InterpolationPolicy::ConstantWithinValidity
            },
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic source-backed DC leakage".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}
fn field_controls() -> ElectrostaticControls {
    ElectrostaticControls {
        linear: CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 128,
            max_corrections: 2,
        },
        charge_tolerance_c: 1e-10,
        energy_tolerance_j: 1e-10,
    }
}
fn controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 2,
        maximum_circuit_trials: 500,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-4,
        cell_energy_tolerance_j: 2e-8,
        thermal_energy_tolerance_j: 2e-8,
        total_energy_tolerance_j: 4e-8,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-10,
            power_tolerance_w: 1e-10,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 4,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn ends(count: usize) -> Vec<(usize, f64)> {
    if count == 1 {
        vec![(0, 1.0), (1, 0.0), (2, 0.0), (3, 1.0)]
    } else {
        vec![(0, 2.0), (3, 2.0), (2, 0.0), (5, 0.0)]
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

fn with_case<R>(
    cards: &[MaterialCard],
    temperatures: &[f64],
    poles: bool,
    linked: bool,
    run: impl FnOnce(
        &ConductionThermalNetwork<'_>,
        DepositionSimulation<'_>,
        PlanarDebyeDielectric,
    ) -> R,
) -> R {
    let curves: Vec<_> = cards.iter().map(caloric).collect();
    let points: Vec<_> = temperatures.iter().copied().map(point).collect();
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells: Vec<_> = cards
        .iter()
        .enumerate()
        .map(|(i, card)| MaterialCell {
            cell: CellId(7 + i as u64),
            material_id: MaterialId(11 + i as u64),
            curve: curves[i].curve(),
            point: ThermalMaterialPoint {
                material: card,
                context: &points[i],
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: linked.then_some((&key, MaterialPropertySelection::SingleClaimOnly)),
        })
        .collect();
    let mut events: Vec<_> = temperatures
        .iter()
        .enumerate()
        .map(|(i, &temperature)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: cells[i].cell,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: cells[i].material_id,
                mass_kg: 1.0,
                specific_h_j_kg: {
                    let knots = curves[i].curve().knots();
                    knots[0].specific_enthalpy_j_kg
                        + (temperature - knots[0].temperature_k)
                            * (knots[1].specific_enthalpy_j_kg - knots[0].specific_enthalpy_j_kg)
                            / (knots[1].temperature_k - knots[0].temperature_k)
                },
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    let links: Vec<_> = if linked {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 19,
                a: cells[0].cell,
                b: cells[1].cell,
                area_m2: 1.0,
                conductance_w_k: 0.0,
            }),
        });
        vec![ConductionLink {
            link: 19,
            distance_a_m: 0.25,
            distance_b_m: 0.25,
        }]
    } else {
        vec![]
    };
    let laws: Vec<_> = cells
        .iter()
        .map(|cell| (cell.material_id, cell.curve as &dyn ThermalLaw))
        .collect();
    let thermal = DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let roles = if poles {
        vec![DebyePoleProperties {
            permittivity_increment: "delta".into(),
            relaxation_time: "tau".into(),
        }]
    } else {
        vec![]
    };
    let materials: Vec<_> = cards
        .iter()
        .zip(&points)
        .map(|(card, point)| {
            ResolvedDebyeMaterial::resolve(
                card,
                point,
                MaterialPropertySelection::SingleClaimOnly,
                "epsilon_inf",
                &roles,
            )
            .unwrap()
        })
        .collect();
    let (vertices, faces, voltages) = if cards.len() == 1 {
        (
            vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
            vec![[0, 1, 2], [0, 2, 3]],
            vec![1.0, 0.0, 0.0, 1.0],
        )
    } else {
        let e0 = materials[0].high_frequency_permittivity_f_m();
        let e1 = materials[1].high_frequency_permittivity_f_m();
        let middle = 2.0 * e0 / (e0 + e1);
        (
            vec![
                [0.0, 0.0],
                [1.0, 0.0],
                [2.0, 0.0],
                [0.0, 1.0],
                [1.0, 1.0],
                [2.0, 1.0],
            ],
            vec![[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]],
            vec![2.0, middle, 0.0, 2.0, middle, 0.0],
        )
    };
    let mesh = TriComplex2::from_indexed_triangles(
        "dielectric-thermal-fixture",
        vertices,
        faces,
        Metric2::planar(1.0).unwrap(),
    )
    .unwrap();
    let regions: Vec<_> = (0..mesh.faces().len())
        .map(|triangle| ConductionRegion {
            triangle,
            cell: cells[triangle / 2].cell,
        })
        .collect();
    let sources: Vec<_> = regions
        .iter()
        .map(|r| materials[r.triangle / 2].clone())
        .collect();
    let polarization = sources
        .iter()
        .map(|law| vec![[0.0; 2]; law.poles().len()])
        .collect();
    let field = PlanarDebyeDielectric::new(
        mesh,
        sources,
        voltages,
        polarization,
        &ends(cards.len())
            .iter()
            .map(|&(vertex, _)| vertex)
            .collect::<Vec<_>>(),
        0.0,
        WeightedAssemblyLimits {
            max_cells: 10,
            max_dofs: 20,
            max_triplets: 200,
        },
        field_controls(),
        &mut || false,
    )
    .unwrap();
    run(
        &ConductionThermalNetwork {
            cells: &cells,
            links: &links,
            regions: &regions,
            volume_tolerance_m3: 1e-12,
        },
        thermal,
        field,
    )
}

#[test]
fn polarization_heat_changes_real_temperature_and_restarts_with_exact_work_account() {
    let cards = [card(2.0, 3.0, 0.2, 10.0, 1.0, 500.0, false)];
    with_case(
        &cards,
        &[300.0],
        true,
        false,
        |network, mut thermal, mut field| {
            let dt = 0.05;
            let prepared = field.prepare_step(dt, &mut || false).unwrap();
            let initial = thermal.energy_j() + field.stored_energy_j();
            let r = (2.0 * 0.2 - dt) / (2.0 * 0.2 + dt);
            let mut work = 0.0;
            for n in 1..=12 {
                let next = network
                    .step_debye_candidate(
                        &thermal,
                        field.time_s() + dt,
                        &field,
                        &prepared,
                        &ends(1),
                        field_controls(),
                        controls(),
                        &[0.0],
                        || false,
                    )
                    .unwrap();
                let p = 3.0 * (1.0 - f64::powi(r, n));
                let heat = p - p * p / 6.0;
                near(next.field().polarization_c_m2()[0][0][0], p, 2e-12);
                near(
                    next.thermal().temperature_k(CellId(7)).unwrap(),
                    300.0 + heat / 10.0,
                    2e-9,
                );
                near(
                    next.report().coupled.deposited_joule_j,
                    next.report().coupled.electrical.heat_j,
                    1e-12,
                );
                let snapshot = &next.report().coupled.source.dielectric[0];
                near(
                    snapshot.query().axes()["T"],
                    next.thermal().temperature_k(CellId(7)).unwrap(),
                    1e-12,
                );
                work += next.report().coupled.electrical.supplied_work_j;
                near(
                    next.thermal().energy_j() + next.field().stored_energy_j(),
                    initial + work,
                    4e-8,
                );
                near(next.report().external_heat_j, 0.0, 2e-10);
                near(next.report().total_energy_residual_j, 0.0, 4e-8);
                (thermal, field, _) = next.into_parts();
            }
        },
    );
}

#[test]
fn dc_leakage_and_polarization_heat_remain_separate_and_revalidate_zero_conductivity() {
    for sigma in [0.0, 4.0] {
        for variable in [false, true] {
            let cards = [with_leakage(
                card(2.0, 3.0, 0.2, 10.0, 1.0, 500.0, false),
                sigma,
                variable,
            )];
            with_case(
                &cards,
                &[300.0],
                true,
                false,
                |network, thermal, reference| {
                    let roles = [DebyePoleProperties {
                        permittivity_increment: "delta".into(),
                        relaxation_time: "tau".into(),
                    }];
                    let source = ResolvedLeakyDielectric::resolve(
                        &cards[0],
                        &point(300.0),
                        MaterialPropertySelection::SingleClaimOnly,
                        "epsilon_inf",
                        &roles,
                        "sigma",
                    )
                    .unwrap();
                    for initial_p in [0.0, 3.0] {
                        let field = PlanarDebyeDielectric::new_leaky(
                            reference.mesh().clone(),
                            vec![source.clone(); 2],
                            reference.voltage_v().to_vec(),
                            vec![vec![[initial_p, 0.0]]; 2],
                            &[0, 1, 2, 3],
                            0.0,
                            WeightedAssemblyLimits {
                                max_cells: 10,
                                max_dofs: 20,
                                max_triplets: 200,
                            },
                            field_controls(),
                            &mut || false,
                        )
                        .unwrap();
                        let dt = 0.1;
                        let prepared = field.prepare_step(dt, &mut || false).unwrap();
                        let native = field
                            .step_candidate(&prepared, &ends(1), field_controls(), &mut || false)
                            .unwrap();
                        let result = network.step_dielectric_field_candidate(
                            &thermal,
                            &native,
                            controls(),
                            &[0.0],
                            || false,
                        );
                        if variable {
                            // Even a sourced zero must not bypass the same constant
                            // conductivity law admission on the thermal trajectory.
                            assert!(result.is_err());
                            assert_eq!(field.time_s(), 0.0);
                            assert_eq!(thermal.time_s(), 0.0);
                            continue;
                        }
                        let next = result.unwrap();
                        let p = initial_p + 2.0 * dt / (2.0 * 0.2 + dt) * (3.0 - initial_p);
                        let polarization_heat = 0.2 * (p - initial_p).powi(2) / (3.0 * dt);
                        let joule_heat = dt * sigma;
                        let report = &next.report().coupled.electrical;
                        near(report.polarization_heat_j, polarization_heat, 3e-12);
                        near(report.joule_heat_j, joule_heat, 3e-12);
                        near(report.heat_j, polarization_heat + joule_heat, 3e-12);
                        near(report.supplied_work_j, p - initial_p + joule_heat, 3e-12);
                        near(report.supplied_displacement_work_j, p - initial_p, 3e-12);
                        near(report.supplied_conduction_work_j, joule_heat, 3e-12);
                        near(
                            next.thermal().temperature_k(CellId(7)).unwrap(),
                            300.0 + (polarization_heat + joule_heat) / 10.0,
                            2e-9,
                        );
                        near(
                            next.report().coupled.deposited_joule_j,
                            report.heat_j,
                            3e-12,
                        );
                        for cell in 0..2 {
                            near(
                                report.cell_heat_j[cell],
                                report.cell_polarization_heat_j[cell]
                                    + report.cell_joule_heat_j[cell],
                                3e-12,
                            );
                            let key = PropertyKey::new("sigma", ELECTRICAL_CONDUCTIVITY_DIMS);
                            let snapshot = &next.report().coupled.source.dielectric[cell];
                            let original = source.material().property_by_key(&key).unwrap();
                            let refreshed = snapshot.property_by_key(&key).unwrap();
                            assert_eq!(
                                refreshed.answer().receipt.selected,
                                original.answer().receipt.selected
                            );
                            near(refreshed.value_si(), sigma, 0.0);
                        }
                        near(next.report().total_energy_residual_j, 0.0, 4e-8);
                        // A borrowed candidate can be retried after refusal, but
                        // its loss cannot be deposited twice into advanced heat.
                        assert!(
                            network
                                .step_dielectric_field_candidate(
                                    next.thermal(),
                                    &native,
                                    controls(),
                                    &[0.0],
                                    || false
                                )
                                .is_err()
                        );
                    }
                },
            );
        }
    }
}

#[test]
fn revolved_dielectric_uses_physical_volumes_in_the_same_thermal_bridge() {
    use fs_em::electrostatic::dispersion::AxisymmetricDebyeDielectric;
    let cards = [with_leakage(
        card(2.0, 3.0, 0.2, 10.0, 1.0, 500.0, false),
        4.0,
        false,
    )];
    with_case(&cards, &[300.0], true, false, |network, thermal, _| {
        // A full cylinder of unit volume. Its meridional area is 1/sqrt(pi),
        // so a planar-area receiver check would reject this valid coupling.
        let radius = core::f64::consts::PI.recip().sqrt();
        let mesh = TriComplex2::from_indexed_triangles(
            "revolved-dielectric-thermal",
            vec![[0.0, 0.0], [radius, 0.0], [radius, 1.0], [0.0, 1.0]],
            vec![[0, 1, 2], [0, 2, 3]],
            Metric2::axisymmetric(core::f64::consts::TAU).unwrap(),
        )
        .unwrap();
        let source = ResolvedLeakyDielectric::resolve(
            &cards[0],
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly,
            "epsilon_inf",
            &[DebyePoleProperties {
                permittivity_increment: "delta".into(),
                relaxation_time: "tau".into(),
            }],
            "sigma",
        )
        .unwrap();
        let field = AxisymmetricDebyeDielectric::new_leaky(
            mesh,
            vec![source; 2],
            vec![1.0, 1.0, 0.0, 0.0],
            vec![vec![[0.0, 0.0]]; 2],
            &[0, 1, 2, 3],
            0.0,
            WeightedAssemblyLimits {
                max_cells: 10,
                max_dofs: 20,
                max_triplets: 200,
            },
            field_controls(),
            &mut || false,
        )
        .unwrap();
        let prepared = field.prepare_step(0.1, &mut || false).unwrap();
        let native = field
            .step_candidate(
                &prepared,
                &[(0, 1.0), (1, 1.0), (2, 0.0), (3, 0.0)],
                field_controls(),
                &mut || false,
            )
            .unwrap();
        let next = network
            .step_dielectric_field_candidate(&thermal, &native, controls(), &[0.0], || false)
            .unwrap();
        let heat = 0.2 * 1.2_f64.powi(2) / (3.0 * 0.1) + 0.1 * 4.0;
        near(next.field().cell_volumes_m3()[0], 2.0 / 3.0, 3e-15);
        near(next.field().cell_volumes_m3()[1], 1.0 / 3.0, 3e-15);
        near(
            next.report().coupled.electrical.cell_heat_j[0],
            heat * 2.0 / 3.0,
            3e-11,
        );
        near(
            next.report().coupled.electrical.cell_heat_j[1],
            heat / 3.0,
            3e-11,
        );
        near(
            next.thermal().temperature_k(CellId(7)).unwrap(),
            300.0 + heat / 10.0,
            2e-9,
        );
        near(next.report().total_energy_residual_j, 0.0, 4e-8);
    });
}

#[test]
fn tetrahedral_dielectric_heats_the_same_native_caloric_owner() {
    use fs_em::spatial::electrostatic::dispersion::SpatialDebyeDielectric;
    let cards = [with_leakage(
        card(2.0, 3.0, 0.2, 10.0, 1.0, 500.0, false),
        4.0,
        false,
    )];
    with_case(&cards, &[300.0], true, false, |network, thermal, _| {
        let source = ResolvedLeakyDielectric::resolve(
            &cards[0],
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly,
            "epsilon_inf",
            &[DebyePoleProperties {
                permittivity_increment: "delta".into(),
                relaxation_time: "tau".into(),
            }],
            "sigma",
        )
        .unwrap();
        // One unit-volume tetrahedron with a fully three-component affine
        // electric field [1,-2,0.5] under its original complete voltage traces.
        let voltage = vec![3.0, 2.0, 7.0, 1.5];
        let field = SpatialDebyeDielectric::new_leaky(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, 0.0, 3.0],
            ],
            vec![[0, 1, 2, 3]],
            vec![source],
            voltage.clone(),
            vec![vec![[0.0; 3]]],
            &[0, 1, 2, 3],
            0.0,
            WeightedAssemblyLimits {
                max_cells: 10,
                max_dofs: 20,
                max_triplets: 200,
            },
            field_controls(),
            &mut || false,
        )
        .unwrap();
        let prepared = field.prepare_step(0.1, &mut || false).unwrap();
        let native = field
            .step_candidate(
                &prepared,
                &voltage.into_iter().enumerate().collect::<Vec<_>>(),
                field_controls(),
                &mut || false,
            )
            .unwrap();
        let regions = [ConductionRegion {
            triangle: 0,
            cell: CellId(7),
        }];
        let layout = ConductionThermalNetwork {
            cells: network.cells,
            links: &[],
            regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        let next = layout
            .step_dielectric_field_candidate(&thermal, &native, controls(), &[0.0], || false)
            .unwrap();
        let norm_squared = 1.0 + 4.0 + 0.25;
        let polarization_heat = 0.2 * 1.2_f64.powi(2) / (3.0 * 0.1) * norm_squared;
        let joule_heat = 0.1 * 4.0 * norm_squared;
        near(
            DielectricFieldState::cell_volumes_m3(next.field())[0],
            1.0,
            1e-14,
        );
        near(
            next.report().coupled.electrical.polarization_heat_j,
            polarization_heat,
            3e-11,
        );
        near(
            next.report().coupled.electrical.joule_heat_j,
            joule_heat,
            3e-11,
        );
        near(
            next.report().coupled.electrical.supplied_work_j,
            1.2 * norm_squared + joule_heat,
            3e-11,
        );
        near(
            next.thermal().temperature_k(CellId(7)).unwrap(),
            300.0 + (polarization_heat + joule_heat) / 10.0,
            2e-9,
        );
        near(next.report().total_energy_residual_j, 0.0, 4e-8);
    });
}

#[test]
fn heterogeneous_polarization_and_fourier_transport_match_independent_two_layer_solution() {
    let cards = [
        card(1.0, 0.8, 0.2, 2.0, 1.0, 500.0, false),
        card(2.0, 1.5, 0.4, 3.0, 3.0, 500.0, false),
    ];
    with_case(
        &cards,
        &[300.0, 320.0],
        true,
        true,
        |network, thermal, field| {
            let dt = 0.1;
            let prepared = field.prepare_step(dt, &mut || false).unwrap();
            let next = network
                .step_debye_candidate(
                    &thermal,
                    dt,
                    &field,
                    &prepared,
                    &ends(2),
                    field_controls(),
                    controls(),
                    &[0.0, 0.0],
                    || false,
                )
                .unwrap();
            let epsilon = [1.0, 2.0];
            let delta = [0.8, 1.5];
            let tau = [0.2, 0.4];
            let old_e = [4.0 / 3.0, 2.0 / 3.0];
            let b = tau.map(|t| dt / (2.0 * t + dt));
            let a = [epsilon[0] + b[0] * delta[0], epsilon[1] + b[1] * delta[1]];
            let history = [b[0] * delta[0] * old_e[0], b[1] * delta[1] * old_e[1]];
            let displacement =
                (2.0 + history[0] / a[0] + history[1] / a[1]) / (1.0 / a[0] + 1.0 / a[1]);
            let e = [
                (displacement - history[0]) / a[0],
                (displacement - history[1]) / a[1],
            ];
            let p = [
                b[0] * delta[0] * (old_e[0] + e[0]),
                b[1] * delta[1] * (old_e[1] + e[1]),
            ];
            let heat = [
                tau[0] * p[0] * p[0] / (delta[0] * dt),
                tau[1] * p[1] * p[1] / (delta[1] * dt),
            ];
            for i in 0..2 {
                near(next.field().electric_field_v_m()[2 * i][0], e[i], 3e-11);
                near(
                    next.report().coupled.electrical.cell_heat_j[2 * i] * 2.0,
                    heat[i],
                    3e-11,
                );
            }
            let g = 1.0 / (0.25 / 1.0 + 0.25 / 3.0);
            let rhs = [heat[0] + dt * g * 20.0, heat[1] - dt * g * 20.0];
            let diagonal = [2.0 + dt * g, 3.0 + dt * g];
            let determinant = diagonal[0] * diagonal[1] - (dt * g) * (dt * g);
            let increments = [
                (diagonal[1] * rhs[0] + dt * g * rhs[1]) / determinant,
                (dt * g * rhs[0] + diagonal[0] * rhs[1]) / determinant,
            ];
            for i in 0..2 {
                near(
                    next.thermal().temperature_k(CellId(7 + i as u64)).unwrap(),
                    [300.0, 320.0][i] + increments[i],
                    2e-9,
                );
            }
            near(
                next.report().coupled.source.thermal.links[0].conductance_w_k,
                g,
                1e-12,
            );
            near(next.report().total_energy_residual_j, 0.0, 4e-8);
            // Equal volume and equal readable material IDs cannot hide a
            // transfer to a different immutable caloric card.
            let swapped: Vec<_> = network
                .regions
                .iter()
                .map(|region| ConductionRegion {
                    triangle: region.triangle,
                    cell: if region.cell == CellId(7) {
                        CellId(8)
                    } else {
                        CellId(7)
                    },
                })
                .collect();
            let wrong_owner = ConductionThermalNetwork {
                cells: network.cells,
                links: network.links,
                regions: &swapped,
                volume_tolerance_m3: network.volume_tolerance_m3,
            };
            assert!(
                wrong_owner
                    .step_debye_candidate(
                        &thermal,
                        dt,
                        &field,
                        &prepared,
                        &ends(2),
                        field_controls(),
                        controls(),
                        &[0.0, 0.0],
                        || false
                    )
                    .is_err()
            );
        },
    );
}

#[test]
fn full_closed_radiation_assembly_keeps_heat_internal_and_field_loss_local() {
    use fs_conduction::{
        SurfaceEmissivity,
        transient::deposition::radiation::{MaterialRadiation, RadiatingCell},
    };
    use fs_exec::{CancelGate, Cx};
    let cards = [
        card(1.0, 0.8, 0.2, 2.0, 1.0, 500.0, false),
        card(2.0, 1.5, 0.4, 3.0, 3.0, 500.0, false),
    ];
    with_case(
        &cards,
        &[300.0, 320.0],
        true,
        false,
        |network, thermal, field| {
            let bulk = MaterialNetwork {
                cells: network.cells,
                links: network.links,
                resistors: &[],
            };
            let electrical = ElectricalContactNetwork {
                bulk: &bulk,
                contacts: &[],
                electrical_bases: &[],
            };
            let surfaces = [0, 1].map(|i| RadiatingCell {
                surface: i,
                cell: network.cells[i].cell,
                material_id: network.cells[i].material_id,
                curve: network.cells[i].curve,
                material: network.cells[i].point,
                selection: MaterialPropertySelection::SingleClaimOnly,
            });
            let emissivities = [0, 1].map(|i| {
                SurfaceEmissivity::from_thermal_state(
                    "Debye thermal enclosure",
                    network.cells[i].point,
                    network.cells[i]
                        .curve
                        .state_at_specific_enthalpy(
                            thermal.cells()[&network.cells[i].cell].specific_h_j_kg,
                        )
                        .unwrap(),
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .unwrap()
                .0
            });
            let enclosure = radiative::enclosure(emissivities);
            let radiation = MaterialRadiation {
                enclosure: &enclosure,
                cells: &surfaces,
                power_tolerance_w: 1e-8,
                radiosity_tolerance_w_m2: 1e-8,
            };
            let heat = MaterialHeatNetwork {
                electrical: &electrical,
                thermal_contacts: &[],
                radiation: Some(&radiation),
            };
            let gate = CancelGate::new();
            let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
            pool.scope(|arena| {
                let cx = Cx::new(
                    &gate,
                    arena,
                    fs_exec::StreamKey {
                        seed: 12,
                        kernel_id: 19,
                        tile: 0,
                        iteration: 0,
                    },
                    fs_exec::Budget::INFINITE,
                    fs_exec::ExecMode::Deterministic,
                );
                let heat = heat.prepare(&cx, &thermal).unwrap();
                let dt = 0.01;
                let prepared = field.prepare_step(dt, &mut || false).unwrap();
                let next = heat
                    .step_debye_candidate(
                        &cx,
                        &thermal,
                        dt,
                        &field,
                        &prepared,
                        &ends(2),
                        network.regions,
                        network.volume_tolerance_m3,
                        field_controls(),
                        controls(),
                        &[0.0, 0.0],
                    )
                    .unwrap();
                let radiation = next
                    .report()
                    .coupled
                    .source
                    .thermal
                    .radiation
                    .as_ref()
                    .unwrap();
                assert!(radiation.sources_w()[&CellId(7)] > 0.0);
                assert!(radiation.sources_w()[&CellId(8)] < 0.0);
                for i in 0..2 {
                    let id = CellId(7 + i as u64);
                    let heat = next.report().coupled.electrical.cell_heat_j[2 * i..2 * i + 2]
                        .iter()
                        .sum::<f64>();
                    let change = next.thermal().cells()[&id].specific_h_j_kg
                        - thermal.cells()[&id].specific_h_j_kg;
                    near(change, heat + dt * radiation.sources_w()[&id], 2e-8);
                }
                near(next.report().external_heat_j, 0.0, 2e-8);
                near(next.report().total_energy_residual_j, 0.0, 4e-8);
            });
        },
    );
}

#[test]
fn zero_loss_changes_storage_without_heating_and_refusals_are_transactional() {
    let cards = [card(2.0, 3.0, 0.2, 10.0, 1.0, 500.0, false)];
    with_case(&cards, &[300.0], false, false, |network, thermal, field| {
        let prepared = field.prepare_step(0.1, &mut || false).unwrap();
        let drive: Vec<_> = ends(1).into_iter().map(|(i, v)| (i, 2.0 * v)).collect();
        let next = network
            .step_debye_candidate(
                &thermal,
                0.1,
                &field,
                &prepared,
                &drive,
                field_controls(),
                controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        near(
            next.field().stored_energy_j(),
            4.0 * field.stored_energy_j(),
            1e-12,
        );
        assert_eq!(next.report().coupled.electrical.heat_j, 0.0);
        assert_eq!(next.thermal().energy_j(), thermal.energy_j());
        assert!(
            network
                .step_debye_candidate(
                    &thermal,
                    0.2,
                    &field,
                    &prepared,
                    &drive,
                    field_controls(),
                    controls(),
                    &[0.0],
                    || false
                )
                .is_err()
        );
    });
    with_case(&cards, &[300.0], true, false, |network, thermal, field| {
        let prepared = field.prepare_step(0.1, &mut || false).unwrap();
        let count = std::cell::Cell::new(0usize);
        let run = |stop: usize| {
            network.step_debye_candidate(
                &thermal,
                0.1,
                &field,
                &prepared,
                &ends(1),
                field_controls(),
                controls(),
                &[0.0],
                || {
                    count.set(count.get() + 1);
                    count.get() >= stop
                },
            )
        };
        let baseline = run(usize::MAX).unwrap();
        let late = count.get() - 1;
        count.set(0);
        assert!(matches!(run(late), Err(CircuitEnthalpyError::Cancelled)));
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(field.time_s(), 0.0);
        assert!(
            field
                .polarization_c_m2()
                .iter()
                .flatten()
                .all(|p| *p == [0.0; 2])
        );
        count.set(0);
        let retry = run(usize::MAX).unwrap();
        assert_eq!(retry.thermal().energy_j(), baseline.thermal().energy_j());
        assert_eq!(
            retry.field().polarization_c_m2(),
            baseline.field().polarization_c_m2()
        );
        let duplicate = [network.regions[0], network.regions[0]];
        let bad = ConductionThermalNetwork {
            cells: network.cells,
            links: network.links,
            regions: &duplicate,
            volume_tolerance_m3: 1e-12,
        };
        assert!(
            bad.step_debye_candidate(
                &thermal,
                0.1,
                &field,
                &prepared,
                &ends(1),
                field_controls(),
                controls(),
                &[0.0],
                || false
            )
            .is_err()
        );
    });
    for source in [
        card(2.0, 3.0, 0.2, 10.0, 1.0, 500.0, true),
        card(2.0, 3.0, 0.2, 10.0, 1.0, 300.001, false),
    ] {
        with_case(
            &[source],
            &[300.0],
            true,
            false,
            |network, thermal, field| {
                let prepared = field.prepare_step(0.1, &mut || false).unwrap();
                assert!(
                    network
                        .step_debye_candidate(
                            &thermal,
                            0.1,
                            &field,
                            &prepared,
                            &ends(1),
                            field_controls(),
                            controls(),
                            &[0.0],
                            || false
                        )
                        .is_err()
                );
                assert_eq!(field.time_s(), 0.0);
                assert_eq!(thermal.time_s(), 0.0);
            },
        );
    }
}
