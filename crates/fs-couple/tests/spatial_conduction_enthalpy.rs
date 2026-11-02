//! G1/G3/G4: actual tetrahedral current closes material and finite heat feedback.
#[path = "circuit_enthalpy/common.rs"]
mod common;
#[path = "spatial_conduction_enthalpy/contacts.rs"]
mod contacts;
#[path = "spatial_conduction_enthalpy/current.rs"]
mod current;
#[path = "spatial_conduction_enthalpy/heat.rs"]
mod heat;
use common::{config, initial, near};
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError,
    conduction::spatial::*,
    material::{ConductionLink, MaterialCell},
};
use fs_em::{conduction::ConductionControls, spatial::conduction::*};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    conductor::{
        ElectricalPropertyBasis, ResolvedConductorState, resolve_conductor_state_with_basis,
        spatial::{ResolvedOrthotropicConductor3, SpatialConductorMaterial},
    },
    phase::{
        EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacitySourceQuery,
        ResolvedHeatCapacityCurve, SolidLiquidPhase,
    },
    state_point::{
        MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility,
        ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point,
    },
    thermal::source::ThermalMaterialPoint,
};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};
use fs_solver::CheckedCgConfig;
use fs_thermal::enthalpy::*;

fn point(temperature_k: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            temperature_k,
        )
        .unwrap()
}

/// Synthetic temperature-only claims: every source coordinate belongs to the
/// thermal owner; no applicability-only material-lot coordinate is invented.
fn material(rho0: f64, cp: f64, alpha: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let basis = ElectricalPropertyBasis::Resistivity;
    let mut claims = ClaimSet::new();
    for (key, lower, upper) in [
        (heat.clone(), cp, cp),
        (density.clone(), 1.0, 1.0),
        (
            PropertyKey::new(basis.property_name(), basis.dimensions()),
            rho0 * (1.0 - 100.0 * alpha),
            rho0 * (1.0 + 1700.0 * alpha),
        ),
        (
            PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            1.0,
            1.0,
        ),
    ] {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, lower), (2000.0, upper)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    200.0,
                    2000.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic tetrahedral electrothermal fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "authored test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(2000.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}

fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 100,
        max_dofs: 300,
        max_triplets: 3600,
    }
}
fn field_controls() -> ConductionControls {
    ConductionControls {
        linear: CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 200,
            max_corrections: 2,
        },
        current_tolerance_a: 1e-9,
        power_tolerance_w: 1e-9,
    }
}
fn slab(xs: &[f64]) -> (Vec<[f64; 3]>, Vec<[u32; 4]>) {
    let mut positions = Vec::new();
    for &x in xs {
        for z in [0.0, 1.0] {
            for y in [0.0, 1.0] {
                positions.push([x, y, z]);
            }
        }
    }
    let mut tets = Vec::new();
    for i in 0..xs.len() - 1 {
        let a = 4 * i as u32;
        let b = a + 4;
        let c = a + 1;
        let d = b + 1;
        let e = a + 2;
        let f = b + 2;
        let g = a + 3;
        let h = b + 3;
        tets.extend([
            [a, b, d, h],
            [a, d, c, h],
            [a, c, g, h],
            [a, g, e, h],
            [a, e, f, h],
            [a, f, b, h],
        ]);
    }
    (positions, tets)
}
fn boundary<C: SpatialConductorMaterial>(
    problem: &SpatialConduction<C>,
    length: f64,
    current: f64,
    voltage: Option<f64>,
) -> SpatialConductionBoundary {
    let mut result = SpatialConductionBoundary::default();
    for face in problem.exterior_faces() {
        let xs = problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0]);
        if xs == [0.0; 3] {
            result.voltage.push(ConductionVoltageFace {
                face,
                voltage_v: [0.0; 3],
            });
        } else if xs == [length; 3] {
            if let Some(voltage) = voltage {
                result.voltage.push(ConductionVoltageFace {
                    face,
                    voltage_v: [voltage; 3],
                });
            } else {
                result.current.push(ConductionCurrentFace {
                    face,
                    inward_current_density_a_m2: current,
                });
            }
        }
    }
    result
}

fn electrodes<C: SpatialConductorMaterial>(
    problem: &SpatialConduction<C>,
    length: f64,
) -> [Vec<usize>; 2] {
    [length, 0.0].map(|x| {
        problem
            .exterior_faces()
            .filter(|&face| {
                problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0]) == [x; 3]
            })
            .collect()
    })
}
fn bindings<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    context: &'a QueryPoint,
    key: &'a PropertyKey,
) -> Vec<MaterialCell<'a>> {
    [CellId(30), CellId(10)]
        .into_iter()
        .enumerate()
        .map(|(i, cell)| MaterialCell {
            cell,
            material_id: MaterialId(i as u64),
            curve: curves[i].curve(),
            point: ThermalMaterialPoint {
                material: &cards[i],
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: Some((key, MaterialPropertySelection::SingleClaimOnly)),
        })
        .collect()
}
const LINKS: [ConductionLink; 1] = [ConductionLink {
    link: 7,
    distance_a_m: 0.5,
    distance_b_m: 0.5,
}];

#[test]
fn native_three_dimensional_material_feedback_heats_local_cells_and_continues() {
    for r0 in [1.0, 2.0] {
        for refinement in [1usize, 2] {
            let (a, ca) = material(r0, 1.0, 0.002);
            let (b, cb) = material(3.0, 2.0, 0.0);
            let cards = [a, b];
            let curves = [ca, cb];
            let context = point(300.0);
            let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
            let cells = bindings(&cards, &curves, &context, &key);
            let mut thermal = initial(&curves, 999.0, None);
            let sources = cards.each_ref().map(|card| {
                resolve_conductor_state_with_basis(
                    card,
                    &context,
                    MaterialPropertySelection::SingleClaimOnly,
                    ElectricalPropertyBasis::Resistivity,
                )
                .unwrap()
            });
            let xs = (0..=2 * refinement)
                .map(|i| i as f64 / refinement as f64)
                .collect::<Vec<_>>();
            let (positions, tets) = slab(&xs);
            let mut problem = SpatialConduction::new(
                positions,
                tets,
                (0..12 * refinement)
                    .map(|i| sources[i / (6 * refinement)].clone())
                    .collect(),
                limits(),
                &mut || false,
            )
            .unwrap();
            let regions = (0..12 * refinement)
                .map(|tetrahedron| SpatialConductionRegion {
                    tetrahedron,
                    cell: if tetrahedron < 6 * refinement {
                        CellId(30)
                    } else {
                        CellId(10)
                    },
                })
                .collect::<Vec<_>>();
            let network = SpatialConductionThermalNetwork {
                cells: &cells,
                links: &LINKS,
                regions: &regions,
                volume_tolerance_m3: 1e-12,
            };
            let dt = 0.5;
            let current = 2.0;
            let mut theta = [0.0; 2];
            for step in 1..=2 {
                let before = thermal.energy_j();
                let drive = boundary(&problem, 2.0, current, None);
                let next = network
                    .step_candidate(
                        &thermal,
                        step as f64 * dt,
                        &problem,
                        &drive,
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0; 2],
                        || false,
                    )
                    .unwrap();
                // Independent two-cell BE elimination, C=[1,2], G=1, I=2:
                // (C1+dtG-dt I² rho0 alpha)theta1-dtG theta2=C1 theta1_old+dt I²rho0.
                let aa = 1.5 - 2.0 * r0 * 0.002;
                let bb = 2.5;
                let det = aa * bb - 0.25;
                let rhs = [theta[0] + 2.0 * r0, 2.0 * theta[1] + 6.0];
                theta = [
                    (bb * rhs[0] + 0.5 * rhs[1]) / det,
                    (0.5 * rhs[0] + aa * rhs[1]) / det,
                ];
                near(
                    next.thermal().temperature_k(CellId(30)).unwrap(),
                    300.0 + theta[0],
                );
                near(
                    next.thermal().temperature_k(CellId(10)).unwrap(),
                    300.0 + theta[1],
                );
                let rho = [r0 * (1.0 + 0.002 * theta[0]), 3.0];
                let voltage = current * (rho[0] + rho[1]);
                let report = next.report();
                for (i, j) in report
                    .electrical
                    .field
                    .current_density_a_m2
                    .iter()
                    .enumerate()
                {
                    near(j[0], -current);
                    near(j[1], 0.0);
                    near(j[2], 0.0);
                    near(
                        report.electrical.field.cell_power_w[i],
                        current * current * rho[i / (6 * refinement)] / (6 * refinement) as f64,
                    );
                }
                near(report.electrical.supplied_work_j, dt * current * voltage);
                near(
                    next.thermal().energy_j() - before,
                    report.electrical.supplied_work_j,
                );
                near(report.deposited_joule_j, report.electrical.heat_j);
                near(report.total_energy_residual_j, 0.0);
                // The independently solved terminal voltage gives the same
                // nonlinear endpoint under voltage drive and reversed current.
                if step == 1 {
                    for drive in [
                        boundary(&problem, 2.0, 0.0, Some(voltage)),
                        boundary(&problem, 2.0, -current, None),
                    ] {
                        let parity = network
                            .step_candidate(
                                &thermal,
                                dt,
                                &problem,
                                &drive,
                                limits(),
                                field_controls(),
                                config(),
                                &[0.0; 2],
                                || false,
                            )
                            .unwrap();
                        near(parity.thermal().energy_j(), next.thermal().energy_j());
                    }
                }
                let (t, p, _) = next.into_parts();
                thermal = t;
                problem = p;
            }
        }
    }
}

fn phase_material(
    basis: ElectricalPropertyBasis,
    slope: f64,
) -> (
    MaterialCard,
    EquilibriumEnthalpyPhaseCurve,
    [PropertyKey; 3],
) {
    let keys =
        std::array::from_fn(|i| PropertyKey::new(format!("principal_{i}"), basis.dimensions()));
    let mut claims = ClaimSet::new();
    for (i, key) in keys.iter().enumerate() {
        let value = if i == 0 {
            PropertyValue::Surface2D {
                axes: ["T".into(), "f".into()],
                axis_dims: [Temperature::DIMS, Dims::NONE],
                samples: [[300.0, 0.0], [600.0, 0.0], [300.0, 1.0], [600.0, 1.0]]
                    .map(|p| (p, 2.0 + slope * p[1] + 0.005 * (p[0] - 400.0)))
                    .to_vec(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: basis.dimensions(),
            }
        } else {
            PropertyValue::Scalar {
                value: 3.0 + i as f64,
                dims: basis.dimensions(),
            }
        };
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value,
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        300.0,
                        600.0,
                    )
                    .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0),
                interpolation: if i == 0 {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic confined-phase 3D transport".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "confined".into(),
            process: "authored".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-100.0, 300.0, 0.0),
            (0.0, 400.0, 0.0),
            (100.0, 400.0, 1.0),
            (300.0, 600.0, 1.0),
        ]
        .map(|(h, t, f)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: 1.0,
        })
        .to_vec(),
    )
    .unwrap();
    (card, curve, keys)
}
fn phase_initial(curve: &EquilibriumEnthalpyPhaseCurve) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(0),
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(0),
                mass_kg: 1.0,
                specific_h_j_kg: 20.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: common::boundary(),
            }),
        }],
        &[(MaterialId(0), curve as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 2,
            max_events: 2,
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
    .0
}

#[test]
fn original_phase_surface_changes_three_dimensional_scalar_and_tensor_current_heat() {
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        for slope in [1.0, 4.0] {
            let (card, curve, keys) = phase_material(basis, slope);
            let context = QueryPoint::new()
                .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 400.0)
                .unwrap()
                .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.2)
                .unwrap();
            let requirements = keys
                .iter()
                .map(|key| {
                    ScalarPropertyRequirement::try_with_key(
                        key,
                        ScalarAdmissibility::StrictlyPositive,
                    )
                    .unwrap()
                })
                .collect::<Vec<_>>();
            let snapshot = resolve_material_state_point(
                &card,
                &context,
                &requirements,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            let scalar = ResolvedConductorState::try_from_material_with_property(
                snapshot.clone(),
                basis,
                &keys[0],
            )
            .unwrap()
            .with_liquid_mass_fraction_axis("f")
            .unwrap();
            let oriented = ResolvedOrthotropicConductor3::try_from_phase_material(
                snapshot,
                basis,
                keys.clone(),
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                "f",
            )
            .unwrap();
            assert!(
                !oriented.has_admitted_static_sources(),
                "induction keeps its original stricter source boundary"
            );
            let thermal = phase_initial(&curve);
            let cells = [MaterialCell {
                cell: CellId(0),
                material_id: MaterialId(0),
                curve: &curve,
                point: ThermalMaterialPoint {
                    material: &card,
                    context: &context,
                    temperature_axis: "T",
                    liquid_fraction_axis: Some("f"),
                },
                conductivity: None,
            }];
            let regions = (0..6)
                .map(|tetrahedron| SpatialConductionRegion {
                    tetrahedron,
                    cell: CellId(0),
                })
                .collect::<Vec<_>>();
            let network = SpatialConductionThermalNetwork {
                cells: &cells,
                links: &[],
                regions: &regions,
                volume_tolerance_m3: 1e-12,
            };
            let (positions, tets) = slab(&[0.0, 1.0]);
            let scalar = SpatialConduction::new(
                positions.clone(),
                tets.clone(),
                vec![scalar; 6],
                limits(),
                &mut || false,
            )
            .unwrap();
            let tensor =
                SpatialConduction::new(positions, tets, vec![oriented; 6], limits(), &mut || false)
                    .unwrap();
            let supply = 0.5 * 4.0;
            let a = 2.0 + 0.2 * slope;
            let heat = if basis == ElectricalPropertyBasis::Resistivity {
                supply * a / (1.0 - supply * slope / 100.0)
            } else {
                2.0 * supply / (a + (a * a + 4.0 * slope * supply / 100.0).sqrt())
            };
            let boundary = boundary(&scalar, 1.0, 2.0, None);
            let [positive, negative] = electrodes(&scalar, 1.0);
            let next = network
                .step_candidate(
                    &thermal,
                    0.5,
                    &scalar,
                    &boundary,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            let tensor_next = network
                .step_candidate(
                    &thermal,
                    0.5,
                    &tensor,
                    &boundary,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            near(next.thermal().energy_j() - thermal.energy_j(), heat);
            near(tensor_next.thermal().energy_j(), next.thermal().energy_j());
            near(next.thermal().temperature_k(CellId(0)).unwrap(), 400.0);
            near(
                tensor_next.problem().conductors()[0]
                    .material()
                    .query()
                    .axes()["f"],
                0.2 + heat / 100.0,
            );
            near(next.report().electrical.heat_j, heat);
            near(next.report().total_energy_residual_j, 0.0);
            // A terminal current must rebuild the resistance at the actual
            // latent fraction, even when temperature stays exactly constant.
            for current_a in [2.0, -2.0] {
                let drive = SpatialConductionDrive::total_current(&positive, &negative, current_a);
                let current_scalar = network
                    .step_candidate(
                        &thermal,
                        0.5,
                        &scalar,
                        drive,
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0],
                        || false,
                    )
                    .unwrap();
                let current_tensor = network
                    .step_candidate(
                        &thermal,
                        0.5,
                        &tensor,
                        drive,
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0],
                        || false,
                    )
                    .unwrap();
                for candidate in [current_scalar.report(), current_tensor.report()] {
                    near(candidate.electrical.heat_j, heat);
                    near(candidate.electrical.supplied_work_j, heat);
                    near(candidate.total_energy_residual_j, 0.0);
                    for density in &candidate.electrical.field.current_density_a_m2 {
                        near(density[0], -current_a);
                        near(density[1], 0.0);
                        near(density[2], 0.0);
                    }
                }
                near(
                    current_scalar.thermal().energy_j(),
                    next.thermal().energy_j(),
                );
                near(
                    current_tensor.thermal().energy_j(),
                    next.thermal().energy_j(),
                );
                near(
                    current_tensor.problem().conductors()[0]
                        .material()
                        .query()
                        .axes()["f"],
                    0.2 + heat / 100.0,
                );
            }
        }
    }
}

#[test]
fn cancelled_final_acceptance_budgets_and_foreign_source_leave_both_histories_intact() {
    let (a, ca) = material(1.0, 1.0, 0.002);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &key);
    let thermal = initial(&curves, 999.0, None);
    let sources = cards.each_ref().map(|c| {
        resolve_conductor_state_with_basis(
            c,
            &context,
            MaterialPropertySelection::SingleClaimOnly,
            ElectricalPropertyBasis::Resistivity,
        )
        .unwrap()
    });
    let (positions, tets) = slab(&[0.0, 1.0, 2.0]);
    let problem = SpatialConduction::new(
        positions,
        tets,
        (0..12).map(|i| sources[i / 6].clone()).collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let regions = (0..12)
        .map(|tetrahedron| SpatialConductionRegion {
            tetrahedron,
            cell: if tetrahedron < 6 {
                CellId(30)
            } else {
                CellId(10)
            },
        })
        .collect::<Vec<_>>();
    let network = SpatialConductionThermalNetwork {
        cells: &cells,
        links: &LINKS,
        regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let drive = boundary(&problem, 2.0, 2.0, None);
    let saved = thermal.cells().clone();
    let mut polls = 0;
    let success = network
        .step_candidate(
            &thermal,
            0.5,
            &problem,
            &drive,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut actual = 0;
    assert!(matches!(
        network.step_candidate(
            &thermal,
            0.5,
            &problem,
            &drive,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || {
                actual += 1;
                actual == polls
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let mut limited = config();
    limited.maximum_circuit_trials = success.report().circuit_trials - 1;
    assert!(matches!(
        network.step_candidate(
            &thermal,
            0.5,
            &problem,
            &drive,
            limits(),
            field_controls(),
            limited,
            &[0.0; 2],
            || false
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let retry = network
        .step_candidate(
            &thermal,
            0.5,
            &problem,
            &drive,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    assert_eq!(retry.thermal().cells(), success.thermal().cells());
    let wrong = problem
        .with_conductors(vec![sources[1].clone(); 12], limits(), &mut || false)
        .unwrap();
    assert!(
        network
            .step_candidate(
                &thermal,
                0.5,
                &wrong,
                &drive,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    assert_eq!(thermal.cells(), &saved);
    assert_eq!(thermal.time_s(), 0.0);
    let mut bad_regions = regions.clone();
    bad_regions[0].cell = CellId(10);
    assert!(
        SpatialConductionThermalNetwork {
            regions: &bad_regions,
            ..network
        }
        .step_candidate(
            &thermal,
            0.5,
            &problem,
            &drive,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false
        )
        .is_err()
    );
}
