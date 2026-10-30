//! G1/G3/G4: physical tensor memory and finite self-heating close together.
mod affine;
mod expansion;
use super::*;
use fs_couple::tensor_material::solid::source::SolidMaterialSource;
use fs_couple::tensor_material::solid::thermal::relaxation::{
    RelaxationThermalCell, RelaxationThermalSolve, RelaxingSolid,
};
use fs_matdb::SelectionPolicy;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use std::collections::BTreeMap;
use std::convert::Infallible;

const FORCE: [[f64; 3]; 4] = [[0.0; 3], [0.0; 3], [0.0; 3], [0.0, 0.0, 0.02]];
const CP: f64 = 0.0001;

fn thermal_card(variable_tau: bool, cp: f64) -> (MaterialCard, TensorPronyProperties) {
    let (mut claims, properties) = fixture(1000.0, 10.0, 0.4);
    claims.push(scalar(
        "specific_heat_capacity",
        Dims([2, 0, -2, -1, 0, 0]),
        cp,
    ));
    if variable_tau {
        let tau = claims
            .iter_mut()
            .find(|c| c.key.name() == "branch_tau")
            .unwrap();
        tau.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: fs_qty::Temperature::DIMS,
            knots: vec![
                (280.0, 0.08),
                (300.0, 0.02),
                (320.0, 0.002),
                (330.0, 0.0005),
            ],
            dims: Time::DIMS,
        };
        tau.interpolation = InterpolationPolicy::LinearInside;
    }
    (card(claims), properties)
}
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0])),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &query().with("T", 280.0).unwrap(),
            upper: &query().with("T", 330.0).unwrap(),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-10,
    )
    .unwrap()
}
fn bind(card: &MaterialCard, properties: TensorPronyProperties) -> ResolvedTensorMaterial {
    ResolvedTensorMaterial::resolve(
        card,
        &query(),
        properties,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn controls() -> RelaxationThermalSolve {
    RelaxationThermalSolve {
        maximum_cells: 4,
        maximum_trials: 512,
        maximum_iterations: 24,
        difference_step_j: 1e-7,
        cell_energy_tolerance_j: 1e-10,
        thermal_energy_tolerance_j: 1e-10,
        total_energy_tolerance_j: 2e-10,
        inventory_relative_tolerance: 1e-12,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-14,
            linear_restart: 4,
            forcing_minimum: 1e-4,
            ..NewtonKrylovConfig::default()
        },
        projection: SolidProjectionBudget {
            max_coordinates: 2,
            max_basis_scalars: 4,
            max_material_points: 2,
            max_scalar_terms: 100_000,
        },
        dynamics: dynamics(),
    }
}
fn elements(source: &ResolvedTensorMaterial, cx: &Cx<'_>) -> ViscoElements {
    ViscoElements::new(
        SolidMeshRef::Tet4 {
            nodes_m: &NODES,
            elements: &[[0, 1, 2, 3]],
        },
        TetViscoField::Uniform(&TetViscoMaterial::from_source(source, IDENTITY).unwrap()),
        &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 128,
        },
        cx,
    )
    .unwrap()
}
fn body<'s>(
    elements: &'s ViscoElements,
    source: ResolvedTensorMaterial,
    cx: &Cx<'_>,
    gate: &CancelGate,
) -> TensorSolidBody<'s> {
    TensorSolidBody::from_material_sources(
        elements
            .project_basis(&[vec![1.0]], controls().projection, cx)
            .unwrap(),
        vec![SolidMaterialSource {
            material: source,
            material_to_mesh: IDENTITY,
        }],
        dynamics(),
        gate,
    )
    .unwrap()
}
fn heat<'a>(
    curves: &[&'a ResolvedHeatCapacityCurve],
    capacity: &[f64],
    conductance: f64,
) -> DepositionSimulation<'a> {
    let mut schedule: Vec<_> = curves
        .iter()
        .enumerate()
        .map(|(i, _)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                material_id: MaterialId(i as u64),
                segment_id: i as u64,
                body_id: i as u64,
                mass_kg: 10.0 / 6.0,
                specific_h_j_kg: 20.0 * capacity[i],
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: 1.0 / 6.0,
                    surface_area_m2: 3.0,
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
    if curves.len() == 2 {
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 55,
                a: CellId(0),
                b: CellId(1),
                area_m2: 0.1,
                conductance_w_k: conductance,
            }),
        });
    }
    let laws: Vec<_> = curves
        .iter()
        .enumerate()
        .map(|(i, c)| (MaterialId(i as u64), c.curve() as &dyn ThermalLaw))
        .collect();
    DepositionSimulation::new(
        schedule,
        &laws,
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 16,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.03,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn cell<'a>(
    id: u64,
    card: &'a MaterialCard,
    caloric: &'a ResolvedHeatCapacityCurve,
    q: &'a QueryPoint,
) -> RelaxationThermalCell<'a> {
    RelaxationThermalCell {
        cell: CellId(id),
        caloric,
        temperature_interval_k: [280.0, 330.0],
        point: ThermalMaterialPoint {
            material: card,
            context: q,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
    }
}
fn no_extra(_: &DepositionSimulation<'_>) -> Result<StepCoefficients<()>, Infallible> {
    Ok(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    })
}
fn tau(t: f64) -> f64 {
    if t <= 300.0 {
        0.08 + (t - 280.0) * (0.02 - 0.08) / 20.0
    } else if t <= 320.0 {
        0.02 + (t - 300.0) * (0.002 - 0.02) / 20.0
    } else {
        0.002 + (t - 320.0) * (0.0005 - 0.002) / 10.0
    }
}

// Independent physical-strain equations: no factors, PHS matrices or source
// projection are used. The native endpoint scheme is closed with bisection.
fn reference(mut y: [f64; 5], relaxation: f64, dt: f64) -> [f64; 5] {
    let rhs = |z: [f64; 5]| {
        [
            z[1],
            (0.02 - (1000.0 * z[0] + 280.0 * z[2]) / 6.0) / (10.0 / 60.0),
            z[1] - z[2] / relaxation,
            280.0 * z[2] * z[2] / (6.0 * relaxation),
            0.02 * z[1],
        ]
    };
    let h = dt / 2048.0;
    for _ in 0..2048 {
        let a = rhs(y);
        let b = rhs(core::array::from_fn(|i| y[i] + 0.5 * h * a[i]));
        let c = rhs(core::array::from_fn(|i| y[i] + 0.5 * h * b[i]));
        let d = rhs(core::array::from_fn(|i| y[i] + h * c[i]));
        y = core::array::from_fn(|i| y[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}
fn coupled_reference(y: [f64; 5], temperature: f64, dt: f64) -> ([f64; 5], f64) {
    let mut lo = temperature;
    let mut hi = 330.0;
    for _ in 0..52 {
        let t = 0.5 * (lo + hi);
        let next = reference(y, tau(t), dt);
        if (10.0 / 6.0) * CP * (t - temperature) - (next[3] - y[3]) > 0.0 {
            hi = t;
        } else {
            lo = t;
        }
    }
    let t = 0.5 * (lo + hi);
    (reference(y, tau(t), dt), t)
}

#[test]
fn finite_temperature_changes_actual_tensor_relaxation_and_continues_full_memory() {
    with_cx(|cx, gate| {
        let (card, properties) = thermal_card(true, CP);
        let curve = caloric(&card);
        let source = bind(&card, properties);
        let elements = elements(&source, cx);
        let body = body(&elements, source, cx, gate);
        let initial = body
            .history(&[0.006], &[0.08], InitialMemory::Instantaneous, gate)
            .unwrap();
        let initial_state = initial.material_history().state().to_vec();
        let context = query();
        let bindings = [cell(0, &card, &curve, &context)];
        let mut thermal = heat(&[&curve], &[CP], 0.0);
        let mut model: Option<RelaxingSolid<'_>> = None;
        let mut expected = [0.006, 0.08, 0.006, 0.0, 0.0];
        let mut temperature = 300.0;
        for target in [0.03, 0.06] {
            let history = model
                .as_ref()
                .map(RelaxingSolid::history)
                .unwrap_or_else(|| initial.clone());
            let (next, t) = coupled_reference(expected, temperature, 0.03);
            let result = history
                .advance_relaxation_thermal_candidate(
                    &thermal,
                    target,
                    &FORCE,
                    &bindings,
                    controls(),
                    &[0.001],
                    no_extra,
                    gate,
                )
                .unwrap();
            let actual = result.model().history();
            near(actual.material_history().displacement_m()[0], next[0], 2e-9);
            near(actual.material_history().velocity_m_s()[0], next[1], 2e-8);
            near(
                actual.material_history().dissipated_heat_j(),
                next[3],
                2e-11,
            );
            near(result.thermal().temperature_k(CellId(0)).unwrap(), t, 2e-7);
            assert_eq!(actual.material_history().parameter_work_j(), 0.0);
            assert_eq!(
                actual.material_history().time_s().to_bits(),
                target.to_bits()
            );
            assert!(result.report().total_energy_residual_j.abs() < 2e-10);
            assert!(result.report().external_thermal_heat_j.abs() < 1e-12);
            let source = &result.model().body().material_sources().unwrap()[0].material;
            near(source.state().state_coordinate("T").unwrap(), t, 2e-7);
            near(
                source.state().property("branch_tau").unwrap().value_si(),
                tau(t),
                1e-9,
            );
            assert!(result.report().mechanical_trials > 1);
            let (next_model, next_thermal, _) = result.into_parts();
            model = Some(next_model);
            thermal = next_thermal;
            expected = next;
            temperature = t;
        }
        assert_eq!(
            model.as_ref().unwrap().history().material_history().steps(),
            2
        );
        assert_eq!(initial.material_history().state(), initial_state);
        let frozen = reference([0.006, 0.08, 0.006, 0.0, 0.0], 0.02, 0.03);
        let changed = coupled_reference([0.006, 0.08, 0.006, 0.0, 0.0], 300.0, 0.03).0;
        assert!((frozen[0] - changed[0]).abs() > 1e-6);
    });
}

#[test]
fn neighboring_material_heat_capacity_changes_the_same_coupled_mechanical_trajectory() {
    with_cx(|cx, gate| {
        let mut results = Vec::new();
        for companion_capacity in [CP, 10.0 * CP] {
            let (a, pa) = thermal_card(true, CP);
            let (b, pb) = thermal_card(true, companion_capacity);
            let ca = caloric(&a);
            let cb = caloric(&b);
            let sources = [bind(&a, pa), bind(&b, pb)];
            let materials = sources
                .iter()
                .map(|s| TetViscoMaterial::from_source(s, IDENTITY).unwrap())
                .collect::<Vec<_>>();
            let mut nodes = NODES.to_vec();
            nodes.extend(NODES.map(|mut p| {
                p[0] += 2.0;
                p
            }));
            let fixed = (0..24).filter(|&n| n != 11 && n != 23).collect::<Vec<_>>();
            let elements = ViscoElements::new(
                SolidMeshRef::Tet4 {
                    nodes_m: &nodes,
                    elements: &[[0, 1, 2, 3], [4, 5, 6, 7]],
                },
                TetViscoField::PerElement(&materials),
                &fixed,
                TetAssemblyBudget::standard(),
                ViscoStateBudget {
                    max_internal_variables: 256,
                },
                cx,
            )
            .unwrap();
            let body = TensorSolidBody::from_material_sources(
                elements
                    .project_basis(&[vec![1.0, 0.0], vec![0.0, 1.0]], controls().projection, cx)
                    .unwrap(),
                sources
                    .into_iter()
                    .map(|material| SolidMaterialSource {
                        material,
                        material_to_mesh: IDENTITY,
                    })
                    .collect(),
                dynamics(),
                gate,
            )
            .unwrap();
            let history = body
                .history(
                    &[0.006, 0.0],
                    &[0.08, 0.0],
                    InitialMemory::Instantaneous,
                    gate,
                )
                .unwrap();
            let thermal = heat(&[&ca, &cb], &[CP, companion_capacity], 0.01);
            let context = query();
            let cells = [cell(0, &a, &ca, &context), cell(1, &b, &cb, &context)];
            let mut force = vec![[0.0; 3]; 8];
            force[3][2] = 0.02;
            let next = history
                .advance_relaxation_thermal_candidate(
                    &thermal,
                    0.03,
                    &force,
                    &cells,
                    controls(),
                    &[0.0004, 0.0004],
                    no_extra,
                    gate,
                )
                .unwrap();
            assert!(next.thermal().temperature_k(CellId(1)).unwrap() > 300.0);
            near(next.report().mechanical.element_heat_j[1], 0.0, 1e-14);
            near(
                next.thermal().energy_j() - thermal.energy_j(),
                next.report().deposited_material_heat_j,
                1e-10,
            );
            assert!(next.report().total_energy_residual_j.abs() < 2e-10);
            results.push((
                next.model().history().material_history().displacement_m()[0],
                next.thermal().temperature_k(CellId(0)).unwrap(),
            ));
        }
        assert!((results[0].0 - results[1].0).abs() > 1e-7);
        assert!((results[0].1 - results[1].1).abs() > 0.01);
    });
}

#[test]
fn unsupported_storage_laws_full_context_and_late_cancellation_refuse_atomically() {
    with_cx(|cx, gate| {
        for key in ["young_modulus", "branch_c00", "density", "branch_tau"] {
            let (mut claims, properties) = fixture(1000.0, 10.0, 0.4);
            let claim = claims.iter_mut().find(|c| c.key.name() == key).unwrap();
            let (value, dims) = match claim.value {
                PropertyValue::Scalar { value, dims } => (value, dims),
                _ => unreachable!(),
            };
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: vec![
                    (280.0, value),
                    (300.0, value),
                    (
                        310.0,
                        if key == "branch_tau" {
                            -value
                        } else {
                            2.0 * value
                        },
                    ),
                    (330.0, value),
                ],
                dims,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
            let card = card(claims);
            let source = bind(&card, properties);
            assert!(
                source
                    .admit_thermal_relaxation(&card, "T", [280.0, 330.0])
                    .is_err()
            );
        }
        let (card, properties) = thermal_card(true, CP);
        let curve = caloric(&card);
        let source = bind(&card, properties);
        let elements = elements(&source, cx);
        let body = body(&elements, source, cx, gate);
        let history = body
            .history(&[0.006], &[0.08], InitialMemory::Instantaneous, gate)
            .unwrap();
        let thermal = heat(&[&curve], &[CP], 0.0);
        let context = query();
        let cells = [cell(0, &card, &curve, &context)];
        let state = history.material_history().state().to_vec();
        let energy = thermal.energy_j();
        let good = history
            .advance_relaxation_thermal_candidate(
                &thermal,
                0.03,
                &FORCE,
                &cells,
                controls(),
                &[0.001],
                no_extra,
                gate,
            )
            .unwrap();
        let wrong_context = context.clone().with("omega", 11.0).unwrap();
        assert!(
            history
                .advance_relaxation_thermal_candidate(
                    &thermal,
                    0.03,
                    &FORCE,
                    &[cell(0, &card, &curve, &wrong_context)],
                    controls(),
                    &[0.001],
                    no_extra,
                    gate
                )
                .is_err()
        );
        assert!(
            history
                .advance_relaxation_thermal_candidate(
                    &thermal,
                    0.03,
                    &FORCE,
                    &cells,
                    controls(),
                    &[0.1],
                    no_extra,
                    gate
                )
                .is_err()
        );
        let later = CancelGate::new_clock_free();
        let failed = history.advance_relaxation_thermal_candidate(
            &thermal,
            0.03,
            &FORCE,
            &cells,
            controls(),
            &[0.001],
            |trial| {
                later.request();
                no_extra(trial)
            },
            &later,
        );
        assert!(matches!(
            failed,
            Err(fs_couple::circuit_enthalpy::CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(history.material_history().state(), state);
        assert_eq!(thermal.energy_j(), energy);
        assert_eq!(history.material_history().time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        let retry = history
            .advance_relaxation_thermal_candidate(
                &thermal,
                0.03,
                &FORCE,
                &cells,
                controls(),
                &[0.001],
                no_extra,
                gate,
            )
            .unwrap();
        assert_eq!(
            retry.model().history().material_history().state(),
            good.model().history().material_history().state()
        );
        assert_eq!(retry.thermal().energy_j(), good.thermal().energy_j());
    });
}
