//! Independent electrical/mechanical and finite-cell heat balances. Authored
//! material data exercise the numerical model, not empirical qualification.
use super::*;
use fs_couple::circuit_enthalpy::material::{
    CellResistor, ConductionLink, MaterialCell, MaterialNetwork,
};
use fs_couple::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpySolve, ResistorCell};
use fs_couple::material_circuit::MaterialCircuitSolve;
use fs_couple::piezoelectric::spatial::transient::load::thermal::{
    LoadThermalError, ThermalElectrodeLoad,
};
use fs_matdb::SelectionPolicy;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Density, QuantitySpec, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

const CELLS: [CellId; 2] = [CellId(40), CellId(10)];
const LOAD: ThermalElectrodeLoad = ThermalElectrodeLoad {
    resistor_branch: 17,
    electrodes: [1, 0],
};
const AREA: f64 = 1e-6;
const VOLUME: f64 = LENGTH * AREA;
const RESISTANCE: f64 = 1e4;
const CAPACITANCE: f64 = EPSILON[2] * LENGTH;

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn thermal_material(
    cp: f64,
    k: f64,
    alpha: f64,
    electrical_upper: f64,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = ClaimSet::new();
    let rho = RESISTANCE * AREA / LENGTH;
    for (name, dims, low, high, upper) in [
        ("density", Density::DIMS, 1.0, 1.0, 400.0),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            cp,
            cp,
            400.0,
        ),
        (
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            k,
            k,
            400.0,
        ),
        (
            ELECTRICAL_RESISTIVITY_PROPERTY,
            ELECTRICAL_RESISTIVITY_DIMS,
            rho * (1.0 - 20.0 * alpha),
            rho * (1.0 + 100.0 * alpha),
            electrical_upper,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, low), (400.0, high)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    280.0,
                    upper,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored finite spatial load fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal load".into(),
            phase: "solid".into(),
            process: "uniform conductor".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(280.0),
            upper: &point(400.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-8,
    )
    .unwrap();
    (card, curve)
}
fn initial(curves: &[ResolvedHeatCapacityCurve; 2]) -> DepositionSimulation<'_> {
    let mut events = Vec::new();
    for (i, id) in CELLS.into_iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: VOLUME,
                specific_h_j_kg: 20.0 * (i + 1) as f64,
                geometry: CellGeometry {
                    start_m: [i as f64 * LENGTH, 0.0, 0.0],
                    end_m: [(i + 1) as f64 * LENGTH, 0.0, 0.0],
                    volume_m3: VOLUME,
                    surface_area_m2: 6.0 * AREA,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 9,
            a: CELLS[0],
            b: CELLS[1],
            area_m2: AREA,
            // The prepared source must replace this deliberately incorrect value.
            conductance_w_k: 999.0,
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1e-7,
            energy_tolerance_j: 1e-16,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn material_cells<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    context: &'a QueryPoint,
    key: &'a PropertyKey,
) -> [MaterialCell<'a>; 2] {
    std::array::from_fn(|i| MaterialCell {
        cell: CELLS[i],
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
}
fn thermal_resistor() -> CellResistor {
    CellResistor {
        destination: ResistorCell {
            branch: LOAD.resistor_branch,
            cell: CELLS[0],
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: LENGTH,
        area_m2: AREA,
    }
}
fn thermal_controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 2,
        maximum_circuit_trials: 256,
        maximum_newton_iterations: 16,
        difference_step_j: 1e-11,
        cell_energy_tolerance_j: 1e-16,
        thermal_energy_tolerance_j: 1e-16,
        total_energy_tolerance_j: 1e-16,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-16,
            power_tolerance_w: 1e-10,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-17,
            relative_tolerance: 1e-14,
            forcing_minimum: 1e-4,
            linear_restart: 2,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn field_controls() -> Controls {
    Controls {
        energy_tolerance_j: 1e-16,
        ..controls()
    }
}

// Independently eliminate one mechanical coordinate and one electrode charge:
// M vdot + K u + G V = 0; Q = C V - G u; Qdot = -V/R(T).
// The fixed fixture uses G=0 and u=v=0, recovering an ordinary RC midpoint step.
fn electrical(old: [f64; 3], mechanical: Option<[f64; 3]>, dt: f64, resistance: f64) -> [f64; 5] {
    let [u, v, q] = old;
    let (a, b, g) = mechanical.map_or((0.0, 0.0, 0.0), |[m, k, g]| {
        let stiffness = 4.0 * m / (dt * dt) + k;
        (
            (4.0 * m * u / (dt * dt) + 2.0 * m * v / dt) / stiffness,
            g / stiffness,
            g,
        )
    });
    let voltage = (q + g * a) / (CAPACITANCE + dt / (2.0 * resistance) + g * b);
    let middle_u = a - b * voltage;
    [
        2.0 * middle_u - u,
        4.0 * (middle_u - u) / dt - v,
        q - dt * voltage / resistance,
        voltage,
        dt * voltage * voltage / resistance,
    ]
}
fn expected(
    old: [f64; 3],
    temperatures: [f64; 2],
    mechanics: Option<[f64; 3]>,
    dt: f64,
    alpha: f64,
) -> ([f64; 2], [f64; 5]) {
    let conductance = AREA / (0.5 * LENGTH / 2.0 + 0.5 * LENGTH / 4.0);
    let transfer = dt * conductance;
    let neighbour = |t| (2.0 * VOLUME * temperatures[1] + transfer * t) / (2.0 * VOLUME + transfer);
    let motion = |t| electrical(old, mechanics, dt, RESISTANCE * (1.0 + alpha * (t - 300.0)));
    let residual =
        |t| VOLUME * (t - temperatures[0]) + transfer * (t - neighbour(t)) - motion(t)[4];
    let (mut low, mut high) = (280.0, 400.0);
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..80 {
        let middle = 0.5 * (low + high);
        if residual(middle) < 0.0 {
            low = middle;
        } else {
            high = middle;
        }
    }
    let t = 0.5 * (low + high);
    ([t, neighbour(t)], motion(t))
}

fn coupled_run(alpha: f64, moving: bool, sign: f64) -> (f64, f64) {
    coupled_run_at(alpha, moving, sign, &[5e-8, 1e-7, 1.5e-7, 2e-7])
}

fn coupled_run_at(alpha: f64, moving: bool, sign: f64, targets: &[f64]) -> (f64, f64) {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let node = *cube.top.last().unwrap();
        let dof = 3 * node + 2;
        let fixed = (0..3 * n)
            .filter(|&i| !moving || i != dof)
            .map(|i| (i, 0.0))
            .collect::<Vec<_>>();
        let zero = vec![[0.0; 3]; n];
        let nodal_charge = vec![0.0; n];
        let mut displacement = zero.clone();
        let mut old = if moving {
            [sign * 1e-7, 0.0, 0.0]
        } else {
            [0.0, 0.0, sign * 100.0 * CAPACITANCE]
        };
        displacement[node][2] = old[0];
        let mut body = Body::new(
            &model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Charge(old[2])),
            &displacement,
            &zero,
            &nodal_charge,
            field_controls(),
            cx,
        )
        .unwrap();
        let mechanics = moving.then(|| {
            let diagonal = |matrix: &fs_sparse::Csr| {
                let (cols, values) = matrix.row(dof);
                values[cols.iter().position(|&i| i == dof).unwrap()]
            };
            let (cols, values) = model.coupling().row(dof);
            [
                diagonal(&model.solid().mass),
                diagonal(&model.solid().stiffness),
                cols.iter()
                    .zip(values)
                    .filter(|(col, _)| cube.top.contains(col))
                    .map(|(_, value)| value)
                    .sum(),
            ]
        });
        let (a, ca) = thermal_material(1.0, 2.0, alpha, 400.0);
        let (b, cb) = thermal_material(2.0, 4.0, 0.0, 400.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let context = point(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = material_cells(&cards, &curves, &context, &key);
        let resistors = [thermal_resistor()];
        let links = [ConductionLink {
            link: 9,
            distance_a_m: 0.5 * LENGTH,
            distance_b_m: 0.5 * LENGTH,
        }];
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &links,
        };
        let mut thermal = initial(&curves);
        let mut temperatures = [300.0; 2];
        let prepared = network.prepare(&thermal, || false).unwrap();
        let initial_energy = body.state().total_energy_j() + thermal.energy_j();
        let mut supplied_heat = 0.0;
        for &target in targets {
            let dt = target - thermal.time_s();
            let (expected_t, next) = expected(old, temperatures, mechanics, dt, alpha);
            let candidate = body
                .step_with_thermal_loads_candidate(
                    &thermal,
                    target,
                    &prepared,
                    &[LOAD],
                    &[Drive::Voltage(0.0), Drive::Current(0.0)],
                    &zero,
                    &nodal_charge,
                    field_controls(),
                    thermal_controls(),
                    &[0.0; 2],
                    cx,
                )
                .unwrap();
            assert_eq!(candidate.thermal().time_s(), target);
            assert_eq!(candidate.body().state().time_s(), target);
            let report = candidate.report();
            assert_eq!(report.electrical.native.time_step_s, dt);
            let load = &report.electrical.loads.loads[0];
            for (cell, t) in CELLS.into_iter().zip(expected_t) {
                near(
                    candidate.thermal().temperature_k(cell).unwrap(),
                    t,
                    2e-6,
                    0.0,
                );
            }
            near(
                load.load.conductor.resistance_ohm(),
                RESISTANCE * (1.0 + alpha * (expected_t[0] - 300.0)),
                1e-3,
                0.0,
            );
            assert_eq!(
                load.load.conductor.material().card_identity(),
                cards[0].content_hash()
            );
            assert_eq!(
                report.electrical.destinations,
                vec![resistors[0].destination]
            );
            near(
                candidate.body().state().displacement_m()[node][2],
                next[0],
                2e-14,
                1e-8,
            );
            near(
                candidate.body().state().velocity_m_s()[node][2],
                next[1],
                2e-7,
                1e-8,
            );
            near(
                candidate.body().state().electrodes()[1].charge_c,
                next[2],
                2e-17,
                1e-8,
            );
            near(load.joule_heat_j, next[4], 2e-16, 1e-8);
            near(report.deposited_joule_j, next[4], 2e-16, 1e-8);
            near(report.non_electrical_heat_j, 0.0, 2e-16, 0.0);
            near(
                candidate.thermal().energy_j() + candidate.body().state().total_energy_j(),
                initial_energy,
                4e-16,
                0.0,
            );
            // Newton probes have not advanced either accepted input.
            near(body.state().electrodes()[1].charge_c, old[2], 2e-17, 1e-8);
            near(
                thermal.temperature_k(CELLS[0]).unwrap(),
                temperatures[0],
                2e-6,
                0.0,
            );
            supplied_heat += report.deposited_joule_j;
            (thermal, body, _) = candidate.into_parts();
            temperatures = expected_t;
            old.copy_from_slice(&next[..3]);
        }
        assert!(supplied_heat > 0.0);
        assert!(thermal.temperature_k(CELLS[1]).unwrap() > 300.0);
        near(
            thermal.energy_j() - 60.0 * VOLUME,
            supplied_heat,
            4e-16,
            0.0,
        );
        (
            body.state().electrodes()[1].voltage_v,
            body.state().displacement_m()[node][2],
        )
    })
}

#[test]
fn g1_spatial_load_heat_flows_to_finite_neighbour_and_changes_rc_charge_decay() {
    let constant = coupled_run(0.0, false, 1.0);
    let heated = coupled_run(0.02, false, 1.0);
    let reverse = coupled_run(0.02, false, -1.0);
    assert!(heated.0 > constant.0 + 1.0);
    near(reverse.0, -heated.0, 2e-6, 0.0);
}

#[test]
fn g1_spatial_motion_and_temperature_dependent_shunt_share_the_same_implicit_endpoint() {
    let constant = coupled_run(0.0, true, 1.0);
    let heated = coupled_run(0.02, true, 1.0);
    assert!((heated.1 - constant.1).abs() > 1e-12);
}

#[test]
fn g1_nonuniform_thermal_targets_keep_exact_shared_clock_and_physical_motion() {
    let targets = [3e-8_f64, 9e-8, 1.2e-7];
    // This pair actually exhibits the one-ULP absolute-time drift. Continue
    // through it while checking the same independent motion and heat balances.
    assert_ne!(targets[0] + (targets[1] - targets[0]), targets[1]);
    coupled_run_at(0.02, true, 1.0, &targets);
}

#[test]
fn g4_late_source_refusal_and_budget_leave_both_histories_retryable() {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let fixed = (0..3 * n).map(|i| (i, 0.0)).collect::<Vec<_>>();
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let body = Body::new(
            &model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Charge(100.0 * CAPACITANCE)),
            &zero,
            &zero,
            &charge,
            field_controls(),
            cx,
        )
        .unwrap();
        let (a, ca) = thermal_material(1.0, 2.0, 0.02, 305.0);
        let (b, cb) = thermal_material(2.0, 4.0, 0.0, 400.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let context = point(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = material_cells(&cards, &curves, &context, &key);
        let resistors = [thermal_resistor()];
        let links = [ConductionLink {
            link: 9,
            distance_a_m: 0.5 * LENGTH,
            distance_b_m: 0.5 * LENGTH,
        }];
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &links,
        };
        let thermal = initial(&curves);
        let prepared = network.prepare(&thermal, || false).unwrap();
        let advance = |target, bindings: &[ThermalElectrodeLoad], solve| {
            body.step_with_thermal_loads_candidate(
                &thermal,
                target,
                &prepared,
                bindings,
                &[Drive::Voltage(0.0), Drive::Current(0.0)],
                &zero,
                &charge,
                field_controls(),
                solve,
                &[0.0; 2],
                cx,
            )
        };
        assert!(matches!(
            advance(5e-8, &[LOAD], thermal_controls()),
            Err(CircuitEnthalpyError::Source(LoadThermalError::Material(_)))
        ));
        assert!(matches!(
            advance(
                1e-9,
                &[LOAD],
                CircuitEnthalpySolve {
                    maximum_circuit_trials: 1,
                    ..thermal_controls()
                }
            ),
            Err(CircuitEnthalpyError::CircuitBudget)
        ));
        for bindings in [
            vec![],
            vec![LOAD, LOAD],
            vec![ThermalElectrodeLoad {
                resistor_branch: 18,
                ..LOAD
            }],
        ] {
            assert!(advance(1e-9, &bindings, thermal_controls()).is_err());
        }
        let retry = advance(1e-9, &[LOAD], thermal_controls()).unwrap();
        assert!(retry.thermal().temperature_k(CELLS[0]).unwrap() > 300.0);
        assert!(
            retry.body().state().electrodes()[1].charge_c < body.state().electrodes()[1].charge_c
        );
        assert!(
            body.step_with_thermal_loads_candidate(
                retry.thermal(),
                2e-9,
                &prepared,
                &[LOAD],
                &[Drive::Voltage(0.0), Drive::Current(0.0)],
                &zero,
                &charge,
                field_controls(),
                thermal_controls(),
                &[0.0; 2],
                cx
            )
            .is_err()
        );
        gate.request();
        assert!(matches!(
            advance(1e-9, &[LOAD], thermal_controls()),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(body.state().time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        near(
            body.state().electrodes()[1].charge_c,
            100.0 * CAPACITANCE,
            1e-20,
            0.0,
        );
        assert_eq!(thermal.temperature_k(CELLS[0]).unwrap(), 300.0);
    });
}
