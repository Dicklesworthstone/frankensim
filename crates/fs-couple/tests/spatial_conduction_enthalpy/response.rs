//! Independent nonlinear 3D current, caloric balance and continuation oracles.
#[path = "response/circuit.rs"]
mod circuit;
#[path = "response/mechanical.rs"]
mod mechanical;
#[path = "response/motion.rs"]
mod motion;
#[path = "response/power.rs"]
mod power;
#[path = "response/terminals.rs"]
mod terminals;
use super::*;
use fs_couple::circuit_enthalpy::conduction::ConductionThermalError;
use fs_em::spatial::conduction::response::{
    SpatialConductionResponse, SpatialConductionResponseControls, SpatialConductionResponseDrive,
    SpatialConductionResponseError,
};
use fs_material::conductor::response::{ConductorMagnitude, ScalarConductorResponse};

const CP: [f64; 2] = [2.0, 3.0];
const ALPHA: [f64; 2] = [0.002, 0.004];

fn query(t: f64, coordinate: ConductorMagnitude) -> QueryPoint {
    point(t)
        .with_quantity("operating", coordinate.quantity(), 0.0)
        .unwrap()
}

pub(super) fn controls() -> SpatialConductionResponseControls {
    SpatialConductionResponseControls {
        maximum_iterations: 100,
        maximum_field_solves: 300,
        absolute_conductivity_tolerance_s_m: 1e-11,
        relative_conductivity_tolerance: 1e-11,
        current_density_tolerance_a_m2: 1e-10,
    }
}

fn nonlinear_material(
    r: f64,
    cp: f64,
    alpha: f64,
    beta: f64,
    coordinate: ConductorMagnitude,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (base, _) = material(r, cp, 0.0);
    let mut claims = base.claims().clone();
    let basis = ElectricalPropertyBasis::Resistivity;
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("original_nonlinear_transport", basis.dimensions()),
            value: PropertyValue::Surface2D {
                axes: ["T".into(), "operating".into()],
                axis_dims: [Temperature::DIMS, coordinate.quantity().dims()],
                samples: [[200.0, 0.0], [2000.0, 0.0], [200.0, 20.0], [2000.0, 20.0]]
                    .into_iter()
                    .map(|p| (p, r + alpha * (p[0] - 300.0) + beta * p[1]))
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: basis.dimensions(),
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    200.0,
                    2000.0,
                )
                .with_quantity("operating", coordinate.quantity(), 0.0, 20.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured nonlinear tetrahedral electrothermal oracle".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &query(200.0, coordinate),
            upper: &query(2000.0, coordinate),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}

fn problem(cards: &[MaterialCard; 2], coordinate: ConductorMagnitude) -> SpatialConductionResponse {
    problem_at(cards, coordinate, &query(300.0, coordinate))
}

fn problem_at(
    cards: &[MaterialCard; 2],
    coordinate: ConductorMagnitude,
    context: &QueryPoint,
) -> SpatialConductionResponse {
    let basis = ElectricalPropertyBasis::Resistivity;
    let key = PropertyKey::new("original_nonlinear_transport", basis.dimensions());
    let sources: Vec<_> = (0..12)
        .map(|tet| {
            let card = &cards[tet / 6];
            let material = resolve_material_state_point(
                card,
                context,
                &[basis.requirement_for(&key).unwrap()],
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            let conductor =
                ResolvedConductorState::try_from_material_with_property(material, basis, &key)
                    .unwrap();
            ScalarConductorResponse::new(card, conductor, &[("operating", coordinate)]).unwrap()
        })
        .collect();
    let (positions, tets) = slab(&[0.0, 1.0, 2.0]);
    SpatialConductionResponse::new(positions, tets, sources, limits(), &mut || false)
        .unwrap()
        .with_relaxation(0.5)
        .unwrap()
}

pub(super) fn boundary_response(
    problem: &SpatialConductionResponse,
    current: f64,
    voltage: Option<f64>,
) -> SpatialConductionBoundary {
    let mut boundary = SpatialConductionBoundary::default();
    for face in problem.exterior_faces() {
        let xs = problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0]);
        if xs == [0.0; 3] {
            boundary.voltage.push(ConductionVoltageFace {
                face,
                voltage_v: [0.0; 3],
            });
        } else if xs == [2.0; 3] {
            if let Some(v) = voltage {
                boundary.voltage.push(ConductionVoltageFace {
                    face,
                    voltage_v: [v; 3],
                });
            } else {
                boundary.current.push(ConductionCurrentFace {
                    face,
                    inward_current_density_a_m2: current,
                });
            }
        }
    }
    boundary
}

fn regions() -> Vec<SpatialConductionRegion> {
    (0..12)
        .rev()
        .map(|tetrahedron| SpatialConductionRegion {
            tetrahedron,
            cell: if tetrahedron < 6 {
                CellId(30)
            } else {
                CellId(10)
            },
        })
        .collect()
}

// Independent two-cell elimination, G=1 W/K, each volume/mass=1.
// rho_i=r_i+alpha_i*theta_i+beta_i*|E_i| gives
// E_i=J*(r_i+alpha_i*theta_i)/(1-beta_i*J).
fn endpoint(
    j: f64,
    dt: f64,
    old: [f64; 2],
    r: [f64; 2],
    beta: [f64; 2],
    field_law: bool,
) -> ([f64; 2], [f64; 2]) {
    let d = beta.map(|b| if field_law { 1.0 - b * j } else { 1.0 });
    let base = std::array::from_fn::<_, 2, _>(|i| r[i] + if field_law { 0.0 } else { beta[i] * j });
    let a = std::array::from_fn::<_, 2, _>(|i| CP[i] + dt - dt * j * j * ALPHA[i] / d[i]);
    let rhs = std::array::from_fn::<_, 2, _>(|i| CP[i] * old[i] + dt * j * j * base[i] / d[i]);
    let det = a[0] * a[1] - dt * dt;
    let theta = [
        (a[1] * rhs[0] + dt * rhs[1]) / det,
        (a[0] * rhs[1] + dt * rhs[0]) / det,
    ];
    (
        theta,
        std::array::from_fn(|i| j * (base[i] + ALPHA[i] * theta[i]) / d[i]),
    )
}

fn voltage_reference(voltage: f64, old: [f64; 2], r: [f64; 2]) -> (f64, [f64; 2], [f64; 2]) {
    let (mut low, mut high) = (0.0, 4.0);
    for _ in 0..100 {
        let j = (low + high) / 2.0;
        let (_, e) = endpoint(j, 0.1, old, r, [0.05, 0.08], true);
        if e[0] + e[1] < voltage {
            low = j;
        } else {
            high = j;
        }
    }
    let j = (low + high) / 2.0;
    let (theta, e) = endpoint(j, 0.1, old, r, [0.05, 0.08], true);
    (j, theta, e)
}

#[test]
fn voltage_dependent_caloric_response_and_applicability_refuse_without_field_free_energy() {
    use fs_material::electrical_interface::response::ElectricalInterfaceCoordinate;
    let coordinate = ConductorMagnitude::ElectricField;
    let voltage = ElectricalInterfaceCoordinate::VoltageJump.quantity();
    let point = |temperature| {
        query(temperature, coordinate)
            .with_quantity("caloric_voltage", voltage, 0.0)
            .unwrap()
    };
    for response in [false, true] {
        let (a, _) = nonlinear_material(1.0, CP[0], ALPHA[0], 0.05, coordinate);
        let (b, _) = nonlinear_material(2.0, CP[1], ALPHA[1], 0.08, coordinate);
        let mut claims = ClaimSet::new();
        for (_, original) in a.claims().claims_ordered() {
            let mut claim = original.clone();
            if claim.key.name() == "specific_heat_capacity" {
                claim.validity =
                    claim
                        .validity
                        .with_quantity("caloric_voltage", voltage, 0.0, 20.0);
                if response {
                    claim.value = PropertyValue::Surface2D {
                        axes: ["T".into(), "caloric_voltage".into()],
                        axis_dims: [Temperature::DIMS, voltage.dims()],
                        samples: [[200.0, 0.0], [2000.0, 0.0], [200.0, 20.0], [2000.0, 20.0]]
                            .into_iter()
                            .map(|p| (p, CP[0] + 0.01 * p[1]))
                            .collect(),
                        triangles: vec![[0, 1, 2], [1, 3, 2]],
                        dims: SPECIFIC_HEAT_CAPACITY_DIMS,
                    };
                }
            }
            claims.insert_claim(claim).unwrap();
        }
        let a = MaterialCard::assemble(a.id().clone(), claims, vec![]).unwrap();
        let cards = [a, b];
        let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
        let density = PropertyKey::new("density", Density::DIMS);
        let curves = cards.each_ref().map(|card| {
            ResolvedHeatCapacityCurve::try_new(
                card,
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
            .unwrap()
        });
        let context = point(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = bindings(&cards, &curves, &context, &key);
        let regions = regions();
        let network = SpatialConductionThermalNetwork {
            cells: &cells,
            links: &LINKS,
            regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        let thermal = initial(&curves, 999.0, None);
        let problem = problem_at(&cards, coordinate, &context);
        let drive = boundary_response(&problem, 0.0, Some(6.0));
        assert!(matches!(
            network.step_response_candidate(
                &thermal,
                0.1,
                &problem,
                &drive,
                limits(),
                field_controls(),
                controls(),
                config(),
                &[0.0; 2],
                || false,
            ),
            Err(CircuitEnthalpyError::Source(
                ConductionThermalError::Invalid(
                    "caloric source response and applicability must be independent of electrical fields"
                )
            ))
        ));
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(
            problem.conductors()[0].material().query().axes()["caloric_voltage"],
            0.0
        );
    }
}

#[test]
fn nonlinear_tetrahedral_current_heats_local_finite_cells_and_continues() {
    let coordinate = ConductorMagnitude::ElectricField;
    let mut heats = Vec::new();
    for first in [1.0, 2.0] {
        let r = [first, 2.0];
        let (a, ca) = nonlinear_material(r[0], CP[0], ALPHA[0], 0.05, coordinate);
        let (b, cb) = nonlinear_material(r[1], CP[1], ALPHA[1], 0.08, coordinate);
        let (cards, curves) = ([a, b], [ca, cb]);
        let context = query(300.0, coordinate);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = bindings(&cards, &curves, &context, &key);
        let regions = regions();
        let layout = SpatialConductionThermalNetwork {
            cells: &cells,
            links: &LINKS,
            regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        let mut thermal = initial(&curves, 999.0, None);
        let mut field = problem(&cards, coordinate);
        let mut old = [0.0; 2];
        for step in 1..=2 {
            let voltage = boundary_response(&field, 0.0, Some(6.0));
            let before = thermal.energy_j();
            let (j, theta, e) = voltage_reference(6.0, old, r);
            let next = layout
                .step_response_candidate(
                    &thermal,
                    step as f64 * 0.1,
                    &field,
                    &voltage,
                    limits(),
                    field_controls(),
                    controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
                near(
                    next.thermal().temperature_k(cell).unwrap(),
                    300.0 + theta[i],
                );
                for tet in 6 * i..6 * i + 6 {
                    let report = &next.report().electrical.field;
                    near(report.field.current_density_a_m2[tet][0], -j);
                    near(report.field.electric_field_v_m[tet][0], -e[i]);
                    near(report.field.cell_power_w[tet], e[i] * j / 6.0);
                    let source = &report.evaluated[tet].conductor;
                    near(source.material().query().axes()["T"], 300.0 + theta[i]);
                    near(source.material().query().axes()["operating"], e[i]);
                    assert_eq!(source.source_key().name(), "original_nonlinear_transport");
                    assert!(source.has_same_electrical_source(&field.conductors()[tet]));
                }
            }
            let report = next.report();
            near(report.electrical.heat_j, 0.1 * 6.0 * j);
            near(report.electrical.supplied_work_j, report.electrical.heat_j);
            near(next.thermal().energy_j() - before, report.deposited_joule_j);
            near(report.total_energy_residual_j, 0.0);
            assert!(report.electrical.field.field_solves > 1);
            assert!(report.circuit_trials > report.electrical.field.field_solves);
            assert_eq!(next.problem().relaxation(), 0.5);
            if step == 1 {
                heats.push(report.deposited_joule_j);
            }
            (thermal, field, _) = next.into_parts();
            old = theta;
        }
        assert_eq!(context.axes()["operating"], 0.0);
    }
    assert!(heats[1] < heats[0]);
}

#[test]
fn mixed_current_electrodes_close_actual_source_heat_and_preserve_shared_budget() {
    let coordinate = ConductorMagnitude::CurrentDensity;
    let r = [1.0, 2.0];
    let beta = [0.2, 0.08];
    let (a, ca) = nonlinear_material(r[0], CP[0], ALPHA[0], beta[0], coordinate);
    let (b, cb) = nonlinear_material(r[1], CP[1], ALPHA[1], beta[1], coordinate);
    let (cards, curves) = ([a, b], [ca, cb]);
    let context = query(300.0, coordinate);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &key);
    let regions = regions();
    let layout = SpatialConductionThermalNetwork {
        cells: &cells,
        links: &LINKS,
        regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let thermal = initial(&curves, 999.0, None);
    let mut solve = config();
    solve.maximum_circuit_trials = 4096;
    let original = problem(&cards, coordinate);
    // Separate volume meshes joined by an original floating equipotential
    // electrode. Three terminals cost TWO basis solves in every iteration.
    let (mut positions, mut tets) = slab(&[0.0, 1.0]);
    let (second, second_tets) = slab(&[1.0, 2.0]);
    positions.extend(second);
    tets.extend(second_tets.into_iter().map(|tet| tet.map(|v| v + 8)));
    let problem = SpatialConductionResponse::new(
        positions,
        tets,
        original.sources().to_vec(),
        limits(),
        &mut || false,
    )
    .unwrap()
    .with_relaxation(0.5)
    .unwrap();
    let faces = [0.0, 1.0, 2.0].map(|x| {
        problem
            .exterior_faces()
            .filter(|&face| {
                problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0]) == [x; 3]
            })
            .collect::<Vec<_>>()
    });
    let terminals = faces.each_ref().map(Vec::as_slice);
    let drive = SpatialConductionResponseDrive::MixedTerminals {
        terminals: &terminals,
        current_terminals: &[1, 2],
        voltages_v: &[13.0],
        currents_a: &[0.0, 2.0],
        response: SpatialConductanceControls {
            maximum_terminals: 3,
            maximum_field_solves: 2,
            maximum_response_values: 1000,
            maximum_dense_work: 1000,
            reference_voltage_v: 1.0,
            coefficient_tolerance_s: 1e-10,
        },
        mixed: SpatialMixedConductanceControls {
            maximum_current_terminals: 2,
            maximum_response_values: 1000,
            maximum_dense_work: 1000,
        },
    };
    let next = layout
        .step_response_candidate(
            &thermal,
            0.1,
            &problem,
            drive,
            limits(),
            field_controls(),
            controls(),
            solve,
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let (theta, e) = endpoint(2.0, 0.1, [0.0; 2], r, beta, false);
    for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
        near(
            next.thermal().temperature_k(cell).unwrap(),
            300.0 + theta[i],
        );
    }
    // The recovered field's nodal voltages are relative to terminal zero (fs-em
    // spatial/conduction/matrix.rs field_at_voltages: "a common potential is removed
    // before all field and power arithmetic"). Terminal zero is the voltage terminal held
    // at 13 V, so the far electrode reads the two slabs' drop, not 13 V plus it.
    near(
        next.report().electrical.field.field.voltage_v[12],
        e[0] + e[1],
    );
    near(
        next.report().electrical.supplied_work_j,
        0.2 * (e[0] + e[1]),
    );
    near(
        next.thermal().energy_j() - thermal.energy_j(),
        next.report().deposited_joule_j,
    );
    assert_eq!(
        next.report().electrical.field.field_solves,
        2 * next.report().electrical.field.iterations
    );
    let mut bounded = solve;
    bounded.maximum_circuit_trials = next.report().circuit_trials - 1;
    assert!(matches!(
        layout.step_response_candidate(
            &thermal,
            0.1,
            &problem,
            drive,
            limits(),
            field_controls(),
            controls(),
            bounded,
            &[0.0; 2],
            || false,
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let unsupported = match drive {
        SpatialConductionResponseDrive::MixedTerminals {
            terminals,
            current_terminals,
            voltages_v,
            response,
            mixed,
            ..
        } => SpatialConductionResponseDrive::MixedTerminals {
            terminals,
            current_terminals,
            voltages_v,
            currents_a: &[0.0, 25.0],
            response,
            mixed,
        },
        _ => unreachable!("fixture uses mixed electrodes"),
    };
    assert!(matches!(
        layout.step_response_candidate(
            &thermal,
            0.1,
            &problem,
            unsupported,
            limits(),
            field_controls(),
            controls(),
            solve,
            &[0.0; 2],
            || false,
        ),
        Err(CircuitEnthalpyError::Source(
            ConductionThermalError::SpatialResponse(SpatialConductionResponseError::Material(_))
        ))
    ));
    assert!(matches!(
        layout.step_response_candidate(
            &thermal,
            0.1,
            &problem,
            drive,
            limits(),
            field_controls(),
            controls(),
            solve,
            &[0.0; 2],
            || true,
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(
        problem.conductors()[0].material().query().axes()["operating"],
        0.0
    );
}
