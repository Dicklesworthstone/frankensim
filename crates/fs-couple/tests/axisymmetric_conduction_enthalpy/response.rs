//! G1/G3/G4: actual meridional E/J laws, thermal feedback and swept work.
//! Manufactured material laws; these fixtures make no empirical-data claim.

#[path = "response/heat.rs"]
mod heat;

#[path = "response/terminals.rs"]
mod terminals;

use super::*;
use fs_couple::circuit_enthalpy::conduction::ConductionThermalError;
use fs_em::conduction::{
    ConductionContact,
    axisymmetric::response::AxisymmetricConductionResponse,
    response::{ConductionResponseControls, contact::ContactResponseControls},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyValue, Provenance, SelectionPolicy, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::{
    conductor::{
        ResolvedConductorState,
        response::{ConductorMagnitude, ScalarConductorResponse},
    },
    electrical_interface::{
        ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
        resolve_electrical_interface_state,
        response::{ElectricalInterfaceCoordinate, ScalarElectricalInterfaceResponse},
    },
    phase::{HeatCapacitySourceQuery, SolidLiquidPhase},
    state_point::{
        SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};
use fs_qty::{
    Density, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const CP: [f64; 2] = [2.0, 3.0];
const ALPHA: [f64; 2] = [0.002, 0.004];
const BETA: [f64; 2] = [0.05, 0.08];
const BASE: [f64; 2] = [1.0, 2.0];

fn response_point(t: f64, coordinate: ConductorMagnitude) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
        .with_quantity("operating", coordinate.quantity(), 0.0)
        .unwrap()
}
fn response_controls() -> ConductionResponseControls {
    ConductionResponseControls {
        maximum_iterations: 100,
        absolute_conductivity_tolerance_s_m: 1e-12,
        relative_conductivity_tolerance: 1e-12,
        current_density_tolerance_a_m2: 1e-10,
    }
}
fn nonlinear_material(
    i: usize,
    coordinate: ConductorMagnitude,
    field_bounded: Option<&str>,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let basis = ElectricalPropertyBasis::Resistivity;
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let electrical = PropertyKey::new("original_transport", basis.dimensions());
    let mut claims = ClaimSet::new();
    for (key, value) in [
        (
            heat.clone(),
            PropertyValue::Scalar {
                value: CP[i],
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            },
        ),
        (
            density.clone(),
            PropertyValue::Scalar {
                value: 1.0,
                dims: Density::DIMS,
            },
        ),
        (
            PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            PropertyValue::Scalar {
                value: 1.0,
                dims: THERMAL_CONDUCTIVITY_DIMS,
            },
        ),
        (
            electrical.clone(),
            PropertyValue::Surface2D {
                axes: ["T".into(), "operating".into()],
                axis_dims: [Temperature::DIMS, coordinate.quantity().dims()],
                samples: [[200.0, 0.0], [2000.0, 0.0], [200.0, 20.0], [2000.0, 20.0]]
                    .into_iter()
                    .map(|p| (p, BASE[i] + ALPHA[i] * (p[0] - 300.0) + BETA[i] * p[1]))
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: basis.dimensions(),
            },
        ),
    ] {
        let mut validity = ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            200.0,
            2000.0,
        );
        if key == electrical || field_bounded.is_some_and(|name| name == key.name()) {
            validity = validity.with_quantity("operating", coordinate.quantity(), 0.0, 20.0);
        }
        let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        };
        claims
            .insert_claim(PropertyClaim {
                key,
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured meridional electrothermal law".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
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
            lower: &response_point(200.0, coordinate),
            upper: &response_point(2000.0, coordinate),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}
fn bulk_sources(
    cards: &[MaterialCard; 2],
    coordinate: ConductorMagnitude,
) -> Vec<ScalarConductorResponse> {
    let basis = ElectricalPropertyBasis::Resistivity;
    let key = PropertyKey::new("original_transport", basis.dimensions());
    [0, 0, 1, 1]
        .into_iter()
        .map(|i| {
            let material = resolve_material_state_point(
                &cards[i],
                &response_point(300.0, coordinate),
                &[
                    basis.requirement_for(&key).unwrap(),
                    ScalarPropertyRequirement::try_new(
                        "density",
                        Density::DIMS,
                        ScalarAdmissibility::StrictlyPositive,
                    )
                    .unwrap(),
                    ScalarPropertyRequirement::try_new(
                        "specific_heat_capacity",
                        SPECIFIC_HEAT_CAPACITY_DIMS,
                        ScalarAdmissibility::StrictlyPositive,
                    )
                    .unwrap(),
                ],
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            let conductor =
                ResolvedConductorState::try_from_material_with_property(material, basis, &key)
                    .unwrap();
            ScalarConductorResponse::new(&cards[i], conductor, &[("operating", coordinate)])
                .unwrap()
        })
        .collect()
}
fn nonlinear(
    cards: &[MaterialCard; 2],
    coordinate: ConductorMagnitude,
    angle: f64,
) -> AxisymmetricConductionResponse {
    AxisymmetricConductionResponse::new(
        mesh(angle),
        bulk_sources(cards, coordinate),
        limits(),
        &mut || false,
    )
    .unwrap()
}

// Divide the two independent BE body balances by A=angle/2. Each body has
// volume A, mass A and G=A. Under fixed |J|, eliminate rho(E,T) algebraically.
fn endpoint(j: f64, dt: f64, old: [f64; 2], field_law: bool) -> ([f64; 2], [f64; 2]) {
    let inverse = BETA.map(|b| if field_law { 1.0 / (1.0 - b * j) } else { 1.0 });
    let base =
        std::array::from_fn::<_, 2, _>(|i| BASE[i] + if field_law { 0.0 } else { BETA[i] * j });
    let a = std::array::from_fn::<_, 2, _>(|i| CP[i] + dt - dt * j * j * ALPHA[i] * inverse[i]);
    let rhs =
        std::array::from_fn::<_, 2, _>(|i| CP[i] * old[i] + dt * j * j * base[i] * inverse[i]);
    let det = a[0] * a[1] - dt * dt;
    let theta = [
        (rhs[0] * a[1] + dt * rhs[1]) / det,
        (rhs[1] * a[0] + dt * rhs[0]) / det,
    ];
    let e = std::array::from_fn(|i| j * (base[i] + ALPHA[i] * theta[i]) * inverse[i]);
    (theta, e)
}

#[test]
fn nonlinear_sources_and_all_drives_close_original_sector_heat_and_work() {
    let j = 2.0;
    let dt = 0.1;
    for coordinate in [
        ConductorMagnitude::ElectricField,
        ConductorMagnitude::CurrentDensity,
    ] {
        let field_law = coordinate == ConductorMagnitude::ElectricField;
        let (a, ca) = nonlinear_material(0, coordinate, None);
        let (b, cb) = nonlinear_material(1, coordinate, None);
        let (cards, curves) = ([a, b], [ca, cb]);
        let context = response_point(300.0, coordinate);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let bindings = cells(&cards, &curves, &context, &key);
        let layout = network(&bindings);
        let (theta, e) = endpoint(j, dt, [0.0; 2], field_law);
        let voltage = e[0] + e[1];
        for angle in [0.7, 2.0, core::f64::consts::TAU] {
            let area = angle / 2.0;
            let thermal = initial(&curves, angle);
            let problem = nonlinear(&cards, coordinate, angle);
            let voltages = [(0, 0.0), (1, 0.0), (4, voltage), (5, voltage)];
            let current_traces = [ConductionCurrentTrace {
                vertices: [4, 5],
                inward_current_density_a_m2: j,
            }];
            for drive in [
                ConductionDrive::Current {
                    positive: &[4, 5],
                    negative: &[0, 1],
                    current_a: j * area,
                },
                ConductionDrive::CurrentDensity {
                    voltage: &[(0, 0.0), (1, 0.0)],
                    currents: &current_traces,
                },
                ConductionDrive::Voltage(&voltages),
                ConductionDrive::Power {
                    positive: &[4, 5],
                    negative: &[0, 1],
                    source: ConductionPowerSource {
                        power_w: area * j * voltage,
                        maximum_voltage_v: 10.0,
                        maximum_current_a: 4.0 * area,
                    },
                },
            ] {
                let mut solve = config();
                solve.maximum_circuit_trials = 4096;
                let next = layout
                    .step_axisymmetric_response_candidate(
                        &thermal,
                        dt,
                        &problem,
                        drive,
                        limits(),
                        field_controls(),
                        response_controls(),
                        solve,
                        &[0.0; 2],
                        || false,
                    )
                    .unwrap();
                let step = &next.report().electrical;
                for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
                    near(
                        next.thermal().temperature_k(cell).unwrap(),
                        300.0 + theta[i],
                    );
                    for triangle in 2 * i..2 * i + 2 {
                        near(step.field.electric_field_v_m[triangle][0], 0.0);
                        near(step.field.electric_field_v_m[triangle][1], -e[i]);
                        near(step.field.current_density_a_m2[triangle][1], -j);
                        near(
                            step.field.cell_power_w[triangle],
                            j * e[i] * problem.mesh().face_measure(triangle).unwrap(),
                        );
                        let selected = &step.field.evaluated[triangle].conductor;
                        near(selected.material().query().axes()["T"], 300.0 + theta[i]);
                        near(
                            selected.material().query().axes()["operating"],
                            if field_law { e[i] } else { j },
                        );
                        assert_eq!(selected.source_key().name(), "original_transport");
                        assert_eq!(
                            selected.source_property().answer().receipt.selected,
                            problem.conductors()[triangle]
                                .source_property()
                                .answer()
                                .receipt
                                .selected
                        );
                    }
                }
                near(step.heat_j, dt * area * j * voltage);
                near(step.supplied_work_j, step.heat_j);
                near(next.report().deposited_joule_j, step.heat_j);
                near(
                    next.thermal().energy_j() - thermal.energy_j(),
                    area * (CP[0] * theta[0] + CP[1] * theta[1]),
                );
                near(next.report().source.links[0].conductance_w_k, area);
                assert!(step.field.field_solves > 1);
                assert!(next.report().circuit_trials > step.field.field_solves);
                assert!(
                    next.report().total_energy_residual_j.abs() < solve.total_energy_tolerance_j
                );
            }
            let current = ConductionDrive::Current {
                positive: &[4, 5],
                negative: &[0, 1],
                current_a: j * area,
            };
            let first = layout
                .step_axisymmetric_response_candidate(
                    &thermal,
                    dt,
                    &problem,
                    current,
                    limits(),
                    field_controls(),
                    response_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            let second = layout
                .step_axisymmetric_response_candidate(
                    first.thermal(),
                    2.0 * dt,
                    first.problem(),
                    current,
                    limits(),
                    field_controls(),
                    response_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            let (theta2, e2) = endpoint(j, dt, theta, field_law);
            near(
                second.thermal().temperature_k(CellId(30)).unwrap(),
                300.0 + theta2[0],
            );
            near(
                second.thermal().temperature_k(CellId(10)).unwrap(),
                300.0 + theta2[1],
            );
            near(second.report().electrical.field.voltage_v[4], e2[0] + e2[1]);
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(
                problem.conductors()[0].material().query().axes()["operating"],
                0.0
            );
        }
    }
}

#[test]
fn nonlinear_field_budgets_source_refusals_and_late_cancel_leave_histories_reusable() {
    let coordinate = ConductorMagnitude::ElectricField;
    let (a, ca) = nonlinear_material(0, coordinate, None);
    let (b, cb) = nonlinear_material(1, coordinate, None);
    let (cards, curves) = ([a, b], [ca, cb]);
    let context = response_point(300.0, coordinate);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let layout = network(&bindings);
    let thermal = initial(&curves, 2.0);
    let original = nonlinear(&cards, coordinate, 2.0);
    let drive = ConductionDrive::Current {
        positive: &[4, 5],
        negative: &[0, 1],
        current_a: 2.0,
    };
    let mut count = 0;
    let good = layout
        .step_axisymmetric_response_candidate(
            &thermal,
            0.1,
            &original,
            drive,
            limits(),
            field_controls(),
            response_controls(),
            config(),
            &[0.0; 2],
            || {
                count += 1;
                false
            },
        )
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        layout.step_axisymmetric_response_candidate(
            &thermal,
            0.1,
            &original,
            drive,
            limits(),
            field_controls(),
            response_controls(),
            config(),
            &[0.0; 2],
            || {
                calls += 1;
                calls >= count - 1
            },
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let mut exhausted = config();
    exhausted.maximum_circuit_trials = 1;
    assert!(matches!(
        layout.step_axisymmetric_response_candidate(
            &thermal,
            0.1,
            &original,
            drive,
            limits(),
            field_controls(),
            response_controls(),
            exhausted,
            &[0.0; 2],
            || false,
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    for property in ["specific_heat_capacity", "thermal_conductivity"] {
        let (a, ca) = nonlinear_material(0, coordinate, Some(property));
        let (b, cb) = nonlinear_material(1, coordinate, Some(property));
        let (cards, curves) = ([a, b], [ca, cb]);
        let thermal = initial(&curves, 2.0);
        let problem = nonlinear(&cards, coordinate, 2.0);
        let bindings = cells(&cards, &curves, &context, &key);
        assert!(matches!(
            network(&bindings).step_axisymmetric_response_candidate(
                &thermal,
                0.1,
                &problem,
                drive,
                limits(),
                field_controls(),
                response_controls(),
                config(),
                &[0.0; 2],
                || false,
            ),
            Err(CircuitEnthalpyError::Source(
                ConductionThermalError::Invalid(_)
            ))
        ));
    }
    let retry = layout
        .step_axisymmetric_response_candidate(
            &thermal,
            0.1,
            &original,
            drive,
            limits(),
            field_controls(),
            response_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    assert_eq!(
        retry.report().electrical.field.voltage_v,
        good.report().electrical.field.voltage_v
    );
    assert_eq!(
        retry.thermal().energy_j().to_bits(),
        good.thermal().energy_j().to_bits()
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(
        original.conductors()[0].material().query().axes()["T"],
        300.0
    );
}

fn contact_mesh(angle: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "two cylinder layers with an original nonlinear disk interface",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ],
        vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap()
}

#[test]
fn nonlinear_contact_jump_and_temperature_follow_the_same_two_body_enthalpy() {
    let coordinate = ConductorMagnitude::ElectricField;
    let jump = ElectricalInterfaceCoordinate::VoltageJumpMagnitude;
    let (a, ca) = nonlinear_material(0, coordinate, None);
    let (b, cb) = nonlinear_material(1, coordinate, None);
    let (cards, curves) = ([a, b], [ca, cb]);
    let context = response_point(300.0, coordinate);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let layout = network(&cells);
    let (r0, gamma, eta, fraction) = (0.5, 0.002, 0.05, 0.25);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
                ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Surface2D {
                axes: ["Tm".into(), "jump".into()],
                axis_dims: [Temperature::DIMS, jump.quantity().dims()],
                samples: [[200.0, 0.0], [1000.0, 0.0], [200.0, 10.0], [1000.0, 10.0]]
                    .into_iter()
                    .map(|p| (p, r0 + gamma * (p[0] - 300.0) + eta * p[1]))
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "Tm",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    200.0,
                    1000.0,
                )
                .with_quantity("jump", jump.quantity(), 0.0, 10.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured nonlinear surface temperature/jump law".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "A".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "B".into(),
        },
        SystemContext {
            medium: "manufactured".into(),
            third_body: None,
            environment: "fixed pressure".into(),
            history: "fixed surface".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let surface_context = super::contacts::interface_context(300.0)
        .with_quantity("jump", jump.quantity(), 0.0)
        .unwrap();
    let source = resolve_electrical_interface_state(
        &card,
        &surface_context,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let law =
        ScalarElectricalInterfaceResponse::new(&card, source.clone(), &[("jump", jump)]).unwrap();
    let contact_controls = ContactResponseControls {
        absolute_resistance_tolerance_ohm_m2: 1e-12,
        relative_resistance_tolerance: 1e-12,
        voltage_tolerance_v: 1e-9,
        current_density_tolerance_a_m2: 1e-9,
    };
    let (j, dt) = (2.0, 0.1);
    // Eliminate both bulk E laws and the actual contact jump:
    // R'' = (r0 + gamma*(theta_a+theta_b)/2) / (1-eta*J).
    let inverse = BETA.map(|b| 1.0 / (1.0 - b * j));
    let rc_base = r0 / (1.0 - eta * j);
    let slope = gamma / (2.0 * (1.0 - eta * j));
    let p = dt * j * j;
    let a = CP[0] + dt - p * (ALPHA[0] * inverse[0] + fraction * slope);
    let b = -dt - p * fraction * slope;
    let c = -dt - p * (1.0 - fraction) * slope;
    let d = CP[1] + dt - p * (ALPHA[1] * inverse[1] + (1.0 - fraction) * slope);
    let x = p * (BASE[0] * inverse[0] + fraction * rc_base);
    let y = p * (BASE[1] * inverse[1] + (1.0 - fraction) * rc_base);
    let det = a * d - b * c;
    let theta = [(x * d - b * y) / det, (a * y - c * x) / det];
    let e: [f64; 2] = std::array::from_fn(|i| j * (BASE[i] + ALPHA[i] * theta[i]) * inverse[i]);
    let rc = rc_base + slope * (theta[0] + theta[1]);
    let voltage = e[0] + e[1] + j * rc;
    let bindings = [super::contacts::binding(
        &cards,
        &card,
        &surface_context,
        fraction,
    )];
    for angle in [0.7, core::f64::consts::TAU] {
        let area = angle / 2.0;
        let thermal = initial(&curves, angle);
        let original = AxisymmetricConductionResponse::new_with_contact_response(
            contact_mesh(angle),
            bulk_sources(&cards, coordinate),
            vec![ConductionContact {
                first_vertices: [2, 3],
                second_vertices: [5, 4],
                source: source.clone(),
            }],
            vec![law.clone()],
            contact_controls,
            limits(),
            &mut || false,
        )
        .unwrap();
        let voltages = [(0, 0.0), (1, 0.0), (6, voltage), (7, voltage)];
        for drive in [
            ConductionDrive::Current {
                positive: &[6, 7],
                negative: &[0, 1],
                current_a: j * area,
            },
            ConductionDrive::Voltage(&voltages),
        ] {
            let next = layout
                .step_axisymmetric_response_contact_candidate(
                    &thermal,
                    dt,
                    &original,
                    drive,
                    &bindings,
                    limits(),
                    field_controls(),
                    response_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            let step = &next.report().electrical;
            near(
                next.thermal().temperature_k(CellId(30)).unwrap(),
                300.0 + theta[0],
            );
            near(
                next.thermal().temperature_k(CellId(10)).unwrap(),
                300.0 + theta[1],
            );
            near(step.field.voltage_v[6], voltage);
            near(step.contacts[0].response.current_a, -j * area);
            near(step.contacts[0].response.power_w, area * j * j * rc);
            near(
                step.field.cell_power_w.iter().sum(),
                area * j * (e[0] + e[1]),
            );
            near(step.heat_j, dt * area * j * voltage);
            near(step.supplied_work_j, step.heat_j);
            near(next.thermal().energy_j() - thermal.energy_j(), step.heat_j);
            near(next.report().deposited_joule_j, step.heat_j);
            assert_eq!(step.contacts[0].cells, [CellId(30), CellId(10)]);
            near(step.contacts[0].heat_fraction_to_first, fraction);
            let accepted = &next.problem().contacts()[0].source;
            near(accepted.resistance_ohm_m2(), rc);
            near(
                accepted.resolved().state_coordinate("jump").unwrap(),
                j * rc,
            );
            near(
                accepted.resolved().state_coordinate("Ta").unwrap(),
                300.0 + theta[0],
            );
            near(
                accepted.resolved().state_coordinate("Tb").unwrap(),
                300.0 + theta[1],
            );
            near(
                accepted.resolved().state_coordinate("pressure").unwrap(),
                1e5,
            );
            assert_eq!(
                accepted.resolved().properties()[0]
                    .answer()
                    .receipt
                    .selected,
                source.resolved().properties()[0].answer().receipt.selected
            );
            assert_eq!(step.field.evaluated_contacts.len(), 1);
            // The same nonlinear interface closes under actual mixed terminal
            // control and a large common reference, with its original heat split.
            let electrodes: [&[usize]; 2] = [&[0, 1], &[6, 7]];
            let current = [j * area];
            let terminal = layout
                .step_axisymmetric_response_terminal_candidate(
                    &thermal,
                    dt,
                    &original,
                    terminals::terminal_drive(
                        &electrodes,
                        fs_em::conduction::terminals::ConductionTerminalControl::Mixed {
                            current_terminals: &[1],
                            currents_a: &current,
                            voltages_v: &[1e250],
                        },
                    ),
                    &bindings,
                    limits(),
                    field_controls(),
                    response_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            for cell in [CellId(30), CellId(10)] {
                near(
                    terminal.thermal().temperature_k(cell).unwrap(),
                    next.thermal().temperature_k(cell).unwrap(),
                );
            }
            let terminal_step = &terminal.report().electrical;
            let operating = terminal_step.terminals.as_ref().unwrap();
            near(operating.voltage_v[1], voltage);
            near(operating.current_a[0], -j * area);
            near(operating.current_a[1], j * area);
            near(terminal_step.supplied_work_j, step.supplied_work_j);
            near(terminal_step.heat_j, step.heat_j);
            near(
                terminal_step.contacts[0].response.power_w,
                step.contacts[0].response.power_w,
            );
            near(terminal_step.contacts[0].heat_fraction_to_first, fraction);
            assert_eq!(terminal_step.contacts[0].cells, step.contacts[0].cells);
            near(
                terminal.problem().contacts()[0].source.resistance_ohm_m2(),
                rc,
            );
            near(
                terminal.problem().contacts()[0]
                    .source
                    .resolved()
                    .state_coordinate("jump")
                    .unwrap(),
                j * rc,
            );
            let continued = layout
                .step_axisymmetric_response_contact_candidate(
                    next.thermal(),
                    2.0 * dt,
                    next.problem(),
                    drive,
                    &bindings,
                    limits(),
                    field_controls(),
                    response_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            assert!(continued.problem().contacts()[0].source.resistance_ohm_m2() > rc);
        }
        assert!(
            layout
                .step_axisymmetric_response_candidate(
                    &thermal,
                    dt,
                    &original,
                    ConductionDrive::Current {
                        positive: &[6, 7],
                        negative: &[0, 1],
                        current_a: j * area
                    },
                    limits(),
                    field_controls(),
                    response_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .is_err()
        );
        near(original.contacts()[0].source.resistance_ohm_m2(), r0);
        assert_eq!(thermal.time_s(), 0.0);
    }
}
