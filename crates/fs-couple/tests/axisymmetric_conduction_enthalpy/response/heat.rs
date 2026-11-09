//! G1/G3: a source thermal interface changes native meridional current and heat.
use super::*;
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
    transient::deposition::contact::{ContactCell, MaterialContactLink},
};
use fs_couple::circuit_enthalpy::{
    conduction::axisymmetric::response::{
        AxisymmetricConductionResponseDrive, step_axisymmetric_conduction_response_candidate,
    },
    contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
    material::MaterialNetwork,
};
use fs_exec::{CancelGate, Cx};
use fs_matdb::ClaimSelection;

const DT: f64 = 0.5;
const VOLTAGE: f64 = 8.0;
const ELECTRICAL_RESISTANCE: f64 = 0.5;
const CONTACT_SLOPE: f64 = 0.002;
const FRACTION: f64 = 0.25;

fn joint(cards: &[MaterialCard; 2], thermal_resistance: f64) -> InterfaceSystemCard {
    let original = crate::contacts::interface(cards, ELECTRICAL_RESISTANCE, CONTACT_SLOPE, 1000.0);
    let mut claims = original.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(RP, RD),
            value: PropertyValue::Scalar {
                value: thermal_resistance,
                dims: RD,
            },
            validity: ValidityDomain::unconstrained()
                .with("Ta", 200.0, 1000.0)
                .with("Tb", 200.0, 1000.0)
                .with("pressure", 1e5, 1e5),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured independent thermal resistance on the original joint".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    InterfaceSystemCard::assemble(
        original.surface_a().clone(),
        original.surface_b().clone(),
        original.context().clone(),
        claims,
        vec![],
    )
    .unwrap()
}

// Independent elimination of both BE body balances, divided by swept area.
// rho_i = BASE_i + ALPHA_i*theta_i + BETA_i*|E_i|,
// R_e'' = r0*(1+beta*(theta_a+theta_b)/2), G/A = 1/(1+R_th'').
fn at_current(j: f64, old: [f64; 2], resistance: f64) -> ([f64; 2], [f64; 2], f64) {
    let inv = BETA.map(|b| 1.0 / (1.0 - b * j));
    let g = 1.0 / (1.0 + resistance);
    let p = DT * j * j;
    let slope = ELECTRICAL_RESISTANCE * CONTACT_SLOPE / 2.0;
    let shares = [FRACTION, 1.0 - FRACTION];
    let a = CP[0] + DT * g - p * (ALPHA[0] * inv[0] + shares[0] * slope);
    let b = -DT * g - p * shares[0] * slope;
    let c = -DT * g - p * shares[1] * slope;
    let d = CP[1] + DT * g - p * (ALPHA[1] * inv[1] + shares[1] * slope);
    let rhs: [f64; 2] = std::array::from_fn(|i| {
        CP[i] * old[i] + p * (BASE[i] * inv[i] + shares[i] * ELECTRICAL_RESISTANCE)
    });
    let determinant = a * d - b * c;
    let theta = [
        (rhs[0] * d - b * rhs[1]) / determinant,
        (a * rhs[1] - c * rhs[0]) / determinant,
    ];
    let e = std::array::from_fn(|i| j * (BASE[i] + ALPHA[i] * theta[i]) * inv[i]);
    let contact = ELECTRICAL_RESISTANCE + slope * (theta[0] + theta[1]);
    (theta, e, contact)
}

fn reference(old: [f64; 2], resistance: f64) -> (f64, [f64; 2], [f64; 2], f64) {
    let (mut lo, mut hi) = (0.0, 3.0);
    let voltage = |j| {
        let (_, e, contact) = at_current(j, old, resistance);
        e[0] + e[1] + j * contact
    };
    assert!(voltage(lo) < VOLTAGE && voltage(hi) > VOLTAGE);
    for _ in 0..64 {
        let mid = (lo + hi) / 2.0;
        if voltage(mid) < VOLTAGE {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let j = (lo + hi) / 2.0;
    let (theta, e, contact) = at_current(j, old, resistance);
    (j, theta, e, contact)
}

#[test]
fn sourced_thermal_joint_controls_current_and_each_original_body_balance() {
    let coordinate = ConductorMagnitude::ElectricField;
    let (a, ca) = nonlinear_material(0, coordinate, None);
    let (b, cb) = nonlinear_material(1, coordinate, None);
    let (cards, curves) = ([a, b], [ca, cb]);
    let context = response_point(300.0, coordinate);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let bulk = MaterialNetwork {
        cells: &bindings,
        resistors: &[],
        links: &LINKS,
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let mut final_currents = Vec::new();
    for resistance in [0.25, 2.0] {
        let card = joint(&cards, resistance);
        let surface_context = crate::contacts::interface_context(300.0);
        let electrical_binding = [crate::contacts::binding(
            &cards,
            &card,
            &surface_context,
            FRACTION,
        )];
        let heat_contacts = [MaterialContactLink {
            link: 7,
            point: electrical_binding[0].point,
            sides: [0, 1].map(|i| ContactCell {
                cell: bindings[i].cell,
                material_id: bindings[i].material_id,
                curve: curves[i].curve(),
            }),
            selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        }];
        let assembly = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &heat_contacts,
            radiation: None,
        };
        let area = 0.7 / 2.0;
        let mut thermal = initial(&curves, 0.7);
        let source = resolve_electrical_interface_state(
            &card,
            &surface_context,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let mut field = AxisymmetricConductionResponse::new_with_phase_contacts(
            contact_mesh(0.7),
            bulk_sources(&cards, coordinate),
            vec![ConductionContact {
                first_vertices: [2, 3],
                second_vertices: [5, 4],
                source: source.clone(),
            }],
            vec![[None, None]],
            limits(),
            &mut || false,
        )
        .unwrap();
        let voltage = [(0, 0.0), (1, 0.0), (6, VOLTAGE), (7, VOLTAGE)];
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 47,
                    kernel_id: 91,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            let prepared = assembly.prepare(&cx, &thermal).unwrap();
            let mut solve = config();
            solve.maximum_circuit_trials = 8192;
            let mut limited = solve;
            limited.maximum_circuit_trials = 1;
            assert!(matches!(
                step_axisymmetric_conduction_response_candidate(
                    &prepared,
                    &cx,
                    &thermal,
                    DT,
                    &field,
                    ConductionDrive::Voltage(&voltage),
                    &electrical_binding,
                    &REGIONS,
                    1e-12,
                    limits(),
                    field_controls(),
                    response_controls(),
                    limited,
                    &[0.0; 2],
                ),
                Err(CircuitEnthalpyError::CircuitBudget)
            ));
            assert_eq!(thermal.time_s(), 0.0);
            let mut old = [0.0; 2];
            for interval in 1..=2 {
                let (j, theta, e, contact_resistance) = reference(old, resistance);
                let electrodes: [&[usize]; 2] = [&[0, 1], &[6, 7]];
                let terminal_voltage = [0.0, VOLTAGE];
                let drive = if interval == 1 {
                    AxisymmetricConductionResponseDrive::from(ConductionDrive::Voltage(&voltage))
                } else {
                    AxisymmetricConductionResponseDrive::from(terminals::terminal_drive(
                        &electrodes,
                        fs_em::conduction::terminals::ConductionTerminalControl::Voltages(
                            &terminal_voltage,
                        ),
                    ))
                };
                let next = step_axisymmetric_conduction_response_candidate(
                    &prepared,
                    &cx,
                    &thermal,
                    f64::from(interval) * DT,
                    &field,
                    drive,
                    &electrical_binding,
                    &REGIONS,
                    1e-12,
                    limits(),
                    field_controls(),
                    response_controls(),
                    solve,
                    &[0.0; 2],
                )
                .unwrap();
                let report = next.report();
                let step = &report.electrical;
                let sample = &report.source.thermal_contacts[0];
                assert_eq!(sample.cells, [CellId(30), CellId(10)]);
                near(sample.interface.value_m2_k_per_w(), resistance);
                near(sample.bulk_conductance_w_k, area);
                near(sample.conductance_w_k, area / (1.0 + resistance));
                let actual = [CellId(30), CellId(10)]
                    .map(|cell| next.thermal().temperature_k(cell).unwrap() - 300.0);
                let contact_power = area * j * j * contact_resistance;
                let exchange = area * (theta[0] - theta[1]) / (1.0 + resistance);
                for i in 0..2 {
                    near(actual[i], theta[i]);
                    let share = if i == 0 { FRACTION } else { 1.0 - FRACTION };
                    let outgoing = if i == 0 { exchange } else { -exchange };
                    near(
                        CP[i] * area * (actual[i] - old[i]),
                        DT * (area * j * e[i] + share * contact_power - outgoing),
                    );
                }
                near(step.contacts[0].response.current_a, -area * j);
                near(step.contacts[0].response.power_w, contact_power);
                near(step.heat_j, DT * area * j * VOLTAGE);
                near(step.supplied_work_j, step.heat_j);
                near(report.deposited_joule_j, step.heat_j);
                near(next.thermal().energy_j() - thermal.energy_j(), step.heat_j);
                near(report.non_electrical_heat_j, 0.0);
                let resolved = &next.problem().contacts()[0].source;
                near(resolved.resistance_ohm_m2(), contact_resistance);
                near(
                    resolved.resolved().state_coordinate("Ta").unwrap(),
                    300.0 + theta[0],
                );
                near(
                    resolved.resolved().state_coordinate("Tb").unwrap(),
                    300.0 + theta[1],
                );
                assert_eq!(
                    resolved.resolved().properties()[0]
                        .answer()
                        .receipt
                        .selected,
                    source.resolved().properties()[0].answer().receipt.selected
                );
                if interval == 2 {
                    let terminals = step.terminals.as_ref().unwrap();
                    near(terminals.current_a[1], area * j);
                    near(terminals.supplied_power_w, area * j * VOLTAGE);
                    final_currents.push(terminals.current_a[1]);
                }
                (thermal, field, _) = next.into_parts();
                old = theta;
            }
        });
    }
    assert!((final_currents[0] - final_currents[1]).abs() > 1e-7);
}
