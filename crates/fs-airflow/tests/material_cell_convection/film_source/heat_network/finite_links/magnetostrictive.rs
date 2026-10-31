//! Finite source coolant changes actual biased current and reciprocal rod motion.
//! The rod's source-constant tangent remains at its independently held temperature.
use super::*;
use fs_couple::magnetostrictive::{MagnetostrictiveLimits, MagnetostrictiveRod};
use fs_couple::material_circuit::magnetostrictive::thermal::{
    MagnetostrictiveThermalError, ThermalMagnetostrictiveCircuit,
};
use fs_material::{
    field_storage::MAGNETIC_PERMEABILITY_DIMS,
    magnetostrictive::{
        COMPLIANCE_33_H_PROPERTY, COMPLIANCE_DIMS, LongitudinalMagnetostrictive,
        MAGNETIC_FIELD_DIMS, PERMEABILITY_33_T_PROPERTY, PIEZOMAGNETIC_33_PROPERTY,
        PIEZOMAGNETIC_DIMS,
    },
};
use fs_qty::QuantitySpec;

const STEP: f64 = 0.02;
const VOLTAGE: f64 = 0.85;
const ROD_LENGTH: f64 = 0.2;
const ROD_AREA: f64 = 0.005;
const ROD_DENSITY: f64 = 3.0;
const TURNS: u32 = 3;
const COMPLIANCE: f64 = 0.1;
const PERMEABILITY: f64 = 2.0;
const BIAS_FIELD: f64 = 10.0;

fn rod(d: f64) -> MagnetostrictiveRod {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (COMPLIANCE_33_H_PROPERTY, COMPLIANCE_DIMS, COMPLIANCE),
        (PIEZOMAGNETIC_33_PROPERTY, PIEZOMAGNETIC_DIMS, d),
        (
            PERMEABILITY_33_T_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            PERMEABILITY,
        ),
        (DENSITY_PROPERTY, Density::DIMS, ROD_DENSITY),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature_quantity(), 300.0, 300.0)
                    .with_quantity(
                        "bias_H",
                        QuantitySpec::dimensional(MAGNETIC_FIELD_DIMS),
                        BIAS_FIELD,
                        BIAS_FIELD,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored finite-coolant biased-rod fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let material = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic magnetostrictive rod".into(),
            phase: "solid".into(),
            process: "source substitution test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let context = QueryPoint::new()
        .with_quantity("T", temperature_quantity(), 300.0)
        .unwrap()
        .with_quantity(
            "bias_H",
            QuantitySpec::dimensional(MAGNETIC_FIELD_DIMS),
            BIAS_FIELD,
        )
        .unwrap();
    let law = LongitudinalMagnetostrictive::resolve(
        &material,
        &context,
        "bias_H",
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    MagnetostrictiveRod::new(
        law,
        ROD_LENGTH,
        ROD_AREA,
        TURNS,
        MagnetostrictiveLimits {
            maximum_abs_strain: 0.05,
            maximum_abs_field_increment_a_m: 3.0,
        },
    )
    .unwrap()
}

// Independent reciprocal storage coefficients derived directly from the
// authored compliance, permeability, piezomagnetic response and geometry.
fn coefficients(d: f64) -> [f64; 5] {
    let turns = f64::from(TURNS);
    [
        ROD_DENSITY * ROD_AREA * ROD_LENGTH / 3.0,
        ROD_AREA / (COMPLIANCE * ROD_LENGTH),
        turns * turns * ROD_AREA * (PERMEABILITY - d * d / COMPLIANCE) / ROD_LENGTH,
        turns * ROD_AREA * d / (COMPLIANCE * ROD_LENGTH),
        ROD_LENGTH * BIAS_FIELD / turns,
    ]
}

struct Reference {
    state: [f64; 5], // winding T, coolant T, extension, velocity, incremental current
    linkage_wb: f64,
    midpoint_current_a: f64,
    joule_j: f64,
    bias_storage_j: f64,
    stored_change_j: f64,
    supplied_j: f64,
    coolant_power_w: f64,
}

// Independent midpoint circuit/rod equations and one scalar coolant balance.
// No production circuit, material coefficient or convection solver supplies
// expected values. The DC-bias storage term is essential to the actual work.
fn reference(old: [f64; 5], viscosity: f64, cp: f64, d: f64) -> Reference {
    let [wall0, fluid0, x0, v0, i0] = old;
    let [mass, stiffness, inductance, coupling, bias] = coefficients(d);
    let linkage0 = inductance * i0 + coupling * x0;
    let evaluate = |wall: f64| {
        let resistance = 1.0 + 0.02 * (wall - 300.0);
        let a = 2.0 * inductance / STEP + resistance;
        let electrical = VOLTAGE - resistance * bias + 2.0 * inductance * i0 / STEP;
        let vm = (-stiffness * x0 + 2.0 * mass * v0 / STEP + coupling * electrical / a)
            / (2.0 * mass / STEP + stiffness * STEP / 2.0 + coupling * coupling / a);
        let im = (electrical - coupling * vm) / a;
        let (x, v, i) = (x0 + STEP * vm, 2.0 * vm - v0, 2.0 * im - i0);
        let linkage = inductance * i + coupling * x;
        let current = bias + im;
        let joule = STEP * resistance * current * current;
        let bias_storage = bias * (linkage - linkage0);
        let stored_change = 0.5 * mass * (v * v - v0 * v0)
            + 0.5 * stiffness * (x * x - x0 * x0)
            + 0.5 * inductance * (i * i - i0 * i0)
            + bias_storage;
        let supplied = STEP * VOLTAGE * current;
        let fluid = fluid0 + (joule - MASS * CP * (wall - wall0)) / (MASS * cp);
        let film = 0.5 * (wall + fluid);
        let k = 0.5 + 0.002 * (film - 300.0);
        let eta = viscosity * 0.001 * (1.0 - 0.01 * (film - 300.0));
        let h = 0.664 * (1000.0 * 0.1 / eta).sqrt() * (eta * cp / k).cbrt() * k / 0.1;
        let power = SURFACE * h * (wall - fluid);
        (
            MASS * CP * (wall - wall0) - joule + STEP * power,
            Reference {
                state: [wall, fluid, x, v, i],
                linkage_wb: linkage,
                midpoint_current_a: current,
                joule_j: joule,
                bias_storage_j: bias_storage,
                stored_change_j: stored_change,
                supplied_j: supplied,
                coolant_power_w: power,
            },
        )
    };
    let (mut low, mut high) = (280.0, 320.0);
    assert!(evaluate(low).0 < 0.0 && evaluate(high).0 > 0.0);
    for _ in 0..100 {
        let wall = 0.5 * (low + high);
        if evaluate(wall).0 > 0.0 {
            high = wall;
        } else {
            low = wall;
        }
    }
    let (residual, result) = evaluate(0.5 * (low + high));
    assert!(residual.abs() < 1e-10);
    near(
        result.stored_change_j + result.joule_j,
        result.supplied_j,
        1e-12,
    );
    assert!(result.state[2].abs() / ROD_LENGTH < 0.05);
    assert!((f64::from(TURNS) * result.state[4] / ROD_LENGTH).abs() < 3.0);
    result
}

#[test]
fn g1_g3_finite_coolant_changes_biased_current_work_and_reciprocal_motion() {
    let run = |viscosity, cp, d| {
        let (winding_card, winding_curve) = body_material();
        let (fluid_card, fluid_curve) = fluid_material(viscosity, cp);
        let winding_context = body_point(310.0);
        let fluid_context = fluid_point(290.0);
        let thermal = initial_pair(&winding_curve, &fluid_curve, cp, false, 0.0);
        let fluid = liquid_law(&fluid_card);
        let cells = [
            (BODY, 0, &winding_card, &winding_curve, &winding_context),
            (FLUID, 1, &fluid_card, &fluid_curve, &fluid_context),
        ]
        .map(|(cell, id, card, curve, context)| MaterialCell {
            cell,
            material_id: MaterialId(id),
            curve: curve.curve(),
            point: ThermalMaterialPoint {
                material: card,
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        });
        let resistors = [CellResistor {
            destination: ResistorCell {
                branch: 0,
                cell: BODY,
            },
            selection: MaterialPropertySelection::SingleClaimOnly,
            length_m: 0.1,
            area_m2: 1e-5,
        }];
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &[],
        };
        let bases = [(0, ElectricalPropertyBasis::Resistivity)];
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &bases,
        };
        let heat = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let links = [binding()];
        let initial = evaluate_source_link_convection(
            &thermal,
            &links,
            1,
            |_, _, t| fluid.at_temperature(t),
            || false,
        )
        .unwrap();
        let prototype = rod(d);
        let [mass, stiffness, inductance, coupling, bias] = coefficients(d);
        near(prototype.modal_mass_kg(), mass, 1e-14);
        near(prototype.fixed_current_stiffness_n_m(), stiffness, 1e-14);
        near(prototype.clamped_inductance_h(), inductance, 1e-14);
        near(prototype.coupling_wb_m(), coupling, 1e-14);
        near(prototype.bias_current_a(), bias, 1e-14);
        let initial_state = [310.0, 290.0, 0.0001, -0.003, 0.01];
        let linkage0 = inductance * initial_state[4] + coupling * initial_state[2];
        // Finite links require their physical owner even during construction.
        assert!(
            ThermalMagnetostrictiveCircuit::new(
                prototype.clone(),
                linkage0,
                initial_state[2],
                initial_state[3],
                VOLTAGE,
                &thermal,
                &bulk,
                ElectricalPropertyBasis::Resistivity,
                || false
            )
            .is_err()
        );
        let driver = ThermalMagnetostrictiveCircuit::new_with_link_conductances(
            prototype.clone(),
            linkage0,
            initial_state[2],
            initial_state[3],
            VOLTAGE,
            &thermal,
            &bulk,
            ElectricalPropertyBasis::Resistivity,
            initial.conductances_w_k(),
            || false,
        )
        .unwrap();
        let initial_cells = thermal.cells().clone();
        let initial_links = thermal.links().clone();
        let gate = fs_exec::CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let mut result = [0.0; 5];
        pool.scope(|arena| {
            let key = fs_exec::StreamKey {
                seed: 197,
                kernel_id: 41,
                tile: 0,
                iteration: 0,
            };
            let cx = fs_exec::Cx::new(
                &gate,
                arena,
                key,
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            assert!(heat.prepare(&cx, &thermal).is_err());
            let prepared = heat
                .prepare_with_link_conductances(&cx, &thermal, initial.conductances_w_k())
                .unwrap();
            let source = SourceConvectionHeatNetwork {
                heat: &prepared,
                cells: &[],
                maximum_cells: 0,
                resolve_fluid: |_: &SourceConvectiveCell,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    fluid.at_temperature(t)
                },
            };
            let assembly = source.with_links(&links, 1, |_, _, t| fluid.at_temperature(t));
            let mut solve = controls();
            solve.maximum_cells = 2;
            solve.maximum_circuit_trials = 1000;
            solve.newton.linear_restart = 2;
            // Initial conductances cannot become an implicitly frozen coolant law.
            assert!(
                driver
                    .step_candidate(&thermal, STEP, VOLTAGE, 0.0, solve, &[0.0; 2], || false)
                    .is_err()
            );
            let (mut accepted, mut actuator) = (thermal.clone(), driver);
            let mut state = initial_state;
            for tick in 1..=2 {
                let expected = reference(state, viscosity, cp, d);
                let old_cells = accepted.cells().clone();
                let old_circuit = actuator.circuit().clone();
                let old_coil = actuator.coil().clone();
                let old_time = actuator.time_s();
                let next = actuator
                    .step_candidate_with_heat_assembly(
                        &cx,
                        &assembly,
                        &accepted,
                        f64::from(tick) * STEP,
                        VOLTAGE,
                        0.0,
                        solve,
                        &[0.0; 2],
                    )
                    .unwrap();
                let report = next.report();
                let step = &report.electrical;
                near(
                    next.thermal().temperature_k(BODY).unwrap(),
                    expected.state[0],
                    3e-7,
                );
                near(
                    next.thermal().temperature_k(FLUID).unwrap(),
                    expected.state[1],
                    3e-7,
                );
                near(
                    next.driver().circuit().extension_m(),
                    expected.state[2],
                    1e-9,
                );
                near(
                    next.driver().circuit().velocity_m_s(),
                    expected.state[3],
                    3e-8,
                );
                near(
                    next.driver().current_a().unwrap(),
                    bias + expected.state[4],
                    3e-8,
                );
                near(
                    next.driver().circuit().flux_linkage_increment_wb(),
                    expected.linkage_wb,
                    3e-9,
                );
                near(step.midpoint_current_a, expected.midpoint_current_a, 3e-8);
                near(step.coil_heat_j, expected.joule_j, 3e-8);
                near(step.bias_storage_change_j, expected.bias_storage_j, 3e-8);
                near(step.stored_change_j, expected.stored_change_j, 3e-8);
                near(step.electrical_supply_j, expected.supplied_j, 3e-8);
                near(step.mechanical_supply_j, 0.0, 1e-12);
                near(step.energy_residual_j, 0.0, 3e-8);
                near(report.deposited_joule_j, expected.joule_j, 3e-8);
                near(
                    report.thermal.stored_energy_change_j,
                    expected.joule_j,
                    3e-8,
                );
                near(report.non_electrical_heat_j, 0.0, 1e-12);
                near(report.total_energy_residual_j, 0.0, 3e-8);
                near(
                    next.thermal().energy_j() - accepted.energy_j() + step.stored_change_j
                        - step.electrical_supply_j,
                    0.0,
                    3e-8,
                );
                near(
                    MASS * CP * (expected.state[0] - state[0]),
                    expected.joule_j - STEP * expected.coolant_power_w,
                    1e-10,
                );
                near(
                    MASS * cp * (expected.state[1] - state[1]),
                    STEP * expected.coolant_power_w,
                    1e-10,
                );
                assert!(
                    step.coil_heat_j > 100.0 * step.incremental.incremental_resistor_dissipation_j
                );
                assert!(step.bias_storage_change_j.abs() > 1e-5);
                let coolant = &report.source.finite_convection.as_ref().unwrap().samples()[0];
                assert_eq!((coolant.wall_cell, coolant.fluid_cell), (BODY, FLUID));
                near(coolant.heat_into_fluid_w, expected.coolant_power_w, 3e-6);
                near(
                    coolant.property_temperature_k,
                    0.5 * (expected.state[0] + expected.state[1]),
                    3e-7,
                );
                assert_eq!(
                    coolant.material.material().card_identity(),
                    fluid_card.content_hash()
                );
                assert_eq!(coolant.material.material().query().axes()["p"], 100000.0);
                assert_eq!(
                    coolant.material.material().query().axes()["solute_mass_fraction"],
                    0.2
                );
                assert_eq!(
                    step.coil.material().card_identity(),
                    winding_card.content_hash()
                );
                assert_eq!(
                    step.coil.material().query(),
                    body_point(next.thermal().temperature_k(BODY).unwrap())
                );
                assert_eq!(next.driver().circuit().rod(), &prototype);
                assert_eq!(next.thermal().links(), &initial_links);
                assert_eq!(accepted.cells(), &old_cells);
                assert_eq!(actuator.circuit(), &old_circuit);
                assert_eq!(actuator.coil(), &old_coil);
                assert_eq!(actuator.time_s(), old_time);
                assert_eq!(next.driver().time_s(), next.thermal().time_s());
                result = [
                    next.thermal().temperature_k(BODY).unwrap(),
                    next.thermal().temperature_k(FLUID).unwrap(),
                    next.driver().circuit().extension_m(),
                    next.driver().circuit().velocity_m_s(),
                    next.driver().current_a().unwrap(),
                ];
                state = expected.state;
                (accepted, actuator, _) = next.into_parts();
            }
            let saved_cells = accepted.cells().clone();
            let saved_circuit = actuator.circuit().clone();
            let saved_coil = actuator.coil().clone();
            let saved_time = actuator.time_s();
            let queries = std::cell::Cell::new(0);
            let refused = SourceConvectionHeatNetwork {
                heat: &prepared,
                cells: &[],
                maximum_cells: 0,
                resolve_fluid: |_: &SourceConvectiveCell,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    Ok::<_, &'static str>(fluid.at_temperature(t).unwrap())
                },
            };
            let failing = refused.with_links(&links, 1, |_, _, t| {
                queries.set(queries.get() + 1);
                if queries.get() > 1 {
                    Err("actuator coolant trial refused")
                } else {
                    Ok(fluid.at_temperature(t).unwrap())
                }
            });
            assert!(matches!(
                actuator.step_candidate_with_heat_assembly(
                    &cx,
                    &failing,
                    &accepted,
                    3.0 * STEP,
                    VOLTAGE,
                    0.0,
                    solve,
                    &[0.0; 2]
                ),
                Err(CircuitEnthalpyError::Source(
                    MagnetostrictiveThermalError::Heat(ConvectionHeatError::Convection(
                        SourceCellConvectionError::Source {
                            cell: FLUID,
                            source: "actuator coolant trial refused",
                            ..
                        }
                    ))
                ))
            ));
            assert!(queries.get() > 1);
            let mut insufficient = solve;
            insufficient.maximum_circuit_trials = 1;
            assert!(matches!(
                actuator.step_candidate_with_heat_assembly(
                    &cx,
                    &assembly,
                    &accepted,
                    3.0 * STEP,
                    VOLTAGE,
                    0.0,
                    insufficient,
                    &[0.0; 2]
                ),
                Err(CircuitEnthalpyError::CircuitBudget)
            ));
            gate.request();
            assert!(matches!(
                actuator.step_candidate_with_heat_assembly(
                    &cx,
                    &assembly,
                    &accepted,
                    3.0 * STEP,
                    VOLTAGE,
                    0.0,
                    solve,
                    &[0.0; 2]
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            let retry_gate = fs_exec::CancelGate::new();
            let retry_cx = fs_exec::Cx::new(
                &retry_gate,
                arena,
                key,
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            let retry = || {
                actuator
                    .step_candidate_with_heat_assembly(
                        &retry_cx,
                        &assembly,
                        &accepted,
                        3.0 * STEP,
                        VOLTAGE,
                        0.0,
                        solve,
                        &[0.0; 2],
                    )
                    .unwrap()
            };
            let a = retry();
            let b = retry();
            assert_eq!(a.thermal().cells(), b.thermal().cells());
            assert_eq!(a.driver().circuit(), b.driver().circuit());
            assert_eq!(a.driver().coil(), b.driver().coil());
            let expected = reference(state, viscosity, cp, d);
            near(
                a.thermal().temperature_k(BODY).unwrap(),
                expected.state[0],
                3e-7,
            );
            near(
                a.thermal().temperature_k(FLUID).unwrap(),
                expected.state[1],
                3e-7,
            );
            near(a.driver().circuit().extension_m(), expected.state[2], 1e-9);
            near(a.driver().circuit().velocity_m_s(), expected.state[3], 3e-8);
            near(
                a.driver().current_a().unwrap(),
                bias + expected.state[4],
                3e-8,
            );
            near(a.report().total_energy_residual_j, 0.0, 3e-8);
            assert_eq!(accepted.cells(), &saved_cells);
            assert_eq!(actuator.circuit(), &saved_circuit);
            assert_eq!(actuator.coil(), &saved_coil);
            assert_eq!(actuator.time_s(), saved_time);
            assert_eq!(thermal.cells(), &initial_cells);
            assert_eq!(thermal.links(), &initial_links);
        });
        result
    };
    let baseline = run(1.0, 4000.0, 0.2);
    let viscous = run(4.0, 4000.0, 0.2);
    let higher_capacity = run(1.0, 8000.0, 0.2);
    let weaker_coupling = run(1.0, 4000.0, 0.1);
    assert!(viscous[0] > baseline[0] && viscous[4] < baseline[4] && viscous[2] < baseline[2]);
    assert!(
        higher_capacity[0] < baseline[0]
            && higher_capacity[4] > baseline[4]
            && higher_capacity[2] > baseline[2]
    );
    assert!(weaker_coupling[2] < 0.6 * baseline[2] && weaker_coupling[4] < baseline[4]);
}
