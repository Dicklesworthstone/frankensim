//! Finite material coolant changes native winding current and magnetic memory.
//! Source-constant play coefficients retain their admitted thermal interval.
use super::*;
use fs_couple::material_circuit::hysteretic::{
    HystereticCircuit, HystereticCircuitInitial, HystereticCircuitSolve,
    thermal::{HystereticThermalError, ThermalHystereticCircuit},
};
use fs_material::field_storage::{
    MAGNETIC_PERMEABILITY_DIMS,
    hysteresis::{
        UniformHystereticCore,
        source::{HysteresisBranchSources, HysteresisSources, ResolvedScalarHysteresis},
    },
    permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
};
use std::sync::Arc;

const CORE: CellId = CellId(70);
const VOLTAGE: f64 = 5.0;

fn roles() -> HysteresisSources {
    HysteresisSources {
        background_permeability: PropertyKey::new("play_background", MAGNETIC_PERMEABILITY_DIMS),
        field_limit: PropertyKey::new("play_field_limit", FIELD_STRENGTH_DIMS),
        branches: vec![HysteresisBranchSources {
            susceptibility: PropertyKey::new("play_slope", MAGNETIC_PERMEABILITY_DIMS),
            pinning_field: PropertyKey::new("play_pinning", FIELD_STRENGTH_DIMS),
            saturation_polarization: PropertyKey::new("play_saturation", FLUX_DENSITY_DIMS),
        }],
    }
}

fn core_material(pinning: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = body_material();
    let mut claims = ClaimSet::new();
    for (_, source) in original.claims().claims_ordered() {
        let mut claim = source.clone();
        if [SPECIFIC_HEAT_CAPACITY_PROPERTY, DENSITY_PROPERTY].contains(&claim.key.name()) {
            claim.value = PropertyValue::Scalar {
                value: 1.0,
                dims: claim.key.dims(),
            };
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        claims.insert_claim(claim).unwrap();
    }
    let keys = roles();
    for (key, value) in [
        (keys.background_permeability, 0.01),
        (keys.field_limit, 30.0),
        (keys.branches[0].susceptibility.clone(), 0.2),
        (keys.branches[0].pinning_field.clone(), pinning),
        (keys.branches[0].saturation_polarization.clone(), 1.0),
    ] {
        let dims = key.dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature_quantity(), 280.0, 320.0)
                    .with("lot", 7.0, 7.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored finite-coolant play-core fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &body_point(280.0),
            upper: &body_point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}

fn initial<'a>(
    core: &'a ResolvedHeatCapacityCurve,
    winding: &'a ResolvedHeatCapacityCurve,
    fluid: &'a ResolvedHeatCapacityCurve,
    fluid_cp: f64,
) -> DepositionSimulation<'a> {
    let mut events = [
        (CORE, 0, 1.0, 1.0, 20.0),
        (BODY, 1, MASS, 1e-6, CP * 30.0),
        (FLUID, 2, MASS, 1e-6, fluid_cp * 10.0),
    ]
    .map(
        |(id, material, mass_kg, volume_m3, specific_h_j_kg)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: material,
                body_id: material,
                material_id: MaterialId(material),
                mass_kg,
                specific_h_j_kg,
                geometry: CellGeometry {
                    start_m: [0.0, material as f64, 0.0],
                    end_m: [1.0, material as f64, 0.0],
                    volume_m3,
                    surface_area_m2: SURFACE,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        },
    )
    .to_vec();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: LINK,
            a: BODY,
            b: FLUID,
            area_m2: SURFACE,
            conductance_w_k: 0.0,
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), core.curve() as &dyn ThermalLaw),
            (MaterialId(1), winding.curve() as &dyn ThermalLaw),
            (MaterialId(2), fluid.curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 3,
            max_links: 1,
            max_events: 4,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn magnetic_controls() -> HystereticCircuitSolve {
    HystereticCircuitSolve {
        maximum_circuit_trials: 256,
        maximum_newton_iterations: 32,
        difference_step_wb: 1e-6,
        terminal_voltage_tolerance_v: 1e-8,
        energy_tolerance_j: 1e-8,
        circuit: controls().circuit,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-13,
            linear_restart: 2,
            forcing_minimum: 1e-4,
            ..NewtonKrylovConfig::default()
        },
    }
}

// Independent play integral: B=mu*H+P, P=chi*(H-hc) while sliding.
// The virgin crossing gives a quadratic; later monotone segments have mean
// current (H0+H1)/2. One scalar winding balance closes the finite coolant.
fn reference(old: [f64; 7], viscosity: f64, cp: f64, pinning: f64) -> [f64; 12] {
    let [core0, winding0, fluid0, field0, flux0, polarization0, _] = old;
    let (mu, chi) = (0.01, 0.2);
    let sum = mu + chi;
    let evaluate = |winding: f64| {
        let resistance = 1.0 + 0.02 * (winding - 300.0);
        let (field, flux, polarization, stored, mean) = if flux0 == 0.0 {
            let a = 1.0 + DT * resistance / (2.0 * sum);
            let b = DT * resistance * chi * pinning / sum - DT * VOLTAGE;
            let c = -DT * resistance * mu * chi * pinning * pinning / (2.0 * sum);
            let flux = (-b + (b * b - 4.0 * a * c).sqrt()) / (2.0 * a);
            let field = (flux + chi * pinning) / sum;
            let polarization = chi * (field - pinning);
            let stored = 0.5 * mu * field * field + 0.5 * polarization * polarization / chi;
            (
                field,
                flux,
                polarization,
                stored,
                (stored + pinning * polarization) / flux,
            )
        } else {
            let field = ((sum - DT * resistance / 2.0) * field0 + DT * VOLTAGE)
                / (sum + DT * resistance / 2.0);
            let polarization = chi * (field - pinning);
            let stored = 0.5 * mu * field * field + 0.5 * polarization * polarization / chi;
            (
                field,
                sum * field - chi * pinning,
                polarization,
                stored,
                0.5 * (field + field0),
            )
        };
        assert!(field > pinning && polarization < 1.0);
        let core_heat = pinning * (polarization - polarization0);
        let joule = DT * resistance * mean * mean;
        let core = core0 + core_heat;
        let fluid = fluid0 + (joule - MASS * CP * (winding - winding0)) / (MASS * cp);
        let film = 0.5 * (winding + fluid);
        let k = 0.5 + 0.002 * (film - 300.0);
        let eta = viscosity * 0.001 * (1.0 - 0.01 * (film - 300.0));
        let h = 0.664 * (1000.0 * 0.1 / eta).sqrt() * (eta * cp / k).cbrt() * k / 0.1;
        let power = SURFACE * h * (winding - fluid);
        (
            MASS * CP * (winding - winding0) - joule + DT * power,
            [
                core,
                winding,
                fluid,
                field,
                flux,
                polarization,
                stored,
                mean,
                core_heat,
                joule,
                power,
                DT * VOLTAGE * mean,
            ],
        )
    };
    let (mut low, mut high) = (280.0, 320.0);
    assert!(evaluate(low).0 < 0.0 && evaluate(high).0 > 0.0);
    for _ in 0..100 {
        let t = 0.5 * (low + high);
        if evaluate(t).0 > 0.0 {
            high = t;
        } else {
            low = t;
        }
    }
    let (residual, result) = evaluate(0.5 * (low + high));
    assert!(residual.abs() < 1e-10);
    result
}

#[test]
fn g1_g3_finite_coolant_changes_native_current_memory_and_pinning_heat() {
    let run = |viscosity, cp, pinning| {
        let (core_card, core_curve) = core_material(pinning);
        let (winding_card, winding_curve) = body_material();
        let (fluid_card, fluid_curve) = fluid_material(viscosity, cp);
        let core_context = body_point(300.0);
        let winding_context = body_point(310.0);
        let fluid_context = fluid_point(290.0);
        let thermal = initial(&core_curve, &winding_curve, &fluid_curve, cp);
        let fluid = liquid_law(&fluid_card);
        let cells = [
            (CORE, 0, &core_card, &core_curve, &core_context),
            (BODY, 1, &winding_card, &winding_curve, &winding_context),
            (FLUID, 2, &fluid_card, &fluid_curve, &fluid_context),
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
                branch: 1,
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
        let bases = [(1, ElectricalPropertyBasis::Resistivity)];
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
        let admitted = evaluate_source_link_convection(
            &thermal,
            &links,
            1,
            |_, _, t| fluid.at_temperature(t),
            || false,
        )
        .unwrap();
        let prepared_bulk = bulk
            .prepare_with_link_conductances(&thermal, &bases, admitted.conductances_w_k(), || false)
            .unwrap();
        let source = ResolvedScalarHysteresis::resolve(
            Arc::new(core_card.clone()),
            &core_context,
            roles(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let core = UniformHystereticCore::from_source(source, 1.0, 1.0, 1.0).unwrap();
        let start = HystereticCircuitInitial {
            core: core
                .candidate_current(&core.initial_state(), 0.0)
                .unwrap()
                .material
                .state,
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![VOLTAGE],
            time_s: 0.0,
        };
        let circuit = HystereticCircuit::new(
            CircuitGraph {
                node_count: 3,
                branches: vec![
                    (1, 0, Branch::VoltageSource { port: 0 }),
                    (1, 2, Branch::Resistor { ohms: f64::NAN }),
                    (2, 0, Branch::Inductor { henries: 999.0 }),
                ],
                transformers: vec![],
            },
            2,
            core,
            prepared_bulk
                .resolve_with_link_conductances(&thermal, admitted.conductances_w_k(), || false)
                .unwrap()
                .resistors,
            start,
        )
        .unwrap();
        assert!(
            ThermalHystereticCircuit::bind(
                circuit.clone(),
                &thermal,
                &bulk,
                &bases,
                CORE,
                [280.0, 320.0],
                || false
            )
            .is_err()
        );
        let driver = ThermalHystereticCircuit::bind_with_link_conductances(
            circuit.clone(),
            &thermal,
            &bulk,
            &bases,
            CORE,
            [280.0, 320.0],
            admitted.conductances_w_k(),
            || false,
        )
        .unwrap();
        let before_cells = thermal.cells().clone();
        let before_circuit = circuit.clone();
        let gate = fs_exec::CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let mut result = [0.0; 3];
        pool.scope(|arena| {
            let cx = fs_exec::Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 191,
                    kernel_id: 37,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            let prepared = heat
                .prepare_with_link_conductances(&cx, &thermal, admitted.conductances_w_k())
                .unwrap();
            let network = SourceConvectionHeatNetwork {
                heat: &prepared,
                cells: &[],
                maximum_cells: 0,
                resolve_fluid: |_: &SourceConvectiveCell,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    fluid.at_temperature(t)
                },
            };
            let assembly = network.with_links(&links, 1, |_, _, t| fluid.at_temperature(t));
            let mut solve = controls();
            solve.maximum_cells = 3;
            solve.maximum_circuit_trials = 4000;
            solve.newton.linear_restart = 3;
            // An initialization conductance is never a silently frozen coolant law.
            assert!(
                driver
                    .step_candidate(
                        &thermal,
                        DT,
                        &[VOLTAGE],
                        solve,
                        magnetic_controls(),
                        &[0.0; 3],
                        || false
                    )
                    .is_err()
            );
            let (mut accepted, mut magnetic) = (thermal.clone(), driver);
            let mut state = [300.0, 310.0, 290.0, 0.0, 0.0, 0.0, 0.0];
            for tick in 1..=2 {
                let expected = reference(state, viscosity, cp, pinning);
                let old_cells = accepted.cells().clone();
                let old_magnetic = magnetic.circuit().clone();
                let next = magnetic
                    .step_candidate_with_heat_assembly(
                        &cx,
                        &assembly,
                        &accepted,
                        f64::from(tick) * DT,
                        &[VOLTAGE],
                        solve,
                        magnetic_controls(),
                        &[0.0; 3],
                    )
                    .unwrap();
                let report = next.report();
                let step = &report.electrical;
                for (cell, wanted) in [CORE, BODY, FLUID].into_iter().zip(expected) {
                    near(next.thermal().temperature_k(cell).unwrap(), wanted, 3e-7);
                }
                near(next.circuit().circuit().current_a(), expected[3], 3e-8);
                near(
                    next.circuit().circuit().flux_linkage_wb(),
                    expected[4],
                    3e-8,
                );
                near(
                    next.circuit().circuit().magnetic_state().polarizations_t()[0],
                    expected[5],
                    3e-8,
                );
                near(step.core.stored_energy_j, expected[6], 3e-8);
                near(step.circuit.resistors[0].current_a, expected[7], 3e-8);
                near(step.core.core_heat_j, expected[8], 3e-8);
                near(step.circuit.heat_j, expected[9], 3e-8);
                near(step.supplied_work_j, expected[11], 3e-8);
                near(report.deposited_joule_j, expected[8] + expected[9], 3e-8);
                near(report.non_electrical_heat_j, 0.0, 1e-12);
                near(
                    report.thermal.stored_energy_change_j,
                    expected[8] + expected[9],
                    3e-8,
                );
                near(
                    next.thermal().energy_j() - accepted.energy_j() + step.stored_change_j
                        - step.supplied_work_j,
                    0.0,
                    3e-8,
                );
                near(report.total_energy_residual_j, 0.0, 3e-8);
                let coolant = &report.source.finite_convection.as_ref().unwrap().samples()[0];
                near(coolant.heat_into_fluid_w, expected[10], 3e-6);
                assert_eq!(
                    coolant.material.material().card_identity(),
                    fluid_card.content_hash()
                );
                assert_eq!(coolant.material.material().query().axes()["p"], 100000.0);
                let source = next.circuit().circuit().core().source().unwrap().material();
                assert_eq!(source.card_identity(), core_card.content_hash());
                assert_eq!(
                    source.query(),
                    body_point(next.thermal().temperature_k(CORE).unwrap())
                );
                let conductor = &step.circuit.resistors[0].material.conductor;
                assert_eq!(
                    conductor.material().card_identity(),
                    winding_card.content_hash()
                );
                near(
                    conductor.material().query().axes()["T"],
                    next.thermal().temperature_k(BODY).unwrap(),
                    1e-12,
                );
                assert_eq!(accepted.cells(), &old_cells);
                assert_eq!(magnetic.circuit(), &old_magnetic);
                assert_eq!(next.thermal().links(), thermal.links());
                assert!(step.core.core_heat_j > 0.0);
                assert!(
                    next.circuit().circuit().magnetic_state().polarizations_t()[0]
                        > old_magnetic.magnetic_state().polarizations_t()[0]
                );
                result = [
                    next.thermal().temperature_k(BODY).unwrap(),
                    next.circuit().circuit().current_a(),
                    next.thermal().temperature_k(CORE).unwrap(),
                ];
                state = expected[..7].try_into().unwrap();
                (accepted, magnetic, _) = next.into_parts();
            }
            let saved_cells = accepted.cells().clone();
            let saved_magnetic = magnetic.circuit().clone();
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
                    Err("magnetic coolant trial refused")
                } else {
                    Ok(fluid.at_temperature(t).unwrap())
                }
            });
            assert!(matches!(
                magnetic.step_candidate_with_heat_assembly(
                    &cx,
                    &failing,
                    &accepted,
                    3.0 * DT,
                    &[VOLTAGE],
                    solve,
                    magnetic_controls(),
                    &[0.0; 3]
                ),
                Err(CircuitEnthalpyError::Source(HystereticThermalError::Heat(
                    ConvectionHeatError::Convection(SourceCellConvectionError::Source {
                        cell: FLUID,
                        source: "magnetic coolant trial refused",
                        ..
                    })
                )))
            ));
            assert!(queries.get() > 1);
            let mut insufficient = solve;
            insufficient.maximum_circuit_trials = 1;
            assert!(matches!(
                magnetic.step_candidate_with_heat_assembly(
                    &cx,
                    &assembly,
                    &accepted,
                    3.0 * DT,
                    &[VOLTAGE],
                    insufficient,
                    magnetic_controls(),
                    &[0.0; 3]
                ),
                Err(CircuitEnthalpyError::CircuitBudget)
            ));
            gate.request();
            assert!(matches!(
                magnetic.step_candidate_with_heat_assembly(
                    &cx,
                    &assembly,
                    &accepted,
                    3.0 * DT,
                    &[VOLTAGE],
                    solve,
                    magnetic_controls(),
                    &[0.0; 3]
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            assert_eq!(accepted.cells(), &saved_cells);
            assert_eq!(magnetic.circuit(), &saved_magnetic);
            assert_eq!(thermal.cells(), &before_cells);
            assert_eq!(circuit, before_circuit);
        });
        result
    };
    let baseline = run(1.0, 4000.0, 2.0);
    let viscous = run(4.0, 4000.0, 2.0);
    let higher_capacity = run(1.0, 8000.0, 2.0);
    let weaker_pinning = run(1.0, 4000.0, 1.0);
    assert!(viscous[0] > baseline[0] && viscous[1] < baseline[1] && viscous[2] < baseline[2]);
    assert!(
        higher_capacity[0] < baseline[0]
            && higher_capacity[1] > baseline[1]
            && higher_capacity[2] > baseline[2]
    );
    assert!(weaker_pinning[1] < baseline[1] && weaker_pinning[2] < baseline[2]);
}
