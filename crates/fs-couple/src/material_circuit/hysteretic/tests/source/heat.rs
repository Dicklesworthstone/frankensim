//! G1/G3/G4: original magnetic memory, finite heat transfer and source-only swaps.
use super::*;
use crate::circuit_enthalpy::{
    contacts::{
        heat::{MaterialHeatNetwork, PreparedHeatNetwork},
        material::ElectricalContactNetwork,
    },
    material::ConductionLink,
};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
    EMISSIVITY_DIMS, GrayDiffuseEnclosure, RadiationSurface, STEFAN_BOLTZMANN_W_M2_K4,
    SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity, ViewFactorMatrix,
    transient::deposition::{contact::*, radiation::*},
};
use fs_exec::{CancelGate, Cx};
use fs_matdb::{ClaimSelection, InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::{
    state_point::THERMAL_CONDUCTIVITY_DIMS, thermal::interface::ThermalInterfacePoint,
};

const CORE: CellId = CellId(10);
const COIL: CellId = CellId(11);
const AREA: f64 = 1.0;
const DT: f64 = 0.05;

fn heat_card(original: MaterialCard, emissivity: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = original.claims().clone();
    for (name, dims, value) in [
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 1.0),
        (SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, emissivity),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured magnetic heat-transfer fixture".into(),
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
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
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

fn with_heat<R>(
    emissivity: f64,
    resistance: f64,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &PreparedHeatNetwork<'_, '_>,
        &DepositionSimulation<'_>,
        &ThermalHystereticCircuit<'_, '_>,
        &MaterialNetwork<'_>,
    ) -> R,
) -> R {
    let (core_card, _) = core_material(2.0, 310.0);
    let (coil_card, _) = material(1.0, 0.005);
    let (a, ca) = heat_card(core_card, emissivity);
    let (b, cb) = heat_card(coil_card, emissivity);
    let cards = [a, b];
    let curves = [ca, cb];
    let original = thermal_initial(&curves);
    let thermal = DepositionSimulation::new(
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(original.cells()[&CORE].clone()),
            },
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(original.cells()[&COIL].clone()),
            },
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 0,
                    a: CORE,
                    b: COIL,
                    area_m2: AREA,
                    conductance_w_k: 999.0,
                }),
            },
        ],
        &[
            (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
        ],
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
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let mut bindings = cells(&cards, &curves, &context);
    for cell in &mut bindings {
        cell.conductivity = Some((&key, MaterialPropertySelection::SingleClaimOnly));
    }
    let resistors = winding();
    let links = [ConductionLink {
        link: 0,
        distance_a_m: 0.05,
        distance_b_m: 0.05,
    }];
    let network = MaterialNetwork {
        cells: &bindings,
        resistors: &resistors,
        links: &links,
    };
    let sources = network.resolve(&thermal, || false).unwrap().resistors;
    let core = sourced_core(&cards[0]);
    // This interval stays on the original sliding branch: mu+chi=0.21 H.
    let start = initial(&core, 3.0, vec![5.0], vec![]);
    let circuit = HystereticCircuit::new(rl_graph(), 2, core, sources, start).unwrap();
    let bases = [(1, ElectricalPropertyBasis::Resistivity)];
    let driver = ThermalHystereticCircuit::bind(
        circuit,
        &thermal,
        &network,
        &bases,
        CORE,
        [290.0, 310.0],
        || false,
    )
    .unwrap();

    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(RP, RD),
            value: PropertyValue::Scalar {
                value: resistance,
                dims: RD,
            },
            validity: ValidityDomain::unconstrained()
                .with("aT", 290.0, 310.0)
                .with("bT", 290.0, 310.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured core/winding thermal interface".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let interface = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "core".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "winding".into(),
        },
        SystemContext {
            medium: "fixed".into(),
            third_body: None,
            environment: "fixture".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let ip = QueryPoint::new()
        .with("aT", 300.0)
        .unwrap()
        .with("bT", 300.0)
        .unwrap();
    let contacts = [MaterialContactLink {
        link: 0,
        sides: [0, 1].map(|i| ContactCell {
            cell: bindings[i].cell,
            material_id: bindings[i].material_id,
            curve: curves[i].curve(),
        }),
        point: ThermalInterfacePoint {
            interface: &interface,
            materials: [&cards[0], &cards[1]],
            context: &ip,
            temperature_axes: ["aT", "bT"],
            liquid_fraction_axes: [None; 2],
            mean_temperature_axis: None,
        },
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }];
    let electrical = ElectricalContactNetwork {
        bulk: &network,
        contacts: &[],
        electrical_bases: &bases,
    };
    let (complex, positions) = fs_conduction::fixtures::box_grid([1; 3], [1.0; 3]);
    let mesh = fs_conduction::mesh::ConductionMesh::new(complex, positions).unwrap();
    let surfaces = [0, 1]
        .map(|i| {
            RadiationSurface::new(
                &mesh,
                format!("face-{i}"),
                |face| fs_conduction::fixtures::on_box_face(face.centroid[0], i as f64),
                SurfaceEmissivity::from_thermal_state(
                    "fixture",
                    bindings[i].point,
                    curves[i].curve().state_at_specific_enthalpy(100.0).unwrap(),
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .unwrap()
                .0,
            )
            .unwrap()
        })
        .to_vec();
    // Algebraic parallel-plate limit; cube faces supply area, not cube view factors.
    let enclosure = GrayDiffuseEnclosure::new(
        surfaces,
        ViewFactorMatrix::infinite_parallel_plates(AREA).unwrap(),
    )
    .unwrap();
    let radiating = [0, 1].map(|i| RadiatingCell {
        surface: i,
        cell: bindings[i].cell,
        material_id: bindings[i].material_id,
        curve: curves[i].curve(),
        material: bindings[i].point,
        selection: MaterialPropertySelection::SingleClaimOnly,
    });
    let radiation = MaterialRadiation {
        enclosure: &enclosure,
        cells: &radiating,
        power_tolerance_w: 1e-9,
        radiosity_tolerance_w_m2: 1e-9,
    };
    let assembly = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &contacts,
        radiation: Some(&radiation),
    };
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 20,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let prepared = assembly.prepare(&cx, &thermal).unwrap();
        run(&cx, &gate, &prepared, &thermal, &driver, &network)
    })
}

// Independent scalar balance: eliminate core temperature with total finite
// caloric energy and eliminate current with the active-branch RL formula.
// No production material evaluator or coupled residual is used here.
fn reference(emissivity: f64, resistance: f64) -> [f64; 5] {
    let evaluate = |coil_t: f64| {
        let r = 1.0 + 0.005 * (coil_t - 300.0);
        let current = (3.0 * (0.21 - 0.5 * DT * r) + DT * 5.0) / (0.21 + 0.5 * DT * r);
        let mean = 0.5 * (3.0 + current);
        let joule = DT * r * mean * mean;
        let pinning = 0.4 * (current - 3.0);
        let core_t = 300.0 + joule + pinning - (coil_t - 300.0);
        let exchange = DT
            * (AREA / (0.1 + resistance) * (coil_t - core_t)
                + STEFAN_BOLTZMANN_W_M2_K4 * AREA * (coil_t.powi(4) - core_t.powi(4))
                    / (2.0 / emissivity - 1.0));
        (
            [core_t, coil_t, current, joule, pinning],
            coil_t - 300.0 - joule + exchange,
        )
    };
    let (mut low, mut high) = (300.0, 301.0);
    assert!(evaluate(low).1 < 0.0 && evaluate(high).1 > 0.0);
    for _ in 0..80 {
        let middle = 0.5 * (low + high);
        if evaluate(middle).1 > 0.0 {
            high = middle;
        } else {
            low = middle;
        }
    }
    evaluate(0.5 * (low + high)).0
}

#[test]
fn sourced_contact_and_radiation_change_current_and_both_finite_heat_stores() {
    let run = |emissivity, resistance| {
        with_heat(
            emissivity,
            resistance,
            |cx, _, prepared, thermal, driver, _| {
                let saved = driver.circuit().clone();
                let cells_before = thermal.cells().clone();
                let next = driver
                    .step_candidate_with_heat_network(
                        cx,
                        prepared,
                        thermal,
                        DT,
                        &[5.0],
                        heat_controls(),
                        controls(),
                        &[0.0; 2],
                    )
                    .unwrap();
                let expected = reference(emissivity, resistance);
                let report = next.report();
                let step = &report.electrical;
                let core_t = next.thermal().temperature_k(CORE).unwrap();
                let coil_t = next.thermal().temperature_k(COIL).unwrap();
                for (actual, expected) in [
                    (core_t, expected[0]),
                    (coil_t, expected[1]),
                    (next.circuit().circuit().current_a(), expected[2]),
                    (step.circuit.heat_j, expected[3]),
                    (step.core.core_heat_j, expected[4]),
                ] {
                    near(actual, expected, 3e-8);
                }
                near(
                    report.deposited_joule_j,
                    step.circuit.heat_j + step.core.core_heat_j,
                    1e-11,
                );
                near(
                    next.thermal().energy_j() - thermal.energy_j(),
                    report.deposited_joule_j,
                    3e-8,
                );
                near(report.non_electrical_heat_j, 0.0, 1e-8);
                near(
                    report.total_energy_residual_j,
                    0.0,
                    heat_controls().total_energy_tolerance_j,
                );
                near(
                    step.stored_change_j + report.deposited_joule_j,
                    step.supplied_work_j,
                    3e-8,
                );
                let source = &report.source;
                assert_eq!(source.thermal_contacts.len(), 1);
                near(
                    source.thermal_contacts[0].conductance_w_k,
                    AREA / (0.1 + resistance),
                    1e-12,
                );
                for surface in source.radiation.as_ref().unwrap().surfaces() {
                    assert_eq!(
                        surface.phase,
                        source
                            .material
                            .cells
                            .iter()
                            .find(|cell| cell.cell == surface.cell)
                            .unwrap()
                            .phase,
                    );
                }
                let accepted_core = next.circuit().circuit().core().source().unwrap().material();
                let original_core = driver.circuit().core().source().unwrap().material();
                assert_eq!(accepted_core.card_identity(), original_core.card_identity());
                near(
                    *accepted_core.query().axes().get("T").unwrap(),
                    core_t,
                    1e-10,
                );
                for property in original_core.properties() {
                    assert_eq!(
                        accepted_core
                            .property(property.requirement().name())
                            .unwrap()
                            .answer()
                            .receipt
                            .selected,
                        property.answer().receipt.selected,
                    );
                }
                assert_eq!(driver.circuit(), &saved);
                assert_eq!(thermal.cells(), &cells_before);
                assert_eq!(thermal.time_s(), 0.0);
                // The endpoint's original memory and source state can continue.
                let continued = next
                    .circuit()
                    .step_candidate_with_heat_network(
                        cx,
                        prepared,
                        next.thermal(),
                        2.0 * DT,
                        &[5.0],
                        heat_controls(),
                        controls(),
                        &[0.0; 2],
                    )
                    .unwrap();
                assert!(
                    continued.circuit().circuit().flux_linkage_wb()
                        > next.circuit().circuit().flux_linkage_wb()
                );
                assert_eq!(
                    continued.thermal().time_s(),
                    continued.circuit().circuit().time_s()
                );
                (core_t, coil_t, next.circuit().circuit().current_a())
            },
        )
    };
    let weak = run(0.1, 0.2);
    let contact = run(0.1, 0.001);
    let both = run(0.9, 0.001);
    assert!(both.0 > contact.0 && contact.0 > weak.0);
    assert!(both.1 < contact.1 && contact.1 < weak.1);
    assert!(both.2 > contact.2 && contact.2 > weak.2);
}

#[test]
fn foreign_assembly_source_exit_and_cancellation_preserve_magnetic_memory() {
    with_heat(
        0.8,
        0.001,
        |cx, gate, prepared, thermal, driver, network| {
            let saved = driver.circuit().clone();
            let cells_before = thermal.cells().clone();
            // Even a newly assembled equivalent network is not the retained owner.
            let foreign = MaterialNetwork {
                cells: network.cells,
                resistors: network.resistors,
                links: network.links,
            };
            let electrical = ElectricalContactNetwork {
                bulk: &foreign,
                contacts: &[],
                electrical_bases: &[(1, ElectricalPropertyBasis::Resistivity)],
            };
            let assembly = MaterialHeatNetwork {
                electrical: &electrical,
                thermal_contacts: &[],
                radiation: None,
            };
            let foreign = assembly.prepare(cx, thermal).unwrap();
            assert!(matches!(
                driver.step_candidate_with_heat_network(
                    cx,
                    &foreign,
                    thermal,
                    DT,
                    &[5.0],
                    heat_controls(),
                    controls(),
                    &[0.0; 2],
                ),
                Err(crate::circuit_enthalpy::CircuitEnthalpyError::Invalid(_)),
            ));
            // A rejected trial outside the original magnetic source interval must
            // not debit any heat or advance any stored polarization.
            assert!(
                driver
                    .step_candidate_with_heat_network(
                        cx,
                        prepared,
                        thermal,
                        DT,
                        &[5.0],
                        heat_controls(),
                        controls(),
                        &[20.0, 0.0],
                    )
                    .is_err()
            );
            assert_eq!(driver.circuit(), &saved);
            assert_eq!(thermal.cells(), &cells_before);
            let retry = || {
                driver
                    .step_candidate_with_heat_network(
                        cx,
                        prepared,
                        thermal,
                        DT,
                        &[5.0],
                        heat_controls(),
                        controls(),
                        &[0.0; 2],
                    )
                    .unwrap()
            };
            let a = retry();
            let b = retry();
            assert_eq!(a.circuit().circuit(), b.circuit().circuit());
            assert_eq!(a.thermal().cells(), b.thermal().cells());
            gate.request();
            assert!(matches!(
                driver.step_candidate_with_heat_network(
                    cx,
                    prepared,
                    thermal,
                    DT,
                    &[5.0],
                    heat_controls(),
                    controls(),
                    &[0.0; 2],
                ),
                Err(crate::circuit_enthalpy::CircuitEnthalpyError::Cancelled),
            ));
            assert_eq!(driver.circuit(), &saved);
            assert_eq!(thermal.cells(), &cells_before);
        },
    );
}
