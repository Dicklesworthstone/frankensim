//! Synthetic chemistry and source cards through the original finite network.
use super::*;
use crate::cavity_phs::piston::reactive::network::ReactiveNetworkError;
use crate::cavity_phs::piston::thermal::HeatExchangeSolve;
use crate::cavity_phs::piston::thermal::boundary::network::material::GasConductionGeometry;
use crate::cavity_phs::piston::thermal::boundary::network::{NetworkGasDrive, NetworkGasSolve};
use crate::circuit_enthalpy::material::{
    CellResistor, ConductionLink, MaterialCell, MaterialNetwork,
};
use crate::circuit_enthalpy::{CircuitEnthalpySolve, ResistorCell};
use crate::material_circuit::{MaterialCircuit, MaterialCircuitSolve};
use fs_matdb::SelectionPolicy;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::SPECIFIC_HEAT_CAPACITY_DIMS;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_ports::{ModuleWindowAudit, PortedModule};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

const BASES: [(usize, ElectricalPropertyBasis); 2] = [
    (2, ElectricalPropertyBasis::Resistivity),
    (1, ElectricalPropertyBasis::Conductivity),
];
fn id(i: usize) -> CellId {
    CellId(if i == 0 { 30 } else { 10 })
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", t)
        .unwrap()
        .with("lot", 7.0)
        .unwrap()
}
struct Fixture {
    cards: [MaterialCard; 2],
    caloric: [ResolvedHeatCapacityCurve; 2],
    context: QueryPoint,
    conductivity: PropertyKey,
}
impl Fixture {
    fn new(k: f64, sigma: f64) -> Self {
        let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cp = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
        let density = PropertyKey::new("density", Density::DIMS);
        let cards = [0, 1].map(|i| {
            let basis = BASES[1 - i].1;
            let mut claims = ClaimSet::new();
            for (key, a, b) in [
                (
                    PropertyKey::new(basis.property_name(), basis.dimensions()),
                    if i == 0 { sigma } else { 0.01 },
                    if i == 0 { 0.5 * sigma } else { 0.02 },
                ),
                (
                    conductivity.clone(),
                    if i == 0 { k } else { 0.4 },
                    if i == 0 { 2.0 * k } else { 0.8 },
                ),
                (cp.clone(), 500.0, 500.0),
                (density.clone(), 4000.0, 4000.0),
            ] {
                let dims = key.quantity().dims();
                claims
                    .insert_claim(PropertyClaim {
                        key,
                        value: PropertyValue::Curve {
                            abscissa: "T".into(),
                            abscissa_dims: Temperature::DIMS,
                            knots: vec![(200.0, a), (1000.0, b)],
                            dims,
                        },
                        validity: ValidityDomain::unconstrained()
                            .with("T", 200.0, 1000.0)
                            .with("lot", 7.0, 7.0),
                        uncertainty: UncertaintyModel::Unstated,
                        interpolation: InterpolationPolicy::LinearInside,
                        observations: vec![],
                        provenance: provenance(),
                    })
                    .unwrap();
            }
            MaterialCard::assemble(
                MaterialStateId {
                    chemistry: format!("synthetic-reactive-wall-{i}"),
                    phase: "solid".into(),
                    process: "test".into(),
                    revision: 0,
                },
                claims,
                vec![],
            )
            .unwrap()
        });
        let caloric = [0, 1].map(|i| {
            ResolvedHeatCapacityCurve::try_new(
                &cards[i],
                HeatCapacitySourceQuery {
                    heat_capacity: &cp,
                    density: &density,
                    temperature_axis: "T",
                    lower: &point(200.0),
                    upper: &point(1000.0),
                    selection: SelectionPolicy::SingleClaimOnly,
                },
                SolidLiquidPhase::Solid,
                0.0,
                1e-5,
            )
            .unwrap()
        });
        Self {
            cards,
            caloric,
            context: point(450.0),
            conductivity,
        }
    }
    fn cells(&self) -> [MaterialCell<'_>; 2] {
        [0, 1].map(|i| MaterialCell {
            cell: id(i),
            material_id: MaterialId(i as u64),
            curve: self.caloric[i].curve(),
            point: ThermalMaterialPoint {
                material: &self.cards[i],
                context: &self.context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: Some((
                &self.conductivity,
                MaterialPropertySelection::SingleClaimOnly,
            )),
        })
    }
    fn state(&self) -> DepositionSimulation<'_> {
        let mut events = [0, 1]
            .map(|i| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: id(i),
                    material_id: MaterialId(i as u64),
                    segment_id: i as u64,
                    body_id: i as u64,
                    mass_kg: 0.001,
                    specific_h_j_kg: if i == 0 { 125000.0 } else { 100000.0 },
                    geometry: CellGeometry {
                        start_m: [i as f64, 0.0, 0.0],
                        end_m: [i as f64 + 0.01, 0.0, 0.0],
                        volume_m3: 0.001 / 4000.0,
                        surface_area_m2: 0.04,
                    },
                    boundary: Boundary {
                        ambient_k: 400.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 2.0,
                    },
                }),
            })
            .to_vec();
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 7,
                a: id(0),
                b: id(1),
                area_m2: 0.002,
                conductance_w_k: 0.1,
            }),
        });
        DepositionSimulation::new(
            events,
            &[
                (MaterialId(0), self.caloric[0].curve() as &dyn ThermalLaw),
                (MaterialId(1), self.caloric[1].curve() as &dyn ThermalLaw),
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
                max_dt_s: 0.001,
                energy_tolerance_j: 1e-7,
            },
            || Ok(()),
        )
        .unwrap()
        .0
    }
    fn piston(&self, reaction: &GasReaction, dilution: f64) -> MaterialGasPiston {
        let body = resolve_body(
            &[RegionBinding {
                key: "piston",
                material: &self.cards[0],
                point: &self.context,
                geometry: VolumeMoments::cuboid(Vec3::new(0.1, 0.1, 0.025)).unwrap(),
                selection: MaterialPropertySelection::SingleClaimOnly,
            }],
            1,
            || false,
        )
        .unwrap();
        MaterialGasPiston::new(
            gas(reaction, [1.0 - dilution, 0.0, dilution]),
            body,
            0.01,
            0.0,
        )
        .unwrap()
    }
}
fn resistors() -> [CellResistor; 2] {
    [1, 0].map(|i| CellResistor {
        destination: ResistorCell {
            branch: i + 1,
            cell: id(i),
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 0.01,
        area_m2: 0.000025,
    })
}
fn links() -> [ConductionLink; 1] {
    [ConductionLink {
        link: 7,
        distance_a_m: 0.001,
        distance_b_m: 0.001,
    }]
}
fn circuit(binding: &MaterialNetwork<'_>, state: &DepositionSimulation<'_>) -> MaterialCircuit {
    MaterialCircuit::new(
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: 0.01 }),
                (1, 2, Branch::Resistor { ohms: 0.0 }),
                (2, 0, Branch::Resistor { ohms: 0.0 }),
            ],
            transformers: vec![],
        },
        binding
            .resolve_with_bases(state, &BASES, || false)
            .unwrap()
            .resistors,
        vec![],
        vec![0.2],
        &[],
    )
    .unwrap()
}
fn controls() -> NetworkGasSolve {
    NetworkGasSolve {
        network: CircuitEnthalpySolve {
            maximum_cells: 2,
            maximum_circuit_trials: 300,
            maximum_newton_iterations: 32,
            difference_step_j: 1e-4,
            cell_energy_tolerance_j: 1e-7,
            thermal_energy_tolerance_j: 1e-7,
            total_energy_tolerance_j: 1e-7,
            circuit: MaterialCircuitSolve {
                energy_tolerance_j: 1e-7,
                power_tolerance_w: 1e-7,
            },
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-8,
                relative_tolerance: 1e-12,
                forcing_minimum: 1e-4,
                linear_restart: 4,
                ..NewtonKrylovConfig::default()
            },
        },
        maximum_boundary_solves: 300,
        total_energy_tolerance_j: 1e-6,
    }
}
fn network_drive(dt: f64) -> NetworkGasDrive {
    NetworkGasDrive {
        cell: id(0),
        input: PistonInput {
            duration_s: dt,
            heat_j: 3.0 * dt,
            external_force_n: 0.2,
            ..drive()
        },
        piston: mechanical_controls(),
        exchange: HeatExchangeSolve {
            heat_bracket_j: [-1.0, 1.0],
            source_tolerance_j: 1e-9,
            energy_tolerance_j: 1e-7,
            maximum_iterations: 128,
        },
    }
}
fn geometry() -> GasConductionGeometry {
    GasConductionGeometry {
        area_m2: 0.01,
        gas_distance_m: 0.0001,
        cell_distance_m: 0.002,
    }
}

#[test]
fn material_heat_resistance_and_reaction_products_share_one_finite_energy_account() {
    let run = |k, sigma, dilution, rate| {
        let fixture = Fixture::new(k, sigma);
        let cells = fixture.cells();
        let resistors = resistors();
        let links = links();
        let binding = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &links,
        };
        let reaction = reaction(8000.0, 4.0, rate);
        let mut thermal = fixture.state();
        let mut electrical = circuit(&binding, &thermal);
        let mut piston = fixture.piston(&reaction, dilution);
        let initial = piston.clone();
        let initial_h = thermal.energy_j();
        let mut supplied = 0.0;
        let mut exchanged = 0.0;
        let mut extent = 0.0;
        for n in 1..=16 {
            let target = f64::from(n) * drive().duration_s;
            let load = network_drive(target - thermal.time_s());
            let single = reaction
                .parcel_candidate(
                    piston.gas(),
                    load.input.duration_s,
                    chemical_controls(),
                    || false,
                )
                .unwrap();
            let candidate = piston
                .reactive_material_network_candidate(
                    &reaction,
                    &binding,
                    &thermal,
                    target,
                    &electrical,
                    &[],
                    &BASES,
                    load,
                    chemical_controls(),
                    controls(),
                    geometry(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            let (next_t, next_c, next_p, report) = candidate.into_parts();
            assert_eq!(
                &report.reaction,
                single.step(),
                "Newton probes cannot react again"
            );
            assert_eq!(next_p.body(), initial.body());
            assert_eq!(next_p.gas().mass_kg(), initial.gas().mass_kg());
            assert!(report.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j);
            let stage = &report.network;
            let cell = stage
                .network
                .source
                .cells
                .iter()
                .find(|c| c.cell == id(0))
                .unwrap();
            let resistor = stage
                .network
                .electrical
                .resistors
                .iter()
                .find(|r| r.material.branch == 1)
                .unwrap();
            assert_eq!(cell.conductivity.as_ref(), Some(&stage.interface));
            assert_eq!(resistor.material.conductor.material(), &stage.interface);
            assert!(stage.interface.query_point().contains(&("lot".into(), 7.0)));
            assert!(stage.network.deposited_joule_j > 0.0);
            let motion = &stage.boundary.piston;
            supplied += stage.coupled.external_network_heat_j
                + load.input.heat_j
                + stage.network.electrical.native.supplied
                + motion.external_force_work_j
                - motion.exterior_pressure_work_j;
            exchanged += stage.boundary.exchanged_heat_j;
            extent += report.reaction.extent_mol;
            let ModuleWindowAudit::Available(window) = next_p.piston().audit_last_window().unwrap()
            else {
                panic!("the closed reacting chamber must retain its complete heat/work account");
            };
            assert!(window.violations().is_empty(), "{:?}", window.violations());
            assert_eq!(
                next_p.piston().ports.next_tick(),
                piston.piston().ports.next_tick() + 1
            );
            thermal = next_t;
            electrical = next_c;
            piston = next_p;
        }
        let gas_change = piston.gas().mass_kg()
            * (piston.gas().state().specific_internal_energy_j_kg()
                - initial.gas().state().specific_internal_energy_j_kg());
        let electric_change = (electrical.charges_c()[0].powi(2) - 0.2_f64.powi(2)) / (2.0 * 0.01);
        let kinetic = 0.5 * piston.piston().mass_kg() * piston.piston().velocity_m_s().powi(2);
        near(
            gas_change + thermal.energy_j() - initial_h + electric_change + kinetic - supplied,
            0.0,
            16.0 * controls().total_energy_tolerance_j,
        );
        [
            electrical.charges_c()[0],
            thermal.temperature_k(id(0)).unwrap(),
            exchanged,
            piston.piston().velocity_m_s(),
            piston.gas().state().state().sound_speed,
            extent,
        ]
    };
    let reference = run(0.2, 50.0, 0.25, 1000.0);
    for changed in [
        run(2.0, 50.0, 0.25, 1000.0),
        run(0.2, 100.0, 0.25, 1000.0),
        run(0.2, 50.0, 0.5, 1000.0),
        run(0.2, 50.0, 0.25, 500.0),
    ] {
        assert!(
            (changed[2] - reference[2]).abs() > 1e-10,
            "wall heat must follow material state"
        );
        assert!(
            (changed[3] - reference[3]).abs() > 1e-12,
            "motion must follow coupled state"
        );
        assert!(
            (changed[4] - reference[4]).abs() > 1e-9,
            "acoustic state must follow coupled state"
        );
    }
    assert!(reference[5] > 0.0);
}

#[test]
fn absent_reagents_reduce_to_existing_network_and_refusals_keep_every_owner() {
    let fixture = Fixture::new(0.2, 50.0);
    let cells = fixture.cells();
    let resistors = resistors();
    let links = links();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &links,
    };
    let reaction = reaction(8000.0, 4.0, 1000.0);
    let thermal = fixture.state();
    let electrical = circuit(&binding, &thermal);
    let piston = fixture.piston(&reaction, 1.0);
    let load = network_drive(drive().duration_s);
    let saved_t = thermal.cells().clone();
    let saved_p = piston.clone();
    let saved_c = electrical.clone();
    let actual = piston
        .reactive_material_network_candidate(
            &reaction,
            &binding,
            &thermal,
            load.input.duration_s,
            &electrical,
            &[],
            &BASES,
            load,
            chemical_controls(),
            controls(),
            geometry(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let expected = piston
        .piston()
        .material_network_candidate_with_bases(
            &binding,
            &thermal,
            load.input.duration_s,
            &electrical,
            &[],
            &BASES,
            load,
            controls(),
            geometry(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let (at, ac, ap, ar) = actual.into_parts();
    let (et, ec, ep, er) = expected.into_parts();
    assert_eq!(ap.gas(), ep.gas());
    assert_eq!(ap.piston().velocity_m_s(), ep.velocity_m_s());
    assert_eq!(ap.piston().mass_kg(), ep.mass_kg());
    assert_eq!(ap.piston().position_m(), ep.position_m());
    assert_eq!(ap.piston().ports.next_tick(), ep.ports.next_tick());
    assert_eq!(at.cells(), et.cells());
    assert_eq!(ac, ec);
    assert_eq!(ar.network.boundary, er.boundary);
    assert_eq!(ar.reaction.extent_mol, 0.0);

    for (target, cancelled, bad_source) in [
        (2.0 * load.input.duration_s, false, false),
        (load.input.duration_s, true, false),
        (load.input.duration_s, false, true),
    ] {
        let bases = if bad_source { &BASES[..1] } else { &BASES[..] };
        let result = piston.reactive_material_network_candidate(
            &reaction,
            &binding,
            &thermal,
            target,
            &electrical,
            &[],
            bases,
            load,
            chemical_controls(),
            controls(),
            geometry(),
            &[0.0; 2],
            || cancelled,
        );
        assert!(result.is_err());
        if cancelled {
            assert!(matches!(result, Err(ReactiveNetworkError::Cancelled)));
        }
        assert_eq!(thermal.cells(), &saved_t);
        assert_eq!(piston, saved_p);
        assert_eq!(electrical, saved_c);
    }
}
