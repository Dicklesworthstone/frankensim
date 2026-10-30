//! Real finite thermal, circuit, material and shared gas owners; synthetic cards.
use super::*;
use crate::cavity_phs::piston::double_acting::electrical::tests as motion_fixture;
use crate::cavity_phs::piston::thermal::double_acting::DoubleHeatExchangeSolve;
use crate::cavity_phs::piston::thermal::{HeatExchangeSolve, source::ThermalMaterialPoint};
use crate::circuit_enthalpy::material::{CellResistor, MaterialCell};
use crate::circuit_enthalpy::{CircuitEnthalpySolve, ResistorCell};
use fs_evidence::ValidityDomain;
use fs_matdb::*;
use fs_material::conductor::ELECTRICAL_RESISTIVITY_DIMS;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
};
use fs_phs::circuit::CircuitGraph;
use fs_qty::{Density, Temperature};
use fs_thermal::enthalpy::*;

pub(in crate::cavity_phs::piston::double_acting) fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap()
}
pub(in crate::cavity_phs::piston::double_acting) fn card(
    k: f64,
    r: f64,
    slope: f64,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = ClaimSet::new();
    for (name, dims, a, b) in [
        (
            "electrical_resistivity",
            ELECTRICAL_RESISTIVITY_DIMS,
            r * 0.0025 * (1.0 - 200.0 * slope),
            r * 0.0025 * (1.0 + 400.0 * slope),
        ),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, k, k),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            500.0,
            500.0,
        ),
        ("density", Density::DIMS, 4000.0, 4000.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, a), (800.0, b)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 200.0, 800.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic electromagnetic network".into(),
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
            process: "authored".into(),
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
            lower: &point(200.0),
            upper: &point(800.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-5,
    )
    .unwrap();
    (card, curve)
}
pub(in crate::cavity_phs::piston::double_acting) fn cells<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    p: &'a QueryPoint,
    key: &'a PropertyKey,
) -> [MaterialCell<'a>; 2] {
    [0, 1].map(|i| MaterialCell {
        cell: CellId(i as u64),
        material_id: MaterialId(i as u64),
        curve: curves[i].curve(),
        point: ThermalMaterialPoint {
            material: &cards[i],
            context: p,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: Some((key, MaterialPropertySelection::SingleClaimOnly)),
    })
}
pub(in crate::cavity_phs::piston::double_acting) fn initial<'a>(
    curves: &'a [ResolvedHeatCapacityCurve; 2],
) -> DepositionSimulation<'a> {
    let events = [0, 1]
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i),
                segment_id: i,
                body_id: i,
                material_id: MaterialId(i),
                mass_kg: 0.001,
                specific_h_j_kg: 100000.0,
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
                    source_w: 0.0,
                },
            }),
        })
        .to_vec();
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 2,
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
pub(in crate::cavity_phs::piston::double_acting) fn resistor() -> [CellResistor; 1] {
    [CellResistor {
        destination: ResistorCell {
            branch: 1,
            cell: CellId(0),
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 0.01,
        area_m2: 0.000025,
    }]
}
pub(in crate::cavity_phs::piston::double_acting) fn circuit(
    binding: &MaterialNetwork<'_>,
    thermal: &DepositionSimulation<'_>,
) -> MaterialCircuit {
    MaterialCircuit::new(
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: 0.01 }),
                (1, 2, Branch::Resistor { ohms: 999.0 }),
                (2, 0, Branch::VoltageSource { port: 0 }),
            ],
            transformers: vec![],
        },
        binding.resolve(thermal, || false).unwrap().resistors,
        vec![],
        vec![0.2],
        &[0.0],
    )
    .unwrap()
}
pub(in crate::cavity_phs::piston::double_acting) fn controls() -> NetworkGasSolve {
    NetworkGasSolve {
        network: CircuitEnthalpySolve {
            maximum_cells: 2,
            maximum_circuit_trials: 3000,
            maximum_newton_iterations: 32,
            difference_step_j: 1e-3,
            cell_energy_tolerance_j: 1e-7,
            thermal_energy_tolerance_j: 1e-7,
            total_energy_tolerance_j: 1e-7,
            circuit: motion_fixture::solve().circuit,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-8,
                relative_tolerance: 1e-14,
                linear_restart: 4,
                forcing_minimum: 1e-4,
                ..NewtonKrylovConfig::default()
            },
        },
        maximum_boundary_solves: 3000,
        total_energy_tolerance_j: 1e-6,
    }
}
pub(in crate::cavity_phs::piston::double_acting) fn drive(
    area: f64,
) -> ElectromagneticNetworkDrive {
    ElectromagneticNetworkDrive {
        gas: DoubleNetworkDrive {
            cells: [CellId(0), CellId(1)],
            input: motion_fixture::input(),
            piston: motion_fixture::mechanical(),
            exchange: DoubleHeatExchangeSolve {
                faces: [HeatExchangeSolve {
                    heat_bracket_j: [-1.0, 1.0],
                    source_tolerance_j: 1e-10,
                    energy_tolerance_j: 1e-7,
                    maximum_iterations: 128,
                }; 2],
                maximum_joint_trials: 20000,
            },
        },
        transducer: TranslationalCircuitPort {
            branch: 2,
            factor_n_a: 2.0,
        },
        motion: motion_fixture::solve(),
        initial_midpoint_velocity_m_s: 0.0,
        geometry: [GasConductionGeometry {
            area_m2: area,
            gas_distance_m: 0.001,
            cell_distance_m: 0.02,
        }; 2],
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-6, "{a} != {b}");
}

#[test]
fn actual_resistor_heat_enters_finite_storage_while_force_drives_the_gases() {
    let (a, ca) = card(0.2, 2.0, 0.0);
    let (b, cb) = card(0.4, 8.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let resistors = resistor();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let t = initial(&curves);
    let c = circuit(&binding, &t);
    let p = motion_fixture::piston(0.0);
    let candidate = p
        .electromagnetic_network_candidate(
            &binding,
            &t,
            0.00005,
            &c,
            &[0.0],
            drive(0.0),
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let direct = p
        .electrical_candidate(
            motion_fixture::electrical(&c, &[0.0], 2.0),
            motion_fixture::input(),
            motion_fixture::mechanical(),
            motion_fixture::solve(),
            0.0,
            || false,
        )
        .unwrap();
    near(
        candidate.piston().primary().velocity_m_s(),
        direct.piston().primary().velocity_m_s(),
    );
    let (tn, cn, _, r) = candidate.into_parts();
    near(tn.energy_j() - t.energy_j(), r.network.electrical.heat_j);
    near(tn.temperature_k(CellId(1)).unwrap(), 400.0);
    near(cn.charges_c()[0], direct.circuit().circuit().charges_c()[0]);
    near(r.coupled.external_network_heat_j, 0.0);
    near(r.coupled.total_energy_residual_j, 0.0);
    assert_eq!(r.boundary.exchanged_heat_j, [0.0; 2]);
    assert!(r.motion.mechanical_work_j > 0.0 && r.network.electrical.heat_j > 0.0);
    assert_eq!(r.network.circuit_trials, r.coupled.boundary_solves);
}

#[test]
fn temperature_dependent_resistivity_and_two_finite_contacts_close_all_storage() {
    let run = |slope| {
        let (a, ca) = card(0.2, 2.0, slope);
        let (b, cb) = card(0.4, 8.0, 0.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let context = point(400.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &context, &key);
        let resistors = resistor();
        let binding = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &[],
        };
        let t = initial(&curves);
        let c = circuit(&binding, &t);
        let p = motion_fixture::piston(0.0);
        let candidate = p
            .electromagnetic_network_candidate(
                &binding,
                &t,
                0.00005,
                &c,
                &[0.0],
                drive(0.01),
                controls(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let (tn, cn, pn, r) = candidate.into_parts();
        let energy = (r.network.electrical.stored_after_j - r.network.electrical.stored_before_j)
            + tn.energy_j()
            - t.energy_j()
            + r.boundary.piston.gas[0].stored_energy_change_j
            + r.boundary.piston.gas[1].stored_energy_change_j
            + r.boundary.piston.kinetic_energy_change_j;
        near(energy, 0.0);
        near(r.motion.transfer_residual_j, 0.0);
        near(r.network.deposited_joule_j, r.network.electrical.heat_j);
        for i in 0..2 {
            assert_eq!(
                r.network
                    .source
                    .cells
                    .iter()
                    .find(|s| s.cell == CellId(i as u64))
                    .unwrap()
                    .conductivity
                    .as_ref(),
                Some(&r.interfaces[i])
            );
        }
        (
            tn.temperature_k(CellId(0)).unwrap(),
            cn.charges_c()[0],
            pn.primary().velocity_m_s(),
        )
    };
    let fixed = run(0.0);
    let variable = run(0.002);
    assert!(variable.0 < fixed.0 && variable.1 > fixed.1 && variable.2 < fixed.2);
}

#[test]
fn common_cell_contacts_are_grouped_and_budget_or_late_cancel_publishes_nothing() {
    let (a, ca) = card(0.2, 2.0, 0.0);
    let (b, cb) = card(0.4, 8.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let resistors = resistor();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let t = initial(&curves);
    let c = circuit(&binding, &t);
    let p = motion_fixture::piston(0.0);
    let old = p.clone();
    let old_c = c.clone();
    let mut drive = drive(0.03);
    drive.gas.cells = [CellId(0); 2];
    assert!(
        p.electromagnetic_network_candidate(
            &binding,
            &t,
            0.00005,
            &c,
            &[0.0],
            drive,
            controls(),
            &[0.0; 2],
            || false
        )
        .is_err()
    );
    drive.geometry[0].area_m2 = 0.001;
    drive.geometry[1].area_m2 = 0.001;
    let ok = p
        .electromagnetic_network_candidate(
            &binding,
            &t,
            0.00005,
            &c,
            &[0.0],
            drive,
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let r = ok.network().report();
    let heat = r.source.boundary.step().exchanged_heat_j;
    near(
        ok.network().thermal().energy_j() - t.energy_j(),
        r.electrical.heat_j - heat[0] - heat[1],
    );
    let mut limited = controls();
    limited.network.maximum_circuit_trials = 1;
    assert!(
        p.electromagnetic_network_candidate(
            &binding,
            &t,
            0.00005,
            &c,
            &[0.0],
            drive,
            limited,
            &[0.0; 2],
            || false
        )
        .is_err()
    );
    let mut polls = 0;
    p.electromagnetic_network_candidate(
        &binding,
        &t,
        0.00005,
        &c,
        &[0.0],
        drive,
        controls(),
        &[0.0; 2],
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut count = 0;
    assert!(matches!(
        p.electromagnetic_network_candidate(
            &binding,
            &t,
            0.00005,
            &c,
            &[0.0],
            drive,
            controls(),
            &[0.0; 2],
            || {
                count += 1;
                count == polls
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(p, old);
    assert_eq!(c, old_c);
    near(t.time_s(), 0.0);
}

#[test]
fn material_inertia_and_source_receipts_survive_direct_electromagnetic_network_evolution() {
    use fs_ga::Vec3;
    use fs_mbd_material::{RegionBinding, VolumeMoments, resolve_body};
    let (a, ca) = card(0.2, 2.0, 0.0);
    let (b, cb) = card(0.4, 8.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let resistors = resistor();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let t = initial(&curves);
    let c = circuit(&binding, &t);
    let numeric = motion_fixture::piston(0.0);
    let body = resolve_body(
        &[RegionBinding {
            key: "moving wall",
            material: &cards[0],
            point: &context,
            geometry: VolumeMoments::cuboid(Vec3::new(0.1, 0.1, 0.025)).unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap();
    let p = MaterialDoubleActingGasPiston::new(
        numeric.primary().gas().clone(),
        numeric.opposing().clone(),
        body.clone(),
        [0.01; 2],
        0.0,
    )
    .unwrap();
    let proposal = p
        .electromagnetic_network_candidate(
            &binding,
            &t,
            0.00005,
            &c,
            &[0.0],
            drive(0.0),
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let direct = p
        .piston()
        .electromagnetic_network_candidate(
            &binding,
            &t,
            0.00005,
            &c,
            &[0.0],
            drive(0.0),
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    assert_eq!(proposal.piston().body(), &body);
    assert_eq!(proposal.piston().piston(), direct.piston());
    assert_eq!(proposal.circuit(), direct.network().circuit());
    near(
        proposal.piston().piston().primary().mass_kg(),
        body.inertia().mass(),
    );
    near(proposal.report().coupled.total_energy_residual_j, 0.0);
    let mut polls = 0;
    p.electromagnetic_network_candidate(
        &binding,
        &t,
        0.00005,
        &c,
        &[0.0],
        drive(0.0),
        controls(),
        &[0.0; 2],
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut count = 0;
    assert!(matches!(
        p.electromagnetic_network_candidate(
            &binding,
            &t,
            0.00005,
            &c,
            &[0.0],
            drive(0.0),
            controls(),
            &[0.0; 2],
            || {
                count += 1;
                count == polls
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(p.body(), &body);
    assert_eq!(p.piston().primary().velocity_m_s(), 0.0);
}
