//! Synthetic source cards through the original flow, circuit, heat and motion owners.
use super::*;
use crate::cavity_phs::piston::material::MaterialGasPiston;
use crate::cavity_phs::piston::tests::{controls as motion_controls, drive as motion_drive, gas};
use crate::cavity_phs::piston::thermal::HeatExchangeSolve;
use crate::cavity_phs::piston::thermal::boundary::network::material::{GasConductionGeometry, resolve_contact};
use crate::circuit_enthalpy::{CircuitEnthalpySolve, ResistorCell};
use crate::circuit_enthalpy::material::{CellResistor, ConductionLink, MaterialCell, MaterialNetwork};
use crate::material_circuit::MaterialCircuitSolve;
use fs_evidence::ValidityDomain;
use fs_ga::Vec3;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferSolve;
use fs_matdb::*;
use fs_mbd_material::{RegionBinding, VolumeMoments, resolve_body};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_ports::{ModuleWindowAudit, PortedModule};
use fs_qty::{Density, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

const BASES: [(usize, ElectricalPropertyBasis); 2] = [
    (2, ElectricalPropertyBasis::Resistivity), (1, ElectricalPropertyBasis::Conductivity),
];
fn id(i: usize) -> CellId { CellId(if i == 0 { 30 } else { 10 }) }
fn point(t: f64) -> QueryPoint { QueryPoint::new().with("T", t).unwrap().with("lot", 7.0).unwrap() }
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
            let basis = if i == 0 { ElectricalPropertyBasis::Conductivity } else { ElectricalPropertyBasis::Resistivity };
            let mut claims = ClaimSet::new();
            for (key, a, b) in [
                (PropertyKey::new(basis.property_name(), basis.dimensions()),
                    if i == 0 { sigma } else { 0.01 }, if i == 0 { 0.5 * sigma } else { 0.02 }),
                (conductivity.clone(), if i == 0 { k } else { 0.4 }, if i == 0 { 2.0 * k } else { 0.8 }),
                (cp.clone(), 500.0, 500.0), (density.clone(), 4000.0, 4000.0),
            ] {
                let dims = key.quantity().dims();
                claims.insert_claim(PropertyClaim { key,
                    value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                        knots: vec![(200.0, a), (800.0, b)], dims },
                    validity: ValidityDomain::unconstrained().with("T", 200.0, 800.0).with("lot", 7.0, 7.0),
                    uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
                    observations: vec![], provenance: Provenance { source: "synthetic open gas/network fixture".into(),
                        license: "CC0-1.0".into(), artifact: None } }).unwrap();
            }
            MaterialCard::assemble(MaterialStateId { chemistry: format!("synthetic-{i}"), phase: "solid".into(),
                process: "test".into(), revision: 0 }, claims, vec![]).unwrap()
        });
        let caloric = [0, 1].map(|i| ResolvedHeatCapacityCurve::try_new(&cards[i], HeatCapacitySourceQuery {
            heat_capacity: &cp, density: &density, temperature_axis: "T", lower: &point(200.0), upper: &point(800.0),
            selection: SelectionPolicy::SingleClaimOnly,
        }, SolidLiquidPhase::Solid, 0.0, 1e-5).unwrap());
        Self { cards, caloric, context: point(400.0), conductivity }
    }
    fn bindings(&self) -> [MaterialCell<'_>; 2] {
        [0, 1].map(|i| MaterialCell { cell: id(i), material_id: MaterialId(i as u64), curve: self.caloric[i].curve(),
            point: ThermalMaterialPoint { material: &self.cards[i], context: &self.context,
                temperature_axis: "T", liquid_fraction_axis: None },
            conductivity: Some((&self.conductivity, MaterialPropertySelection::SingleClaimOnly)) })
    }
    fn state(&self) -> DepositionSimulation<'_> {
        let mut events = [0, 1].map(|i| TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id: id(i), material_id: MaterialId(i as u64), segment_id: i as u64, body_id: i as u64,
            mass_kg: 0.001, specific_h_j_kg: if i == 0 { 125000.0 } else { 110000.0 },
            geometry: CellGeometry { start_m: [i as f64, 0.0, 0.0], end_m: [i as f64 + 0.01, 0.0, 0.0],
                volume_m3: 0.001 / 4000.0, surface_area_m2: 0.04 },
            boundary: Boundary { ambient_k: 400.0, htc_w_m2_k: 0.0, reservoir: None,
                source_w: if i == 0 { 0.0 } else { 2.0 } },
        }) }).to_vec();
        events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
            id: 7, a: id(0), b: id(1), area_m2: 0.002, conductance_w_k: 999.0,
        }) });
        DepositionSimulation::new(events, &[(MaterialId(0), self.caloric[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), self.caloric[1].curve() as &dyn ThermalLaw)],
            DepositionBudget { max_cells: 2, max_links: 1, max_events: 3, max_substeps: 100 }).unwrap()
            .advance_to(0.0, StepConfig { max_dt_s: 0.001, energy_tolerance_j: 1e-7 }, || Ok(())).unwrap().0
    }
    fn piston(&self) -> MaterialGasPiston {
        let body = resolve_body(&[RegionBinding { key: "piston", material: &self.cards[0], point: &self.context,
            geometry: VolumeMoments::cuboid(Vec3::new(0.1, 0.1, 0.025)).unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly }], 1, || false).unwrap();
        MaterialGasPiston::new(gas(0.75, 0.001), body, 0.01, 0.0).unwrap()
    }
}
fn resistors() -> [CellResistor; 2] {
    [1, 0].map(|i| CellResistor { destination: ResistorCell { branch: i + 1, cell: id(i) },
        selection: MaterialPropertySelection::SingleClaimOnly, length_m: 0.01, area_m2: 0.000025 })
}
fn links() -> [ConductionLink; 1] { [ConductionLink { link: 7, distance_a_m: 0.001, distance_b_m: 0.001 }] }
fn circuit(binding: &MaterialNetwork<'_>, state: &DepositionSimulation<'_>) -> MaterialCircuit {
    MaterialCircuit::new(CircuitGraph { node_count: 3, branches: vec![
        (1, 0, Branch::Capacitor { farads: 0.01 }), (1, 2, Branch::Resistor { ohms: 0.0 }),
        (2, 0, Branch::Resistor { ohms: 0.0 }),
    ], transformers: vec![] }, binding.resolve_with_bases(state, &BASES, || false).unwrap().resistors,
        vec![], vec![1.0], &[]).unwrap()
}
fn solve() -> NetworkGasSolve {
    NetworkGasSolve { network: CircuitEnthalpySolve {
        maximum_cells: 2, maximum_circuit_trials: 300, maximum_newton_iterations: 32,
        difference_step_j: 1e-4, cell_energy_tolerance_j: 1e-7, thermal_energy_tolerance_j: 1e-7,
        total_energy_tolerance_j: 1e-7, circuit: MaterialCircuitSolve { energy_tolerance_j: 1e-7, power_tolerance_w: 1e-7 },
        newton: NewtonKrylovConfig { absolute_tolerance: 1e-8, relative_tolerance: 1e-12,
            forcing_minimum: 1e-4, linear_restart: 4, ..NewtonKrylovConfig::default() },
    }, maximum_boundary_solves: 300, total_energy_tolerance_j: 1e-6 }
}
fn drive(dt: f64) -> NetworkGasDrive {
    NetworkGasDrive { cell: id(0), input: PistonInput { duration_s: dt, heat_j: 3.0 * dt,
        external_force_n: 0.2, ..motion_drive(0.0) }, piston: motion_controls(),
        exchange: HeatExchangeSolve { heat_bracket_j: [-1.0, 1.0], source_tolerance_j: 1e-9,
            energy_tolerance_j: 1e-7, maximum_iterations: 128 } }
}
fn flow() -> GasNozzleStepSolve {
    GasNozzleStepSolve { nozzle: motion_controls().gas, maximum_donor_fraction: 0.1,
        transfer: GasMassTransferSolve { donor: motion_controls().gas, receiver: motion_controls().gas,
            mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 } }
}
fn nozzle() -> GasNozzle { GasNozzle { area_m2: 1e-5, discharge_coefficient: 0.8 } }
fn reservoir(x: f64) -> ClosedGasParcel {
    let mut source = gas(x, 0.001);
    source.add_heat(100.0, motion_controls().gas, || false).unwrap(); source
}
fn geometry() -> GasConductionGeometry {
    GasConductionGeometry { area_m2: 0.01, gas_distance_m: 0.001, cell_distance_m: 0.02 }
}

#[test]
fn data_only_mixture_conductivity_and_resistance_changes_reach_every_owner() {
    let run = |k, sigma, fraction| {
        let fixture = Fixture::new(k, sigma); let cells = fixture.bindings();
        let resistors = resistors(); let links = links();
        let binding = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
        let mut thermal = fixture.state(); let mut circuit = circuit(&binding, &thermal);
        let mut piston = fixture.piston(); let mut reservoir = reservoir(fraction);
        let initial_body = piston.body().clone(); let initial_gas = piston.gas().clone();
        let initial_reservoir = reservoir.clone(); let initial_h = thermal.energy_j();
        let initial_mass = piston.gas().mass_kg() + reservoir.mass_kg();
        let mut supplied = 0.0; let mut wall_heat = 0.0;
        for n in 1..=32 {
            let target = f64::from(n) * 0.00005; let drive = drive(target - thermal.time_s());
            let candidate = piston.nozzle_material_network_candidate(&reservoir, nozzle(), &binding,
                &thermal, target, &circuit, &[], &BASES, drive, flow(), solve(), geometry(), &[0.0; 2], || false).unwrap();
            let (next_t, next_c, next_p, next_r, report) = candidate.into_parts();
            assert_eq!(next_p.body(), &initial_body);
            assert_eq!(next_p.piston().mass_kg(), initial_body.inertia().mass());
            assert_eq!(report.direction, Some(NozzleDirection::IntoPiston));
            assert!(report.total_mass_residual_kg.abs() <= flow().transfer.mass_tolerance_kg);
            assert!(report.total_energy_residual_j.abs() <= solve().total_energy_tolerance_j);
            let stage = &report.network;
            let original = stage.network.source.cells.iter().find(|c| c.cell == id(0)).unwrap();
            assert_eq!(original.conductivity.as_ref(), Some(&stage.interface));
            let branch = stage.network.electrical.resistors.iter().find(|r| r.material.branch == 1).unwrap();
            assert_eq!(branch.material.conductor.material(), &stage.interface);
            assert_eq!(branch.material.conductor.source_basis(), ElectricalPropertyBasis::Conductivity);
            assert!(stage.interface.query_point().contains(&("lot".into(), 7.0)));
            assert!(stage.network.deposited_joule_j > 0.0);
            assert!((report.reservoir_energy_change_j + report.transfer.transferred_enthalpy_j).abs() < 1e-8);
            let motion = &stage.boundary.piston;
            supplied += stage.coupled.external_network_heat_j + drive.input.heat_j
                + stage.network.electrical.native.supplied + motion.external_force_work_j - motion.exterior_pressure_work_j;
            wall_heat += stage.boundary.exchanged_heat_j;
            thermal = next_t; circuit = next_c; piston = next_p; reservoir = next_r;
        }
        // Independently reconstruct all finite stores; no sum of step residuals.
        let gas_change = super::super::stored_change(&initial_gas, piston.gas())
            + super::super::stored_change(&initial_reservoir, &reservoir);
        let electrical_change = (circuit.charges_c()[0].powi(2) - 1.0) / (2.0 * 0.01);
        let mechanical_change = 0.5 * piston.piston().mass_kg() * piston.piston().velocity_m_s().powi(2);
        let total = electrical_change + thermal.energy_j() - initial_h + gas_change + mechanical_change - supplied;
        assert!(total.abs() < 32.0 * solve().total_energy_tolerance_j);
        assert!((piston.gas().mass_kg() + reservoir.mass_kg() - initial_mass).abs() < 1e-13);
        assert!(matches!(piston.piston().audit_last_window().unwrap(), ModuleWindowAudit::Unavailable { .. }));
        [circuit.charges_c()[0], thermal.temperature_k(id(0)).unwrap(), wall_heat,
            piston.piston().velocity_m_s(), piston.gas().state().state().sound_speed]
    };
    let reference = run(0.2, 200.0, 0.01);
    for changed in [run(2.0, 200.0, 0.01), run(0.2, 400.0, 0.01), run(0.2, 200.0, 0.9)] {
        assert!((changed[2] - reference[2]).abs() > 1e-10, "wall heat must react to shared state");
        assert!((changed[3] - reference[3]).abs() > 1e-12, "actual momentum must react");
        assert!((changed[4] - reference[4]).abs() > 1e-9, "accepted acoustic material state must react");
    }
}

#[test]
fn shut_path_is_exactly_the_existing_closed_network_including_history() {
    let fixture = Fixture::new(0.2, 200.0); let cells = fixture.bindings();
    let resistors = resistors(); let links = links();
    let binding = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let thermal = fixture.state(); let circuit = circuit(&binding, &thermal);
    let piston = fixture.piston(); let source = reservoir(0.25); let drive = drive(0.00005);
    let expected = piston.piston().material_network_candidate_with_bases(&binding, &thermal,
        drive.input.duration_s, &circuit, &[], &BASES, drive, solve(), geometry(), &[0.0; 2], || false).unwrap();
    let actual = piston.nozzle_material_network_candidate(&source, GasNozzle { area_m2: 0.0, ..nozzle() },
        &binding, &thermal, drive.input.duration_s, &circuit, &[], &BASES, drive, flow(), solve(), geometry(),
        &[0.0; 2], || false).unwrap();
    let (t, c, p, r, report) = actual.into_parts();
    let (et, ec, ep, er) = expected.into_parts();
    assert_eq!(p.piston(), &ep); assert_eq!(r, source);
    assert_eq!(t.cells(), et.cells()); assert_eq!(c.charges_c(), ec.charges_c());
    assert_eq!(report.network.boundary, er.boundary); assert_eq!(report.network.coupled, er.coupled);
    assert_eq!(report.network.interface, er.interface);
    assert_eq!(report.direction, None); assert_eq!(report.transfer.transferred_mass_kg, 0.0);
}

#[test]
fn multiple_network_trials_withdraw_the_reservoir_only_once_and_reverse_flow_works() {
    let fixture = Fixture::new(0.2, 200.0); let cells = fixture.bindings();
    let resistors = resistors(); let links = links();
    let binding = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let thermal = fixture.state(); let circuit = circuit(&binding, &thermal); let drive = drive(0.00005);
    let destinations = resistors.iter().map(|r| r.destination).collect::<Vec<_>>();
    let prepared = binding.prepare_with_bases(&thermal, &BASES, || false).unwrap();
    for outward in [false, true] {
        let (chamber, source) = if outward { (reservoir(0.75), gas(0.25, 0.001)) }
            else { (gas(0.75, 0.001), reservoir(0.25)) };
        let piston = GasPiston::new(chamber, 1.0, 0.01, 0.0).unwrap();
        let (transfer, expected_flow) = if outward {
            piston.gas().through_nozzle_to(&source, nozzle(), drive.input.duration_s, flow(), || false).unwrap()
        } else {
            source.through_nozzle_to(piston.gas(), nozzle(), drive.input.duration_s, flow(), || false).unwrap()
        };
        let (donor, receiver, expected_transfer) = transfer.into_parts();
        let candidate = piston.nozzle_network_candidate(&source, nozzle(), &thermal, drive.input.duration_s,
            CircuitEnthalpyInput { circuit: &circuit, ports: &[], destinations: &destinations }, drive, flow(), solve(),
            &[0.0; 2], |state| prepared.resolve(state, || false),
            |state, source, gas| resolve_contact(&binding, id(0), state, source, gas, geometry()), || false).unwrap();
        assert_eq!(candidate.reservoir(), if outward { &receiver } else { &donor });
        assert_eq!(candidate.report().transfer, expected_transfer);
        assert_eq!(candidate.report().flow, expected_flow);
        assert!(candidate.report().network.coupled.boundary_solves > 1);
        assert_eq!(candidate.report().direction, Some(if outward { NozzleDirection::OutOfPiston } else { NozzleDirection::IntoPiston }));
        assert!(candidate.report().total_energy_residual_j.abs() <= solve().total_energy_tolerance_j);
    }
}

#[test]
fn late_failure_or_cancellation_preserves_all_four_accepted_owners() {
    let fixture = Fixture::new(0.2, 200.0); let cells = fixture.bindings();
    let resistors = resistors(); let links = links();
    let binding = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let thermal = fixture.state(); let circuit = circuit(&binding, &thermal);
    let piston = fixture.piston(); let source = reservoir(0.25); let drive = drive(0.00005);
    let saved_piston = piston.clone(); let saved_source = source.clone();
    let saved_cells = thermal.cells().clone(); let saved_charges = circuit.charges_c().to_vec();
    let run = |geometry: GasConductionGeometry, target: f64, cancelled: &mut dyn FnMut() -> bool|
        piston.nozzle_material_network_candidate(&source, nozzle(),
        &binding, &thermal, target, &circuit, &[], &BASES, drive, flow(), solve(), geometry,
        &[0.0; 2], cancelled);
    let mut polls = 0;
    run(geometry(), drive.input.duration_s, &mut || { polls += 1; false }).unwrap();
    assert!(polls > 10);
    let mut seen = 0;
    assert!(matches!(run(geometry(), drive.input.duration_s,
        &mut || { seen += 1; seen == polls }), Err(CircuitEnthalpyError::Cancelled)));
    // The geometry refusal occurs after the valid nozzle transfer, inside the
    // actual material wall evaluation. Its tentative species cannot escape.
    assert!(run(GasConductionGeometry { area_m2: 1.0, ..geometry() }, drive.input.duration_s,
        &mut || false).is_err());
    assert!(matches!(run(geometry(), 2.0 * drive.input.duration_s,
        &mut || false), Err(CircuitEnthalpyError::Invalid(_))));
    assert_eq!(piston, saved_piston); assert_eq!(source, saved_source);
    assert_eq!(thermal.cells(), &saved_cells); assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(circuit.charges_c(), saved_charges); assert!(circuit.fluxes_wb().is_empty());
}
