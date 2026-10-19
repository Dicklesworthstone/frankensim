//! Synthetic source data through the real material, circuit, thermal and gas paths.
use super::*;
use crate::cavity_phs::piston::tests::{gas, controls as motion_controls, drive as motion_drive};
use crate::cavity_phs::piston::thermal::boundary::network::surfaces::{
    NetworkThermalPort, material::{MaterialGasSurface, MaterialSurfacesDrive},
};
use crate::circuit_enthalpy::material::{MaterialCell, CellResistor, ConductionLink};
use crate::circuit_enthalpy::ResistorCell;
use crate::material_circuit::MaterialCircuitSolve;
use fs_evidence::ValidityDomain;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_matdb::*;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{QuantitySpec, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

const BASES: [(usize, ElectricalPropertyBasis); 2] = [
    (2, ElectricalPropertyBasis::Resistivity), (1, ElectricalPropertyBasis::Conductivity),
];
fn id(i: usize) -> CellId { CellId(if i == 0 { 30 } else { 10 }) }
struct Fixture {
    cards: [MaterialCard; 2],
    curves: [EquilibriumEnthalpyPhaseCurve; 2],
    point: QueryPoint,
    key: PropertyKey,
}
impl Fixture {
    fn new(sigma: f64, split_property: Option<&str>) -> Self {
        let cards = [0, 1].map(|i| {
            let basis = if i == 0 { ElectricalPropertyBasis::Conductivity } else { ElectricalPropertyBasis::Resistivity };
            let mut claims = ClaimSet::new();
            for (name, dims, value) in [
                (basis.property_name(), basis.dimensions(), if i == 0 { sigma } else { 0.4 }),
                ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, if i == 0 { 0.2 } else { 0.4 }),
            ] {
                let intervals = if i == 0 && split_property == Some(name) {
                    vec![(200.0, 430.0), (431.0, 800.0)]
                } else { vec![(200.0, 800.0)] };
                for (lo, hi) in intervals {
                    claims.insert_claim(PropertyClaim {
                        key: PropertyKey::new(name, dims),
                        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                            knots: vec![(lo, value), (hi, 2.0 * value)], dims },
                        validity: ValidityDomain::unconstrained()
                            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), lo, hi).with("lot", 7.0, 7.0),
                        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
                        observations: vec![], provenance: Provenance { source: "synthetic snapshot wiring fixture".into(),
                            license: "CC0-1.0".into(), artifact: None },
                    }).unwrap();
                }
            }
            MaterialCard::assemble(MaterialStateId { chemistry: format!("synthetic-{i}"), phase: "solid".into(),
                process: "test".into(), revision: 0 }, claims, vec![]).unwrap()
        });
        let curves = [0, 1].map(|i| EquilibriumEnthalpyPhaseCurve::try_single_phase(cards[i].content_hash(),
            SolidLiquidPhase::Solid, vec![
                EnthalpyPhaseKnot { specific_enthalpy_j_kg: 0.0, temperature_k: 200.0, liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1.0 },
                EnthalpyPhaseKnot { specific_enthalpy_j_kg: 600000.0, temperature_k: 800.0, liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1.0 },
            ]).unwrap());
        Self { cards, curves, point: QueryPoint::new().with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 400.0)
            .unwrap().with("lot", 7.0).unwrap(), key: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS) }
    }
    fn cells(&self) -> [MaterialCell<'_>; 2] {
        [0, 1].map(|i| MaterialCell { cell: id(i), material_id: MaterialId(i as u64), curve: &self.curves[i],
            point: ThermalMaterialPoint { material: &self.cards[i], context: &self.point,
                temperature_axis: "T", liquid_fraction_axis: None },
            conductivity: Some((&self.key, MaterialPropertySelection::SingleClaimOnly)) })
    }
    fn state(&self, temperatures: [f64; 2]) -> DepositionSimulation<'_> {
        let mut events = [0, 1].map(|i| TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id: id(i), material_id: MaterialId(i as u64), segment_id: i as u64, body_id: i as u64,
            mass_kg: 0.001, specific_h_j_kg: 1000.0 * (temperatures[i] - 200.0),
            geometry: CellGeometry { start_m: [i as f64, 0.0, 0.0], end_m: [i as f64 + 0.1, 0.0, 0.0],
                volume_m3: 0.001, surface_area_m2: 0.04 },
            boundary: Boundary { ambient_k: 400.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 },
        }) }).to_vec();
        events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
            id: 7, a: id(0), b: id(1), area_m2: 0.002, conductance_w_k: 999.0,
        }) });
        DepositionSimulation::new(events, &[(MaterialId(0), &self.curves[0] as &dyn ThermalLaw),
            (MaterialId(1), &self.curves[1] as &dyn ThermalLaw)],
            DepositionBudget { max_cells: 2, max_links: 1, max_events: 3, max_substeps: 100 }).unwrap()
            .advance_to(0.0, StepConfig { max_dt_s: 0.001, energy_tolerance_j: 1e-7 }, || Ok(())).unwrap().0
    }
}
fn resistors() -> [CellResistor; 2] {
    // Reverse both branch and cell ordering to exercise the existing identity mapping.
    [1, 0].map(|i| CellResistor { destination: ResistorCell { branch: i + 1, cell: id(i) },
        selection: MaterialPropertySelection::SingleClaimOnly, length_m: 0.1, area_m2: 0.01 })
}
fn links() -> [ConductionLink; 1] { [ConductionLink { link: 7, distance_a_m: 0.001, distance_b_m: 0.001 }] }
fn controls() -> NetworkGasSolve {
    NetworkGasSolve { network: CircuitEnthalpySolve {
        maximum_cells: 2, maximum_circuit_trials: 300, maximum_newton_iterations: 32,
        difference_step_j: 1e-4, cell_energy_tolerance_j: 1e-7, thermal_energy_tolerance_j: 1e-7,
        total_energy_tolerance_j: 1e-7, circuit: MaterialCircuitSolve { energy_tolerance_j: 1e-7, power_tolerance_w: 1e-7 },
        newton: NewtonKrylovConfig { absolute_tolerance: 1e-8, relative_tolerance: 1e-12,
            forcing_minimum: 1e-4, linear_restart: 4, ..NewtonKrylovConfig::default() },
    }, maximum_boundary_solves: 300, total_energy_tolerance_j: 1e-6 }
}
fn drive(dt: f64) -> NetworkGasDrive {
    NetworkGasDrive { cell: id(0), input: PistonInput { duration_s: dt, ..motion_drive(0.0) }, piston: motion_controls(),
        exchange: HeatExchangeSolve { heat_bracket_j: [-1.0, 1.0], source_tolerance_j: 1e-9,
            energy_tolerance_j: 1e-7, maximum_iterations: 128 } }
}
fn geometry() -> GasConductionGeometry { GasConductionGeometry { area_m2: 0.01, gas_distance_m: 0.001, cell_distance_m: 0.02 } }
fn circuit(binding: &MaterialNetwork<'_>, state: &DepositionSimulation<'_>) -> MaterialCircuit {
    let initial = binding.resolve_with_bases(state, &BASES, || false).unwrap();
    MaterialCircuit::new(CircuitGraph { node_count: 3, branches: vec![
        (1, 0, Branch::Capacitor { farads: 0.01 }), (1, 2, Branch::Resistor { ohms: 0.0 }),
        (2, 0, Branch::Resistor { ohms: 0.0 }),
    ], transformers: vec![] }, initial.resistors, vec![], vec![1.0], &[]).unwrap()
}

#[test]
fn prepared_sources_refresh_values_and_refuse_switching_either_transport_claim() {
    for split in ["electrical_conductivity", "thermal_conductivity"] {
        let fixture = Fixture::new(5.0, Some(split));
        let cells = fixture.cells(); let resistors = resistors(); let links = links();
        let binding = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
        let cold = fixture.state([400.0, 400.0]); let warm = fixture.state([420.0, 400.0]);
        let hot = fixture.state([500.0, 400.0]);
        let prepared = binding.prepare_with_bases(&cold, &BASES, || false).unwrap();
        let first = prepared.resolve(&cold, || false).unwrap();
        assert_eq!(first.source, binding.resolve_with_bases(&cold, &BASES, || false).unwrap().source);
        let second = prepared.resolve(&warm, || false).unwrap();
        let property = |sample: &MaterialNetworkSample| sample.cells.iter().find(|c| c.cell == id(0))
            .unwrap().conductivity.as_ref().unwrap().property(split).unwrap().clone();
        assert_ne!(property(&first.source).value_si(), property(&second.source).value_si());
        assert_eq!(property(&first.source).answer().receipt.selected, property(&second.source).answer().receipt.selected);
        assert!(binding.resolve_with_bases(&hot, &BASES, || false).is_ok());
        assert!(prepared.resolve(&hot, || false).is_err(), "a new in-domain claim is not the selected model");
        assert!(matches!(prepared.resolve(&warm, || true), Err(MaterialNetworkError::Thermal(DepositionError::Cancelled))));
        assert_eq!(prepared.resolve(&cold, || false).unwrap().source, first.source);
    }
}

#[test]
fn preparation_cancellation_and_incomplete_basis_do_not_publish_or_mutate_history() {
    let fixture = Fixture::new(5.0, None); let cells = fixture.cells(); let resistors = resistors(); let links = links();
    let binding = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let state = fixture.state([400.0, 400.0]); let saved = state.cells().clone();
    assert!(binding.prepare(&state, || false).is_err(), "legacy basis remains resistivity-only");
    assert!(binding.prepare_with_bases(&state, &BASES[..1], || false).is_err());
    assert!(binding.prepare_with_bases(&state, &[BASES[0]; 2], || false).is_err());
    let mut polls = 0;
    binding.prepare_with_bases(&state, &BASES, || { polls += 1; false }).unwrap();
    let mut seen = 0;
    assert!(matches!(binding.prepare_with_bases(&state, &BASES, || { seen += 1; seen == polls }),
        Err(MaterialNetworkError::Thermal(DepositionError::Cancelled))));
    assert_eq!(state.cells(), &saved);
    assert!(binding.prepare_with_bases(&state, &BASES, || false).is_ok());
}

#[test]
fn conductivity_only_electrical_data_reaches_single_and_multiple_surface_motion() {
    let run = |sigma, multiple| {
        let fixture = Fixture::new(sigma, None); let cells = fixture.cells(); let resistors = resistors(); let links = links();
        let binding = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
        let mut thermal = fixture.state([450.0, 420.0]); let mut electrical = circuit(&binding, &thermal);
        let mut piston = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap();
        assert!(fixture.cards[0].claims_for("electrical_resistivity").is_empty());
        for n in 1..=40 {
            let target = f64::from(n) * 0.00005; let drive = drive(target - thermal.time_s());
            let (next, circuit, gas, report, sources, residual) = if multiple {
                let surfaces = [MaterialGasSurface { port: NetworkThermalPort { id: 9, cell: id(0) }, geometry: geometry() }];
                let candidate = piston.material_network_surfaces_candidate_with_bases(&binding, &thermal, target, &electrical,
                    &[], &BASES, MaterialSurfacesDrive { surfaces: &surfaces, input: drive.input, piston: drive.piston,
                        exchange: drive.exchange }, controls(), &[0.0; 2], || false).unwrap();
                let (t, c, g, r) = candidate.into_parts();
                (t, c, g, r.network, r.ports.into_iter().map(|p| p.source).collect::<Vec<_>>(), r.coupled.total_energy_residual_j)
            } else {
                let candidate = piston.material_network_candidate_with_bases(&binding, &thermal, target, &electrical, &[],
                    &BASES, drive, controls(), geometry(), &[0.0; 2], || false).unwrap();
                let (t, c, g, r) = candidate.into_parts();
                (t, c, g, r.network, vec![r.interface], r.coupled.total_energy_residual_j)
            };
            let body = report.source.cells.iter().find(|c| c.cell == id(0)).unwrap();
            let conductor = &report.electrical.resistors.iter().find(|r| r.material.branch == 1).unwrap().material.conductor;
            assert_eq!(conductor.source_basis(), ElectricalPropertyBasis::Conductivity);
            assert_eq!(body.conductivity.as_ref(), Some(&sources[0]));
            assert_eq!(conductor.material(), &sources[0]);
            assert_eq!(sources[0].properties().len(), 2);
            assert!(residual.abs() <= controls().total_energy_tolerance_j);
            thermal = next; electrical = circuit; piston = gas;
        }
        (electrical.charges_c()[0], piston.position_m(), piston.gas().state().state().sound_speed)
    };
    for multiple in [false, true] {
        let a = run(5.0, multiple); let b = run(10.0, multiple);
        assert!((a.0 - b.0).abs() > 1e-7);
        assert!((a.1 - b.1).abs() > 1e-12);
        assert!((a.2 - b.2).abs() > 1e-8);
    }
}
