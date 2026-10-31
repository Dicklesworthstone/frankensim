//! Natural and forced patches share one finite-cell/circuit/chamber transaction.
use super::{controls, exchange, input, near, piston, scope, source, surface, DT};
use std::collections::BTreeMap;
use fs_airflow::graph::thermal::medium::convection::finite_gas::{
    NaturalGasSurface, network::exchange::{self as coupled, ExchangeContact, GasHeatLaw, GasHeatSample},
};
use fs_conduction::transient::deposition::*;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::{Chamber, Contact, NetworkSolve};
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::{
    exchange::MaterialExchangeNetwork, state::{ChamberDrive, MaterialExchangeState},
};
use fs_couple::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
use fs_couple::circuit_enthalpy::{CircuitEnthalpySolve, contacts::material::ElectricalContactNetwork,
    material::{MaterialCell, MaterialNetwork}};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve};
use fs_matdb::{ClaimSet, MaterialCard, MaterialStateId, QueryPoint};
use fs_matdb::{InterpolationPolicy, PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel};
use fs_evidence::ValidityDomain;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::{GasNozzle, network::snapshot::SnapshotNozzleLink},
    species::{network::{SpeciesNetworkSolve, transport::{self as shared, GasTransportSolve}},
        diffusion::{BinaryDiffusionMaterial, GasDiffusionGeometry, BINARY_DIFFUSIVITY_DIMS,
            BINARY_DIFFUSIVITY_PROPERTY, network::{DiffusionFaceSolve, DiffusionNetworkLink}}},
};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::network::GasNetworkSolve;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::Length;
use fs_qty::{QuantitySpec, chemistry::SpeciesId, semantic::{QuantityKind, SemanticType, ValueForm}};

const HOT: CellId = CellId(10);
const COLD: CellId = CellId(20);
const MASS: f64 = 0.001;

fn fixture<R>(run: impl FnOnce(&MaterialExchangeNetwork<'_>, DepositionSimulation<'_>, MaterialCircuit) -> R) -> R {
    // The same deliberately synthetic card and caloric chart as the other
    // finite-body tests. Neither cell claims an unsupplied conductivity.
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic body".into(),
        phase: "equilibrium phase".into(), process: "manufactured".into(), revision: 0 },
        ClaimSet::new(), vec![]).unwrap();
    let chart = source::body_curve();
    let context = QueryPoint::new().with("T", 400.0).unwrap();
    let bindings = [HOT, COLD].map(|cell| MaterialCell { cell, material_id: MaterialId(0),
        curve: &chart, point: ThermalMaterialPoint { material: &card, context: &context,
            temperature_axis: "T", liquid_fraction_axis: None }, conductivity: None });
    let bulk = MaterialNetwork { cells: &bindings, resistors: &[], links: &[] };
    let electrical = ElectricalContactNetwork { bulk: &bulk, contacts: &[], electrical_bases: &[] };
    let assembly = MaterialExchangeNetwork { electrical: &electrical, thermal_contacts: &[], radiation: None };
    let events = [(HOT, 260000.0), (COLD, 50000.0)].map(|(id, enthalpy)| TimedEvent {
        time_s: 0.0, event: Event::Birth(DepositedCell { id, segment_id: id.0, body_id: id.0,
            material_id: MaterialId(0), mass_kg: MASS, specific_h_j_kg: enthalpy,
            geometry: CellGeometry { start_m: [0.0; 3], end_m: [0.1, 0.0, 0.0],
                volume_m3: 1e-6, surface_area_m2: 0.02 },
            boundary: Boundary { ambient_k: 400.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 },
        })
    }).to_vec();
    let thermal = DepositionSimulation::new(events, &[(MaterialId(0), &chart as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 2, max_links: 0, max_events: 2, max_substeps: 16 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: DT, energy_tolerance_j: 1e-7 }, || Ok(())).unwrap().0;
    // Passive electrical storage remains part of the same accepted state.
    let circuit = MaterialCircuit::new(CircuitGraph { node_count: 2,
        branches: vec![(1, 0, Branch::Capacitor { farads: 0.1 })], transformers: vec![] },
        vec![], vec![], vec![0.02], &[]).unwrap();
    run(&assembly, thermal, circuit)
}

fn patches() -> [ExchangeContact; 3] {
    let natural = |area_m2| NaturalGasSurface {
        geometry: fs_convection::NaturalConvectionGeometry::VerticalPlate,
        area_m2, characteristic_length: Length::new(0.1), gravity_m_s2: 9.81 };
    [
        ExchangeContact { contact: Contact { id: 1, chamber: 7, cell: HOT }, law: GasHeatLaw::NaturalConvection(natural(0.01)) },
        ExchangeContact { contact: Contact { id: 2, chamber: 7, cell: HOT }, law: GasHeatLaw::ForcedConvection(
            super::ForcedGasSurface { area_m2: 0.01, ..surface() }) },
        ExchangeContact { contact: Contact { id: 3, chamber: 7, cell: COLD }, law: GasHeatLaw::NaturalConvection(natural(0.02)) },
    ]
}

fn network_controls() -> NetworkSolve {
    NetworkSolve { network: CircuitEnthalpySolve {
        maximum_cells: 2, maximum_circuit_trials: 400, maximum_newton_iterations: 24,
        difference_step_j: 1e-5, cell_energy_tolerance_j: 1e-7,
        thermal_energy_tolerance_j: 1e-7, total_energy_tolerance_j: 1e-7,
        circuit: MaterialCircuitSolve { energy_tolerance_j: 1e-7, power_tolerance_w: 1e-7 },
        newton: fs_solver::NewtonKrylovConfig { absolute_tolerance: 1e-10, relative_tolerance: 1e-12,
            // The heat residual is admitted at 1e-10 J. Differencing it across
            // 1e-5 J gives an O(2e-5) relative Jacobian noise floor. Do not ask
            // FGMRES for analytic-Jacobian accuracy after the first Newton step;
            // the outer residual and independent physical gates remain strict.
            forcing_minimum: 1e-4, linear_restart: 4, ..fs_solver::NewtonKrylovConfig::default() },
    }, maximum_chambers: 1, maximum_contacts: 3, maximum_boundary_solves: 1000, total_energy_tolerance_j: 1e-6 }
}

#[test]
fn mixed_patches_share_actual_sources_and_conserve_finite_stores_under_material_swap() {
    let run = |fraction| fixture(|assembly, thermal, circuit| {
        let mut result = None;
        scope(|cx, _| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let gas = piston(fraction, 400.0, 100000.0);
            let chambers = [Chamber { id: 7, piston: &gas, input: input(), piston_solve: controls(), exchange: exchange() }];
            let next = coupled::step_candidate(&prepared, cx, &thermal, DT, &circuit, &[], &chambers,
                &[], &patches(), network_controls(), &[0.0; 2]).unwrap();
            let accepted = &next.chambers()[0];
            let mut heat = BTreeMap::<CellId, f64>::new();
            let mut natural_count = 0;
            let mut forced_count = 0;
            for patch in &accepted.surfaces {
                let (contact, sample) = &patch.source;
                let sample = match sample {
                    GasHeatSample::NaturalConvection(sample) => { natural_count += 1; sample },
                    GasHeatSample::ForcedConvection(sample) => { forced_count += 1; sample },
                    _ => panic!("unexpected heat law"),
                };
                assert_eq!(sample.gas, accepted.boundary.piston.midpoint);
                near(sample.body.temperature_k(), next.thermal().temperature_k(contact.cell).unwrap(), 1e-12);
                let h = sample.evaluation.as_ref().unwrap().coefficient.value.value();
                let q = DT * h * sample.area_m2 * (sample.body.temperature_k() - sample.gas.state().temperature);
                *heat.entry(contact.cell).or_default() += q;
            }
            assert_eq!((natural_count, forced_count), (2, 1));
            assert!(heat[&HOT] > 0.0 && heat[&COLD] < 0.0);
            for cell in [HOT, COLD] {
                near(MASS * (next.thermal().cells()[&cell].specific_h_j_kg
                    - thermal.cells()[&cell].specific_h_j_kg), -heat[&cell], 1e-7);
            }
            let du = gas.gas().mass_kg() * (accepted.piston.gas().state().specific_internal_energy_j_kg()
                - gas.gas().state().specific_internal_energy_j_kg());
            let kinetic = 0.5 * gas.mass_kg() * accepted.piston.velocity_m_s().powi(2);
            let exterior = input().exterior_pressure_pa * (accepted.piston.gas().volume_m3() - gas.gas().volume_m3());
            near(next.thermal().energy_j() - thermal.energy_j() + du + kinetic + exterior, 0.0, 1e-6);
            assert_eq!(next.circuit().charges_c(), circuit.charges_c());
            assert_eq!(thermal.time_s(), 0.0);
            result = Some([next.thermal().temperature_k(HOT).unwrap(), next.thermal().temperature_k(COLD).unwrap()]);
        });
        result.unwrap()
    });
    let a = run(0.75);
    let b = run(0.25);
    assert!((a[0] - b[0]).abs() > 1e-7);
    assert!((a[1] - b[1]).abs() > 1e-7);
}

#[test]
fn live_flow_path_preserves_patch_ownership_and_publishes_all_histories_together() {
    fixture(|assembly, thermal, circuit| scope(|cx, gate| {
        let prepared = assembly.prepare(cx, &thermal).unwrap();
        let gas = piston(0.75, 400.0, 100000.0);
        let mut state = MaterialExchangeState::from_parts(thermal, circuit, BTreeMap::from([(7, gas)]));
        let saved_cells = state.thermal().cells().clone();
        let saved_circuit = state.circuit().clone();
        let saved_gases = state.chambers().clone();
        let drives = [ChamberDrive { chamber: 7, heat_j: 0.0, exterior_pressure_pa: 100000.0,
            external_force_n: 0.0, piston_solve: controls(), exchange: exchange() }];
        let flow = GasNetworkSolve { maximum_parcels: 1, maximum_links: 0,
            maximum_withdrawn_fraction: 0.1, mass_tolerance_kg: 1e-12, energy_tolerance_j: 1e-7 };
        let mut overlap = patches();
        overlap[0].law = GasHeatLaw::NaturalConvection(NaturalGasSurface {
            geometry: fs_convection::NaturalConvectionGeometry::VerticalPlate,
            area_m2: 0.02, characteristic_length: Length::new(0.1), gravity_m_s2: 9.81,
        });
        assert!(coupled::advance_with_flow(&mut state, &prepared, cx, DT, &[], &drives, &[], &overlap,
            &[], flow, network_controls(), &[0.0; 2]).is_err());
        assert_eq!(state.thermal().cells(), &saved_cells);
        assert_eq!(state.circuit(), &saved_circuit);
        assert_eq!(state.chambers(), &saved_gases);
        let (preview, expected) = coupled::flow_step_candidate(&state, &prepared, cx, DT, &[], &drives, &[],
            &patches(), &[], flow, network_controls(), &[0.0; 2]).unwrap();
        assert_eq!(state.time_s(), 0.0);
        let actual = coupled::advance_with_flow(&mut state, &prepared, cx, DT, &[], &drives, &[], &patches(),
            &[], flow, network_controls(), &[0.0; 2]).unwrap();
        assert_eq!(state.time_s(), DT);
        assert_eq!(state.thermal().cells(), preview.thermal().cells());
        assert_eq!(state.circuit(), preview.circuit());
        assert_eq!(state.chambers(), preview.chambers());
        assert_eq!(actual.total_energy_residual_j, expected.total_energy_residual_j);
        assert!(actual.total_energy_residual_j.abs() <= 1e-6);
        let accepted_cells = state.thermal().cells().clone();
        let accepted_gases = state.chambers().clone();
        let accepted_circuit = state.circuit().clone();
        gate.request();
        assert!(coupled::advance_with_flow(&mut state, &prepared, cx, 2.0*DT, &[], &drives, &[], &patches(),
            &[], flow, network_controls(), &[0.0; 2]).is_err());
        assert_eq!(state.time_s(), DT);
        assert_eq!(state.thermal().cells(), &accepted_cells);
        assert_eq!(state.circuit(), &accepted_circuit);
        assert_eq!(state.chambers(), &accepted_gases);
    }));
}

fn diffusion_card(value: f64) -> MaterialCard {
    let typed = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS),
        value: PropertyValue::Scalar { value, dims: BINARY_DIFFUSIVITY_DIMS },
        validity: ValidityDomain::unconstrained()
            .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), 200.0, 1000.0)
            .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e7),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![], provenance: Provenance { source: "manufactured binary pair; not empirical".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic A/B".into(), phase: "gas".into(),
        process: "test".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn diffusion_face<'a>(pairs: &'a [BinaryDiffusionMaterial<'a>]) -> DiffusionNetworkLink<'a> {
    DiffusionNetworkLink { id: 27, parcels: [7, 8], pairs,
        geometry: GasDiffusionGeometry { area_m2: 0.1, distance_m: 0.0001 },
        solve: DiffusionFaceSolve { maximum_temperature_difference_k: 5.0,
            maximum_pressure_difference_pa: 2000.0, gradient_tolerance_per_m: 1e-8,
            molar_flux_tolerance_mol_m2_s: 1e-10 } }
}
fn transport_controls() -> GasTransportSolve {
    GasTransportSolve { exchange: SpeciesNetworkSolve { maximum_species_fraction: 0.05,
        mass_tolerance_kg: 1e-13, energy_tolerance_j: 1e-7 }, entropy_tolerance_j_k: 1e-10 }
}
fn transport_chambers() -> BTreeMap<u64, super::GasPiston> {
    BTreeMap::from([(7, piston(0.75, 400.0, 100000.0)),
        (8, piston(0.25, 400.0, 100000.0)), (9, piston(0.75, 400.0, 60000.0))])
}
fn transport_drives(state: &MaterialExchangeState<'_>) -> Vec<ChamberDrive> {
    state.chambers().iter().map(|(&id, gas)| ChamberDrive { chamber: id, heat_j: 0.0,
        exterior_pressure_pa: gas.gas().state().state().pressure, external_force_n: 0.0,
        piston_solve: controls(), exchange: exchange() }).collect()
}
fn nozzle() -> SnapshotNozzleLink {
    SnapshotNozzleLink { id: 901, parcels: [7, 9], nozzle: GasNozzle {
        area_m2: 2e-6, discharge_coefficient: 0.8 }, nozzle_solve: controls().gas }
}

#[test]
fn combined_transport_composition_supplies_every_actual_convection_trial() {
    fixture(|assembly, thermal, circuit| scope(|cx, _| {
        let card = diffusion_card(2e-5);
        let context = QueryPoint::new();
        let pairs = [BinaryDiffusionMaterial { species: [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
            material: &card, context: &context, selection: MaterialPropertySelection::SingleClaimOnly }];
        let faces = [diffusion_face(&pairs)];
        let prepared = assembly.prepare(cx, &thermal).unwrap();
        let state = MaterialExchangeState::from_parts(thermal, circuit, transport_chambers());
        let drives = transport_drives(&state);
        let coupled_controls = NetworkSolve { maximum_chambers: 3, ..network_controls() };
        let parcels = state.chambers().iter().map(|(&id, g)| (id, g.gas().clone())).collect();
        let endpoints = drives.iter().map(|d| (d.chamber, d.piston_solve.gas)).collect();
        let (moved, flow) = shared::step_candidate(&parcels, &[nozzle()], &faces, DT,
            &endpoints, transport_controls(), || false).unwrap().into_parts();
        let staged = MaterialExchangeState::from_parts(state.thermal().clone(), state.circuit().clone(),
            moved.iter().map(|(&id, gas)| (id, super::GasPiston::new(gas.clone(), 1.0, 0.01, 0.0).unwrap())).collect());
        let empty_flow = GasNetworkSolve { maximum_parcels: 3, maximum_links: 0,
            maximum_withdrawn_fraction: 0.05, mass_tolerance_kg: 1e-13, energy_tolerance_j: 1e-7 };
        let (expected, _) = coupled::flow_step_candidate(&staged, &prepared, cx, DT, &[], &drives,
            &[], &patches(), &[], empty_flow, coupled_controls, &[0.0; 2]).unwrap();
        let (actual, report) = coupled::transport_step_candidate(&state, &prepared, cx, DT, &[], &drives,
            &[], &patches(), &[nozzle()], &faces, transport_controls(), coupled_controls, &[0.0; 2]).unwrap();
        assert_eq!(report.flow, flow);
        assert!(report.flow.nozzles[0].flow.mass_flow_kg_s > 0.0);
        assert!(report.flow.exchange.parcels.iter().any(|p| p.incoming_mass_kg > 0.0));
        assert_ne!(moved[&7].model(), state.chambers()[&7].gas().model());
        assert_eq!(actual.thermal().cells(), expected.thermal().cells());
        assert_eq!(actual.circuit(), expected.circuit());
        for (&id, gas) in actual.chambers() {
            assert_eq!(gas.gas(), expected.chambers()[&id].gas());
            assert_eq!(gas.velocity_m_s(), expected.chambers()[&id].velocity_m_s());
        }
        let mut wall_heat = BTreeMap::<CellId, f64>::new();
        for chamber in &report.exchange.report.source.chambers {
            for patch in &chamber.surfaces {
                let (contact, sample) = &patch.source;
                let sample = match sample {
                    GasHeatSample::NaturalConvection(sample) | GasHeatSample::ForcedConvection(sample) => sample,
                    _ => panic!("expected convective patch"),
                };
                assert_eq!(sample.gas, chamber.boundary.piston.midpoint);
                assert_eq!(sample.gas.thermodynamics().receipt().mole_fraction_bits(),
                    moved[&chamber.id].state().thermodynamics().receipt().mole_fraction_bits());
                near(sample.body.temperature_k(), actual.thermal().temperature_k(contact.cell).unwrap(), 1e-12);
                let h = sample.evaluation.as_ref().unwrap().coefficient.value.value();
                *wall_heat.entry(contact.cell).or_default() += DT * h * sample.area_m2
                    * (sample.body.temperature_k() - sample.gas.state().temperature);
            }
        }
        for cell in [HOT, COLD] {
            near(MASS * (actual.thermal().cells()[&cell].specific_h_j_kg
                - state.thermal().cells()[&cell].specific_h_j_kg), -wall_heat[&cell], 1e-7);
        }
        // Independent physical endpoint accounting, including the outside
        // pressure work; carried species enthalpy and wall heat cancel inside.
        let mut balance = actual.thermal().energy_j() - state.thermal().energy_j()
            + (actual.circuit().charges_c()[0].powi(2) - state.circuit().charges_c()[0].powi(2)) / (2.0 * 0.1);
        for (&id, before) in state.chambers() {
            let after = &actual.chambers()[&id];
            balance += after.gas().mass_kg() * after.gas().state().specific_internal_energy_j_kg()
                - before.gas().mass_kg() * before.gas().state().specific_internal_energy_j_kg()
                + 0.5 * before.mass_kg() * (after.velocity_m_s().powi(2) - before.velocity_m_s().powi(2))
                + before.gas().state().state().pressure * (after.gas().volume_m3() - before.gas().volume_m3());
        }
        near(balance, 0.0, 1e-6);
        assert!(report.total_energy_residual_j.abs() <= coupled_controls.total_energy_tolerance_j);
        assert!(report.total_mass_residual_kg.abs() <= transport_controls().exchange.mass_tolerance_kg);
        assert_eq!(state.time_s(), 0.0);
    }));
}

#[test]
fn diffusion_data_changes_mixed_convection_and_motion_without_replacing_physics() {
    let run = |value| fixture(|assembly, thermal, circuit| {
        let mut observed = None;
        scope(|cx, _| {
            let card = diffusion_card(value); let context = QueryPoint::new();
            let pairs = [BinaryDiffusionMaterial { species: [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
                material: &card, context: &context, selection: MaterialPropertySelection::SingleClaimOnly }];
            let faces = [diffusion_face(&pairs)];
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut state = MaterialExchangeState::from_parts(thermal, circuit, transport_chambers());
            let drives = transport_drives(&state);
            let controls = NetworkSolve { maximum_chambers: 3, ..network_controls() };
            for _ in 0..8 {
                let target = state.time_s() + DT;
                let report = coupled::advance_with_transport(&mut state, &prepared, cx, target, &[], &drives,
                    &[], &patches(), &[nozzle()], &faces, transport_controls(), controls, &[0.0; 2]).unwrap();
                assert!(report.total_energy_residual_j.abs() <= controls.total_energy_tolerance_j);
            }
            observed = Some([state.thermal().temperature_k(HOT).unwrap(),
                state.chambers()[&7].position_m(), state.chambers()[&7].gas().state().state().sound_speed]);
        });
        observed.unwrap()
    });
    let a = run(2e-5); let b = run(4e-5);
    assert_ne!(a[0], b[0]); assert_ne!(a[1], b[1]); assert!((a[2] - b[2]).abs() > 0.01);
}

#[test]
fn combined_transport_late_heat_refusal_preserves_every_accepted_owner() {
    fixture(|assembly, thermal, circuit| scope(|cx, gate| {
        let card = diffusion_card(2e-5); let context = QueryPoint::new();
        let pairs = [BinaryDiffusionMaterial { species: [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
            material: &card, context: &context, selection: MaterialPropertySelection::SingleClaimOnly }];
        let faces = [diffusion_face(&pairs)];
        let prepared = assembly.prepare(cx, &thermal).unwrap();
        let mut state = MaterialExchangeState::from_parts(thermal, circuit, transport_chambers());
        let drives = transport_drives(&state);
        let controls = NetworkSolve { maximum_chambers: 3, ..network_controls() };
        let saved_cells = state.thermal().cells().clone(); let saved_gases = state.chambers().clone();
        let saved_circuit = state.circuit().clone();
        let mut invalid = patches();
        invalid[1].law = GasHeatLaw::ForcedConvection(super::ForcedGasSurface { speed_m_s: 1e20, ..surface() });
        assert!(coupled::advance_with_transport(&mut state, &prepared, cx, DT, &[], &drives,
            &[], &invalid, &[nozzle()], &faces, transport_controls(), controls, &[0.0; 2]).is_err());
        assert_eq!(state.time_s(), 0.0); assert_eq!(state.thermal().cells(), &saved_cells);
        assert_eq!(state.chambers(), &saved_gases); assert_eq!(state.circuit(), &saved_circuit);
        gate.request();
        assert!(matches!(coupled::advance_with_transport(&mut state, &prepared, cx, DT, &[], &drives,
            &[], &patches(), &[nozzle()], &faces, transport_controls(), controls, &[0.0; 2]),
            Err(fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::state::flow::MaterialGasLawFlowError::Cancelled)));
        assert_eq!(state.time_s(), 0.0); assert_eq!(state.thermal().cells(), &saved_cells);
        assert_eq!(state.chambers(), &saved_gases); assert_eq!(state.circuit(), &saved_circuit);
    }));
}
