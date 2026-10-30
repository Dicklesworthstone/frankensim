//! G1/G3/G4: material heat transport changes the same thermally powered actuator.
//! Coefficients are manufactured; these checks do not qualify a measured alloy.
use super::*;
use fs_couple::material_circuit::motion::MaterialFieldMotion;
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_ga::Vec3;
use fs_material::field_storage::{MAGNETIC_PERMEABILITY_DIMS,
    MAGNETIC_PERMEABILITY_PROPERTY, ResolvedInductor};
use fs_mbd_material::{RegionBinding, VolumeMoments, resolve_body};

#[path = "../../electrical_contacts/source.rs"]
mod field_source;

fn actuator(circuit: &MaterialCircuit, linkage: f64, velocity: f64) -> MaterialFieldMotion {
    // Replace the fixture's zero-voltage external load branch with a physical
    // uniform-gap inductor. Its separate load resistor and contact are retained.
    let mut graph = circuit.graph().clone();
    graph.branches[4].2 = Branch::Inductor { henries: f64::NAN };
    let field = field_source::scalar_card(MAGNETIC_PERMEABILITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS, 1.0);
    let storage = MaterialStorage { branch: 4,
        element: MaterialStorageElement::Inductor(ResolvedInductor::resolve(&field,
            &QueryPoint::new(), MaterialPropertySelection::SingleClaimOnly, 1.0, 1.0, 1).unwrap()) };
    let circuit = MaterialCircuit::new_with_storage(graph, circuit.resistors().to_vec(),
        vec![storage], vec![linkage], vec![], &[0.0; 2]).unwrap()
        .rebind_contacts_candidate(circuit.contacts().to_vec(), &[0.0; 2]).unwrap();
    let density = field_source::scalar_card("density", Density::DIMS, 2.0);
    let body = resolve_body(&[RegionBinding { key: "guided pole", material: &density,
        point: &QueryPoint::new(), selection: MaterialPropertySelection::SingleClaimOnly,
        geometry: VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap() }], 1, || false).unwrap();
    MaterialFieldMotion::new(circuit, 4, body, velocity, [0.5, 1.5]).unwrap()
}

fn stored(model: &MaterialFieldMotion) -> f64 {
    // Independent uniform-gap energy: mu=N=A=1, guided mass=2 kg.
    0.5 * model.circuit().fluxes_wb()[0].powi(2) * model.gap_m()
        + model.velocity_m_s().powi(2)
}

#[test]
fn g1_material_heat_paths_power_the_original_actuator_and_close_every_cell() {
    let run = |emissivity, resistance, layers| with_heat(emissivity, resistance,
        layers, None, 4, |cx, _, prepared, initial, circuit| {
        let mut moving = actuator(circuit, 0.0, 0.0);
        let initial_energy = initial.energy_j();
        let mut history = None;
        let mut last_current = 0.0;
        for index in 1..=4 {
            let accepted = history.as_ref().unwrap_or(initial);
            let before = [0, 1, 2, 3].map(|id| accepted.temperature_k(CellId(id)).unwrap());
            let target = f64::from(index) * 0.02;
            let dt = target - accepted.time_s();
            let next = prepared.step_thermoelectric_motion_candidate(cx, accepted, target,
                &moving, &[0.0; 2], 0.0, &legs(), config(), &[0.0; 4]).unwrap();
            let report = next.report();
            let [a, b, c, d] = [0, 1, 2, 3].map(|id| next.thermal().temperature_k(CellId(id)).unwrap());
            let current = report.coupled.source.legs[0].transfer.current_a;
            near(report.coupled.source.legs[1].transfer.current_a, current);
            // Independent KVL with actual linkage change and the 1+0.5 ohm load.
            near((next.circuit().fluxes_wb()[0] - moving.circuit().fluxes_wb()[0]) / dt,
                6.5 * current + 0.03 * (b - a));
            near(report.coupled.deposited_joule_j, dt * 6.5 * current * current);
            near(report.coupled.electrical.contacts[0].heat_j, dt * 0.5 * current * current);
            near(report.coupled.source.reversible_electrical_work_j, dt * 0.03 * (a - b) * current);
            near(report.coupled.source.reversible_heat_j,
                -report.coupled.source.reversible_electrical_work_j);

            // Independent finite-volume balances: Fourier resistance or exact
            // integral of affine k(T), gray parallel-plate radiation, Kelvin
            // junction heat, bulk Joule loss and the authored 1/4:3/4 partition.
            let ka = 0.3 + 0.001 * a;
            let kb = 0.3 + 0.001 * b;
            let ga = 1.0 / (0.5 / ka + 0.5);
            let gb = 1.0 / (0.5 / kb + 0.5);
            let conduction = if layers {
                assert_eq!(resistance, 0.0);
                (0.3 + 0.0005 * (a + b)) * (a - b)
            } else {
                (a - b) / (0.5 / ka + resistance + 0.5 / kb)
            };
            let radiation = 5.670_374_419e-8 * (a.powi(4) - b.powi(4))
                / (2.0 / emissivity - 1.0);
            let losses = current * current;
            let power = [
                ga * (c + d - 2.0 * a) - conduction - radiation - 0.03 * a * current + 1.125 * losses,
                gb * (c + d - 2.0 * b) + conduction + radiation + 0.03 * b * current + 0.375 * losses,
                ga * (a - c) + gb * (b - c) + 2.0 * losses,
                ga * (a - d) + gb * (b - d) + 3.0 * losses,
            ];
            for ((old, new), power) in before.into_iter().zip([a, b, c, d]).zip(power) {
                near(10.0 * (new - old), dt * power);
            }
            near(next.thermal().energy_j() + stored(next.motion()), initial_energy);
            near(report.external_electrical_work_j, 0.0);
            near(report.external_mechanical_work_j, 0.0);
            near(report.external_heat_j, 0.0);
            near(report.closed_radiation_net_j, 0.0);
            near(report.total_energy_residual_j, 0.0);
            assert_eq!(report.coupled.source.network.thermal_contacts.len(), usize::from(!layers));
            assert_eq!(report.coupled.source.network.integrated_conduction.len(), usize::from(layers));
            assert!(report.coupled.source.network.radiation.is_some());
            assert!(current > 0.0);
            let (thermal, model, _) = next.into_parts();
            history = Some(thermal); moving = model; last_current = current;
        }
        assert!(moving.gap_m() < 1.0 && moving.velocity_m_s() < 0.0);
        (last_current, moving.gap_m(), moving.circuit().fluxes_wb()[0])
    });
    let baseline = run(0.2, 0.5, false);
    for changed in [run(0.8, 0.5, false), run(0.2, 3.0, false), run(0.2, 0.0, true)] {
        assert!((changed.0 - baseline.0).abs() > 1e-7);
        assert!((changed.1 - baseline.1).abs() > 1e-12);
        assert!((changed.2 - baseline.2).abs() > 1e-8);
    }
}

#[test]
fn g3_empty_heat_assembly_matches_original_motion_with_actual_external_work() {
    let (cards, curves) = fixture(0.02, -0.01, 3.0, 1.0);
    let context = point(350.0); let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key); let rs = resistors(); let links = links();
    let bulk = MaterialNetwork { cells: &cells, resistors: &rs, links: &links };
    let thermal = initial(&curves, [400.0, 300.0, 350.0, 350.0]);
    let fixed = circuit(&bulk, &thermal, false, 0.0);
    let moving = actuator(&fixed, 0.4, 0.02);
    let electrical = ElectricalContactNetwork { bulk: &bulk, contacts: &[], electrical_bases: &BASES };
    let heat = MaterialHeatNetwork { electrical: &electrical, thermal_contacts: &[], radiation: None };
    scoped(|cx, _| {
        let prepared = heat.prepare(cx, &thermal).unwrap();
        let original = thermoelectric::motion::step_candidate(&bulk, &thermal, 0.1, &moving,
            &[0.0; 2], 0.2, &legs(), &BASES, config(), &[0.0; 4], || false).unwrap();
        let actual = prepared.step_thermoelectric_motion_candidate(cx, &thermal, 0.1, &moving,
            &[0.0; 2], 0.2, &legs(), config(), &[0.0; 4]).unwrap();
        assert_eq!(actual.thermal().cells(), original.thermal().cells());
        assert_eq!(actual.motion(), original.motion());
        near(actual.report().external_mechanical_work_j, 0.2 * (actual.motion().gap_m() - moving.gap_m()));
        near(actual.thermal().energy_j() - thermal.energy_j() + stored(actual.motion()) - stored(&moving),
            actual.report().external_mechanical_work_j);
    });
}

#[test]
fn g4_unresolved_kelvin_edge_and_cancellation_preserve_all_histories() {
    with_heat(0.2, 0.5, false, None, 0, |cx, _, prepared, thermal, circuit| {
        let moving = actuator(circuit, 0.4, 0.0);
        let saved = thermal.cells().clone(); let saved_motion = moving.clone();
        assert!(matches!(prepared.step_thermoelectric_motion_candidate(cx, thermal, 0.01,
            &moving, &[0.0; 2], 0.0, &legs(), config(), &[0.0; 4]),
            Err(CircuitEnthalpyError::Invalid("a Kelvin edge requires its original junction-to-bulk temperature trace"))));
        assert_eq!(thermal.cells(), &saved); assert_eq!(moving, saved_motion);
    });
    with_heat(0.2, 0.0, true, None, 4, |cx, gate, prepared, thermal, circuit| {
        let moving = actuator(circuit, 0.4, 0.0);
        let first = prepared.step_thermoelectric_motion_candidate(cx, thermal, 0.01,
            &moving, &[0.0; 2], 0.0, &legs(), config(), &[0.0; 4]).unwrap();
        let replay = prepared.step_thermoelectric_motion_candidate(cx, thermal, 0.01,
            &moving, &[0.0; 2], 0.0, &legs(), config(), &[0.0; 4]).unwrap();
        assert_eq!(first.motion(), replay.motion());
        assert_eq!(first.thermal().cells(), replay.thermal().cells());
        let saved = thermal.cells().clone(); let saved_motion = moving.clone();
        gate.request();
        assert!(matches!(prepared.step_thermoelectric_motion_candidate(cx, thermal, 0.01,
            &moving, &[0.0; 2], 0.0, &legs(), config(), &[0.0; 4]),
            Err(CircuitEnthalpyError::Cancelled)));
        assert_eq!(thermal.cells(), &saved); assert_eq!(moving, saved_motion);
    });
}
