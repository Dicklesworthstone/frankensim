//! Native circuit, finite cells, two gases, shared motion and original DG loads.
//! Synthetic source cards exercise coupling, not experimental material accuracy.

use super::*;
use crate::cavity_phs::piston::double_acting::acoustic::tests as acoustic_fixture;
use crate::cavity_phs::piston::double_acting::electrical::network::tests as fixture;
use crate::cavity_phs::piston::double_acting::electrical::tests as motion_fixture;
use fs_matdb::PropertyKey;
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;
use fs_thermal::enthalpy::CellId;

// A dyadic interval preserves the exact common DG clock under repeated steps.
const DT: f64 = 1.0 / 32768.0;
const AREAS: [f64; 2] = [0.003, -0.002];

fn drive() -> ElectromagneticNetworkDrive {
    let mut drive = fixture::drive(0.01);
    drive.gas.input.duration_s = DT;
    drive.gas.input.heat_j = [0.002, 0.001];
    drive.gas.input.external_force_n = 0.3;
    drive.gas.piston.volume_tolerance_m3 = 1e-18;
    drive.transducer.factor_n_a = 20.0;
    drive
}

fn acoustic() -> AcousticMotionSolve {
    AcousticMotionSolve {
        energy_tolerance_j: 1e-6,
        ..acoustic_fixture::motion()
    }
}

fn faces(loads: &[DgPortLoad], areas: [f64; 2]) -> [AcousticMotionFace<'_>; 2] {
    [0, 1].map(|i| AcousticMotionFace {
        load: &loads[i],
        signed_area_m2: areas[i],
    })
}

fn thermal_energy(thermal: &DepositionSimulation<'_>) -> f64 {
    thermal
        .cells()
        .values()
        .map(|cell| cell.mass_kg * cell.specific_h_j_kg)
        .sum()
}

fn physical_storage(
    piston: &DoubleActingGasPiston,
    thermal: &DepositionSimulation<'_>,
    circuit: &MaterialCircuit,
    loads: &[DgPortLoad],
) -> f64 {
    // The fixture contains one 0.01 F capacitor and no inductors. Reconstruct
    // that energy, cell enthalpy and gas caloric integrals without reports.
    circuit.charges_c()[0].powi(2) / (2.0 * 0.01)
        + thermal_energy(thermal)
        + acoustic_fixture::caloric_energy(piston.primary().gas())
        + acoustic_fixture::caloric_energy(piston.opposing())
        + 0.5 * piston.primary().mass_kg() * piston.primary().velocity_m_s().powi(2)
        + loads.iter().map(DgPortLoad::stored_energy_j).sum::<f64>()
}

#[test]
fn g1_acoustic_reaction_changes_current_and_finite_heat_with_one_closed_energy_account() {
    let (a, ca) = fixture::card(0.2, 2.0, 0.002);
    let (b, cb) = fixture::card(0.4, 8.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = fixture::point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = fixture::cells(&cards, &curves, &context, &key);
    let resistors = fixture::resistor();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let mut thermal = fixture::initial(&curves);
    let mut circuit = fixture::circuit(&binding, &thermal);
    let mut piston = motion_fixture::piston(0.02);
    let mut loads = vec![
        acoustic_fixture::load(0.15, DT),
        acoustic_fixture::load(0.85, DT),
    ];
    let initial_energy = physical_storage(&piston, &thermal, &circuit, &loads);
    let mut external_work = 0.0;
    let mut terminal_work = 0.0;
    for index in 0..4 {
        let original_loads = loads.clone();
        let heat_before = thermal_energy(&thermal);
        let target = thermal.time_s() + DT;
        let candidate = piston
            .electroacoustic_network_candidate(
                &faces(&loads, AREAS),
                &binding,
                &thermal,
                target,
                &circuit,
                &[0.0],
                drive(),
                fixture::controls(),
                acoustic(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let source = &candidate.native().network().report().source;
        let step = &source.boundary.step().piston;
        let displacement = step.displacement_m;
        acoustic_fixture::near(source.motion.transfer_residual_j, 0.0, 1e-7);
        acoustic_fixture::near(candidate.motion().transfer_residual_j, 0.0, 1e-6);
        acoustic_fixture::near(candidate.total_energy_residual_j(), 0.0, 1e-6);
        acoustic_fixture::near(
            source.motion.back_emf_v,
            drive().transducer.factor_n_a * displacement / DT,
            1e-7,
        );
        acoustic_fixture::near(
            source.motion.force_n,
            drive().transducer.factor_n_a * source.motion.current_a,
            1e-12,
        );
        if index == 0 {
            let unloaded = piston
                .electromagnetic_network_candidate(
                    &binding,
                    &thermal,
                    target,
                    &circuit,
                    &[0.0],
                    drive(),
                    fixture::controls(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            let baseline = &unloaded.network().report().source.motion;
            assert!(source.motion.current_a > baseline.current_a);
            assert!(
                candidate.native().piston().primary().velocity_m_s()
                    < unloaded.piston().primary().velocity_m_s()
            );
            assert!(
                candidate.native().network().report().electrical.heat_j
                    > unloaded.network().report().electrical.heat_j
            );
            assert_ne!(
                candidate
                    .native()
                    .network()
                    .thermal()
                    .temperature_k(CellId(0))
                    .unwrap(),
                unloaded
                    .network()
                    .thermal()
                    .temperature_k(CellId(0))
                    .unwrap()
            );
        }
        external_work += drive().gas.input.heat_j.iter().sum::<f64>()
            + drive().gas.input.external_force_n * displacement;
        terminal_work += candidate.motion().terminal_work_j;
        for i in 0..2 {
            let mut reference = original_loads[i].clone();
            reference
                .preview_step(AREAS[i] * (displacement / DT))
                .unwrap()
                .accept();
            acoustic_fixture::assert_load_same(&candidate.loads()[i], &reference);
            assert_eq!(candidate.loads()[i].state().step_count(), index + 1);
        }
        let (native, next_loads, _, _) = candidate.into_parts();
        let (next_thermal, next_circuit, next_piston, report) = native.into_parts();
        acoustic_fixture::near(
            thermal_energy(&next_thermal) - heat_before,
            report.network.electrical.heat_j - report.boundary.exchanged_heat_j.iter().sum::<f64>(),
            1e-6,
        );
        acoustic_fixture::near(
            report.network.deposited_joule_j,
            report.network.electrical.heat_j,
            1e-12,
        );
        acoustic_fixture::near(report.motion.external_electrical_work_j, 0.0, 1e-12);
        acoustic_fixture::near(report.coupled.external_network_heat_j, 0.0, 1e-12);
        thermal = next_thermal;
        circuit = next_circuit;
        piston = next_piston;
        loads = next_loads;
    }
    let acoustic_loss = loads
        .iter()
        .map(DgPortLoad::cumulative_dissipated_j)
        .sum::<f64>();
    acoustic_fixture::near(
        physical_storage(&piston, &thermal, &circuit, &loads) - initial_energy
            + acoustic_loss
            + terminal_work
            - external_work,
        0.0,
        4e-6,
    );
}

#[test]
fn g3_zero_faces_reduce_to_the_original_finite_network_and_continue_acoustic_memory() {
    let (a, ca) = fixture::card(0.2, 2.0, 0.002);
    let (b, cb) = fixture::card(0.4, 8.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = fixture::point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = fixture::cells(&cards, &curves, &context, &key);
    let resistors = fixture::resistor();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let thermal = fixture::initial(&curves);
    let circuit = fixture::circuit(&binding, &thermal);
    let piston = motion_fixture::piston(0.02);
    let loads = vec![
        acoustic_fixture::load(0.15, DT),
        acoustic_fixture::load(0.85, DT),
    ];
    let prefix = piston
        .electroacoustic_network_candidate(
            &faces(&loads, AREAS),
            &binding,
            &thermal,
            DT,
            &circuit,
            &[0.0],
            drive(),
            fixture::controls(),
            acoustic(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let (native, loads, _, _) = prefix.into_parts();
    let (thermal, circuit, piston, _) = native.into_parts();
    assert!(loads.iter().all(|load| load.stored_energy_j() > 0.0));
    let coupled = piston
        .electroacoustic_network_candidate(
            &faces(&loads, [0.0; 2]),
            &binding,
            &thermal,
            2.0 * DT,
            &circuit,
            &[0.0],
            drive(),
            fixture::controls(),
            acoustic(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let original = piston
        .electromagnetic_network_candidate(
            &binding,
            &thermal,
            2.0 * DT,
            &circuit,
            &[0.0],
            drive(),
            fixture::controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    acoustic_fixture::assert_piston_same(coupled.native().piston(), original.piston());
    assert_eq!(
        coupled.native().network().circuit(),
        original.network().circuit()
    );
    assert_eq!(
        coupled.native().network().thermal().cells(),
        original.network().thermal().cells()
    );
    assert_eq!(
        coupled.native().network().thermal().time_s(),
        original.network().thermal().time_s()
    );
    assert_eq!(
        coupled.native().network().report().source.boundary.step(),
        original.network().report().source.boundary.step()
    );
    assert_eq!(
        coupled.native().network().report().electrical.heat_j,
        original.network().report().electrical.heat_j
    );
    assert_eq!(
        coupled.native().coupled().total_energy_residual_j,
        original.coupled().total_energy_residual_j
    );
    assert_eq!(coupled.motion().force_n, 0.0);
    assert_eq!(coupled.motion().mechanical_work_j, 0.0);
    for (i, load) in loads.iter().enumerate() {
        let mut reference = load.clone();
        reference.preview_step(0.0).unwrap().accept();
        acoustic_fixture::assert_load_same(&coupled.loads()[i], &reference);
        assert_eq!(coupled.loads()[i].state().step_count(), 2);
    }
}

#[test]
fn g4_failed_joint_work_and_final_cancel_preserve_all_owners_and_exact_retry() {
    let (a, ca) = fixture::card(0.2, 2.0, 0.002);
    let (b, cb) = fixture::card(0.4, 8.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = fixture::point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = fixture::cells(&cards, &curves, &context, &key);
    let resistors = fixture::resistor();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let thermal = fixture::initial(&curves);
    let circuit = fixture::circuit(&binding, &thermal);
    let piston = motion_fixture::piston(0.02);
    let loads = vec![
        acoustic_fixture::load(0.15, DT),
        acoustic_fixture::load(0.85, DT),
    ];
    let old_piston = piston.clone();
    let old_circuit = circuit.clone();
    let old_cells = thermal.cells().clone();
    let old_loads = loads.clone();
    let mut polls = 0;
    let reference = piston
        .electroacoustic_network_candidate(
            &faces(&loads, AREAS),
            &binding,
            &thermal,
            DT,
            &circuit,
            &[0.0],
            drive(),
            fixture::controls(),
            acoustic(),
            &[0.0; 2],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut limited = fixture::controls();
    limited.maximum_boundary_solves = 1;
    assert!(
        piston
            .electroacoustic_network_candidate(
                &faces(&loads, AREAS),
                &binding,
                &thermal,
                DT,
                &circuit,
                &[0.0],
                drive(),
                limited,
                acoustic(),
                &[0.0; 2],
                || false,
            )
            .is_err()
    );
    limited = fixture::controls();
    limited.network.maximum_circuit_trials = 1;
    assert!(
        piston
            .electroacoustic_network_candidate(
                &faces(&loads, AREAS),
                &binding,
                &thermal,
                DT,
                &circuit,
                &[0.0],
                drive(),
                limited,
                acoustic(),
                &[0.0; 2],
                || false,
            )
            .is_err()
    );
    let mut missing = fixture::cells(&cards, &curves, &context, &key);
    missing[0].conductivity = None;
    let missing_source = MaterialNetwork {
        cells: &missing,
        resistors: &resistors,
        links: &[],
    };
    assert!(
        piston
            .electroacoustic_network_candidate(
                &faces(&loads, AREAS),
                &missing_source,
                &thermal,
                DT,
                &circuit,
                &[0.0],
                drive(),
                fixture::controls(),
                acoustic(),
                &[0.0; 2],
                || false,
            )
            .is_err()
    );
    let mut count = 0;
    assert!(matches!(
        piston.electroacoustic_network_candidate(
            &faces(&loads, AREAS),
            &binding,
            &thermal,
            DT,
            &circuit,
            &[0.0],
            drive(),
            fixture::controls(),
            acoustic(),
            &[0.0; 2],
            || {
                count += 1;
                count == polls
            },
        ),
        Err(AcousticMotionError::Cancelled)
    ));
    acoustic_fixture::assert_piston_same(&piston, &old_piston);
    assert_eq!(circuit, old_circuit);
    assert_eq!(thermal.cells(), &old_cells);
    assert_eq!(thermal.time_s(), 0.0);
    for (load, old) in loads.iter().zip(&old_loads) {
        acoustic_fixture::assert_load_same(load, old);
    }
    let retry = piston
        .electroacoustic_network_candidate(
            &faces(&loads, AREAS),
            &binding,
            &thermal,
            DT,
            &circuit,
            &[0.0],
            drive(),
            fixture::controls(),
            acoustic(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    assert_eq!(
        retry.native().network().circuit(),
        reference.native().network().circuit()
    );
    assert_eq!(
        retry.native().network().thermal().cells(),
        reference.native().network().thermal().cells()
    );
    acoustic_fixture::assert_piston_same(retry.native().piston(), reference.native().piston());
    assert_eq!(retry.motion(), reference.motion());
    assert_eq!(
        retry.total_energy_residual_j(),
        reference.total_energy_residual_j()
    );
    for (load, expected) in retry.loads().iter().zip(reference.loads()) {
        acoustic_fixture::assert_load_same(load, expected);
    }
}

#[test]
fn g3_material_body_and_original_port_clocks_survive_two_accepted_network_steps() {
    use crate::cavity_phs::piston::material::{MaterialDoubleActingGasPiston, tests::body};

    let (a, ca) = fixture::card(0.2, 2.0, 0.002);
    let (b, cb) = fixture::card(0.4, 8.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = fixture::point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = fixture::cells(&cards, &curves, &context, &key);
    let resistors = fixture::resistor();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let mut thermal = fixture::initial(&curves);
    let mut circuit = fixture::circuit(&binding, &thermal);
    let specimen = body(8000.0);
    let gas_fixture = motion_fixture::piston(0.02);
    let mut piston = MaterialDoubleActingGasPiston::new(
        gas_fixture.primary().gas().clone(),
        gas_fixture.opposing().clone(),
        specimen.clone(),
        [0.01; 2],
        0.02,
    )
    .unwrap();
    let mut loads = vec![
        acoustic_fixture::load(0.15, DT),
        acoustic_fixture::load(0.85, DT),
    ];
    for sample in 0..2 {
        let target = thermal.time_s() + DT;
        let original = piston
            .piston()
            .electroacoustic_network_candidate(
                &faces(&loads, AREAS),
                &binding,
                &thermal,
                target,
                &circuit,
                &[0.0],
                drive(),
                fixture::controls(),
                acoustic(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let sourced = piston
            .electroacoustic_network_candidate(
                &faces(&loads, AREAS),
                &binding,
                &thermal,
                target,
                &circuit,
                &[0.0],
                drive(),
                fixture::controls(),
                acoustic(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        assert_eq!(sourced.native().piston().body(), &specimen);
        acoustic_fixture::assert_piston_same(
            sourced.native().piston().piston(),
            original.native().piston(),
        );
        assert_eq!(
            sourced.native().piston().piston().ports.next_tick(),
            sample + 2
        );
        assert_eq!(
            sourced.native().piston().piston().primary().mass_kg(),
            specimen.inertia().mass()
        );
        assert_eq!(
            sourced.native().circuit(),
            original.native().network().circuit()
        );
        assert_eq!(
            sourced.native().thermal().cells(),
            original.native().network().thermal().cells()
        );
        assert_eq!(sourced.motion(), original.motion());
        assert_eq!(
            sourced.total_energy_residual_j(),
            original.total_energy_residual_j()
        );
        for (load, expected) in sourced.loads().iter().zip(original.loads()) {
            acoustic_fixture::assert_load_same(load, expected);
            assert_eq!(load.state().step_count(), sample + 1);
        }
        let (native, next_loads, _, _) = sourced.into_parts();
        let (next_thermal, next_circuit, next_piston, _) = native.into_parts();
        thermal = next_thermal;
        circuit = next_circuit;
        piston = next_piston;
        loads = next_loads;
    }
    assert_eq!(piston.body(), &specimen);
    assert_eq!(piston.piston().ports.next_tick(), 3);
}
