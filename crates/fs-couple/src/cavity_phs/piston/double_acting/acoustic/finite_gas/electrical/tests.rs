//! Original circuit, motion, finite regional gas and complete acoustic owners.
//! Synthetic source laws exercise conservative composition, not measurements.

use super::super::{
    FiniteGasAcousticCandidate, FiniteGasAcousticError,
    tests::{
        Inputs, Loads, Parcels, assert_loads_same, faces, fixture, nozzles, original_faces,
        transport_controls,
    },
};
use crate::cavity_phs::piston::PistonError;
use crate::cavity_phs::piston::double_acting::acoustic::tests::{
    assert_piston_same, caloric_energy, motion, near,
};
use crate::cavity_phs::piston::double_acting::electrical::tests::{
    circuit, electrical, input, mechanical, piston, solve,
};
use crate::cavity_phs::piston::double_acting::electrical::{
    ElectricalMotionError, ElectricalPistonCandidate,
};
use crate::cavity_phs::piston::double_acting::{DoubleActingGasPiston, DoublePistonInput};
use crate::dg_port::DgPortLoad;
use crate::material_circuit::MaterialCircuit;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::transport::GasTransportPaths;

type Candidate = FiniteGasAcousticCandidate<ElectricalPistonCandidate>;
type Error = FiniteGasAcousticError<ElectricalMotionError<PistonError>>;

fn run(
    piston: &DoubleActingGasPiston,
    circuit: &MaterialCircuit,
    parcels: &Parcels,
    loads: &Loads,
    inputs: &Inputs,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Candidate, Error> {
    let links = nozzles();
    piston.transported_electroacoustic_candidate(
        &faces(loads),
        parcels,
        inputs,
        GasTransportPaths {
            nozzles: &links,
            ..Default::default()
        },
        electrical(circuit, &[0.0], 2.0),
        input(),
        mechanical(),
        solve(),
        motion(),
        transport_controls(),
        piston.primary().velocity_m_s(),
        cancelled,
    )
}

fn physical_storage(
    piston: &DoubleActingGasPiston,
    circuit: &MaterialCircuit,
    parcels: &Parcels,
    loads: &Loads,
) -> f64 {
    // The fixture has one 0.01 F capacitor and no inductors. Integrate the
    // synthetic Cp/R polynomials and actual finite masses independently of
    // reported gas, electrical, transport or enclosure residuals.
    caloric_energy(piston.primary().gas())
        + caloric_energy(piston.opposing())
        + 0.5 * piston.primary().mass_kg() * piston.primary().velocity_m_s().powi(2)
        + circuit.charges_c()[0].powi(2) / (2.0 * 0.01)
        + parcels.values().map(caloric_energy).sum::<f64>()
        + loads.values().map(DgPortLoad::stored_energy_j).sum::<f64>()
}

fn assert_same(actual: &Candidate, expected: &Candidate) {
    let a = actual.native();
    let b = expected.native();
    assert_piston_same(a.piston(), b.piston());
    assert_eq!(a.piston_step(), b.piston_step());
    assert_eq!(a.circuit().circuit(), b.circuit().circuit());
    assert_eq!(a.circuit().step().native.x, b.circuit().step().native.x);
    assert_eq!(a.circuit().step().native.y, b.circuit().step().native.y);
    assert_eq!(a.circuit().step().resistors, b.circuit().step().resistors);
    assert_eq!(a.circuit().step().heat_j, b.circuit().step().heat_j);
    assert_eq!(
        [
            a.motion().current_a,
            a.motion().back_emf_v,
            a.motion().mechanical_work_j,
            a.motion().electrical_work_j,
            a.total_energy_residual_j(),
        ],
        [
            b.motion().current_a,
            b.motion().back_emf_v,
            b.motion().mechanical_work_j,
            b.motion().electrical_work_j,
            b.total_energy_residual_j(),
        ]
    );
    assert_eq!(actual.parcels(), expected.parcels());
    assert_loads_same(actual.loads(), expected.loads());
    assert_eq!(actual.step().motion, expected.step().motion);
    assert_eq!(actual.step().face_parcels, expected.step().face_parcels);
    assert_eq!(
        actual.step().transport.parameter_work_j,
        expected.step().transport.parameter_work_j
    );
    assert_eq!(
        actual.step().transport.combined_energy_residual_j,
        expected.step().transport.combined_energy_residual_j
    );
    assert_eq!(
        actual.total_energy_residual_j(),
        expected.total_energy_residual_j()
    );
}

#[test]
fn g1_regional_gas_changes_the_next_current_without_replaying_the_current_interval() {
    let p = piston(0.03);
    let c = circuit(2.0, 0.2, false);
    let (parcels, loads, mut inputs) = fixture(input().duration_s);
    inputs.get_mut(&10).unwrap().supplied_heat_j = 1e-3;
    inputs.get_mut(&20).unwrap().supplied_heat_j = 5e-4;
    let frozen = p
        .electroacoustic_candidate(
            &original_faces(&loads),
            electrical(&c, &[0.0], 2.0),
            input(),
            mechanical(),
            solve(),
            motion(),
            p.primary().velocity_m_s(),
            || false,
        )
        .unwrap();
    let candidate = run(&p, &c, &parcels, &loads, &inputs, &mut || false).unwrap();

    // Gas feedback follows exactly the original native frozen-medium sample.
    assert_piston_same(candidate.native().piston(), frozen.native().piston());
    assert_eq!(
        candidate.native().piston_step(),
        frozen.native().piston_step()
    );
    assert_eq!(
        candidate.native().circuit().circuit(),
        frozen.native().circuit().circuit()
    );
    assert_eq!(
        candidate.native().motion().current_a,
        frozen.native().motion().current_a
    );
    assert_eq!(
        candidate.native().motion().back_emf_v,
        frozen.native().motion().back_emf_v
    );
    for (index, id) in [10, 20].into_iter().enumerate() {
        let updated = &candidate.loads()[&id];
        let baseline = &frozen.loads()[index];
        let actual = &candidate.step().motion.faces[index];
        let original = &frozen.motion().faces[index];
        assert_eq!(updated.state().step_count(), 1);
        assert_eq!(updated.state().step_count(), baseline.state().step_count());
        assert_eq!(actual.inlet_flow_m3_s, original.inlet_flow_m3_s);
        assert_eq!(actual.mean_pressure_pa, original.mean_pressure_pa);
        assert_eq!(actual.inlet_work_j, original.inlet_work_j);
        assert_eq!(
            actual.duct.ledger.port_work_step_j,
            original.duct.ledger.port_work_step_j
        );
        assert_eq!(
            updated.terminal().unwrap().pole_states(),
            baseline.terminal().unwrap().pole_states()
        );
        assert_eq!(
            updated.last_terminal_frames(),
            baseline.last_terminal_frames()
        );
        assert_eq!(
            updated.duct().gas_field().unwrap().uniform_state(),
            Some(candidate.parcels()[&id].state())
        );
        assert_ne!(
            candidate.parcels()[&id]
                .model()
                .thermodynamics()
                .mole_composition(),
            parcels[&id].model().thermodynamics().mole_composition()
        );
    }
    assert!(
        candidate.parcels()[&10].state().state().temperature
            > parcels[&10].state().state().temperature
    );

    // Hold the newly accepted circuit and piston fixed between two comparisons.
    // Only the acoustic endpoint medium differs. Its reaction must reach back
    // through the original velocity/back-emf residual into actual current.
    let frozen_loads: Loads = [
        (10, frozen.loads()[0].clone()),
        (20, frozen.loads()[1].clone()),
    ]
    .into_iter()
    .collect();
    let next_p = candidate.native().piston();
    let next_c = candidate.native().circuit().circuit();
    let refreshed = next_p
        .electroacoustic_candidate(
            &original_faces(candidate.loads()),
            electrical(next_c, &[0.0], 2.0),
            input(),
            mechanical(),
            solve(),
            motion(),
            next_p.primary().velocity_m_s(),
            || false,
        )
        .unwrap();
    let unchanged = next_p
        .electroacoustic_candidate(
            &original_faces(&frozen_loads),
            electrical(next_c, &[0.0], 2.0),
            input(),
            mechanical(),
            solve(),
            motion(),
            next_p.primary().velocity_m_s(),
            || false,
        )
        .unwrap();
    let a = refreshed.native().motion();
    let b = unchanged.native().motion();
    assert!((a.current_a - b.current_a).abs() > 1e-12);
    assert_ne!(a.back_emf_v, b.back_emf_v);
    near(
        a.current_a - b.current_a,
        (b.back_emf_v - a.back_emf_v) / (2.0 + input().duration_s / (2.0 * 0.01)),
        1e-10,
    );
}

#[test]
fn g1_transport_electrical_and_acoustic_storage_close_one_external_energy_window() {
    let mut p = piston(0.03);
    let mut c = circuit(2.0, 0.2, true);
    let (mut parcels, mut loads, inputs) = fixture(input().duration_s);
    let initial_energy = physical_storage(&p, &c, &parcels, &loads);
    let initial_mass: f64 = parcels.values().map(|gas| gas.mass_kg()).sum();
    let initial_loss: f64 = loads
        .values()
        .map(DgPortLoad::cumulative_dissipated_j)
        .sum();
    let mut external = 0.0;
    let mut resistor_heat = 0.0;
    let mut terminal_work = 0.0;
    let links = nozzles();
    for sample in 0..4_u64 {
        let drive = DoublePistonInput {
            heat_j: if sample < 2 { [0.002, 0.001] } else { [0.0; 2] },
            external_force_n: if sample < 2 { 0.3 } else { 0.0 },
            ..input()
        };
        let applied = [0.0, if sample < 2 { 3.0 } else { 0.0 }];
        let candidate = p
            .transported_electroacoustic_candidate(
                &faces(&loads),
                &parcels,
                &inputs,
                GasTransportPaths {
                    nozzles: &links,
                    ..Default::default()
                },
                electrical(&c, &applied, 2.0),
                drive,
                mechanical(),
                solve(),
                motion(),
                transport_controls(),
                p.primary().velocity_m_s(),
                &mut || false,
            )
            .unwrap();
        let native = candidate.native();
        let sound = &candidate.step().motion;
        let dx = native.piston_step().displacement_m;

        // The reserved source is internal. Port 1 is the only imposed voltage.
        external += drive.duration_s * applied[1] * native.circuit().step().native.y[1]
            + drive.heat_j.iter().sum::<f64>()
            + drive.external_force_n * dx
            + inputs
                .values()
                .map(|input| input.supplied_heat_j)
                .sum::<f64>();
        resistor_heat += native.circuit().step().heat_j;
        terminal_work += sound.terminal_work_j;
        near(
            sound.mechanical_work_j
                + sound
                    .faces
                    .iter()
                    .map(|face| face.inlet_work_j)
                    .sum::<f64>(),
            0.0,
            1e-10,
        );
        near(
            candidate
                .parcels()
                .values()
                .map(|gas| gas.mass_kg())
                .sum::<f64>(),
            initial_mass,
            2e-16,
        );
        for load in candidate.loads().values() {
            assert_eq!(load.state().step_count(), sample + 1);
        }
        assert_eq!(native.piston().ports.next_tick() - 1, sample + 1);
        let acoustic_loss: f64 = candidate
            .loads()
            .values()
            .map(DgPortLoad::cumulative_dissipated_j)
            .sum();
        near(
            physical_storage(
                native.piston(),
                native.circuit().circuit(),
                candidate.parcels(),
                candidate.loads(),
            ) - initial_energy
                + acoustic_loss
                - initial_loss
                + resistor_heat
                + terminal_work
                - external,
            0.0,
            (sample + 1) as f64 * motion().energy_tolerance_j,
        );
        let (native, next_parcels, next_loads, _) = candidate.into_parts();
        let (next_c, next_p, _, _, _, _) = native.into_parts();
        p = next_p;
        c = next_c;
        parcels = next_parcels;
        loads = next_loads;
    }
    assert!(loads.values().all(|load| load.stored_energy_j() > 0.0));
}

#[test]
fn g4_late_gas_refusal_and_cancellation_preserve_all_owners_and_the_retry_suffix() {
    let p = piston(0.03);
    let c = circuit(2.0, 0.2, false);
    let (parcels, loads, inputs) = fixture(input().duration_s);
    let saved_p = p.clone();
    let saved_c = c.clone();
    let saved_parcels = parcels.clone();
    let saved_loads = loads.clone();
    let mut refused_inputs = inputs.clone();
    refused_inputs.get_mut(&10).unwrap().supplied_heat_j = 1e6;
    assert!(run(&p, &c, &parcels, &loads, &refused_inputs, &mut || false,).is_err());

    let mut polls = 0;
    let expected = run(&p, &c, &parcels, &loads, &inputs, &mut || {
        polls += 1;
        false
    })
    .unwrap();
    let mut count = 0;
    let cancelled = run(&p, &c, &parcels, &loads, &inputs, &mut || {
        count += 1;
        count == polls
    });
    assert!(matches!(cancelled, Err(FiniteGasAcousticError::Cancelled)));
    assert_eq!(count, polls);
    assert_piston_same(&p, &saved_p);
    assert_eq!(c, saved_c);
    assert_eq!(parcels, saved_parcels);
    assert_loads_same(&loads, &saved_loads);

    let retried = run(&p, &c, &parcels, &loads, &inputs, &mut || false).unwrap();
    assert_same(&retried, &expected);
    let (a, parcels_a, loads_a, _) = retried.into_parts();
    let (b, parcels_b, loads_b, _) = expected.into_parts();
    let (c_a, p_a, _, _, _, _) = a.into_parts();
    let (c_b, p_b, _, _, _, _) = b.into_parts();
    let next_a = run(&p_a, &c_a, &parcels_a, &loads_a, &inputs, &mut || false).unwrap();
    let next_b = run(&p_b, &c_b, &parcels_b, &loads_b, &inputs, &mut || false).unwrap();
    assert_same(&next_a, &next_b);
}

#[test]
fn g1_finite_heat_network_and_transport_retain_one_native_thermal_history() {
    use crate::cavity_phs::piston::double_acting::electrical::network::tests as network_fixture;
    use crate::circuit_enthalpy::material::MaterialNetwork;
    use fs_matdb::PropertyKey;
    use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;
    use fs_thermal::enthalpy::DepositionSimulation;

    fn thermal_energy(thermal: &DepositionSimulation<'_>) -> f64 {
        thermal
            .cells()
            .values()
            .map(|cell| cell.mass_kg * cell.specific_h_j_kg)
            .sum()
    }

    // Dyadic duration preserves the exact common DG/thermal clock.
    let dt = 1.0 / 32768.0;
    let (a, ca) = network_fixture::card(0.2, 2.0, 0.002);
    let (b, cb) = network_fixture::card(0.4, 8.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = network_fixture::point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = network_fixture::cells(&cards, &curves, &context, &key);
    let resistors = network_fixture::resistor();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let mut thermal = network_fixture::initial(&curves);
    let mut c = network_fixture::circuit(&binding, &thermal);
    let mut p = piston(0.02);
    let (mut parcels, mut loads, inputs) = fixture(dt);
    let mut drive = network_fixture::drive(0.01);
    drive.gas.input.duration_s = dt;
    drive.gas.input.heat_j = [0.002, 0.001];
    drive.gas.input.external_force_n = 0.3;
    drive.gas.piston.volume_tolerance_m3 = 1e-18;
    drive.transducer.factor_n_a = 20.0;
    let mut acoustic = motion();
    acoustic.energy_tolerance_j = 1e-6;
    let mut transport = transport_controls();
    transport.energy_tolerance_j = 1e-6;
    let initial_energy = physical_storage(&p, &c, &parcels, &loads) + thermal_energy(&thermal);
    let initial_loss: f64 = loads
        .values()
        .map(DgPortLoad::cumulative_dissipated_j)
        .sum();
    let mut external = 0.0;
    let mut terminal_work = 0.0;
    let links = nozzles();

    for sample in 0..2_u64 {
        let target = thermal.time_s() + dt;
        let frozen = p
            .electroacoustic_network_candidate(
                &original_faces(&loads),
                &binding,
                &thermal,
                target,
                &c,
                &[0.0],
                drive,
                network_fixture::controls(),
                acoustic,
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let candidate = p
            .transported_electroacoustic_network_candidate(
                &faces(&loads),
                &parcels,
                &inputs,
                GasTransportPaths {
                    nozzles: &links,
                    ..Default::default()
                },
                &binding,
                &thermal,
                target,
                &c,
                &[0.0],
                drive,
                network_fixture::controls(),
                acoustic,
                transport,
                &[0.0; 2],
                &mut || false,
            )
            .unwrap();
        let native = candidate.native();
        assert_piston_same(native.piston(), frozen.native().piston());
        assert_eq!(
            native.network().circuit(),
            frozen.native().network().circuit()
        );
        assert_eq!(
            native.network().thermal().cells(),
            frozen.native().network().thermal().cells()
        );
        assert_eq!(
            native.network().report().source.boundary.step(),
            frozen.native().network().report().source.boundary.step()
        );
        assert_eq!(
            native.network().report().electrical.heat_j,
            frozen.native().network().report().electrical.heat_j
        );
        near(
            native
                .network()
                .report()
                .source
                .motion
                .external_electrical_work_j,
            0.0,
            1e-12,
        );
        near(native.coupled().external_network_heat_j, 0.0, 1e-12);
        external += drive.gas.input.heat_j.iter().sum::<f64>()
            + drive.gas.input.external_force_n
                * native
                    .network()
                    .report()
                    .source
                    .boundary
                    .step()
                    .piston
                    .displacement_m
            + inputs
                .values()
                .map(|input| input.supplied_heat_j)
                .sum::<f64>();
        terminal_work += candidate.step().motion.terminal_work_j;
        let loss: f64 = candidate
            .loads()
            .values()
            .map(DgPortLoad::cumulative_dissipated_j)
            .sum();

        // Native resistor heat is already in finite cell enthalpy. It is not
        // exported a second time when the complete finite network is enclosed.
        near(
            physical_storage(
                native.piston(),
                native.network().circuit(),
                candidate.parcels(),
                candidate.loads(),
            ) + thermal_energy(native.network().thermal())
                - initial_energy
                + loss
                - initial_loss
                + terminal_work
                - external,
            0.0,
            (sample + 1) as f64 * 1e-6,
        );
        for load in candidate.loads().values() {
            assert_eq!(load.state().step_count(), sample + 1);
        }
        assert_eq!(native.piston().ports.next_tick() - 1, sample + 1);
        assert_eq!(native.network().thermal().time_s(), target);
        let (native, next_parcels, next_loads, _) = candidate.into_parts();
        let (next_thermal, next_c, next_p, _) = native.into_parts();
        p = next_p;
        c = next_c;
        thermal = next_thermal;
        parcels = next_parcels;
        loads = next_loads;
    }
}
