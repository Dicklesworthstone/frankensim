//! G1/G3/G4: the original elastic pressure port consumes complete real-gas state.
//! Synthetic PR/NASA sources exercise the coupling, not empirical fluid accuracy.
use super::*;
use fs_material::thermochemical_gas::energy::GasEnergyError;
use fs_material::thermochemical_gas::nonideal::energy::parcel::IsochoricRealGasParcel;
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
};
use fs_qty::Temperature;
use fs_thermochem::SpeciesId;
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseV1;

const PRESSURE: f64 = 5.0e6;
const RATE: u32 = 20_000;

fn parcel() -> IsochoricRealGasParcel {
    let source = gas(0.75);
    let model = ThermochemicalPengRobinsonModel::new(
        source.model().clone(),
        vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A").unwrap(),
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B").unwrap(),
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6),
                acentric_factor: -0.00219,
            },
        ],
        vec![0.0, 0.03, 0.03, 0.0],
    )
    .unwrap();
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            300.0,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), PRESSURE)
        .unwrap()
        .with("specimen-lot", 42.0)
        .unwrap();
    IsochoricRealGasParcel::new(model, point, VOLUME, PengRobinsonPhaseV1::Vapor).unwrap()
}

fn solve() -> ModalGasSolve {
    ModalGasSolve {
        gas: GasEnergySolve {
            minimum_temperature_k: 280.0,
            maximum_temperature_k: 450.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
        energy_tolerance_j: 5e-7,
        ..controls()
    }
}

#[test]
fn dense_gas_drives_the_same_elastic_basis_with_native_pressure_work_and_calorics() {
    with_cx(|cx| {
        let gas = parcel();
        let ideal = ClosedGasParcel::new(
            gas.model().source().clone(),
            gas.point().clone(),
            gas.volume_m3(),
        )
        .unwrap();
        assert!((gas.state().density_kg_m3() / ideal.state().state().density - 1.0).abs() > 0.05);
        let wall = motion(1000.0, 0.0, 0.1, AREA, RATE, cx);
        let identity = wall.valve().identity();
        let initial_energy = energy(&wall);
        let initial_u = gas.state().specific_internal_energy_j_kg();
        let mut real_run = GasValveMotion::new(gas.clone(), wall.clone(), cx).unwrap();
        let mut ideal_run = GasValveMotion::new(ideal, wall, cx).unwrap();
        let (mut supplied_heat, mut loss, mut exterior_work) = (0.0, 0.0, 0.0);
        for step in 1..=8 {
            let old_volume = real_run.gas().volume_m3();
            let report = real_run.step(25.0, PRESSURE, solve(), cx).unwrap();
            ideal_run.step(25.0, PRESSURE, solve(), cx).unwrap();
            let physical = real_run.motion().nodal_motion(cx).unwrap();
            let dv = real_run.gas().volume_m3() - old_volume;
            near(
                report.gas.work_by_gas_j,
                report.midpoint.pressure_pa() * dv,
                1e-12,
            );
            near(
                report.motion.modal.input_work_j,
                (report.midpoint.pressure_pa() - PRESSURE) * report.motion.swept_volume_change_m3,
                1e-11,
            );
            near(
                real_run.gas().volume_m3() - VOLUME,
                AREA * physical.displacement_m[7][0],
                solve().volume_tolerance_m3,
            );
            assert!(report.pressure_work_residual_j.abs() <= solve().energy_tolerance_j);
            assert!(report.total_energy_residual_j.abs() <= solve().energy_tolerance_j);
            assert_eq!(real_run.gas().mass_kg(), gas.mass_kg());
            assert_eq!(real_run.gas().model(), gas.model());
            assert_eq!(
                report.midpoint.constitutive_law(),
                gas.state().constitutive_law()
            );
            assert_eq!(
                real_run.gas().point().axes().get("specimen-lot"),
                Some(&42.0)
            );
            assert_eq!(real_run.motion().valve().identity(), identity);
            assert_eq!(real_run.motion().accepted_steps(), step);
            supplied_heat += report.gas.supplied_heat_j;
            loss += report.motion.modal.viscous_dissipation_j;
            exterior_work += report.exterior_pressure_work_j;
        }
        let stored =
            gas.mass_kg() * (real_run.gas().state().specific_internal_energy_j_kg() - initial_u);
        near(
            stored + energy(real_run.motion()) - initial_energy + loss + exterior_work,
            supplied_heat,
            8.0 * solve().energy_tolerance_j,
        );
        let real_wall = real_run.motion().nodal_motion(cx).unwrap();
        let ideal_wall = ideal_run.motion().nodal_motion(cx).unwrap();
        assert!(real_wall.displacement_m[7][0] > 0.0);
        assert!(real_wall.velocity_m_s[7][0] > 0.0);
        assert!((real_wall.velocity_m_s[7][0] - ideal_wall.velocity_m_s[7][0]).abs() > 1e-8);
        assert!(
            (real_run.gas().state().pressure_pa() - ideal_run.gas().state().state().pressure).abs()
                > 1.0
        );
    });
}

#[test]
fn real_gas_equilibrium_keeps_inventory_and_elastic_state_on_the_accepted_clock() {
    with_cx(|cx| {
        let gas = parcel();
        let wall = motion(1000.0, 0.0, 0.0, AREA, RATE, cx);
        let states = wall.model().states().to_vec();
        let mut run = GasValveMotion::new(gas.clone(), wall, cx).unwrap();
        for step in 1..=4 {
            let report = run
                .step(0.0, gas.state().pressure_pa(), solve(), cx)
                .unwrap();
            assert_eq!(run.gas(), &gas);
            assert_eq!(report.midpoint, *gas.state());
            assert_eq!(run.motion().model().states(), states);
            assert_eq!(run.motion().accepted_steps(), step);
            assert_eq!(report.gas.work_by_gas_j, 0.0);
            assert_eq!(report.total_energy_residual_j, 0.0);
            assert_eq!(report.motion.time_s, step as f64 / RATE as f64);
        }
    });
}

#[test]
fn real_eos_bulk_modulus_controls_the_coupled_clock() {
    with_cx(|cx| {
        let gas = parcel();
        let wall = motion(1000.0, 0.0, 0.0, AREA, 2000, cx);
        let shape = wall.valve().valve().source_modes()[0].phi[0];
        let elastic = wall.model().modes()[0].angular_frequency_rad_s;
        // One free nodal coordinate: the physical pressure force is A*p and
        // its unit-mass overlap is A*phi. The EOS supplies the adiabatic bulk.
        let bulk = gas.state().density_kg_m3() * gas.state().sound_speed_m_s().powi(2);
        let expected = (elastic.powi(2) + bulk * (AREA * shape).powi(2) / VOLUME).sqrt();
        let reference = gas.state().source().state();
        let wrong = (elastic.powi(2)
            + reference.density * reference.sound_speed.powi(2) * (AREA * shape).powi(2) / VOLUME)
            .sqrt();
        assert!((expected / wrong - 1.0).abs() > 1e-3);
        match GasValveMotion::new(gas, wall, cx).unwrap_err() {
            ModalGasError::TimeResolution(fs_phs::linear_exact::ResolvedBandError::Outside {
                rate_rad_s,
                ..
            }) => near(rate_rad_s, expected, 1e-12 * expected),
            other => panic!("expected the real-gas clock refusal, got {other:?}"),
        }
    });
}

#[test]
fn real_gas_source_refusal_and_cancellation_preserve_both_accepted_histories() {
    with_cx(|cx| {
        let mut run =
            GasValveMotion::new(parcel(), motion(1000.0, 0.0, 0.1, AREA, RATE, cx), cx).unwrap();
        run.step(10.0, PRESSURE, solve(), cx).unwrap();
        let before = run.clone();
        let unchanged = |candidate: &GasValveMotion<IsochoricRealGasParcel>| {
            assert_eq!(candidate.gas(), before.gas());
            assert_eq!(
                candidate.motion().model().states(),
                before.motion().model().states()
            );
            assert_eq!(
                candidate.motion().accepted_steps(),
                before.motion().accepted_steps()
            );
            assert_eq!(
                candidate.motion().valve().identity(),
                before.motion().valve().identity()
            );
        };
        let mut unsupported = solve();
        unsupported.gas.maximum_temperature_k = 2100.0; // Original NASA data stops at 2000 K.
        assert!(matches!(
            run.step(25.0, PRESSURE, unsupported, cx),
            Err(ModalGasError::Gas(GasEnergyError::Material(_)))
        ));
        unchanged(&run);
        let mut limited = solve();
        limited.maximum_trials = 1;
        assert!(matches!(
            run.step(25.0, PRESSURE, limited, cx),
            Err(ModalGasError::TrialBudget)
        ));
        unchanged(&run);
        let gate = CancelGate::new_clock_free();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(matches!(
                run.step(25.0, PRESSURE, solve(), &cancelled),
                Err(ModalGasError::Cancelled)
            ));
        });
        unchanged(&run);
        let mut clean = before.clone();
        let retried = run.step(25.0, PRESSURE, solve(), cx).unwrap();
        let expected = clean.step(25.0, PRESSURE, solve(), cx).unwrap();
        assert_eq!(retried.midpoint, expected.midpoint);
        assert_eq!(retried.gas, expected.gas);
        assert_eq!(run.gas(), clean.gas());
        assert_eq!(
            run.motion().model().states(),
            clean.motion().model().states()
        );
        assert_eq!(run.motion().accepted_steps(), 2);
    });
}
