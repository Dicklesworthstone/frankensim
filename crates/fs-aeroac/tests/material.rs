//! G0/G3: one material state drives actual flow, SI force and radiation.
//! These short transient runs verify coupling, not a converged edge-tone or SPL.

use fs_aeroac::SCOPE_STATEMENT;
use fs_aeroac::curle2d::dipole_pressure;
use fs_aeroac::jetlab::run_jet_labium;
use fs_aeroac::material::{
    JetLabiumSetup, JetPhysicalState, MaterialJetLabium, MaterialJetLabiumRun,
};
use fs_lbm::core2::CollisionModel2;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::GasState;
use fs_material::gas_species::input::read_gas_state;
use fs_math::{c64::C64, det};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;
use fs_thermochem::water::transport::If97TransportState;

fn setup() -> JetLabiumSetup {
    JetLabiumSetup {
        nx: 64,
        ny: 32,
        slot_half: 3.0,
        slot_smoothing: 1.0,
        edge_distance: 12,
        plate_length: 16,
        fringe_width: 12,
        fringe_sigma: 0.3,
        steps_settle: 32,
        steps_record: 64,
        seed_amplitude: 0.02,
        nozzle_thickness: 2,
        collision: CollisionModel2::Bgk,
    }
}

fn physical() -> JetPhysicalState {
    JetPhysicalState {
        cell_size_m: 1e-4,
        timestep_s: 1e-5,
        jet_speed_m_s: 0.5,
        max_speed_m_s: 1.4,
        max_cells: 64 * 32,
    }
}

fn gas(record: &str) -> GasState {
    let mut words = record.split_ascii_whitespace();
    let state = read_gas_state(&mut words, None).unwrap();
    assert!(words.next().is_none());
    state
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-30),
        "{actual:e} != {expected:e}"
    );
}

// Independent Hann-window DFT of a real force history. Positive
// exponent returns amplitudes for the workspace's exp(-i omega t) convention.
fn force_line(run: &MaterialJetLabiumRun<'_>, bin: usize) -> ([C64; 2], f64) {
    let n = run.force_series_n_m().len() as f64;
    let mut force = [C64::new(0.0, 0.0); 2];
    let means = [0, 1].map(|axis| run.force_series_n_m().iter().map(|f| f[axis]).sum::<f64>() / n);
    let window_sum = 0.5 * (n - 1.0);
    for (i, value) in run.force_series_n_m().iter().enumerate() {
        let phase = core::f64::consts::TAU * bin as f64 * i as f64 / n;
        let window = 0.5 - 0.5 * (core::f64::consts::TAU * i as f64 / (n - 1.0)).cos();
        for axis in 0..2 {
            force[axis] = force[axis]
                + C64::new(det::cos(phase), det::sin(phase))
                    .scale(2.0 * (value[axis] - means[axis]) * window / window_sum);
        }
    }
    (force, bin as f64 / (n * run.timestep_s()))
}

#[test]
fn g0_actual_lattice_force_and_frequency_have_correct_physical_units() {
    let gas = gas("300 100000 0");
    let physical = physical();
    let model = MaterialJetLabium::try_new(&gas, &setup(), physical).unwrap();
    let lattice = run_jet_labium(model.lattice_config()).unwrap();
    let run = model.run().unwrap();
    let dx = physical.cell_size_m;
    let dt = physical.timestep_s;
    let force_unit = gas.density * dx.powi(3) / dt.powi(2);
    close(
        model.lattice_config().u_jet,
        physical.jet_speed_m_s * dt / dx,
    );
    close(
        model.lattice_config().tau,
        0.5 + 3.0 * gas.dynamic_viscosity / gas.density * dt / dx.powi(2),
    );
    for (actual, original) in run.force_series_n_m().iter().zip(&lattice.force_series) {
        for axis in 0..2 {
            close(actual[axis], original[axis] * force_unit);
        }
    }
    close(
        run.lattice_diagnostics.reynolds,
        gas.density * physical.jet_speed_m_s * 6.0 * dx / gas.dynamic_viscosity,
    );
    close(
        run.sampled_physical_mach_max,
        lattice.diagnostics.mach_max_lattice * dx / dt / gas.sound_speed,
    );
    assert_eq!(run.gas(), &gas);
    assert_eq!(run.scope(), SCOPE_STATEMENT);
    let (force, frequency) = force_line(&run, 4);
    let observer = [0.2, 0.1];
    let source = [0.0, 0.0];
    let observed = run.observe_spectral_line(4, observer, source).unwrap();
    close(observed.frequency_hz, frequency);
    for axis in 0..2 {
        assert!((observed.force_n_m[axis] - force[axis]).abs() < 2e-11 * force[axis].abs());
    }
    let actual = run
        .radiate_force_line(force, frequency, observer, source)
        .unwrap();
    let expected = dipole_pressure(
        force,
        core::f64::consts::TAU * frequency / gas.sound_speed,
        observer,
        source,
    )
    .unwrap();
    close(actual.pressure.re, expected.pressure.re);
    close(actual.pressure.im, expected.pressure.im);
    assert!((observed.field.pressure - expected.pressure).abs() < 2e-11 * expected.pressure.abs());
    assert!(actual.pressure.abs() > 0.0);
    assert_eq!(actual.scope, SCOPE_STATEMENT);
    let prepared = run
        .prepare_coherent_observer(frequency, observer, &[source, source], 2)
        .unwrap();
    assert_eq!(prepared.acoustic(), run.acoustic());
    assert_eq!(prepared.frequency_hz(), frequency);
    let summed = prepared.evaluate(&[force, force]).unwrap();
    close(summed.pressure_pa.re, 2.0 * expected.pressure.re);
    close(summed.pressure_pa.im, 2.0 * expected.pressure.im);
    let independent = prepared
        .evaluate_independent_sources(&[force, force])
        .unwrap();
    close(
        independent.mean_square_pressure_pa2,
        expected.pressure.norm_sq(),
    );
    for axis in 0..2 {
        close(
            independent.active_intensity_w_m2[axis],
            0.5 * summed.active_intensity_w_m2[axis],
        );
    }
    assert!(
        run.prepare_coherent_observer(0.5 / dt, observer, &[source], 1)
            .is_err()
    );
    assert!(
        run.prepare_coherent_observer(frequency, observer, &[source, source], 1)
            .is_err()
    );
    assert!(run.observe_spectral_line(0, observer, source).is_err());
    assert!(run.observe_spectral_line(32, observer, source).is_err());
    assert!(
        run.observe_spectral_line(usize::MAX, observer, source)
            .is_err()
    );
    assert!(
        run.radiate_force_line(force, 0.5 / dt, observer, source)
            .is_err()
    );
    assert!(
        run.radiate_force_line(force, f64::NAN, observer, source)
            .is_err()
    );
}

#[test]
fn g3_data_only_pressure_temperature_and_composition_swaps_reach_flow_and_sound() {
    let states = [
        gas("mixture 300 100000 mole 2 N2 0.8 O2 0.2"),
        gas("mixture 300 200000 mole 2 N2 0.8 O2 0.2"),
        gas("mixture 300 100000 mole 3 N2 0.4 O2 0.1 CO2 0.5"),
        gas("mixture 320 100000 mole 2 N2 0.8 O2 0.2"),
    ];
    let setup = setup();
    let models: Vec<_> = states
        .iter()
        .map(|gas| MaterialJetLabium::try_new(gas, &setup, physical()).unwrap())
        .collect();
    let runs: Vec<_> = models.iter().map(|model| model.run().unwrap()).collect();
    let observed: Vec<_> = runs
        .iter()
        .map(|run| {
            run.observe_spectral_line(4, [0.2, 0.1], [0.0, 0.0])
                .unwrap()
        })
        .collect();
    close(states[1].density, 2.0 * states[0].density);
    close(states[1].sound_speed, states[0].sound_speed);
    close(
        runs[1].lattice_diagnostics.reynolds,
        2.0 * runs[0].lattice_diagnostics.reynolds,
    );
    for i in 1..states.len() {
        // Geometry, imposed speed and timestep remain exactly fixed: changing
        // the material record changes the real solver's viscosity and forces.
        assert_eq!(
            models[i].lattice_config().u_jet,
            models[0].lattice_config().u_jet
        );
        assert_ne!(
            models[i].lattice_config().tau,
            models[0].lattice_config().tau
        );
        assert_ne!(runs[i].force_series_n_m(), runs[0].force_series_n_m());
        assert!(
            (observed[i].field.pressure - observed[0].field.pressure).abs()
                > 1e-5 * observed[0].field.pressure.abs()
        );
    }
    let (force, frequency) = force_line(&runs[0], 4);
    let pressure: Vec<_> = runs
        .iter()
        .map(|run| {
            run.radiate_force_line(force, frequency, [0.2, 0.1], [0.0, 0.0])
                .unwrap()
                .pressure
        })
        .collect();
    // Holding the source force fixed isolates propagation. Ideal-gas pressure
    // alone leaves c unchanged; changing molecular composition or T changes c.
    assert_eq!(pressure[0], pressure[1]);
    assert!(states[2].sound_speed < states[0].sound_speed);
    assert!(states[3].sound_speed > states[0].sound_speed);
    assert!((pressure[2] - pressure[0]).abs() > 1e-3 * pressure[0].abs());
    assert!((pressure[3] - pressure[0]).abs() > 1e-3 * pressure[0].abs());
}

#[test]
fn g0_material_and_physical_applicability_refuse_before_running() {
    let good = gas("300 100000 0");
    let setup = setup();
    let mut stale = good;
    stale.density *= 1.01;
    assert!(MaterialJetLabium::try_new(&stale, &setup, physical()).is_err());
    let high_physical_mach = JetPhysicalState {
        jet_speed_m_s: 200.0,
        max_speed_m_s: 400.0,
        timestep_s: 1e-8,
        ..physical()
    };
    // Lattice speed is only 0.04, yet the real gas is outside low-Mach physics.
    assert!(MaterialJetLabium::try_new(&good, &setup, high_physical_mach).is_err());
    let high_lattice_mach = JetPhysicalState {
        timestep_s: 0.01,
        ..physical()
    };
    assert!(MaterialJetLabium::try_new(&good, &setup, high_lattice_mach).is_err());
    let small_budget = JetPhysicalState {
        max_cells: 64 * 32 - 1,
        ..physical()
    };
    assert!(MaterialJetLabium::try_new(&good, &setup, small_budget).is_err());
    let small_envelope = JetPhysicalState {
        max_speed_m_s: 0.5,
        ..physical()
    };
    assert!(MaterialJetLabium::try_new(&good, &setup, small_envelope).is_err());
    let overflow_geometry = JetLabiumSetup {
        plate_length: usize::MAX,
        ..setup
    };
    assert!(MaterialJetLabium::try_new(&good, &overflow_geometry, physical()).is_err());
}

#[test]
fn g3_original_water_eos_drives_flow_force_and_same_state_radiation() {
    let physical = JetPhysicalState {
        cell_size_m: 1e-5,
        timestep_s: 1e-6,
        ..physical()
    };
    let states = [300.0, 330.0].map(|t| {
        If97TransportState::at_pressure_temperature(
            Pressure::new(2e6),
            Temperature::new(t),
            DensePressureSolve {
                pressure_tolerance_pa: 1e-4,
                density_tolerance_kg_m3: 1e-8,
                maximum_iterations: 128,
            },
            || false,
        )
        .unwrap()
    });
    let models: Vec<_> = states
        .iter()
        .map(|water| MaterialJetLabium::try_new(water, &setup(), physical).unwrap())
        .collect();
    let runs: Vec<_> = models.iter().map(|model| model.run().unwrap()).collect();
    for (i, water) in states.iter().enumerate() {
        let thermo = water.thermodynamics();
        let rho = thermo.density().value();
        let speed = thermo.sound_speed().value();
        let nu = water.dynamic_viscosity().value() / rho;
        let dx = physical.cell_size_m;
        let dt = physical.timestep_s;
        close(
            models[i].lattice_config().tau,
            0.5 + 3.0 * nu * dt / (dx * dx),
        );
        assert_eq!(models[i].acoustic(), water.acoustic_state().unwrap());
        assert!(core::ptr::eq(models[i].scaling().material(), water));
        assert!(core::ptr::eq(runs[i].material(), water));
        assert_eq!(runs[i].acoustic().sound_speed_m_s(), speed);
        let lattice = run_jet_labium(models[i].lattice_config()).unwrap();
        for (actual, original) in runs[i].force_series_n_m().iter().zip(lattice.force_series) {
            for axis in 0..2 {
                close(actual[axis], original[axis] * rho * dx.powi(3) / dt.powi(2));
            }
        }
        let line = runs[i]
            .observe_spectral_line(4, [0.2, 0.1], [0.0, 0.0])
            .unwrap();
        let expected = dipole_pressure(
            line.force_n_m,
            core::f64::consts::TAU * line.frequency_hz / speed,
            [0.2, 0.1],
            [0.0, 0.0],
        )
        .unwrap();
        assert_eq!(line.field.pressure, expected.pressure);
        assert!(line.field.pressure.abs() > 0.0);
    }
    assert_ne!(
        models[0].lattice_config().tau,
        models[1].lattice_config().tau
    );
    assert_ne!(runs[0].force_series_n_m(), runs[1].force_series_n_m());
    assert_ne!(
        runs[0].acoustic().sound_speed_m_s(),
        runs[1].acoustic().sound_speed_m_s()
    );
    // Hold the physical source fixed to isolate the original water EOS in propagation.
    let force = [C64::new(0.0, 0.0), C64::new(1.0, 0.5)];
    let pressures: Vec<_> = runs
        .iter()
        .map(|run| {
            run.radiate_force_line(force, 5000.0, [0.2, 0.1], [0.0, 0.0])
                .unwrap()
                .pressure
        })
        .collect();
    assert_ne!(pressures[0], pressures[1]);
}

#[test]
fn g0_mixed_transport_and_acoustic_density_refuses_both_native_rigs() {
    struct IncoherentFluid;
    impl FluidTransport for IncoherentFluid {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            FluidTransportState::try_new(1000.0, 4200.0, 0.6, 0.001)
        }
    }
    impl FluidAcoustics for IncoherentFluid {
        fn acoustic_state(&self) -> Result<FluidAcousticState, fs_material::MaterialError> {
            FluidAcousticState::try_new(900.0, 1500.0)
        }
    }
    let error = MaterialJetLabium::try_new(&IncoherentFluid, &setup(), physical())
        .err()
        .unwrap();
    assert!(error.to_string().contains("same density"));
    let setup3 = fs_aeroac::material::SlotJet3dSetup {
        nx: 64,
        ny: 32,
        nz: 8,
        slot_half: 3.0,
        higher_order_rate: 1.2,
        nozzle_thickness: 2,
        edge_distance: 12,
        plate_length: 16,
        fringe_width: 12,
        fringe_sigma: 0.3,
        seed_amplitude: 0.02,
        steps_settle: 16,
        steps_record: 64,
    };
    let error = fs_aeroac::material::MaterialSlotJet3d::try_new(
        &IncoherentFluid,
        &setup3,
        JetPhysicalState {
            max_cells: 64 * 32 * 8,
            ..physical()
        },
    )
    .err()
    .unwrap();
    assert!(error.to_string().contains("same density"));
}
