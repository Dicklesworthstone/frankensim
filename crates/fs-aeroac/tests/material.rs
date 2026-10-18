//! G0/G3: one material state drives actual flow, SI force and radiation.
//! These short transient runs verify coupling, not a converged edge-tone or SPL.

use fs_aeroac::SCOPE_STATEMENT;
use fs_aeroac::curle2d::dipole_pressure;
use fs_aeroac::jetlab::run_jet_labium;
use fs_aeroac::material::{
    JetLabiumSetup, JetPhysicalState, MaterialJetLabium, MaterialJetLabiumRun,
};
use fs_lbm::core2::CollisionModel2;
use fs_material::gas::GasState;
use fs_material::gas_species::input::read_gas_state;
use fs_math::{c64::C64, det};

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
