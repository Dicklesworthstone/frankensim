//! G0/G3: the same gas drives actual D3Q19 flow and physical sound-source data.
//! Short transients test integration, not a converged broadband spectrum or SPL.

use fs_aeroac::SCOPE_STATEMENT;
use fs_aeroac::material::{JetPhysicalState, MaterialSlotJet3d, SlotJet3dSetup};
use fs_aeroac::slot_jet_3d::run_slot_jet_3d;
use fs_lbm::d3q19::CollisionModel3;
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::gas::GasState;
use fs_material::gas_species::input::read_gas_state;
use fs_math::c64::C64;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;
use fs_thermochem::water::transport::If97TransportState;

fn setup() -> SlotJet3dSetup {
    SlotJet3dSetup {
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
    }
}

fn physical() -> JetPhysicalState {
    JetPhysicalState {
        cell_size_m: 1e-4,
        timestep_s: 1e-5,
        jet_speed_m_s: 0.5,
        max_speed_m_s: 1.4,
        max_cells: 64 * 32 * 8,
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

fn second_order(model: &MaterialSlotJet3d<'_>) -> f64 {
    match model.lattice_config().collision {
        CollisionModel3::CentralMoment {
            second_order_rate,
            higher_order_rate,
        } => {
            assert_eq!(higher_order_rate, 1.2);
            second_order_rate
        }
        _ => panic!("material-bound 3-D rig must keep its central-moment owner"),
    }
}

#[test]
fn g0_actual_3d_force_uses_volume_mass_and_shared_physical_spectrum() {
    let gas = gas("300 100000 0");
    let physical = physical();
    let model = MaterialSlotJet3d::try_new(&gas, &setup(), physical).unwrap();
    let lattice = run_slot_jet_3d(model.lattice_config()).unwrap();
    let run = model.run().unwrap();
    let dx = physical.cell_size_m;
    let dt = physical.timestep_s;
    // A 3-D impulse belongs to a volume, not to unit span: N, not N/m.
    let force_unit_n = gas.density * dx.powi(4) / dt.powi(2);
    close(
        model.lattice_config().u_jet,
        physical.jet_speed_m_s * dt / dx,
    );
    close(
        1.0 / second_order(&model),
        0.5 + 3.0 * gas.dynamic_viscosity / gas.density * dt / dx.powi(2),
    );
    for (actual, original) in run.force_series_n().iter().zip(&lattice.force_series) {
        for axis in 0..2 {
            close(actual[axis], original[axis] * force_unit_n);
        }
    }
    close(run.span_m(), 8.0 * dx);
    close(
        run.lattice_diagnostics.reynolds,
        gas.density * physical.jet_speed_m_s * 6.0 * dx / gas.dynamic_viscosity,
    );
    close(
        run.final_plane_physical_mach_max,
        lattice.diagnostics.mach_max_lattice * dx / dt / gas.sound_speed,
    );
    assert_eq!(run.gas(), &gas);
    assert_eq!(run.scope(), SCOPE_STATEMENT);

    let bin = 4;
    let line = run.force_spectral_line(bin).unwrap();
    let n = run.force_series_n().len() as f64;
    close(line.frequency_hz, bin as f64 / (n * dt));
    close(
        line.acoustic_wavenumber_rad_m,
        core::f64::consts::TAU * line.frequency_hz / gas.sound_speed,
    );
    // Independent direct DFT, including physical amplitude and phase.
    let means = [0, 1].map(|axis| run.force_series_n().iter().map(|f| f[axis]).sum::<f64>() / n);
    let mut expected = [C64::ZERO; 2];
    for (i, force) in run.force_series_n().iter().enumerate() {
        let phase = core::f64::consts::TAU * bin as f64 * i as f64 / n;
        let window = 0.5 - 0.5 * (core::f64::consts::TAU * i as f64 / (n - 1.0)).cos();
        for axis in 0..2 {
            expected[axis] = expected[axis]
                + C64::new(phase.cos(), phase.sin())
                    .scale(2.0 * (force[axis] - means[axis]) * window / (0.5 * (n - 1.0)));
        }
    }
    for axis in 0..2 {
        assert!((line.force_n[axis] - expected[axis]).abs() < 2e-10 * expected[axis].abs());
    }
    assert!(line.force_n[1].abs() > 0.0);
    for bin in [0, 32, usize::MAX] {
        assert!(run.force_spectral_line(bin).is_err());
    }
}

#[test]
fn g3_pressure_temperature_and_composition_reach_the_3d_solver() {
    let states = [
        gas("mixture 300 100000 mole 2 N2 0.8 O2 0.2"),
        gas("mixture 300 200000 mole 2 N2 0.8 O2 0.2"),
        gas("mixture 300 100000 mole 3 N2 0.4 O2 0.1 CO2 0.5"),
        gas("mixture 320 100000 mole 2 N2 0.8 O2 0.2"),
    ];
    let models: Vec<_> = states
        .iter()
        .map(|gas| MaterialSlotJet3d::try_new(gas, &setup(), physical()).unwrap())
        .collect();
    let runs: Vec<_> = models.iter().map(|model| model.run().unwrap()).collect();
    let lines: Vec<_> = runs
        .iter()
        .map(|run| run.force_spectral_line(4).unwrap())
        .collect();
    close(
        runs[1].lattice_diagnostics.reynolds,
        2.0 * runs[0].lattice_diagnostics.reynolds,
    );
    // Ideal-gas pressure changes inertia and viscosity/Re, while c stays fixed.
    assert_eq!(
        lines[0].acoustic_wavenumber_rad_m,
        lines[1].acoustic_wavenumber_rad_m
    );
    assert!(lines[2].acoustic_wavenumber_rad_m > lines[0].acoustic_wavenumber_rad_m);
    assert!(lines[3].acoustic_wavenumber_rad_m < lines[0].acoustic_wavenumber_rad_m);
    for i in 1..states.len() {
        assert_eq!(
            models[i].lattice_config().u_jet,
            models[0].lattice_config().u_jet
        );
        assert_ne!(second_order(&models[i]), second_order(&models[0]));
        assert_ne!(runs[i].force_series_n(), runs[0].force_series_n());
        assert!(
            (lines[i].force_n[1] - lines[0].force_n[1]).abs() > 1e-5 * lines[0].force_n[1].abs()
        );
        assert_eq!(runs[i].gas(), &states[i]);
    }
}

#[test]
fn g0_incoherent_material_mach_and_allocation_refuse_before_execution() {
    let good = gas("300 100000 0");
    let mut stale = good;
    stale.dynamic_viscosity *= 2.0;
    assert!(MaterialSlotJet3d::try_new(&stale, &setup(), physical()).is_err());
    let high_physical_mach = JetPhysicalState {
        jet_speed_m_s: 200.0,
        max_speed_m_s: 400.0,
        timestep_s: 1e-8,
        ..physical()
    };
    assert!(MaterialSlotJet3d::try_new(&good, &setup(), high_physical_mach).is_err());
    let high_lattice_mach = JetPhysicalState {
        timestep_s: 0.01,
        ..physical()
    };
    assert!(MaterialSlotJet3d::try_new(&good, &setup(), high_lattice_mach).is_err());
    let small_budget = JetPhysicalState {
        max_cells: 64 * 32 * 8 - 1,
        ..physical()
    };
    assert!(MaterialSlotJet3d::try_new(&good, &setup(), small_budget).is_err());
    let overflow = SlotJet3dSetup {
        nz: usize::MAX,
        ..setup()
    };
    assert!(MaterialSlotJet3d::try_new(&good, &overflow, physical()).is_err());
    // The real seed has a 1.5 spanwise multiplier. A 2-D bound would admit this.
    let seeded = SlotJet3dSetup {
        seed_amplitude: 0.19,
        ..setup()
    };
    let short_envelope = JetPhysicalState {
        max_speed_m_s: 0.51,
        ..physical()
    };
    assert!(MaterialSlotJet3d::try_new(&good, &seeded, short_envelope).is_err());
    let invalid_rate = SlotJet3dSetup {
        higher_order_rate: 2.1,
        ..setup()
    };
    assert!(MaterialSlotJet3d::try_new(&good, &invalid_rate, physical()).is_err());
}

#[test]
fn g1_original_liquid_source_reaches_native_3d_flow_and_force_spectrum() {
    let water = If97TransportState::at_pressure_temperature(
        Pressure::new(2e6),
        Temperature::new(300.0),
        DensePressureSolve {
            pressure_tolerance_pa: 1e-4,
            density_tolerance_kg_m3: 1e-8,
            maximum_iterations: 128,
        },
        || false,
    )
    .unwrap();
    let physical = JetPhysicalState {
        cell_size_m: 1e-5,
        timestep_s: 1e-6,
        ..physical()
    };
    let model = MaterialSlotJet3d::try_new(&water, &setup(), physical).unwrap();
    let run = model.run().unwrap();
    let lattice = run_slot_jet_3d(model.lattice_config()).unwrap();
    let rho = water.thermodynamics().density().value();
    let speed = water.thermodynamics().sound_speed().value();
    let dx = physical.cell_size_m;
    let dt = physical.timestep_s;
    assert!(core::ptr::eq(model.scaling().material(), &water));
    assert!(core::ptr::eq(run.material(), &water));
    assert_eq!(model.acoustic(), water.acoustic_state().unwrap());
    assert_eq!(run.acoustic(), model.acoustic());
    let rate = match model.lattice_config().collision {
        CollisionModel3::CentralMoment {
            second_order_rate, ..
        } => second_order_rate,
        _ => panic!("original central-moment solver required"),
    };
    close(
        1.0 / rate,
        0.5 + 3.0 * water.dynamic_viscosity().value() / rho * dt / (dx * dx),
    );
    for (actual, original) in run.force_series_n().iter().zip(lattice.force_series) {
        for axis in 0..2 {
            close(actual[axis], original[axis] * rho * dx.powi(4) / dt.powi(2));
        }
    }
    close(
        run.final_plane_physical_mach_max,
        lattice.diagnostics.mach_max_lattice * dx / dt / speed,
    );
    let line = run.force_spectral_line(4).unwrap();
    close(
        line.acoustic_wavenumber_rad_m,
        core::f64::consts::TAU * line.frequency_hz / speed,
    );
    assert!(line.force_n[1].abs() > 0.0);
}
