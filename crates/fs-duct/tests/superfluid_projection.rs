//! G1: explicit in-phase superfluid density waves use the original acoustic
//! frequency/time-domain owners. These authored properties are synthetic.

use fs_duct::compose::ComposeOptions;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh};
use fs_duct::{ConeSlices, Duct, Segment, WallModel, fluid};
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::superfluid::{TwoFluidProperties, TwoFluidState};
use fs_math::c64::C64;

fn near(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn in_phase_projection_reaches_native_impedance_and_energy_preserving_dg() {
    let properties = TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: 0.0,
        normal_viscosity_pa_s: 0.03,
    };
    let material = TwoFluidState::try_new(properties).unwrap();
    let first_sound = material.first_sound_projection().unwrap();
    let acoustic = first_sound.acoustic_state().unwrap();
    let radius: f64 = 0.01;
    let length = 1.0;
    let area = core::f64::consts::PI * radius * radius;
    let duct = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    let omega = 4.0;
    let wave = fluid::segment_wave(&first_sound, radius, omega).unwrap();
    near(wave.wavenumber.re, omega / 10.0, 1.0, 1e-14);
    near(
        wave.characteristic_impedance.re,
        40.0 / area,
        40.0 / area,
        1e-14,
    );
    assert_eq!(wave.wavenumber.im, 0.0);
    assert_eq!(wave.characteristic_impedance.im, 0.0);

    let zc = 40.0 / area;
    let response = fluid::input_impedance_load(
        &duct,
        &first_sound,
        omega,
        Some(C64::from_re(zc)),
        WallModel::Rigid,
        ConeSlices::Auto,
        ComposeOptions::default(),
        || false,
    )
    .unwrap();
    near(response.z_in.re, zc, zc, 3e-12);
    near(response.z_in.im, 0.0, zc, 3e-12);
    let phase = omega * length / 10.0;
    near(response.u_out_over_u_in.re, phase.cos(), 1.0, 3e-12);
    near(response.u_out_over_u_in.im, phase.sin(), 1.0, 3e-12);

    let dg = DgDuct::from_material(
        &duct,
        &first_sound,
        DgMesh {
            elements: 4,
            order: 4,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap();
    assert_eq!(dg.material().source(), &material);
    assert_eq!(dg.acoustic_material(), acoustic);
    let amplitude = 0.01;
    let mut state = dg.project_initial_state(
        |x| amplitude * (core::f64::consts::PI * x / length).cos(),
        |_| 0.0,
    );
    let energy = state.stored_energy_j();
    let expected_energy = area * length * amplitude * amplitude / (4.0 * 400.0);
    near(energy, expected_energy, expected_energy, 2e-5);
    let dt = length / (8.0 * 10.0);
    let propagator = dg.propagator(dt, 1e-10).unwrap();
    for step in 1..=8 {
        let report = propagator.step(&mut state, &[]).unwrap();
        assert!(report.ledger.balance_checked);
        let phase = (core::f64::consts::PI * step as f64 / 8.0).cos();
        let pressure_error = dg
            .l2_pressure_error(&state, |x| {
                amplitude * phase * (core::f64::consts::PI * x / length).cos()
            })
            .unwrap();
        assert!(pressure_error / amplitude < 1e-4);
        near(state.stored_energy_j(), energy, energy, 3e-10);
    }

    // Thermal expansion couples pressure to entropy; that state must use the
    // full two-fluid channel rather than silently dropping its second mode.
    let coupled = TwoFluidState::try_new(TwoFluidProperties {
        thermal_expansion_per_k: -0.01,
        ..properties
    })
    .unwrap();
    assert!(coupled.first_sound_projection().is_err());
}
