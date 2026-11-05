//! G1/G3 compact radiation, independent of the coupled pHS state equations.

use super::*;
use fs_material::gas::GasSpec;
use fs_math::c64::C64;
use fs_phs::linear_exact::{LinearPropagator, split_quadratic, step_split};

fn panel() -> PlateCavitySpec {
    PlateCavitySpec {
        omegas: vec![600.0],
        zetas: vec![0.08],
        drive: vec![1.0],
        areas: vec![0.08],
        volume_m3: 0.008,
        neck_radius_m: 0.02,
        neck_length_m: 0.03,
        inner_end: fs_scenario::NeckEnd::Flanged,
        outer_end: fs_scenario::NeckEnd::Flanged,
        gas: GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap(),
        relative_humidity: 0.0,
    }
}

fn geometry(spec: &PlateCavitySpec) -> fs_scenario::HelmholtzCavity {
    fs_scenario::HelmholtzCavity {
        volume_m3: spec.volume_m3,
        neck_radius_m: spec.neck_radius_m,
        neck_length_m: spec.neck_length_m,
        inner_end: spec.inner_end,
        outer_end: spec.outer_end,
    }
}

// Continuous physical equations with the e^(-i omega t) convention:
// (s² + 2 zeta omega_s s + omega_s²) q = drive F + A p,
// (M s + R) U_out = p, C s p = -A s q - U_out.
// R is the declared loss coefficient, not an independently fitted oscillator.
fn harmonic_pressure(spec: &PlateCavitySpec, frequency: f64, distance: f64) -> (C64, C64) {
    let real = C64::from_re;
    let s = C64::new(0.0, -2.0 * core::f64::consts::PI * frequency);
    let compliance = spec.volume_m3 / (spec.gas.density * spec.gas.sound_speed.powi(2));
    let neck_area = core::f64::consts::PI * spec.neck_radius_m.powi(2);
    // Independent compact-flanged end correction, one for each neck end.
    let mass =
        spec.gas.density * (spec.neck_length_m + 2.0 * 0.8216 * spec.neck_radius_m) / neck_area;
    let resistance =
        crate::cavity_modal::resonance_neck_resistance(&geometry(spec), &spec.gas).unwrap();
    let neck_impedance = real(mass) * s + real(resistance);
    let pressure_per_displacement =
        -s * real(spec.areas[0]) / (s * real(compliance) + neck_impedance.recip());
    let mechanical = s * s
        + real(2.0 * spec.zetas[0] * spec.omegas[0]) * s
        + real(spec.omegas[0].powi(2));
    let displacement = real(spec.drive[0])
        / (mechanical - real(spec.areas[0]) * pressure_per_displacement);
    let plate_flow = real(spec.areas[0]) * s * displacement;
    let neck_outflow = pressure_per_displacement * displacement / neck_impedance;
    let observation = real(spec.gas.density / (2.0 * core::f64::consts::PI * distance)) * s;
    (
        observation * plate_flow,
        observation * (plate_flow + neck_outflow),
    )
}

fn measured_response(samples: &[f64], frequency: f64, rate: u32) -> C64 {
    // Whole periods in the last second; first two seconds let transients decay.
    let mut amplitude = C64::ZERO;
    let start = samples.len() - rate as usize;
    for (index, &sample) in samples.iter().enumerate().skip(start) {
        let phase = 2.0 * core::f64::consts::PI * frequency * index as f64 / f64::from(rate);
        amplitude.re += sample * phase.cos();
        amplitude.im += sample * phase.sin();
    }
    amplitude * C64::from_re(2.0 / f64::from(rate))
}

#[test]
fn g1_compact_pressure_matches_independent_coupled_harmonic_equations() {
    let spec = panel();
    let rate = 16_000;
    let distance = 2.0;
    for frequency in [20.0, 80.0, 160.0] {
        let forces: Vec<_> = (0..3 * rate)
            .map(|i| {
                (2.0 * core::f64::consts::PI * frequency * f64::from(i) / f64::from(rate)).cos()
            })
            .collect();
        let samples = realize_plate_cavity_compact_with_observer(
            &spec,
            &forces,
            rate,
            distance,
            GasAbsorptionModel::Lossless,
            None,
        )
        .unwrap();
        let observed = measured_response(&samples, frequency, rate);
        let (plate, expected) = harmonic_pressure(&spec, frequency, distance);
        let relative = (observed - expected).abs() / expected.abs();
        assert!(
            relative < 0.003,
            "{frequency} Hz: observed {observed:?}, oracle {expected:?}, relative {relative}"
        );
        if frequency == 20.0 {
            assert!(
                expected.abs() < 0.1 * plate.abs(),
                "outward plate and neck inflow must cancel below resonance"
            );
        }
    }
}

#[test]
fn g3_density_and_mechanical_mass_similarity_preserves_flow_and_scales_pressure() {
    let mut base = panel();
    // An aperture with zero cylindrical length removes the Bessel wall loss.
    // Remaining acoustic inertia and radiation resistance both scale with rho.
    base.neck_length_m = 0.0;
    let mut dense = base.clone();
    dense.gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 202_650.0).unwrap();
    let ratio = dense.gas.density / base.gas.density;
    dense.areas[0] /= ratio.sqrt();
    dense.drive[0] /= ratio.sqrt();
    // Scale physical plate mass, stiffness, damping and applied force together.
    // Its mass-normalized displacement scales by sqrt(ratio), leaving flow fixed.
    let forces: Vec<_> = (0..1_000)
        .map(|i| {
            if i < 32 {
                (0.1 * f64::from(i)).sin()
            } else {
                0.0
            }
        })
        .collect();
    let dense_forces: Vec<_> = forces.iter().map(|force| ratio * force).collect();
    let a = realize_plate_cavity_compact_with_observer(
        &base,
        &forces,
        8_000,
        1.5,
        GasAbsorptionModel::Lossless,
        None,
    )
    .unwrap();
    let b = realize_plate_cavity_compact_with_observer(
        &dense,
        &dense_forces,
        8_000,
        1.5,
        GasAbsorptionModel::Lossless,
        None,
    )
    .unwrap();
    let scale = b.iter().map(|v| v.abs()).fold(0.0, f64::max);
    assert!(scale > 0.0);
    for (first, second) in a.iter().zip(b) {
        assert!((second - ratio * first).abs() < 2.0e-10 * scale);
    }
}

#[test]
fn prepared_plate_observer_preserves_original_exact_zoh_trajectory() {
    let spec = panel();
    let geometry = geometry(&spec);
    let plate = plate_force_and_flow(&spec.omegas, &spec.zetas, &spec.drive, &spec.areas).unwrap();
    let cavity = helmholtz_resonator_flow(
        spec.volume_m3,
        spec.neck_radius_m,
        spec.neck_length_m,
        spec.gas.density,
        spec.gas.sound_speed,
        crate::cavity_modal::resonance_neck_resistance(&geometry, &spec.gas).unwrap(),
        crate::cavity_modal::neck_ends(&geometry),
    )
    .unwrap();
    let system = transformer(plate, cavity, 1, 0, 1.0).unwrap();
    let split = split_quadratic(
        &system,
        LinearSplitBudget {
            max_modes: system.state_dim(),
            max_omega_dt: core::f64::consts::PI,
        },
    )
    .unwrap();
    let forces: Vec<_> = (0..24).map(|i| (0.7 * f64::from(i)).cos()).collect();
    let pressure = realize_plate_cavity_with_observer(
        &spec,
        &forces,
        8_000,
        1.0,
        GasAbsorptionModel::Lossless,
        None,
    )
    .unwrap();
    let mut state = vec![0.0; system.state_dim()];
    for (force, observed) in forces.into_iter().zip(pressure) {
        let step = step_split(
            &system,
            &split,
            LinearPropagator::ExactZoh,
            &state,
            &[force],
            1.0 / 8_000.0,
        )
        .unwrap();
        let expected = spec.gas.density * spec.areas[0] * (step.x[1] - state[1]) * 8_000.0
            / (2.0 * core::f64::consts::PI);
        assert!((observed - expected).abs() < 1.0e-11 * expected.abs().max(1.0e-10));
        state = step.x;
    }
}

#[test]
fn compact_observation_keeps_causality_and_refuses_incompatible_baffle() {
    let mut spec = panel();
    let mut force = vec![0.0; 256];
    force[100] = 1.0;
    let first = realize_plate_cavity_compact_with_observer(
        &spec,
        &force,
        8_000,
        1.0,
        GasAbsorptionModel::Lossless,
        None,
    )
    .unwrap();
    assert!(first[..100].iter().all(|x| *x == 0.0));
    assert!(first[100..].iter().any(|x| *x != 0.0));
    assert_eq!(
        first,
        realize_plate_cavity_compact_with_observer(
            &spec,
            &force,
            8_000,
            1.0,
            GasAbsorptionModel::Lossless,
            None,
        )
        .unwrap()
    );
    spec.outer_end = fs_scenario::NeckEnd::Unflanged;
    assert!(matches!(
        realize_plate_cavity_compact_with_observer(
            &spec,
            &force,
            8_000,
            1.0,
            GasAbsorptionModel::Lossless,
            None,
        ),
        Err(CavityPhsError::Invalid { what })
            if what == "compact plate and neck observation requires a common exterior baffle"
    ));
    assert!(
        realize_plate_cavity_with_observer(
            &spec,
            &force,
            8_000,
            1.0,
            GasAbsorptionModel::Lossless,
            None,
        )
        .is_ok(),
        "the plate-only observer retains its previous geometry domain"
    );
}
