//! G1/G3: original linearized-Euler equations, boundary data and laboratory clock.
use super::*;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.3, 340.0).unwrap()
}

fn field(wave: &HalfPlaneDiffraction, point: [f64; 3]) -> [C64; 4] {
    wave.amplitudes(point, C64::new(0.4, -0.7), || false)
        .unwrap()
}

fn derivative(wave: &HalfPlaneDiffraction, point: [f64; 3], axis: usize) -> [C64; 4] {
    let h = 2e-5;
    let samples: [[C64; 4]; 4] = [-2.0, -1.0, 1.0, 2.0].map(|offset| {
        let mut at = point;
        at[axis] += offset * h;
        field(wave, at)
    });
    core::array::from_fn(|i| {
        (samples[0][i] - samples[1][i].scale(8.0) + samples[2][i].scale(8.0)
            - samples[3][i])
            .scale(1.0 / (12.0 * h))
    })
}

#[test]
fn convected_diffraction_satisfies_the_original_pressure_and_momentum_equations() {
    let rho = medium().density_kg_m3();
    let bulk = medium().adiabatic_bulk_modulus_pa();
    let frequency = 380.0;
    let minus_i_omega = C64::new(0.0, -core::f64::consts::TAU * frequency);
    for flow in [-180.0, 0.0, 120.0] {
        for boundary in [IdealPlaneBoundary::Rigid, IdealPlaneBoundary::PressureRelease] {
            let wave = HalfPlaneDiffraction::with_edge_parallel_flow(
                medium(), frequency, [-0.48, -0.64, 0.6], boundary, flow,
            ).unwrap();
            for point in [[0.3, 0.2, 0.1], [0.4, -0.2, -0.3], [-0.15, -0.25, 0.2]] {
                let value = field(&wave, point);
                let derivatives: [[C64; 4]; 3] = core::array::from_fn(|axis| derivative(&wave, point, axis));
                // Dp + K div(u) = 0 and rho D(u) + grad(p) = 0,
                // using the original LABORATORY frequency, not the stored Omega.
                let material_p = minus_i_omega * value[0] + derivatives[2][0].scale(flow);
                let divergence = derivatives[0][1] + derivatives[1][2] + derivatives[2][3];
                assert!((material_p + divergence.scale(bulk)).abs() < 2e-6);
                for axis in 0..3 {
                    let material_u = minus_i_omega * value[axis + 1]
                        + derivatives[2][axis + 1].scale(flow);
                    assert!((material_u.scale(rho) + derivatives[axis][0]).abs() < 2e-8);
                }
            }
        }
    }
}

#[test]
fn flow_retains_the_incident_dispersion_lab_clock_and_both_sheet_faces() {
    let frequency = 380.0;
    for flow in [-180.0, 120.0] {
        for boundary in [IdealPlaneBoundary::Rigid, IdealPlaneBoundary::PressureRelease] {
            let wave = HalfPlaneDiffraction::with_edge_parallel_flow(
                medium(), frequency, [-0.48, -0.64, 0.6], boundary, flow,
            ).unwrap();
            let expected_k = core::f64::consts::TAU * frequency / (340.0 + flow * 0.6);
            assert!((wave.wavenumber_rad_m() / expected_k - 1.0).abs() < 3e-16);
            assert_eq!(wave.laboratory_frequency_hz(), frequency);
            assert_eq!(wave.mean_flow_velocity_m_s(), [0.0, 0.0, flow]);
            assert!((wave.intrinsic_angular_frequency_rad_s()
                + flow * expected_k * 0.6 - core::f64::consts::TAU * frequency).abs() < 1e-12);
            for y in [0.0, -0.0] {
                let line = wave.evaluate_harmonic_line([0.3, y, 0.2], C64::ONE, || false).unwrap();
                assert_eq!(line.frequency_hz, frequency);
                match boundary {
                    IdealPlaneBoundary::Rigid => assert_eq!(line.particle_velocity_m_s[1], C64::ZERO),
                    IdealPlaneBoundary::PressureRelease => assert_eq!(line.pressure_pa, C64::ZERO),
                }
            }
        }
    }
}

#[test]
fn rest_limit_is_identical_and_transverse_incidence_still_advects_energy() {
    let direction = [-0.6, -0.8, 0.0];
    let point = [0.1, -0.2, 0.3];
    let amplitude = C64::new(0.4, -0.7);
    let stationary = HalfPlaneDiffraction::new(medium(), 380.0, direction, IdealPlaneBoundary::Rigid).unwrap();
    for flow in [0.0, 120.0] {
        let moving = HalfPlaneDiffraction::with_edge_parallel_flow(
            medium(), 380.0, direction, IdealPlaneBoundary::Rigid, flow,
        ).unwrap();
        // kz=0: the entire p/u solution is unchanged even with nonzero flow.
        assert_eq!(field(&stationary, point), field(&moving, point));
        let line = moving.evaluate_harmonic_line(point, amplitude, || false).unwrap();
        let result = moving.evaluate_energy_transport(point, amplitude, || false).unwrap();
        let expected_energy = line.pressure_pa.norm_sq() / (4.0 * medium().adiabatic_bulk_modulus_pa())
            + 0.25 * medium().density_kg_m3()
                * line.particle_velocity_m_s.iter().map(|v| v.norm_sq()).sum::<f64>();
        assert!((result.energy_density_j_m3 / expected_energy - 1.0).abs() < 4e-16);
        assert_eq!(result.pressure_work_flux_w_m2[2], 0.0);
        assert!((result.energy_flux_w_m2[2] - flow * expected_energy).abs() < 1e-18);
    }
}

#[test]
fn unresolved_flow_and_cancelled_observations_refuse() {
    for speed in [f64::NAN, f64::INFINITY, -340.0, 340.0, 400.0] {
        assert!(HalfPlaneDiffraction::with_edge_parallel_flow(
            medium(), 380.0, [-0.48, -0.64, 0.6], IdealPlaneBoundary::Rigid, speed,
        ).is_err());
    }
    let wave = HalfPlaneDiffraction::with_edge_parallel_flow(
        medium(), 380.0, [-0.48, -0.64, 0.6], IdealPlaneBoundary::Rigid, 120.0,
    ).unwrap();
    assert!(matches!(wave.evaluate_energy_transport([0.1, 0.2, 0.3], C64::ONE, || true),
        Err(AeroacError::Cancelled { .. })));
}

#[test]
fn coherent_flow_spectrum_conserves_transport_including_interference_advection() {
    let flow = 120.0;
    let frame = HalfPlaneFrame::new([0.0; 3], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).unwrap();
    let waves = [
        HalfPlaneIncidentWave {
            propagation_direction: [-0.48, -0.64, 0.6],
            pressure_pa: C64::new(0.4, -0.7),
        },
        HalfPlaneIncidentWave {
            propagation_direction: [0.36, -0.48, -0.8],
            pressure_pa: C64::new(0.3, 0.2),
        },
    ];
    let spectrum = PreparedHalfPlaneSpectrum::with_edge_parallel_flow(
        medium(), 380.0, frame, IdealPlaneBoundary::Rigid, flow, [0.0; 3],
        &waves, HalfPlaneSpectrumBudget { max_components: 2 }, || false,
    ).unwrap();
    let point = [0.3, 0.2, 0.1];
    let h = 2e-5;
    let mut transported_divergence = 0.0;
    let mut pressure_work_divergence = 0.0;
    for axis in 0..3 {
        let mut lo = point;
        let mut hi = point;
        lo[axis] -= h;
        hi[axis] += h;
        let a = spectrum.evaluate_energy_transport(HalfPlaneReceiver::WorldPosition(lo), || false).unwrap();
        let b = spectrum.evaluate_energy_transport(HalfPlaneReceiver::WorldPosition(hi), || false).unwrap();
        transported_divergence += (b.energy_flux_w_m2[axis] - a.energy_flux_w_m2[axis]) / (2.0 * h);
        pressure_work_divergence += (b.pressure_work_flux_w_m2[axis] - a.pressure_work_flux_w_m2[axis]) / (2.0 * h);
    }
    // The same laboratory frequency has DIFFERENT kz/Omega for each component.
    // Their interference energy varies along the flowing edge; omitting U*E
    // breaks the original steady linearized-Euler energy equation.
    assert!(transported_divergence.abs() < 2e-9, "{transported_divergence:e}");
    assert!(pressure_work_divergence.abs() > 1e-4);
    let line = spectrum.evaluate_harmonic_line(HalfPlaneReceiver::WorldPosition(point), || false).unwrap();
    assert_eq!(line.frequency_hz, 380.0);
}

#[test]
fn placed_flow_spectrum_retains_each_wave_number_and_absolute_source_phase() {
    let frequency = 380.0;
    let flow = 120.0;
    let origin = [1.25, -2.0, 0.5];
    let reference = [-0.125, 0.5, 0.25];
    let frame = HalfPlaneFrame::new(origin, [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]).unwrap();
    let local_point = [0.25, -0.375, 0.5];
    let world_point = [origin[0] + local_point[2], origin[1] + local_point[0], origin[2] + local_point[1]];
    let local_directions = [[-0.48, -0.64, 0.6], [0.36, -0.48, -0.8]];
    let amplitudes = [C64::new(0.4, -0.7), C64::new(0.3, 0.2)];
    let waves: [HalfPlaneIncidentWave; 2] = core::array::from_fn(|i| {
        let [x, y, z] = local_directions[i];
        HalfPlaneIncidentWave { propagation_direction: [z, x, y], pressure_pa: amplitudes[i] }
    });
    let spectrum = PreparedHalfPlaneSpectrum::with_edge_parallel_flow(
        medium(), frequency, frame, IdealPlaneBoundary::Rigid, flow, reference,
        &waves, HalfPlaneSpectrumBudget { max_components: 2 }, || false,
    ).unwrap();
    let observed = spectrum.evaluate_harmonic_line(HalfPlaneReceiver::WorldPosition(world_point), || false).unwrap();
    let mut expected_p = C64::ZERO;
    let mut expected_u = [C64::ZERO; 3];
    for i in 0..2 {
        let direction = local_directions[i];
        let intrinsic_hz = frequency / (1.0 + flow * direction[2] / 340.0);
        let k = core::f64::consts::TAU * frequency / (340.0 + flow * direction[2]);
        let phase: f64 = k * (0..3).map(|j| waves[i].propagation_direction[j] * (origin[j] - reference[j])).sum::<f64>();
        let amplitude = amplitudes[i] * C64::new(phase.cos(), phase.sin());
        let original = HalfPlaneDiffraction::new(medium(), intrinsic_hz, direction, IdealPlaneBoundary::Rigid).unwrap();
        let field = original.evaluate_harmonic_line(local_point, amplitude, || false).unwrap();
        expected_p = expected_p + field.pressure_pa;
        for (world, local) in [(0, 2), (1, 0), (2, 1)] {
            expected_u[world] = expected_u[world] + field.particle_velocity_m_s[local];
        }
    }
    assert!((observed.pressure_pa - expected_p).abs() < 1e-13);
    for axis in 0..3 {
        assert!((observed.particle_velocity_m_s[axis] - expected_u[axis]).abs() < 1e-15);
    }
    assert_eq!(observed.frequency_hz, frequency);
    assert_eq!(spectrum.mean_flow_velocity_m_s(), [flow, 0.0, 0.0]);
}
