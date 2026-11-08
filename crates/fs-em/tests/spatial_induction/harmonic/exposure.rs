//! Independent finite-time heat oracles on the original solved 3D field.
use super::*;
use fs_em::spatial::vector_potential::induction::harmonic::{
    SpatialHarmonicInductionCell,
    exposure::{SpatialHarmonicHeatCell, SpatialHarmonicHeatExposure, SpatialHarmonicHeatLimits},
};

fn limits() -> SpatialHarmonicHeatLimits {
    SpatialHarmonicHeatLimits {
        maximum_cells: 6,
        maximum_output_bytes: 4096,
    }
}

fn relative(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual - expected).abs() <= allowance * expected.abs(),
        "actual={actual:e} expected={expected:e} relative={:e}",
        (actual - expected).abs() / expected.abs()
    );
}

fn analytic_square(z: C64, omega: f64, phase: f64, duration: f64) -> f64 {
    // Independent direct antiderivatives, used only away from cancellation.
    let end = phase + omega * duration;
    let cosine = (2.0 * end).sin() - (2.0 * phase).sin();
    let sine = (2.0 * phase).cos() - (2.0 * end).cos();
    (z.re * z.re + z.im * z.im) * duration / 2.0
        + (z.re * z.re - z.im * z.im) * cosine / (4.0 * omega)
        + z.re * z.im * sine / (2.0 * omega)
}

fn integrated_physical_cell(
    cell: &SpatialHarmonicInductionCell,
    omega: f64,
    phase: f64,
    duration: f64,
) -> f64 {
    // Independent positive four-point degree-two tetrahedral quadrature,
    // followed by Simpson integration of real E(t).J(t). Neither the native
    // P1 mass formula nor the analytic harmonic moment kernel is called.
    let high = (5.0 + 3.0 * 5.0_f64.sqrt()) / 20.0;
    let low = (5.0 - 5.0_f64.sqrt()) / 20.0;
    let intervals = 2048;
    let mut integral = 0.0;
    for point in 0..4 {
        let shape: [f64; 4] = std::array::from_fn(|i| if i == point { high } else { low });
        let electric: [C64; 3] = std::array::from_fn(|axis| {
            (0..4).fold(C64::ZERO, |sum, i| {
                sum + cell.conductor_electric_field_v_m[i][axis].scale(shape[i])
            })
        });
        let current: [C64; 3] = std::array::from_fn(|axis| {
            (0..4).fold(C64::ZERO, |sum, i| {
                sum + cell.current_density_a_m2[i][axis].scale(shape[i])
            })
        });
        let mut sum = 0.0;
        for step in 0..=intervals {
            let time = duration * f64::from(step) / f64::from(intervals);
            let (sine, cosine) = (phase + omega * time).sin_cos();
            let physical = |z: C64| z.re * cosine + z.im * sine;
            let density: f64 = (0..3)
                .map(|axis| physical(electric[axis]) * physical(current[axis]))
                .sum();
            let weight = if step == 0 || step == intervals {
                1.0
            } else if step % 2 == 0 {
                2.0
            } else {
                4.0
            };
            sum += weight * density;
        }
        integral += sum * (duration / (3.0 * f64::from(intervals))) * (cell.volume_m3 / 4.0);
    }
    integral
}

#[test]
fn finite_interval_heat_matches_closed_cube_mode_and_independent_spacetime_quadrature() {
    let (problem, sources) = cube(1, None);
    let initial = problem
        .solve(
            &[[0.0; 3]; 6],
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let (m, k, omega) = (3.0 * 5.0 / 40.0, 2.0, 3.0);
    let drive_phase = C64::new(0.8, -0.6);
    for alpha in [0.0, 0.35] {
        let velocity: Vec<_> = problem
            .positions_m()
            .iter()
            .map(|p| p.map(|x| alpha * x))
            .collect();
        let response = state
            .solve_harmonic(
                &[[drive_phase.scale(4.0); 3]; 6],
                &velocity,
                omega,
                harmonic_controls(),
                &mut || false,
            )
            .unwrap();
        let amplitude = drive_phase.scale(k) / C64::new(k + 2.0 * alpha * m, -omega * m);
        let conductor_rate = C64::new(2.0 * alpha, -omega) * amplitude;
        for (phase, duration) in [(-0.7, 0.13), (0.2, 0.37), (1.4, 1.7)] {
            let heat = response
                .observe_joule_heat(duration, phase, limits(), &mut || false)
                .unwrap();
            let expected = m * analytic_square(conductor_rate, omega, phase, duration);
            relative(heat.heat_j(), expected, 3e-9);
            relative(heat.mean_power_w(), expected / duration, 3e-9);
            assert!(core::ptr::eq(heat.response(), &response));
            assert_eq!(heat.duration_s(), duration);
            assert_eq!(heat.phase_at_start_rad(), phase);
            assert_eq!(heat.cells().len(), response.cells().len());
            for (observed, original) in heat.cells().iter().zip(response.cells()) {
                assert_eq!(observed.cell, original.cell);
                assert_eq!(observed.component, original.component);
                assert_eq!(observed.volume_m3, original.volume_m3);
                assert!(core::ptr::eq(
                    observed.source,
                    state.conductors()[observed.cell].as_ref().unwrap(),
                ));
                let independently_integrated =
                    integrated_physical_cell(original, omega, phase, duration);
                relative(observed.heat_j, independently_integrated, 3e-10);
                assert!(observed.heat_j >= 0.0);
            }
            near(
                heat.cells().iter().map(|cell| cell.heat_j).sum(),
                heat.heat_j(),
                1e-13,
            );
        }
        let cycles = response
            .observe_joule_heat(
                8.0 * std::f64::consts::TAU / omega,
                -0.3,
                limits(),
                &mut || false,
            )
            .unwrap();
        relative(
            cycles.mean_power_w(),
            response.balance().mean_ohmic_power_w,
            2e-13,
        );
        if alpha > 0.0 {
            assert!(cycles.mean_power_w() > response.balance().mean_supplied_power_w);
        }
    }
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb(), initial.potential_wb());
}

#[test]
fn tiny_windows_retain_positive_heat_at_a_field_zero() {
    let (problem, sources) = cube(1, None);
    let initial = problem
        .solve(
            &[[0.0; 3]; 6],
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    // K=2, M=3/8 and omega=4: this source gives a real A amplitude of one.
    // The physical conducting E begins at zero and grows as omega*sin(omega*t).
    let omega = 4.0;
    let response = state
        .solve_harmonic(
            &[[C64::new(4.0, -3.0); 3]; 6],
            &[[0.0; 3]; 8],
            omega,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    for duration in [1e-6, 1e-9, 1e-12] {
        let heat = response
            .observe_joule_heat(duration, 0.0, limits(), &mut || false)
            .unwrap();
        let expected: f64 = response
            .cells()
            .iter()
            .map(|cell| integrated_physical_cell(cell, omega, 0.0, duration))
            .sum();
        assert!(heat.heat_j() > 0.0);
        relative(heat.heat_j(), expected, 2e-11);
        // A rounded 1-sinc(2*omega*T) already vanishes at T=1e-9.
        // Below that scale, compare to the original field quadrature so the
        // oracle includes the solve's finite numerical residue at the zero.
        if duration >= 1e-9 {
            let leading = (3.0 / 8.0) * omega.powi(4) * duration.powi(3) / 3.0;
            relative(heat.heat_j(), leading, 3e-6);
        }
        assert!(heat.heat_j() < response.balance().mean_ohmic_power_w * duration * 1e-9);
    }
}

#[test]
fn phase_covariance_interval_composition_and_zero_field_preserve_original_response() {
    let (problem, sources) = cube(1, None);
    let initial = problem
        .solve(
            &[[0.0; 3]; 6],
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let omega = 3.0;
    let response = state
        .solve_harmonic(
            &[[C64::new(4.0, 0.0); 3]; 6],
            &[[0.0; 3]; 8],
            omega,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let (phase, duration, first) = (0.4, 0.73, 0.21);
    let whole = response
        .observe_joule_heat(duration, phase, limits(), &mut || false)
        .unwrap();
    let left = response
        .observe_joule_heat(first, phase, limits(), &mut || false)
        .unwrap();
    let right = response
        .observe_joule_heat(
            duration - first,
            phase + omega * first,
            limits(),
            &mut || false,
        )
        .unwrap();
    relative(left.heat_j() + right.heat_j(), whole.heat_j(), 2e-14);
    let shift: f64 = -0.9;
    let rotated = state
        .solve_harmonic(
            &[[C64::new(shift.cos(), shift.sin()).scale(4.0); 3]; 6],
            &[[0.0; 3]; 8],
            omega,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let same_physical = rotated
        .observe_joule_heat(duration, phase + shift, limits(), &mut || false)
        .unwrap();
    relative(same_physical.heat_j(), whole.heat_j(), 3e-10);
    let repeated = response
        .observe_joule_heat(duration, phase, limits(), &mut || false)
        .unwrap();
    assert_eq!(repeated.heat_j().to_bits(), whole.heat_j().to_bits());
    assert!(core::ptr::eq(whole.response(), repeated.response()));
    let zero = state
        .solve_harmonic(
            &[[C64::ZERO; 3]; 6],
            &[[0.0; 3]; 8],
            omega,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let heat = zero
        .observe_joule_heat(1e-200, 1e200, limits(), &mut || false)
        .unwrap();
    assert_eq!(heat.heat_j(), 0.0);
    assert_eq!(heat.mean_power_w(), 0.0);
}

#[test]
fn heat_bounds_invalid_windows_and_late_cancellation_leave_the_field_reusable() {
    let (problem, sources) = cube(1, None);
    let initial = problem
        .solve(
            &[[0.0; 3]; 6],
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let response = state
        .solve_harmonic(
            &[[C64::new(4.0, 0.0); 3]; 6],
            &[[0.0; 3]; 8],
            3.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let bytes = core::mem::size_of::<SpatialHarmonicHeatExposure<'_, '_, '_>>()
        + 6 * core::mem::size_of::<SpatialHarmonicHeatCell<'_>>();
    let exact = SpatialHarmonicHeatLimits {
        maximum_output_bytes: bytes,
        ..limits()
    };
    let baseline = response
        .observe_joule_heat(0.4, 0.2, exact, &mut || false)
        .unwrap();
    for settings in [
        SpatialHarmonicHeatLimits {
            maximum_cells: 5,
            ..exact
        },
        SpatialHarmonicHeatLimits {
            maximum_output_bytes: bytes - 1,
            ..exact
        },
    ] {
        assert!(matches!(
            response.observe_joule_heat(0.4, 0.2, settings, &mut || false),
            Err(SpatialMagneticError::Invalid(_))
        ));
    }
    for duration in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            response
                .observe_joule_heat(duration, 0.0, exact, &mut || false)
                .is_err()
        );
    }
    for phase in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            response
                .observe_joule_heat(0.4, phase, exact, &mut || false)
                .is_err()
        );
    }
    for stop_at in [1, 5, 8] {
        let mut polls = 0;
        assert!(matches!(
            response.observe_joule_heat(0.4, 0.2, exact, &mut || {
                polls += 1;
                polls >= stop_at
            }),
            Err(SpatialMagneticError::Cancelled)
        ));
    }
    let repeated = response
        .observe_joule_heat(0.4, 0.2, exact, &mut || false)
        .unwrap();
    assert_eq!(baseline.heat_j().to_bits(), repeated.heat_j().to_bits());
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb(), initial.potential_wb());
}
