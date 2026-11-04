//! G1/G3/G4: native loaded-mouth observation, coherent radiation and refusals.
use fs_aeroac::observation::{HarmonicSampleBudget, HarmonicSampleClock};
use fs_couple::acoustic_realize::duct_harmonic::{
    DuctHarmonicBudget, DuctHarmonicError, DuctOpeningFlows, DuctOpeningId, DuctOpeningReceiver,
    observe_duct_harmonics,
};
use fs_couple::pcm_wav::{
    baffled::CircularOutletReceiver,
    observation::{PressureRenderer, harmonic::HarmonicPressureRenderer},
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd};
use fs_orbit::HbOrbit;
use fs_phs::acoustic_outlet::{
    AcousticOutletFlow, AcousticTapFlow, acoustic_chain_with_opening_flows,
};
use fs_phs::{AcousticSection, AcousticTap, ChainOutlet, MouthFlange, PortHamiltonian};
use std::f64::consts::{PI, TAU};

const RADIUS: f64 = 0.012;
const FREQUENCY: f64 = 120.0;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.2, 340.0).unwrap()
}
fn budget() -> DuctHarmonicBudget {
    DuctHarmonicBudget {
        max_openings: 8,
        max_harmonics: 8,
        max_projection_terms: 100_000,
        max_kernel_terms: 100_000,
        max_workspace_bytes: 8 << 20,
        minimum_clearance_m: 0.05,
    }
}
fn receiver(opening: DuctOpeningId) -> DuctOpeningReceiver {
    DuctOpeningReceiver {
        opening,
        outward_normal: [0.0, 0.0, 1.0],
        geometry: CircularOutletReceiver {
            position_m: [0.0, 0.0, 0.6],
            radial_rings: 32,
            angular_points: 32,
            maximum_frequency_hz: 1000.0,
        },
    }
}
fn native(
    outlet: ChainOutlet,
    taps: &[AcousticTap],
) -> (PortHamiltonian, AcousticOutletFlow, Vec<AcousticTapFlow>) {
    acoustic_chain_with_opening_flows(
        &[AcousticSection {
            length: 0.24,
            radius: RADIUS,
            outlet_radius: RADIUS,
            cells: 4,
        }],
        1.2,
        340.0,
        outlet,
        1,
        taps,
        None,
        None,
    )
    .unwrap()
}
fn waveform(states: usize, harmonics: usize) -> HbOrbit {
    HbOrbit {
        omega: TAU * FREQUENCY,
        coeffs: vec![vec![C64::ZERO; harmonics + 1]; states],
        residual_trace: Vec::new(),
        residual: f64::NAN,
        aft_samples: 0,
        sampling_defect: f64::NAN,
    }
}
fn flows<'a>(bell: &'a AcousticOutletFlow, taps: &'a [AcousticTapFlow]) -> DuctOpeningFlows<'a> {
    DuctOpeningFlows {
        bell,
        bell_radius_m: RADIUS,
        taps,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}
fn phasor(angle: f64) -> C64 {
    C64::new(angle.cos(), angle.sin())
}
fn set_flow(orbit: &mut HbOrbit, row: &AcousticOutletFlow, dc: f64, peaks: &[C64]) {
    let (coordinate, &weight) = row
        .coefficients()
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.abs().total_cmp(&b.abs()))
        .unwrap();
    assert_ne!(weight, 0.0);
    orbit.coeffs[coordinate][0] = C64::from_re(dc / weight);
    for (h, &peak) in peaks.iter().enumerate() {
        orbit.coeffs[coordinate][h + 1] = peak.conj().scale(0.5 / weight);
    }
}

#[test]
fn loaded_bell_flow_not_inlet_flow_reaches_the_physical_pressure_renderer() {
    let (system, bell, taps) = native(ChainOutlet::Radiating(MouthFlange::Flanged), &[]);
    let mut orbit = waveform(system.state_dim(), 3);
    let peaks = [
        C64::new(2e-5, -3e-5),
        C64::new(-4e-6, 2e-6),
        C64::new(3e-9, 1e-9),
    ];
    set_flow(&mut orbit, &bell, 7e-6, &peaks);
    let initial: Vec<_> = (0..system.state_dim())
        .map(|i| orbit.sample(i, 0.0))
        .collect();
    let native_flow = bell.volume_flow_m3_s(&initial).unwrap();
    let inlet_flow = system.output(&initial)[0];
    assert!((native_flow - inlet_flow).abs() > 1e-6);
    // The original radiation inertia is retained in the native outlet row.
    // Replacing it by the bare last-cell rho*dx/A would give a different flow.
    let largest = bell
        .coefficients()
        .iter()
        .map(|v| v.abs())
        .fold(0.0_f64, f64::max);
    let bare_cell_mass = 1.2 * (0.24 / 4.0) / (PI * RADIUS * RADIUS);
    assert!(largest < 1.0 / bare_cell_mass);
    let original = orbit.coeffs.clone();
    let result = observe_duct_harmonics(
        &orbit,
        flows(&bell, &taps),
        &[receiver(DuctOpeningId::Bell)],
        medium(),
        budget(),
        || false,
    )
    .unwrap();
    assert_eq!(orbit.coeffs, original);
    near(result.openings[0].dc_volume_flow_m3_s, 7e-6, 1e-20);
    for (actual, expected) in result.openings[0].peak_volume_flow_m3_s.iter().zip(peaks) {
        near(actual.re, expected.re, 1e-20);
        near(actual.im, expected.im, 1e-20);
    }
    // Independent on-axis uniform-disk Rayleigh integral:
    // p/Q = rho*c/A * (exp(i*k*z) - exp(i*k*sqrt(z^2+a^2))).
    // This includes finite aperture and physical delay, without a sample filter.
    let z: f64 = 0.6;
    let edge = (z * z + RADIUS * RADIUS).sqrt();
    let expected: Vec<_> = peaks
        .iter()
        .enumerate()
        .map(|(h, &q)| {
            let k = (h + 1) as f64 * orbit.omega / 340.0;
            (phasor(k * z) - phasor(k * edge)).scale(1.2 * 340.0 / (PI * RADIUS * RADIUS)) * q
        })
        .collect();
    let scale = expected.iter().map(|p| p.abs()).sum::<f64>();
    let clock = HarmonicSampleClock::new(Dd::from_f64(0.013), 8000).unwrap();
    let mut renderer = HarmonicPressureRenderer::new(
        &result.field,
        clock,
        HarmonicSampleBudget {
            max_samples: 64,
            max_tone_samples: 192,
        },
    )
    .unwrap();
    let mut pressure = [0.0; 64];
    renderer.block(&mut pressure).unwrap();
    for (j, actual) in pressure.into_iter().enumerate() {
        let time = 0.013 + j as f64 / 8000.0;
        let reference = expected
            .iter()
            .enumerate()
            .map(|(h, p)| {
                let angle = (h + 1) as f64 * orbit.omega * time;
                p.re * angle.cos() + p.im * angle.sin()
            })
            .sum::<f64>();
        near(actual, reference, 1e-6 * scale);
    }
}

#[test]
fn active_tap_ids_relative_phase_cancellation_and_clock_rotation_survive() {
    let tap = |station, open_fraction| AcousticTap {
        station,
        neck_length: 0.002,
        neck_radius: 0.003,
        open_fraction,
    };
    // The closed entry remains in native IDs, while active stations are unsorted.
    let (system, bell, taps) = native(
        ChainOutlet::Rigid,
        &[tap(0.5, 0.0), tap(0.8, 1.0), tap(0.2, 1.0)],
    );
    assert_eq!(taps.iter().map(|t| t.tap_index).collect::<Vec<_>>(), [1, 2]);
    let receivers = [
        receiver(DuctOpeningId::Bell),
        receiver(DuctOpeningId::Tap(1)),
        receiver(DuctOpeningId::Tap(2)),
    ];
    let mut orbit = waveform(system.state_dim(), 1);
    let q = C64::new(1e-5, -2e-5);
    set_flow(&mut orbit, &taps[0].flow, 0.0, &[q]);
    set_flow(&mut orbit, &taps[1].flow, 0.0, &[-q]);
    let result = observe_duct_harmonics(
        &orbit,
        flows(&bell, &taps),
        &receivers,
        medium(),
        budget(),
        || false,
    )
    .unwrap();
    assert_eq!(
        result
            .openings
            .iter()
            .map(|v| v.opening)
            .collect::<Vec<_>>(),
        [
            DuctOpeningId::Bell,
            DuctOpeningId::Tap(1),
            DuctOpeningId::Tap(2)
        ]
    );
    assert_eq!(result.openings[1].station, Some(0.8));
    assert_eq!(result.openings[2].station, Some(0.2));
    for time in [0.0, 0.0013, 0.0037] {
        near(
            result
                .field
                .pressure_at_clock_dd(Dd::from_f64(time))
                .unwrap(),
            0.0,
            1e-14,
        );
    }
    set_flow(&mut orbit, &taps[1].flow, 0.0, &[q.scale(-0.3)]);
    let base = observe_duct_harmonics(
        &orbit,
        flows(&bell, &taps),
        &receivers,
        medium(),
        budget(),
        || false,
    )
    .unwrap();
    let shift = 0.47;
    for row in &mut orbit.coeffs {
        row[1] = row[1] * phasor(shift);
    }
    let rotated = observe_duct_harmonics(
        &orbit,
        flows(&bell, &taps),
        &receivers,
        medium(),
        budget(),
        || false,
    )
    .unwrap();
    for time in [0.0, 0.0013, 0.0037] {
        near(
            rotated
                .field
                .pressure_at_clock_dd(Dd::from_f64(time))
                .unwrap(),
            base.field
                .pressure_at_clock_dd(Dd::from_f64(time + shift / orbit.omega))
                .unwrap(),
            1e-12,
        );
    }
}

#[test]
fn lateral_tone_hole_and_bell_fields_rotate_in_one_physical_frame() {
    let (system, bell, taps) = native(
        ChainOutlet::Radiating(MouthFlange::Flanged),
        &[AcousticTap {
            station: 0.4,
            neck_length: 0.002,
            neck_radius: 0.003,
            open_fraction: 1.0,
        }],
    );
    let mut orbit = waveform(system.state_dim(), 1);
    set_flow(&mut orbit, &bell, 0.0, &[C64::new(1e-5, 2e-5)]);
    set_flow(&mut orbit, &taps[0].flow, 0.0, &[C64::new(3e-6, -7e-6)]);
    // A forward-facing bell and lateral hole, with one common receiver. The
    // hole is 0.1 m behind the bell in this explicitly supplied geometry.
    let mut receivers = [
        receiver(DuctOpeningId::Bell),
        receiver(DuctOpeningId::Tap(0)),
    ];
    receivers[0].geometry.position_m = [0.4, 0.2, 0.6];
    receivers[1].geometry.position_m = [0.4, 0.2, 0.7];
    receivers[1].outward_normal = [5.0, 0.0, 0.0];
    let original = observe_duct_harmonics(
        &orbit,
        flows(&bell, &taps),
        &receivers,
        medium(),
        budget(),
        || false,
    )
    .unwrap();
    assert_eq!(original.openings[1].outward_normal, [1.0, 0.0, 0.0]);
    // Rigidly rotate the entire physical configuration by +pi/2 about y.
    // Pressure is invariant; the velocity must undergo the same rotation.
    let rotate = |v: [f64; 3]| [v[2], v[1], -v[0]];
    for receiver in &mut receivers {
        receiver.geometry.position_m = rotate(receiver.geometry.position_m);
        receiver.outward_normal = rotate(receiver.outward_normal);
    }
    let rotated = observe_duct_harmonics(
        &orbit,
        flows(&bell, &taps),
        &receivers,
        medium(),
        budget(),
        || false,
    )
    .unwrap();
    for time in [0.0, 0.0013, 0.0037] {
        let a = original
            .field
            .field_at_clock_dd(Dd::from_f64(time))
            .unwrap();
        let b = rotated.field.field_at_clock_dd(Dd::from_f64(time)).unwrap();
        assert!(a.particle_velocity_m_s.iter().any(|u| u.abs() > 1e-8));
        near(b.pressure_pa, a.pressure_pa, 1e-12);
        for (actual, expected) in b
            .particle_velocity_m_s
            .into_iter()
            .zip(rotate(a.particle_velocity_m_s))
        {
            near(actual, expected, 1e-15);
        }
    }
}

#[test]
fn a_rigid_bell_and_closed_tap_retain_state_but_emit_no_harmonics() {
    let closed = AcousticTap {
        station: 0.4,
        neck_length: 0.002,
        neck_radius: 0.003,
        open_fraction: 0.0,
    };
    let (system, bell, taps) = native(ChainOutlet::Rigid, &[closed]);
    assert!(taps.is_empty());
    let mut orbit = waveform(system.state_dim(), 2);
    for (i, row) in orbit.coeffs.iter_mut().enumerate() {
        row[0] = C64::from_re((i + 1) as f64 * 1e-7);
        row[1] = C64::new(1e-4, 2e-4);
        row[2] = C64::new(-3e-5, 4e-5);
    }
    let original = orbit.coeffs.clone();
    let result = observe_duct_harmonics(
        &orbit,
        flows(&bell, &taps),
        &[receiver(DuctOpeningId::Bell)],
        medium(),
        budget(),
        || false,
    )
    .unwrap();
    assert_eq!(result.openings[0].dc_volume_flow_m3_s, 0.0);
    assert_eq!(result.openings[0].peak_volume_flow_m3_s, [C64::ZERO; 2]);
    assert_eq!(orbit.coeffs, original);
    for time in [0.0, 0.002, 0.017] {
        assert_eq!(
            result
                .field
                .pressure_at_clock_dd(Dd::from_f64(time))
                .unwrap(),
            0.0
        );
    }
}

#[test]
fn wrong_rows_receiver_band_budgets_and_cancellation_refuse_atomically() {
    let (system, bell, taps) = native(ChainOutlet::Radiating(MouthFlange::Flanged), &[]);
    let mut orbit = waveform(system.state_dim(), 2);
    let mut limits = budget();
    limits.max_projection_terms = 1;
    assert!(matches!(
        observe_duct_harmonics(
            &orbit,
            flows(&bell, &taps),
            &[receiver(DuctOpeningId::Bell)],
            medium(),
            limits,
            || false
        ),
        Err(DuctHarmonicError::Budget(_))
    ));
    limits = budget();
    limits.max_kernel_terms = 1;
    assert!(matches!(
        observe_duct_harmonics(
            &orbit,
            flows(&bell, &taps),
            &[receiver(DuctOpeningId::Bell)],
            medium(),
            limits,
            || false
        ),
        Err(DuctHarmonicError::Budget(_))
    ));
    limits = budget();
    limits.max_workspace_bytes = 1;
    assert!(matches!(
        observe_duct_harmonics(
            &orbit,
            flows(&bell, &taps),
            &[receiver(DuctOpeningId::Bell)],
            medium(),
            limits,
            || false
        ),
        Err(DuctHarmonicError::Budget(_))
    ));
    for at in [0, 2, 7, 15, 30] {
        let mut polls = 0;
        assert!(matches!(
            observe_duct_harmonics(
                &orbit,
                flows(&bell, &taps),
                &[receiver(DuctOpeningId::Bell)],
                medium(),
                budget(),
                || {
                    let stop = polls == at;
                    polls += 1;
                    stop
                }
            ),
            Err(DuctHarmonicError::Cancelled)
        ));
    }
    let mut low_band = receiver(DuctOpeningId::Bell);
    low_band.geometry.maximum_frequency_hz = 1.99 * FREQUENCY;
    assert!(matches!(
        observe_duct_harmonics(
            &orbit,
            flows(&bell, &taps),
            &[low_band],
            medium(),
            budget(),
            || false
        ),
        Err(DuctHarmonicError::Invalid(_))
    ));
    assert!(matches!(
        observe_duct_harmonics(
            &orbit,
            flows(&bell, &taps),
            &[receiver(DuctOpeningId::Tap(0))],
            medium(),
            budget(),
            || false
        ),
        Err(DuctHarmonicError::Invalid(_))
    ));
    for normal in [[0.0; 3], [f64::NAN, 0.0, 1.0], [0.0, 0.0, -1.0]] {
        let mut geometry = receiver(DuctOpeningId::Bell);
        geometry.outward_normal = normal;
        assert!(matches!(
            observe_duct_harmonics(
                &orbit,
                flows(&bell, &taps),
                &[geometry],
                medium(),
                budget(),
                || false
            ),
            Err(DuctHarmonicError::Invalid(_))
        ));
    }
    orbit.coeffs[0][0].im = 1.0;
    assert!(matches!(
        observe_duct_harmonics(
            &orbit,
            flows(&bell, &taps),
            &[receiver(DuctOpeningId::Bell)],
            medium(),
            budget(),
            || false
        ),
        Err(DuctHarmonicError::Invalid(_))
    ));
    orbit.coeffs[0][0].im = 0.0;
    orbit.coeffs[0].pop();
    assert!(matches!(
        observe_duct_harmonics(
            &orbit,
            flows(&bell, &taps),
            &[receiver(DuctOpeningId::Bell)],
            medium(),
            budget(),
            || false
        ),
        Err(DuctHarmonicError::Invalid(_))
    ));
}
