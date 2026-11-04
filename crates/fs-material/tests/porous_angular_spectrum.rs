//! G0/G3 coherent refraction, interference and source-retention controls.
use core::cell::Cell;
use fs_material::MaterialError;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::layers::{
    AcousticSpectrumError, AngularSpectrumBudget, IncidentPlaneWave, Layer,
    MAX_ANGULAR_SPECTRUM_COMPONENTS, TransmittingLayers,
};
use fs_math::c64::C64;

const OMEGA: f64 = core::f64::consts::TAU * 1000.0;

fn fluid() -> ThermoviscousFluidState {
    // Explicit analytic fixture; no measured material claim.
    ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::try_new(1.2, 340.0).unwrap(),
        FluidTransportState::try_new(1.2, 1005.0, 0.0257, 1.81e-5).unwrap(),
        1.4,
    )
    .unwrap()
}
fn budget(max_components: usize) -> AngularSpectrumBudget {
    AngularSpectrumBudget { max_components }
}
fn phase(angle: f64) -> C64 {
    C64::new(angle.cos(), angle.sin())
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}"
    );
}
fn near_real(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn single_component_matches_independent_snell_fields_and_interface_work() {
    let incident = fluid();
    let receiver = FluidAcousticState::try_new(2.4, 510.0).unwrap();
    let stack = TransmittingLayers::new(&incident, &receiver, OMEGA, &[]).unwrap();
    let amplitude = C64::new(0.7, -0.3);
    let direction = [0.36, 0.48, 0.8];
    let beam = stack
        .prepare_angular_spectrum(
            &[IncidentPlaneWave {
                direction,
                pressure_pa: amplitude,
            }],
            budget(1),
            &mut || false,
        )
        .unwrap();
    let zi = 1.2 * 340.0;
    let zt = 2.4 * 510.0;
    let sin_t = (510.0 / 340.0) * 0.6;
    let cos_t = f64::sqrt(1.0 - sin_t * sin_t);
    let reflection = (zt * 0.8 - zi * cos_t) / (zt * 0.8 + zi * cos_t);
    let transmission = 1.0 + reflection;
    let position = [0.03, -0.015, 0.12];
    let argument = -OMEGA
        * ((direction[0] * position[0] + direction[1] * position[1]) / 340.0
            + cos_t * position[2] / 510.0);
    let p = amplitude.scale(transmission) * phase(argument);
    let sample = beam.transmitted_field(position, &mut || false).unwrap();
    near(sample.pressure_pa, p, 4e-15);
    for (actual, factor) in sample.particle_velocity_m_s.into_iter().zip([
        direction[0] / (340.0 * 2.4),
        direction[1] / (340.0 * 2.4),
        cos_t / zt,
    ]) {
        near(actual, p.scale(factor), 2e-17);
    }
    let entering = beam
        .incident_field([0.03, -0.015, 0.0], &mut || false)
        .unwrap();
    let exiting = beam
        .transmitted_field([0.03, -0.015, 0.0], &mut || false)
        .unwrap();
    near(entering.pressure_pa, exiting.pressure_pa, 3e-15);
    near(
        entering.particle_velocity_m_s[2],
        exiting.particle_velocity_m_s[2],
        2e-17,
    );
    let incident_flux = 0.5 * amplitude.norm_sq() * direction[2] / zi;
    near_real(
        exiting.active_intensity_w_m2[2],
        incident_flux * (1.0 - reflection * reflection),
        4e-18,
    );
    near_real(
        entering.active_intensity_w_m2[2],
        exiting.active_intensity_w_m2[2],
        4e-18,
    );
    assert!(
        sample
            .reactive_intensity_w_m2
            .iter()
            .all(|v| v.abs() < 2e-18)
    );
    assert!(core::ptr::eq(beam.stack(), &stack));
    assert_eq!(beam.component_count(), 1);
}

#[test]
fn matched_fluids_keep_entrance_and_exit_phase_origins() {
    let fluid = fluid();
    let direction = [0.6, 0.0, 0.8];
    let amplitude = C64::new(0.4, 0.9);
    let k = OMEGA / 340.0;
    for thickness in [0.0, 0.17] {
        let stack = TransmittingLayers::new(
            &fluid,
            &fluid,
            OMEGA,
            &[Layer::FluidGap {
                thickness_m: thickness,
            }],
        )
        .unwrap();
        let beam = stack
            .prepare_angular_spectrum(
                &[IncidentPlaneWave {
                    direction,
                    pressure_pa: amplitude,
                }],
                budget(1),
                &mut || false,
            )
            .unwrap();
        let exit_point = [0.04, 0.02, 0.2];
        let expected =
            amplitude * phase(-k * (0.6 * exit_point[0] + 0.8 * (exit_point[2] + thickness)));
        let transmitted = beam.transmitted_field(exit_point, &mut || false).unwrap();
        near(transmitted.pressure_pa, expected, 6e-15);
        for (actual, d) in transmitted.particle_velocity_m_s.into_iter().zip(direction) {
            near(actual, expected.scale(d / (1.2 * 340.0)), 3e-17);
        }
        let entry_point = [0.04, 0.02, -0.2];
        let incident = beam.incident_field(entry_point, &mut || false).unwrap();
        near(
            incident.pressure_pa,
            amplitude * phase(-k * (0.6 * entry_point[0] + 0.8 * entry_point[2])),
            6e-15,
        );
    }
}

#[test]
fn coherent_crossing_waves_keep_interference_and_reactive_transverse_flow() {
    let fluid = fluid();
    let stack = TransmittingLayers::new(&fluid, &fluid, OMEGA, &[]).unwrap();
    let components = [
        IncidentPlaneWave {
            direction: [0.6, 0.0, 0.8],
            pressure_pa: C64::from_re(0.5),
        },
        IncidentPlaneWave {
            direction: [-0.6, 0.0, 0.8],
            pressure_pa: C64::from_re(0.5),
        },
    ];
    let beam = stack
        .prepare_angular_spectrum(&components, budget(2), &mut || false)
        .unwrap();
    let z = 0.13;
    let k = OMEGA / 340.0;
    let impedance = 1.2 * 340.0;
    for angle in [
        0.0,
        core::f64::consts::FRAC_PI_4,
        core::f64::consts::FRAC_PI_2,
    ] {
        let x = angle / (k * 0.6);
        let carrier = phase(-k * 0.8 * z);
        let expected_p = carrier.scale(angle.cos());
        let expected_vx = carrier * C64::new(0.0, -0.6 / impedance * angle.sin());
        let sample = beam.transmitted_field([x, 0.0, z], &mut || false).unwrap();
        near(sample.pressure_pa, expected_p, 3e-15);
        near(sample.particle_velocity_m_s[0], expected_vx, 1e-17);
        near(
            sample.particle_velocity_m_s[2],
            expected_p.scale(0.8 / impedance),
            1e-17,
        );
        near_real(sample.active_intensity_w_m2[0], 0.0, 3e-18);
        near_real(
            sample.active_intensity_w_m2[2],
            0.5 * 0.8 / impedance * angle.cos().powi(2),
            4e-18,
        );
        near_real(
            sample.reactive_intensity_w_m2[0],
            0.5 * 0.6 / impedance * angle.cos() * angle.sin(),
            3e-18,
        );
        assert_eq!(
            sample,
            beam.transmitted_field([x, 0.0, z], &mut || false).unwrap()
        );
    }
    // At the central antinode coherent normal intensity is twice the sum of
    // the two isolated component intensities. Incoherent angular_power would
    // erase both this enhancement and the off-axis nodes.
    let center = beam
        .transmitted_field([0.0, 0.0, z], &mut || false)
        .unwrap();
    let incoherent = 2.0 * (0.5 * 0.5_f64.powi(2) * 0.8 / impedance);
    near_real(center.active_intensity_w_m2[2], 2.0 * incoherent, 4e-18);
}

#[test]
fn total_internal_reflection_retains_decay_and_reactive_normal_intensity() {
    let incident = fluid();
    let receiver = FluidAcousticState::try_new(0.3, 680.0).unwrap();
    let stack = TransmittingLayers::new(&incident, &receiver, OMEGA, &[]).unwrap();
    let beam = stack
        .prepare_angular_spectrum(
            &[IncidentPlaneWave {
                direction: [0.8, 0.0, 0.6],
                pressure_pa: C64::ONE,
            }],
            budget(1),
            &mut || false,
        )
        .unwrap();
    let kt = OMEGA / 340.0 * 0.8;
    let kr = OMEGA / 680.0;
    let alpha = (kt * kt - kr * kr).sqrt();
    let yi = 0.6 / (1.2 * 340.0);
    let yr = C64::new(0.0, -alpha / (0.3 * OMEGA));
    let transmission = C64::from_re(2.0 * yi) / (C64::from_re(yi) + yr);
    for depth in [0.0, 1.0 / alpha, 5.0 / alpha] {
        let expected = transmission.scale((-alpha * depth).exp());
        let sample = beam
            .transmitted_field([0.0, 0.0, depth], &mut || false)
            .unwrap();
        near(sample.pressure_pa, expected, 4e-15);
        near(sample.particle_velocity_m_s[2], expected * yr, 3e-17);
        near_real(sample.active_intensity_w_m2[2], 0.0, 3e-18);
        near_real(
            sample.reactive_intensity_w_m2[2],
            0.5 * alpha / (0.3 * OMEGA) * expected.norm_sq(),
            2e-17,
        );
        near_real(
            sample.active_intensity_w_m2[0],
            0.5 * kt / (0.3 * OMEGA) * expected.norm_sq(),
            2e-17,
        );
        assert!(sample.reactive_intensity_w_m2[2] > 0.0);
    }
}

#[test]
fn source_warm_air_layer_preserves_evanescent_tunnelling_amplitude_and_phase() {
    let model = GasSpec::dry_air_ussa1976();
    let incident = GasState::try_new(&model, 293.15, 101_325.0).unwrap();
    let warm = GasState::try_new(&model, 353.15, 101_325.0).unwrap();
    let cosine = 0.2_f64;
    let sine = (1.0 - cosine * cosine).sqrt();
    let kt = OMEGA / incident.sound_speed * sine;
    let kw = OMEGA / warm.sound_speed;
    let alpha = (kt * kt - kw * kw).sqrt();
    let eta = OMEGA * warm.density * cosine / (incident.density * incident.sound_speed * alpha);
    let mut previous = 1.0_f64;
    for thickness in [0.02, 0.08, 0.3] {
        let stack = TransmittingLayers::new(
            &incident,
            &incident,
            OMEGA,
            &[Layer::FluidLayer {
                fluid: &warm,
                thickness_m: thickness,
            }],
        )
        .unwrap();
        let beam = stack
            .prepare_angular_spectrum(
                &[IncidentPlaneWave {
                    direction: [sine, 0.0, cosine],
                    pressure_pa: C64::ONE,
                }],
                budget(1),
                &mut || false,
            )
            .unwrap();
        // Independent conventional transfer matrix, safe for these modest
        // optical depths: T = 2 / (2 cosh(a d) + i(eta-1/eta) sinh(a d)).
        let depth = alpha * thickness;
        let expected =
            C64::from_re(2.0) / C64::new(2.0 * depth.cosh(), (eta - eta.recip()) * depth.sinh());
        let sample = beam.transmitted_field([0.0; 3], &mut || false).unwrap();
        near(sample.pressure_pa, expected, 2e-14);
        assert!(sample.active_intensity_w_m2[2] > 0.0);
        let magnitude = sample.pressure_pa.abs();
        assert!(magnitude < previous);
        previous = magnitude;
    }
}

struct CountedFluid {
    state: ThermoviscousFluidState,
    samples: Cell<usize>,
}
impl FluidAcoustics for CountedFluid {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.samples.set(self.samples.get() + 1);
        Ok(self.state.acoustic())
    }
}
impl FluidTransport for CountedFluid {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.state.transport())
    }
}
impl ThermoviscousFluid for CountedFluid {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        self.samples.set(self.samples.get() + 1);
        Ok(self.state)
    }
}
struct MutableReceiver {
    state: Cell<FluidAcousticState>,
    samples: Cell<usize>,
}
impl FluidAcoustics for MutableReceiver {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.samples.set(self.samples.get() + 1);
        Ok(self.state.get())
    }
}

#[test]
fn beam_preparation_and_observation_share_original_frozen_sources() {
    let incident = CountedFluid {
        state: fluid(),
        samples: Cell::new(0),
    };
    let receiver = MutableReceiver {
        state: Cell::new(FluidAcousticState::try_new(1.8, 420.0).unwrap()),
        samples: Cell::new(0),
    };
    let stack = TransmittingLayers::new(
        &incident,
        &receiver,
        OMEGA,
        &[Layer::FluidLayer {
            fluid: &receiver,
            thickness_m: 0.04,
        }],
    )
    .unwrap();
    assert_eq!(incident.samples.get(), 1);
    assert_eq!(receiver.samples.get(), 1);
    let components = [IncidentPlaneWave {
        direction: [0.6, 0.0, 0.8],
        pressure_pa: C64::ONE,
    }; 3];
    let beam = stack
        .prepare_angular_spectrum(&components, budget(3), &mut || false)
        .unwrap();
    let first = beam
        .transmitted_field([0.02, 0.0, 0.13], &mut || false)
        .unwrap();
    receiver
        .state
        .set(FluidAcousticState::try_new(50.0, 900.0).unwrap());
    for _ in 0..10 {
        assert_eq!(
            beam.transmitted_field([0.02, 0.0, 0.13], &mut || false)
                .unwrap(),
            first
        );
        beam.incident_field([0.0, 0.0, -0.1], &mut || false)
            .unwrap();
    }
    let another = stack
        .prepare_angular_spectrum(&components, budget(3), &mut || false)
        .unwrap();
    assert_eq!(
        another
            .transmitted_field([0.02, 0.0, 0.13], &mut || false)
            .unwrap(),
        first
    );
    assert_eq!(incident.samples.get(), 1);
    assert_eq!(receiver.samples.get(), 1);
}

#[test]
fn explicit_budgets_domains_and_cancellation_refuse_without_partial_fields() {
    let fluid = fluid();
    let stack = TransmittingLayers::new(&fluid, &fluid, OMEGA, &[]).unwrap();
    let wave = IncidentPlaneWave {
        direction: [0.6, 0.0, 0.8],
        pressure_pa: C64::ONE,
    };
    let components = [wave; 3];
    for max in [0, 2, MAX_ANGULAR_SPECTRUM_COMPONENTS + 1] {
        assert!(matches!(
            stack.prepare_angular_spectrum(&components, budget(max), &mut || false),
            Err(AcousticSpectrumError::Invalid(_))
        ));
    }
    assert!(
        stack
            .prepare_angular_spectrum(&[], budget(1), &mut || false)
            .is_err()
    );
    for direction in [
        [0.0; 3],
        [0.6, 0.0, -0.8],
        [1.0, 0.0, 0.0],
        [0.1, 0.0, 0.8],
        [f64::NAN, 0.0, 1.0],
    ] {
        assert!(
            stack
                .prepare_angular_spectrum(
                    &[IncidentPlaneWave { direction, ..wave }],
                    budget(1),
                    &mut || false
                )
                .is_err()
        );
    }
    assert!(
        stack
            .prepare_angular_spectrum(
                &[IncidentPlaneWave {
                    pressure_pa: C64::new(f64::INFINITY, 0.0),
                    ..wave
                }],
                budget(1),
                &mut || false
            )
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        stack.prepare_angular_spectrum(&components, budget(3), &mut || {
            polls += 1;
            polls == 3
        }),
        Err(AcousticSpectrumError::Cancelled)
    ));
    let beam = stack
        .prepare_angular_spectrum(&components, budget(3), &mut || false)
        .unwrap();
    let expected = beam
        .transmitted_field([0.02, 0.0, 0.1], &mut || false)
        .unwrap();
    for stop in [1, 3, 5] {
        polls = 0;
        assert_eq!(
            beam.transmitted_field([0.02, 0.0, 0.1], &mut || {
                polls += 1;
                polls == stop
            }),
            Err(AcousticSpectrumError::Cancelled)
        );
        assert_eq!(
            beam.transmitted_field([0.02, 0.0, 0.1], &mut || false)
                .unwrap(),
            expected
        );
    }
    assert!(beam.incident_field([0.0, 0.0, 0.1], &mut || false).is_err());
    assert!(
        beam.transmitted_field([0.0, 0.0, -0.1], &mut || false)
            .is_err()
    );
    assert!(
        beam.transmitted_field([f64::NAN, 0.0, 0.1], &mut || false)
            .is_err()
    );
}
