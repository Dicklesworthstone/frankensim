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
use fs_material::porous::{PorousCard, PorousFluidLaw};
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

// Independent unscaled cosine/sine transfer, used only at modest depths.
// Returns cos(k_n d), sin(k_n d)/k_n and k_n sin(k_n d), including k_n=0.
fn normal_matrix(kn_squared: f64, depth: f64) -> [f64; 3] {
    if kn_squared > 0.0 {
        let kn = kn_squared.sqrt();
        let phase = kn * depth;
        [phase.cos(), phase.sin() / kn, kn * phase.sin()]
    } else if kn_squared < 0.0 {
        let alpha = (-kn_squared).sqrt();
        let decay = alpha * depth;
        [decay.cosh(), decay.sinh() / alpha, -alpha * decay.sinh()]
    } else {
        [1.0, depth, 0.0]
    }
}

#[test]
fn interior_refraction_matches_independent_fields_across_critical_incidence() {
    let incident = fluid();
    let layer = FluidAcousticState::try_new(0.7, 425.0).unwrap();
    let thickness = 0.17;
    let stack = TransmittingLayers::new(
        &incident,
        &incident,
        OMEGA,
        &[Layer::FluidLayer {
            fluid: &layer,
            thickness_m: thickness,
        }],
    )
    .unwrap();
    let amplitude = C64::new(0.6, -0.8);
    for cosine in [1.0_f64, 0.9, 0.6, 0.3] {
        let sine = (1.0 - cosine * cosine).sqrt();
        let direction = [0.6 * sine, 0.8 * sine, cosine];
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
        let interior = beam.prepare_layer_field(0, &mut || false).unwrap();
        let kt = OMEGA / 340.0 * sine;
        let kn_squared = (OMEGA / 425.0).powi(2) - kt * kt;
        let yi = cosine / (1.2 * 340.0);
        let rho_omega = 0.7 * OMEGA;
        let [c, s_over_k, k_s] = normal_matrix(kn_squared, thickness);
        let t = C64::from_re(2.0)
            / C64::new(2.0 * c, rho_omega * yi * s_over_k + k_s / (rho_omega * yi));
        for depth in [0.0, 0.2 * thickness, 0.5 * thickness, thickness] {
            let [c, s_over_k, k_s] = normal_matrix(kn_squared, thickness - depth);
            let position = [0.04, -0.03, depth];
            let carrier = amplitude * phase(-kt * (0.6 * position[0] + 0.8 * position[1]));
            let p = carrier * t * C64::new(c, rho_omega * yi * s_over_k);
            let vz = carrier * t * C64::new(yi * c, k_s / rho_omega);
            let sample = interior.field(position, &mut || false).unwrap();
            near(sample.pressure_pa, p, 2e-14);
            near(
                sample.particle_velocity_m_s[0],
                p.scale(0.6 * kt / rho_omega),
                2e-16,
            );
            near(
                sample.particle_velocity_m_s[1],
                p.scale(0.8 * kt / rho_omega),
                2e-16,
            );
            near(sample.particle_velocity_m_s[2], vz, 2e-16);
            near_real(
                sample.active_intensity_w_m2[2],
                0.5 * yi * t.norm_sq(),
                3e-17,
            );
            assert_eq!(sample, interior.field(position, &mut || false).unwrap());
        }
        assert_eq!(interior.layer_index(), 0);
        assert!(core::ptr::eq(interior.spectrum(), &beam));
    }
}

#[test]
fn coherent_interior_layers_keep_interface_work_and_zero_thickness_identity() {
    let incident = fluid();
    let other = FluidAcousticState::try_new(0.4, 500.0).unwrap();
    let layers = [
        Layer::FluidGap { thickness_m: 0.03 },
        Layer::FluidLayer {
            fluid: &other,
            thickness_m: 0.08,
        },
        Layer::FluidGap { thickness_m: 0.0 },
    ];
    let stack = TransmittingLayers::new(&incident, &incident, OMEGA, &layers).unwrap();
    let components = [
        IncidentPlaneWave {
            direction: [0.6, 0.0, 0.8],
            pressure_pa: C64::new(0.7, 0.2),
        },
        IncidentPlaneWave {
            direction: [-0.8, 0.0, 0.6],
            pressure_pa: C64::new(-0.3, 0.5),
        },
    ];
    let beam = stack
        .prepare_angular_spectrum(&components, budget(2), &mut || false)
        .unwrap();
    let interiors: Vec<_> = (0..3)
        .map(|index| beam.prepare_layer_field(index, &mut || false).unwrap())
        .collect();
    let x = 0.013;
    let y = -0.017;
    let mut previous = beam.incident_field([x, y, 0.0], &mut || false).unwrap();
    for (index, thickness) in [0.03, 0.08, 0.0].into_iter().enumerate() {
        let front = interiors[index].field([x, y, 0.0], &mut || false).unwrap();
        near(previous.pressure_pa, front.pressure_pa, 2e-14);
        near(
            previous.particle_velocity_m_s[2],
            front.particle_velocity_m_s[2],
            2e-16,
        );
        near_real(
            previous.active_intensity_w_m2[2],
            front.active_intensity_w_m2[2],
            3e-17,
        );
        previous = interiors[index]
            .field([x, y, thickness], &mut || false)
            .unwrap();
    }
    let transmitted = beam.transmitted_field([x, y, 0.0], &mut || false).unwrap();
    near(previous.pressure_pa, transmitted.pressure_pa, 2e-14);
    near(
        previous.particle_velocity_m_s[2],
        transmitted.particle_velocity_m_s[2],
        2e-16,
    );
    let interior = &interiors[1];
    for point in [[x, y, -0.1], [x, y, 0.09], [f64::NAN, y, 0.01]] {
        assert!(interior.field(point, &mut || false).is_err());
    }
    assert!(beam.prepare_layer_field(3, &mut || false).is_err());
    assert!(matches!(
        beam.prepare_layer_field(0, &mut || true),
        Err(AcousticSpectrumError::Cancelled)
    ));
    let mut polls = 0;
    assert_eq!(
        interior.field([x, y, 0.02], &mut || {
            polls += 1;
            polls == 3
        }),
        Err(AcousticSpectrumError::Cancelled)
    );
}

#[test]
fn opaque_evanescent_layer_preserves_finite_interior_near_entrance() {
    let incident = fluid();
    let layer = FluidAcousticState::try_new(0.6, 680.0).unwrap();
    let thickness = 1.0e100;
    let stack = TransmittingLayers::new(
        &incident,
        &incident,
        OMEGA,
        &[Layer::FluidLayer {
            fluid: &layer,
            thickness_m: thickness,
        }],
    )
    .unwrap();
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
    let interior = beam.prepare_layer_field(0, &mut || false).unwrap();
    let kt = OMEGA * 0.8 / 340.0;
    let alpha = (kt * kt - (OMEGA / 680.0).powi(2)).sqrt();
    let yi = 0.6 / (1.2 * 340.0);
    let yl = C64::new(0.0, -alpha / (0.6 * OMEGA));
    let p0 = C64::from_re(2.0 * yi) / (C64::from_re(yi) + yl);
    for z in [0.0, 0.1 / alpha, 1.0 / alpha, 4.0 / alpha] {
        let p = p0.scale((-alpha * z).exp());
        let sample = interior.field([0.0, 0.0, z], &mut || false).unwrap();
        near(sample.pressure_pa, p, 5e-14);
        near(sample.particle_velocity_m_s[2], p * yl, 2e-16);
        near_real(sample.active_intensity_w_m2[2], 0.0, 3e-17);
        assert!(sample.reactive_intensity_w_m2[2] > 0.0);
    }
    let deep = interior
        .field([0.0, 0.0, 0.5 * thickness], &mut || false)
        .unwrap();
    assert_eq!(deep.pressure_pa, C64::ZERO);
    assert_eq!(deep.particle_velocity_m_s, [C64::ZERO; 3]);
}

#[test]
fn porous_interior_flux_divergence_matches_original_local_constitutive_loss() {
    let incident = fluid();
    let law = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let material = law.effective_fluid(&incident, OMEGA).unwrap();
    let density = material.dynamic_density_kg_m3();
    let compliance = material.dynamic_bulk_modulus_pa().recip();
    for thickness in [0.06, 1000.0] {
        let stack = TransmittingLayers::new(
            &incident,
            &incident,
            OMEGA,
            &[Layer::Porous {
                law: &law,
                thickness_m: thickness,
            }],
        )
        .unwrap();
        let beam = stack
            .prepare_angular_spectrum(
                &[IncidentPlaneWave {
                    direction: [0.36, 0.48, 0.8],
                    pressure_pa: C64::new(0.9, 0.2),
                }],
                budget(1),
                &mut || false,
            )
            .unwrap();
        let interior = beam.prepare_layer_field(0, &mut || false).unwrap();
        let entrance = beam.incident_field([0.0; 3], &mut || false).unwrap();
        let front = interior.field([0.0; 3], &mut || false).unwrap();
        near(front.pressure_pa, entrance.pressure_pa, 3e-14);
        near(
            front.particle_velocity_m_s[2],
            entrance.particle_velocity_m_s[2],
            2e-16,
        );
        let exit = beam.transmitted_field([0.0; 3], &mut || false).unwrap();
        let back = interior
            .field([0.0, 0.0, thickness], &mut || false)
            .unwrap();
        near(back.pressure_pa, exit.pressure_pa, 3e-14);
        near(
            back.particle_velocity_m_s[2],
            exit.particle_velocity_m_s[2],
            2e-16,
        );
        let z = 0.02;
        let h = 1e-6;
        let sample = interior.field([0.0, 0.0, z], &mut || false).unwrap();
        let front_flux = interior
            .field([0.0, 0.0, z - h], &mut || false)
            .unwrap()
            .active_intensity_w_m2[2];
        let back_flux = interior
            .field([0.0, 0.0, z + h], &mut || false)
            .unwrap()
            .active_intensity_w_m2[2];
        let divergence = (back_flux - front_flux) / (2.0 * h);
        let loss = 0.5
            * OMEGA
            * (-density.im
                * sample
                    .particle_velocity_m_s
                    .iter()
                    .map(|v| v.norm_sq())
                    .sum::<f64>()
                - compliance.im * sample.pressure_pa.norm_sq());
        assert!(loss > 0.0);
        near_real(divergence, -loss, 2e-8 * loss);
    }
}

#[test]
fn small_interior_depth_is_not_lost_against_a_long_lossless_layer() {
    let incident = fluid();
    let thickness = 1.0e16;
    let z = 0.013;
    assert_eq!(thickness - z, thickness);
    let stack = TransmittingLayers::new(
        &incident,
        &incident,
        OMEGA,
        &[Layer::FluidGap {
            thickness_m: thickness,
        }],
    )
    .unwrap();
    let beam = stack
        .prepare_angular_spectrum(
            &[IncidentPlaneWave {
                direction: [0.6, 0.0, 0.8],
                pressure_pa: C64::ONE,
            }],
            budget(1),
            &mut || false,
        )
        .unwrap();
    let interior = beam.prepare_layer_field(0, &mut || false).unwrap();
    let sample = interior.field([0.0, 0.0, z], &mut || false).unwrap();
    let expected = phase(-OMEGA / 340.0 * 0.8 * z);
    near(sample.pressure_pa, expected, 2e-14);
    near(
        sample.particle_velocity_m_s[2],
        expected.scale(0.8 / (1.2 * 340.0)),
        2e-16,
    );
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
    let interior = beam.prepare_layer_field(0, &mut || false).unwrap();
    let original_interior = interior.field([0.02, 0.0, 0.03], &mut || false).unwrap();
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
        assert_eq!(
            interior.field([0.02, 0.0, 0.03], &mut || false).unwrap(),
            original_interior
        );
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

// Independent conventional transfer from the receiving face. Growing cos/sin
// are used only for this well-conditioned reference, not for production fields.
fn porous_heat_reference(
    fluid: &ThermoviscousFluidState,
    law: &PorousFluidLaw,
    direction: [f64; 3],
    amplitude: C64,
    thickness: f64,
    depth: f64,
) -> ([C64; 4], [f64; 2]) {
    let response = law.effective_fluid(fluid, OMEGA).unwrap();
    let rho = response.dynamic_density_kg_m3();
    let compliance = response.dynamic_bulk_modulus_pa().recip();
    let kt = [direction[0] * OMEGA / 340.0, direction[1] * OMEGA / 340.0];
    let mut kn = (rho * compliance.scale(OMEGA * OMEGA)
        - C64::from_re(kt[0] * kt[0] + kt[1] * kt[1]))
    .sqrt();
    if kn.im > 0.0 {
        kn = -kn;
    }
    let zn = rho.scale(OMEGA) / kn;
    let y0 = direction[2] / (1.2 * 340.0);
    let from_back = |distance: f64| {
        let z = kn.scale(distance);
        let cos = C64::new(z.re.cos() * z.im.cosh(), -z.re.sin() * z.im.sinh());
        let i_sin = C64::new(-z.re.cos() * z.im.sinh(), z.re.sin() * z.im.cosh());
        [cos + (zn * i_sin).scale(y0), i_sin / zn + cos.scale(y0)]
    };
    let front = from_back(thickness);
    let transmission = amplitude.scale(2.0) / (front[0] + front[1].scale(y0.recip()));
    let at = from_back(thickness - depth);
    let p = transmission * at[0];
    let v = [
        p * C64::from_re(kt[0] / OMEGA) / rho,
        p * C64::from_re(kt[1] / OMEGA) / rho,
        transmission * at[1],
    ];
    let heat = [
        -0.5 * OMEGA * rho.im * v.iter().map(|v| v.norm_sq()).sum::<f64>(),
        -0.5 * OMEGA * compliance.im * p.norm_sq(),
    ];
    ([p, v[0], v[1], v[2]], heat)
}

#[test]
fn g1_coherent_heat_matches_independent_original_boundary_equations() {
    let fluid = fluid();
    let law = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let direction = [0.36, 0.48, 0.8];
    for thickness in [1e-12, 1e-5, 0.06] {
        let stack = TransmittingLayers::new(
            &fluid,
            &fluid,
            OMEGA,
            &[Layer::Porous {
                law: &law,
                thickness_m: thickness,
            }],
        )
        .unwrap();
        // Numerical homogeneity of the linear model, not finite-amplitude physics.
        for scale in [1e-40, 1.0, 1e40] {
            let amplitude = C64::new(0.7, -0.3).scale(scale);
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
            let inside = beam.prepare_layer_field(0, &mut || false).unwrap();
            for fraction in [0.0, 0.2, 0.7, 1.0] {
                let depth = fraction * thickness;
                let sample = inside
                    .field_with_heat([0.0, 0.0, depth], &mut || false)
                    .unwrap();
                let (reference, heat) =
                    porous_heat_reference(&fluid, &law, direction, amplitude, thickness, depth);
                near(sample.field.pressure_pa, reference[0], 8e-13 * scale);
                for (actual, expected) in sample
                    .field
                    .particle_velocity_m_s
                    .into_iter()
                    .zip(&reference[1..])
                {
                    near(actual, *expected, 5e-15 * scale);
                }
                near_real(sample.viscous_heating_w_m3, heat[0], 3e-12 * heat[0]);
                near_real(sample.thermal_heating_w_m3, heat[1], 3e-12 * heat[1]);
            }
        }
    }
}

#[test]
fn g3_coherent_local_heat_integrates_to_original_layer_absorption() {
    use fs_material::porous::layers::LayerAbsorptionBudget;
    let fluid = fluid();
    let law = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let other = FluidAcousticState::try_new(2.1, 410.0).unwrap();
    let thicknesses = [0.04, 0.017, 0.03, 0.025];
    let stack = TransmittingLayers::new(
        &fluid,
        &fluid,
        OMEGA,
        &[
            Layer::Porous {
                law: &law,
                thickness_m: thicknesses[0],
            },
            Layer::FluidGap {
                thickness_m: thicknesses[1],
            },
            Layer::Porous {
                law: &law,
                thickness_m: thicknesses[2],
            },
            Layer::FluidLayer {
                fluid: &other,
                thickness_m: thicknesses[3],
            },
        ],
    )
    .unwrap();
    let wave = IncidentPlaneWave {
        direction: [0.6, 0.0, 0.8],
        pressure_pa: C64::new(0.7, 0.2),
    };
    let beam = stack
        .prepare_angular_spectrum(&[wave], budget(1), &mut || false)
        .unwrap();
    // Unit physical face area: incident watts equal normal W/m² numerically.
    let incident = 0.5 * wave.pressure_pa.norm_sq() * 0.8 / (1.2 * 340.0);
    let work = stack
        .layer_absorption(
            0.8,
            incident,
            LayerAbsorptionBudget { max_layers: 4 },
            || false,
        )
        .unwrap();
    for (i, d) in thicknesses.into_iter().enumerate() {
        let inside = beam.prepare_layer_field(i, &mut || false).unwrap();
        let mut integrated = [0.0; 2];
        // Independent composite Simpson integration, not the analytic loss owner.
        const N: usize = 128;
        for j in 0..=N {
            let sample = inside
                .field_with_heat([0.0, 0.0, d * (j as f64 / N as f64)], &mut || false)
                .unwrap();
            let weight = if j == 0 || j == N {
                1.0
            } else if j % 2 == 0 {
                2.0
            } else {
                4.0
            };
            integrated[0] += weight * sample.viscous_heating_w_m3;
            integrated[1] += weight * sample.thermal_heating_w_m3;
        }
        let expected = work.layers()[i];
        for (value, fraction) in integrated.into_iter().zip([
            expected.viscous_absorbed_fraction,
            expected.thermal_absorbed_fraction,
        ]) {
            near_real(
                value * d / (3.0 * N as f64),
                incident * fraction,
                2e-9 * incident,
            );
        }
    }
}

#[test]
fn g1_coherent_heat_retains_velocity_nodes_and_exact_source_cancellation() {
    let fluid = fluid();
    let law = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let stack = TransmittingLayers::new(
        &fluid,
        &fluid,
        OMEGA,
        &[Layer::Porous {
            law: &law,
            thickness_m: 0.06,
        }],
    )
    .unwrap();
    let wave = IncidentPlaneWave {
        direction: [0.6, 0.0, 0.8],
        pressure_pa: C64::new(0.7, -0.3),
    };
    let single = stack
        .prepare_angular_spectrum(&[wave], budget(1), &mut || false)
        .unwrap();
    let crossing = stack
        .prepare_angular_spectrum(
            &[
                wave,
                IncidentPlaneWave {
                    direction: [-0.6, 0.0, 0.8],
                    ..wave
                },
            ],
            budget(2),
            &mut || false,
        )
        .unwrap();
    let single_layer = single.prepare_layer_field(0, &mut || false).unwrap();
    let crossing_layer = crossing.prepare_layer_field(0, &mut || false).unwrap();
    let depth = 0.021;
    let one = single_layer
        .field_with_heat([0.0, 0.0, depth], &mut || false)
        .unwrap();
    let antinode = crossing_layer
        .field_with_heat([0.0, 0.0, depth], &mut || false)
        .unwrap();
    near(
        antinode.field.pressure_pa,
        one.field.pressure_pa.scale(2.0),
        1e-14,
    );
    assert_eq!(antinode.field.particle_velocity_m_s[0], C64::ZERO);
    near_real(
        antinode.thermal_heating_w_m3,
        4.0 * one.thermal_heating_w_m3,
        1e-14,
    );
    let x = core::f64::consts::FRAC_PI_2 / (0.6 * OMEGA / 340.0);
    let node = crossing_layer
        .field_with_heat([x, 0.0, depth], &mut || false)
        .unwrap();
    near(node.field.pressure_pa, C64::ZERO, 2e-15);
    near(
        node.field.particle_velocity_m_s[0],
        C64::new(0.0, -2.0) * one.field.particle_velocity_m_s[0],
        2e-17,
    );
    assert!(node.thermal_heating_w_m3 < 1e-28 * one.thermal_heating_w_m3);
    let rho = stack
        .prepared_layers()
        .material_response(0)
        .unwrap()
        .dynamic_density_kg_m3();
    let expected_viscous = -2.0 * OMEGA * rho.im * one.field.particle_velocity_m_s[0].norm_sq();
    near_real(
        node.viscous_heating_w_m3,
        expected_viscous,
        2e-13 * expected_viscous,
    );
    assert!(node.viscous_heating_w_m3 > 0.0);
    let cancelled = stack
        .prepare_angular_spectrum(
            &[
                wave,
                IncidentPlaneWave {
                    pressure_pa: -wave.pressure_pa,
                    ..wave
                },
            ],
            budget(2),
            &mut || false,
        )
        .unwrap();
    let cancelled_layer = cancelled.prepare_layer_field(0, &mut || false).unwrap();
    let zero = cancelled_layer
        .field_with_heat([0.03, 0.0, depth], &mut || false)
        .unwrap();
    assert_eq!(zero.field.pressure_pa, C64::ZERO);
    assert_eq!(zero.field.particle_velocity_m_s, [C64::ZERO; 3]);
    assert_eq!(zero.viscous_heating_w_m3, 0.0);
    assert_eq!(zero.thermal_heating_w_m3, 0.0);
}

#[test]
fn g0_coherent_heat_is_a_pure_observation_and_keeps_lossless_identity_semantics() {
    let incident = CountedFluid {
        state: fluid(),
        samples: Cell::new(0),
    };
    let stack = TransmittingLayers::new(
        &incident,
        &incident,
        OMEGA,
        &[
            Layer::FluidGap { thickness_m: 0.1 },
            Layer::FluidGap { thickness_m: 0.0 },
        ],
    )
    .unwrap();
    // The original constructor owns source sampling. This observation must
    // introduce no additional queries, independent of its preparation count.
    let constructor_queries = incident.samples.get();
    let waves = [IncidentPlaneWave {
        direction: [0.6, 0.0, 0.8],
        pressure_pa: C64::ONE,
    }; 2];
    let beam = stack
        .prepare_angular_spectrum(&waves, budget(2), &mut || false)
        .unwrap();
    let inside = beam.prepare_layer_field(0, &mut || false).unwrap();
    assert_eq!(incident.samples.get(), constructor_queries);
    let point = [0.02, 0.03, 0.04];
    let expected = inside.field_with_heat(point, &mut || false).unwrap();
    assert_eq!(expected.field, inside.field(point, &mut || false).unwrap());
    assert_eq!(expected.viscous_heating_w_m3, 0.0);
    assert_eq!(expected.thermal_heating_w_m3, 0.0);
    assert_eq!(incident.samples.get(), constructor_queries);
    let identity = beam.prepare_layer_field(1, &mut || false).unwrap();
    assert!(identity.field([0.0; 3], &mut || false).is_ok());
    assert!(identity.field_with_heat([0.0; 3], &mut || false).is_err());
    assert_eq!(incident.samples.get(), constructor_queries);
    for bad in [[0.0, 0.0, -1e-9], [0.0, 0.0, 0.11], [f64::NAN, 0.0, 0.05]] {
        assert!(inside.field_with_heat(bad, &mut || false).is_err());
        assert_eq!(incident.samples.get(), constructor_queries);
    }
    let mut polls = 0;
    assert!(matches!(
        inside.field_with_heat(point, &mut || {
            polls += 1;
            polls == 3
        }),
        Err(AcousticSpectrumError::Cancelled)
    ));
    assert_eq!(
        inside.field_with_heat(point, &mut || false).unwrap(),
        expected
    );
    assert_eq!(incident.samples.get(), constructor_queries);
}
