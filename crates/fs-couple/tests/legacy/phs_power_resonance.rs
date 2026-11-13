//! End-to-end original physical heat/flux peaks, widths and parameter motion.
use fs_couple::phs_power_resonance::*;
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{PowerExtremumKind, ResonanceError};
use fs_phs::harmonic::scattering::{ScatteringBudget, ScatteringModel};
use fs_phs::harmonic::{HarmonicBudget, HarmonicError, HarmonicModel, HarmonicParameterDirection};
use fs_phs::mechanical::{
    LinearMechanicalSystem, MechanicalHarmonicDirection, MechanicalHarmonicDirectionControls,
    MechanicalRelaxationDirection,
};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch};

fn native_budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 32,
        max_dense_bytes: 4_000_000,
        max_dense_work: 8_000_000,
        relative_residual_tolerance: 2e-11,
    }
}
fn controls() -> ResonanceControls {
    ResonanceControls {
        absolute_angular_tolerance_rad_s: 2e-12,
        relative_angular_tolerance: 2e-12,
        maximum_evaluations: 400,
    }
}
fn oscillator(m: f64, c: f64, k: f64, g: f64) -> HarmonicModel {
    let system = PortHamiltonian::new(
        2,
        1,
        vec![0., 1., -1., 0.],
        vec![0., 0., 0., c],
        vec![0., g],
        Box::new(QuadraticStorage::new(vec![k, 0., 0., 1. / m], 2).unwrap()),
    )
    .unwrap();
    HarmonicModel::from_ode(&system, native_budget(), || false).unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-9 * b.abs().max(1e-4), "{a:e} != {b:e}");
}

/// Frequency-dependent coherent phasors are lent from reusable source arrays.
/// Phase and phase-parameter laws must cancel from this scalar physical heat.
struct Source {
    arrays: [[C64; 1]; 5],
    calls: [usize; 3],
}
impl Source {
    fn new() -> Self {
        Self {
            arrays: [[C64::ZERO; 1]; 5],
            calls: [0; 3],
        }
    }
}
impl PhsPowerSource for Source {
    fn jet(
        &mut self,
        w: f64,
        order: PhsPowerSourceOrder,
    ) -> Result<PhsPowerSourceJet<'_>, HarmonicError> {
        let (delay, delay_theta, gain_theta) = (0.31, 0.07, 0.12);
        let a = C64::new(0.9, -0.2) * C64::new((delay * w).cos(), (delay * w).sin());
        self.arrays[0][0] = a;
        self.arrays[1][0] = a * C64::new(0., delay);
        self.arrays[2][0] = a.scale(-delay * delay);
        self.arrays[3][0] = a * C64::new(gain_theta, w * delay_theta);
        self.arrays[4][0] = self.arrays[3][0] * C64::new(0., delay) + a * C64::new(0., delay_theta);
        self.calls[match order {
            PhsPowerSourceOrder::Frequency => 0,
            PhsPowerSourceOrder::Curvature => 1,
            PhsPowerSourceOrder::Parameter => 2,
        }] += 1;
        Ok(PhsPowerSourceJet {
            input: &self.arrays[0],
            angular_derivative: Some(&self.arrays[1]),
            angular_second_derivative: Some(&self.arrays[2]),
            parameter_derivative: Some(&self.arrays[3]),
            angular_parameter_derivative: Some(&self.arrays[4]),
        })
    }
}

#[test]
fn original_mechanical_heat_width_and_material_sensitivity_keep_source_phase() {
    let (m, c, k, g) = (1.3, 0.17, 4.2, 0.9);
    let (dm, dc, dk, dg) = (0.2, -0.03, 0.31, -0.07);
    let model = oscillator(m, c, k, g);
    let query = PhsPowerQuery::new(
        PhsPowerModel::Original(&model),
        PhsPowerObservation::Dissipation,
    )
    .unwrap();
    let dq = [dk, 0., 0., -dm / (m * m)];
    let dr = [0., 0., 0., dc];
    let dp = [0., dg];
    let direction = PhsPowerParameterDirection {
        physical: HarmonicParameterDirection {
            storage: Some(&dq),
            dissipation: Some(&dr),
            ports: Some(&dp),
            ..Default::default()
        },
        references: None,
    };
    let mut source = Source::new();
    let peak = query
        .extremum_with_parameter_derivative((1.2, 2.3), &mut source, direction, controls(), || {
            false
        })
        .unwrap();
    let omega = (k / m).sqrt();
    let omega_theta = 0.5 * omega * (dk / k - dm / m);
    let height = 0.5 * g * g * (0.9 * 0.9 + 0.2 * 0.2) / c;
    close(peak.extremum.angular_frequency_rad_s, omega);
    close(peak.angular_frequency_derivative_rad_s, omega_theta);
    close(peak.extremum.sample.average_power_w, height);
    close(
        peak.average_power_derivative,
        height * (2. * dg / g + 2. * 0.12 - dc / c),
    );
    assert_eq!(source.calls[1..], [0, 1]);
    let band = query
        .bandwidth_with_parameter_derivative(
            &peak,
            (0.8, peak.extremum.angular_frequency_rad_s),
            (peak.extremum.angular_frequency_rad_s, 3.),
            &mut source,
            direction,
            controls(),
            || false,
        )
        .unwrap();
    let width = c / m;
    let width_theta = dc / m - c * dm / (m * m);
    close(band.bandwidth.angular_bandwidth_rad_s, width);
    close(band.angular_bandwidth_derivative_rad_s, width_theta);
    close(band.bandwidth.response_quality_factor, omega / width);
    close(
        band.response_quality_factor_derivative,
        omega_theta / width - omega * width_theta / (width * width),
    );
    assert_eq!(source.calls[1..], [0, 3]);
    assert!(source.calls[0] > 3);
}

#[test]
fn loaded_material_heat_and_outgoing_flux_keep_reference_feedback_and_its_direction() {
    let (m, c, k, g, z) = (1.3, 0.17, 4.2, 0.9, 2.4);
    let (dm, dc, dk, dg, dz) = (0.2, -0.03, 0.31, -0.07, 0.23);
    let model = ScatteringModel::new(
        oscillator(m, c, k, g),
        &[(0, z)],
        ScatteringBudget {
            maximum_ports: 2,
            maximum_payload_bytes: 1_000_000,
            maximum_work: 2_000_000,
        },
        || false,
    )
    .unwrap();
    let heat = PhsPowerQuery::new(
        PhsPowerModel::Scattering(&model),
        PhsPowerObservation::Dissipation,
    )
    .unwrap();
    let outgoing = PhsPowerQuery::new(
        PhsPowerModel::Scattering(&model),
        PhsPowerObservation::OutgoingWave(0),
    )
    .unwrap();
    let dq = [dk, 0., 0., -dm / (m * m)];
    let dr = [0., 0., 0., dc];
    let dp = [0., dg];
    let dzs = [dz];
    let direction = PhsPowerParameterDirection {
        physical: HarmonicParameterDirection {
            storage: Some(&dq),
            dissipation: Some(&dr),
            ports: Some(&dp),
            ..Default::default()
        },
        references: Some(&dzs),
    };
    let mut source = Source::new();
    let peak = heat
        .extremum_with_parameter_derivative((1., 2.6), &mut source, direction, controls(), || false)
        .unwrap();
    let omega = (k / m).sqrt();
    let omega_theta = 0.5 * omega * (dk / k - dm / m);
    let total_c = c + g * g / z;
    let total_c_theta = dc + 2. * g * dg / z - g * g * dz / (z * z);
    let a2 = 0.9 * 0.9 + 0.2 * 0.2;
    let height = 2. * c * g * g * a2 / (z * z * total_c * total_c);
    let height_theta =
        height * (dc / c + 2. * dg / g + 2. * 0.12 - 2. * dz / z - 2. * total_c_theta / total_c);
    close(peak.extremum.angular_frequency_rad_s, omega);
    close(peak.angular_frequency_derivative_rad_s, omega_theta);
    close(peak.extremum.sample.average_power_w, height);
    close(peak.average_power_derivative, height_theta);
    let band = heat
        .bandwidth_with_parameter_derivative(
            &peak,
            (0.7, peak.extremum.angular_frequency_rad_s),
            (peak.extremum.angular_frequency_rad_s, 3.),
            &mut source,
            direction,
            controls(),
            || false,
        )
        .unwrap();
    close(band.bandwidth.angular_bandwidth_rad_s, total_c / m);
    close(
        band.angular_bandwidth_derivative_rad_s,
        total_c_theta / m - total_c * dm / (m * m),
    );
    let notch = outgoing
        .extremum_with_parameter_derivative((1., 2.6), &mut source, direction, controls(), || false)
        .unwrap();
    assert_eq!(notch.extremum.kind, PowerExtremumKind::Minimum);
    close(notch.extremum.angular_frequency_rad_s, omega);
    close(notch.angular_frequency_derivative_rad_s, omega_theta);
    close(
        notch.extremum.sample.average_power_w,
        a2 / (2. * z) - height,
    );
    close(
        notch.average_power_derivative,
        a2 / (2. * z) * (2. * 0.12 - dz / z) - height_theta,
    );
    // A complete value-only survey requests one three-solve curvature, not a
    // dummy material direction. Grid calls and final calls share one cap.
    let before = source.calls;
    let survey = outgoing
        .survey(
            &[0.4, 1.2, 1.6, 2.1, 3.],
            &mut source,
            controls(),
            2,
            || false,
        )
        .unwrap();
    assert_eq!(survey.extrema.len(), 1);
    assert_eq!(survey.extrema[0].kind, PowerExtremumKind::Minimum);
    assert_eq!(source.calls[1] - before[1], 1);
    assert_eq!(source.calls[2] - before[2], 0);
    let before = source.calls;
    assert!(matches!(
        heat.extremum((1., 2.6), &mut source, controls(), || true),
        Err(ResonanceError::Cancelled)
    ));
    assert_eq!(source.calls, before);
    assert!(
        PhsPowerQuery::new(
            PhsPowerModel::Scattering(&model),
            PhsPowerObservation::OutgoingWave(1)
        )
        .is_err()
    );
}

#[test]
fn causal_maxwell_heat_peak_and_asymmetric_width_move_with_original_relaxation_time() {
    let mechanical = LinearMechanicalSystem::new(
        &[1.],
        &[0.],
        &[1.],
        &[1.],
        1,
        1,
        vec![RelaxationBranch {
            projection: vec![1.],
            stiffness: 11. / 8.,
            relaxation_time_s: 1.,
        }],
    )
    .unwrap();
    let model = HarmonicModel::from_ode(mechanical.system(), native_budget(), || false).unwrap();
    let arm = [MechanicalRelaxationDirection {
        relaxation_time_s: 1.,
        ..Default::default()
    }];
    let prepared = mechanical
        .prepare_harmonic_parameter_direction(
            MechanicalHarmonicDirection {
                branches: Some(&arm),
                ..Default::default()
            },
            MechanicalHarmonicDirectionControls {
                max_state_dimension: 16,
                max_dense_bytes: 1_000_000,
                max_multiply_terms: 1_000_000,
            },
            || false,
        )
        .unwrap();
    let direction = PhsPowerParameterDirection {
        physical: prepared.as_harmonic_direction(),
        references: None,
    };
    let query = PhsPowerQuery::new(
        PhsPowerModel::Original(&model),
        PhsPowerObservation::Dissipation,
    )
    .unwrap();
    let input = [C64::ONE];
    let mut source = input.as_slice();
    let peak = query
        .extremum_with_parameter_derivative((1.2, 1.6), &mut source, direction, controls(), || {
            false
        })
        .unwrap();
    let omega = 2.0_f64.sqrt();
    close(peak.extremum.angular_frequency_rad_s, omega);
    close(peak.angular_frequency_derivative_rad_s, 1. / (3. * omega));
    close(peak.extremum.sample.average_power_w, 44. / 41.);
    close(peak.average_power_derivative, 1012. / 1681.);
    let band = query
        .bandwidth_with_parameter_derivative(
            &peak,
            (0.5, peak.extremum.angular_frequency_rad_s),
            (peak.extremum.angular_frequency_rad_s, 2.5),
            &mut source,
            direction,
            controls(),
            || false,
        )
        .unwrap();
    // Independent cleared physical denominator, not a fitted Lorentzian:
    // P=(11/16)*x / [(1-x)^2+x*(19/8-x)^2], x=omega².
    // These positive half-height roots were independently evaluated at90digits.
    let lower = 1.126_258_885_627_684_5;
    let upper = 1.663_281_746_433_824;
    close(band.bandwidth.lower_angular_frequency_rad_s, lower);
    close(band.bandwidth.upper_angular_frequency_rad_s, upper);
    assert!(((omega - lower) - (upper - omega)).abs() > 0.03);
    let edge_direction = |w: f64| {
        let x = w * w;
        let d = (1. - x) * (1. - x) + x * (19. / 8. - x) * (19. / 8. - x);
        let dx = 233. / 64. - (15. / 2.) * x + 3. * x * x;
        let dt = 2. * x * (19. / 8. - x) * (19. / 8. - x);
        let p = (11. / 16.) * x / d;
        let pt = p * (1. - dt / d);
        let pw = 2. * w * (11. / 16.) * (d - x * dx) / (d * d);
        (0.5 * (1012. / 1681.) - pt) / pw
    };
    close(
        band.lower_angular_frequency_derivative_rad_s,
        edge_direction(lower),
    );
    close(
        band.upper_angular_frequency_derivative_rad_s,
        edge_direction(upper),
    );
    let width = upper - lower;
    let width_theta = edge_direction(upper) - edge_direction(lower);
    close(band.angular_bandwidth_derivative_rad_s, width_theta);
    close(
        band.response_quality_factor_derivative,
        1. / (3. * omega * width) - omega * width_theta / (width * width),
    );
    let survey = query
        .survey_with_parameter_derivative(
            &[0.4, 1.2, 1.6, 2.2],
            &mut source,
            direction,
            controls(),
            2,
            || false,
        )
        .unwrap();
    assert_eq!(survey.extrema.len(), 1);
    close(
        survey.extrema[0].angular_frequency_derivative_rad_s,
        1. / (3. * omega),
    );
}
