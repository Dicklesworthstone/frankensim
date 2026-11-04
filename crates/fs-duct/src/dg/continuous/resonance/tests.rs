use super::*;
use crate::dg::{DgBoundary, DgMesh};
use crate::{Duct, Segment};
use fs_material::gas::{GasSpec, GasState};

const PORT: PortId = PortId(7);
const WORK: usize = 10_000_000;
const LENGTH: f64 = 0.6;
const RADIUS: f64 = 0.008;

fn duct() -> (DgDuct, f64) {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * RADIUS * RADIUS);
    let fundamental = core::f64::consts::PI * gas.sound_speed / LENGTH;
    let duct = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: LENGTH,
            }],
        },
        &gas,
        DgMesh {
            elements: 4,
            order: 4,
        },
        DgBoundary::Port(PORT),
        DgBoundary::resistive(20.0 * zc).unwrap(),
    )
    .unwrap();
    (duct, fundamental)
}

fn controls() -> ResonanceSurveyControls {
    ResonanceSurveyControls {
        absolute_frequency_tolerance_rad_s: 1e-5,
        relative_frequency_tolerance: 1e-9,
        maximum_evaluations: 256,
        maximum_peaks: 8,
    }
}

fn grid(fundamental: f64) -> Vec<f64> {
    // Three half-wave peaks, each with resolved rising and falling flanks.
    // No point is deliberately placed at an expected peak or minimum.
    (0..=24)
        .map(|i| fundamental * (0.67 + f64::from(i) * 0.11))
        .collect()
}

fn pressure(duct: &DgDuct, omega: f64, inflow: C64) -> C64 {
    duct.harmonic_response(omega, &[(PORT, inflow)], WORK, 1e-10, &mut || false)
        .unwrap()
        .ports[0]
        .pressure_pa
}

fn near(actual: C64, expected: C64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * actual.abs().max(expected.abs()) + 1e-20,
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn survey_resolves_three_physical_duct_pressure_peaks_and_preserves_complex_gain() {
    let (duct, fundamental) = duct();
    let grid = grid(fundamental);
    let flow = C64::new(1e-7, 4e-8);
    let identity = duct.operator_identity();
    let survey = duct
        .survey_port_pressure_resonances(
            &grid,
            &[(PORT, flow)],
            PORT,
            controls(),
            WORK,
            1e-10,
            &mut || false,
        )
        .unwrap();
    assert_eq!(survey.peaks.len(), 3);
    assert_eq!(survey.grid_points, grid.len());
    assert!(survey.evaluations <= controls().maximum_evaluations);
    assert!(survey.maximum_relative_residual <= 1e-10);
    assert_eq!(duct.operator_identity(), identity);
    for (i, peak) in survey.peaks.iter().enumerate() {
        let omega = peak.peak.angular_frequency_rad_s;
        // Independently, a uniform lossless line terminated by R > Zc has
        // pressure-response maxima at kL = n*pi, regardless of finite R.
        // The tolerance admits the original spatial discretization error.
        let expected = fundamental * (i + 1) as f64;
        assert!((omega / expected - 1.0).abs() < 0.003);
        assert!(peak.half_power.is_none());
        near(peak.observation, pressure(&duct, omega, flow), 1e-12);
        for neighbor in [omega - 0.01 * fundamental, omega + 0.01 * fundamental] {
            assert!(peak.observation.abs() > pressure(&duct, neighbor, flow).abs());
        }
    }

    let gain = C64::new(-2.0, 0.75);
    let scaled = duct
        .survey_port_pressure_resonances(
            &grid,
            &[(PORT, flow * gain)],
            PORT,
            controls(),
            WORK,
            1e-10,
            &mut || false,
        )
        .unwrap();
    assert_eq!(scaled.peaks.len(), survey.peaks.len());
    for (original, scaled) in survey.peaks.iter().zip(&scaled.peaks) {
        assert!(
            (original.peak.angular_frequency_rad_s - scaled.peak.angular_frequency_rad_s).abs()
                <= 2e-5
        );
        near(scaled.observation, original.observation * gain, 2e-6);
    }
}

#[test]
fn survey_preserves_native_port_work_band_and_cancellation_refusals() {
    let (duct, fundamental) = duct();
    let grid = grid(fundamental);
    let flows = [(PORT, C64::from_re(1e-7))];
    let unknown = [(PortId(99), flows[0].1)];
    for (inflows, observed, work) in [
        (&flows[..], PortId(99), WORK),
        (&unknown[..], PORT, WORK),
        (&flows[..], PORT, 1),
    ] {
        assert!(matches!(
            duct.survey_port_pressure_resonances(
                &grid,
                inflows,
                observed,
                controls(),
                work,
                1e-10,
                &mut || false,
            ),
            Err(ResonanceSurveyError::Search(ResonanceSearchError::Source(
                DuctError::BadParameter { .. }
            )))
        ));
    }
    assert!(matches!(
        duct.survey_port_pressure_resonances(
            &grid,
            &flows,
            PORT,
            controls(),
            WORK,
            1e-10,
            &mut || true,
        ),
        Err(ResonanceSurveyError::Search(ResonanceSearchError::Source(
            DuctError::Cancelled { .. }
        )))
    ));
    // Wall diffusion depends on omega * rho * radius² / viscosity. This
    // low-frequency band is resolved with the declared eight memory poles.
    let bounded = duct
        .with_viscothermal(8, [50.0 / 64.0, 2000.0 / 64.0], 0.02)
        .unwrap();
    let band = bounded.viscothermal_band_rad_s().unwrap();
    let outside = [band[1] * 1.01, band[1] * 1.02, band[1] * 1.03];
    let native = bounded
        .harmonic_response_with_frequency_derivative(
            outside[0],
            &flows,
            None,
            WORK,
            1e-10,
            &mut || false,
        )
        .unwrap_err();
    assert_eq!(
        bounded.survey_port_pressure_resonances(
            &outside,
            &flows,
            PORT,
            controls(),
            WORK,
            1e-10,
            &mut || false,
        ),
        Err(ResonanceSurveyError::Search(ResonanceSearchError::Source(native)))
    );
}
