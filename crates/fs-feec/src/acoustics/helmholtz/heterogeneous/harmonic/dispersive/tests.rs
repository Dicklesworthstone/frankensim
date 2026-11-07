//! Original FE integration regressions. Native leaf tests independently compare
//! the numerical condensation/tangent with the full wall/cavity state equations.

use super::*;
use crate::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_phs::{WallPin, acoustic_outlet::AcousticWallRelaxation};
use fs_rep_mesh::TetComplex;

const POSITIONS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const MEDIUM: Medium = Medium {
    density_kg_m3: 1.0,
    sound_speed_m_s: 1.0,
};
const WALL: [u32; 3] = [0, 1, 2];
const RELEASE: [u32; 3] = [1, 2, 3];

fn wall(resistance: f64) -> MechanicalBoundary {
    MechanicalBoundary {
        faces: vec![WALL],
        wall: WallPin {
            surface_density: 0.4,
            stiffness_per_area: 2.0,
            resistance,
        },
        relaxation: Vec::new(),
    }
}

#[test]
fn original_face_projection_keeps_signed_slopes_and_frequency_binding() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 2, MEDIUM).unwrap();
    let model = HeterogeneousHelmholtzFe::new(
        &geometry,
        4,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap();
    let mut boundary = wall(0.3);
    boundary.relaxation.push(AcousticWallRelaxation {
        stiffness_per_area: 2.0,
        relaxation_time_s: 0.2,
    });
    let hz = 3.0 / core::f64::consts::TAU;
    let sample =
        BoundaryAdmittanceFrequencySample::mechanical(hz, &boundary, &mut || false).unwrap();
    assert!(sample.derivative_m_pa_s_per_rad_s.re < 0.0);
    let prepared = model
        .prepare_with_dispersive_admittances(
            hz,
            core::slice::from_ref(&sample),
            &[RELEASE],
            HarmonicLimits::default(),
            &mut || false,
        )
        .unwrap();
    let (g, b) = prepared.boundary_frequency_derivative.as_ref().unwrap();
    for (matrix, coefficient) in [
        (g, sample.derivative_m_pa_s_per_rad_s.re),
        (b, sample.derivative_m_pa_s_per_rad_s.im),
    ] {
        let (columns, values) = matrix.row(0);
        let diagonal = values[columns.iter().position(|&col| col == 0).unwrap()];
        // Original hierarchical vertex trace is lambda_0; its squared
        // integral on this area-1/2 triangle is exactly 1/12.
        assert!((diagonal - coefficient / 12.0).abs() < 2e-14 * coefficient.abs());
    }
    let traces = [(WALL, [C64::new(0.2, -0.1); 3])];
    let response = prepared
        .solve_with_frequency_derivative(&traces, None, &mut || false)
        .unwrap();
    let frozen = model
        .prepare_with_pressure_release(
            hz,
            core::slice::from_ref(&sample.boundary),
            &[RELEASE],
            HarmonicLimits::default(),
            &mut || false,
        )
        .unwrap()
        .solve_with_frequency_derivative(&traces, None, &mut || false)
        .unwrap();
    assert_eq!(response.response.pressure_pa, frozen.response.pressure_pa);
    assert!(
        (response.pressure_derivative_pa_per_rad_s[0] - frozen.pressure_derivative_pa_per_rad_s[0])
            .abs()
            > 1e-5
    );
    for vertex in [1, 2, 3] {
        assert_eq!(response.response.pressure_pa[vertex], C64::ZERO);
        assert_eq!(response.pressure_derivative_pa_per_rad_s[vertex], C64::ZERO);
    }
    assert!(
        model
            .prepare_with_dispersive_admittances(
                hz * 1.01,
                core::slice::from_ref(&sample),
                &[RELEASE],
                HarmonicLimits::default(),
                &mut || false,
            )
            .is_err()
    );
    assert!(
        model
            .prepare_with_dispersive_admittances(
                hz,
                &[sample.clone(), sample],
                &[RELEASE],
                HarmonicLimits::default(),
                &mut || false,
            )
            .is_err()
    );
}

#[test]
fn mechanical_wall_survey_finds_the_two_coupled_discrete_pressure_peaks() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    // One free P1 pressure coordinate gives an independently solvable
    // discrete two-oscillator regression, not a continuum accuracy claim.
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = HeterogeneousHelmholtzFe::new(
        &geometry,
        3,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap();
    let boundary = wall(0.02);
    let drive = C64::new(0.2, -0.1);
    let grid = [1.0, 1.6, 1.9, 2.1, 2.5, 5.0, 6.0, 6.7, 7.5, 8.0];
    let survey = model
        .survey_pressure_resonances_with_mechanical_walls(
            &grid,
            &[],
            &[boundary],
            &[RELEASE],
            &[(WALL, [drive; 3])],
            (0, [1.0, 0.0, 0.0, 0.0]),
            HarmonicLimits::default(),
            ResonanceSurveyControls {
                absolute_frequency_tolerance_rad_s: 1e-8,
                relative_frequency_tolerance: 1e-9,
                maximum_evaluations: 256,
                maximum_peaks: 4,
            },
            &mut || false,
        )
        .unwrap();
    assert_eq!(survey.peaks.len(), 2);
    let analytic = |omega: f64| {
        // Exact original P1 forms: K00=1/2, M00=1/60, W00=1/12,
        // integral(lambda_0 v_in)=v_in/6. Use wall dynamic stiffness,
        // independently of the production impedance condensation.
        let d = C64::new(2.0 - 0.4 * omega * omega, 0.02 * omega);
        let a = C64::from_re(0.5 - omega * omega / 60.0) - C64::from_re(omega * omega / 12.0) / d;
        C64::new(0.0, omega / 6.0) * drive / a
    };
    for peak in survey.peaks {
        let omega = peak.peak.angular_frequency_rad_s;
        let expected = analytic(omega);
        assert!((peak.observation - expected).abs() < 2e-9 * expected.abs());
        assert!(expected.abs() > analytic(omega - 0.01).abs());
        assert!(expected.abs() > analytic(omega + 0.01).abs());
        assert!(peak.half_power.is_none());
    }
}
