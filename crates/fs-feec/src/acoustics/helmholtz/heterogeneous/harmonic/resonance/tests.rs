//! G1/G3 original-field survey wiring. Full FE execution requires the complete
//! constellation dependency environment; the numerical tangent has native tests.

use super::*;
use crate::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use crate::fixtures::masked_cube_grid;
use fs_math::c64::C64;

const MEDIUM: Medium = Medium {
    density_kg_m3: 1.3,
    sound_speed_m_s: 2.0,
};

fn controls() -> ResonanceSurveyControls {
    ResonanceSurveyControls {
        absolute_frequency_tolerance_rad_s: 1e-7,
        relative_frequency_tolerance: 1e-8,
        maximum_evaluations: 128,
        maximum_peaks: 4,
    }
}

fn end(fe: &HelmholtzFe<'_>, x: f64) -> Vec<[u32; 3]> {
    fe.boundary_faces()
        .into_iter()
        .filter(|face| {
            face.iter()
                .all(|&v| (fe.positions()[v as usize][0] - x).abs() < 1e-12)
        })
        .collect()
}

fn observation(fe: &HelmholtzFe<'_>, x: f64) -> (usize, [f64; 4]) {
    for (element, vertices) in fe.tetrahedra().iter().enumerate() {
        for (local, &vertex) in vertices.iter().enumerate() {
            if (fe.positions()[vertex as usize][0] - x).abs() < 1e-12 {
                let mut barycentric = [0.0; 4];
                barycentric[local] = 1.0;
                return (element, barycentric);
            }
        }
    }
    panic!("fixture end vertex is missing")
}

#[test]
fn loaded_duct_survey_locates_the_original_complex_pressure_peak() {
    let (complex, mut points) = masked_cube_grid(4, 1, 1, &|_, _, _| true);
    for p in &mut points {
        p[0] /= 4.0;
        p[1] *= 0.2;
        p[2] *= 0.3;
    }
    let geometry = HelmholtzFe::new(&complex, &points, 2, MEDIUM).unwrap();
    let model = HeterogeneousHelmholtzFe::new(
        &geometry,
        4,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap();
    let drive = C64::new(0.02, -0.01);
    let traces: Vec<_> = end(&geometry, 0.0)
        .into_iter()
        .map(|face| (face, [drive; 3]))
        .collect();
    let boundaries = [BoundaryAdmittance {
        faces: end(&geometry, 1.0),
        admittance_m_pa_s: C64::from_re(0.2 / (MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s)),
    }];
    let sample = observation(&geometry, 0.0);
    let survey = model
        .survey_pressure_resonances(
            &[4.8, 5.6, 6.4, 7.2, 8.0],
            &boundaries,
            &[],
            &traces,
            sample,
            HarmonicLimits::default(),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(survey.peaks.len(), 1);
    let peak = &survey.peaks[0];
    let omega = peak.peak.angular_frequency_rad_s;
    // Independently, a lossless line with finite R > rho*c peaks at kL=pi.
    // The tolerance admits the P2 spatial dispersion of this original mesh.
    assert!((omega / (core::f64::consts::PI * MEDIUM.sound_speed_m_s) - 1.0).abs() < 0.005);
    assert!(peak.half_power.is_none());
    let at = |w| {
        let response = model
            .prepare(
                w / core::f64::consts::TAU,
                &boundaries,
                HarmonicLimits::default(),
                &mut || false,
            )
            .unwrap()
            .solve(&traces, &mut || false)
            .unwrap();
        model
            .pressure_at(&response.pressure_pa, sample.0, sample.1)
            .unwrap()
    };
    let actual = at(omega);
    assert!((actual - peak.observation).abs() < 1e-11 * actual.abs());
    assert!(actual.abs() > at(omega - 0.05).abs());
    assert!(actual.abs() > at(omega + 0.05).abs());
    assert!(survey.maximum_relative_residual <= HarmonicLimits::default().residual_tolerance);
}

#[test]
fn exact_pressure_release_observation_is_silent_and_refusals_remain_native() {
    let (complex, points) = masked_cube_grid(1, 1, 1, &|_, _, _| true);
    let geometry = HelmholtzFe::new(&complex, &points, 2, MEDIUM).unwrap();
    let model = HeterogeneousHelmholtzFe::new(
        &geometry,
        4,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap();
    let inlet = end(&geometry, 0.0);
    let release = end(&geometry, 1.0);
    let traces: Vec<_> = inlet.iter().map(|&face| (face, [C64::ONE; 3])).collect();
    let boundaries = [BoundaryAdmittance {
        faces: inlet,
        admittance_m_pa_s: C64::from_re(0.1),
    }];
    let sample = observation(&geometry, 1.0);
    let grid = [1.8, 2.5, 3.1, 3.8, 4.5];
    let survey = model
        .survey_pressure_resonances(
            &grid,
            &boundaries,
            &release,
            &traces,
            sample,
            HarmonicLimits::default(),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(survey.peaks.is_empty());
    assert_eq!(survey.evaluations, grid.len());
    let invalid = model.survey_pressure_resonances(
        &grid,
        &boundaries,
        &release,
        &traces,
        (usize::MAX, [0.25; 4]),
        HarmonicLimits::default(),
        controls(),
        &mut || false,
    );
    assert!(matches!(
        invalid,
        Err(ResonanceSurveyError::Search(ResonanceSearchError::Source(
            _
        )))
    ));
    let cancelled = model.survey_pressure_resonances(
        &grid,
        &boundaries,
        &release,
        &traces,
        sample,
        HarmonicLimits::default(),
        controls(),
        &mut || true,
    );
    assert!(matches!(
        cancelled,
        Err(ResonanceSurveyError::Search(ResonanceSearchError::Source(
            HeterogeneousError::Cancelled
        )))
    ));
}
