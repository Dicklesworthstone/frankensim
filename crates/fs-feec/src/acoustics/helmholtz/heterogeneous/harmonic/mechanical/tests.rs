//! Original FE integration regressions. Full-crate execution requires its
//! constellation dependencies; the complete numerical owner has native tests.

use super::super::BoundaryAdmittanceFrequencySample;
use super::*;
use crate::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_phs::{WallPin, acoustic_outlet::AcousticWallRelaxation};
use fs_rep_mesh::TetComplex;

const POINTS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const MEDIUM: Medium = Medium {
    density_kg_m3: 1.0,
    sound_speed_m_s: 1.0,
};
const WALL_FACE: [u32; 3] = [0, 1, 2];
const ROBIN_FACE: [u32; 3] = [0, 1, 3];
const RELEASE: [u32; 3] = [1, 2, 3];

#[test]
fn lossless_wall_pole_retains_motion_storage_and_is_crossed_by_the_full_survey() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POINTS, 1, MEDIUM).unwrap();
    let model = HeterogeneousHelmholtzFe::new(
        &geometry,
        3,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap();
    let wall = MechanicalBoundary {
        faces: vec![WALL_FACE],
        wall: WallPin {
            surface_density: 1.0,
            stiffness_per_area: 1.0,
            resistance: 0.0,
        },
        relaxation: vec![],
    };
    let robin = ResistiveBoundary {
        faces: vec![ROBIN_FACE],
        admittance_m_pa_s: 0.15,
    };
    let hz = 1.0 / core::f64::consts::TAU;
    let drive = C64::new(0.2, -0.1);
    let traces = [(WALL_FACE, [drive; 3])];
    assert!(BoundaryAdmittanceFrequencySample::mechanical(hz, &wall, &mut || false).is_err());
    let prepared = model
        .prepare_with_coupled_mechanical_walls(
            hz,
            core::slice::from_ref(&robin),
            core::slice::from_ref(&wall),
            &[],
            MechanicalHarmonicLimits::default(),
            &mut || false,
        )
        .unwrap();
    let result = prepared
        .solve_with_frequency_derivative(&traces, None, &mut || false)
        .unwrap();
    assert!(
        result
            .response
            .acoustic
            .pressure_pa
            .iter()
            .all(|p| p.abs() < 1e-12)
    );
    for (i, &global) in prepared.wall_dofs(0).unwrap().iter().enumerate() {
        assert!((result.response.walls[0].normal_velocity_m_s[i] - drive).abs() < 1e-12);
        assert!(
            (result.response.walls[0].displacement_m[i] - C64::new(0.0, -1.0) * drive).abs()
                < 1e-12
        );
        assert!(
            (result.pressure_derivative_pa_per_rad_s[global] - C64::new(0.0, 2.0) * drive).abs()
                < 1e-11
        );
        let kinematic = C64::new(0.0, 1.0)
            * (result.response.walls[0].displacement_m[i]
                + result.wall_derivatives[0].displacement_m_per_rad_s[i]);
        assert!(
            (result.wall_derivatives[0].normal_velocity_m_s_per_rad_s[i] - kinematic).abs() < 1e-11
        );
    }
    // Original wall area=1/2; kinetic and elastic mean stores each equal
    // |v|^2/8. No pressure-only energy account can recover this finite motion.
    assert!((result.response.total_energy_j() - 0.25 * drive.norm_sq()).abs() < 1e-12);
    assert_eq!(result.response.wall_accounting.dissipation_w(), 0.0);
    assert!(result.response.power_balance_residual_w().abs() < 1e-12);
    let grid = [0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0];
    let controls = ResonanceSurveyControls {
        absolute_frequency_tolerance_rad_s: 1e-7,
        relative_frequency_tolerance: 1e-8,
        maximum_evaluations: 256,
        maximum_peaks: 8,
    };
    let observation = (0, [1.0, 0.0, 0.0, 0.0]);
    let survey = model
        .survey_pressure_resonances_with_coupled_mechanical_walls(
            &grid,
            core::slice::from_ref(&robin),
            core::slice::from_ref(&wall),
            &[],
            &traces,
            observation,
            MechanicalHarmonicLimits::default(),
            controls,
            &mut || false,
        )
        .unwrap();
    assert!(!survey.peaks.is_empty());
    for peak in &survey.peaks {
        let at = model
            .prepare_with_coupled_mechanical_walls(
                peak.peak.angular_frequency_rad_s / core::f64::consts::TAU,
                core::slice::from_ref(&robin),
                core::slice::from_ref(&wall),
                &[],
                MechanicalHarmonicLimits::default(),
                &mut || false,
            )
            .unwrap()
            .solve(&traces, &mut || false)
            .unwrap();
        let observed = model
            .pressure_at(&at.acoustic.pressure_pa, observation.0, observation.1)
            .unwrap();
        assert!((observed - peak.observation).abs() < 1e-10 * observed.abs().max(1.0));
        assert!(peak.half_power.is_none());
    }
    assert!(
        model
            .survey_pressure_resonances_with_mechanical_walls(
                &grid,
                &[BoundaryAdmittance {
                    faces: robin.faces.clone(),
                    admittance_m_pa_s: C64::from_re(robin.admittance_m_pa_s),
                }],
                &[wall],
                &[],
                &traces,
                observation,
                super::super::HarmonicLimits::default(),
                controls,
                &mut || false
            )
            .is_err()
    );
}

#[test]
fn original_p2_trace_and_maxwell_heat_match_the_admitted_condensed_equations() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POINTS, 2, MEDIUM).unwrap();
    let model = HeterogeneousHelmholtzFe::new(
        &geometry,
        4,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap();
    let wall = MechanicalBoundary {
        faces: vec![WALL_FACE],
        wall: WallPin {
            surface_density: 0.4,
            stiffness_per_area: 2.0,
            resistance: 0.3,
        },
        relaxation: vec![AcousticWallRelaxation {
            stiffness_per_area: 2.0,
            relaxation_time_s: 0.2,
        }],
    };
    let robin = ResistiveBoundary {
        faces: vec![ROBIN_FACE],
        admittance_m_pa_s: 0.15,
    };
    let hz = 3.0 / core::f64::consts::TAU;
    let traces = [(WALL_FACE, [C64::new(0.2, -0.1); 3])];
    let slopes = [(WALL_FACE, [C64::new(-0.02, 0.03); 3])];
    let prepared = model
        .prepare_with_coupled_mechanical_walls(
            hz,
            core::slice::from_ref(&robin),
            core::slice::from_ref(&wall),
            &[RELEASE],
            MechanicalHarmonicLimits::default(),
            &mut || false,
        )
        .unwrap();
    let full = prepared
        .solve_with_frequency_derivative(&traces, Some(&slopes), &mut || false)
        .unwrap();
    assert_eq!(prepared.wall_dofs(0).unwrap().len(), 6); // Including constrained shared-edge coordinates.
    assert_eq!(full.response.walls[0].relaxation_memory[0].len(), 6);
    let condensed = model
        .prepare_with_dispersive_admittances(
            hz,
            &[
                BoundaryAdmittanceFrequencySample::mechanical(hz, &wall, &mut || false).unwrap(),
                BoundaryAdmittanceFrequencySample {
                    frequency_hz: hz,
                    boundary: BoundaryAdmittance {
                        faces: robin.faces,
                        admittance_m_pa_s: C64::from_re(robin.admittance_m_pa_s),
                    },
                    derivative_m_pa_s_per_rad_s: C64::ZERO,
                },
            ],
            &[RELEASE],
            super::super::HarmonicLimits::default(),
            &mut || false,
        )
        .unwrap()
        .solve_with_frequency_derivative(&traces, Some(&slopes), &mut || false)
        .unwrap();
    for (actual, expected) in full
        .response
        .acoustic
        .pressure_pa
        .iter()
        .zip(&condensed.response.pressure_pa)
        .chain(
            full.pressure_derivative_pa_per_rad_s
                .iter()
                .zip(&condensed.pressure_derivative_pa_per_rad_s),
        )
    {
        assert!((*actual - *expected).abs() < 1e-10 * expected.abs().max(1.0));
    }
    let heat = full.response.acoustic.boundary_dissipation_w
        + full.response.wall_accounting.dissipation_w();
    assert!((heat - condensed.response.boundary_dissipation_w).abs() < 1e-11);
    assert!(full.response.wall_accounting.relaxation_energy_j > 0.0);
    assert!(full.response.wall_accounting.relaxation_dissipation_w > 0.0);
    assert!(
        full.response
            .wall_accounting
            .power_balance_residual_w()
            .abs()
            < 1e-11
    );
    assert!(full.response.fluid_power_balance_residual_w().abs() < 1e-11);
    assert!(full.response.power_balance_residual_w().abs() < 1e-11);
    for i in [1, 2, 3] {
        assert_eq!(full.response.acoustic.pressure_pa[i], C64::ZERO);
        assert_eq!(full.pressure_derivative_pa_per_rad_s[i], C64::ZERO);
    }
}
