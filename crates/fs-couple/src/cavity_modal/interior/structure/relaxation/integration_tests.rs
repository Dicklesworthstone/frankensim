//! Actual FE/Cx binding of the same original full-state response and receivers.
//! Native numerical tests run separately when the complete Cargo graph is absent.

use super::*;
use crate::cavity_modal::interior::structure::relaxation::{
    MaxwellDampingBranch, ResonanceControls, ResonanceSurveyControls,
    StructuralInteriorRelaxationHarmonicLimits,
};
use crate::cavity_modal::interior::structure::resonance::{
    StructuralInteriorReceiver, StructuralInteriorSources,
};

#[test]
fn original_fe_binding_retains_complete_forced_fields_and_actual_receiver_peaks() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 0.2,
        sound_speed_m_s: 1.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let mass = dense(&[1.0], 1);
    let stiffness = dense(&[0.7], 1);
    let trace = geometry
        .boundary_faces()
        .into_iter()
        .map(|face| {
            let (vertices, normal, _) = geometry.boundary_geometry(face).unwrap();
            (
                face,
                vertices.map(|p| (0..3).map(|i| (p[i] - 0.25) * normal[i]).sum()),
            )
        })
        .collect();
    let branches = [MaxwellDampingBranch {
        strain: vec![1.0, 0.0, 0.0, 0.0, 0.0],
        stiffness: 1.25,
        relaxation_time_s: 1.0,
    }];
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: None,
            },
            &[trace],
            &[],
            controls(),
            cx,
        )
        .unwrap();
        let prepared = model
            .prepare_relaxation_response(
                &branches,
                StructuralInteriorRelaxationHarmonicLimits::default(),
                cx,
            )
            .unwrap();
        let f = [C64::new(0.7, -0.2)];
        let q = [C64::new(-0.005, 0.0125); 4];
        let omega = 1.37;
        let s = C64::new(0.0, -omega);
        let u = (f[0] + q[0].scale(2.4) / s) / (s * s + C64::ONE + s.scale(1.25) / (s + C64::ONE));
        let p = u.scale(-0.6) + q[0].scale(4.8) / s;
        let factor = prepared.factor(omega, cx).unwrap();
        let result = factor
            .solve_with_frequency_derivative(&f, &q, None, cx)
            .unwrap();
        assert!((result.response.structural_displacement[0] - u).abs() < 3e-10);
        for &value in &result.response.pressure_pa {
            assert!((value - p).abs() < 3e-10);
        }
        assert!(result.response.power_balance_residual_w.abs() < 3e-10);
        assert!(result.relative_derivative_residual < 1e-9);
        assert!(factor.solve(&[], &q, cx).is_err());
        let replay = factor.solve(&f, &q, cx).unwrap();
        assert_eq!(
            replay.structural_displacement,
            result.response.structural_displacement
        );

        let zero_flow = [C64::ZERO; 4];
        let force = [C64::ONE];
        let sources = StructuralInteriorSources {
            structural_force: &force,
            inward_flow_load: &zero_flow,
        };
        let control = ResonanceControls {
            absolute_frequency_tolerance_rad_s: 1e-9,
            relative_frequency_tolerance: 1e-9,
            maximum_evaluations: 256,
            measure_half_power_bandwidth: true,
        };
        let displacement = prepared
            .driven_resonance(
                (0.8, 2.4),
                sources,
                StructuralInteriorReceiver::StructuralDisplacement(&[1.0]),
                control,
                cx,
            )
            .unwrap();
        let pressure = prepared
            .driven_resonance(
                (0.8, 2.4),
                sources,
                StructuralInteriorReceiver::PressurePoint {
                    tetrahedron: 0,
                    barycentric: [0.25; 4],
                },
                control,
                cx,
            )
            .unwrap();
        assert!(
            (pressure.peak.angular_frequency_rad_s - displacement.peak.angular_frequency_rad_s)
                .abs()
                < 1e-8
        );
        assert!((pressure.observation - displacement.observation.scale(-0.6)).abs() < 1e-8);
        assert!(displacement.half_power.unwrap().response_quality_factor > 0.0);
        let survey = prepared
            .survey_resonances(
                &[0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4],
                sources,
                StructuralInteriorReceiver::StructuralVelocity(&[1.0]),
                ResonanceSurveyControls {
                    absolute_frequency_tolerance_rad_s: 1e-9,
                    relative_frequency_tolerance: 1e-9,
                    maximum_evaluations: 256,
                    maximum_peaks: 2,
                },
                cx,
            )
            .unwrap();
        assert_eq!(survey.peaks.len(), 1);
        assert!(
            (survey.peaks[0].peak.angular_frequency_rad_s
                - displacement.peak.angular_frequency_rad_s)
                .abs()
                > 0.01
        );
    });
}
