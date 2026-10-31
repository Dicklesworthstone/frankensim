//! G1/G3/G4: differential inductance of the actual nonlinear spatial field.
//! Manufactured H(B) sources test implementation, not measured core materials.

#[path = "nonlinear_winding_circuit/common.rs"]
mod common;

use fs_em::{
    planar::{
        MagneticError,
        nonlinear::{NonlinearMagneticError, PlanarNonlinearMagnetostatics},
    },
    winding::{
        WindingConfig, WindingError,
        incremental::{
            IncrementalWindingConfig, IncrementalWindingError, PlanarIncrementalWindingBasis,
        },
    },
};

fn config() -> IncrementalWindingConfig {
    IncrementalWindingConfig {
        winding: WindingConfig {
            max_windings: 4,
            max_basis_values: 1000,
            pairing_tolerance_h: 1e-10,
            current_balance_tolerance_turns: 1e-10,
            inverse_residual_tolerance: 1e-10,
        },
        minimum_smooth_field_margin_t: 1e-8,
    }
}

fn basis<'a>(
    problem: &'a PlanarNonlinearMagnetostatics,
    current: f64,
    end: f64,
) -> PlanarIncrementalWindingBasis<'a> {
    PlanarIncrementalWindingBasis::new(
        problem,
        common::profiles(1),
        &common::boundary(),
        vec![current],
        vec![end],
        &[0.0; 8],
        config(),
        common::material::controls().nonlinear,
        &mut || false,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-9 * (1.0 + expected.abs()),
        "{actual} != {expected}"
    );
}

#[test]
fn saturated_inductance_is_the_flux_derivative_not_secant_or_energy_equivalent() {
    let depth = 0.7;
    let source = common::source(true, 2.0);
    let problem = common::problem(&source, depth);
    for current in [-2.0, -0.5, 0.0, 0.5, 2.0] {
        let model = basis(&problem, current, 0.03);
        // Each square's single free potential gives H(2A)=i. The signed
        // outgoing/return pairing is lambda_active=2*depth*B(i).
        let slope = if current.abs() > 1.0 { 4.0 } else { 1.0 };
        let active = 2.0 * depth / slope;
        near(model.active_inductance_h()[0], active);
        near(model.inductance_h()[0], active + 0.03);
        near(
            model.operating_flux_linkages_wb()[0],
            common::linkage(current, depth, true) + 0.03 * current,
        );
        near(
            model.operating_field().energy_j(),
            2.0 * depth * common::density_energy(common::b_from_current(current, true), true),
        );
        near(
            model.voltage_from_current_rate_v(&[3.0]).unwrap()[0],
            3.0 * (active + 0.03),
        );
        let delta = model.linearized_flux_increment_wb(&[0.0002]).unwrap();
        near(
            model.linearized_current_increment_a(&delta).unwrap()[0],
            0.0002,
        );
        assert!(model.maximum_pairing_defect_h() < 1e-10);
        assert!(model.inverse_residual() < 1e-10);
        assert!(model.unit_responses()[0].report().converged_euclidean());
        let report = model.operating_field().report();
        assert!(report.residual_norm_a_m <= report.residual_target_a_m);
        if current.abs() > 1.0 {
            let secant = (model.operating_flux_linkages_wb()[0] - 0.03 * current) / current;
            let energy_equivalent = 2.0 * model.operating_field().energy_j() / current.powi(2);
            assert!(secant > 2.0 * active);
            assert!((energy_equivalent - active).abs() > 0.1);
        }
    }
}

#[test]
fn signed_mutual_tangent_depth_and_end_storage_follow_the_same_nonlinear_field() {
    let source = common::source(true, 2.0);
    for depth in [0.3, 0.9] {
        let problem = common::problem(&source, depth);
        let primary = common::profiles(1).remove(0);
        let profiles = vec![primary.clone(), primary.iter().map(|j| -3.0 * j).collect()];
        let currents = [2.6, 0.2]; // i_effective = 2 A, in the saturated segment.
        let ends = [0.02, 0.04];
        let model = PlanarIncrementalWindingBasis::new(
            &problem,
            profiles.clone(),
            &common::boundary(),
            currents.to_vec(),
            ends.to_vec(),
            &[0.0; 8],
            config(),
            common::material::controls().nonlinear,
            &mut || false,
        )
        .unwrap();
        let active = 2.0 * depth / 4.0;
        for (&a, b) in model.active_inductance_h().iter().zip([
            active,
            -3.0 * active,
            -3.0 * active,
            9.0 * active,
        ]) {
            near(a, b);
        }
        assert_eq!(
            model.inductance_h()[1].to_bits(),
            model.inductance_h()[2].to_bits()
        );
        near(model.inductance_h()[0], active + ends[0]);
        near(model.inductance_h()[3], 9.0 * active + ends[1]);
        let flux = common::linkage(2.0, depth, true);
        near(
            model.operating_flux_linkages_wb()[0],
            flux + ends[0] * currents[0],
        );
        near(
            model.operating_flux_linkages_wb()[1],
            -3.0 * flux + ends[1] * currents[1],
        );
        let voltage = model.voltage_from_current_rate_v(&[0.7, -0.2]).unwrap();
        near(
            voltage[0],
            (active + ends[0]) * 0.7 + (-3.0 * active) * (-0.2),
        );
        near(
            voltage[1],
            (-3.0 * active) * 0.7 + (9.0 * active + ends[1]) * (-0.2),
        );
        // Independently perturb the full nonlinear solve, not the extracted L.
        for j in 0..2 {
            let epsilon = 1e-5;
            let mut linkages = [[0.0; 2]; 2];
            for (side, sign) in [-1.0, 1.0].iter().enumerate() {
                let mut i = currents;
                i[j] += sign * epsilon;
                let density: Vec<_> = (0..8)
                    .map(|cell| profiles[0][cell] * i[0] + profiles[1][cell] * i[1])
                    .collect();
                let field = problem
                    .solve(
                        &density,
                        &common::boundary(),
                        model.operating_field().potential_wb_m(),
                        common::material::controls().nonlinear,
                        &mut || false,
                    )
                    .unwrap();
                for winding in 0..2 {
                    let mut linkage = ends[winding] * i[winding];
                    for (cell, face) in problem.mesh().faces().iter().enumerate() {
                        let sum: f64 = face
                            .iter()
                            .map(|&v| field.potential_wb_m()[v as usize])
                            .sum();
                        linkage += sum / 3.0
                            * profiles[winding][cell]
                            * problem.mesh().face_measure(cell).unwrap();
                    }
                    linkages[side][winding] = linkage;
                }
            }
            for receiver in 0..2 {
                near(
                    (linkages[1][receiver] - linkages[0][receiver]) / (2.0 * epsilon),
                    model.inductance_h()[receiver * 2 + j],
                );
            }
        }
        let delta_i = [0.001, -0.002];
        let delta_lambda = model.linearized_flux_increment_wb(&delta_i).unwrap();
        let recovered = model.linearized_current_increment_a(&delta_lambda).unwrap();
        for i in 0..2 {
            near(recovered[i], delta_i[i]);
        }
    }
}

#[test]
fn biased_tangent_retains_transverse_reluctivity_and_homogeneous_incremental_traces() {
    let depth = 0.8;
    let source = common::source(true, 2.0);
    let problem = common::problem(&source, depth);
    let bias_t = 1.5;
    for gauge in [0.0, 7.0] {
        let initial: Vec<_> = problem
            .mesh()
            .vertices()
            .iter()
            .map(|p| gauge + bias_t * p[0])
            .collect();
        let boundary: Vec<_> = (0..6).map(|v| (v, initial[v])).collect();
        let model = PlanarIncrementalWindingBasis::new(
            &problem,
            common::profiles(1),
            &boundary,
            vec![0.0],
            vec![0.0],
            &initial,
            config(),
            common::material::controls().nonlinear,
            &mut || false,
        )
        .unwrap();
        // The base uniform B has H/B=2 and dH/dB=4. Each lobe's center
        // couples two radial and two transverse triangles: K_t=12*depth,
        // load=2*depth, dA_center/di=1/6, L=4*depth/6.
        near(model.inductance_h()[0], 2.0 * depth / 3.0);
        near(model.operating_flux_linkages_wb()[0], -6.0 * depth * bias_t);
        let response = &model.unit_responses()[0];
        for v in 0..6 {
            assert_eq!(response.potential_per_current_h_m()[v], 0.0);
        }
        near(response.potential_per_current_h_m()[6], 1.0 / 6.0);
        near(response.potential_per_current_h_m()[7], -1.0 / 6.0);
        for b in model.operating_field().magnetic_flux_density_t() {
            near(b[0], 0.0);
            near(b[1], -bias_t);
        }
        assert_eq!(model.boundary_potential_wb_m(), boundary);
        assert_eq!(model.net_turns_by_component(), &[vec![0.0]]);
    }
}

#[test]
fn source_corners_and_range_edges_cannot_masquerade_as_unique_differential_inductance() {
    let source = common::source(true, 2.0);
    let problem = common::problem(&source, 0.7);
    for magnitude in [1.0, 4.0] {
        // An exact uniform field meets its boundary and has zero volume residual.
        let initial: Vec<_> = problem
            .mesh()
            .vertices()
            .iter()
            .map(|p| magnitude * p[0])
            .collect();
        let boundary: Vec<_> = (0..6).map(|v| (v, initial[v])).collect();
        let result = PlanarIncrementalWindingBasis::new(
            &problem,
            common::profiles(1),
            &boundary,
            vec![0.0],
            vec![0.01],
            &initial,
            config(),
            common::material::controls().nonlinear,
            &mut || false,
        );
        if magnitude == 1.0 {
            assert!(matches!(
                result,
                Err(IncrementalWindingError::Nondifferentiable { .. })
            ));
        } else {
            assert!(matches!(
                result,
                Err(IncrementalWindingError::SourceBoundary { .. })
            ));
        }
    }
    // A source data knot with identical adjacent slopes is smooth and valid.
    let smooth = common::material::source(&[(0.0, 0.0), (1.0, 1.0), (4.0, 4.0)], 2.0);
    let linear = common::problem(&smooth, 0.7);
    near(basis(&linear, 1.0, 0.0).inductance_h()[0], 1.4);
}

#[test]
fn extraction_refuses_unclosed_profiles_singular_modes_unconverged_bias_and_cancelled_work() {
    let source = common::source(true, 2.0);
    let problem = common::problem(&source, 0.7);
    let run = |profiles, ends, cfg, controls, cancelled: &mut dyn FnMut() -> bool| {
        PlanarIncrementalWindingBasis::new(
            &problem,
            profiles,
            &common::boundary(),
            vec![2.0],
            ends,
            &[0.0; 8],
            cfg,
            controls,
            &mut || cancelled(),
        )
    };
    let controls = common::material::controls().nonlinear;
    assert!(matches!(
        run(
            vec![vec![6.0; 8]],
            vec![0.01],
            config(),
            controls,
            &mut || false
        ),
        Err(IncrementalWindingError::Winding(WindingError::Invalid(_)))
    ));
    let mut tiny = config();
    tiny.winding.max_basis_values = 1;
    assert!(matches!(
        run(common::profiles(1), vec![0.01], tiny, controls, &mut || {
            false
        }),
        Err(IncrementalWindingError::Winding(WindingError::Invalid(_)))
    ));
    let mut no_newton = controls;
    no_newton.max_newton_iterations = 0;
    assert!(matches!(
        run(
            common::profiles(1),
            vec![0.01],
            config(),
            no_newton,
            &mut || false
        ),
        Err(IncrementalWindingError::Magnetic(
            NonlinearMagneticError::Nonconvergence(_)
        ))
    ));
    let duplicate = common::profiles(1).remove(0);
    let singular = PlanarIncrementalWindingBasis::new(
        &problem,
        vec![duplicate.clone(), duplicate],
        &common::boundary(),
        vec![1.0, 1.0],
        vec![0.0; 2],
        &[0.0; 8],
        config(),
        controls,
        &mut || false,
    );
    assert!(matches!(
        singular,
        Err(IncrementalWindingError::Winding(
            WindingError::SingularInductance
        ))
    ));
    let mut checks = 0;
    run(
        common::profiles(1),
        vec![0.01],
        config(),
        controls,
        &mut || {
            checks += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    let interrupted = run(
        common::profiles(1),
        vec![0.01],
        config(),
        controls,
        &mut || {
            seen += 1;
            seen >= checks - 1
        },
    );
    assert!(matches!(
        interrupted,
        Err(IncrementalWindingError::Magnetic(
            NonlinearMagneticError::Magnetic(MagneticError::Cancelled)
        )) | Err(IncrementalWindingError::Winding(WindingError::Magnetic(
            MagneticError::Cancelled
        )))
    ));
    let accepted = basis(&problem, 2.0, 0.01);
    assert!(accepted.voltage_from_current_rate_v(&[f64::NAN]).is_err());
    assert!(accepted.linearized_current_increment_a(&[]).is_err());
    near(accepted.inductance_h()[0], 0.36);
}
