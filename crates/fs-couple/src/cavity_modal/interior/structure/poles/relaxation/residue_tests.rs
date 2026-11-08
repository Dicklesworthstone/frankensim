//! Full FE adapter checks against an independently factored coupled cubic.

use super::*;
use crate::cavity_modal::interior::structure::poles::relaxation::{
    MaxwellDampingBranch, StructuralInteriorRelaxationError,
    residue::{
        MaxwellBranchDirection, RelaxationResidueOptions, StructuralInteriorRelaxationDirection,
        StructuralInteriorRelaxationLoad,
    },
};
use fs_modal::relaxation_modes::MaxwellModeError;

fn fixture(f: impl FnOnce(&StructuralInterior<'_, '_, '_>, &[MaxwellDampingBranch], &Cx<'_>)) {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 0.2,
        sound_speed_m_s: 1.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let mass = dense(&[1.0], 1);
    let stiffness = dense(&[0.7], 1);
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
            &[dilation(&geometry)],
            &[],
            controls(),
            cx,
        )
        .unwrap();
        f(&model, &branches, cx);
    });
}

fn roots() -> [C64; 3] {
    [
        C64::from_re(-0.5),
        C64::new(-0.25, -31.0_f64.sqrt() / 4.0),
        C64::new(-0.25, 31.0_f64.sqrt() / 4.0),
    ]
}

#[test]
fn force_and_inward_flow_residues_recover_every_original_field() {
    // Elimination is used ONLY for this independent oracle, never by the
    // implementation. T(s)=s²+1+1.25*s/(1+s), and source F+0.6*Q/s.
    fixture(|model, branches, cx| {
        let prepared = model
            .prepare_relaxation_modes(branches, options(), cx)
            .unwrap();
        let force = [C64::new(0.7, -0.3)];
        let inward = [C64::new(0.2, 0.1); 4];
        let zero_force = [C64::ZERO];
        let zero_flow = [C64::ZERO; 4];
        let loads = [
            StructuralInteriorRelaxationLoad {
                structural_force: &force,
                inward_flow_m3_s: &zero_flow,
            },
            StructuralInteriorRelaxationLoad {
                structural_force: &zero_force,
                inward_flow_m3_s: &inward,
            },
            StructuralInteriorRelaxationLoad {
                structural_force: &force,
                inward_flow_m3_s: &inward,
            },
        ];
        for s in roots() {
            let right = C64::from_re(-0.6) / s;
            let pair = prepared
                .refine_residue(
                    s.scale(1.02),
                    &[C64::ONE],
                    &[right; 4],
                    &[C64::ONE],
                    &[-right; 4],
                    RelaxationResidueOptions::default(),
                    cx,
                )
                .unwrap();
            close(pair.mode().laplace_pole_per_s, s, 3e-11);
            let transpose = pair.transpose_displacement_potential();
            close(transpose[1] / transpose[0], -right, 3e-11);
            let residuals = pair.transpose_source_residuals();
            assert!(residuals.0 < 2e-12 && residuals.1 < 2e-12);
            let results = pair.force_residues(&loads, cx).unwrap();
            let slope = s.scale(2.0) + C64::from_re(1.25) / ((C64::ONE + s) * (C64::ONE + s));
            for (load, residue) in loads.iter().zip(&results) {
                let flow = load
                    .inward_flow_m3_s
                    .iter()
                    .copied()
                    .fold(C64::ZERO, |a, b| a + b);
                let expected = (load.structural_force[0] + flow.scale(0.6) / s) / slope;
                close(residue.structural_displacement_residue[0], expected, 3e-11);
                close(residue.structural_velocity_residue[0], s * expected, 3e-11);
                for a in 0..4 {
                    close(
                        residue.pressure_potential_residue_pa[a],
                        right * expected,
                        3e-11,
                    );
                    close(
                        residue.pressure_residue_pa_per_s[a],
                        expected.scale(-0.6),
                        3e-11,
                    );
                    close(
                        residue.displaced_flow_residue_m3_s2[a],
                        (s * expected).scale(0.125),
                        3e-11,
                    );
                }
                close(residue.fluid_force_residue[0], expected.scale(-0.3), 3e-11);
                close(
                    residue.branch_stress_residue[0],
                    expected * s.scale(1.25) / (C64::ONE + s),
                    3e-11,
                );
                close(
                    pair.mode().structural_displacement[0] * residue.participation,
                    expected,
                    3e-11,
                );
            }
            close(
                results[0].structural_displacement_residue[0]
                    + results[1].structural_displacement_residue[0],
                results[2].structural_displacement_residue[0],
                3e-12,
            );
        }
    });
}

#[test]
fn original_matrix_and_material_directions_give_exact_cubic_pole_slopes() {
    fixture(|model, branches, cx| {
        let prepared = model
            .prepare_relaxation_modes(branches, options(), cx)
            .unwrap();
        let zero = [0.0];
        let one = [1.0];
        let zero_acoustic = [0.0; 16];
        let zero_coupling = [0.0; 4];
        let zero_strain = [0.0; 5];
        let unit_strain = [1.0, 0.0, 0.0, 0.0, 0.0];
        let zero_branch = [MaxwellBranchDirection {
            stiffness: 0.0,
            log_relaxation_time: 0.0,
            strain: &zero_strain,
        }];
        let stiffness_branch = [MaxwellBranchDirection {
            stiffness: 1.0,
            ..zero_branch[0]
        }];
        let time_branch = [MaxwellBranchDirection {
            log_relaxation_time: 1.0,
            ..zero_branch[0]
        }];
        let strain_branch = [MaxwellBranchDirection {
            strain: &unit_strain,
            ..zero_branch[0]
        }];
        let spatial = StructuralInteriorPoleDirection {
            structural_mass: &zero,
            structural_damping: &zero,
            structural_stiffness: &zero,
            acoustic_mass: &zero_acoustic,
            acoustic_stiffness: &zero_acoustic,
            wall_damping: &zero_acoustic,
            coupling: &zero_coupling,
        };
        let base = StructuralInteriorRelaxationDirection {
            spatial,
            branches: &zero_branch,
        };
        let directions = [
            StructuralInteriorRelaxationDirection {
                spatial: StructuralInteriorPoleDirection {
                    structural_mass: &one,
                    ..spatial
                },
                ..base
            },
            StructuralInteriorRelaxationDirection {
                spatial: StructuralInteriorPoleDirection {
                    structural_damping: &one,
                    ..spatial
                },
                ..base
            },
            StructuralInteriorRelaxationDirection {
                spatial: StructuralInteriorPoleDirection {
                    structural_stiffness: &one,
                    ..spatial
                },
                ..base
            },
            StructuralInteriorRelaxationDirection {
                spatial: StructuralInteriorPoleDirection {
                    coupling: &model.coupling,
                    ..spatial
                },
                ..base
            },
            StructuralInteriorRelaxationDirection {
                branches: &stiffness_branch,
                ..base
            },
            StructuralInteriorRelaxationDirection {
                branches: &time_branch,
                ..base
            },
            StructuralInteriorRelaxationDirection {
                branches: &strain_branch,
                ..base
            },
            base,
        ];
        for s in roots() {
            let right = C64::from_re(-0.6) / s;
            let pair = prepared
                .refine_residue(
                    s,
                    &[C64::ONE],
                    &[right; 4],
                    &[C64::ONE],
                    &[-right; 4],
                    RelaxationResidueOptions::default(),
                    cx,
                )
                .unwrap();
            let slopes = pair.pole_derivatives(&directions, cx).unwrap();
            let denominator = C64::ONE + s;
            let ts = s.scale(2.0) + C64::from_re(1.25) / (denominator * denominator);
            let changes = [
                s * s,
                s,
                C64::ONE,
                C64::from_re(0.6),
                s / denominator,
                s.scale(1.25) / (denominator * denominator),
                s.scale(2.5) / denominator,
                C64::ZERO,
            ];
            for (actual, change) in slopes.iter().zip(changes) {
                let expected = -change / ts;
                close(actual.laplace_pole_derivative_per_s, expected, 3e-11);
                assert!((actual.amplitude_decay_rate_derivative_per_s + expected.re).abs() < 3e-11);
                if s.im == 0.0 {
                    assert_eq!(actual.oscillation_frequency_derivative_rad_s, None);
                } else {
                    assert!(
                        (actual.oscillation_frequency_derivative_rad_s.unwrap()
                            - s.im.signum() * expected.im)
                            .abs()
                            < 3e-11
                    );
                }
            }
        }
    });
}

#[test]
fn residue_budget_bad_batch_and_cancellation_preserve_the_reusable_pair() {
    fixture(|model, branches, cx| {
        let prepared = model
            .prepare_relaxation_modes(branches, options(), cx)
            .unwrap();
        let s = roots()[0];
        let right = C64::from_re(-0.6) / s;
        let low = RelaxationResidueOptions {
            max_physical_bytes: 1,
            ..RelaxationResidueOptions::default()
        };
        assert!(matches!(
            prepared.refine_residue(
                s,
                &[C64::ONE],
                &[right; 4],
                &[C64::ONE],
                &[-right; 4],
                low,
                cx
            ),
            Err(StructuralInteriorRelaxationError::Mode(
                MaxwellModeError::Budget { .. }
            ))
        ));
        let limits = RelaxationResidueOptions {
            max_loads: 1,
            max_directions: 1,
            ..RelaxationResidueOptions::default()
        };
        let pair = prepared
            .refine_residue(
                s,
                &[C64::ONE],
                &[right; 4],
                &[C64::ONE],
                &[-right; 4],
                limits,
                cx,
            )
            .unwrap();
        let force = [C64::ONE];
        let flow = [C64::ZERO; 4];
        let load = StructuralInteriorRelaxationLoad {
            structural_force: &force,
            inward_flow_m3_s: &flow,
        };
        let baseline = pair.force_residues(&[load], cx).unwrap();
        assert!(pair.force_residues(&[load, load], cx).is_err());
        assert!(pair.force_residues(&[], cx).is_err());
        assert!(
            pair.force_residues(
                &[StructuralInteriorRelaxationLoad {
                    structural_force: &[],
                    ..load
                }],
                cx
            )
            .is_err()
        );
        assert!(pair.pole_derivatives(&[], cx).is_err());
        let gate = CancelGate::new();
        gate.request();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0,
                    kernel_id: 0,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(pair.force_residues(&[load], &cancelled).is_err());
        });
        let retry = pair.force_residues(&[load], cx).unwrap();
        assert_eq!(
            baseline[0].structural_displacement_residue,
            retry[0].structural_displacement_residue
        );
        assert_eq!(
            baseline[0].branch_stress_residue,
            retry[0].branch_stress_residue
        );
    });
}
