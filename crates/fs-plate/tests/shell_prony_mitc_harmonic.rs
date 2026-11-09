//! Full source-owned membrane/bending/shear response on the native MITC field.

use fs_exec::{CancelGate, Cx};
use fs_material::elastic::OrthotropicElastic;
use fs_material::laminate::harmonic::{HarmonicLaminate, HarmonicLaminateBudget};
use fs_material::laminate::{
    LaminaLayer, LaminateSection, LaminateShearRecovery, LayerLaw, ReferenceSurface,
};
use fs_material::visco::{Mandel6, PronyBranch, TensorProny};
use fs_math::{c64::C64, det};
use fs_modal::harmonic::HarmonicResponseControls;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::connections::{
    ShellConnection, ShellConnectionBudget, ShellConnectionTerm,
};
use fs_plate::shell::mitc::harmonic::{HarmonicMitcShell, MitcHarmonicBudget, MitcHarmonicError};
use fs_plate::shell::mitc::strain::{PreparedElementStrain, ShellStrainPoint, integration_rule};
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, MitcShellModel, NodalDirectors, ShearShellSection,
    ShellConstraint, ShellMass, assemble_mitc3p,
};
use fs_sparse::Csr;
use std::sync::Arc;

fn with_cx(run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 59,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ));
    });
}
fn source_budget() -> HarmonicLaminateBudget {
    HarmonicLaminateBudget {
        maximum_layers: 8,
        maximum_payload_bytes: 1 << 20,
        maximum_work: 20_000_000,
    }
}
fn budget() -> MitcHarmonicBudget {
    MitcHarmonicBudget {
        maximum_elements: 4,
        maximum_layers: 16,
        maximum_payload_bytes: 1 << 22,
        maximum_work: 40_000_000,
    }
}
fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 36,
        maximum_payload_bytes: 1 << 22,
        maximum_relative_residual: 2e-11,
    }
}
fn elastic(scale: f64) -> OrthotropicElastic {
    OrthotropicElastic::new(
        [14e6 * scale, 7e6 * scale, 10e6 * scale],
        [0.21, 0.27, 0.18],
        [3e6 * scale, 2e6 * scale, 4e6 * scale],
        0.03,
    )
    .unwrap()
}
fn reference(angle: f64) -> Arc<LaminateSection> {
    Arc::new(
        LaminateSection::try_new_equilibrated_shear(
            vec![
                LaminaLayer {
                    law: LayerLaw::Elastic(elastic(1.0)),
                    thickness_m: 0.026,
                    density_kg_m3: 450.0,
                    orientation_rad: 0.31 + angle,
                    material: None,
                },
                LaminaLayer {
                    law: LayerLaw::Elastic(elastic(0.6)),
                    thickness_m: 0.002,
                    density_kg_m3: 1100.0,
                    orientation_rad: -0.47 + angle,
                    material: None,
                },
            ],
            ReferenceSurface::Offset(0.008),
            0.27 + angle,
        )
        .unwrap(),
    )
}
fn matrix(diagonal: [f64; 6], vector: [f64; 6], scale: f64) -> Mandel6 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            scale * (vector[i] * vector[j] + if i == j { diagonal[i] } else { 0.0 })
        })
    })
}
fn prony(scale: f64) -> TensorProny {
    let equilibrium = matrix(
        [17.0, 13.0, 11.0, 7.0, 5.0, 9.0],
        [1.3, -0.7, 0.9, 0.4, -0.5, 0.6],
        1e6 * scale,
    );
    let branch = matrix(
        [1.1, 0.9, 0.7, 0.5, 0.4, 0.3],
        [0.7, 0.4, -0.8, 0.9, -0.3, 0.5],
        1e6 * scale,
    );
    let kelvin = matrix(
        [0.7, 0.8, 0.9, 0.5, 0.3, 0.4],
        [0.3, 0.9, -0.4, 0.7, 0.8, -0.6],
        100.0 * scale,
    );
    TensorProny::new_mandel(
        equilibrium,
        vec![PronyBranch::new(branch, 0.014).unwrap()],
        0.03,
    )
    .unwrap()
    .with_kelvin_viscosity_mandel(kelvin)
    .unwrap()
    .with_band((0.0, 100.0))
    .unwrap()
}
fn actual(
    reference: &LaminateSection,
    density_scale: f64,
    elastic_scale: Option<f64>,
) -> Arc<LaminateSection> {
    Arc::new(
        LaminateSection::try_new(
            reference
                .layers()
                .iter()
                .enumerate()
                .map(|(i, layer)| LaminaLayer {
                    thickness_m: layer.thickness_m,
                    orientation_rad: layer.orientation_rad,
                    density_kg_m3: density_scale * layer.density_kg_m3,
                    material: None,
                    law: if let Some(scale) = elastic_scale {
                        LayerLaw::Elastic(elastic(scale * if i == 0 { 1.0 } else { 0.6 }))
                    } else {
                        LayerLaw::Prony(prony(if i == 0 { 1.0 } else { 0.6 }))
                    },
                })
                .collect(),
            ReferenceSurface::Offset(0.008),
        )
        .unwrap(),
    )
}
fn source(reference: Arc<LaminateSection>, actual: Arc<LaminateSection>) -> Arc<HarmonicLaminate> {
    let shear = Arc::new(LaminateShearRecovery::try_new(reference).unwrap());
    Arc::new(HarmonicLaminate::try_new(actual, shear, source_budget(), || false).unwrap())
}
fn model(reference: &LaminateSection, angle: f64) -> MitcShellModel {
    // Tilted facet, signed laminate first moments and oblique free-node axes.
    let mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [0.4, 0.03, 0.05], [0.07, 0.3, 0.04]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let normal = mesh.facet(0).unwrap().frame[2];
    let directors = NodalDirectors::new(&mesh, vec![normal; 3], 0.9).unwrap();
    let (s, c) = (det::sin(0.43), det::cos(0.43));
    let constraints = [
        ShellConstraint {
            node: 0,
            mask: [true; 6],
            frame: ConstraintFrame::Global,
        },
        ShellConstraint {
            node: 1,
            mask: [true; 6],
            frame: ConstraintFrame::Global,
        },
        ShellConstraint {
            node: 2,
            mask: [false; 6],
            frame: ConstraintFrame::Local([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]]),
        },
    ];
    assemble_mitc3p(
        &mesh,
        &directors,
        &[ShearShellSection::from_laminate(reference, angle).unwrap()],
        &constraints,
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: 3,
            max_triangles: 1,
        },
        0.08,
        || false,
    )
    .unwrap()
}
fn near(actual: C64, expected: C64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-14),
        "{actual:?} != {expected:?}"
    );
}
fn real_near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-18),
        "{actual:e} != {expected:e}"
    );
}
fn work(k: &Csr, c: &Csr, omega: f64, q: &[C64]) -> C64 {
    (0..q.len()).fold(C64::ZERO, |sum, i| {
        sum + q[i].conj()
            * (0..q.len()).fold(C64::ZERO, |row, j| {
                row + C64::new(k.get(i, j), -omega * c.get(i, j)) * q[j]
            })
    })
}
fn real_work(k: &Csr, q: &[C64]) -> f64 {
    (0..q.len())
        .map(|i| {
            (q[i].conj() * (0..q.len()).fold(C64::ZERO, |row, j| row + q[j].scale(k.get(i, j)))).re
        })
        .sum()
}
fn strain_in_section(point: ShellStrainPoint, angle: f64) -> [f64; 8] {
    let (s, c) = (det::sin(angle), det::cos(angle));
    let rotate = |[x, y, xy]: [f64; 3]| {
        [
            c * c * x + s * s * y + c * s * xy,
            s * s * x + c * c * y - c * s * xy,
            -2.0 * c * s * x + 2.0 * c * s * y + (c * c - s * s) * xy,
        ]
    };
    let e = rotate(point.membrane);
    let k = rotate(point.gradient_per_m);
    let [xz, yz] = point.transverse_shear;
    [
        e[0],
        e[1],
        e[2],
        k[0],
        k[1],
        k[2],
        c * xz + s * yz,
        -s * xz + c * yz,
    ]
}

// Independent virtual-work observation of an ACTUAL coherent nodal field,
// not the implementation's cached unit-column matrix contraction.
fn source_work(
    model: &MitcShellModel,
    material: &HarmonicLaminate,
    angle: f64,
    omega: f64,
    q: &[C64],
) -> (C64, C64) {
    let (mesh, directors, sections) = model.original_geometry();
    let native =
        PreparedElementStrain::new(mesh, 0, directors, &sections[0], &mut || false).unwrap();
    let real = model
        .expand(&q.iter().map(|v| v.re).collect::<Vec<_>>())
        .unwrap();
    let imaginary = model
        .expand(&q.iter().map(|v| v.im).collect::<Vec<_>>())
        .unwrap();
    let field = native
        .harmonic_field(
            &core::array::from_fn(|i| real[i]),
            &core::array::from_fn(|i| imaginary[i]),
            &mut || false,
        )
        .unwrap();
    let section = material
        .prepare_harmonic(omega, source_budget(), || false)
        .unwrap();
    let area = mesh.facet(0).unwrap().area_m2;
    let mut complete = C64::ZERO;
    let mut cross = C64::ZERO;
    for (lambda, weight) in integration_rule() {
        let point = field.at_barycentric(lambda).unwrap();
        let re = strain_in_section(point.real, angle);
        let im = strain_in_section(point.imaginary, angle);
        let e: [C64; 8] = core::array::from_fn(|i| C64::new(re[i], im[i]));
        let resultant = section.resultants(&e, || false).unwrap();
        for i in 0..8 {
            complete = complete + (e[i].conj() * resultant[i]).scale(area * weight);
            for j in 0..8 {
                if (i < 6) != (j < 6) {
                    cross = cross
                        + (e[i].conj()
                            * C64::new(section.storage()[8 * i + j], -section.loss()[8 * i + j])
                            * e[j])
                            .scale(area * weight);
                }
            }
        }
    }
    (complete, cross)
}
fn motion() -> Vec<C64> {
    vec![
        C64::new(8e-7, -3e-7),
        C64::new(-3e-7, 4e-7),
        C64::new(7e-7, 2e-7),
        C64::new(2e-5, -4e-6),
        C64::new(-8e-6, 3e-6),
        C64::new(4e-6, 2e-6),
    ]
}

#[test]
fn lossless_reference_preserves_original_mass_storage_and_physical_springs() {
    with_cx(|cx| {
        let reference = reference(0.0);
        let laminate = source(Arc::clone(&reference), Arc::clone(&reference));
        let base = model(&reference, 0.37);
        let term = ShellConnectionTerm {
            node: 2,
            translation: [0.4, -0.6, 0.2],
            rotation_m: [0.03, 0.01, -0.02],
        };
        let spring = 13_500.0;
        let attached = base
            .with_connections(
                &[ShellConnection {
                    terms: vec![term.clone()],
                    stiffness_n_m: spring,
                    damping_n_s_m: 0.0,
                }],
                ShellConnectionBudget {
                    max_connections: 1,
                    max_terms: 1,
                    max_work: 1_000_000,
                },
                || false,
            )
            .unwrap();
        let bare =
            HarmonicMitcShell::try_new(&base, vec![Arc::clone(&laminate)], budget(), cx).unwrap();
        let coupled = HarmonicMitcShell::try_new(&attached, vec![laminate], budget(), cx).unwrap();
        let mut physical_row = vec![0.0; 18];
        physical_row[12..15].copy_from_slice(&term.translation);
        physical_row[15..18].copy_from_slice(&term.rotation_m);
        let row = base.project_load(&physical_row).unwrap();
        let scale = (0..6)
            .flat_map(|i| (0..6).map(move |j| (i, j)))
            .map(|(i, j)| attached.pencil().k.get(i, j).abs())
            .fold(0.0_f64, f64::max);
        for omega in [17.0, 61.0] {
            let original = bare
                .prepare_harmonic_operator(omega, controls(), cx)
                .unwrap();
            let assembled = coupled
                .prepare_harmonic_operator(omega, controls(), cx)
                .unwrap();
            assert!(std::ptr::eq(assembled.mass(), &attached.pencil().m));
            assert_eq!(assembled.mass(), original.mass());
            assert_eq!(assembled.damping().nnz(), 0);
            for i in 0..6 {
                for j in 0..6 {
                    assert!(
                        (assembled.stiffness().get(i, j) - attached.pencil().k.get(i, j)).abs()
                            < 2e-10 * scale
                    );
                    assert!(
                        (assembled.stiffness().get(i, j)
                            - original.stiffness().get(i, j)
                            - spring * row[i] * row[j])
                            .abs()
                            < 1e-10 * spring
                    );
                }
            }
        }
    });
}

#[test]
fn full_anisotropic_shear_cross_work_survives_oblique_support_and_section_rotations() {
    with_cx(|cx| {
        let reference = reference(0.0);
        let reference_source = source(Arc::clone(&reference), Arc::clone(&reference));
        let actual = actual(&reference, 1.0, None);
        let laminate = source(Arc::clone(&reference), actual);
        let model = model(&reference, 0.37);
        let shell =
            HarmonicMitcShell::try_new(&model, vec![Arc::clone(&laminate)], budget(), cx).unwrap();
        assert!(core::ptr::eq(shell.laminate(0).unwrap(), laminate.as_ref()));
        let q = motion();
        let original_material = source_work(&model, &reference_source, 0.37, 37.0, &q).0.re;
        let fixed = real_work(&model.pencil().k, &q) - original_material;
        let mut works = Vec::new();
        for omega in [17.0, 61.0] {
            let assembled = shell
                .prepare_harmonic_operator(omega, controls(), cx)
                .unwrap();
            assembled.admit_displacement(&q, cx).unwrap();
            let (physical, cross) = source_work(&model, &laminate, 0.37, omega, &q);
            assert!(
                cross.abs() > 1e-12,
                "fixture must exercise full transverse cross work"
            );
            let actual = work(assembled.stiffness(), assembled.damping(), omega, &q);
            near(actual, physical + C64::from_re(fixed), 3e-10);
            assert!(actual.re > 0.0 && actual.im < 0.0);
            works.push(actual);
        }
        assert!((works[1] - works[0]).abs() > 1e-3 * works[0].abs());
    });
}

#[test]
fn full_native_factor_matches_force_work_and_refuses_out_of_source_motion() {
    with_cx(|cx| {
        let reference = reference(0.0);
        let laminate = source(Arc::clone(&reference), actual(&reference, 1.0, None));
        let model = model(&reference, 0.37);
        let shell = HarmonicMitcShell::try_new(&model, vec![laminate], budget(), cx).unwrap();
        let omega = 37.0;
        let operator = shell
            .prepare_harmonic_operator(omega, controls(), cx)
            .unwrap();
        let factor = shell
            .prepare_harmonic_response(omega, controls(), cx)
            .unwrap();
        let force = [
            C64::new(1e-5, -3e-6),
            C64::new(-2e-5, 4e-6),
            C64::new(3e-5, 5e-6),
            C64::new(-2e-7, 3e-8),
            C64::new(5e-7, -7e-8),
            C64::new(4e-8, 2e-8),
        ];
        let response = factor.solve(&force, cx).unwrap();
        assert!(response.relative_residual < controls().maximum_relative_residual);
        assert!(response.average_viscous_power_w > 0.0);
        let material = work(
            operator.stiffness(),
            operator.damping(),
            omega,
            &response.displacement,
        );
        real_near(
            response.average_viscous_power_w,
            -0.5 * omega * material.im,
            2e-9,
        );
        real_near(
            response.average_input_power_w,
            response.average_viscous_power_w,
            2e-8,
        );
        for i in 0..6 {
            let left = (0..6).fold(C64::ZERO, |sum, j| {
                sum + C64::new(
                    operator.stiffness().get(i, j) - omega * omega * operator.mass().get(i, j),
                    -omega * operator.damping().get(i, j),
                ) * response.displacement[j]
            });
            near(left, force[i], 2e-8);
        }
        assert!(matches!(
            operator.admit_displacement(&[C64::new(1.0, 1.0); 6], cx),
            Err(MitcHarmonicError::Material { .. })
        ));
        assert!(matches!(
            factor.solve(&force.map(|value| value.scale(1e12)), cx),
            Err(MitcHarmonicError::Material { .. })
        ));
        let zero = factor.solve(&[C64::ZERO; 6], cx).unwrap();
        assert!(zero.displacement.iter().all(|q| *q == C64::ZERO));
        assert_eq!(zero.average_viscous_power_w, 0.0);
    });
}

#[test]
fn original_material_bands_density_and_fixed_trial_binding_remain_mandatory() {
    with_cx(|cx| {
        let reference = reference(0.0);
        let laminate = source(Arc::clone(&reference), actual(&reference, 1.0, None));
        let model = model(&reference, 0.37);
        let shell =
            HarmonicMitcShell::try_new(&model, vec![Arc::clone(&laminate)], budget(), cx).unwrap();
        assert!(matches!(
            shell.prepare_harmonic_operator(101.0, controls(), cx),
            Err(MitcHarmonicError::Material { .. })
        ));
        assert!(
            shell
                .prepare_harmonic_operator(0.0, controls(), cx)
                .is_err()
        );
        let wrong_density = source(Arc::clone(&reference), actual(&reference, 2.0, None));
        assert!(matches!(
            HarmonicMitcShell::try_new(&model, vec![wrong_density], budget(), cx),
            Err(MitcHarmonicError::SourceMismatch { .. })
        ));
        let other = crate::reference(0.12);
        let wrong_trial = source(Arc::clone(&other), other);
        assert!(matches!(
            HarmonicMitcShell::try_new(&model, vec![wrong_trial], budget(), cx),
            Err(MitcHarmonicError::SourceMismatch { .. })
        ));
        assert!(matches!(
            HarmonicMitcShell::try_new(
                &model,
                vec![laminate],
                MitcHarmonicBudget {
                    maximum_work: 1,
                    ..budget()
                },
                cx
            ),
            Err(MitcHarmonicError::Budget { .. })
        ));
        assert!(matches!(
            shell.prepare_harmonic_operator(
                37.0,
                HarmonicResponseControls {
                    maximum_payload_bytes: 1,
                    ..controls()
                },
                cx
            ),
            Err(MitcHarmonicError::Budget { .. })
        ));
    });
}

#[test]
fn very_soft_actual_material_is_retained_without_elastic_reference_subtraction() {
    with_cx(|cx| {
        let reference = reference(0.0);
        let scale = 1e-20;
        let laminate = source(Arc::clone(&reference), actual(&reference, 1.0, Some(scale)));
        let mesh = ShellMesh::new(
            vec![[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [0.0, 0.3, 0.0]],
            vec![[0, 1, 2]],
        )
        .unwrap();
        let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 3], 0.9).unwrap();
        // One axial coordinate and fixed drilling rotations: no numerical
        // stiffness can conceal the twenty-decade actual material reduction.
        let constraints = [
            ShellConstraint {
                node: 0,
                mask: [true; 6],
                frame: ConstraintFrame::Global,
            },
            ShellConstraint {
                node: 1,
                mask: [true; 6],
                frame: ConstraintFrame::Global,
            },
            ShellConstraint {
                node: 2,
                mask: [true, false, true, true, true, true],
                frame: ConstraintFrame::Global,
            },
        ];
        let model = assemble_mitc3p(
            &mesh,
            &directors,
            &[ShearShellSection::from_laminate(&reference, 0.0).unwrap()],
            &constraints,
            ShellMass::Consistent,
            MitcBudget {
                max_nodes: 3,
                max_triangles: 1,
            },
            0.08,
            || false,
        )
        .unwrap();
        assert_eq!(model.pencil().free, 1);
        let shell = HarmonicMitcShell::try_new(&model, vec![laminate], budget(), cx).unwrap();
        let operator = shell
            .prepare_harmonic_operator(17.0, controls(), cx)
            .unwrap();
        let expected = scale * model.pencil().k.get(0, 0);
        assert!(operator.stiffness().get(0, 0) > 0.0);
        real_near(operator.stiffness().get(0, 0), expected, 5e-9);
        assert_eq!(operator.mass(), &model.pencil().m);
        assert_eq!(operator.damping().nnz(), 0);
    });
}
