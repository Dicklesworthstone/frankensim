//! Offset and oriented layer physics in the original full-coordinate shell.

use fs_exec::{CancelGate, Cx};
use fs_material::elastic::OrthotropicElastic;
use fs_material::laminate::{
    LaminaLayer, LaminateError, LaminateSection, LayerLaw, ReferenceSurface,
};
use fs_material::visco::tensor::{PronyBranch, TensorProny};
use fs_math::c64::C64;
use fs_modal::harmonic::{HarmonicParameterControls, HarmonicResponseControls};
use fs_plate::dkt_curvature_quadrature;
use fs_plate::shell::viscosity::ShellViscosityBudget;
use fs_plate::shell::viscosity::harmonic::{
    TensorPronyShell, TensorPronyShellSection, TensorShellHarmonicError,
};
use fs_plate::shell::{CstDktMass, ShellMesh, ShellSupport};

fn with_cx(run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 59,
                kernel_id: 23,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ));
    });
}

fn elastic(limit: f64) -> OrthotropicElastic {
    OrthotropicElastic::new(
        [14.0, 7.0, 10.0],
        [0.21, 0.27, 0.18],
        [3.0, 2.0, 4.0],
        limit,
    )
    .unwrap()
}
fn prony() -> TensorProny {
    let c = elastic(0.03).stiffness();
    // Different shear/normal branch proportions make the reduced tensor
    // frequency dependent beyond one scalar factor, even before layer rotation.
    let branch = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            if i == j {
                [6.0, 1.0, 4.0, 2.0, 3.0, 5.0][i]
            } else {
                0.0
            }
        })
    });
    TensorProny::new_mandel(c, vec![PronyBranch::new(branch, 0.3).unwrap()], 0.03)
        .unwrap()
        .with_band((0.0, 20.0))
        .unwrap()
}
fn layer(law: LayerLaw, h: f64, rho: f64, angle: f64) -> LaminaLayer {
    LaminaLayer {
        thickness_m: h,
        orientation_rad: angle,
        density_kg_m3: rho,
        material: None,
        law,
    }
}
fn laminate() -> LaminateSection {
    LaminateSection::try_new(
        vec![
            layer(LayerLaw::Elastic(elastic(0.03)), 0.08, 4.0, 0.31),
            layer(LayerLaw::Prony(prony()), 0.02, 9.0, -0.47),
        ],
        ReferenceSurface::Offset(0.03),
    )
    .unwrap()
}
fn budget() -> ShellViscosityBudget {
    ShellViscosityBudget {
        max_elements: 4,
        max_full_dofs: 36,
        max_dofs: 36,
        max_modes: 18,
        max_work: 20_000_000,
    }
}
fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 36,
        maximum_payload_bytes: 1 << 22,
        maximum_relative_residual: 2e-11,
    }
}
fn source(
    section: LaminateSection,
    mass: CstDktMass,
    cx: &Cx<'_>,
) -> Result<TensorPronyShell, TensorShellHarmonicError> {
    TensorPronyShell::assemble(
        ShellMesh::new(
            vec![[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [0.0, 0.3, 0.0]],
            vec![[0, 1, 2]],
        )
        .unwrap(),
        vec![TensorPronyShellSection::from_laminate(section).unwrap()],
        &[0, 1],
        ShellSupport::Clamped,
        mass,
        0.0,
        budget(),
        cx,
    )
}
fn real_near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-14),
        "{actual:e} != {expected:e}"
    );
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-8),
        "{actual:?} != {expected:?}"
    );
}
fn multiply(matrix: &fs_sparse::Csr, q: &[C64]) -> Vec<C64> {
    (0..q.len())
        .map(|i| (0..q.len()).fold(C64::ZERO, |sum, j| sum + q[j].scale(matrix.get(i, j))))
        .collect()
}

// Independent virtual-work evaluation of the retained CLT section on the
// public DKT curvature rule. The free node is (0, 0.3); its CST gradients are
// (0, 1/0.3). Physical theta_x -> beta_y; theta_y -> -beta_x.
fn physical_stiffness(section: &LaminateSection, omega: f64) -> [[C64; 6]; 6] {
    let (a, b, d) = section.complex_abd(omega).unwrap();
    let (curvature, _) = dkt_curvature_quadrature(&[0.0, 0.4, 0.0], &[0.0, 0.0, 0.3], 0).unwrap();
    let mut k = [[C64::ZERO; 6]; 6];
    for point in curvature {
        let mut strain = [[0.0; 6]; 6];
        strain[1][1] = 1.0 / 0.3;
        strain[2][0] = 1.0 / 0.3;
        for r in 0..3 {
            strain[r + 3][2] = point[r][6];
            strain[r + 3][3] = point[r][8];
            strain[r + 3][4] = -point[r][7];
        }
        for i in 0..6 {
            for j in 0..6 {
                for r in 0..6 {
                    for s in 0..6 {
                        let constitutive = match (r < 3, s < 3) {
                            (true, true) => a[3 * r + s],
                            (false, false) => d[3 * (r - 3) + s - 3],
                            _ => b[3 * (r % 3) + s % 3],
                        };
                        k[i][j] = k[i][j] + constitutive.scale(0.02 * strain[r][i] * strain[s][j]);
                    }
                }
            }
        }
    }
    k
}

#[test]
fn reference_pencil_retains_rotated_abd_and_signed_physical_mass_moments() {
    with_cx(|cx| {
        let section = laminate();
        let original = section.clone();
        let shell = source(section, CstDktMass::Consistent, cx).unwrap();
        assert_eq!(shell.sections()[0].laminate(), Some(&original));
        assert!(shell.sections()[0].law().is_none());
        let k = physical_stiffness(&original, 0.0);
        // The v row has zero numerical drilling strain. Its rotation coupling
        // is entirely physical B; centered homogeneous substitution would zero it.
        for j in 0..6 {
            real_near(shell.model().k.get(1, j), k[1][j].re, 3e-12);
        }
        assert!(k[1][3].abs() > 1e-3);
        let i0 = 4.0 * 0.08 + 9.0 * 0.02;
        let i1 = 4.0 * 0.08 * 0.01 + 9.0 * 0.02 * 0.06;
        let i2 = (4.0 * (0.05_f64.powi(3) + 0.03_f64.powi(3))
            + 9.0 * (0.07_f64.powi(3) - 0.05_f64.powi(3)))
            / 3.0;
        // Integral N_2^2 dA = area/6 = 0.01. The first moment couples
        // u to +theta_y and v to -theta_x under the original director motion.
        real_near(shell.model().m.get(2, 2), 0.01 * i0, 3e-13);
        real_near(shell.model().m.get(0, 4), 0.01 * i1, 3e-13);
        real_near(shell.model().m.get(1, 3), -0.01 * i1, 3e-13);
        real_near(shell.model().m.get(3, 3), 0.01 * i2, 3e-13);
        assert!(source(original, CstDktMass::Lumped, cx).is_err());
    });
}

#[test]
fn mixed_oriented_layer_motion_heat_and_analytic_frequency_tangent_agree() {
    with_cx(|cx| {
        let section = laminate();
        let shell = source(section.clone(), CstDktMass::Consistent, cx).unwrap();
        let q = [
            C64::ZERO,
            C64::new(2e-6, -1e-6),
            C64::new(3e-6, 8e-7),
            C64::new(-4e-6, 2e-6),
            C64::new(1e-6, -7e-7),
            C64::ZERO,
        ];
        for omega in [0.7, 2.4, 7.0] {
            let k = physical_stiffness(&section, omega);
            let mq = multiply(&shell.model().m, &q);
            let force: Vec<_> = (0..6)
                .map(|i| {
                    (0..6).fold(C64::ZERO, |sum, j| sum + k[i][j].conj() * q[j])
                        - mq[i].scale(omega * omega)
                })
                .collect();
            // This displacement has zero drilling strain: all material work
            // below is independently integrated layer loss, with no gauge work.
            let result = shell
                .prepare_harmonic_response(omega, controls(), cx)
                .unwrap()
                .with_frequency_derivative(1 << 20, cx)
                .unwrap()
                .solve(
                    &force,
                    None,
                    None,
                    HarmonicParameterControls {
                        maximum_additional_payload_bytes: 1 << 20,
                        maximum_work: 10_000_000,
                    },
                    cx,
                )
                .unwrap();
            for i in 0..6 {
                near(result.response.displacement[i], q[i], 2e-8);
            }
            let loss_work: f64 = (0..6)
                .map(|i| {
                    (0..6)
                        .map(|j| (q[i].conj() * q[j]).re * k[i][j].im)
                        .sum::<f64>()
                })
                .sum();
            assert!(loss_work > 0.0);
            real_near(
                result.response.average_viscous_power_w,
                0.5 * omega * loss_work,
                2e-10,
            );
            let h = 1e-5 * omega;
            let sample = |w| {
                shell
                    .prepare_harmonic_response(w, controls(), cx)
                    .unwrap()
                    .solve(&force, cx)
                    .unwrap()
            };
            let plus = sample(omega + h);
            let minus = sample(omega - h);
            for i in 0..6 {
                near(
                    result.displacement_derivative[i],
                    (plus.displacement[i] - minus.displacement[i]).scale(0.5 / h),
                    4e-6,
                );
            }
            real_near(
                result.average_viscous_power_derivative,
                (plus.average_viscous_power_w - minus.average_viscous_power_w) / (2.0 * h),
                4e-6,
            );
        }
        assert!(matches!(
            shell.prepare_harmonic_response(21.0, controls(), cx),
            Err(TensorShellHarmonicError::Laminate {
                element: 0,
                source: LaminateError::Layer { layer: 1, .. }
            })
        ));
    });
}

#[test]
fn original_interior_layer_strain_limit_is_enforced_before_response_publication() {
    with_cx(|cx| {
        let section = LaminateSection::try_new(
            vec![
                layer(LayerLaw::Elastic(elastic(0.02)), 0.04, 4.0, 0.0),
                layer(LayerLaw::Elastic(elastic(1e-5)), 0.02, 4.0, 0.41),
                layer(LayerLaw::Elastic(elastic(0.02)), 0.04, 4.0, 0.0),
            ],
            ReferenceSurface::Offset(0.0),
        )
        .unwrap();
        let shell = source(section, CstDktMass::Consistent, cx).unwrap();
        let q = [
            C64::ZERO,
            C64::new(3e-5, 2e-5),
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
        ];
        let kq = multiply(&shell.model().k, &q);
        let mq = multiply(&shell.model().m, &q);
        let force: Vec<_> = kq.iter().zip(mq).map(|(&k, m)| k - m.scale(4.0)).collect();
        assert!(matches!(
            shell
                .prepare_harmonic_response(2.0, controls(), cx)
                .unwrap()
                .solve(&force, cx),
            Err(TensorShellHarmonicError::Laminate {
                element: 0,
                source: LaminateError::Layer { layer: 1, .. }
            })
        ));
    });
}

#[test]
fn reduced_only_layer_cannot_invent_out_of_plane_strain_constants() {
    let section = LaminateSection::try_new(
        vec![layer(
            LayerLaw::PlaneStress {
                stiffness_pa: [4.0, 1.0, 0.0, 1.0, 4.0, 0.0, 0.0, 0.0, 2.0],
                transverse_shear_pa: [2.0, 3.0],
            },
            0.1,
            4.0,
            0.3,
        )],
        ReferenceSurface::Midsurface,
    )
    .unwrap();
    assert!(TensorPronyShellSection::from_laminate(section).is_err());
}
