//! Full source tensor through the native plane continuum and its reductions.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::{OrthotropicElastic, SmallStrainLaw, tensor};
use fs_solid::{
    Formulation, IncrementSettings, Mesh2, PlaneAssemblyBudget, PlaneKind, SmallStrainLawField,
    SmallStrainPlaneProblem,
};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 121,
                kernel_id: 121,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn principal() -> OrthotropicElastic {
    OrthotropicElastic::new([12e9, 4e9, 2e9], [0.2, 0.1, 0.25], [3e9, 1e9, 0.8e9], 0.02).unwrap()
}
fn material(limit: f64) -> ElasticAcousticState {
    let q = tensor::rotation([1.0, 1.0, 2.0].map(|v| v / 6.0_f64.sqrt()), 0.63);
    ElasticAcousticState::try_new_mandel(
        tensor::rotate_mandel(principal().stiffness(), q),
        1800.0,
        ContentHash([121; 32]),
        limit,
    )
    .unwrap()
}
fn close(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() <= 2e-10 * scale.max(a.abs()).max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn mandel(e: [f64; 6]) -> [f64; 6] {
    std::array::from_fn(|i| {
        if i < 3 {
            e[i]
        } else {
            core::f64::consts::SQRT_2 * e[i]
        }
    })
}

#[test]
fn g1_full_tensor_stress_tangent_and_work_share_the_original_mandel_law() {
    let law = material(0.02);
    let e = [0.001, -0.0004, 0.0002, 0.0003, -0.0002, 0.0001];
    let s = law.stress(&e, &());
    let c = law.tangent(&e, &());
    let response = law.strain_response(mandel(e)).unwrap();
    assert_eq!(
        response.material_state_identity,
        *law.material_state_identity()
    );
    close(
        response.stored_energy_density_j_m3,
        0.5 * tensor::contract(&e, &s),
        1.0,
    );
    for j in 0..6 {
        let mut plus = e;
        let mut minus = e;
        let h = 1e-6;
        plus[j] += h;
        minus[j] -= h;
        let sp = law.stress(&plus, &());
        let sm = law.stress(&minus, &());
        for i in 0..6 {
            close((sp[i] - sm[i]) / (2.0 * h), c[i][j], 1e9);
            let wi = if i < 3 { 1.0 } else { 2.0 };
            let wj = if j < 3 { 1.0 } else { 2.0 };
            close(wi * c[i][j], wj * c[j][i], 1e9);
        }
        let energy = |x| {
            law.strain_response(mandel(x))
                .unwrap()
                .stored_energy_density_j_m3
        };
        let weight = if j < 3 { 1.0 } else { 2.0 };
        close(
            (energy(plus) - energy(minus)) / (2.0 * h),
            weight * s[j],
            1e6,
        );
    }
}

#[test]
fn g3_full_tensor_material_and_strain_co_rotate_without_changing_elastic_work() {
    let original =
        ElasticAcousticState::from_orthotropic(&principal(), 1800.0, ContentHash([120; 32]))
            .unwrap();
    let q = tensor::rotation([1.0, 1.0, 2.0].map(|v| v / 6.0_f64.sqrt()), 0.63);
    let rotated = material(0.02);
    let e = [0.001, -0.0004, 0.0002, 0.0003, -0.0002, 0.0001];
    let expected = tensor::rotate(&original.stress(&e, &()), &q);
    let actual = rotated.stress(&tensor::rotate(&e, &q), &());
    for i in 0..6 {
        close(actual[i], expected[i], 1e6);
    }
    close(
        original
            .strain_response(mandel(e))
            .unwrap()
            .stored_energy_density_j_m3,
        rotated
            .strain_response(mandel(tensor::rotate(&e, &q)))
            .unwrap()
            .stored_energy_density_j_m3,
        1.0,
    );
}

fn reference(law: &ElasticAcousticState, plane: PlaneKind) -> ([f64; 6], [f64; 6], f64) {
    let c = law.stiffness_mandel_pa();
    let mut e = mandel([0.002, -0.0004, 0.0, 0.0003, 0.0, 0.0]);
    if plane == PlaneKind::Stress {
        // Independent SPD block solution in Mandel coordinates. The native
        // reduction instead performs local Newton/LU in tensor-Voigt space.
        let transverse = [2, 4, 5];
        let block: Vec<_> = transverse
            .iter()
            .flat_map(|&i| transverse.map(|j| c[i][j]))
            .collect();
        let factor = fs_la::factor::cholesky(&block, 3).unwrap();
        let mut rhs = transverse.map(|i| -(c[i][0] * e[0] + c[i][1] * e[1] + c[i][3] * e[3]));
        factor.solve(&mut rhs);
        for i in 0..3 {
            e[transverse[i]] = rhs[i];
        }
        assert!(
            e[4].abs() > 1e-6 && e[5].abs() > 1e-6,
            "exercise transverse shear relaxation"
        );
    }
    let s: [f64; 6] = std::array::from_fn(|i| c[i].iter().zip(e).map(|(a, b)| a * b).sum());
    let energy = 0.5 * e.iter().zip(s).map(|(a, b)| a * b).sum::<f64>();
    let to_voigt = |v: [f64; 6]| {
        std::array::from_fn(|i| {
            if i < 3 {
                v[i]
            } else {
                core::f64::consts::FRAC_1_SQRT_2 * v[i]
            }
        })
    };
    (to_voigt(e), to_voigt(s), energy)
}
fn affine(p: [f64; 2]) -> [f64; 2] {
    [0.002 * p[0] + 0.0003 * p[1], 0.0003 * p[0] - 0.0004 * p[1]]
}

#[test]
fn g1_rotated_anisotropic_tensor_drives_native_p1_q1_plane_stress_and_strain() {
    let law = material(0.02);
    with_cx(|cx| {
        for mesh in [
            Mesh2::quads(2.0, 1.0, 2, 2),
            Mesh2::triangles(2.0, 1.0, 2, 2),
        ] {
            let mut fixed = Vec::new();
            for (node, &p) in mesh.nodes.iter().enumerate() {
                if p[0] == 0.0 || p[0] == 2.0 || p[1] == 0.0 || p[1] == 1.0 {
                    let u = affine(p);
                    fixed.extend([(2 * node, u[0]), (2 * node + 1, u[1])]);
                }
            }
            assert!(
                fixed.len() < 2 * mesh.nodes.len(),
                "leave internal DOFs to solve"
            );
            for plane in [PlaneKind::Stress, PlaneKind::Strain] {
                let problem = SmallStrainPlaneProblem {
                    mesh: &mesh,
                    laws: SmallStrainLawField::Uniform(&law),
                    plane,
                    formulation: Formulation::Standard,
                    thickness_m: 0.2,
                    reduction: ReductionConfig {
                        max_iterations: 8,
                        stress_tolerance_pa: 1e-6,
                        strain_limit: 0.02,
                    },
                    rotation_limit_rad: 0.02,
                    prescribed_m: &fixed,
                    nodal_forces_n: &[],
                    load_path: &[1.0],
                    budget: PlaneAssemblyBudget::default(),
                    settings: IncrementSettings::default(),
                };
                let solved = problem.solve(cx).unwrap();
                let (expected_e, expected_s, energy_density) = reference(&law, plane);
                for (&p, u) in mesh.nodes.iter().zip(&solved.displacement_m) {
                    let expected = affine(p);
                    for i in 0..2 {
                        close(u[i], expected[i], 0.01);
                    }
                }
                for (state, stress) in solved
                    .quadrature_states
                    .iter()
                    .zip(&solved.quadrature_stress_pa)
                {
                    for i in 0..6 {
                        close(state.strain[i], expected_e[i], 0.01);
                        close(stress[i], expected_s[i], 1e7);
                    }
                }
                let u: Vec<f64> = solved.displacement_m.iter().flatten().copied().collect();
                let k = problem
                    .tangent_matrix(&u, &solved.quadrature_states, cx)
                    .unwrap();
                let dense = k.to_dense();
                let n = u.len();
                let energy: f64 = (0..n)
                    .map(|i| {
                        (0..n)
                            .map(|j| 0.5 * u[i] * dense[n * i + j] * u[j])
                            .sum::<f64>()
                    })
                    .sum();
                close(energy, 2.0 * 1.0 * 0.2 * energy_density, 1.0);
                for i in 0..n {
                    for j in 0..n {
                        close(dense[n * i + j], dense[n * j + i], 1e9);
                    }
                }
            }
        }
    });
}

#[test]
fn g0_original_tensor_strain_domain_survives_the_native_law_adapter() {
    let law = material(1e-3);
    assert!(
        law.strain_response([0.002, 0.0, 0.0, 0.0, 0.0, 0.0])
            .is_err()
    );
    let e = [0.0, 0.0, 0.0, 0.0008, 0.0, 0.0]; // sqrt(2)*epsilon_xy exceeds ceiling
    assert!(law.stress(&e, &()).iter().all(|v| v.is_nan()));
    assert!(law.tangent(&e, &()).iter().flatten().all(|v| v.is_nan()));
    assert!(law.strain_response([f64::NAN; 6]).is_err());
    assert_eq!(
        law.strain_response([0.0; 6])
            .unwrap()
            .stored_energy_density_j_m3,
        0.0
    );
}
