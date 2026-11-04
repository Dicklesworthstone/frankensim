//! Original quadrature thermal gradients, compatible bending and memory work.
//! Every constitutive number is an explicit synthetic software fixture.

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::Cx;
use fs_feec::highorder::quad1d::gauss_legendre;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, MaterialPropertySelection,
    ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES,
};
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_solid::linear::Formulation;
use fs_solid::linear3::hex::{
    HEX20_REFERENCE_NODES, HEX27_REFERENCE_NODES, HexConnectivity, HexElasticProblem, HexFamily,
    HexQuadrature,
};
use fs_solid::linear3::source::{
    NodalThermalTemperatures, QuadraticThermalSource, QuadraticThermalSources,
    ThermalExpansionModel,
};
use fs_solid::linear3::tet10::{QuadraticTetProblem, TetQuadrature, promote_p1_to_p2};
use fs_solid::linear3::visco::*;
use fs_solid::linear3::{
    QuadraticThermalStrainField, TetAssemblyBudget, TetElasticError, TetMaterialField,
    TetStaticSolveConfig, TetThermalStrainField, TetThermalStrainState,
};

const FRAME: ContentHash = ContentHash([97; 32]);
const HEX_NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
];
const HEXES: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];

fn with_cx(f: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate)) {
    let gate = fs_exec::CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 927,
                    kernel_id: 927,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance={tolerance:e}"
    );
}
fn relative(a: f64, b: f64) {
    near(a, b, 3e-10 * a.abs().max(b.abs()).max(1e-18));
}
fn diagonal(d: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { d[i] } else { 0.0 }))
}
fn quadratic(a: [f64; 6], c: &Mandel6, b: [f64; 6]) -> f64 {
    (0..6)
        .map(|i| a[i] * (0..6).map(|j| c[i][j] * b[j]).sum::<f64>())
        .sum()
}
fn tensors() -> (Mandel6, Mandel6, Mandel6) {
    let mut c = diagonal([90.0, 70.0, 50.0, 60.0, 30.0, 40.0]);
    c[0][3] = -8.0;
    c[3][0] = -8.0;
    c[1][3] = 5.0;
    c[3][1] = 5.0;
    let v = [2.0, -1.0, 1.0, 3.0, 0.0, -1.0];
    let branch = core::array::from_fn(|i| core::array::from_fn(|j| 4.0 * v[i] * v[j]));
    (c, branch, diagonal([0.3, 0.2, 0.1, 0.5, 0.2, 0.3]))
}
fn material(memory: bool) -> TetViscoMaterial {
    let (c, branch, eta) = tensors();
    let mut law = TensorProny::new_mandel(
        c,
        if memory {
            vec![PronyBranch::new(branch, 0.2).unwrap()]
        } else {
            vec![]
        },
        0.025,
    )
    .unwrap();
    if memory {
        law = law.with_kelvin_viscosity_mandel(eta).unwrap();
    }
    TetViscoMaterial {
        law,
        density_kg_m3: 27.0,
        material_state_identity: [98; 32],
    }
}
fn config() -> TetStaticSolveConfig {
    TetStaticSolveConfig {
        maximum_iterations: 400,
        relative_residual_tolerance: 2e-11,
    }
}
fn prepare(
    mesh: SolidMeshRef<'_>,
    m: &TetViscoMaterial,
    fixed: &[usize],
    cx: &Cx<'_>,
) -> ViscoElements {
    ViscoElements::new(
        mesh,
        TetViscoField::Uniform(m),
        fixed,
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 10_000,
        },
        cx,
    )
    .unwrap()
}
fn brick<'a>(nodes: &'a [[f64; 3]], formulation: Formulation) -> SolidMeshRef<'a> {
    SolidMeshRef::Hex8 {
        nodes_m: nodes,
        elements: &HEXES,
        formulation,
        stiffness_points_per_dir: 2,
    }
}
fn thermal(m: &TetViscoMaterial, free: [f64; 6], temperature: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        temperature,
        300.0,
        free,
        m.material_state_identity,
        ContentHash([99; 32]),
        0.025,
    )
    .unwrap()
}
fn field(states: &[TetThermalStrainState]) -> QuadraticThermalStrainField<'_> {
    QuadraticThermalStrainField::AtQuadrature {
        states,
        points_per_element: 8,
    }
}
fn brick_samples(
    m: &TetViscoMaterial,
    strain: impl Fn([f64; 3]) -> [f64; 6],
) -> Vec<TetThermalStrainState> {
    let (q, _) = gauss_legendre(2);
    let mut samples = Vec::new();
    for &x in &q {
        for &y in &q {
            for &z in &q {
                let p = [x, y, z].map(|v| 0.5 * (1.0 + v));
                samples.push(thermal(m, strain(p), 300.0 + 100.0 * p[2]));
            }
        }
    }
    samples
}
fn expansion_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, value) in ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES
        .into_iter()
        .zip([1e-5, 0.0, 0.0])
    {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, INVERSE_TEMPERATURE_DIMS),
                value: PropertyValue::Scalar {
                    value,
                    dims: INVERSE_TEMPERATURE_DIMS,
                },
                validity: ValidityDomain::unconstrained().with("T", 250.0, 450.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic directional expansion".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "gradient fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn original_nodal_temperature_sources_drive_compatible_quadratic_bending_and_cooling() {
    with_cx(|cx, _| {
        let p1 = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (tet_nodes, tet_cells) = promote_p1_to_p2(&p1, &[[0, 1, 2, 3]], None, cx).unwrap();
        let physical = |p: [i8; 3]| p.map(|x| 0.5 * (1.0 + f64::from(x)));
        let h20 = HEX20_REFERENCE_NODES.map(physical);
        let h27 = HEX27_REFERENCE_NODES.map(physical);
        let c20 = [core::array::from_fn(|i| i)];
        let c27 = [core::array::from_fn(|i| i)];
        let m = material(true);
        let elastic = m.equilibrium_material().unwrap();
        let card = expansion_card();
        let model = ThermalExpansionModel::OrientedOrthotropic {
            principal_to_mesh: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            principal_frame: FRAME,
        };
        let reference = QueryPoint::new().with("T", 300.0).unwrap();
        let original = model
            .resolve(
                &card,
                &reference,
                &reference,
                MaterialPropertySelection::SingleClaimOnly,
                &elastic,
                FRAME,
                0.025,
            )
            .unwrap();
        let source = QuadraticThermalSource {
            card: &card,
            reference_expansion: original.reference(),
            model,
        };
        for (mesh, nodes) in [
            (
                SolidMeshRef::Tet10 {
                    nodes_m: &tet_nodes,
                    elements: &tet_cells,
                    quadrature: TetQuadrature {
                        stiffness_points_per_dir: 3,
                        mass_points_per_dir: 4,
                    },
                },
                &tet_nodes[..],
            ),
            (
                SolidMeshRef::Hex {
                    nodes_m: &h20,
                    elements: HexConnectivity::Hex20(&c20),
                    quadrature: HexQuadrature { points_per_dir: 3 },
                },
                &h20[..],
            ),
            (
                SolidMeshRef::Hex {
                    nodes_m: &h27,
                    elements: HexConnectivity::Hex27(&c27),
                    quadrature: HexQuadrature { points_per_dir: 3 },
                },
                &h27[..],
            ),
        ] {
            let fixed: Vec<_> = nodes
                .iter()
                .enumerate()
                .filter(|(_, p)| p[0] == 0.0)
                .flat_map(|(n, _)| [3 * n, 3 * n + 1, 3 * n + 2])
                .collect();
            let solid = prepare(mesh, &m, &fixed, cx);
            let mass = solid.equilibrium().total_mass_kg;
            let operator = solid.quasistatic_operator(0.05, cx).unwrap();
            let mut history = solid.state(cx).unwrap();
            let zero = vec![0.0; 3 * nodes.len()];
            let temperatures0 = vec![300.0; nodes.len()];
            for scale in [1.0, 0.4, -0.2] {
                let temperatures1: Vec<_> =
                    nodes.iter().map(|p| 300.0 + 100.0 * scale * p[2]).collect();
                let temperatures = NodalThermalTemperatures {
                    reference_k: &temperatures0,
                    current_k: &temperatures1,
                };
                let sources = QuadraticThermalSources::Uniform(&source);
                let samples = match mesh {
                    SolidMeshRef::Tet10 {
                        nodes_m,
                        elements,
                        quadrature,
                    } => QuadraticTetProblem {
                        nodes_m,
                        elements,
                        quadrature,
                        materials: TetMaterialField::Uniform(&elastic),
                        fixed_dofs: &fixed,
                        budget: TetAssemblyBudget::standard(),
                    }
                    .resolve_nodal_thermal_strains(
                        temperatures,
                        sources,
                        FRAME,
                        0.025,
                        cx,
                    ),
                    SolidMeshRef::Hex {
                        nodes_m,
                        elements,
                        quadrature,
                    } => HexElasticProblem {
                        nodes_m,
                        elements,
                        quadrature,
                        family: match elements {
                            HexConnectivity::Hex20(_) => HexFamily::Serendipity20,
                            HexConnectivity::Hex27(_) => HexFamily::Lagrange27,
                        },
                        materials: TetMaterialField::Uniform(&elastic),
                        fixed_dofs: &fixed,
                        budget: TetAssemblyBudget::standard(),
                    }
                    .resolve_nodal_thermal_strains(
                        temperatures,
                        sources,
                        FRAME,
                        0.025,
                        cx,
                    ),
                    _ => unreachable!(),
                }
                .unwrap();
                let step = history
                    .step_quasistatic_with_thermal_strain(
                        &operator,
                        &zero,
                        &zero,
                        &samples,
                        config(),
                        cx,
                    )
                    .unwrap();
                for (u, &[x, _, z]) in step.displacement_m.chunks_exact(3).zip(nodes) {
                    let b = 0.001 * scale;
                    for (&actual, expected) in u.iter().zip([b * x * z, 0.0, -0.5 * b * x * x]) {
                        near(actual, expected, 3e-11);
                    }
                }
                assert!(step.constitutive.energy.stored_j < 1e-18);
                assert!(step.constitutive.step_dissipated_j < 1e-18);
                assert_eq!(history.thermal_strains().unwrap(), samples.strains());
                assert!(matches!(
                    history.thermal_field(),
                    Some(QuadraticThermalStrainField::AtQuadrature {
                        points_per_element: 27,
                        ..
                    })
                ));
                assert_eq!(solid.equilibrium().total_mass_kg, mass);
            }
        }
    });
}

#[test]
fn zero_mean_thermal_gradient_keeps_relaxation_energy_and_drives_native_release() {
    with_cx(|cx, _| {
        let m = material(true);
        let (c, branch, eta) = tensors();
        let (a, h, tau): (f64, f64, f64) = (0.003, 0.08, 0.2);
        let phi = tau / h * (1.0 - (-h / tau).exp());
        let fixed: Vec<_> = HEX_NODES
            .iter()
            .enumerate()
            .filter(|(_, p)| p[0] == 0.0)
            .flat_map(|(n, _)| [3 * n, 3 * n + 1, 3 * n + 2])
            .collect();
        for form in [Formulation::Standard, Formulation::BBar] {
            let solid = prepare(brick(&HEX_NODES, form), &m, &fixed, cx);
            let samples = brick_samples(&m, |p| {
                let d = a * (p[2] - 0.5);
                [d, 0.0, 0.0, 0.6 * d, 0.0, 0.0]
            });
            let direction = if form == Formulation::Standard {
                [1.0, 0.0, 0.0, 0.6, 0.0, 0.0]
            } else {
                [2.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0, 0.6, 0.0, 0.0]
            };
            let eq = quadratic(direction, &c, direction);
            let memory = quadratic(direction, &branch, direction);
            let viscosity = quadratic(direction, &eta, direction);
            // Integral (z-1/2)^2 dV=1/12 on the unit cube.
            let stored = a * a * (eq + memory * phi * phi) / 24.0;
            let work = a * a * (0.5 * eq + memory * tau / h * (1.0 - phi) + viscosity / h) / 12.0;
            let mut history = solid.state(cx).unwrap();
            let heated = history
                .step_internal_force_with_thermal_strain(&[0.0; 24], h, field(&samples), cx)
                .unwrap();
            relative(heated.energy.stored_j, stored);
            relative(heated.step_eigenstrain_work_j, work);
            relative(heated.step_dissipated_j, work - stored);
            assert_eq!(heated.step_mechanical_work_j, 0.0);
            let hold: f64 = 0.12;
            let held = history.step_internal_force(&[0.0; 24], hold, cx).unwrap();
            let relaxed = a * a * (eq + memory * (phi * (-hold / tau).exp()).powi(2)) / 24.0;
            relative(held.energy.stored_j, relaxed);
            relative(held.step_dissipated_j, stored - relaxed);
            assert_eq!(held.step_eigenstrain_work_j, 0.0);
            assert_eq!(history.thermal_strains().unwrap(), samples);
            let previous_heat = history.energy().dissipated_j;
            let operator = solid
                .dynamic_operator(
                    0.005,
                    ViscoDynamicConfig {
                        maximum_steps: 30,
                        maximum_free_dofs: 30,
                        solve: config(),
                        absolute_energy_tolerance_j: 1e-13,
                        relative_energy_tolerance: 2e-9,
                    },
                    cx,
                )
                .unwrap();
            let mut motion = operator.start(history, &[0.0; 24], cx).unwrap();
            for _ in 0..20 {
                let step = motion.advance(&[0.0; 24], cx).unwrap();
                assert_eq!(step.constitutive.step_eigenstrain_work_j, 0.0);
                near(
                    step.kinetic_energy_j
                        + step.constitutive.energy.stored_j
                        + step.constitutive.energy.dissipated_j
                        - previous_heat,
                    relaxed,
                    3e-12,
                );
            }
            assert!(motion.kinetic_energy_j() > 1e-10);
            assert_eq!(
                motion.material_history().thermal_strains().unwrap(),
                samples
            );
        }
    });
}

#[test]
fn warped_bbar_thermal_projection_uses_physical_volume_and_full_tensor_conjugacy() {
    with_cx(|cx, _| {
        let m = material(false);
        let (c, _, _) = tensors();
        let k = 0.6;
        let nodes = HEX_NODES.map(|[x, y, z]| [x, y, z * (1.0 + k * x)]);
        // x=s, y=t, z=r(1+k*s). detJ=1+k*x in unit reference coordinates.
        let moments = [1.0 + k / 2.0, 0.5 + k / 3.0, 1.0 / 3.0 + k / 4.0];
        let mean_x = moments[1] / moments[0];
        let (a, b, g) = (0.002, 0.001, 0.0006);
        let samples = brick_samples(&m, |p| {
            [p[0] * (a + b), p[0] * (a - b), p[0] * a, p[0] * g, 0.0, 0.0]
        });
        let c0 = [a * mean_x, a * mean_x, a * mean_x, 0.0, 0.0, 0.0];
        let c1 = [b, -b, 0.0, g, 0.0, 0.0];
        let expected = 0.5
            * (moments[0] * quadratic(c0, &c, c0)
                + 2.0 * moments[1] * quadratic(c0, &c, c1)
                + moments[2] * quadratic(c1, &c, c1));
        let solid = prepare(brick(&nodes, Formulation::BBar), &m, &[], cx);
        let mut history = solid.state(cx).unwrap();
        let response = history
            .step_internal_force_with_thermal_strain(&[0.0; 24], 0.1, field(&samples), cx)
            .unwrap();
        relative(response.energy.stored_j, expected);
        relative(response.step_eigenstrain_work_j, expected);
        assert_eq!(response.step_dissipated_j, 0.0);
        relative(
            solid.equilibrium().total_mass_kg,
            m.density_kg_m3 * moments[0],
        );
        // Virtual v=(0,x*y,0) has projected normals (mean_x-x,mean_x+2x,mean_x-x)/3
        // and Mandel xy=y/sqrt(2). Its mean over y supplies these two polynomials.
        let d0 = [
            mean_x / 3.0,
            mean_x / 3.0,
            mean_x / 3.0,
            0.5 * core::f64::consts::FRAC_1_SQRT_2,
            0.0,
            0.0,
        ];
        let d1 = [-1.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0, 0.0, 0.0, 0.0];
        let virtual_work = -(moments[0] * quadratic(d0, &c, c0)
            + moments[1] * (quadratic(d0, &c, c1) + quadratic(d1, &c, c0))
            + moments[2] * quadratic(d1, &c, c1));
        let nodal: f64 = nodes
            .iter()
            .enumerate()
            .map(|(n, p)| response.full_force_n[3 * n + 1] * p[0] * p[1])
            .sum();
        relative(nodal, virtual_work);
        let unweighted = [a * 0.5, a * 0.5, a * 0.5, 0.0, 0.0, 0.0];
        let wrong = 0.5
            * (moments[0] * quadratic(unweighted, &c, unweighted)
                + 2.0 * moments[1] * quadratic(unweighted, &c, c1)
                + moments[2] * quadratic(c1, &c, c1));
        assert!((wrong - expected).abs() > 1e-6);
        assert_eq!(history.thermal_strains().unwrap(), samples);

        let rotation = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let rotated = TetViscoMaterial {
            law: m.law.rotated(rotation).unwrap(),
            ..m.clone()
        };
        let rotated_nodes = nodes.map(|[x, y, z]| [-y, x, z]);
        let map = fs_material::tensor::mandel_rotation(rotation);
        let rotated_samples: Vec<_> = samples
            .iter()
            .map(|s| {
                let old = s.free_strain_mandel();
                thermal(
                    &rotated,
                    core::array::from_fn(|i| (0..6).map(|j| map[i][j] * old[j]).sum()),
                    s.temperature_k(),
                )
            })
            .collect();
        let rotated_solid = prepare(brick(&rotated_nodes, Formulation::BBar), &rotated, &[], cx);
        let rotated_response = rotated_solid
            .state(cx)
            .unwrap()
            .step_internal_force_with_thermal_strain(&[0.0; 24], 0.1, field(&rotated_samples), cx)
            .unwrap();
        relative(rotated_response.energy.stored_j, expected);
        for (old, new) in response
            .full_force_n
            .chunks_exact(3)
            .zip(rotated_response.full_force_n.chunks_exact(3))
        {
            for (&actual, exact) in new.iter().zip([-old[1], old[0], old[2]]) {
                near(actual, exact, 2e-13);
            }
        }
    });
}

#[test]
fn original_constant_field_is_exact_and_point_source_domain_failures_are_transactional() {
    with_cx(|cx, gate| {
        let mut m = material(true);
        let (c, branch, _) = tensors();
        m.law = TensorProny::new_mandel(c, vec![PronyBranch::new(branch, 0.2).unwrap()], 0.012)
            .unwrap();
        let solid = prepare(brick(&HEX_NODES, Formulation::BBar), &m, &[], cx);
        let constant = thermal(&m, [0.001, 0.0002, -0.0001, 0.0003, 0.0, 0.0], 330.0);
        let uniform = vec![constant; 8];
        let displacement: Vec<_> = HEX_NODES
            .iter()
            .flat_map(|[x, y, _]| [0.0005 * x * y, 0.0, 0.0])
            .collect();
        let mut old_path = solid.state(cx).unwrap();
        let mut point_path = solid.state(cx).unwrap();
        let a = old_path
            .step_internal_force_with_thermal_strain(
                &displacement,
                0.1,
                TetThermalStrainField::Uniform(&constant),
                cx,
            )
            .unwrap();
        let b = point_path
            .step_internal_force_with_thermal_strain(&displacement, 0.1, field(&uniform), cx)
            .unwrap();
        assert_eq!(a.full_force_n, b.full_force_n);
        assert_eq!(a.full_mean_force_n, b.full_mean_force_n);
        assert_eq!(a.energy, b.energy);
        assert_eq!(a.step_eigenstrain_work_j, b.step_eigenstrain_work_j);
        let saved = point_path.clone();
        let before = point_path.energy();
        let work = point_path.eigenstrain_work_j();
        let mut wrong_id = uniform.clone();
        wrong_id[7] = TetThermalStrainState::try_new(
            330.0,
            300.0,
            constant.free_strain_mandel(),
            [255; 32],
            ContentHash([99; 32]),
            0.025,
        )
        .unwrap();
        assert!(matches!(
            point_path.step_internal_force_with_thermal_strain(
                &displacement,
                0.1,
                field(&wrong_id),
                cx
            ),
            Err(ViscoElementError::Elastic(
                TetElasticError::ThermalElasticStateMismatch { element: 0 }
            ))
        ));
        let mut wrong_reference = uniform.clone();
        wrong_reference[7] = TetThermalStrainState::try_new(
            330.0,
            301.0,
            constant.free_strain_mandel(),
            m.material_state_identity,
            ContentHash([99; 32]),
            0.025,
        )
        .unwrap();
        assert!(
            point_path
                .step_internal_force_with_thermal_strain(
                    &displacement,
                    0.1,
                    field(&wrong_reference),
                    cx
                )
                .is_err()
        );
        for bad in [
            QuadraticThermalStrainField::AtQuadrature {
                states: &uniform,
                points_per_element: 7,
            },
            field(&uniform[..7]),
        ] {
            assert!(
                point_path
                    .step_internal_force_with_thermal_strain(&displacement, 0.1, bad, cx)
                    .is_err()
            );
        }
        let materials = vec![m.equilibrium_material().unwrap(); 8];
        assert!(
            point_path
                .step_internal_force_with_thermal_strain(
                    &displacement,
                    0.1,
                    QuadraticThermalStrainField::ThermomechanicalAtQuadrature {
                        states: &uniform,
                        materials: &materials,
                        points_per_element: 8,
                    },
                    cx
                )
                .is_err()
        );
        let hidden = brick_samples(&m, |p| {
            let s = if p[2] < 0.5 { -0.008 } else { 0.008 };
            [s, s, s, 0.0, 0.0, 0.0]
        });
        // Its projected thermal volume strain is zero, but the original raw
        // mechanical norm sqrt(3)*0.008 exceeds the source's 0.012 limit.
        assert!(matches!(
            point_path.step_internal_force_with_thermal_strain(&[0.0; 24], 0.1, field(&hidden), cx),
            Err(ViscoElementError::Material(_))
        ));
        assert_eq!(point_path.energy(), before);
        assert_eq!(point_path.eigenstrain_work_j(), work);
        assert_eq!(point_path.displacement_m(), displacement);
        assert_eq!(point_path.thermal_strains().unwrap(), uniform);
        let expected = saved
            .clone()
            .step_internal_force(&displacement, 0.1, cx)
            .unwrap();
        let actual = point_path
            .step_internal_force(&displacement, 0.1, cx)
            .unwrap();
        assert_eq!(actual.full_force_n, expected.full_force_n);
        assert_eq!(actual.energy, expected.energy);
        let before_cancel = point_path.energy();
        gate.request();
        assert!(matches!(
            point_path.step_internal_force_with_thermal_strain(
                &displacement,
                0.1,
                field(&uniform),
                cx
            ),
            Err(ViscoElementError::Elastic(TetElasticError::Cancelled))
        ));
        assert_eq!(point_path.energy(), before_cancel);
        assert_eq!(point_path.thermal_strains().unwrap(), uniform);
    });
}
