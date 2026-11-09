//! Consistent graded reference inertia must also own body loads and reactions.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear::Formulation;
use fs_solid::linear3::hex::{
    HEX20_REFERENCE_NODES, HEX27_REFERENCE_NODES, HexConnectivity, HexElasticProblem, HexFamily,
    HexQuadrature,
};
use fs_solid::linear3::hex8::Hex8ElasticProblem;
use fs_solid::linear3::tet10::{QuadraticTetProblem, TetQuadrature};
use fs_solid::linear3::{
    QuadraticThermalResponse, ReferenceAccelerationField, TetAssemblyBudget, TetElasticAssembly,
    TetElasticError, TetElasticMaterial, TetMaterialField, TetStaticSolveConfig,
    TetThermalStrainField, TetThermalStrainState,
};

#[derive(Debug, Clone, Copy)]
enum Kind {
    Tet10,
    Hex8(Formulation),
    Hex20,
    Hex27,
}
const KINDS: [Kind; 5] = [
    Kind::Tet10,
    Kind::Hex8(Formulation::Standard),
    Kind::Hex8(Formulation::BBar),
    Kind::Hex20,
    Kind::Hex27,
];
const TET_NODES: [[f64; 3]; 10] = [
    [0., 0., 0.],
    [1., 0., 0.],
    [0., 1., 0.],
    [0., 0., 1.],
    [0.5, 0., 0.],
    [0.5, 0.5, 0.],
    [0., 0.5, 0.],
    [0., 0., 0.5],
    [0.5, 0., 0.5],
    [0., 0.5, 0.5],
];
const TET: [[usize; 10]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]];
const HEX8: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];
const HEX20: [[usize; 20]; 1] = [[
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
]];
const HEX27: [[usize; 27]; 1] = [[
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26,
]];

struct Profiles {
    tet: [[f64; 4]; 1],
    hex: [[f64; 8]; 1],
}
impl Profiles {
    fn graded() -> Self {
        Self {
            tet: [[2., 3., 4., 5.]],
            hex: [core::array::from_fn(|i| {
                let [r, s, t] = HEX27_REFERENCE_NODES[i].map(f64::from);
                2. + 0.5 * r + 0.25 * s + 0.125 * t + 0.0625 * r * s * t
            })],
        }
    }
    fn uniform(value: f64) -> Self {
        Self {
            tet: [[value; 4]],
            hex: [[value; 8]],
        }
    }
}

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 67,
                    kernel_id: 67,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}

fn material(kind: Kind) -> TetElasticMaterial {
    TetElasticMaterial::try_new(
        if matches!(kind, Kind::Tet10) { 3.5 } else { 2. },
        1e6,
        0.25,
        [67; 32],
    )
    .unwrap()
}
fn tet<'a>(
    nodes: &'a [[f64; 3]],
    law: &'a TetElasticMaterial,
    fixed: &'a [usize],
    points: usize,
) -> QuadraticTetProblem<'a> {
    QuadraticTetProblem {
        nodes_m: nodes,
        elements: &TET,
        materials: TetMaterialField::Uniform(law),
        fixed_dofs: fixed,
        budget: TetAssemblyBudget::standard(),
        quadrature: TetQuadrature {
            stiffness_points_per_dir: 3,
            mass_points_per_dir: points,
        },
    }
}
fn hex8<'a>(
    nodes: &'a [[f64; 3]],
    law: &'a TetElasticMaterial,
    fixed: &'a [usize],
    form: Formulation,
) -> Hex8ElasticProblem<'a> {
    Hex8ElasticProblem {
        nodes_m: nodes,
        hexahedra: &HEX8,
        materials: TetMaterialField::Uniform(law),
        fixed_dofs: fixed,
        budget: TetAssemblyBudget::standard(),
        formulation: form,
        stiffness_points_per_dir: 2,
    }
}
fn hex<'a>(
    kind: Kind,
    nodes: &'a [[f64; 3]],
    law: &'a TetElasticMaterial,
    fixed: &'a [usize],
) -> HexElasticProblem<'a> {
    let (family, elements) = match kind {
        Kind::Hex20 => (HexFamily::Serendipity20, HexConnectivity::Hex20(&HEX20)),
        Kind::Hex27 => (HexFamily::Lagrange27, HexConnectivity::Hex27(&HEX27)),
        _ => unreachable!(),
    };
    HexElasticProblem {
        nodes_m: nodes,
        family,
        elements,
        materials: TetMaterialField::Uniform(law),
        fixed_dofs: fixed,
        budget: TetAssemblyBudget::standard(),
        quadrature: HexQuadrature { points_per_dir: 3 },
    }
}
#[allow(clippy::too_many_arguments)]
fn load(
    kind: Kind,
    nodes: &[[f64; 3]],
    law: &TetElasticMaterial,
    fixed: &[usize],
    profile: &Profiles,
    acceleration: ReferenceAccelerationField<'_>,
    points: usize,
    cx: &Cx<'_>,
) -> Result<Vec<(usize, f64)>, TetElasticError> {
    match kind {
        Kind::Tet10 => tet(nodes, law, fixed, points)
            .assemble_acceleration_load_with_reference_density(&profile.tet, acceleration, cx),
        Kind::Hex8(form) => hex8(nodes, law, fixed, form)
            .assemble_acceleration_load_with_reference_density(&profile.hex, acceleration, cx),
        _ => hex(kind, nodes, law, fixed).assemble_acceleration_load_with_reference_density(
            &profile.hex,
            HexQuadrature {
                points_per_dir: points,
            },
            acceleration,
            cx,
        ),
    }
}
fn mass(
    kind: Kind,
    nodes: &[[f64; 3]],
    law: &TetElasticMaterial,
    profile: &Profiles,
    cx: &Cx<'_>,
) -> TetElasticAssembly {
    match kind {
        Kind::Tet10 => tet(nodes, law, &[], 6).assemble_with_reference_density(&profile.tet, cx),
        Kind::Hex8(form) => {
            hex8(nodes, law, &[], form).assemble_with_reference_density(&profile.hex, cx)
        }
        _ => hex(kind, nodes, law, &[]).assemble_with_reference_density(
            &profile.hex,
            HexQuadrature { points_per_dir: 6 },
            cx,
        ),
    }
    .unwrap()
}
fn response(
    kind: Kind,
    nodes: &[[f64; 3]],
    law: &TetElasticMaterial,
    fixed: &[usize],
    forces: &[(usize, f64)],
    cx: &Cx<'_>,
) -> QuadraticThermalResponse {
    let state = TetThermalStrainState::try_new(
        300.,
        300.,
        [0.; 6],
        law.material_state_identity,
        ContentHash([68; 32]),
        0.02,
    )
    .unwrap();
    let field = TetThermalStrainField::Uniform(&state);
    let config = TetStaticSolveConfig::default();
    let frame = ContentHash([69; 32]);
    match kind {
        Kind::Tet10 => tet(nodes, law, fixed, 6)
            .solve_thermal_response_with_nodal_forces(field, forces, config, frame, 0.02, cx),
        Kind::Hex8(form) => hex8(nodes, law, fixed, form)
            .solve_thermal_response_with_nodal_forces(field, forces, config, frame, 0.02, cx),
        _ => hex(kind, nodes, law, fixed)
            .solve_thermal_response_with_nodal_forces(field, forces, config, frame, 0.02, cx),
    }
    .unwrap()
}

type Poly = Vec<([usize; 3], f64)>;
fn multiply(a: &Poly, b: &Poly) -> Poly {
    a.iter()
        .flat_map(|(p, c)| {
            b.iter()
                .map(move |(q, d)| (core::array::from_fn(|i| p[i] + q[i]), c * d))
        })
        .collect()
}
fn integrate(kind: Kind, p: &Poly) -> f64 {
    let factorial = |n: usize| (1..=n).map(|n| n as f64).product::<f64>();
    p.iter()
        .map(|(power, c)| {
            c * if matches!(kind, Kind::Tet10) {
                power.iter().map(|&n| factorial(n)).product::<f64>()
                    / factorial(3 + power.iter().sum::<usize>())
            } else {
                power
                    .iter()
                    .map(|n| if n % 2 == 0 { 2. / (n + 1) as f64 } else { 0. })
                    .product::<f64>()
            }
        })
        .sum()
}
fn evaluate(p: &Poly, x: [f64; 3]) -> f64 {
    p.iter()
        .map(|(power, c)| c * (0..3).map(|i| x[i].powi(power[i] as i32)).product::<f64>())
        .sum()
}
fn geometry(kind: Kind, curved: bool) -> (Vec<[f64; 3]>, [Poly; 3], Poly) {
    let mut coordinates: [Poly; 3] = core::array::from_fn(|axis| {
        let mut p = [0; 3];
        p[axis] = 1;
        vec![(p, 1.)]
    });
    let mut jacobian = vec![([0; 3], 1.)];
    if curved {
        match kind {
            Kind::Tet10 => {
                for (axis, a) in [0.25, 0.125, 0.375].into_iter().enumerate() {
                    let mut p = [0; 3];
                    p[axis] = 2;
                    coordinates[axis].push((p, a));
                    p[axis] = 1;
                    jacobian = multiply(&jacobian, &vec![([0; 3], 1.), (p, 2. * a)]);
                }
            }
            Kind::Hex8(_) => {
                coordinates[0].push(([0, 1, 1], 0.125));
                coordinates[1].push(([1, 0, 1], 0.25));
                coordinates[2].push(([1, 1, 0], 0.375));
                jacobian = vec![
                    ([0; 3], 1.),
                    ([2, 0, 0], -0.25 * 0.375),
                    ([0, 2, 0], -0.125 * 0.375),
                    ([0, 0, 2], -0.125 * 0.25),
                    ([1, 1, 1], 2. * 0.125 * 0.25 * 0.375),
                ];
            }
            _ => {
                coordinates[0].push(([2, 0, 0], 0.125));
                coordinates[1].push(([2, 1, 0], 0.25));
                coordinates[2].push(([2, 0, 1], 0.375));
                for factor in [
                    vec![([0; 3], 1.), ([1, 0, 0], 0.25)],
                    vec![([0; 3], 1.), ([2, 0, 0], 0.25)],
                    vec![([0; 3], 1.), ([2, 0, 0], 0.375)],
                ] {
                    jacobian = multiply(&jacobian, &factor);
                }
            }
        }
    }
    let reference: Vec<_> = match kind {
        Kind::Tet10 => TET_NODES.to_vec(),
        Kind::Hex8(_) => HEX27_REFERENCE_NODES[..8]
            .iter()
            .map(|p| p.map(f64::from))
            .collect(),
        Kind::Hex20 => HEX20_REFERENCE_NODES
            .iter()
            .map(|p| p.map(f64::from))
            .collect(),
        Kind::Hex27 => HEX27_REFERENCE_NODES
            .iter()
            .map(|p| p.map(f64::from))
            .collect(),
    };
    let density = if matches!(kind, Kind::Tet10) {
        vec![
            ([0; 3], 2.),
            ([1, 0, 0], 1.),
            ([0, 1, 0], 2.),
            ([0, 0, 1], 3.),
        ]
    } else {
        vec![
            ([0; 3], 2.),
            ([1, 0, 0], 0.5),
            ([0, 1, 0], 0.25),
            ([0, 0, 1], 0.125),
            ([1, 1, 1], 0.0625),
        ]
    };
    let nodes = reference
        .into_iter()
        .map(|p| core::array::from_fn(|axis| evaluate(&coordinates[axis], p)))
        .collect();
    (nodes, coordinates, multiply(&density, &jacobian))
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 4e-11 * b.abs().max(1.),
        "{a:.17e} != {b:.17e}"
    );
}
fn resultant(nodes: &[[f64; 3]], force: &[(usize, f64)]) -> ([f64; 3], [f64; 3]) {
    assert_eq!(force.len(), 3 * nodes.len());
    let mut total = [0.; 3];
    let mut torque = [0.; 3];
    for (i, (point, rows)) in nodes.iter().zip(force.chunks_exact(3)).enumerate() {
        for (axis, &(dof, value)) in rows.iter().enumerate() {
            assert_eq!(dof, 3 * i + axis);
            total[axis] += value;
        }
        let moment = cross(*point, [rows[0].1, rows[1].1, rows[2].1]);
        for axis in 0..3 {
            torque[axis] += moment[axis];
        }
    }
    (total, torque)
}

#[test]
fn uniform_body_load_preserves_exact_mass_centroid_torque_and_rigid_virtual_work() {
    with_cx(|cx, _| {
        for kind in KINDS {
            for curved in [false, true] {
                let (nodes, coordinates, measure) = geometry(kind, curved);
                let law = material(kind);
                let profile = Profiles::graded();
                let g = [1.3, -2.4, 0.7];
                let fixed: Vec<_> = (0..3 * nodes.len()).collect();
                let force = load(
                    kind,
                    &nodes,
                    &law,
                    &fixed,
                    &profile,
                    ReferenceAccelerationField::Uniform(g),
                    6,
                    cx,
                )
                .unwrap();
                let exact_mass = integrate(kind, &measure);
                let first = core::array::from_fn(|axis| {
                    integrate(kind, &multiply(&coordinates[axis], &measure))
                });
                let (total, torque) = resultant(&nodes, &force);
                for axis in 0..3 {
                    close(total[axis], exact_mass * g[axis]);
                    close(torque[axis], cross(first, g)[axis]);
                }
                let translation = [0.2, -0.1, 0.3];
                let rotation = [0.1, 0.2, -0.15];
                let work: f64 = nodes
                    .iter()
                    .zip(force.chunks_exact(3))
                    .map(|(p, f)| {
                        let r = cross(rotation, *p);
                        (0..3)
                            .map(|axis| (translation[axis] + r[axis]) * f[axis].1)
                            .sum::<f64>()
                    })
                    .sum();
                close(
                    work,
                    dot(translation, g) * exact_mass + dot(rotation, cross(first, g)),
                );
                let nodal = vec![g; nodes.len()];
                assert_eq!(
                    force,
                    load(
                        kind,
                        &nodes,
                        &law,
                        &fixed,
                        &profile,
                        ReferenceAccelerationField::PerNode(&nodal),
                        6,
                        cx
                    )
                    .unwrap()
                );
                if matches!(kind, Kind::Tet10 | Kind::Hex20) && !curved {
                    assert!(force[0].1 < 0., "signed quadratic corner load was lumped");
                }
            }
        }
    });
}

#[test]
fn full_consistent_mass_action_keeps_supported_acceleration_columns_and_reciprocity() {
    with_cx(|cx, _| {
        for kind in KINDS {
            let (nodes, _, _) = geometry(kind, true);
            let law = material(kind);
            let profile = Profiles::graded();
            let native = mass(kind, &nodes, &law, &profile, cx);
            let fixed: Vec<_> = (1..3 * nodes.len()).collect();
            let acceleration: Vec<_> = (0..nodes.len())
                .map(|i| {
                    [
                        0.02 * i as f64 - 0.3,
                        0.1 * (i % 3) as f64 - 0.2,
                        0.05 * (i % 5) as f64,
                    ]
                })
                .collect();
            let force = load(
                kind,
                &nodes,
                &law,
                &fixed,
                &profile,
                ReferenceAccelerationField::PerNode(&acceleration),
                6,
                cx,
            )
            .unwrap();
            let mut expected = vec![0.; 3 * nodes.len()];
            native.mass.spmv(
                &acceleration.iter().flatten().copied().collect::<Vec<_>>(),
                &mut expected,
            );
            for ((_, actual), expected) in force.iter().zip(expected) {
                close(*actual, expected);
            }
            let mut support = vec![[0.; 3]; nodes.len()];
            support[1][0] = 1.;
            let action = load(
                kind,
                &nodes,
                &law,
                &fixed,
                &profile,
                ReferenceAccelerationField::PerNode(&support),
                6,
                cx,
            )
            .unwrap();
            close(action[0].1, native.mass.get(0, 3));
            assert!(action[0].1.abs() > 1e-8, "supported column was discarded");
            let all_fixed: Vec<_> = (0..3 * nodes.len()).collect();
            assert_eq!(
                action,
                load(
                    kind,
                    &nodes,
                    &law,
                    &all_fixed,
                    &profile,
                    ReferenceAccelerationField::PerNode(&support),
                    6,
                    cx
                )
                .unwrap()
            );
            let left: f64 = acceleration
                .iter()
                .flatten()
                .zip(&action)
                .map(|(a, (_, f))| a * f)
                .sum();
            close(left, force[3].1);
            let resampled = TetElasticMaterial::try_new(0.5, 1e6, 0.25, [67; 32]).unwrap();
            assert_eq!(
                force,
                load(
                    kind,
                    &nodes,
                    &resampled,
                    &fixed,
                    &profile,
                    ReferenceAccelerationField::PerNode(&acceleration),
                    6,
                    cx
                )
                .unwrap()
            );
        }
    });
}

#[test]
fn self_weight_changes_equilibrium_and_support_reactions_at_equal_total_mass() {
    with_cx(|cx, _| {
        for kind in KINDS {
            let (nodes, coordinates, measure) = geometry(kind, false);
            let law = material(kind);
            let fixed: Vec<_> = (1..3 * nodes.len()).collect();
            let all_fixed: Vec<_> = (0..3 * nodes.len()).collect();
            let g = [9.81, 0., 0.];
            let acceleration = ReferenceAccelerationField::Uniform(g);
            let graded = load(
                kind,
                &nodes,
                &law,
                &fixed,
                &Profiles::graded(),
                acceleration,
                6,
                cx,
            )
            .unwrap();
            let uniform = load(
                kind,
                &nodes,
                &law,
                &fixed,
                &Profiles::uniform(law.density_kg_m3),
                acceleration,
                6,
                cx,
            )
            .unwrap();
            let (total, torque) = resultant(&nodes, &graded);
            let (uniform_total, uniform_torque) = resultant(&nodes, &uniform);
            for axis in 0..3 {
                close(total[axis], uniform_total[axis]);
            }
            assert!((torque[1] - uniform_torque[1]).abs() > 1e-3);
            let a = response(kind, &nodes, &law, &fixed, &graded, cx);
            let b = response(kind, &nodes, &law, &fixed, &uniform, cx);
            assert!(a.iterations() > 0);
            close(
                a.displacement_m()[0][0] / b.displacement_m()[0][0],
                graded[0].1 / uniform[0].1,
            );
            assert!((a.displacement_m()[0][0] / b.displacement_m()[0][0] - 1.).abs() > 0.01);
            let reactions: Vec<_> = a
                .nodal_force_balance_n()
                .iter()
                .flatten()
                .copied()
                .enumerate()
                .collect();
            let (support, moment) = resultant(&nodes, &reactions);
            for axis in 0..3 {
                close(support[axis], -total[axis]);
                close(moment[axis], -torque[axis]);
            }
            let clamped = response(kind, &nodes, &law, &all_fixed, &graded, cx);
            assert!(clamped.displacement_m().iter().flatten().all(|v| *v == 0.));
            for (reaction, (_, force)) in clamped
                .nodal_force_balance_n()
                .iter()
                .flatten()
                .zip(&graded)
            {
                close(*reaction, -force);
            }
            let first = core::array::from_fn(|axis| {
                integrate(kind, &multiply(&coordinates[axis], &measure))
            });
            for axis in 0..3 {
                close(torque[axis], cross(first, g)[axis]);
            }
        }
    });
}

#[test]
fn curved_profiles_retain_native_quadrature_and_acceleration_refusals() {
    for kind in KINDS {
        with_cx(|cx, gate| {
            let (nodes, _, _) = geometry(kind, true);
            let law = material(kind);
            let profile = Profiles::graded();
            let fixed: Vec<_> = (0..3 * nodes.len()).collect();
            let call = |profile: &Profiles, a: ReferenceAccelerationField<'_>, points| {
                load(kind, &nodes, &law, &fixed, profile, a, points, cx)
            };
            if !matches!(kind, Kind::Hex8(_)) {
                assert!(matches!(
                    call(&profile, ReferenceAccelerationField::Uniform([0.; 3]), 5),
                    Err(TetElasticError::QuadratureInadmissible { .. })
                ));
            }
            let a = call(
                &profile,
                ReferenceAccelerationField::Uniform([1., 2., 3.]),
                6,
            )
            .unwrap();
            let b = call(
                &profile,
                ReferenceAccelerationField::Uniform([1., 2., 3.]),
                7,
            )
            .unwrap();
            for ((_, a), (_, b)) in a.iter().zip(b) {
                close(*a, b);
            }
            assert!(matches!(
                call(&profile, ReferenceAccelerationField::PerNode(&[]), 6),
                Err(TetElasticError::InvalidBodyAcceleration { .. })
            ));
            for value in [f64::NAN, f64::INFINITY] {
                assert!(matches!(
                    call(
                        &profile,
                        ReferenceAccelerationField::Uniform([value, 0., 0.]),
                        6
                    ),
                    Err(TetElasticError::InvalidBodyAcceleration { .. })
                ));
                let mut nodal = vec![[0.; 3]; nodes.len()];
                nodal[1][1] = value;
                assert!(matches!(
                    call(&profile, ReferenceAccelerationField::PerNode(&nodal), 6),
                    Err(TetElasticError::InvalidBodyAcceleration { .. })
                ));
            }
            for density in [0., -1., f64::NAN, f64::INFINITY] {
                assert!(matches!(
                    call(
                        &Profiles::uniform(density),
                        ReferenceAccelerationField::Uniform([0.; 3]),
                        6
                    ),
                    Err(TetElasticError::InvalidMaterial { .. })
                ));
            }
            assert!(
                call(
                    &Profiles::uniform(1e300),
                    ReferenceAccelerationField::Uniform([1e300; 3]),
                    6
                )
                .is_err()
            );
            gate.request();
            assert!(matches!(
                call(&profile, ReferenceAccelerationField::Uniform([1.; 3]), 6),
                Err(TetElasticError::Cancelled)
            ));
        });
    }
}
