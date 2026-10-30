//! G1/G3/G4: same source coefficients drive real inertia, motion and material heat.
//! Manufactured constitutive cards; no measured-material validation is claimed.

use fs_couple::tensor_material::{
    InitialMemory, TensorMaterialBudget,
    membrane::{
        TensorMembrane,
        source::{MembraneMass, MembraneMaterialSource},
    },
    solid::TensorSolidBody,
};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    state_point::MaterialPropertySelection,
    visco_tensor::{
        TensorProjectionBudget,
        source::{
            ResolvedTensorMaterial, TensorBranchProperties, TensorElasticity, TensorPronyProperties,
        },
    },
};
use fs_plate::shell::ShellMesh;
use fs_qty::{Density, Dims, Pressure, Time};
use fs_solid::linear3::{
    TetAssemblyBudget,
    visco::{
        SolidMeshRef, SolidProjectionBudget, TetViscoField, TetViscoMaterial, ViscoElements,
        ViscoStateBudget,
    },
};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TURN: [[f64; 3]; 3] = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 719,
                    kernel_id: 719,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}
fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained()
            .with("T", 280.0, 330.0)
            .with("omega", 0.0, if name == "density" { 300.0 } else { 1000.0 }),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured same-card tensor specimen".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn fixture(young: f64, density: f64, fraction: f64) -> (Vec<PropertyClaim>, TensorPronyProperties) {
    let mut claims = vec![
        scalar("density", Density::DIMS, density),
        scalar("young_modulus", Pressure::DIMS, young),
        scalar("poisson_ratio", Dims::NONE, 0.0),
        scalar("strain_limit", Dims::NONE, 0.1),
    ];
    let mut branches = Vec::new();
    if fraction > 0.0 {
        let branch = TensorBranchProperties {
            stiffness_upper: core::array::from_fn(|i| format!("branch_c{i:02}")),
            relaxation_time: "branch_tau".into(),
        };
        claims.push(scalar("branch_tau", Time::DIMS, 0.02));
        let mut entry = 0;
        for row in 0..6 {
            for column in row..6 {
                let coefficient = if row == column {
                    [1.0, 0.1, 0.7, 0.3, 0.8, 0.2][row]
                } else if row == 0 && column == 1 {
                    0.05
                } else {
                    0.0
                };
                claims.push(scalar(
                    &branch.stiffness_upper[entry],
                    Pressure::DIMS,
                    fraction * young * coefficient,
                ));
                entry += 1;
            }
        }
        branches.push(branch);
    }
    (
        claims,
        TensorPronyProperties {
            elasticity: TensorElasticity::Isotropic,
            strain_limit: "strain_limit".into(),
            branches,
        },
    )
}
fn card(claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "same material name".into(),
            phase: "solid".into(),
            process: "manufactured tensor test".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap()
}
fn query() -> QueryPoint {
    QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with("omega", 10.0)
        .unwrap()
}
fn material(young: f64, density: f64, fraction: f64) -> ResolvedTensorMaterial {
    let (claims, properties) = fixture(young, density, fraction);
    ResolvedTensorMaterial::resolve(
        &card(claims),
        &query(),
        properties,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn dynamics() -> TensorMaterialBudget {
    TensorMaterialBudget {
        max_state_coordinates: 128,
        max_material_points: 16,
    }
}
fn projection() -> TensorProjectionBudget {
    TensorProjectionBudget {
        max_coordinates: 16,
        max_scalar_terms: 1_000_000,
    }
}

#[test]
fn tetrahedral_inertia_restoring_force_motion_and_heat_follow_card_substitutions() {
    with_cx(|cx, gate| {
        let mut motions = Vec::new();
        for (young, density, fraction) in [
            (1000.0, 10.0, 0.0),
            (1000.0, 40.0, 0.0),
            (4000.0, 10.0, 0.0),
            (1000.0, 10.0, 0.4),
        ] {
            let source = material(young, density, fraction);
            let field = TetViscoMaterial::from_source(&source, IDENTITY).unwrap();
            assert_eq!(field.material_state_identity, source.state().identity().0);
            assert_eq!(field.law.band_rad_s(), (0.0, 300.0));
            let solid = ViscoElements::new(
                SolidMeshRef::Tet4 {
                    nodes_m: &NODES,
                    elements: &[[0, 1, 2, 3]],
                },
                TetViscoField::Uniform(&field),
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                TetAssemblyBudget::standard(),
                ViscoStateBudget {
                    max_internal_variables: 100,
                },
                cx,
            )
            .unwrap();
            let basis = solid
                .project_basis(
                    &[vec![1.0]],
                    SolidProjectionBudget {
                        max_coordinates: 1,
                        max_basis_scalars: 1,
                        max_material_points: 1,
                        max_scalar_terms: 100_000,
                    },
                    cx,
                )
                .unwrap();
            let body = TensorSolidBody::new(basis, dynamics(), gate).unwrap();
            // On a unit tetrahedron, integral z² dV = 1/60 and volume = 1/6.
            // Only node 3's z displacement is free; nu=0 makes K=E/6 exactly.
            near(body.body().mass_matrix_kg()[0], density / 60.0, 1e-13);
            let relaxed = body
                .history(&[0.002], &[0.0], InitialMemory::Relaxed, gate)
                .unwrap();
            near(
                relaxed.material_history().stored_energy_j(),
                0.5 * young / 6.0 * 0.002_f64.powi(2),
                1e-13,
            );
            let mut history = body
                .history(&[0.002], &[0.0], InitialMemory::Instantaneous, gate)
                .unwrap();
            let record = history.advance(&[[0.0; 3]; 4], 0.003, 1e-10, gate).unwrap();
            let motion = history.nodal_displacement_m(gate).unwrap()[3][2];
            motions.push(motion);
            near(history.material_history().balance_residual_j(), 0.0, 1e-11);
            near(
                record.element_heat_j.iter().sum(),
                record.dynamics.dissipated,
                1e-11,
            );
            if fraction == 0.0 {
                let omega = fs_math::det::sqrt(10.0 * young / density);
                near(motion, 0.002 * fs_math::det::cos(omega * 0.003), 1e-12);
                near(record.dynamics.dissipated, 0.0, 1e-14);
            } else {
                assert!(record.element_heat_j[0] > 0.0);
            }
        }
        assert!(motions[1] > motions[0] && motions[0] > motions[2]);
        assert!((motions[3] - motions[0]).abs() > 1e-7);
        assert_eq!(
            material(1000.0, 10.0, 0.0).state().material(),
            material(4000.0, 10.0, 0.0).state().material()
        );
        assert_ne!(
            material(1000.0, 10.0, 0.0).state().card_identity(),
            material(4000.0, 10.0, 0.0).state().card_identity()
        );
    });
}

fn membrane(
    source: ResolvedTensorMaterial,
    axes: [[f64; 3]; 3],
    rotated: bool,
    reordered: bool,
    mass: MembraneMass,
    gate: &CancelGate,
) -> TensorMembrane {
    let turn = |v: [f64; 3]| [-v[1], v[0], v[2]];
    let nodes = NODES[..3]
        .iter()
        .copied()
        .map(|p| if rotated { turn(p) } else { p })
        .collect();
    let mesh = ShellMesh::new(nodes, vec![if reordered { [1, 2, 0] } else { [0, 1, 2] }]).unwrap();
    let fixed: &[usize] = if rotated {
        &[0, 1, 2, 3, 5, 8]
    } else {
        &[0, 1, 2, 4, 5, 8]
    };
    TensorMembrane::from_material_sources(
        mesh,
        vec![MembraneMaterialSource {
            material: source,
            thickness_m: 0.1,
            material_to_global: axes,
        }],
        fixed,
        mass,
        projection(),
        dynamics(),
        gate,
    )
    .unwrap()
}

#[test]
fn complete_anisotropic_memory_follows_material_axes_not_triangle_order() {
    let gate = CancelGate::new_clock_free();
    let source = material(1000.0, 10.0, 0.4);
    for mass in [MembraneMass::Consistent, MembraneMass::Lumped] {
        let a = membrane(source.clone(), IDENTITY, false, false, mass, &gate);
        let b = membrane(source.clone(), TURN, true, false, mass, &gate);
        let c = membrane(source.clone(), IDENTITY, false, true, mass, &gate);
        assert_eq!(a.material_sources().unwrap()[0].material, source);
        for node_mass in a.nodal_masses_kg() {
            near(*node_mass, 1.0 / 6.0, 1e-14);
        }
        let zero = [[0.0; 3]; 3];
        let mut x = a
            .history(&zero, &zero, InitialMemory::Relaxed, &gate)
            .unwrap();
        let mut y = b
            .history(&zero, &zero, InitialMemory::Relaxed, &gate)
            .unwrap();
        let mut z = c
            .history(&zero, &zero, InitialMemory::Relaxed, &gate)
            .unwrap();
        let turn = |v: [f64; 3]| [-v[1], v[0], v[2]];
        let load = [[0.0; 3], [0.01, 0.0, 0.0], [0.003, -0.002, 0.0]];
        let f = a.project_nodal_forces(&load).unwrap();
        let g = b.project_nodal_forces(&load.map(turn)).unwrap();
        let h = c.project_nodal_forces(&load).unwrap();
        for _ in 0..8 {
            x.advance(&f, 0.002, &gate).unwrap();
            y.advance(&g, 0.002, &gate).unwrap();
            z.advance(&h, 0.002, &gate).unwrap();
        }
        let original = a.expand_coordinates(x.displacement_m()).unwrap();
        let rotated = b.expand_coordinates(y.displacement_m()).unwrap();
        let reordered = c.expand_coordinates(z.displacement_m()).unwrap();
        for ((a, b), c) in original.into_iter().zip(rotated).zip(reordered) {
            for (a, c) in a.into_iter().zip(c) {
                near(a, c, 1e-11);
            }
            for (a, b) in turn(a).into_iter().zip(b) {
                near(a, b, 1e-11);
            }
        }
        near(x.dissipated_heat_j(), y.dissipated_heat_j(), 1e-11);
        near(x.dissipated_heat_j(), z.dissipated_heat_j(), 1e-11);
        assert!(x.dissipated_heat_j() > 0.0);
        near(x.balance_residual_j(), 0.0, 1e-11);
        // Rotate only the anisotropic material while holding geometry and load
        // fixed: this changes physical branch storage and must remain visible
        // to the mechanical operator, rather than being erased by facet axes.
        let turned = membrane(source.clone(), TURN, false, false, mass, &gate);
        let u = [[0.0; 3], [0.001, 0.0, 0.0], [0.0; 3]];
        let e0 = a
            .history(&u, &zero, InitialMemory::Instantaneous, &gate)
            .unwrap()
            .stored_energy_j();
        let e1 = turned
            .history(&u, &zero, InitialMemory::Instantaneous, &gate)
            .unwrap()
            .stored_energy_j();
        assert!((e0 - e1).abs() > 1e-6);
    }
}

#[test]
fn unsupported_density_sources_incomplete_assignments_and_frames_refuse() {
    let (claims, properties) = fixture(1000.0, 10.0, 0.4);
    let bind = |claims, properties| {
        ResolvedTensorMaterial::resolve(
            &card(claims),
            &query(),
            properties,
            MaterialPropertySelection::SingleClaimOnly,
        )
    };
    let mut missing = claims.clone();
    missing.remove(0);
    assert!(bind(missing, properties.clone()).is_err());
    let mut wrong_units = claims.clone();
    wrong_units[0] = scalar("density", Pressure::DIMS, 10.0);
    assert!(bind(wrong_units, properties.clone()).is_err());
    let mut dispersive = claims.clone();
    dispersive[0].value = PropertyValue::Curve {
        abscissa: "omega".into(),
        abscissa_dims: Dims([0, 0, -1, 0, 0, 0]),
        knots: vec![(0.0, 10.0), (300.0, 12.0)],
        dims: Density::DIMS,
    };
    dispersive[0].interpolation = InterpolationPolicy::LinearInside;
    assert!(bind(dispersive, properties.clone()).is_err());
    let mut incomplete = claims;
    incomplete.pop();
    assert!(bind(incomplete, properties).is_err());
    let source = material(1000.0, 10.0, 0.4);
    let reflected = [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    assert!(TetViscoMaterial::from_source(&source, reflected).is_err());
    let gate = CancelGate::new_clock_free();
    let mesh = || ShellMesh::new(NODES[..3].to_vec(), vec![[0, 1, 2]]).unwrap();
    let sources = || {
        vec![MembraneMaterialSource {
            material: source.clone(),
            thickness_m: 0.1,
            material_to_global: IDENTITY,
        }]
    };
    assert!(
        TensorMembrane::from_material_sources(
            mesh(),
            vec![],
            &[],
            MembraneMass::Consistent,
            projection(),
            dynamics(),
            &gate
        )
        .is_err()
    );
    let mut bad_frame = sources();
    bad_frame[0].material_to_global = reflected;
    assert!(
        TensorMembrane::from_material_sources(
            mesh(),
            bad_frame,
            &[0, 1, 2, 4, 5, 8],
            MembraneMass::Consistent,
            projection(),
            dynamics(),
            &gate
        )
        .is_err()
    );
    gate.request();
    assert!(
        TensorMembrane::from_material_sources(
            mesh(),
            sources(),
            &[0, 1, 2, 4, 5, 8],
            MembraneMass::Consistent,
            projection(),
            dynamics(),
            &gate
        )
        .is_err()
    );
}
