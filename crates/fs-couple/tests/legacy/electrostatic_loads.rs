//! G1/G3/G4: accepted dielectric fields feed the native solid volume-load owner.
//! Materials are synthetic SI fixtures, not experimentally validated specimens.

#[path = "electrostatic_loads/nonlinear.rs"]
mod nonlinear;

use fs_blake3::ContentHash;
use fs_couple::electrostatic::spatial_coulomb_volume_load;
use fs_em::{
    electrostatic::ElectrostaticControls,
    spatial::electrostatic::{
        ElectrostaticVoltageFace, SpatialElectrostaticBoundary, SpatialElectrostatics,
    },
};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY},
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};
use fs_solid::{
    linear3::{
        TetElasticMaterial, TetLinearElasticProblem, TetMaterialField, TetThermalStrainField,
        TetThermalStrainState,
    },
    loads3::{LoadError3, TetLoadGeometry3},
};
use fs_solver::CheckedCgConfig;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 73,
                kernel_id: 73,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}
fn problem(nodes: &[[f64; 3]], tets: &[[usize; 4]]) -> SpatialElectrostatics {
    let dims = ELECTRIC_PERMITTIVITY_DIMS;
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(ELECTRIC_PERMITTIVITY_PROPERTY, dims),
            value: PropertyValue::Scalar { value: 2e-8, dims },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic electrostatic load fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "authored uniform permittivity".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let material = resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &[ScalarPropertyRequirement::try_new(
            ELECTRIC_PERMITTIVITY_PROPERTY,
            dims,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    SpatialElectrostatics::new(
        nodes.to_vec(),
        tets.iter().map(|t| t.map(|i| i as u32)).collect(),
        vec![material; tets.len()],
        WeightedAssemblyLimits {
            max_cells: 128,
            max_dofs: 1024,
            max_triplets: 8192,
        },
        &mut || false,
    )
    .unwrap()
}
fn controls() -> ElectrostaticControls {
    ElectrostaticControls {
        linear: CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 1000,
            max_corrections: 2,
        },
        charge_tolerance_c: 1e-15,
        energy_tolerance_j: 1e-12,
    }
}
fn affine_voltage(p: &SpatialElectrostatics, electric: [f64; 3]) -> SpatialElectrostaticBoundary {
    SpatialElectrostaticBoundary {
        voltage: p
            .exterior_faces()
            .map(|face| ElectrostaticVoltageFace {
                face,
                voltage_v: p.mesh().faces[face]
                    .map(|node| -dot(electric, p.positions_m()[node as usize])),
            })
            .collect(),
        displacement: vec![],
    }
}
fn cube(n: usize, length: f64) -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let (mesh, nodes) = fs_feec::fixtures::kuhn_cube(n);
    (
        nodes.into_iter().map(|p| p.map(|x| length * x)).collect(),
        mesh.tets.iter().map(|t| t.map(|i| i as usize)).collect(),
    )
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn g1_constant_p1_charge_force_preserves_si_resultant_moment_and_virtual_work() {
    let nodes = [
        [1.0, -2.0, 0.5],
        [3.0, -2.0, 0.5],
        [1.0, 1.0, 0.5],
        [1.0, -2.0, 4.5],
    ];
    let tets = [[0, 1, 2, 3]];
    let p = problem(&nodes, &tets);
    let electric = [2e3, -3e3, 5e3];
    let origin = [-0.5, 0.75, -1.0];
    // On a single P1 element all nodal voltages prescribe a constant field.
    // The volume source changes electrode reactions. This is an exact load-
    // transfer oracle for that discrete field, not a continuum Poisson claim.
    for charge_sign in [-1.0, 1.0] {
        for field_sign in [-1.0, 1.0] {
            let rho = charge_sign * 1e-9;
            let e = electric.map(|v| field_sign * v);
            let field = p
                .solve(&[rho], &affine_voltage(&p, e), controls(), &mut || false)
                .unwrap();
            with_cx(|cx, _| {
                let geometry =
                    TetLoadGeometry3::new(&nodes, &tets, Default::default(), 4, cx).unwrap();
                let load = spatial_coulomb_volume_load(&field, &geometry, &[0], cx).unwrap();
                let expected = e.map(|v| 4.0 * rho * v); // V = 2*3*4/6 = 4 m3.
                let (force, moment) = load.resultant(origin, cx).unwrap();
                let centroid = [1.5, -1.25, 1.5];
                let expected_moment =
                    cross(std::array::from_fn(|i| centroid[i] - origin[i]), expected);
                for axis in 0..3 {
                    near(force[axis], expected[axis], 1e-17);
                    near(moment[axis], expected_moment[axis], 2e-17);
                }
                let displacement =
                    |p: [f64; 3]| [0.01 * p[0] - 0.02 * p[1], 0.003 + p[2] * 0.04, -p[0] * 0.02];
                let work: f64 = load
                    .nodal_forces_n()
                    .iter()
                    .map(|&(dof, f)| f * displacement(nodes[dof / 3])[dof % 3])
                    .sum();
                near(work, dot(expected, displacement(centroid)), 1e-18);
            });
        }
    }
}

#[test]
fn g1_computed_space_charge_load_drives_supported_native_solid() {
    let length = 0.01;
    let (nodes, tets) = cube(2, length);
    let p = problem(&nodes, &tets);
    let voltage = 2000.0;
    let boundary = SpatialElectrostaticBoundary {
        voltage: p
            .exterior_faces()
            .filter_map(|face| {
                let xs = p.mesh().faces[face].map(|n| p.positions_m()[n as usize][0]);
                if xs == [0.0; 3] {
                    Some(ElectrostaticVoltageFace {
                        face,
                        voltage_v: [voltage; 3],
                    })
                } else if xs == [length; 3] {
                    Some(ElectrostaticVoltageFace {
                        face,
                        voltage_v: [0.0; 3],
                    })
                } else {
                    None
                }
            })
            .collect(),
        displacement: vec![],
    };
    let rho = 0.05;
    let field = p
        .solve(&vec![rho; tets.len()], &boundary, controls(), &mut || false)
        .unwrap();
    assert!(
        field
            .electric_field_v_m()
            .iter()
            .any(|e| (e[0] - voltage / length).abs() > 100.0)
    );
    let fixed: Vec<_> = nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| p[0] == 0.0)
        .flat_map(|(node, _)| (0..3).map(move |axis| 3 * node + axis))
        .collect();
    let material = TetElasticMaterial::try_new(1200.0, 1e7, 0.25, [73; 32]).unwrap();
    let thermal = TetThermalStrainState::try_new(
        300.0,
        300.0,
        [0.0; 6],
        [73; 32],
        ContentHash([74; 32]),
        0.01,
    )
    .unwrap();
    with_cx(|cx, _| {
        let geometry =
            TetLoadGeometry3::new(&nodes, &tets, Default::default(), 4 * tets.len(), cx).unwrap();
        let selected: Vec<_> = (0..tets.len()).collect();
        let load = spatial_coulomb_volume_load(&field, &geometry, &selected, cx).unwrap();
        let solid = TetLinearElasticProblem {
            nodes_m: &nodes,
            tetrahedra: &tets,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &fixed,
            budget: Default::default(),
        };
        let response = load
            .solve(
                &solid,
                TetThermalStrainField::Uniform(&thermal),
                Default::default(),
                ContentHash([75; 32]),
                0.01,
                cx,
            )
            .unwrap();
        // Integral E_x dV = voltage * area, including the solved space field.
        let total_force = rho * voltage * length * length;
        let (force, _) = load.resultant([0.0; 3], cx).unwrap();
        near(force[0], total_force, 1e-12);
        for axis in 0..3 {
            let reaction: f64 = response.support_reactions_n.iter().map(|r| r[axis]).sum();
            near(reaction, -force[axis], 1e-10);
        }
        assert!(response.free_force_residual_n < 1e-10);
        assert!(
            nodes
                .iter()
                .zip(response.displacement.displacement_m())
                .filter(|(p, _)| p[0] == length)
                .all(|(_, u)| u[0] > 0.0)
        );
        let work: f64 = load
            .nodal_forces_n()
            .iter()
            .map(|&(dof, force)| force * response.displacement.displacement_m()[dof / 3][dof % 3])
            .sum();
        assert!(work > 0.0);
    });
}

#[test]
fn g3_cell_selection_does_not_turn_electrode_charge_into_body_charge() {
    let (nodes, tets) = cube(1, 1.0);
    let p = problem(&nodes, &tets);
    let electric = [2.0, -3.0, 5.0];
    let boundary = affine_voltage(&p, electric);
    let rho = [1e-9, 2e-9, 0.0, -3e-9, -2e-9, 4e-9];
    let field = p.solve(&rho, &boundary, controls(), &mut || false).unwrap();
    let uncharged = p
        .solve(&[0.0; 6], &boundary, controls(), &mut || false)
        .unwrap();
    assert!(
        uncharged
            .nodal_boundary_charge_c()
            .iter()
            .flatten()
            .any(|q| *q != 0.0)
    );
    with_cx(|cx, gate| {
        let geometry = TetLoadGeometry3::new(&nodes, &tets, Default::default(), 24, cx).unwrap();
        let load = spatial_coulomb_volume_load(&field, &geometry, &[1, 5], cx).unwrap();
        let (force, _) = load.resultant([0.0; 3], cx).unwrap();
        for axis in 0..3 {
            near(force[axis], (rho[1] + rho[5]) * electric[axis] / 6.0, 1e-22);
        }
        for (field, selected) in [(&field, &[][..]), (&uncharged, &[0, 1, 2, 3, 4, 5][..])] {
            let load = spatial_coulomb_volume_load(field, &geometry, selected, cx).unwrap();
            assert!(load.nodal_forces_n().iter().all(|(_, f)| *f == 0.0));
        }
        for selected in [&[1, 1][..], &[5, 1][..], &[6][..]] {
            assert!(matches!(
                spatial_coulomb_volume_load(&field, &geometry, selected, cx),
                Err(LoadError3::Invalid(_))
            ));
        }
        gate.request();
        assert!(matches!(
            spatial_coulomb_volume_load(&field, &geometry, &[1], cx),
            Err(LoadError3::Cancelled)
        ));
    });
}

#[test]
fn g3_equal_size_different_frame_or_element_order_refuses_load_transfer() {
    let (nodes, tets) = cube(1, 1.0);
    let p = problem(&nodes, &tets);
    let field = p
        .solve(
            &[1e-9; 6],
            &affine_voltage(&p, [2.0, 0.0, 0.0]),
            controls(),
            &mut || false,
        )
        .unwrap();
    let mut reordered = tets.clone();
    reordered.reverse();
    let shifted: Vec<_> = nodes.iter().map(|p| [p[0] + 1.0, p[1], p[2]]).collect();
    with_cx(|cx, _| {
        for (n, t) in [(&nodes, &reordered), (&shifted, &tets)] {
            let geometry = TetLoadGeometry3::new(n, t, Default::default(), 24, cx).unwrap();
            assert!(matches!(
                spatial_coulomb_volume_load(&field, &geometry, &[0], cx),
                Err(LoadError3::Invalid(_))
            ));
        }
    });
}
