//! G1/G3/G4: independent 3D series, anisotropic and current-trace solutions.
//! Manufactured source laws test implementation, not empirical material data.
#[path = "spatial_conduction/contacts.rs"]
mod contacts;
#[path = "spatial_conduction/matrix.rs"]
mod matrix;
#[path = "spatial_conduction/motion.rs"]
mod motion;
#[path = "spatial_conduction/terminals.rs"]
mod terminals;

use fs_em::{conduction::ConductionControls, spatial::conduction::*};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{
        ElectricalPropertyBasis, ResolvedConductorState,
        spatial::{ResolvedOrthotropicConductor3, SpatialConductorMaterial},
    },
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};
use fs_solver::CheckedCgConfig;

fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 100,
        max_dofs: 300,
        max_triplets: 3600,
    }
}
fn controls() -> ConductionControls {
    ConductionControls {
        linear: CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 200,
            max_corrections: 2,
        },
        current_tolerance_a: 1e-10,
        power_tolerance_w: 1e-9,
    }
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-9 * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}
fn material(
    basis: ElectricalPropertyBasis,
    values: &[f64],
) -> fs_material::state_point::ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (i, &value) in values.iter().enumerate() {
        let name = if values.len() == 1 {
            basis.property_name().to_owned()
        } else {
            format!("principal_{i}")
        };
        let key = PropertyKey::new(&name, basis.dimensions());
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value: PropertyValue::Scalar {
                    value,
                    dims: basis.dimensions(),
                },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured 3D conductor".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_with_key(&key, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "authored".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn slab(xs: &[f64]) -> (Vec<[f64; 3]>, Vec<[u32; 4]>) {
    let mut positions = Vec::new();
    for &x in xs {
        for z in [0.0, 1.0] {
            for y in [0.0, 1.0] {
                positions.push([x, y, z]);
            }
        }
    }
    let mut tets = Vec::new();
    for i in 0..xs.len() - 1 {
        let a = 4 * i as u32;
        let b = a + 4;
        let c = a + 1;
        let d = b + 1;
        let e = a + 2;
        let f = b + 2;
        let g = a + 3;
        let h = b + 3;
        tets.extend([
            [a, b, d, h],
            [a, d, c, h],
            [a, c, g, h],
            [a, g, e, h],
            [a, e, f, h],
            [a, f, b, h],
        ]);
    }
    (positions, tets)
}
fn ends<C: SpatialConductorMaterial>(
    problem: &SpatialConduction<C>,
    right: f64,
    voltage: Option<f64>,
    current: f64,
) -> SpatialConductionBoundary {
    let mut boundary = SpatialConductionBoundary::default();
    for face in problem.exterior_faces() {
        let xs = problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0]);
        if xs == [0.0; 3] {
            if voltage.is_some() {
                boundary.voltage.push(ConductionVoltageFace {
                    face,
                    voltage_v: [0.0; 3],
                });
            } else {
                boundary.current.push(ConductionCurrentFace {
                    face,
                    inward_current_density_a_m2: -current,
                });
            }
        } else if xs == [right; 3] {
            if current == 0.0 {
                boundary.voltage.push(ConductionVoltageFace {
                    face,
                    voltage_v: [voltage.unwrap(); 3],
                });
            } else {
                boundary.current.push(ConductionCurrentFace {
                    face,
                    inward_current_density_a_m2: current,
                });
            }
        }
    }
    boundary
}

#[test]
fn heterogeneous_volume_current_heat_and_drive_follow_data() {
    for (basis, a, b) in [
        (ElectricalPropertyBasis::Resistivity, 2.0, 5.0),
        (ElectricalPropertyBasis::Conductivity, 0.25, 0.5),
    ] {
        let sources = [a, b].map(|value| {
            ResolvedConductorState::try_from_material_with_basis(material(basis, &[value]), basis)
                .unwrap()
        });
        let (positions, tets) = slab(&[0.0, 1.0, 2.0]);
        let problem = SpatialConduction::new(
            positions,
            tets,
            (0..12).map(|i| sources[i / 6].clone()).collect(),
            limits(),
            &mut || false,
        )
        .unwrap();
        let rho = if basis == ElectricalPropertyBasis::Resistivity {
            [a, b]
        } else {
            [a.recip(), b.recip()]
        };
        let current = 1.25;
        for boundary in [
            ends(&problem, 2.0, Some(current * (rho[0] + rho[1])), 0.0),
            ends(&problem, 2.0, Some(0.0), current),
            ends(&problem, 2.0, None, current),
        ] {
            let field = problem.solve(&boundary, controls(), &mut || false).unwrap();
            for (v, p) in problem.positions_m().iter().enumerate() {
                let expected = current * (rho[0] * p[0].min(1.0) + rho[1] * (p[0] - 1.0).max(0.0));
                near(field.voltage_v[v], expected);
            }
            for i in 0..12 {
                near(field.current_density_a_m2[i][0], -current);
                near(field.current_density_a_m2[i][1], 0.0);
                near(field.current_density_a_m2[i][2], 0.0);
                near(field.electric_field_v_m[i][0], -current * rho[i / 6]);
                near(field.cell_power_w[i], current * current * rho[i / 6] / 6.0);
            }
            near(field.power_w, current * current * (rho[0] + rho[1]));
            near(field.supplied_power_w, field.power_w);
            near(field.current_residual_a, 0.0);
            assert!(field.maximum_free_current_residual_a < 1e-10);
        }
    }
}

#[test]
fn rotated_orthotropic_material_and_geometry_scaling_use_full_tensor() {
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let axes = [[s, s, 0.0], [-s, s, 0.0], [0.0, 0.0, 1.0]];
    for scale in [1.0, 0.25] {
        for basis in [
            ElectricalPropertyBasis::Conductivity,
            ElectricalPropertyBasis::Resistivity,
        ] {
            let values = if basis == ElectricalPropertyBasis::Conductivity {
                [2.0, 8.0, 3.0]
            } else {
                [0.5, 0.125, 1.0 / 3.0]
            };
            let source = ResolvedOrthotropicConductor3::try_from_material(
                material(basis, &values),
                basis,
                std::array::from_fn(|i| {
                    PropertyKey::new(format!("principal_{i}"), basis.dimensions())
                }),
                axes,
            )
            .unwrap();
            let (mut positions, tets) = slab(&[0.0, 1.0]);
            for p in &mut positions {
                for x in p {
                    *x *= scale;
                }
            }
            let problem =
                SpatialConduction::new(positions, tets, vec![source; 6], limits(), &mut || false)
                    .unwrap();
            // E=(-1,2,-3) and sigma=[[5,-3,0],[-3,5,0],[0,0,3]].
            // J=(-11,13,-9), E.J=64 W/m3, independent of mesh scaling.
            let phi = |p: [f64; 3]| 17.0 + p[0] - 2.0 * p[1] + 3.0 * p[2];
            let boundary = SpatialConductionBoundary {
                voltage: problem
                    .exterior_faces()
                    .map(|face| ConductionVoltageFace {
                        face,
                        voltage_v: problem.mesh().faces[face]
                            .map(|v| phi(problem.positions_m()[v as usize])),
                    })
                    .collect(),
                current: vec![],
            };
            let field = problem.solve(&boundary, controls(), &mut || false).unwrap();
            for i in 0..6 {
                for d in 0..3 {
                    near(field.electric_field_v_m[i][d], [-1.0, 2.0, -3.0][d]);
                    near(field.current_density_a_m2[i][d], [-11.0, 13.0, -9.0][d]);
                }
            }
            near(field.power_w, 64.0 * scale.powi(3));
            near(field.supplied_power_w, field.power_w);
        }
    }
}

#[test]
fn source_swap_gauge_refusals_and_cancellation_preserve_input() {
    let source = ResolvedConductorState::try_from_material(material(
        ElectricalPropertyBasis::Conductivity,
        &[2.0],
    ))
    .unwrap();
    let (positions, tets) = slab(&[0.0, 0.5, 1.0]);
    let problem =
        SpatialConduction::new(positions, tets, vec![source; 12], limits(), &mut || false).unwrap();
    let boundary = ends(&problem, 1.0, Some(3.0), 0.0);
    let first = problem.solve(&boundary, controls(), &mut || false).unwrap();
    let replacement = ResolvedConductorState::try_from_material(material(
        ElectricalPropertyBasis::Conductivity,
        &[5.0],
    ))
    .unwrap();
    let changed = problem
        .with_conductors(vec![replacement; 12], limits(), &mut || false)
        .unwrap();
    near(
        changed
            .solve(&boundary, controls(), &mut || false)
            .unwrap()
            .power_w,
        45.0,
    );
    near(first.power_w, 18.0);
    let mut count = 0;
    assert!(matches!(
        problem.solve(&boundary, controls(), &mut || {
            count += 1;
            count > 15
        }),
        Err(SpatialConductionError::Cancelled)
    ));
    assert_eq!(
        problem
            .solve(&boundary, controls(), &mut || false)
            .unwrap()
            .voltage_v,
        first.voltage_v
    );
    let mut unbalanced = ends(&problem, 1.0, None, 1.0);
    unbalanced.current.pop();
    assert!(
        problem
            .solve(&unbalanced, controls(), &mut || false)
            .is_err()
    );
    let mut duplicate = boundary.clone();
    duplicate.voltage.push(duplicate.voltage[0]);
    assert!(
        problem
            .solve(&duplicate, controls(), &mut || false)
            .is_err()
    );
    let internal = problem
        .mesh()
        .faces
        .iter()
        .enumerate()
        .find(|(id, _)| !problem.exterior_faces().any(|x| x == *id))
        .unwrap()
        .0;
    assert!(
        problem
            .solve(
                &SpatialConductionBoundary {
                    voltage: vec![ConductionVoltageFace {
                        face: internal,
                        voltage_v: [0.0; 3]
                    }],
                    current: vec![]
                },
                controls(),
                &mut || false
            )
            .is_err()
    );
    assert!(
        problem
            .with_conductors(
                problem.conductors().to_vec(),
                WeightedAssemblyLimits {
                    max_triplets: 1,
                    ..limits()
                },
                &mut || false
            )
            .is_err()
    );
}
