//! Original nonlinear Gauss fields feed the same native solid volume integral.
//! Independent P1 source and affine-work oracles; no continuum Poisson claim.
use super::*;
use fs_couple::electrostatic::spatial_nonlinear_coulomb_volume_load;
use fs_em::electrostatic::nonlinear::{
    ElectrodeControl, NonlinearElectrostaticControls,
    spatial::{SpatialDielectricMaterial, SpatialNonlinearElectrostatics},
};
use fs_material::field_storage::dielectric::reversible::{
    ELECTRIC_DISPLACEMENT_DIMS, ELECTRIC_FIELD_DIMS, REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
    ReversibleDielectricMaterial,
};
use fs_qty::QuantitySpec;

const GAP: f64 = 1.0 / 16.0;
const WIDTH: f64 = 1.0 / 8.0;
const DEPTH: f64 = 1.0 / 4.0;
const E0: f64 = 1e4;
const D0: f64 = 1e-6;
const ORIGIN: [f64; 3] = [0.125, -0.25, 0.5];

fn nonlinear_material() -> ReversibleDielectricMaterial {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
                ELECTRIC_DISPLACEMENT_DIMS,
            ),
            value: PropertyValue::Curve {
                abscissa: "E".into(),
                abscissa_dims: ELECTRIC_FIELD_DIMS,
                knots: [(0.0, 0.0), (1.0, 1.0), (3.0, 5.0), (4.0, 9.0)]
                    .map(|(e, d)| (e * E0, d * D0))
                    .to_vec(),
                dims: ELECTRIC_DISPLACEMENT_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "E",
                QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS),
                0.0,
                4.0 * E0,
            ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "authored nonlinear dielectric for solid volume-load oracles".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic nonlinear dielectric".into(),
            phase: "solid".into(),
            process: "authored D(E) load fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new()
        .with_quantity("E", QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS), 0.0)
        .unwrap();
    ReversibleDielectricMaterial::resolve(
        &card,
        &query,
        "E",
        4.0 * E0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn nonlinear_problem(nodes: &[[f64; 3]], tets: &[[usize; 4]]) -> SpatialNonlinearElectrostatics {
    SpatialNonlinearElectrostatics::new(
        nodes.to_vec(),
        tets.iter().map(|tet| tet.map(|v| v as u32)).collect(),
        vec![SpatialDielectricMaterial::Reversible(nonlinear_material()); tets.len()],
        WeightedAssemblyLimits {
            max_cells: 8,
            max_dofs: 64,
            max_triplets: 288,
        },
        &mut || false,
    )
    .unwrap()
}

fn nonlinear_controls() -> NonlinearElectrostaticControls {
    NonlinearElectrostaticControls {
        absolute_residual_tolerance_c: 1e-20,
        relative_residual_tolerance: 1e-12,
        charge_tolerance_c: 1e-18,
        energy_tolerance_j: 1e-13,
        max_newton_iterations: 32,
        max_backtracks: 24,
        max_operator_applications: 2048,
        linear: controls().linear,
    }
}

fn tetrahedron(axis: usize) -> [[f64; 3]; 4] {
    let mut nodes = [ORIGIN; 4];
    nodes[1][axis] += GAP;
    nodes[2][(axis + 1) % 3] += WIDTH;
    nodes[3][(axis + 2) % 3] += DEPTH;
    nodes
}

fn faces_at(problem: &SpatialNonlinearElectrostatics, axis: usize, value: f64) -> Vec<usize> {
    problem
        .exterior_faces()
        .into_iter()
        .filter(|&face| {
            problem.mesh().faces[face]
                .iter()
                .all(|&node| problem.positions_m()[node as usize][axis] == value)
        })
        .collect()
}

#[test]
fn nonlinear_charged_tetrahedron_has_exact_body_force_moment_and_affine_work() {
    let tets = [[0, 1, 2, 3]];
    let volume = GAP * WIDTH * DEPTH / 6.0;
    let moment_origin = [-0.3, 0.25, -0.1];
    for axis in 0..3 {
        let nodes = tetrahedron(axis);
        let problem = nonlinear_problem(&nodes, &tets);
        let ground = faces_at(&problem, axis, ORIGIN[axis]);
        assert_eq!(ground.len(), 1);
        let electrodes = problem
            .prepare_electrodes(&[&ground], &mut || false)
            .unwrap();
        for charge_sign in [-1.0, 1.0] {
            // The only free P1 shape function is x/G in the rotated frame.
            // Its weak Gauss equation gives D_axis=-rho*G/4. The authored
            // curve then gives |E|=2 E0 at |D|=3 D0, crossing its first knot.
            // A fixed initial permittivity would incorrectly give |E|=3 E0.
            let rho = charge_sign * 12.0 * D0 / GAP;
            let field = electrodes
                .solve(
                    &[rho],
                    &[ElectrodeControl::Voltage(0.0)],
                    &[0.0; 4],
                    nonlinear_controls(),
                    &mut || false,
                )
                .unwrap();
            assert!(field.report().newton_iterations > 0);
            let mut electric = [0.0; 3];
            electric[axis] = -charge_sign * 2.0 * E0;
            for component in 0..3 {
                near(
                    field.electric_field_v_m()[0][component],
                    electric[component],
                    1e-6,
                );
            }
            near(field.terminal_charge_c()[0], -rho * volume, 1e-18);
            near(field.energy_j(), 3.5 * volume * E0 * D0, 1e-14);
            near(field.coenergy_j(), 2.5 * volume * E0 * D0, 1e-14);
            let expected_force = electric.map(|e| rho * volume * e);
            let force_scale = rho.abs() * volume * 2.0 * E0;
            let centroid: [f64; 3] =
                std::array::from_fn(|i| nodes.iter().map(|x| x[i] / 4.0).sum());
            let expected_moment = cross(
                std::array::from_fn(|i| centroid[i] - moment_origin[i]),
                expected_force,
            );
            let motion = |x: [f64; 3]| {
                [
                    0.001 + 0.02 * x[0] - 0.01 * x[1],
                    -0.002 + 0.03 * x[1] + 0.004 * x[2],
                    0.003 - 0.01 * x[0] + 0.025 * x[2],
                ]
            };
            with_cx(|cx, _| {
                let geometry =
                    TetLoadGeometry3::new(&nodes, &tets, Default::default(), 4, cx).unwrap();
                let load =
                    spatial_nonlinear_coulomb_volume_load(&field, &geometry, &[0], cx).unwrap();
                let (force, moment) = load.resultant(moment_origin, cx).unwrap();
                for component in 0..3 {
                    near(
                        force[component],
                        expected_force[component],
                        2e-9 * force_scale,
                    );
                    near(
                        moment[component],
                        expected_moment[component],
                        2e-9 * force_scale,
                    );
                }
                let mut nodal = [[0.0; 3]; 4];
                for &(dof, force) in load.nodal_forces_n() {
                    nodal[dof / 3][dof % 3] = force;
                }
                // The three grounded-face nodes retain their mechanical load;
                // the adapter does not drop reactions at any electrical DOF.
                for force in nodal {
                    for component in 0..3 {
                        near(
                            force[component],
                            expected_force[component] / 4.0,
                            2e-9 * force_scale,
                        );
                    }
                }
                let work: f64 = load
                    .nodal_forces_n()
                    .iter()
                    .map(|&(dof, force)| force * motion(nodes[dof / 3])[dof % 3])
                    .sum();
                near(
                    work,
                    dot(expected_force, motion(centroid)),
                    2e-10 * force_scale,
                );
                let empty =
                    spatial_nonlinear_coulomb_volume_load(&field, &geometry, &[], cx).unwrap();
                assert!(empty.nodal_forces_n().is_empty());
            });
        }
    }
}

#[test]
fn nonlinear_electrode_and_constitutive_charge_do_not_become_body_force() {
    let (nodes, tets) = cube(1, GAP);
    let problem = nonlinear_problem(&nodes, &tets);
    let left = faces_at(&problem, 0, 0.0);
    let right = faces_at(&problem, 0, GAP);
    let electrodes = problem
        .prepare_electrodes(&[&left, &right], &mut || false)
        .unwrap();
    let initial: Vec<_> = nodes.iter().map(|x| 2.0 * E0 * x[0]).collect();
    let field = electrodes
        .solve(
            &[0.0; 6],
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Voltage(2.0 * E0 * GAP),
            ],
            &initial,
            nonlinear_controls(),
            &mut || false,
        )
        .unwrap();
    assert!(field.energy_j() > 0.0);
    assert!(
        field
            .terminal_charge_c()
            .iter()
            .all(|charge| charge.abs() > D0 * GAP * GAP)
    );
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&nodes, &tets, Default::default(), 24, cx).unwrap();
        let load =
            spatial_nonlinear_coulomb_volume_load(&field, &geometry, &[0, 1, 2, 3, 4, 5], cx)
                .unwrap();
        assert!(load.nodal_forces_n().is_empty());
    });
}

#[test]
fn nonlinear_load_mismatch_and_cancellation_preserve_the_accepted_field_for_reuse() {
    let nodes = tetrahedron(0);
    let tets = [[0, 1, 2, 3]];
    let problem = nonlinear_problem(&nodes, &tets);
    let ground = faces_at(&problem, 0, ORIGIN[0]);
    let electrodes = problem
        .prepare_electrodes(&[&ground], &mut || false)
        .unwrap();
    let field = electrodes
        .solve(
            &[12.0 * D0 / GAP],
            &[ElectrodeControl::Voltage(0.0)],
            &[0.0; 4],
            nonlinear_controls(),
            &mut || false,
        )
        .unwrap();
    let original_electric = field.electric_field_v_m().to_vec();
    let original_voltage = field.relative_voltage_v().to_vec();
    let original_charge = field.terminal_charge_c().to_vec();
    let original_energy = field.energy_j();
    let shifted = nodes.map(|x| [x[0] + GAP, x[1], x[2]]);
    let reversed_winding = [[0, 2, 1, 3]];
    let mut original_load = Vec::new();
    with_cx(|cx, gate| {
        for (positions, connectivity) in [(&shifted, &tets), (&nodes, &reversed_winding)] {
            let mismatch =
                TetLoadGeometry3::new(positions, connectivity, Default::default(), 4, cx).unwrap();
            assert!(matches!(
                spatial_nonlinear_coulomb_volume_load(&field, &mismatch, &[0], cx),
                Err(LoadError3::Invalid(_))
            ));
        }
        let geometry = TetLoadGeometry3::new(&nodes, &tets, Default::default(), 4, cx).unwrap();
        assert!(matches!(
            spatial_nonlinear_coulomb_volume_load(&field, &geometry, &[1], cx),
            Err(LoadError3::Invalid(_))
        ));
        assert!(matches!(
            spatial_nonlinear_coulomb_volume_load(&field, &geometry, &[0, 0], cx),
            Err(LoadError3::Budget)
        ));
        let load = spatial_nonlinear_coulomb_volume_load(&field, &geometry, &[0], cx).unwrap();
        original_load = load.nodal_forces_n().to_vec();
        gate.request();
        assert!(matches!(
            spatial_nonlinear_coulomb_volume_load(&field, &geometry, &[0], cx),
            Err(LoadError3::Cancelled)
        ));
    });
    assert_eq!(field.electric_field_v_m(), original_electric.as_slice());
    assert_eq!(field.relative_voltage_v(), original_voltage.as_slice());
    assert_eq!(field.terminal_charge_c(), original_charge.as_slice());
    assert_eq!(field.energy_j(), original_energy);
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&nodes, &tets, Default::default(), 4, cx).unwrap();
        let repeated = spatial_nonlinear_coulomb_volume_load(&field, &geometry, &[0], cx).unwrap();
        assert_eq!(repeated.nodal_forces_n(), original_load.as_slice());
    });
}
