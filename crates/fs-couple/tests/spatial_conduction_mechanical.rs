//! G1/G3/G4: original sourced inertia and moving-current feedback, not a motor fit.
#[path = "spatial_conduction_mechanical/response.rs"]
mod response;
use fs_couple::circuit_enthalpy::conduction::spatial::mechanical::{
    PreparedSpatialConductionKick, SpatialConductionKickControls, SpatialConductionKickError,
};
use fs_em::{
    conduction::ConductionControls,
    spatial::conduction::{
        ConductionCurrentFace, ConductionVoltageFace, SpatialConduction, SpatialConductionBoundary,
        motion::SpatialLorentzLoadLimits,
    },
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_ga::{Se3, Twist, Vec3, Wrench};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{ElectricalPropertyBasis, resolve_conductor_state_with_basis},
    state_point::MaterialPropertySelection,
};
use fs_mbd::articulated::FreeFloatingBaseState;
use fs_mbd_material::{
    RegionBinding, ResolvedBody, VolumeMoments, deformation::TetrahedralMass, resolve_body,
};
use fs_qty::Density;
use fs_solver::CheckedCgConfig;

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 2e-8 * a.abs().max(b.abs()).max(1.0),
        "{a:.16e} != {b:.16e}"
    );
}
fn controls() -> SpatialConductionKickControls {
    SpatialConductionKickControls {
        field: ConductionControls {
            linear: CheckedCgConfig {
                tolerance: 1e-12,
                max_iterations: 200,
                max_corrections: 2,
            },
            current_tolerance_a: 1e-9,
            power_tolerance_w: 1e-9,
        },
        loads: SpatialLorentzLoadLimits {
            max_vertices: 100,
            max_cells: 100,
        },
        maximum_field_solves: 7,
        angular_velocity_tolerance_rad_s: 1e-8,
        linear_velocity_tolerance_m_s: 1e-8,
        energy_tolerance_j: 1e-8,
    }
}
fn assembly() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 100,
        max_dofs: 300,
        max_triplets: 3600,
    }
}
fn fixture(sigma: f64, density: f64, xs: &[f64]) -> (SpatialConduction, ResolvedBody) {
    let basis = ElectricalPropertyBasis::Conductivity;
    let mut claims = ClaimSet::new();
    for (key, value) in [
        (
            PropertyKey::new(basis.property_name(), basis.dimensions()),
            sigma,
        ),
        (PropertyKey::new("density", Density::DIMS), density),
    ] {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic electromechanical law".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured conductor".into(),
            phase: "solid".into(),
            process: "test law".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new();
    let source = resolve_conductor_state_with_basis(
        &card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
        basis,
    )
    .unwrap();
    let mut positions = Vec::new();
    for &x in xs {
        for z in [-0.5, 0.5] {
            for y in [-0.5, 0.5] {
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
    let count = tets.len();
    let problem = SpatialConduction::new(
        positions,
        tets,
        vec![source; count],
        assembly(),
        &mut || false,
    )
    .unwrap();
    let keys: Vec<_> = (0..count).map(|i| format!("tet-{i:04}")).collect();
    let regions: Vec<_> = problem
        .mesh()
        .tets
        .iter()
        .enumerate()
        .map(|(i, tet)| RegionBinding {
            key: &keys[i],
            material: &card,
            point: &point,
            selection: MaterialPropertySelection::SingleClaimOnly,
            geometry: VolumeMoments::tetrahedron(tet.map(|v| {
                let p = problem.positions_m()[v as usize];
                Vec3::new(p[0], p[1], p[2])
            }))
            .unwrap(),
        })
        .collect();
    let body = resolve_body(&regions, count, || false).unwrap();
    (problem, body)
}
fn bind_mass<'a>(problem: &SpatialConduction, body: &'a ResolvedBody) -> TetrahedralMass<'a> {
    let keys: Vec<_> = (0..problem.mesh().tets.len())
        .map(|i| format!("tet-{i:04}"))
        .collect();
    let names: Vec<_> = keys.iter().map(String::as_str).collect();
    let tets: Vec<_> = problem
        .mesh()
        .tets
        .iter()
        .map(|tet| tet.map(|v| v as usize))
        .collect();
    TetrahedralMass::bind(body, problem.positions_m(), &tets, &names, 100, || false).unwrap()
}
fn ends(
    problem: &SpatialConduction,
    voltage: Option<f64>,
    current: f64,
) -> SpatialConductionBoundary {
    let mut boundary = SpatialConductionBoundary::default();
    for face in problem.exterior_faces() {
        let xs = problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0]);
        if xs == [-0.5; 3] {
            boundary.voltage.push(ConductionVoltageFace {
                face,
                voltage_v: [0.0; 3],
            });
        } else if xs == [0.5; 3] {
            if let Some(value) = voltage {
                boundary.voltage.push(ConductionVoltageFace {
                    face,
                    voltage_v: [value; 3],
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
fn translating(speed: f64) -> FreeFloatingBaseState {
    FreeFloatingBaseState::new(
        Se3::identity(),
        Twist::new(Vec3::new(0., 0., 0.), Vec3::new(0., speed, 0.)),
    )
}

#[test]
fn uniform_bar_brakes_generates_and_motors_with_actual_sourced_mass_feedback() {
    for xs in [&[-0.5, 0.5][..], &[-0.5, -0.25, 0.25, 0.5][..]] {
        for density in [1.0, 4.0] {
            let (problem, body) = fixture(2.0, density, xs);
            let mass = bind_mass(&problem, &body);
            let flux = vec![[0., 0., 1.5]; problem.mesh().tets.len()];
            for (voltage, speed, external_y) in [
                (0., 1.2, 0.),
                (0.4, 1.2, 0.),
                (2.4, 0.2, 0.),
                (0.4, 1.2, 0.7),
            ] {
                let boundary = ends(&problem, Some(voltage), 0.0);
                let dt = 0.2;
                let mut prepared = PreparedSpatialConductionKick::new(
                    &problem,
                    &mass,
                    &boundary,
                    &flux,
                    Se3::identity(),
                    dt,
                    controls(),
                    &mut || false,
                )
                .unwrap();
                assert_eq!(prepared.mass().reference(), &body);
                let kick = prepared
                    .kick_candidate(
                        translating(speed),
                        Wrench::new(Vec3::new(0., 0., 0.), Vec3::new(0., external_y, 0.)),
                        &mut || false,
                    )
                    .unwrap();
                // Unit-volume bar: Jx=sigma*(B*v_mid-V), Fy=-B*Jx.
                // Original sourced mass is density*1 m3, not an authored gain.
                let drag = 2.0 * 1.5 * 1.5;
                let drive = 2.0 * 1.5 * voltage + external_y;
                let after = ((density - 0.5 * dt * drag) * speed + dt * drive)
                    / (density + 0.5 * dt * drag);
                let mid = 0.5 * (speed + after);
                let current = 2.0 * (1.5 * mid - voltage);
                near(kick.mechanical.base_after.twist_body.linear.y, after);
                near(kick.midpoint_twist_body.linear.y, mid);
                let angular = kick.mechanical.base_after.twist_body.angular;
                for value in [angular.x, angular.y, angular.z] {
                    near(value, 0.0);
                }
                near(kick.lorentz_wrench_body.force.y, -1.5 * current);
                near(
                    kick.mechanical.linear_momentum_change_world.y,
                    density * (after - speed),
                );
                near(kick.joule_heat_j, dt * current * current / 2.0);
                near(kick.terminal_work_j, -dt * voltage * current);
                near(kick.mechanical_work_j, -dt * 1.5 * current * mid);
                near(kick.external_work_j, dt * external_y * mid);
                near(
                    kick.mechanical.kinetic_energy_change,
                    0.5 * density * (after - speed) * (after + speed),
                );
                near(kick.total_energy_residual_j, 0.0);
                near(kick.mechanical_work_residual_j, 0.0);
                assert_eq!(kick.mechanical.base_after.world_from_base, Se3::identity());
                assert_eq!(kick.field_solves, 1);
                for j in kick.field.current_density_a_m2() {
                    near(j[0], current);
                }
                if voltage == 0.0 {
                    assert!(after < speed && kick.mechanical.kinetic_energy_change < 0.);
                }
                if voltage == 0.4 {
                    assert!(kick.terminal_work_j < 0. && kick.joule_heat_j > 0.);
                }
                if voltage == 2.4 {
                    assert!(after > speed && kick.mechanical_work_j > 0.);
                }
            }
        }
    }
}

#[test]
fn original_current_source_adjusts_voltage_as_the_native_body_accelerates() {
    let (problem, body) = fixture(2.0, 3.0, &[-0.5, 0.0, 0.5]);
    let mass = bind_mass(&problem, &body);
    let boundary = ends(&problem, None, 1.25);
    let flux = vec![[0., 0., 2.]; problem.mesh().tets.len()];
    let dt = 0.3;
    let mut prepared = PreparedSpatialConductionKick::new(
        &problem,
        &mass,
        &boundary,
        &flux,
        Se3::identity(),
        dt,
        controls(),
        &mut || false,
    )
    .unwrap();
    let mut state = translating(-0.4);
    for _ in 0..3 {
        let old = state.twist_body.linear.y;
        let kick = prepared
            .kick_candidate(state, Wrench::default(), &mut || false)
            .unwrap();
        let after = old + dt * 2.0 * 1.25 / 3.0;
        let midpoint = 0.5 * (old + after);
        near(kick.mechanical.base_after.twist_body.linear.y, after);
        for j in kick.field.current_density_a_m2() {
            near(j[0], -1.25);
        }
        let voltage = 1.25 / 2.0 + 2.0 * midpoint;
        near(kick.terminal_work_j, dt * 1.25 * voltage);
        near(kick.joule_heat_j, dt * 1.25 * 1.25 / 2.0);
        near(kick.total_energy_residual_j, 0.0);
        state = kick.mechanical.base_after;
    }
}

#[test]
fn non_diagonal_rotational_feedback_uses_the_full_original_moment_and_same_stage_work() {
    let (problem, body) = fixture(6.0, 6.0, &[-0.5, 0.5]);
    let mass = bind_mass(&problem, &body);
    let boundary = SpatialConductionBoundary {
        voltage: problem
            .exterior_faces()
            .map(|face| ConductionVoltageFace {
                face,
                voltage_v: [0.; 3],
            })
            .collect(),
        current: vec![],
    };
    // On these six original tets all voltage DOFs are on the zero-potential
    // boundary. The represented cell current is exactly sigma*(u_centroid x B).
    // Direct centroid moments give D=24*[[1/24,0,1/48],
    // [0,1/24,1/48],[1/48,1/48,1/12]], while sourced cube I_com=identity.
    // Solving (I+0.1D)w1=(I-0.1D)*(1,0,0) independently gives this rational w1.
    let expected = [2377.0 / 2893.0, 10.0 / 2893.0, -220.0 / 2893.0];
    for pose in [
        Se3::identity(),
        Se3::exp(Twist::new(
            Vec3::new(0.2, -0.3, 0.4),
            Vec3::new(2., -1., 3.),
        ))
        .unwrap(),
    ] {
        let mut prepared = PreparedSpatialConductionKick::new(
            &problem,
            &mass,
            &boundary,
            &[[0., 0., 2.]; 6],
            pose,
            0.2,
            controls(),
            &mut || false,
        )
        .unwrap();
        let old =
            FreeFloatingBaseState::new(pose, Twist::new(Vec3::new(1., 0., 0.), Vec3::default()));
        let kick = prepared
            .kick_candidate(old, Wrench::default(), &mut || false)
            .unwrap();
        let spin = kick.mechanical.base_after.twist_body.angular;
        for (a, b) in [spin.x, spin.y, spin.z].into_iter().zip(expected) {
            near(a, b);
        }
        let linear = kick.mechanical.base_after.twist_body.linear;
        for v in [linear.x, linear.y, linear.z] {
            near(v, 0.0);
        }
        near(kick.terminal_work_j, 0.0);
        near(
            kick.joule_heat_j,
            0.5 * (1.0 - expected.iter().map(|w| w * w).sum::<f64>()),
        );
        near(kick.mechanical_work_j, -kick.joule_heat_j);
        near(kick.total_energy_residual_j, 0.0);
        let angular_delta = pose
            .transform_vector(Vec3::new(expected[0] - 1., expected[1], expected[2]))
            .unwrap();
        for (a, b) in [
            kick.mechanical.angular_momentum_change_world.x,
            kick.mechanical.angular_momentum_change_world.y,
            kick.mechanical.angular_momentum_change_world.z,
        ]
        .into_iter()
        .zip([angular_delta.x, angular_delta.y, angular_delta.z])
        {
            near(a, b);
        }
    }
}

#[test]
fn open_circuit_screens_translation_and_fixed_operator_time_refinement_is_second_order() {
    let (problem, body) = fixture(2.0, 3.0, &[-0.5, 0.0, 0.5]);
    let mass = bind_mass(&problem, &body);
    let flux = vec![[0., 0., 1.5]; 12];
    let mut open = PreparedSpatialConductionKick::new(
        &problem,
        &mass,
        &SpatialConductionBoundary::default(),
        &flux,
        Se3::identity(),
        0.2,
        controls(),
        &mut || false,
    )
    .unwrap();
    let unchanged = open
        .kick_candidate(translating(1.2), Wrench::default(), &mut || false)
        .unwrap();
    near(unchanged.mechanical.base_after.twist_body.linear.y, 1.2);
    near(unchanged.joule_heat_j, 0.0);
    for e in unchanged.field.electric_field_v_m() {
        near(e[0], -1.8);
    }
    let short = ends(&problem, Some(0.0), 0.0);
    let mut errors = Vec::new();
    let total_time: f64 = 0.6;
    for count in [4, 8, 16] {
        let mut prepared = PreparedSpatialConductionKick::new(
            &problem,
            &mass,
            &short,
            &flux,
            Se3::identity(),
            total_time / count as f64,
            controls(),
            &mut || false,
        )
        .unwrap();
        let mut state = translating(1.2);
        let mut heat = 0.0;
        for _ in 0..count {
            let kick = prepared
                .kick_candidate(state, Wrench::default(), &mut || false)
                .unwrap();
            heat += kick.joule_heat_j;
            state = kick.mechanical.base_after;
        }
        let speed = state.twist_body.linear.y;
        near(heat, 0.5 * 3.0 * (1.2 - speed) * (1.2 + speed));
        errors.push((speed - 1.2 * (-4.5 * total_time / 3.0).exp()).abs());
    }
    assert!(errors[0] / errors[1] > 3.9 && errors[1] / errors[2] > 3.9);
}

#[test]
fn wrong_sources_pose_budgets_and_cancelled_kicks_refuse_without_poisoning_retry() {
    let (problem, body) = fixture(2.0, 3.0, &[-0.5, 0.5]);
    let mass = bind_mass(&problem, &body);
    let boundary = ends(&problem, Some(0.0), 0.0);
    let flux = [[0., 0., 1.5]; 6];
    let prepare = |control, cancel: &mut dyn FnMut() -> bool| {
        PreparedSpatialConductionKick::new(
            &problem,
            &mass,
            &boundary,
            &flux,
            Se3::identity(),
            0.2,
            control,
            &mut || cancel(),
        )
    };
    assert!(
        prepare(
            SpatialConductionKickControls {
                maximum_field_solves: 6,
                ..controls()
            },
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        prepare(controls(), &mut || true),
        Err(SpatialConductionKickError::Cancelled)
    ));
    let (_, foreign_body) = fixture(2.0, 7.0, &[-0.5, 0.5]);
    let foreign_mass = bind_mass(&problem, &foreign_body);
    assert!(matches!(
        PreparedSpatialConductionKick::new(
            &problem,
            &foreign_mass,
            &boundary,
            &flux,
            Se3::identity(),
            0.2,
            controls(),
            &mut || false,
        ),
        Err(SpatialConductionKickError::Invalid(_))
    ));
    let mut prepared = prepare(controls(), &mut || false).unwrap();
    let mut calls = 0;
    let expected = prepared
        .kick_candidate(translating(1.2), Wrench::default(), &mut || {
            calls += 1;
            false
        })
        .unwrap();
    for cutoff in [1, calls / 2, calls] {
        let mut count = 0;
        assert!(matches!(
            prepared.kick_candidate(translating(1.2), Wrench::default(), &mut || {
                count += 1;
                count == cutoff
            }),
            Err(SpatialConductionKickError::Cancelled)
        ));
    }
    let mut foreign_pose = translating(1.2);
    foreign_pose.world_from_base =
        Se3::exp(Twist::new(Vec3::default(), Vec3::new(1., 0., 0.))).unwrap();
    assert!(
        prepared
            .kick_candidate(foreign_pose, Wrench::default(), &mut || false)
            .is_err()
    );
    assert!(
        prepared
            .kick_candidate(translating(f64::NAN), Wrench::default(), &mut || false)
            .is_err()
    );
    let retry = prepared
        .kick_candidate(translating(1.2), Wrench::default(), &mut || false)
        .unwrap();
    assert_eq!(retry.mechanical, expected.mechanical);
    assert_eq!(retry.field.voltage_v(), expected.field.voltage_v());
    assert_eq!(
        retry.field.current_density_a_m2(),
        expected.field.current_density_a_m2()
    );
}
