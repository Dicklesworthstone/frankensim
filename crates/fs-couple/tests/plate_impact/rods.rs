//! G0/G1/G3/G4: original spatial rod inertia coupled to consistent plate mass.

use super::*;
use fs_contact::multicontact::{MultiImpactConfig, MultiImpactDeterminacy};
use fs_couple::plate_impact::multiple::{
    PlateImpactSite, PlateRigidImpactSetBudget, PlateRigidImpactSetControls,
};
use fs_couple::plate_impact::plates::PlateImpactSurface;
use fs_couple::plate_impact::rods::{
    RodPlateImpact, RodPlateImpactBudget, RodPlateImpactError, RodPlateImpactReport, RodPlatePoint,
};
use fs_couple::rod_motion::{RodMotion, RodMotionControls, RodNodalWrench, SpatialRodMotion};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_ga::{So3, So3Tangent, Vec3 as GaVec3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_plate_transient::PlateLoad;
use fs_qty::{Density, Pressure};
use fs_solid::rod::material_field::SpatialMaterialRod;
use fs_solid::rod::material_field::motion::RodRootBoundary;
use fs_solid::rod::material_field::tension::{MaterialTensionLink, TensionLinkAttachment};
use std::sync::Arc;

fn material(density: f64) -> ConcentricShearSection {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [
        ("density", density, Density::DIMS),
        ("young_modulus", 1e4, Pressure::DIMS),
        ("shear_modulus", 4e3, Pressure::DIMS),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic rod/plate fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let source = resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ConcentricElasticSection::try_new(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: 0.05,
            material: source,
        }],
    )
    .unwrap()
    .with_shear_modulus("shear_modulus", 0.8)
    .unwrap()
}
fn eccentric_source() -> SpatialMaterialRod {
    let section = material(1000.0);
    let link = MaterialTensionLink::try_new(
        section.elastic().clone(),
        [
            TensionLinkAttachment {
                node: 0,
                offset_body_m: [0.05, 0.08, -0.03],
            },
            TensionLinkAttachment {
                node: 1,
                offset_body_m: [0.05, 0.08, -0.03],
            },
        ],
        0.9,
        1,
        0.2,
    )
    .unwrap();
    SpatialMaterialRod::try_straight(0.8, vec![section])
        .unwrap()
        .with_tension_links(vec![link], &CancelGate::new())
        .unwrap()
}
fn centered(
    velocity: [[f64; 3]; 2],
    translation: [f64; 3],
    clamped: bool,
    density: f64,
) -> RodMotion {
    let gate = CancelGate::new();
    let model = Arc::new(
        SpatialMaterialRod::try_straight(0.8, vec![material(density)])
            .unwrap()
            .into_motion_model(
                if clamped {
                    RodRootBoundary::Clamped
                } else {
                    RodRootBoundary::Free
                },
                &gate,
            )
            .unwrap(),
    );
    let mut pose = model.initial_pose();
    for p in &mut pose.positions_m {
        *p = add(*p, translation);
    }
    RodMotion::new(model, pose, &velocity, &[[0.0; 3]; 2], 0.7, &gate).unwrap()
}
fn bridge_controls() -> PlateRigidImpactSetControls {
    let single = controls(0.0);
    let mut contact = MultiImpactConfig::default();
    contact.solver.gradient_tolerance = 1e-12;
    PlateRigidImpactSetControls {
        impulse: single.impulse,
        contact,
        point_tolerance_m: single.point_tolerance_m,
        speed_tolerance_m_s: single.speed_tolerance_m_s,
        linear_impulse_tolerance_ns: single.linear_impulse_tolerance_ns,
        angular_impulse_tolerance_nms: single.angular_impulse_tolerance_nms,
        energy_tolerance_j: single.energy_tolerance_j,
    }
}
fn bridge_budget() -> RodPlateImpactBudget {
    RodPlateImpactBudget {
        plate: PlateRigidImpactSetBudget {
            plate: budget(),
            max_contacts: 8,
            max_mobility_products: 10000,
        },
        maximum_rod_nodes: 2,
        max_workspace_components: 20000,
    }
}
fn sites(points: &[[f64; 2]], restitution: f64) -> Vec<PlateImpactSite> {
    points
        .iter()
        .map(|&point_xy| PlateImpactSite {
            point_xy,
            activation_tolerance_m: 1e-9,
            restitution,
        })
        .collect()
}
fn surface<'a>(
    mesh: &'a PlateMesh,
    model: &'a PlateModel,
    full: &'a Csr,
    mass_scale: f64,
    damping: f64,
    reference_height_m: f64,
) -> PlateImpactSurface<'a, NoPlateLoad> {
    PlateImpactSurface {
        dynamics: dynamics(model, mass_scale, damping),
        mesh,
        full_mass: full,
        reference_height_m,
    }
}
fn plate_state(n: usize) -> SecondOrderState {
    SecondOrderState::new(0.7, &vec![0.0; n], &vec![0.0; n], &vec![0.0; n])
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|c| a[c] + b[c])
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn rotate(q: So3, a: [f64; 3]) -> [f64; 3] {
    let r = q.rotate(GaVec3::new(a[0], a[1], a[2])).unwrap();
    [r.x, r.y, r.z]
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn row(point: [f64; 2]) -> Vec<f64> {
    vec![
        1.0 - point[0] - point[1],
        0.0,
        0.0,
        point[0],
        0.0,
        0.0,
        point[1],
        0.0,
        0.0,
    ]
}
fn audit(report: &RodPlateImpactReport) {
    let b = report.balance.as_ref().unwrap();
    close_vec(b.linear_balance_residual, Vec3::ZERO);
    close_vec(b.angular_balance_residual, Vec3::ZERO);
    close(b.work_energy_residual_j, 0.0);
    close(b.energy_residual_j, 0.0);
    close(b.kinetic_energy_change_j, -report.contact.dissipated_energy);
    for (&a, &b) in report
        .relative_normal_speeds_after_m_s
        .iter()
        .zip(&report.contact.relative_normal_speeds_after)
    {
        close(a, b);
    }
}

#[test]
fn full_spatial_mass_and_consistent_plate_mass_share_one_event_and_continue_native_motion() {
    let gate = CancelGate::new();
    let model = Arc::new(
        eccentric_source()
            .into_spatial_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let inertia = model.inertia()[1];
    let center = inertia.center_of_mass_body_m();
    let tensor = inertia.inertia_com_body_kg_m2();
    let spatial = fs_mbd::articulated::SpatialInertia::new(
        inertia.mass_kg(),
        GaVec3::new(center[0], center[1], center[2]),
        Mat3 {
            m: std::array::from_fn(|c| tensor[c / 3][c % 3]),
        },
    )
    .unwrap()
    .matrix()
    .unwrap();
    let arms = [[-0.1, 0.0, 0.0], [0.1, 0.0, 0.0], [0.0, 0.15, 0.0]];
    let points = arms.map(|arm| [0.4 + arm[0], 0.3 + arm[1]]);
    let rows: Vec<_> = points.map(row).into_iter().collect();
    let rod_responses: Vec<_> = arms
        .iter()
        .map(|&arm| {
            let t = cross(arm, [0.0, 0.0, 1.0]);
            let mut response = [t[0], t[1], t[2], 0.0, 0.0, 1.0];
            lu(&spatial.m, 6).unwrap().solve(&mut response);
            response
        })
        .collect();
    let contacts: Vec<_> = arms
        .map(|point_body_m| RodPlatePoint {
            node: 1,
            point_body_m,
        })
        .into_iter()
        .collect();
    let rod_v = [-0.2, 0.1, -0.4];
    let rod_w = [0.1, -0.2, 0.15];
    for restitution in [0.0, 0.4, 1.0] {
        let mut reference = None;
        for scene in 0..2 {
            let rotation = So3::exp(So3Tangent::new(GaVec3::new(
                0.0,
                0.0,
                if scene == 0 { 0.0 } else { 0.7 },
            )))
            .unwrap();
            let translation = if scene == 0 {
                [0.0; 3]
            } else {
                [1.2, -0.7, 2.0]
            };
            let order = if scene == 0 { [0, 1, 2] } else { [2, 0, 1] };
            let (mut mesh, plate_model, full) = fixture(false);
            mesh.nodes = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
                .map(|p| {
                    let p = add(rotate(rotation, p), translation);
                    (p[0], p[1])
                })
                .into_iter()
                .collect();
            let mut pose = model.initial_pose();
            for n in 0..2 {
                pose.positions_m[n] = add(
                    rotate(rotation, add(pose.positions_m[n], [-0.4, 0.3, 0.0])),
                    translation,
                );
                pose.frames[n] = rotation;
            }
            let mut rod = SpatialRodMotion::new_spatial(
                model.clone(),
                pose,
                &[rotate(rotation, rod_v); 2],
                &[rod_w; 2],
                0.7,
                1e-12,
                &gate,
            )
            .unwrap();
            let native_before = rod.clone();
            let ordered: Vec<_> = order
                .map(|i| {
                    let p = add(rod.pose().positions_m[1], rotate(rotation, arms[i]));
                    [p[0], p[1]]
                })
                .into_iter()
                .collect();
            let bridge = RodPlateImpact::prepare(
                surface(&mesh, &plate_model, &full, 1.3, 0.17, translation[2]),
                &sites(&ordered, restitution),
                bridge_controls(),
                bridge_budget(),
                &mut || false,
            )
            .unwrap();
            let dense: Vec<_> = full.to_dense().iter().map(|v| 1.3 * v).collect();
            let mut plate_responses = rows.clone();
            for response in &mut plate_responses {
                lu(&dense, 9).unwrap().solve(response);
            }
            let mut matrix = vec![0.0; 9];
            for i in 0..3 {
                for j in 0..3 {
                    let r = rod_responses[j];
                    matrix[3 * i + j] = r[5]
                        + cross([r[0], r[1], r[2]], arms[i])[2]
                        + dot(&rows[i], &plate_responses[j]);
                }
            }
            // Independent dense simultaneous equation with every row approaching.
            // A positive impulse vector alone does not guarantee approaching
            // velocities when off-diagonal mobility entries have either sign.
            let mut exact = [0.3 * (1.0 + restitution); 3];
            lu(&matrix, 3).unwrap().solve(&mut exact);
            assert!(exact.iter().all(|j| *j > 0.0));
            let target: Vec<_> = (0..3)
                .map(|i| {
                    let before = -dot(&matrix[3 * i..3 * i + 3], &exact) / (1.0 + restitution);
                    assert!(before < 0.0);
                    rod_v[2] + cross(rod_w, arms[i])[2] - before
                })
                .collect();
            let interpolation: Vec<_> = rows.iter().flat_map(|r| [r[0], r[3], r[6]]).collect();
            let mut nodal = target;
            lu(&interpolation, 3).unwrap().solve(&mut nodal);
            let mut plate = plate_state(9);
            for n in 0..3 {
                plate.v[3 * n] = nodal[n];
                plate.a[3 * n] = -0.17 * nodal[n];
            }
            let plate_before = plate.clone();
            let report = bridge
                .apply(&mut rod, &mut plate, &order.map(|i| contacts[i]), &gate)
                .unwrap();
            audit(&report);
            for i in 0..3 {
                for j in 0..3 {
                    close(
                        report.delassus_kg_inverse[3 * i + j],
                        matrix[3 * order[i] + order[j]],
                    );
                }
                close(
                    report.relative_normal_speeds_before_m_s[i],
                    -dot(&matrix[3 * order[i]..3 * order[i] + 3], &exact) / (1.0 + restitution),
                );
                close(report.contact.normal_impulses[i], exact[order[i]]);
            }
            for n in 0..9 {
                let change: f64 = (0..3).map(|j| -exact[j] * plate_responses[j][n]).sum();
                close(plate.v[n], plate_before.v[n] + change);
                close(plate.a[n], -0.17 * plate.v[n]);
            }
            let change: [f64; 6] =
                std::array::from_fn(|c| (0..3).map(|j| exact[j] * rod_responses[j][c]).sum());
            let expected_v = rotate(rotation, add(rod_v, [change[3], change[4], change[5]]));
            let expected_w = add(rod_w, [change[0], change[1], change[2]]);
            for c in 0..3 {
                close(rod.velocities_m_s().unwrap()[1][c], expected_v[c]);
                close(
                    rod.angular_velocities_body_rad_s().unwrap()[1][c],
                    expected_w[c],
                );
            }
            assert_eq!(rod.pose(), native_before.pose());
            assert_eq!(rod.time_s(), native_before.time_s());
            assert_eq!(rod.accepted_steps(), native_before.accepted_steps());
            assert!(std::ptr::eq(rod.model(), native_before.model()));
            let mut only_jump = plate_before.clone();
            only_jump.v = plate.v.clone();
            only_jump.a = plate.a.clone();
            assert_eq!(plate, only_jump);
            if let Some(loss) = reference {
                close(loss, report.contact.dissipated_energy);
            }
            reference = Some(report.contact.dissipated_energy);
            let dt = 1e-6;
            let continuing = rod
                .step_candidate(
                    dt,
                    &[RodNodalWrench::default(); 2],
                    RodMotionControls {
                        maximum_step_s: dt,
                        maximum_steps: 10,
                        energy_defect_limit_j: 1e-6,
                        maximum_rotation_rad: 0.1,
                        balance_origin_m: [0.0; 3],
                    },
                    &gate,
                )
                .unwrap()
                .0;
            OperatorGeneralizedAlpha::new(9, dt, 1.0, ImplicitSolveConfig::default())
                .step(&mut plate, &dynamics(&plate_model, 1.3, 0.17), &[0.0; 9])
                .unwrap();
            close(continuing.time_s(), plate.t);
            assert_eq!(continuing.accepted_steps(), 1);
        }
    }
}

#[test]
fn fixed_rod_attachment_and_plate_supports_take_original_reactions_without_work() {
    let gate = CancelGate::new();
    let (mut mesh, model, full) = fixture(true);
    mesh.nodes = vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
    let points = [[0.25, 0.25]];
    let bridge = RodPlateImpact::prepare(
        surface(&mesh, &model, &full, 1.7, 0.0, 0.0),
        &sites(&points, 0.4),
        bridge_controls(),
        bridge_budget(),
        &mut || false,
    )
    .unwrap();
    let mut rod = centered([[0.0; 3]; 2], [0.0; 3], true, 1000.0);
    let before = rod.clone();
    let mut plate = plate_state(6);
    plate.v[0] = 0.5;
    plate.v[3] = 0.5;
    let initial = plate.clone();
    let report = bridge
        .apply(
            &mut rod,
            &mut plate,
            &[RodPlatePoint {
                node: 0,
                point_body_m: [0.25, 0.25, 0.0],
            }],
            &gate,
        )
        .unwrap();
    audit(&report);
    assert_eq!(rod.linear_momenta_kg_m_s(), before.linear_momenta_kg_m_s());
    assert_eq!(rod.body_spins_kg_m2_s(), before.body_spins_kg_m2_s());
    close(report.rod.as_ref().unwrap().kinetic_energy_change_j, 0.0);
    let j = report.contact.normal_impulses[0];
    let support = &report.rod.as_ref().unwrap().support_impulses[0];
    close(support.linear_n_s[2], -j);
    close(support.angular_nm_s[0], -0.25 * j);
    close(support.angular_nm_s[1], 0.25 * j);
    let mut full_delta = [0.0; 9];
    for (r, &free) in model.dof_map.iter().enumerate() {
        if let Some(free) = free {
            full_delta[r] = plate.v[free] - initial.v[free];
        }
    }
    let dense = full.to_dense();
    for &(r, reaction) in &report.balance.as_ref().unwrap().plate_support_impulses {
        let expected = 1.7 * dot(&dense[9 * r..9 * r + 9], &full_delta) + j * row(points[0])[r];
        close(reaction, expected);
    }
    assert!(
        report
            .balance
            .as_ref()
            .unwrap()
            .plate_support_impulses
            .iter()
            .any(|(_, v)| v.abs() > 1e-3)
    );
}

#[test]
fn distinct_rod_nodes_couple_through_the_plate_and_redundant_rows_retain_determinacy() {
    let gate = CancelGate::new();
    let (mut mesh, model, full) = fixture(false);
    mesh.nodes = vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
    let base_points = [[0.05, 0.1], [0.85, 0.1]];
    let base_rows = base_points.map(row);
    let mut responses = base_rows.clone();
    for r in &mut responses {
        lu(&full.to_dense(), 9).unwrap().solve(r);
    }
    let initial = centered(
        [[0.0, 0.0, -1.0], [0.0, 0.0, -0.5]],
        [0.05, 0.1, 0.0],
        false,
        1000.0,
    );
    let m = initial.model().inertia()[0].mass_kg();
    let matrix = [
        1.0 / m + dot(&base_rows[0], &responses[0]),
        dot(&base_rows[0], &responses[1]),
        dot(&base_rows[1], &responses[0]),
        1.0 / m + dot(&base_rows[1], &responses[1]),
    ];
    let mut exact = [1.0, 0.5];
    lu(&matrix, 2).unwrap().solve(&mut exact);
    assert!(exact.iter().all(|j| *j > 0.0));
    let mut reference: Option<(Vec<f64>, Vec<[f64; 3]>)> = None;
    for duplicate in [false, true] {
        let mapping = if duplicate { vec![0, 1, 1] } else { vec![0, 1] };
        let points: Vec<_> = mapping.iter().map(|&i| base_points[i]).collect();
        let contacts: Vec<_> = mapping
            .iter()
            .map(|&node| RodPlatePoint {
                node,
                point_body_m: [0.0; 3],
            })
            .collect();
        let bridge = RodPlateImpact::prepare(
            surface(&mesh, &model, &full, 1.0, 0.0, 0.0),
            &sites(&points, 0.0),
            bridge_controls(),
            bridge_budget(),
            &mut || false,
        )
        .unwrap();
        let mut rod = initial.clone();
        let mut plate = plate_state(9);
        let report = bridge
            .apply(&mut rod, &mut plate, &contacts, &gate)
            .unwrap();
        audit(&report);
        close(report.contact.normal_impulses[0], exact[0]);
        close(report.contact.normal_impulses[1..].iter().sum(), exact[1]);
        if duplicate {
            assert!(matches!(
                report.contact.determinacy,
                MultiImpactDeterminacy::UniqueVelocity(_)
            ));
        }
        close(report.delassus_kg_inverse[1], matrix[1]);
        if let Some((v, p)) = reference {
            for (a, b) in plate.v.iter().zip(v) {
                close(*a, b);
            }
            for (a, b) in rod.linear_momenta_kg_m_s().iter().zip(p) {
                for c in 0..3 {
                    close(a[c], b[c]);
                }
            }
        }
        reference = Some((plate.v.clone(), rod.linear_momenta_kg_m_s().to_vec()));
    }
}

#[test]
fn raw_native_work_retains_small_loss_below_large_common_translation_energy() {
    let gate = CancelGate::new();
    let (mut mesh, model, full) = fixture(false);
    mesh.nodes = vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
    // Keep the native mass-solve work defect below its original 1e-10 J
    // allowance while the total kinetic-energy ulp still exceeds the loss.
    let common = 4_194_304.0;
    let kick = 3.0 / 64.0;
    // The density and original section integration give unit nodal masses.
    let density = 1.0 / (0.4 * (core::f64::consts::PI * 0.05 * 0.05));
    let mut rod = centered(
        [[0.0, 0.0, common], [0.0, 0.0, common - kick]],
        [-0.8, 0.0, 1e12],
        false,
        density,
    );
    assert_eq!(rod.model().inertia()[1].mass_kg(), 1.0);
    let mut plate = plate_state(9);
    for n in [0, 3, 6] {
        plate.v[n] = common;
    }
    plate.v[0] += kick;
    let before = (rod.clone(), plate.clone());
    let bridge = RodPlateImpact::prepare(
        surface(&mesh, &model, &full, 3.0, 0.0, 1e12),
        &sites(&[[0.0, 0.0]], 0.0),
        bridge_controls(),
        bridge_budget(),
        &mut || false,
    )
    .unwrap();
    let report = bridge
        .apply(
            &mut rod,
            &mut plate,
            &[RodPlatePoint {
                node: 1,
                point_body_m: [0.0; 3],
            }],
            &gate,
        )
        .unwrap();
    audit(&report);
    close(report.contact.normal_impulses[0], kick);
    close(report.contact.dissipated_energy, kick * kick);
    let energy = |state: &SecondOrderState| {
        let mut mv = vec![0.0; 9];
        full.spmv(&state.v, &mut mv);
        1.5 * dot(&state.v, &mv)
    };
    let naive = rod.observe([0.0; 3]).unwrap().kinetic_energy_j + energy(&plate)
        - before.0.observe([0.0; 3]).unwrap().kinetic_energy_j
        - energy(&before.1);
    assert!(
        (naive + report.contact.dissipated_energy).abs() > 0.5 * report.contact.dissipated_energy
    );
}

#[test]
fn native_load_refusal_cancellation_geometry_clock_and_budgets_publish_neither_participant() {
    use std::cell::Cell;
    struct Load<'a> {
        mode: Cell<u8>,
        gate: &'a CancelGate,
    }
    impl PlateLoad for Load<'_> {
        fn parameter_count(&self) -> usize {
            0
        }
        fn forcing(&self, _: f64, out: &mut [f64]) -> Result<(), String> {
            match self.mode.get() {
                1 => return Err("native plate load refusal".into()),
                2 => self.gate.request(),
                _ => {}
            }
            out.fill(0.0);
            Ok(())
        }
        fn forcing_vjp(&self, _: f64, _: &[f64], _: &mut [f64]) -> Result<(), String> {
            Ok(())
        }
    }
    let gate = CancelGate::new();
    let load = Load {
        mode: Cell::new(0),
        gate: &gate,
    };
    let (mut mesh, model, full) = fixture(false);
    mesh.nodes = vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
    let dynamics = PlateDynamics::new(
        &model,
        PlateDynamicsParameters {
            stiffness_scale: 1.0,
            mass_scale: 1.0,
            mass_damping_per_s: 0.2,
            stiffness_damping_s: 0.0,
        },
        &load,
        PlateDynamicsBudget {
            max_dofs: 9,
            max_nonzeros: 100,
            max_full_dofs: 9,
            max_load_parameters: 0,
        },
        &mut || false,
    )
    .unwrap();
    let make = |budget| {
        RodPlateImpact::prepare(
            PlateImpactSurface {
                dynamics,
                mesh: &mesh,
                full_mass: &full,
                reference_height_m: 0.0,
            },
            &sites(&[[0.0, 0.0]], 0.3),
            bridge_controls(),
            budget,
            &mut || false,
        )
    };
    let mut short = bridge_budget();
    short.max_workspace_components = 1;
    assert!(matches!(
        make(short),
        Err(RodPlateImpactError::Plate(PlateRigidImpactError::Budget(_)))
    ));
    let mut overflow = bridge_budget();
    overflow.maximum_rod_nodes = usize::MAX;
    assert!(make(overflow).is_err());
    let bridge = make(bridge_budget()).unwrap();
    let initial = centered([[0.0, 0.0, -1.0]; 2], [-0.8, 0.0, 0.0], false, 1000.0);
    let contacts = [RodPlatePoint {
        node: 1,
        point_body_m: [0.0; 3],
    }];
    let mut accepted = initial.clone();
    let mut accepted_plate = plate_state(9);
    bridge
        .apply(&mut accepted, &mut accepted_plate, &contacts, &gate)
        .unwrap();
    for case in 0..6 {
        let mut rod = initial.clone();
        let mut plate = plate_state(9);
        let mut points = contacts;
        match case {
            0 => plate.t += 0.1,
            1 => points[0].node = 2,
            2 => points[0].point_body_m[0] = 0.01,
            3 => points[0].point_body_m[2] = -0.01,
            4 => load.mode.set(1),
            _ => load.mode.set(2),
        }
        let before = (format!("{rod:?}"), plate.clone());
        assert!(bridge.apply(&mut rod, &mut plate, &points, &gate).is_err());
        assert_eq!((format!("{rod:?}"), plate), before);
        load.mode.set(0);
    }
    // A fresh gate retries the same native states after the late cancellation.
    let mut rod = initial.clone();
    let mut plate = plate_state(9);
    bridge
        .apply(&mut rod, &mut plate, &contacts, &CancelGate::new())
        .unwrap();
    assert_eq!(
        rod.linear_momenta_kg_m_s(),
        accepted.linear_momenta_kg_m_s()
    );
    assert_eq!(plate, accepted_plate);
    let mut separating = centered([[0.0, 0.0, 1.0]; 2], [-0.8, 0.0, 0.0], false, 1000.0);
    let before = format!("{separating:?}");
    let mut plate = plate_state(9);
    let original = plate.clone();
    let report = bridge
        .apply(&mut separating, &mut plate, &contacts, &CancelGate::new())
        .unwrap();
    assert!(report.rod.is_none() && report.plate.is_none() && report.balance.is_none());
    assert_eq!(format!("{separating:?}"), before);
    assert_eq!(plate, original);

    // A rounded world point would turn this -0.25 m penetration into a
    // +0.25 m inactive gap at 2^53. Subtract original origins/offsets first.
    let height = 9_007_199_254_740_992.0;
    let shifted = RodPlateImpact::prepare(
        PlateImpactSurface {
            dynamics,
            mesh: &mesh,
            full_mass: &full,
            reference_height_m: height,
        },
        &sites(&[[0.0, 0.0]], 0.3),
        bridge_controls(),
        bridge_budget(),
        &mut || false,
    )
    .unwrap();
    let mut rod = centered([[0.0, 0.0, -1.0]; 2], [-0.8, 0.0, height], false, 1000.0);
    let mut plate = plate_state(9);
    for n in [0, 3, 6] {
        plate.q[n] = 1.75;
    }
    let before = (format!("{rod:?}"), plate.clone());
    assert!(matches!(
        shifted.apply(
            &mut rod,
            &mut plate,
            &[RodPlatePoint {
                node: 1,
                point_body_m: [0.0, 0.0, 1.5],
            }],
            &CancelGate::new(),
        ),
        Err(RodPlateImpactError::Plate(
            PlateRigidImpactError::MultipleContact(_)
        ))
    ));
    assert_eq!((format!("{rod:?}"), plate), before);
}
