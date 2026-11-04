//! E4.6a-ii battery (bead wf-root-guzez.5.13.2): force build-up
//! cross-checks (per-axis oracles, never totals-only), fixed-control trim
//! within physical bands, determinism, the untrimmable falsifier (canard
//! authority removed -> typed trim-not-found), state-domain refusals,
//! golden trim state. Repro: cargo test -p fs-flyer --test aircraft_battery

use fs_flyer::aircraft::{MAX_TRIM_ITERATIONS, wright_openloop_v1};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-e46aii\",\"case\":\"{case}\",{payload}}}");
}

const RHO: f64 = 1.294;

#[test]
fn g3_resolved_components_feed_pitch_roll_and_yaw_from_their_own_incident_flow() {
    use fs_flyer::dragledger::{
        ComponentDragGeometry, ComponentDragShape, DragComponent, DragLedger,
    };
    let mut design = wright_openloop_v1();
    design.coupling = fs_flyer::propcoupling::CANDIDATE_B;
    design.ledger = DragLedger {
        components: [("left", 1.2), ("right", 0.8)]
            .map(|(component_id, cd)| DragComponent {
                component_id,
                area_m2: 1.0,
                cd,
                cd_source: "synthetic component-force fixture, not historical data",
                uncertainty_frac: 0.2,
            })
            .to_vec(),
        unresolved_interference_frac: 0.2,
    };
    let aggregate_identity = design.digest();
    let ledger_identity = design.ledger.digest();
    let arms = [[-0.7, 0.9, -0.25], [0.8, -1.3, 0.35]];
    let areas = [[0.1, 0.3, 0.2], [0.15, 0.25, 0.4]];
    design.drag_geometry = design
        .ledger
        .components
        .iter()
        .enumerate()
        .map(|(i, c)| ComponentDragGeometry {
            component_id: c.component_id,
            point_body_m: std::array::from_fn(|j| design.cg_m[j] + arms[i][j]),
            shape: ComponentDragShape::OrthogonalFaces {
                normals_body: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                areas_m2: areas[i],
            },
        })
        .collect();
    assert_ne!(aggregate_identity, design.digest());
    assert_eq!(ledger_identity, design.ledger.digest());
    let mut legacy = design.clone();
    legacy.drag_geometry.clear();
    assert_eq!(aggregate_identity, legacy.digest());
    let (speed, alpha, pitch) = (13.86, 0.06, 0.12);
    let velocity = [
        speed * fs_math::det::cos(alpha),
        0.0,
        speed * fs_math::det::sin(alpha),
    ];
    let gusts = [[0.8, -1.5, 0.2], [-0.2, 0.5, -0.1]];
    let run = |co_moving| {
        design
            .force_buildup_with_incident(speed, alpha, 0.03, 40.0, pitch, RHO, None, &|point| {
                Ok(
                    match design
                        .drag_geometry
                        .iter()
                        .position(|geometry| geometry.point_body_m == point)
                    {
                        Some(_) if co_moving => [
                            -pitch * (point[2] - design.cg_m[2]),
                            0.0,
                            pitch * (point[0] - design.cg_m[0]),
                        ],
                        Some(i) => std::array::from_fn(|j| velocity[j] + gusts[i][j]),
                        None => velocity,
                    },
                )
            })
            .unwrap()
    };
    let unloaded = run(true);
    let loaded = run(false);
    assert_eq!(unloaded.parasite_wrench, Default::default());
    assert_eq!(unloaded.coupled.gamma, loaded.coupled.gamma);
    assert_eq!(unloaded.thrust_n, loaded.thrust_n);
    let (mut expected_force, mut expected_moment, mut expected_power) = ([0.0; 3], [0.0; 3], 0.0);
    for i in 0..2 {
        let local = [
            velocity[0] + gusts[i][0] + pitch * arms[i][2],
            gusts[i][1],
            velocity[2] + gusts[i][2] - pitch * arms[i][0],
        ];
        let force: [f64; 3] = std::array::from_fn(|j| {
            -0.5 * RHO
                * design.ledger.components[i].cd
                * areas[i][j]
                * local[j].abs()
                * local[j]
                * 1.2
        });
        let moment = [
            arms[i][1] * force[2] - arms[i][2] * force[1],
            arms[i][2] * force[0] - arms[i][0] * force[2],
            arms[i][0] * force[1] - arms[i][1] * force[0],
        ];
        for j in 0..3 {
            expected_force[j] += force[j];
            expected_moment[j] += moment[j];
            expected_power -= force[j] * local[j];
        }
    }
    for j in 0..3 {
        assert!((loaded.force_n[j] - unloaded.force_n[j] - expected_force[j]).abs() < 1e-9);
        assert!((loaded.parasite_wrench.moment_nm[j] - expected_moment[j]).abs() < 1e-9);
    }
    assert!((loaded.moment_y_nm - unloaded.moment_y_nm - expected_moment[1]).abs() < 1e-9);
    assert!((loaded.parasite_wrench.dissipated_power_w - expected_power).abs() < 1e-9);
    let before = design
        .distributed_lateral_loads(&unloaded, 0.03, RHO)
        .unwrap()
        .net();
    let after = design
        .distributed_lateral_loads(&loaded, 0.03, RHO)
        .unwrap()
        .net();
    assert!(expected_moment[0].abs() > 0.01 && expected_moment[2].abs() > 0.1);
    assert!((after.roll_nm - before.roll_nm - expected_moment[0]).abs() < 1e-9);
    assert!((after.yaw_nm - before.yaw_nm - expected_moment[2]).abs() < 1e-9);
}

#[test]
fn g3_parasite_line_samples_local_gust_and_rigid_pitch_velocity() {
    let mut design = wright_openloop_v1();
    design.coupling = fs_flyer::propcoupling::CANDIDATE_B;
    design.drag_z_m = design.cg_m[2] - 0.35;
    let alpha = 0.06;
    let speed = 13.86;
    let pitch_rate = 0.12;
    let velocity = [
        speed * fs_math::det::cos(alpha),
        0.0,
        speed * fs_math::det::sin(alpha),
    ];
    let point = [design.cg_m[0], design.cg_m[1], design.drag_z_m];
    let arm = design.drag_z_m - design.cg_m[2];
    let gust = [2.0, -1.5, 0.8];
    let run = |increment: [f64; 3]| {
        design
            .force_buildup_with_incident(speed, alpha, 0.03, 40.0, pitch_rate, RHO, None, &|p| {
                Ok(std::array::from_fn(|i| {
                    velocity[i] + if p == point { increment[i] } else { 0.0 }
                }))
            })
            .unwrap()
    };
    let calm = run([0.0; 3]);
    let local = run(gust);
    assert_eq!(calm.coupled.gamma, local.coupled.gamma);
    assert_eq!(calm.thrust_n, local.thrust_n);
    let expected = |increment: [f64; 3]| {
        let flow = [
            velocity[0] + pitch_rate * arm + increment[0],
            increment[1],
            velocity[2] + increment[2],
        ];
        let norm = flow.iter().map(|v| v * v).sum::<f64>().sqrt();
        let magnitude = design.ledger.evaluate(RHO, norm).unwrap().total_parasite_n;
        (flow.map(|v| -magnitude * v / norm), magnitude)
    };
    let (before, before_magnitude) = expected([0.0; 3]);
    let (after, after_magnitude) = expected(gust);
    for i in 0..3 {
        assert!((local.force_n[i] - calm.force_n[i] - (after[i] - before[i])).abs() < 1e-9);
    }
    assert!((local.moment_y_nm - calm.moment_y_nm - arm * (after[0] - before[0])).abs() < 1e-9);
    assert!((local.drag_n - calm.drag_n - (after_magnitude - before_magnitude)).abs() < 1e-9);
    let before_lateral = design.distributed_lateral_loads(&calm, 0.03, RHO).unwrap();
    let after_lateral = design.distributed_lateral_loads(&local, 0.03, RHO).unwrap();
    assert_eq!(before_lateral.circulation, after_lateral.circulation);
    assert_eq!(before_lateral.profile_drag, after_lateral.profile_drag);
    assert_eq!(before_lateral.propellers, after_lateral.propellers);
    let roll_delta = -arm * (after[1] - before[1]);
    assert!(roll_delta.abs() > 1.0);
    assert!((after_lateral.net().roll_nm - before_lateral.net().roll_nm - roll_delta).abs() < 1e-9);
    assert_eq!(after_lateral.net().yaw_nm, before_lateral.net().yaw_nm);
    let result =
        design.force_buildup_with_incident(speed, alpha, 0.03, 40.0, pitch_rate, RHO, None, &|p| {
            if p == point {
                Err(fs_flyer::Refusal {
                    code: "drag-sample-refused",
                    message: "local test field unavailable".into(),
                    ranked_repairs: vec![],
                })
            } else {
                Ok(velocity)
            }
        });
    assert_eq!(result.unwrap_err().code, "drag-sample-refused");
}

#[test]
fn g0_separated_profile_and_induced_drag_reach_the_actual_aircraft_force_and_work() {
    let mut design = wright_openloop_v1();
    design.coupling = fs_flyer::propcoupling::CANDIDATE_B;
    let alpha = 0.50_f64;
    let canard = 0.10;
    let speed = 13.86;
    let loads = design
        .force_buildup(speed, alpha, canard, 40.0, 0.0, RHO)
        .unwrap();
    let (panels, strips) = design.layout(canard).unwrap();
    let local = loads.coupled.strip_velocity_mps.as_ref().unwrap();
    let cg_velocity = [speed * alpha.cos(), speed * alpha.sin()];
    let mut expected = 0.0;
    let mut attached_only = 0.0;
    let mut profile_power = 0.0;
    let mut kj_power = 0.0;
    let mut kj_force = [0.0; 2];
    let mut separated_count = 0;
    for (s, strip) in strips.iter().enumerate() {
        let a = loads.coupled.strip_alpha_eff_rad[s];
        let weight = if a.abs() <= 0.30 {
            0.0
        } else if a.abs() >= 0.45 {
            1.0
        } else {
            let t = (a.abs() - 0.30) / 0.15;
            t * t * (3.0 - 2.0 * t)
        };
        separated_count += usize::from(weight > 0.0);
        // Independent normal-force -> drag projection of the flat plate,
        // using the final induced incidence rather than geometric alpha.
        let pressure_drag = 1.98 * a.sin().powi(2);
        let cd = design.cd0_wing * (1.0 - weight) + pressure_drag * weight;
        let panel = panels[strip.panel_indices[0]];
        let area = strip.chord_m * (panel.b[1] - panel.a[1]).abs();
        let [u, _, w] = local[s];
        let local_speed = (u * u + w * w).sqrt();
        let q = 0.5 * RHO * local_speed * local_speed;
        let drag = q * area * cd;
        expected += drag;
        attached_only += q * area * design.cd0_wing;
        profile_power += drag * (u * cg_velocity[0] + w * cg_velocity[1]) / local_speed;
        // Lift is normal to the local washed flow, so its work at the CG
        // is generally nonzero. Include that resolved momentum transfer
        // instead of assuming the CG velocity is every strip's velocity.
        for &j in &strip.panel_indices {
            let panel = panels[j];
            let scale = RHO * loads.coupled.gamma[j] * (panel.b[1] - panel.a[1]);
            kj_force[0] -= scale * w;
            kj_force[1] += scale * u;
            kj_power += scale * (-w * cg_velocity[0] + u * cg_velocity[1]);
        }
    }
    // Reconstruct the signed transverse force from accepted local flow and
    // circulation, then project using reference velocity components. This
    // does not call the production projection helper or infer lift from -Fz.
    let transverse_lift = (kj_force[0] * cg_velocity[1] - kj_force[1] * cg_velocity[0]) / speed;
    let q_reference = 0.5 * RHO * speed * speed;
    // Existing biplane convention: S=2bc and AR=b/c, so S*AR=2b^2.
    let drag_denominator =
        2.0 * q_reference * core::f64::consts::PI * design.wing_span_m.powi(2) * design.oswald_e;
    let expected_induced = transverse_lift.powi(2) / drag_denominator;
    let body_normal_induced = kj_force[1].powi(2) / drag_denominator;
    assert!((loads.induced_drag_n - expected_induced).abs() < 1e-9);
    assert!((loads.lift_n + kj_force[1]).abs() < 1e-9);
    assert!(
        loads.induced_drag_n > 1.2 * body_normal_induced,
        "high incidence must expose missing streamwise momentum loss"
    );
    assert!(separated_count > 0);
    assert!(
        expected > 2.0 * attached_only,
        "missing separated drag must be measurable"
    );
    let parasite = design.ledger.evaluate(RHO, speed).unwrap().total_parasite_n;
    let profile = loads.drag_n - parasite - loads.induced_drag_n;
    assert!((profile - expected).abs() < 1e-9);
    let weight = design.gross_mass_kg * 9.80665;
    let aero_x = loads.force_n[0] - loads.thrust_n.iter().sum::<f64>() + weight * alpha.sin();
    let aero_z = loads.force_n[2] - weight * alpha.cos();
    let dissipated_power = -(aero_x * cg_velocity[0] + aero_z * cg_velocity[1]);
    let expected_power = speed * (parasite + loads.induced_drag_n) + profile_power - kj_power;
    assert!(profile_power > 0.0 && dissipated_power > 0.0);
    assert!((dissipated_power - expected_power).abs() < 1e-7);
}

#[test]
fn local_canard_gust_changes_pitch_and_uses_local_force_velocity() {
    let mut design = wright_openloop_v1();
    design.coupling = fs_flyer::propcoupling::CANDIDATE_B;
    let velocity = [
        13.86 * fs_math::det::cos(0.06),
        0.0,
        13.86 * fs_math::det::sin(0.06),
    ];
    let baseline = design
        .force_buildup_with_incident(13.86, 0.06, 0.05, 40.0, 0.0, RHO, None, &|_| Ok(velocity))
        .unwrap();
    let gust = design
        .force_buildup_with_incident(13.86, 0.06, 0.05, 40.0, 0.0, RHO, None, &|point| {
            let mut v = velocity;
            if point[0] > 1.0 {
                v[2] += 0.2;
            }
            Ok(v)
        })
        .unwrap();
    assert!(
        gust.moment_y_nm > baseline.moment_y_nm + 5.0,
        "canard upgust must pitch nose-up"
    );
    assert_ne!(gust.force_n[2].to_bits(), baseline.force_n[2].to_bits());
    let (panels, strips) = design.layout(0.05).unwrap();
    let local = gust.coupled.strip_velocity_mps.as_ref().unwrap();
    let mut lift = 0.0;
    for (s, strip) in strips.iter().enumerate() {
        for &j in &strip.panel_indices {
            let panel = panels[j];
            let cross_z =
                local[s][0] * (panel.b[1] - panel.a[1]) - local[s][1] * (panel.b[0] - panel.a[0]);
            lift -= RHO * gust.coupled.gamma[j] * cross_z;
        }
    }
    assert!((gust.lift_n - lift).abs() < 1.0e-9);
    assert!((gust.lift_n - gust.coupled.wing_lift_n).abs() < 1.0e-9);
}

#[test]
fn distributed_pitch_motion_enters_relative_flow_once() {
    let mut design = wright_openloop_v1();
    design.coupling = fs_flyer::propcoupling::CANDIDATE_B;
    let velocity = [
        13.86 * fs_math::det::cos(0.06),
        0.0,
        13.86 * fs_math::det::sin(0.06),
    ];
    let q = 0.1;
    let rotating = design
        .force_buildup_with_incident(13.86, 0.06, 0.05, 40.0, q, RHO, None, &|_| Ok(velocity))
        .unwrap();
    let explicit_motion = design
        .force_buildup_with_incident(13.86, 0.06, 0.05, 40.0, 0.0, RHO, None, &|point| {
            Ok([
                velocity[0] + q * (point[2] - design.cg_m[2]),
                0.0,
                velocity[2] - q * (point[0] - design.cg_m[0]),
            ])
        })
        .unwrap();
    assert_eq!(
        rotating, explicit_motion,
        "no second pitch-rate twist may accompany q cross r"
    );
}

#[test]
fn uniform_flight_uses_the_same_material_velocity_at_zero_and_finite_pitch_rate() {
    let mut design = wright_openloop_v1();
    design.coupling = fs_flyer::propcoupling::CANDIDATE_B;
    let velocity = [
        13.86 * fs_math::det::cos(0.06),
        0.0,
        13.86 * fs_math::det::sin(0.06),
    ];
    for q in [-0.2, 0.0, 0.2] {
        let uniform = design
            .force_buildup(13.86, 0.06, 0.05, 40.0, q, RHO)
            .unwrap();
        let distributed = design
            .force_buildup_with_incident(13.86, 0.06, 0.05, 40.0, q, RHO, None, &|_| Ok(velocity))
            .unwrap();
        assert_eq!(uniform, distributed);
        assert!(
            (uniform.lift_n - uniform.coupled.wing_lift_n).abs() < 1e-9,
            "circulation and assembled lift must use the same washed velocity at q={q}"
        );
    }
    let stationary = design
        .force_buildup(13.86, 0.06, 0.05, 40.0, 0.0, RHO)
        .unwrap();
    for q in [-1e-8, 1e-8] {
        let moving = design
            .force_buildup(13.86, 0.06, 0.05, 40.0, q, RHO)
            .unwrap();
        assert!((moving.moment_y_nm - stationary.moment_y_nm).abs() < 1e-3);
        for axis in 0..3 {
            assert!(
                (moving.force_n[axis] - stationary.force_n[axis]).abs() < 1e-3,
                "the force law must be continuous through zero pitch rate"
            );
        }
    }
}

#[test]
fn pitch_velocity_matches_moving_geometry_and_reaches_propeller_hubs() {
    let mut design = wright_openloop_v1();
    design.coupling = fs_flyer::propcoupling::CANDIDATE_B;
    // Different hub heights expose the axial part of omega cross r in both
    // rotor solves; a strip-only incidence correction cannot reproduce it.
    design.disks[0].center_m[2] = -1.4;
    design.disks[1].center_m[2] = -0.1;
    let velocity = [
        13.86 * fs_math::det::cos(0.06),
        0.0,
        13.86 * fs_math::det::sin(0.06),
    ];
    let (panels, strips) = design.layout(0.05).unwrap();
    let q = 0.2;
    let loads = design
        .force_buildup(13.86, 0.06, 0.05, 40.0, q, RHO)
        .unwrap();
    let local = loads.coupled.strip_velocity_mps.as_ref().unwrap();
    // Independent central difference of the material point's finite rigid
    // rotation about body y. This checks axial and vertical signs together.
    let point_velocity = |point: [f64; 3]| {
        let arm = [point[0] - design.cg_m[0], point[2] - design.cg_m[2]];
        let epsilon = 1e-5;
        let (s, c) = (q * epsilon).sin_cos();
        let forward = [c * arm[0] + s * arm[1], -s * arm[0] + c * arm[1]];
        let backward = [c * arm[0] - s * arm[1], s * arm[0] + c * arm[1]];
        [
            velocity[0] + (forward[0] - backward[0]) / (2.0 * epsilon),
            0.0,
            velocity[2] + (forward[1] - backward[1]) / (2.0 * epsilon),
        ]
    };
    let mut unwashed = 0;
    for (s, strip) in strips.iter().enumerate() {
        let p = panels[strip.panel_indices[0]];
        let point = core::array::from_fn(|axis| 0.5 * (p.a[axis] + p.b[axis]));
        let expected = point_velocity(point);
        assert!((local[s][2] - expected[2]).abs() < 1e-10);
        // Outside both finite wake windows, the axial velocity is exactly
        // the derivative of body geometry, with no induction to subtract.
        if design.disks.iter().all(|disk| {
            (point[1] - disk.center_m[1]).abs() >= design.rotor.radius_m
                || (point[0] - disk.center_m[0]).abs() >= 2.0 * design.rotor.radius_m
        }) {
            assert!((local[s][0] - expected[0]).abs() < 1e-10);
            unwashed += 1;
        }
    }
    assert!(unwashed >= 4);
    for (k, disk) in design.disks.iter().enumerate() {
        let induced =
            fs_wing::induced_velocity_free(disk.center_m, &panels, &loads.coupled.gamma, velocity);
        let axial = point_velocity(disk.center_m)[0] + induced[0];
        let rotor = fs_airscrew::bemt_solve(&design.rotor, RHO, axial, 40.0).unwrap();
        assert!((loads.thrust_n[k] - rotor.thrust_n).abs() < 1e-7);
        assert!((loads.coupled.torque_nm[k] - rotor.torque_nm).abs() < 1e-7);
        let without_rotation =
            fs_airscrew::bemt_solve(&design.rotor, RHO, velocity[0] + induced[0], 40.0).unwrap();
        assert!(
            (rotor.thrust_n - without_rotation.thrust_n).abs() > 0.05,
            "the regression must expose omitted hub velocity"
        );
    }
}

#[test]
fn force_buildup_per_axis_oracles() {
    let d = wright_openloop_v1();
    let b = d.force_buildup(13.86, 0.06, 0.05, 40.0, 0.0, RHO).unwrap();
    // Each strip's lift/drag is normal/parallel to its own incident velocity.
    // With prop wash the ratios w/u differ across strips, so rotating the
    // summed lift/drag by the CG incidence is not a physical force oracle.
    let (panels, strips) = d.layout(0.05).unwrap();
    let local = b.coupled.strip_velocity_mps.as_ref().unwrap();
    let mut suction_x = 0.0;
    let mut profile_drag = [0.0; 2];
    for (s, strip) in strips.iter().enumerate() {
        let [u, _, z] = local[s];
        let speed2 = u * u + z * z;
        let speed = speed2.sqrt();
        for &j in &strip.panel_indices {
            let panel = panels[j];
            let lift = -RHO * b.coupled.gamma[j] * u * (panel.b[1] - panel.a[1]);
            suction_x += lift * z / u;
        }
        let front = panels[strip.panel_indices[0]];
        let drag = 0.5 * RHO * speed2 * strip.chord_m * front.width_m * d.cd0_wing;
        profile_drag[0] += drag * u / speed;
        profile_drag[1] += drag * z / speed;
    }
    let w = d.gross_mass_kg * 9.80665;
    let alpha = 0.06f64;
    let cg_drag = b.induced_drag_n + d.ledger.evaluate(RHO, 13.86).unwrap().total_parasite_n;
    let fz_ref = -b.lift_n + w * alpha.cos() - profile_drag[1] - cg_drag * alpha.sin();
    assert!(
        (b.force_n[2] - fz_ref).abs() < 1e-9,
        "Fz {} vs reconstructed {fz_ref}",
        b.force_n[2]
    );
    let fx_ref = suction_x + b.thrust_n[0] + b.thrust_n[1]
        - profile_drag[0]
        - cg_drag * alpha.cos()
        - w * alpha.sin();
    assert!(
        (b.force_n[0] - fx_ref).abs() < 1e-9,
        "Fx {} vs reconstructed {fx_ref}",
        b.force_n[0]
    );
    // Physical sanity: lift in the weight class at this state.
    assert!(b.lift_n > 2000.0 && b.lift_n < 6000.0, "lift {}", b.lift_n);
    assert!(b.thrust_n[0] > 20.0 && b.thrust_n[0] < 400.0);
    assert!(b.induced_drag_n > 10.0 && b.induced_drag_n < 300.0);
    jlog(
        "buildup",
        &format!(
            "\"lift\":{},\"drag\":{},\"thrust\":{},\"my\":{}",
            b.lift_n,
            b.drag_n,
            b.thrust_n[0] + b.thrust_n[1],
            b.moment_y_nm
        ),
    );
}

#[test]
fn canard_command_sign_matches_the_frozen_convention() {
    // control-signs-v1: positive canard command -> POSITIVE pitch moment.
    let d = wright_openloop_v1();
    let lo = d.force_buildup(13.86, 0.06, 0.00, 40.0, 0.0, RHO).unwrap();
    let hi = d.force_buildup(13.86, 0.06, 0.10, 40.0, 0.0, RHO).unwrap();
    assert!(
        hi.moment_y_nm > lo.moment_y_nm + 10.0,
        "positive dc must raise M_y: {} -> {}",
        lo.moment_y_nm,
        hi.moment_y_nm
    );
    jlog(
        "control-sign",
        &format!("\"dM\":{}", hi.moment_y_nm - lo.moment_y_nm),
    );
}

#[test]
fn trim_converges_inside_physical_bands() {
    let d = wright_openloop_v1();
    let t = d.trim(RHO, [13.0, 0.06, 0.1, 45.0]).unwrap();
    jlog(
        "trim-state",
        &format!(
            "\"v\":{},\"alpha\":{},\"dc\":{},\"omega\":{},\"iters\":{}",
            t.v_mps, t.alpha_rad, t.delta_canard_rad, t.omega_prop_rad_s, t.iterations
        ),
    );
    assert!(t.iterations <= MAX_TRIM_ITERATIONS);
    for (i, r) in t.residuals.iter().enumerate() {
        assert!(r.abs() < 0.5, "residual {i} = {r}");
    }
    // The torque-balance state must land at a physical prop speed.
    assert!(
        t.omega_prop_rad_s > 25.0 && t.omega_prop_rad_s < 80.0,
        "trim prop speed {} rad/s outside the physical band",
        t.omega_prop_rad_s
    );
    // Physical bands: the Flyer trims slow, nose-up, canard loaded.
    assert!(
        t.v_mps > 9.0 && t.v_mps < 20.0,
        "trim speed {} outside the plausible band",
        t.v_mps
    );
    // MODEL-plausibility band, not a historical claim: the 2*pi
    // thin-airfoil camber closure over-lifts vs the 1901-anchored
    // sections (E4.1: modern-equiv cl 0.659 at 5 deg vs 1.18 predicted),
    // so this model trims ~1 deg lower and ~11% faster than the Dec-17
    // state; measured trim alpha is -0.0149 rad. The historical-match
    // refinement is section-data territory (E4.1 tables), not trim math.
    assert!(
        t.alpha_rad > -0.05 && t.alpha_rad < 0.30,
        "trim alpha {} outside the plausible band",
        t.alpha_rad
    );
    assert!(
        t.delta_canard_rad.abs() < 0.5,
        "canard {}",
        t.delta_canard_rad
    );
    // Lift ~ weight at trim (residual closure already enforces it; this
    // is the human-readable receipt line).
    let w = d.gross_mass_kg * 9.80665;
    assert!((t.buildup.lift_n / w - 1.0).abs() < 0.15);
    jlog(
        "trim",
        &format!(
            "\"v\":{},\"alpha\":{},\"dc\":{},\"omega\":{},\"iters\":{},\"lift\":{},\"thrust\":{}",
            t.v_mps,
            t.alpha_rad,
            t.delta_canard_rad,
            t.omega_prop_rad_s,
            t.iterations,
            t.buildup.lift_n,
            t.buildup.thrust_n[0] + t.buildup.thrust_n[1]
        ),
    );
}

#[test]
fn trim_is_deterministic() {
    let d = wright_openloop_v1();
    let a = d.trim(RHO, [13.0, 0.06, 0.1, 45.0]).unwrap();
    let b = d.trim(RHO, [13.0, 0.06, 0.1, 45.0]).unwrap();
    assert_eq!(a.v_mps.to_bits(), b.v_mps.to_bits());
    assert_eq!(a.alpha_rad.to_bits(), b.alpha_rad.to_bits());
    assert_eq!(a.delta_canard_rad.to_bits(), b.delta_canard_rad.to_bits());
    assert_eq!(a.omega_prop_rad_s.to_bits(), b.omega_prop_rad_s.to_bits());
    jlog("determinism", "\"bitwise\":true");
}

#[test]
fn untrimmable_configuration_refuses_typed() {
    // FALSIFIER: shrink the canard to a sliver — no pitch authority, no
    // equilibrium. The trim must refuse with the typed code, never hand
    // back a best-effort state.
    let mut d = wright_openloop_v1();
    d.canard_span_m = 0.2;
    d.canard_chord_m = 0.05;
    let err = d.trim(RHO, [13.0, 0.06, 0.1, 45.0]).unwrap_err();
    assert_eq!(err.code, "trim-not-found");
    assert!(
        err.message.contains("residual trail") || err.message.contains("|r|"),
        "refusal must carry the residual trajectory: {}",
        err.message
    );
    jlog("falsifier", "\"untrimmable_refused\":true");
}

#[test]
fn state_domain_refusals() {
    let d = wright_openloop_v1();
    assert_eq!(
        d.force_buildup(0.5, 0.06, 0.0, 40.0, 0.0, RHO)
            .unwrap_err()
            .code,
        "state-invalid"
    );
    assert_eq!(
        d.force_buildup(13.0, 0.7, 0.0, 40.0, 0.0, RHO)
            .unwrap_err()
            .code,
        "state-invalid"
    );
    assert_eq!(
        d.force_buildup(13.0, 0.06, 0.0, 200.0, 0.0, RHO)
            .unwrap_err()
            .code,
        "state-invalid"
    );
    assert_eq!(
        d.force_buildup(13.0, 0.06, 0.0, 40.0, 0.0, f64::NAN)
            .unwrap_err()
            .code,
        "state-invalid"
    );
    // Design digest sensitivity: moving the CG moves the digest.
    let mut d2 = wright_openloop_v1();
    d2.cg_m[0] -= 0.1;
    assert_ne!(d2.digest(), d.digest());
    jlog("domain", "\"refusals_typed\":true");
}

#[test]
fn trim_golden_digest() {
    let d = wright_openloop_v1();
    let t = d.trim(RHO, [13.0, 0.06, 0.1, 45.0]).unwrap();
    let mut payload = Vec::new();
    for v in [
        t.v_mps,
        t.alpha_rad,
        t.delta_canard_rad,
        t.buildup.lift_n,
        t.buildup.drag_n,
        t.buildup.thrust_n[0],
    ] {
        payload.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    payload.extend_from_slice(t.design_digest.as_bytes());
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer.e46aii-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    // Golden bump 2026-08-21 (bead guzez.7.2.1): fs-airscrew Prandtl
    // acos + fs-flyer contact tanh routed through det::, unifying the
    // native and wasm lanes; the trim point moved by ulps.
    // Golden bump 2026-10-28: ordinary prop wash now applies disk induction
    // once, matching its momentum root. Aircraft parameters and all physical
    // trim bands remain unchanged.
    // Golden bump 2026-11-03: ordinary loads now use the accepted local
    // washed velocity, matching their circulation solve. The original trim
    // digest was reproduced before this correction; physical bands are intact.
    // Golden bump 2026-11-03: induced drag now uses transverse circulatory
    // lift. The measured trim changes; its physical acceptance bands do not.
    // Golden bump 2026-11-03: the coupling spec v3 binds accepted-circulation
    // continuation, with changed final solver ulps and unchanged trim bands.
    assert_eq!(
        digest, "c4e842d466e5d9a77b5026acbbb94d684e18dd5479f601c5d5d3a5608799f4e0",
        "trim golden moved — determinism regression or an intentional \
         model change requiring the golden-bump protocol"
    );
}
