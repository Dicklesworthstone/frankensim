//! G0/G1/G3 spatial disk geometry, physical shaft signs and native momentum.
use fs_flyer::aircraft::wright_openloop_v1;
use fs_flyer::propcoupling::spatial::{self, DiskAxes};
use fs_flyer::propcoupling::{
    CANDIDATE_B, IncidentInflow, PropCouplingSolverSpec, PropDisk,
    coupled_prop_airframe_step_with_incident, wash_factor,
};
use fs_flyer::spine::SixDofState;
use fs_flyer::spine::rotors::geared::{GearedAssembly, GearedRotor, GearedState};
use fs_flyer::spine::rotors::propulsion::{
    GearedPropulsion, PropellerHandedness, PropellerMount, PropulsionInputs,
};
use fs_ga::{Mat3, Quat, So3, So3Tangent, Vec3};
use fs_wing::nonlinear::{InfluenceOperator, StripRegime, StripSpec};
use fs_wing::{SurfaceId, flat_surface, induced_velocity_free};

const ZERO: [f64; 3] = [0.0; 3];
fn v(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn a(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn unit(x: [f64; 3]) -> [f64; 3] {
    a(v(x).scale(1.0 / v(x).norm()))
}
fn closure(_: usize, alpha: f64) -> (f64, StripRegime) {
    (
        core::f64::consts::TAU * (alpha + 0.08),
        StripRegime::Attached,
    )
}
fn setup() -> (Vec<fs_wing::Panel>, Vec<StripSpec>, [PropDisk; 2]) {
    let panels = flat_surface(SurfaceId::WingLower, 8.0, 1.5, 0.0, 0.0, 8, 2).unwrap();
    let strips = (0..8)
        .map(|i| StripSpec {
            panel_indices: vec![i, i + 8],
            chord_m: 1.5,
            twist_rad: 0.0,
        })
        .collect();
    let disks = [
        PropDisk {
            center_m: [-2.0, -1.5, 0.0],
            omega_rad_s: 35.0,
        },
        PropDisk {
            center_m: [-2.0, 1.5, 0.0],
            omega_rad_s: 38.0,
        },
    ];
    (panels, strips, disks)
}
fn precise() -> PropCouplingSolverSpec {
    PropCouplingSolverSpec {
        tol: 1e-11,
        cap: 64,
        ..CANDIDATE_B
    }
}

#[test]
fn g3_true_disk_plane_overlap_and_rotation_covariance() {
    let axes = DiskAxes::new([[1.0, 0.0, 0.0], unit([1.0, 0.3, -0.2])]).unwrap();
    let centers = [[0.0, 0.0, 0.0], [-1.0, 0.2, -0.1]];
    assert_eq!(
        axes.wash_factors([0.2, 0.0, 1.01], centers, 1.0).unwrap()[0],
        0.0
    );
    assert!(wash_factor(0.2, 0.0, 1.0) > 0.5); // old span-only gate washed this point
    let point = [0.1, 0.25, 0.3];
    let factors = axes.wash_factors(point, centers, 1.0).unwrap();
    let expected = 1.0 - 0.1 / f64::sqrt(1.01);
    assert!((factors[0] - expected).abs() < 2e-16);
    let r = So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.4))).unwrap();
    let transformed = DiskAxes::new(axes.directions().map(|d| a(r.rotate(v(d)).unwrap()))).unwrap();
    let rotated = transformed
        .wash_factors(
            a(r.rotate(v(point)).unwrap()),
            centers.map(|c| a(r.rotate(v(c)).unwrap())),
            1.0,
        )
        .unwrap();
    for i in 0..2 {
        assert!((rotated[i] - factors[i]).abs() < 8e-16);
    }
    for bad in [ZERO, [2.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0]] {
        assert!(DiskAxes::new([bad, [1.0, 0.0, 0.0]]).is_err());
    }
    assert!(axes.wash_factors(point, centers, 0.0).is_err());
}

#[test]
fn g0_coplanar_axial_limit_matches_native_legacy_but_identifies_spatial_model() {
    let (panels, strips, disks) = setup();
    let velocity = [12.0, 0.0, 0.8];
    let local = vec![velocity; strips.len()];
    let flow = IncidentInflow {
        strip_velocity_mps: &local,
        disk_velocity_mps: [velocity; 2],
    };
    let op = InfluenceOperator::build(&panels, velocity, 1.225).unwrap();
    let rotor = wright_openloop_v1().rotor;
    let spec = precise();
    let old = coupled_prop_airframe_step_with_incident(
        &op, &panels, &strips, &closure, &rotor, &disks, velocity, 1.225, &spec, None, flow,
    )
    .unwrap();
    let axes = DiskAxes::new([[1.0, 0.0, 0.0]; 2]).unwrap();
    let spatial = spatial::coupled_step(
        &axes, &op, &panels, &strips, &closure, &rotor, &disks, velocity, 1.225, &spec, None, flow,
    )
    .unwrap();
    for i in 0..2 {
        assert!((old.thrust_n[i] - spatial.thrust_n[i]).abs() < 1e-9);
        assert!((old.torque_nm[i] - spatial.torque_nm[i]).abs() < 1e-10);
    }
    for (x, y) in old.gamma.iter().zip(&spatial.gamma) {
        assert!((x - y).abs() < 1e-11);
    }
    assert_ne!(old.spec_digest, spatial.spec_digest);
    assert_eq!(old.spec_digest, spec.digest());
}

#[test]
fn g0_canted_disk_inflow_thrust_and_vector_wash_reconstruct_independently() {
    let (panels, strips, mut disks) = setup();
    disks[0].center_m[2] = -0.3;
    disks[1].center_m[2] = 0.2;
    let axes = DiskAxes::new([unit([1.0, 0.08, 0.12]), unit([1.0, -0.04, -0.06])]).unwrap();
    let velocity = [12.0, 0.1, 0.8];
    let local = vec![velocity; strips.len()];
    let disk_flow = [[12.2, 0.3, 0.7], [11.8, -0.2, 0.9]];
    let flow = IncidentInflow {
        strip_velocity_mps: &local,
        disk_velocity_mps: disk_flow,
    };
    let op = InfluenceOperator::build(&panels, velocity, 1.225).unwrap();
    let rotor = wright_openloop_v1().rotor;
    let result = spatial::coupled_step(
        &axes,
        &op,
        &panels,
        &strips,
        &closure,
        &rotor,
        &disks,
        velocity,
        1.225,
        &precise(),
        None,
        flow,
    )
    .unwrap();
    let area = core::f64::consts::PI * rotor.radius_m * rotor.radius_m;
    for k in 0..2 {
        let wing = induced_velocity_free(disks[k].center_m, &panels, &result.gamma, velocity);
        let axial = (v(wing) + v(disk_flow[k])).dot(v(axes.directions()[k]));
        let independent =
            fs_airscrew::bemt_solve(&rotor, 1.225, axial, disks[k].omega_rad_s).unwrap();
        assert!((independent.thrust_n - result.thrust_n[k]).abs() < 1e-9);
        assert!((independent.torque_nm - result.torque_nm[k]).abs() < 1e-10);
        let vi = result.w_slip[k];
        assert!((2.0 * 1.225 * area * vi * (axial + vi) - result.thrust_n[k]).abs() < 1e-9);
    }
    let mut transverse = 0.0_f64;
    for (i, strip) in strips.iter().enumerate() {
        let p = panels[strip.panel_indices[0]];
        let point = (v(p.a) + v(p.b)).scale(0.5);
        let mut expected = v(velocity);
        for k in 0..2 {
            let direction = v(axes.directions()[k]);
            let offset = point - v(disks[k].center_m);
            let axial = offset.dot(direction);
            let radial = (offset - direction.scale(axial)).norm();
            let factor = if radial < rotor.radius_m && axial.abs() < 2.0 * rotor.radius_m {
                1.0 - axial / f64::sqrt(axial * axial + rotor.radius_m * rotor.radius_m)
            } else {
                0.0
            };
            expected = expected + direction.scale(factor * result.w_slip[k]);
        }
        let accepted = v(result.strip_velocity_mps.as_ref().unwrap()[i]);
        assert!((accepted - expected).norm() < 2e-10);
        transverse = transverse.max((accepted.z - velocity[2]).abs());
    }
    assert!(transverse > 1e-3);
}

fn assembly(mass: f64, axes: [[f64; 3]; 2], ratios: [f64; 2]) -> GearedAssembly<2> {
    GearedAssembly::new(
        mass,
        Mat3 {
            m: [1800.0, 10.0, 5.0, 10.0, 420.0, 2.0, 5.0, 2.0, 1850.0],
        },
        core::array::from_fn(|i| GearedRotor {
            axis_body: axes[i],
            ratio: ratios[i],
            axial_inertia_kgm2: [1.5, 2.5][i],
        }),
    )
    .unwrap()
}
fn state(assembly: &GearedAssembly<2>) -> GearedState {
    assembly
        .state_from_speed(
            SixDofState {
                pos_m: [1.0, -2.0, -10.0],
                vel_mps: [12.0, 0.1, 0.8],
                quat: [1.0, 0.0, 0.0, 0.0],
                omega_body: [0.04, -0.02, 0.03],
            },
            0.0,
            35.0,
        )
        .unwrap()
}
fn controls() -> PropulsionInputs {
    PropulsionInputs {
        density_kg_m3: 1.225,
        canard_rad: 0.03,
        motor_nm: 80.0,
        gravity_world_mps2: [0.0, 0.0, 9.80665],
    }
}
fn mounts(axes: [[f64; 3]; 2]) -> [PropellerMount; 2] {
    core::array::from_fn(|i| PropellerMount {
        rotor_index: i,
        thrust_direction_body: axes[i],
        handedness: if i == 0 {
            PropellerHandedness::Positive
        } else {
            PropellerHandedness::Negative
        },
    })
}
fn rotation(s: &GearedState) -> So3 {
    let q = s.rigid.quat;
    So3::try_from_quat(Quat {
        w: q[0],
        x: q[1],
        y: q[2],
        z: q[3],
    })
    .unwrap()
}
fn physical_momentum(assembly: &GearedAssembly<2>, s: &GearedState) -> Vec3 {
    let mut body = assembly.locked_inertia_kgm2().apply(v(s.rigid.omega_body));
    for rotor in assembly.rotors() {
        body = body
            + v(rotor.axis_body)
                .scale(rotor.axial_inertia_kgm2 * rotor.ratio * assembly.speed(s).unwrap());
    }
    rotation(s).rotate(body).unwrap()
        + v(s.rigid.pos_m).cross(v(s.rigid.vel_mps).scale(assembly.mass_kg()))
}

#[test]
fn g0_canted_native_flight_closes_spatial_momentum_and_is_invariant_to_shaft_sign_choice() {
    let mut design = wright_openloop_v1();
    design.coupling = precise();
    let axes = [unit([1.0, 0.08, 0.12]), unit([1.0, -0.04, -0.06])];
    let assembly = assembly(design.gross_mass_kg, axes, [1.0, -1.0]);
    let state = state(&assembly);
    assert!(GearedPropulsion::new(&design, &assembly, mounts(axes)).is_err());
    let model = GearedPropulsion::new_spatial(&design, &assembly, mounts(axes)).unwrap();
    let dt = 0.002;
    let result = model
        .step(
            &state,
            0.0,
            dt,
            ZERO,
            None,
            |_, _| Ok(controls()),
            |_, _| Ok(ZERO),
        )
        .unwrap();
    assert!(
        (physical_momentum(&assembly, &result.motion.state)
            - physical_momentum(&assembly, &state)
            - v(result.motion.external_angular_impulse_world))
        .norm()
            < 1e-9
    );
    assert!(v(result.motion.linear_momentum_residual).norm() < 1e-10);
    let midpoint = &result.midpoint;
    let torque = v(midpoint.loads.carrier.moment_nm)
        + v(axes[0]).scale(-midpoint.coupled.torque_nm[0])
        + v(axes[1]).scale(midpoint.coupled.torque_nm[1]);
    let xmid = (v(state.rigid.pos_m) + v(result.motion.state.rigid.pos_m)).scale(0.5);
    let expected = (rotation(&state).rotate(torque).unwrap()
        + rotation(&result.motion.state).rotate(torque).unwrap())
    .scale(0.5 * dt)
        + xmid.cross(v(midpoint.loads.carrier.force_n)).scale(dt);
    assert!((expected - v(result.motion.external_angular_impulse_world)).norm() < 1e-10);
    assert!(
        midpoint
            .loads
            .carrier
            .moment_nm
            .iter()
            .all(|m| m.abs() > 0.01)
    );
    // Change only the coordinate sign of the right mechanical axis. The same
    // physical blade, spin and torque must give the same six-DOF trajectory.
    let reversed = [axes[0], a(v(axes[1]).scale(-1.0))];
    let other_assembly = self::assembly(design.gross_mass_kg, reversed, [1.0, 1.0]);
    let mut other_mounts = mounts(axes);
    other_mounts[1].handedness = PropellerHandedness::Positive;
    let other = GearedPropulsion::new_spatial(&design, &other_assembly, other_mounts).unwrap();
    let other_state = self::state(&other_assembly);
    let same = other
        .step(
            &other_state,
            0.0,
            dt,
            ZERO,
            None,
            |_, _| Ok(controls()),
            |_, _| Ok(ZERO),
        )
        .unwrap();
    assert_eq!(result.motion.state, same.motion.state);
    assert_eq!(
        result.motion.external_angular_impulse_world,
        same.motion.external_angular_impulse_world
    );
    assert_eq!(result.midpoint.coupled, same.midpoint.coupled);
    assert_eq!(
        result.midpoint.loads.aerodynamic_nm[1],
        -same.midpoint.loads.aerodynamic_nm[1]
    );
}

#[test]
fn g4_invalid_mount_and_late_stage_refusal_return_no_motion() {
    let mut design = wright_openloop_v1();
    design.coupling = precise();
    let axes = [unit([1.0, 0.08, 0.12]), unit([1.0, -0.04, -0.06])];
    let assembly = assembly(design.gross_mass_kg, axes, [1.0, -1.0]);
    let mut wrong = mounts(axes);
    wrong[0].thrust_direction_body = [1.0, 0.0, 0.0];
    assert!(GearedPropulsion::new_spatial(&design, &assembly, wrong).is_err());
    let model = GearedPropulsion::new_spatial(&design, &assembly, mounts(axes)).unwrap();
    let original = state(&assembly);
    let state = original;
    let result = model.step(
        &state,
        0.0,
        0.002,
        ZERO,
        None,
        |time, _| {
            if time > 0.0 {
                Err(fs_flyer::Refusal {
                    code: "cancelled",
                    message: "late stop".into(),
                    ranked_repairs: vec![],
                })
            } else {
                Ok(controls())
            }
        },
        |_, _| Ok(ZERO),
    );
    assert_eq!(result.unwrap_err().code, "cancelled");
    assert_eq!(state, original);
}
