//! Actual native rotor streamtubes, accepted wing forces and carrier motion.
use super::*;
use fs_flyer::propcoupling::swirl::wake::WakeTransport;

fn material<'p, 'a>(base: &'p GearedPropulsion<'a, 2>) -> WakePropulsion<'p, 'a, 2> {
    WakePropulsion::new(base, Controls::default(), WakeControls::default())
        .with_transport(WakeTransport::MaterialStreamtubes)
}

#[test]
fn native_annuli_preserve_mass_torque_and_material_surface_continuity() {
    let (d, panels, strips, properties, groups) = fixture();
    let surfaces = LiftingSurfaces::new(&panels, &strips, &properties, &groups).unwrap();
    let a = assembly(d.gross_mass_kg, AXES, [1.0; 2]);
    let base = GearedPropulsion::new_spatial(&d, &a, mounts(AXES, [1.0; 2]))
        .unwrap()
        .with_lifting_surfaces(&surfaces)
        .unwrap();
    let result = material(&base)
        .evaluate(&state(&a, 12.0, 0.0, 35.0), inputs(), &|_| Ok(ZERO), None)
        .unwrap();
    let wake = &result.coupled.wakes[0];
    let native = &result.coupled.solution.rotors[0];
    // Independent native power identity supplies V only for this nonzero-thrust fixture.
    let incident = native.power.thrust_w / native.thrust_n;
    let rho = inputs().density_kg_m3;
    let disk_flux = wake.cross_section(0.0).unwrap().mass_flux_kg_s;
    let mut powers = Vec::new();
    for s in [0.2, 0.8, 1.8] {
        let section = wake.material_annuli(s).unwrap();
        let total = wake.cross_section(s).unwrap();
        near(total.mass_flux_kg_s, disk_flux, 2e-13);
        near(total.angular_momentum_flux_nm, native.torque_nm, 2e-13);
        assert!(total.absolute_cell_power_defect_w <= wake.maximum_power_defect_w() + 1e-12);
        assert!(section.last().unwrap().outer_radius_m < d.rotor.radius_m);
        let mut radial_power = 0.0;
        for (index, cell) in section.iter().enumerate() {
            let source = wake.annuli()[index];
            let expected_mass = rho
                * (incident + native.stations[index].induction.axial_mps)
                * std::f64::consts::PI
                * (source.outer_radius_m.powi(2) - source.inner_radius_m.powi(2));
            near(cell.mass_flux_kg_s, expected_mass, 3e-13);
            let mut mass = 0.0;
            let mut angular = 0.0;
            let mut swirl_power = 0.0;
            let mut radial_integral = 0.0;
            let nr = 96;
            let nt = 12;
            let dr = (cell.outer_radius_m - cell.inner_radius_m) / nr as f64;
            for ir in 0..nr {
                let r = cell.inner_radius_m + (ir as f64 + 0.5) * dr;
                for it in 0..nt {
                    let theta = std::f64::consts::TAU * (it as f64 + 0.5) / nt as f64;
                    let er = Vec3::new(0.0, theta.cos(), theta.sin());
                    let et = Vec3::new(0.0, -theta.sin(), theta.cos());
                    let point =
                        array(vec(d.disks[0].center_m) + Vec3::new(-s, 0.0, 0.0) + er.scale(r));
                    let relative = vec(wake.relative_velocity_at(point).unwrap());
                    let u = incident + relative.x;
                    let vr = -relative.dot(er);
                    let vt = -relative.dot(et);
                    let dm = rho * u * r * dr * std::f64::consts::TAU / nt as f64;
                    mass += dm;
                    angular += dm * r * vt;
                    swirl_power += 0.5 * dm * vt * vt;
                    radial_integral += 0.5 * dm * vr * vr;
                }
            }
            near(mass, cell.mass_flux_kg_s, 3e-11);
            near(angular, source.torque_nm, 2e-11);
            near(swirl_power, cell.resolved_swirl_power_w, 2e-6);
            near(radial_integral, cell.radial_kinetic_power_w, 2e-6);
            radial_power += cell.radial_kinetic_power_w;
            // Finite-difference geometry and the physical fluid velocity follow
            // the same material boundary; interface tangential jumps do not mix.
            let h = 1e-5;
            let before = wake.material_annuli(s - h).unwrap();
            let after = wake.material_annuli(s + h).unwrap();
            let boundary_slope =
                (after[index].outer_radius_m - before[index].outer_radius_m) / (2.0 * h);
            near(boundary_slope, cell.radius_derivatives[1], 2e-9);
            let r = cell.outer_radius_m * (1.0 - 1e-10);
            let relative = wake
                .relative_velocity_at([-s, d.disks[0].center_m[1] + r, 0.0])
                .unwrap();
            near(-relative[1] / cell.axial_speed_mps, boundary_slope, 2e-9);
        }
        assert!(radial_power > 0.0);
        powers.push(radial_power);
    }
    // Central differences of the source field independently test cylindrical
    // divergence: dU/ds + (1/r)*d(r*v_r)/dr = 0. It converges at order two.
    let s = 0.8;
    let cell = wake.material_annuli(s).unwrap()[2];
    let r = 0.5 * (cell.inner_radius_m + cell.outer_radius_m);
    let velocity = |s: f64, r: f64| {
        let v = wake
            .relative_velocity_at([-s, d.disks[0].center_m[1] + r, 0.0])
            .unwrap();
        [incident + v[0], -v[1]]
    };
    let mut defects = Vec::new();
    for h in [0.004, 0.002, 0.001] {
        let axial = (velocity(s + h, r)[0] - velocity(s - h, r)[0]) / (2.0 * h);
        let radial =
            ((r + h) * velocity(s, r + h)[1] - (r - h) * velocity(s, r - h)[1]) / (2.0 * h * r);
        defects.push((axial + radial).abs());
    }
    assert!(defects[2] < 1e-5, "{defects:?}");
    for pair in defects.windows(2) {
        assert!((3.8..4.2).contains(&(pair[0] / pair[1])), "{defects:?}");
    }
}

#[test]
fn full_streamtube_field_replaces_axial_wash_and_drives_actual_midpoint_momentum() {
    let (d, mut panels, strips, properties, groups) = fixture();
    let a = assembly(d.gross_mass_kg, AXES, [1.0; 2]);
    let surfaces = LiftingSurfaces::new(&panels, &strips, &properties, &groups).unwrap();
    let base = GearedPropulsion::new_spatial(&d, &a, mounts(AXES, [1.0; 2]))
        .unwrap()
        .with_lifting_surfaces(&surfaces)
        .unwrap();
    let model = material(&base);
    let start = state(&a, 12.0, 0.0, 35.0);
    let result = model
        .step(
            &start,
            0.0,
            0.001,
            [3.0, -1.0, 2.0],
            None,
            |_, _| Ok(inputs()),
            |_, _| Ok(ZERO),
        )
        .unwrap();
    let stage = &result.initial.coupled;
    for k in 0..2 {
        let point = array((vec(panels[k].a) + vec(panels[k].b)).scale(0.5));
        let full = vec(stage.wakes[0].relative_velocity_at(point).unwrap())
            + vec(stage.wakes[1].relative_velocity_at(point).unwrap());
        let actual = vec(stage.solution.coupled.strip_velocity_mps.as_ref().unwrap()[k]);
        assert!(
            (actual - (Vec3::new(12.0, 0.0, 0.0) + full)).norm()
                <= WakeControls::default().velocity_tolerance_mps * 1.01
        );
        assert!(full.x > 0.0 && full.y > 0.0 && full.z < 0.0);
    }
    assert!(result.midpoint.aerodynamic_force_body_n[2] > 0.1);
    let dp = (vec(result.motion.state.rigid.vel_mps) - vec(start.rigid.vel_mps)).scale(a.mass_kg());
    assert!((dp - vec(result.motion.external_linear_impulse_world)).norm() < 2e-10);
    let dl = angular_momentum(&a, &result.motion.state, [3.0, -1.0, 2.0])
        - angular_momentum(&a, &start, [3.0, -1.0, 2.0]);
    assert!((dl - vec(result.motion.external_angular_impulse_world)).norm() < 2e-9);

    // Entirely upstream surfaces retain the original axial-coupling arithmetic.
    let mut upstream_panels = panels.clone();
    for panel in &mut upstream_panels {
        panel.a[0] = 0.8;
        panel.b[0] = 0.8;
        panel.ctrl[0] = 0.7;
    }
    let upstream_surfaces =
        LiftingSurfaces::new(&upstream_panels, &strips, &properties, &groups).unwrap();
    let upstream_base = GearedPropulsion::new_spatial(&d, &a, mounts(AXES, [1.0; 2]))
        .unwrap()
        .with_lifting_surfaces(&upstream_surfaces)
        .unwrap();
    let retained = material(&upstream_base)
        .evaluate(&start, inputs(), &|_| Ok(ZERO), None)
        .unwrap();
    let axial = SwirlPropulsion::new(&upstream_base, Controls::default())
        .evaluate(&start, inputs(), &|_| Ok(ZERO), None)
        .unwrap();
    assert_eq!(
        retained.aerodynamic_force_body_n,
        axial.aerodynamic_force_body_n
    );
    assert_eq!(
        retained.coupled.solution.coupled.strip_velocity_mps,
        axial.coupled.coupled.strip_velocity_mps
    );

    // A point inside the old cylinder but outside the material wake receives
    // neither its old axial increment nor a retained angular force.
    for (k, panel) in panels.iter_mut().enumerate() {
        let y = d.disks[k].center_m[1] + 0.9999 * d.rotor.radius_m;
        panel.a[1] = y - 0.15;
        panel.b[1] = y + 0.15;
        panel.ctrl[1] = y;
    }
    let outer_surfaces = LiftingSurfaces::new(&panels, &strips, &properties, &groups).unwrap();
    let outer_base = GearedPropulsion::new_spatial(&d, &a, mounts(AXES, [1.0; 2]))
        .unwrap()
        .with_lifting_surfaces(&outer_surfaces)
        .unwrap();
    let contracted = material(&outer_base)
        .evaluate(&start, inputs(), &|_| Ok(ZERO), None)
        .unwrap();
    let cylinder = WakePropulsion::new(&outer_base, Controls::default(), WakeControls::default())
        .evaluate(&start, inputs(), &|_| Ok(ZERO), None)
        .unwrap();
    assert!(cylinder.aerodynamic_force_body_n[2] > 0.1);
    assert!(contracted.aerodynamic_force_body_n[2].abs() < 1e-6);
    for velocity in contracted
        .coupled
        .solution
        .coupled
        .strip_velocity_mps
        .unwrap()
    {
        near(velocity[0], 12.0, 1.1e-7);
    }
}

#[test]
fn windmilling_expands_material_wake_and_retains_signed_torque_and_refusals() {
    let (mut d, panels, strips, properties, groups) = fixture();
    d.rotor = windmill_rotor();
    let surfaces = LiftingSurfaces::new(&panels, &strips, &properties, &groups).unwrap();
    let a = assembly(d.gross_mass_kg, AXES, [1.0, -1.0]);
    let base = GearedPropulsion::new_spatial(&d, &a, mounts(AXES, [1.0, -1.0]))
        .unwrap()
        .with_lifting_surfaces(&surfaces)
        .unwrap();
    let start = state(&a, 10.0, 0.0, 60.0);
    let model = material(&base);
    let result = model
        .step(
            &start,
            0.0,
            0.001,
            ZERO,
            None,
            |_, _| {
                Ok(PropulsionInputs {
                    motor_nm: 0.0,
                    canard_rad: 0.0,
                    ..inputs()
                })
            },
            |_, _| Ok(ZERO),
        )
        .unwrap();
    assert!(result.motion.speed_rad_s > a.speed(&start).unwrap());
    for (k, wake) in result.initial.coupled.wakes.iter().enumerate() {
        let disk = wake.cross_section(0.0).unwrap();
        let far = wake.cross_section(1.8).unwrap();
        near(disk.mass_flux_kg_s, far.mass_flux_kg_s, 1e-13);
        assert!(
            wake.material_annuli(1.8)
                .unwrap()
                .last()
                .unwrap()
                .outer_radius_m
                > 1.0
        );
        assert!(far.angular_momentum_flux_nm * [1.0, -1.0][k] < 0.0);
        let flow = wake
            .relative_velocity_at([-0.8, d.disks[k].center_m[1] + 0.7, 0.0])
            .unwrap();
        assert!(flow[0] < 0.0 && flow[1] < 0.0 && flow[2] * [1.0, -1.0][k] > 0.0);
        assert!(far.absolute_cell_power_defect_w <= wake.maximum_power_defect_w() + 1e-12);
    }
    // The mechanical transaction returns no partial motion on a late-stage refusal.
    let original = start;
    let mut calls = 0;
    let failure = model.step(
        &start,
        0.0,
        0.001,
        ZERO,
        None,
        |_, _| {
            calls += 1;
            Ok(PropulsionInputs {
                density_kg_m3: if calls == 2 { f64::NAN } else { 1.225 },
                ..inputs()
            })
        },
        |_, _| Ok(ZERO),
    );
    assert!(failure.is_err() && calls == 2);
    assert_eq!(start, original);

    // Coplanar disk footprints are disjoint, but the admitted expanded wake
    // can reach the other receiver. Its centerline sample cannot hide that.
    d.disks[1].center_m = [0.0, d.disks[0].center_m[1] + 2.05, 0.0];
    let close_base = GearedPropulsion::new_spatial(&d, &a, mounts(AXES, [1.0, -1.0]))
        .unwrap()
        .with_lifting_surfaces(&surfaces)
        .unwrap();
    let failure = material(&close_base)
        .evaluate(&start, inputs(), &|_| Ok(ZERO), None)
        .unwrap_err();
    assert_eq!(failure.code, "swirl-wake-disk-interaction");
    assert!(failure.message.contains("expanded material streamtube"));
}
