//! Independent physical loads for explicit swept, canted and vertical sections.
use fs_flyer::Refusal;
use fs_flyer::aircraft::{OpenLoopDesign, wright_openloop_v1};
use fs_flyer::propcoupling::{CANDIDATE_B, PropCouplingSolverSpec};
use fs_flyer::spine::SixDofState;
use fs_flyer::spine::rotors::geared::{GearedAssembly, GearedRotor};
use fs_flyer::spine::rotors::propulsion::surfaces::{
    InducedDragGroup, LiftingSurfaces, SectionProperties,
};
use fs_flyer::spine::rotors::propulsion::{
    GearedPropulsion, PropellerHandedness, PropellerMount, PropulsionEvaluation, PropulsionInputs,
};
use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_wing::nonlinear::StripSpec;
use fs_wing::{Panel, SurfaceId, flat_surface};

const RHO: f64 = 1.225;
const ZERO: [f64; 3] = [0.0; 3];

fn v(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}
fn a(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
fn rotation(axis_angle: [f64; 3]) -> So3 {
    So3::exp(So3Tangent::new(v(axis_angle))).unwrap()
}
fn rotate(r: So3, value: [f64; 3]) -> [f64; 3] {
    a(r.rotate(v(value)).unwrap())
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "{actual:e} != {expected:e}; tolerance {tolerance:e}"
    );
}
fn near_vector(actual: Vec3, expected: Vec3, tolerance: f64) {
    assert!(
        (actual - expected).norm() < tolerance,
        "{actual:?} != {expected:?}; tolerance {tolerance:e}"
    );
}
fn refused<T>(result: Result<T, Refusal>, code: &str) {
    match result {
        Ok(_) => panic!("expected {code}"),
        Err(error) => assert_eq!(error.code, code),
    }
}
fn inputs() -> PropulsionInputs {
    PropulsionInputs {
        density_kg_m3: RHO,
        canard_rad: 0.0,
        motor_nm: 0.0,
        gravity_world_mps2: ZERO,
    }
}

#[derive(Clone)]
struct Scene {
    design: OpenLoopDesign,
    panels: Vec<Panel>,
    strips: Vec<StripSpec>,
    properties: Vec<SectionProperties>,
    groups: Vec<InducedDragGroup>,
    // Analytic fixture axes, independent of SectionFrame::from_panel/flow.
    chord: Vec<[f64; 3]>,
    span: Vec<[f64; 3]>,
    axes: [[f64; 3]; 2],
    velocity: [f64; 3],
    omega: [f64; 3],
}

impl Scene {
    fn swept(cant: f64) -> Self {
        let mut design = wright_openloop_v1();
        design.coupling = PropCouplingSolverSpec {
            tol: 1e-10,
            cap: 64,
            ..CANDIDATE_B
        };
        design.drag_z_m = design.cg_m[2]; // Parasite line at CG rotates objectively.
        let sweep: f64 = 0.25;
        let cant = rotation([cant, 0.0, 0.0]);
        let mut panels = flat_surface(SurfaceId::WingLower, 6.0, 1.2, 0.0, 0.0, 6, 2).unwrap();
        let place = |mut point: [f64; 3]| {
            point[0] += sweep.tan() * point[1];
            a(v(rotate(cant, point)) + Vec3::new(0.2, 0.0, -0.8))
        };
        for panel in &mut panels {
            panel.a = place(panel.a);
            panel.b = place(panel.b);
            panel.ctrl = place(panel.ctrl);
            panel.normal = rotate(cant, panel.normal);
        }
        let strips = (0..6)
            .map(|s| StripSpec {
                panel_indices: vec![s, s + 6],
                chord_m: 1.2 * sweep.cos(),
                twist_rad: 0.01 + s as f64 * 0.001,
            })
            .collect();
        Self {
            design,
            panels,
            strips,
            properties: (0..6)
                .map(|s| SectionProperties {
                    camber_ratio: 0.012 + 0.002 * s as f64,
                    attached_cd: 0.01 + 0.002 * s as f64,
                    canard_gain: 0.0,
                    induced_drag_group: 0,
                })
                .collect(),
            groups: vec![InducedDragGroup {
                span_m: 6.0 / sweep.cos(),
                oswald_e: 0.82,
            }],
            chord: vec![rotate(cant, [sweep.cos(), -sweep.sin(), 0.0]); 6],
            span: vec![rotate(cant, [sweep.sin(), sweep.cos(), 0.0]); 6],
            axes: [[1.0, 0.0, 0.0]; 2],
            velocity: [12.0, 2.3, 0.7],
            omega: [0.015, -0.01, 0.02],
        }
    }

    fn opposing() -> Self {
        let mut scene = Self::swept(0.0);
        scene.panels.clear();
        scene.strips.clear();
        scene.properties.clear();
        scene.chord.clear();
        scene.span.clear();
        scene.design.cg_m = [-0.27, 0.0, 0.0];
        scene.design.drag_z_m = 0.0;
        for disk in &mut scene.design.disks {
            disk.center_m[2] = 0.0;
        }
        for (g, id) in [SurfaceId::WingUpper, SurfaceId::WingLower]
            .into_iter()
            .enumerate()
        {
            let start = scene.panels.len();
            scene
                .panels
                .extend(flat_surface(id, 4.0, 0.9, 0.0, [-1.5, 1.5][g], 4, 2).unwrap());
            for s in 0..4 {
                scene.strips.push(StripSpec {
                    panel_indices: vec![start + s, start + s + 4],
                    chord_m: 0.9,
                    twist_rad: [0.12, -0.12][g],
                });
                scene.properties.push(SectionProperties {
                    camber_ratio: 0.0,
                    attached_cd: 0.01,
                    canard_gain: 0.0,
                    induced_drag_group: g,
                });
                scene.chord.push([1.0, 0.0, 0.0]);
                scene.span.push([0.0, 1.0, 0.0]);
            }
        }
        scene.groups = vec![
            InducedDragGroup {
                span_m: 4.0,
                oswald_e: 0.8,
            };
            2
        ];
        scene.velocity = [12.0, 0.0, 0.0];
        scene.omega = ZERO;
        scene
    }

    fn rotated(&self, r: So3) -> Self {
        let mut moved = self.clone();
        for panel in &mut moved.panels {
            panel.a = rotate(r, panel.a);
            panel.b = rotate(r, panel.b);
            panel.ctrl = rotate(r, panel.ctrl);
            panel.normal = rotate(r, panel.normal);
        }
        moved.design.cg_m = rotate(r, self.design.cg_m);
        moved.design.drag_z_m = moved.design.cg_m[2];
        for disk in &mut moved.design.disks {
            disk.center_m = rotate(r, disk.center_m);
        }
        moved.chord = self.chord.iter().map(|&axis| rotate(r, axis)).collect();
        moved.span = self.span.iter().map(|&axis| rotate(r, axis)).collect();
        moved.axes = self.axes.map(|axis| rotate(r, axis));
        moved.velocity = rotate(r, self.velocity);
        moved.omega = rotate(r, self.omega);
        moved
    }

    fn assembly(&self) -> GearedAssembly<2> {
        GearedAssembly::new(
            self.design.gross_mass_kg,
            Mat3 {
                m: [1000.0, 0.0, 0.0, 0.0, 1000.0, 0.0, 0.0, 0.0, 1000.0],
            },
            core::array::from_fn(|k| GearedRotor {
                axis_body: self.axes[k],
                ratio: [1.0, -1.0][k],
                axial_inertia_kgm2: [1.5, 2.5][k],
            }),
        )
        .unwrap()
    }

    fn mounts(&self) -> [PropellerMount; 2] {
        core::array::from_fn(|k| PropellerMount {
            rotor_index: k,
            thrust_direction_body: self.axes[k],
            handedness: [PropellerHandedness::Positive, PropellerHandedness::Negative][k],
        })
    }

    fn solve(&self) -> PropulsionEvaluation<2> {
        let assembly = self.assembly();
        let surfaces =
            LiftingSurfaces::new(&self.panels, &self.strips, &self.properties, &self.groups)
                .unwrap();
        let model = GearedPropulsion::new_spatial(&self.design, &assembly, self.mounts())
            .unwrap()
            .with_lifting_surfaces(&surfaces)
            .unwrap();
        let state = assembly
            .state_from_speed(
                SixDofState {
                    pos_m: ZERO,
                    vel_mps: self.velocity,
                    quat: [1.0, 0.0, 0.0, 0.0],
                    omega_body: self.omega,
                },
                0.0,
                35.0,
            )
            .unwrap();
        model
            .evaluate(&state, inputs(), &|_| Ok(ZERO), None)
            .unwrap()
    }
}

struct Oracle {
    force: Vec3,
    moment: Vec3,
    group_force: Vec<Vec3>,
    group_drag: Vec<f64>,
    profile_force: Vec3,
    pitch_moment: Vec3,
}

// Independent panel quadrature and dimensional section/drag equations. No
// production SectionFlow, CirculatoryWrenchAccumulator or load helper is used.
fn oracle(scene: &Scene, result: &PropulsionEvaluation<2>) -> Oracle {
    let mut force = Vec3::default();
    let mut moment = Vec3::default();
    let mut groups = vec![Vec3::default(); scene.groups.len()];
    let mut profile = Vec3::default();
    let mut pitch = Vec3::default();
    let local = result.coupled.strip_velocity_mps.as_ref().unwrap();
    for (s, strip) in scene.strips.iter().enumerate() {
        assert!(result.coupled.strip_alpha_eff_rad[s].abs() < 0.3);
        let incident = v(local[s]);
        let section = scene.properties[s];
        for &j in &strip.panel_indices {
            let panel = scene.panels[j];
            let bound = v(panel.b) - v(panel.a);
            let panel_force = incident.cross(bound).scale(RHO * result.coupled.gamma[j]);
            let arm = (v(panel.a) + v(panel.b)).scale(0.5) - v(scene.design.cg_m);
            near(panel_force.dot(incident), 0.0, 2e-9);
            force = force + panel_force;
            moment = moment + arm.cross(panel_force);
            groups[section.induced_drag_group] = groups[section.induced_drag_group] + panel_force;
        }
        let chord = v(scene.chord[s]);
        let span = v(scene.span[s]);
        let normal = chord.cross(span);
        let uc = incident.dot(chord);
        let un = incident.dot(normal);
        let speed2 = uc * uc + un * un;
        let speed = speed2.sqrt();
        let projected = incident - span.scale(incident.dot(span));
        let q = 0.5 * RHO * speed2;
        let first = scene.panels[strip.panel_indices[0]];
        let area = strip.chord_m * (v(first.b) - v(first.a)).norm();
        let drag = q * area * section.attached_cd;
        let profile_force = projected.scale(-drag / speed);
        near(profile_force.dot(incident), -drag * speed, 2e-9);
        let couple =
            span.scale(-core::f64::consts::PI * section.camber_ratio * q * area * strip.chord_m);
        let arm = (v(first.a) + v(first.b)).scale(0.5) - v(scene.design.cg_m);
        force = force + profile_force;
        moment = moment + arm.cross(profile_force) + couple;
        profile = profile + profile_force;
        pitch = pitch + couple;
    }
    let speed2 = v(scene.velocity).dot(v(scene.velocity));
    let direction = v(scene.velocity).scale(1.0 / speed2.sqrt());
    let q = 0.5 * RHO * speed2;
    let group_drag: Vec<f64> = scene
        .groups
        .iter()
        .zip(&groups)
        .map(|(group, load)| {
            let transverse = *load - direction.scale(load.dot(direction));
            let drag = transverse.dot(transverse)
                / (q * core::f64::consts::PI * group.oswald_e * group.span_m.powi(2));
            force = force - direction.scale(drag);
            drag
        })
        .collect();
    let parasite = q
        * scene
            .design
            .ledger
            .components
            .iter()
            .map(|c| c.area_m2 * c.cd)
            .sum::<f64>()
        * (1.0 + scene.design.ledger.unresolved_interference_frac);
    let parasite_force = direction.scale(-parasite);
    force = force + parasite_force;
    moment = moment
        + Vec3::new(0.0, 0.0, scene.design.drag_z_m - scene.design.cg_m[2]).cross(parasite_force);
    for k in 0..2 {
        let thrust = v(scene.axes[k]).scale(result.coupled.thrust_n[k]);
        force = force + thrust;
        moment = moment + (v(scene.design.disks[k].center_m) - v(scene.design.cg_m)).cross(thrust);
        // Shaft torques remain separate external ports, never duplicated in
        // the carrier-only moment assembled above.
        near(
            result.loads.aerodynamic_nm[k],
            [-1.0, 1.0][k] * result.coupled.torque_nm[k],
            1e-12,
        );
    }
    Oracle {
        force,
        moment,
        group_force: groups,
        group_drag,
        profile_force: profile,
        pitch_moment: pitch,
    }
}

#[test]
fn actual_swept_canted_and_vertical_loads_match_independent_full_wrench() {
    for cant in [0.65, core::f64::consts::FRAC_PI_2] {
        let scene = Scene::swept(cant);
        let result = scene.solve();
        let expected = oracle(&scene, &result);
        near_vector(v(result.aerodynamic_force_body_n), expected.force, 3e-8);
        near_vector(v(result.loads.carrier.force_n), expected.force, 3e-8);
        near_vector(v(result.loads.carrier.moment_nm), expected.moment, 6e-8);
        assert!(expected.profile_force.norm() > 5.0);
        assert!(expected.pitch_moment.x.abs() > 1.0);
        assert!(expected.pitch_moment.z.abs() > 1.0);
        let panel = scene.panels[0];
        assert!((v(panel.b) - v(panel.a)).norm() > (panel.b[1] - panel.a[1]).abs() + 0.1);
        if cant == core::f64::consts::FRAC_PI_2 {
            assert!((panel.b[1] - panel.a[1]).abs() < 1e-12);
            assert!(
                expected.group_force[0].y < -10.0,
                "vertical sideforce must oppose positive lateral flow"
            );
        }
    }
}

#[test]
fn separate_drag_groups_preserve_drag_of_opposing_lifting_surfaces() {
    let scene = Scene::opposing();
    let separated = scene.solve();
    let physical = oracle(&scene, &separated);
    assert!(physical.group_force[0].z * physical.group_force[1].z < 0.0);
    assert!((physical.group_force[0] + physical.group_force[1]).norm() < 1e-8);
    assert!(physical.group_drag.iter().sum::<f64>() > 1.0);
    near_vector(v(separated.loads.carrier.moment_nm), physical.moment, 6e-8);
    near_vector(v(separated.aerodynamic_force_body_n), physical.force, 3e-8);

    let mut merged = scene.clone();
    merged.groups.truncate(1);
    for section in &mut merged.properties {
        section.induced_drag_group = 0;
    }
    let cancellation = merged.solve();
    let merged_oracle = oracle(&merged, &cancellation);
    assert_eq!(cancellation.coupled.gamma, separated.coupled.gamma);
    assert_eq!(cancellation.coupled.thrust_n, separated.coupled.thrust_n);
    assert!(merged_oracle.group_drag[0] < 1e-15);
    let drag = physical.group_drag.iter().sum::<f64>();
    near_vector(
        v(separated.aerodynamic_force_body_n) - v(cancellation.aerodynamic_force_body_n),
        Vec3::new(-drag, 0.0, 0.0),
        3e-8,
    );
    // Group drag acts at CG; changing only group membership adds no moment.
    assert_eq!(
        separated.loads.carrier.moment_nm,
        cancellation.loads.carrier.moment_nm
    );
}

#[test]
fn native_custom_scene_rotates_objectively_including_negative_global_forward_axis() {
    let scene = Scene::swept(0.55);
    let original = scene.solve();
    for r in [
        rotation([0.35, -0.25, 0.4]),
        rotation([0.0, core::f64::consts::PI, 0.0]),
    ] {
        let moved = scene.rotated(r);
        let result = moved.solve();
        let independent = oracle(&moved, &result);
        near_vector(v(result.aerodynamic_force_body_n), independent.force, 3e-8);
        near_vector(v(result.loads.carrier.moment_nm), independent.moment, 6e-8);
        near_vector(
            v(result.aerodynamic_force_body_n),
            r.rotate(v(original.aerodynamic_force_body_n)).unwrap(),
            8e-7,
        );
        near_vector(
            v(result.loads.carrier.moment_nm),
            r.rotate(v(original.loads.carrier.moment_nm)).unwrap(),
            1e-6,
        );
        for (g, reference) in result.coupled.gamma.iter().zip(&original.coupled.gamma) {
            near(*g, *reference, 4e-9);
        }
        for k in 0..2 {
            near(
                result.coupled.thrust_n[k],
                original.coupled.thrust_n[k],
                1e-7,
            );
            near(
                result.loads.aerodynamic_nm[k],
                original.loads.aerodynamic_nm[k],
                1e-8,
            );
        }
        if moved.velocity[0] < 0.0 {
            assert!(moved.velocity[0] < -10.0);
            for chord in &moved.chord {
                assert!(v(moved.velocity).dot(v(*chord)) > 10.0);
            }
            for axis in &moved.axes {
                assert!(v(moved.velocity).dot(v(*axis)) > 10.0);
            }
        }
    }
}

#[test]
fn custom_layout_rejects_incomplete_or_incompatible_physical_ownership() {
    let scene = Scene::swept(0.6);
    let admit = |panels: &[Panel],
                 strips: &[StripSpec],
                 properties: &[SectionProperties],
                 groups: &[InducedDragGroup]| {
        LiftingSurfaces::new(panels, strips, properties, groups).map(|_| ())
    };
    admit(
        &scene.panels,
        &scene.strips,
        &scene.properties,
        &scene.groups,
    )
    .unwrap();
    refused(
        admit(
            &scene.panels,
            &scene.strips,
            &scene.properties[..5],
            &scene.groups,
        ),
        "lifting-surfaces-invalid",
    );
    let mut groups = scene.groups.clone();
    groups.push(groups[0]);
    refused(
        admit(&scene.panels, &scene.strips, &scene.properties, &groups),
        "lifting-surfaces-invalid",
    );
    groups[0].span_m = 0.0;
    refused(
        admit(
            &scene.panels,
            &scene.strips,
            &scene.properties,
            &groups[..1],
        ),
        "lifting-surfaces-invalid",
    );
    let mut strips = scene.strips.clone();
    strips[1].panel_indices[0] = strips[0].panel_indices[0];
    refused(
        admit(&scene.panels, &strips, &scene.properties, &scene.groups),
        "lifting-surfaces-invalid",
    );
    let mut properties = scene.properties.clone();
    properties[0].canard_gain = f64::NAN;
    refused(
        admit(&scene.panels, &scene.strips, &properties, &scene.groups),
        "lifting-surfaces-invalid",
    );
    for shift in [
        v(scene.span[0]).scale(0.1),
        v(scene.chord[0]).cross(v(scene.span[0])).scale(0.1),
    ] {
        let mut panels = scene.panels.clone();
        let row = &mut panels[scene.strips[0].panel_indices[1]];
        row.a = a(v(row.a) + shift);
        row.b = a(v(row.b) + shift);
        row.ctrl = a(v(row.ctrl) + shift);
        refused(
            admit(&panels, &scene.strips, &scene.properties, &scene.groups),
            "lifting-surfaces-invalid",
        );
    }
    let mut panels = scene.panels.clone();
    panels[scene.strips[0].panel_indices[1]].surface = SurfaceId::WingUpper;
    refused(
        admit(&panels, &scene.strips, &scene.properties, &scene.groups),
        "lifting-surfaces-invalid",
    );
    let mut panels = scene.panels.clone();
    let row = &mut panels[scene.strips[0].panel_indices[1]];
    core::mem::swap(&mut row.a, &mut row.b);
    refused(
        admit(&panels, &scene.strips, &scene.properties, &scene.groups),
        "section-frame-invalid",
    );

    let assembly = scene.assembly();
    let surfaces = LiftingSurfaces::new(
        &scene.panels,
        &scene.strips,
        &scene.properties,
        &scene.groups,
    )
    .unwrap();
    let legacy = GearedPropulsion::new(&scene.design, &assembly, scene.mounts()).unwrap();
    refused(
        legacy.with_lifting_surfaces(&surfaces),
        "lifting-surfaces-invalid",
    );
}
