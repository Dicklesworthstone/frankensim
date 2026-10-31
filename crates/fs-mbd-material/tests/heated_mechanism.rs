//! Material-defined link dynamics driven by actual inter-link heat transport.
#[path = "heating/common.rs"]
mod common;
use common::*;
use fs_ga::{Se3, So3, Vec3, Wrench};
use fs_matdb::MaterialCard;
use fs_material::phase::{ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd::articulated::{BaseState, JointModel, forward_dynamics, forward_kinematics};
use fs_mbd_material::articulated::{MaterialLink, ResolvedArticulated, resolve_articulated};
use fs_mbd_material::thermal::JointMaterialAttachment;
use fs_mbd_material::thermal::heating::{HeatedArticulated, HeatingError, LinkHeating};
use fs_mbd_material::{RegionBinding, VolumeMoments};
use fs_thermal::enthalpy::*;
use fs_thermal::observers::ObserverCouplingError;
use std::collections::BTreeMap;

fn mechanism(a: &MaterialCard, b: &MaterialCard, temperatures: [f64; 2]) -> ResolvedArticulated {
    let pa = point(temperatures[0], 1e5);
    let pb = point(temperatures[1], 1e5);
    let cube = VolumeMoments::cuboid(Vec3::new(1., 1., 1.)).unwrap();
    let bind = |card, p, x| RegionBinding { key: "a", material: card, point: p,
        geometry: VolumeMoments::new(1., Vec3::new(x, 0., 0.), cube.inertia_m5()).unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly };
    let mount = [bind(a, &pa, 0.)];
    let rotor = [bind(b, &pb, 1.)];
    resolve_articulated(&[
        MaterialLink { name: "mount", parent: None, parent_from_child_zero: Se3::identity(),
            joint: JointModel::FIXED, regions: &mount },
        MaterialLink { name: "rotor", parent: Some(0), parent_from_child_zero: Se3::identity(),
            joint: JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(), regions: &rotor },
    ], 2, || false).unwrap()
}

fn bindings<'a>(a: &'a MaterialCard, b: &'a MaterialCard,
    ca: &'a ResolvedHeatCapacityCurve, cb: &'a ResolvedHeatCapacityCurve,
    watts: [f64; 2],
) -> Vec<LinkHeating<'a>> {
    vec![
        LinkHeating { link: "mount", regions: vec![region(a, ca, "a", 0, watts[0])] },
        LinkHeating { link: "rotor", regions: vec![region(b, cb, "a", 1, watts[1])] },
    ]
}
fn contact() -> TimedEvent {
    TimedEvent { time_s: 0., event: Event::Contact(ThermalLink {
        id: 77, a: CellId(0), b: CellId(1), area_m2: 1., conductance_w_k: 20.,
    }) }
}
fn alpha(body: &HeatedArticulated<'_>, torque: f64, gravity: f64) -> f64 {
    forward_dynamics(body.configuration().model(), BaseState::stationary(Se3::identity()),
        &[0.], &[0.], &[torque], Vec3::new(0., -gravity, 0.), &[Wrench::default(); 2])
        .unwrap().generalized_acceleration[0]
}

#[test]
fn cross_link_heat_uses_one_timeline_and_reconfigures_every_link() {
    let a = card(10., 1e-3, 500.);
    let b = card(20., 2e-3, 500.);
    let reference = mechanism(&a, &b, [300.; 2]);
    let ca = caloric(&a, 1e5, SolidLiquidPhase::Solid);
    let cb = caloric(&b, 1e5, SolidLiquidPhase::Solid);
    let mut input = bindings(&a, &b, &ca, &cb, [40., 0.]);
    input.reverse(); // local region key "a" may appear on both links
    let body = HeatedArticulated::new(&reference, input, vec![contact()], budget(), config(), || Ok(())).unwrap();
    let cfg = StepConfig { max_dt_s: 10., energy_tolerance_j: 1e-9 };
    let (plain, expected) = body.thermal().advance_to(1.3, cfg, || Ok(())).unwrap();
    let (next, actual) = body.advance_to(1.3, cfg, || Ok(())).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(next.thermal().cells(), plain.cells());
    assert!(actual.substeps > 1);
    close(actual.external_heat_j, 52.);
    close(next.thermal().energy_j(), 52.);
    assert!(next.thermal().temperature_k(CellId(1)).unwrap() > 300.);
    for (i, coefficient) in [1e-3, 2e-3].into_iter().enumerate() {
        let region = &next.configuration().bodies()[i].regions()[0];
        let t = next.thermal().temperature_k(CellId(i as u64)).unwrap();
        close(region.expansion().free_linear_strain(), coefficient * (t - 300.));
        assert_eq!(region.inertia().mass().to_bits(), reference.bodies()[i].inertia().mass().to_bits());
        assert_eq!(next.thermal().cells()[&CellId(i as u64)].mass_kg.to_bits(), region.inertia().mass().to_bits());
    }
    // Uniform density cube with COM one metre from the z hinge: Izz = 7/3.
    let s = 1. + 2e-3 * (next.thermal().temperature_k(CellId(1)).unwrap() - 300.);
    close(alpha(&next, 1., 0.), 3. / (7. * s * s));
    close(alpha(&next, 0., 9.81), -6. * 9.81 / (7. * s));
    let base = BaseState::stationary(Se3::identity());
    assert_eq!(forward_kinematics(reference.model(), base, &[0.3], &[0.2]).unwrap(),
        forward_kinematics(next.configuration().model(), base, &[0.3], &[0.2]).unwrap());
    assert_eq!(body.thermal().time_s(), 0.);
}

#[test]
fn heat_capacity_only_swap_changes_contact_history_and_actual_mechanical_response() {
    let mount = card(10., 1e-3, 500.);
    let rotor_a = card(10., 2e-3, 500.);
    let rotor_b = card(20., 2e-3, 500.);
    let cm = caloric(&mount, 1e5, SolidLiquidPhase::Solid);
    let ca = caloric(&rotor_a, 1e5, SolidLiquidPhase::Solid);
    let cb = caloric(&rotor_b, 1e5, SolidLiquidPhase::Solid);
    let ra = mechanism(&mount, &rotor_a, [300.; 2]);
    let rb = mechanism(&mount, &rotor_b, [300.; 2]);
    assert_eq!(ra.model(), rb.model()); // no imposed mass, torque or geometry change
    let make = |r, c, caloric| HeatedArticulated::new(r,
        bindings(&mount, c, &cm, caloric, [40., 0.]), vec![contact()], budget(), config(), || Ok(())).unwrap();
    let (a, aa) = make(&ra, &rotor_a, &ca).advance_to(2., config(), || Ok(())).unwrap();
    let (b, ab) = make(&rb, &rotor_b, &cb).advance_to(2., config(), || Ok(())).unwrap();
    close(aa.external_heat_j, ab.external_heat_j);
    assert!(a.thermal().temperature_k(CellId(1)).unwrap() > b.thermal().temperature_k(CellId(1)).unwrap());
    assert!(alpha(&a, 1., 0.) < alpha(&b, 1., 0.));
    assert_eq!(a.thermal().mass_kg(), b.thermal().mass_kg());
}

#[test]
fn different_reference_temperatures_equilibrate_without_replacing_matter() {
    let a = card(10., 1e-3, 500.);
    let b = card(20., 2e-3, 500.);
    let reference = mechanism(&a, &b, [350., 300.]);
    let ca = caloric(&a, 1e5, SolidLiquidPhase::Solid);
    let cb = caloric(&b, 1e5, SolidLiquidPhase::Solid);
    let mut input = bindings(&a, &b, &ca, &cb, [0.; 2]);
    input[0].regions[0].initial_specific_enthalpy_j_kg = 500.;
    let old = HeatedArticulated::new(&reference, input, vec![contact()], budget(), config(), || Ok(())).unwrap();
    let (next, audit) = old.advance_to(12., config(), || Ok(())).unwrap();
    let equilibrium = (20. * 350. + 40. * 300.) / 60.;
    for id in [CellId(0), CellId(1)] {
        assert!((next.thermal().temperature_k(id).unwrap() - equilibrium).abs() < 1e-7);
    }
    close(audit.external_heat_j, 0.);
    close(next.thermal().energy_j(), old.thermal().energy_j());
    close(next.thermal().mass_kg(), 4.);
    assert!(next.configuration().bodies()[0].regions()[0].geometry().volume_m3() < 1.);
    assert!(next.configuration().bodies()[1].regions()[0].geometry().volume_m3() > 1.);
}

#[derive(Clone)]
struct Heat { calls: usize, work: f64, fail: bool }
impl DepositionObserver for Heat {
    type Error = DepositionError;
    fn sources(&mut self, _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        self.calls += 1;
        Ok(BTreeMap::from([(CellId(0), 30.), (CellId(1), 10.)]))
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> { Ok(()) }
    fn step(&mut self, _: f64, dt: f64, _: &BTreeMap<CellId, f64>, _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.work += 40. * dt;
        if self.fail { return Err(DepositionError::InvalidInput("generated heat refuses")); }
        Ok(())
    }
}

#[test]
fn all_link_sources_and_late_failures_share_one_transaction() {
    let a = card(10., 1e-3, 500.);
    let reference = mechanism(&a, &a, [300.; 2]);
    let c = caloric(&a, 1e5, SolidLiquidPhase::Solid);
    let body = HeatedArticulated::new(&reference, bindings(&a, &a, &c, &c, [5., 0.]),
        vec![contact()], budget(), config(), || Ok(())).unwrap();
    let source = Heat { calls: 0, work: 0., fail: false };
    let (next, accepted, audit) = body.advance_with(&source, 1.3, config(), || Ok(())).unwrap();
    close(audit.external_heat_j, 45. * 1.3);
    close(next.thermal().energy_j(), 45. * 1.3);
    close(accepted.work, 40. * 1.3);
    assert_eq!(accepted.calls, audit.substeps);
    let bad = Heat { fail: true, ..source.clone() };
    assert!(matches!(body.advance_with(&bad, 1., config(), || Ok(())), Err(ObserverCouplingError::First(_))));
    let mut calls = 0;
    body.advance_with(&source, 1., config(), || { calls += 1; Ok(()) }).unwrap();
    let mut counter = 0;
    assert!(body.advance_with(&source, 1., config(), || {
        counter += 1;
        if counter == calls { Err(DepositionError::Cancelled) } else { Ok(()) }
    }).is_err());
    assert_eq!(body.thermal().time_s(), 0.);
    assert_eq!(body.thermal().energy_j(), 0.);
    assert_eq!(source.calls, 0);
    assert_eq!(source.work, 0.);
    assert_eq!(body.configuration().model(), reference.model());
}

#[test]
fn intermediate_material_refusal_rolls_back_other_links_too() {
    let a = card(10., 1e-3, 500.);
    let b = card(10., 1e-3, 301.);
    let reference = mechanism(&a, &b, [300.; 2]);
    let ca = caloric(&a, 1e5, SolidLiquidPhase::Solid);
    let cb = caloric(&b, 1e5, SolidLiquidPhase::Solid);
    let events = vec![TimedEvent { time_s: 1., event: Event::Boundary(CellId(1), boundary(-40.)) }];
    let body = HeatedArticulated::new(&reference, bindings(&a, &b, &ca, &cb, [5., 40.]),
        events, budget(), config(), || Ok(())).unwrap();
    let (plain, _) = body.thermal().advance_to(2., config(), || Ok(())).unwrap();
    close(plain.temperature_k(CellId(1)).unwrap(), 300.); // endpoint-only check would miss the excursion
    assert!(matches!(body.advance_to(2., config(), || Ok(())), Err(HeatingError::Mechanical(_))));
    assert_eq!(body.thermal().time_s(), 0.);
    assert_eq!(body.configuration().model(), reference.model());
}

#[test]
fn global_admission_and_aligned_restarts_preserve_the_complete_mechanism() {
    let a = card(10., 1e-3, 500.);
    let reference = mechanism(&a, &a, [300.; 2]);
    let c = caloric(&a, 1e5, SolidLiquidPhase::Solid);
    let input = bindings(&a, &a, &c, &c, [40., 0.]);
    let build = |input, budget| HeatedArticulated::new(&reference, input, vec![contact(),
        TimedEvent { time_s: 1., event: Event::Boundary(CellId(0), boundary(0.)) }], budget, config(), || Ok(()));
    for which in 0..5 {
        let mut invalid = input.clone();
        match which {
            0 => invalid[1].link = "mount",
            1 => invalid[1].link = "foreign",
            2 => invalid[1].regions[0].cell_id = CellId(0),
            3 => invalid[1].regions[0].material_id = MaterialId(10),
            _ => invalid[1].regions.clear(),
        }
        assert!(matches!(build(invalid, budget()), Err(HeatingError::Binding(_))));
    }
    assert!(build(input[..1].to_vec(), budget()).is_err());
    for b in [DepositionBudget { max_cells: 1, ..budget() },
        DepositionBudget { max_events: 3, ..budget() }] {
        assert!(matches!(build(input.clone(), b), Err(HeatingError::Thermal(DepositionError::Budget(_)))));
    }
    let old = build(input.clone(), budget()).unwrap();
    let mut reverse = input.clone(); reverse.reverse();
    let other = build(reverse, budget()).unwrap();
    assert_eq!(old.configuration(), other.configuration());
    assert_eq!(old.thermal().cells(), other.thermal().cells());
    let (whole, audit) = old.advance_to(3., config(), || Ok(())).unwrap();
    let (first, a1) = old.advance_to(1., config(), || Ok(())).unwrap();
    let (rest, a2) = first.advance_to(3., config(), || Ok(())).unwrap();
    assert_eq!(whole.thermal().cells(), rest.thermal().cells());
    assert_eq!(whole.configuration(), rest.configuration());
    close(audit.external_heat_j, a1.external_heat_j + a2.external_heat_j);
    let small = build(input, DepositionBudget { max_substeps: 1, ..budget() }).unwrap();
    assert!(matches!(small.advance_to(3., config(), || Ok(())), Err(HeatingError::Thermal(DepositionError::Budget(_)))));
    assert_eq!(small.thermal().time_s(), 0.);
}

#[test]
fn sourced_parent_expansion_moves_joint_geometry_and_native_lever_response() {
    for expansion_per_k in [1e-3, 2e-3] {
        let a = card(10., expansion_per_k, 500.);
        let b = card(20., 0., 500.);
        let ca = caloric(&a, 1e5, SolidLiquidPhase::Solid);
        let cb = caloric(&b, 1e5, SolidLiquidPhase::Solid);
        let p = point(300., 1e5);
        let geometry = VolumeMoments::cuboid(Vec3::new(1., 1., 1.)).unwrap();
        let parent = [RegionBinding {
            key: "a",
            material: &a,
            point: &p,
            geometry,
            selection: MaterialPropertySelection::SingleClaimOnly,
        }];
        let child = [RegionBinding {
            key: "a",
            material: &b,
            point: &p,
            geometry,
            selection: MaterialPropertySelection::SingleClaimOnly,
        }];
        let hinge = JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap();
        let reference = resolve_articulated(
            &[
                MaterialLink {
                    name: "mount",
                    parent: None,
                    parent_from_child_zero: Se3::identity(),
                    joint: JointModel::FIXED,
                    regions: &parent,
                },
                MaterialLink {
                    name: "lever",
                    parent: Some(0),
                    parent_from_child_zero: Se3::identity(),
                    joint: hinge,
                    regions: &parent,
                },
                MaterialLink {
                    name: "tip",
                    parent: Some(1),
                    parent_from_child_zero: Se3::from_parts(
                        So3::identity(),
                        Vec3::new(0.5, 0., 0.),
                    )
                    .unwrap(),
                    joint: hinge,
                    regions: &child,
                },
            ],
            3,
            || false,
        )
        .unwrap();
        let original = HeatedArticulated::new(
            &reference,
            vec![
                LinkHeating {
                    link: "mount",
                    regions: vec![region(&a, &ca, "a", 0, 0.)],
                },
                LinkHeating {
                    link: "lever",
                    regions: vec![region(&a, &ca, "a", 1, 40.)],
                },
                LinkHeating {
                    link: "tip",
                    regions: vec![region(&b, &cb, "a", 2, 0.)],
                },
            ],
            vec![],
            budget(),
            config(),
            || Ok(()),
        )
        .unwrap();
        let attachments = ["lever", "tip"].map(|name| JointMaterialAttachment {
            link: name.into(),
            parent_region: "a".into(),
            child_region: "a".into(),
        });
        let old = original
            .with_joint_attachments(&attachments, || Ok(()))
            .unwrap();
        assert_eq!(old.configuration().model(), reference.model());
        let reverse = [attachments[1].clone(), attachments[0].clone()];
        assert_eq!(
            old.configuration(),
            original
                .with_joint_attachments(&reverse, || Ok(()))
                .unwrap()
                .configuration()
        );
        let (next, audit) = old.advance_to(2., config(), || Ok(())).unwrap();
        // Only the parent's source alpha changes between these two material cases.
        // rho=2, V=1, Cp=10: 80 J heats the parent by 4 K.
        let stretch = 1. + 4. * expansion_per_k;
        close(audit.external_heat_j, 80.);
        let native = next.configuration().model();
        close(
            native.links()[2].parent_from_child_zero().translation().x,
            0.5 * stretch,
        );
        let kinematics = forward_kinematics(
            native,
            BaseState::stationary(Se3::identity()),
            &[0., 0.],
            &[0., 0.],
        )
        .unwrap();
        close(kinematics.world_from_link[2].translation().x, 0.5 * stretch);
        // H11-H12 = I_parent + m_child L^2 = (1/3 + 1/2) s^2.
        let response = forward_dynamics(
            native,
            BaseState::stationary(Se3::identity()),
            &[0., 0.],
            &[0., 0.],
            &[1., 0.],
            Vec3::default(),
            &[Wrench::default(); 3],
        )
        .unwrap();
        close(
            response.generalized_acceleration[0],
            6. / (5. * stretch * stretch),
        );
        close(
            response.generalized_acceleration[1],
            -response.generalized_acceleration[0],
        );
        for (body, reference) in next.configuration().bodies().iter().zip(reference.bodies()) {
            assert_eq!(
                body.inertia().mass().to_bits(),
                reference.inertia().mass().to_bits()
            );
        }
        let (first, _) = old.advance_to(1., config(), || Ok(())).unwrap();
        let (rest, _) = first.advance_to(2., config(), || Ok(())).unwrap();
        assert_eq!(next.configuration(), rest.configuration());
        // Rebinding an already hot model still starts at ORIGINAL material geometry.
        assert_eq!(
            next.configuration(),
            next.with_joint_attachments(&attachments, || Ok(()))
                .unwrap()
                .configuration()
        );
        assert_eq!(original.configuration().model(), reference.model());
    }
}

#[test]
fn material_joint_assignments_refuse_unknown_roles_and_unsupported_child_motion_atomically() {
    let a = card(10., 1e-3, 500.);
    let c = caloric(&a, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&a, &a, [300.; 2]);
    let body = HeatedArticulated::new(
        &reference,
        bindings(&a, &a, &c, &c, [40., 0.]),
        vec![],
        budget(),
        config(),
        || Ok(()),
    )
    .unwrap();
    let valid = JointMaterialAttachment {
        link: "rotor".into(),
        parent_region: "a".into(),
        child_region: "a".into(),
    };
    for which in 0..4 {
        let mut bad = valid.clone();
        match which {
            0 => bad.link = "mount".into(),
            1 => bad.link = "unknown".into(),
            2 => bad.parent_region = "unknown".into(),
            _ => bad.child_region = "unknown".into(),
        }
        assert!(body.with_joint_attachments(&[bad], || Ok(())).is_err());
    }
    assert!(
        body.with_joint_attachments(&[valid.clone(), valid.clone()], || Ok(()))
            .is_err()
    );
    let mut input = bindings(&a, &a, &c, &c, [40., 0.]);
    input[1].regions[0].fixed_point_m = Vec3::new(0.1, 0., 0.);
    let moving_origin =
        HeatedArticulated::new(&reference, input, vec![], budget(), config(), || Ok(())).unwrap();
    assert!(
        moving_origin
            .with_joint_attachments(&[valid.clone()], || Ok(()))
            .is_err()
    );
    let mut count = 0;
    body.with_joint_attachments(&[valid.clone()], || {
        count += 1;
        Ok(())
    })
    .unwrap();
    let mut calls = 0;
    let refused = body.with_joint_attachments(&[valid], || {
        calls += 1;
        if calls == count {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        refused,
        Err(HeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert!(body.configuration().joint_attachments().is_empty());
    assert_eq!(body.configuration().model(), reference.model());
    assert_eq!(body.thermal().time_s(), 0.);
}
