//! G1/G3: coupled source heating changes a continuing joint trajectory.
#[path = "heating/common.rs"]
mod common;
use common::*;
use fs_ga::{Se3, So3, Vec3};
use fs_matdb::MaterialCard;
use fs_material::phase::{ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd::articulated::JointModel;
use fs_mbd_material::articulated::{MaterialLink, ResolvedArticulated, resolve_articulated};
use fs_mbd_material::thermal::JointMaterialAttachment;
use fs_mbd_material::thermal::heating::articulated_motion::{
    Controls, Error, HeatedArticulatedMotion,
};
use fs_mbd_material::thermal::heating::{HeatedArticulated, LinkHeating};
use fs_mbd_material::{RegionBinding, VolumeMoments};
use fs_thermal::enthalpy::*;
use fs_thermal::observers::ObserverCouplingError;
use std::collections::BTreeMap;

fn mechanism(a: &MaterialCard, b: &MaterialCard) -> ResolvedArticulated {
    mechanism_at(a, b, Vec3::default())
}
fn mechanism_at(a: &MaterialCard, b: &MaterialCard, second_origin: Vec3) -> ResolvedArticulated {
    let p = point(300., 1e5);
    let bind = |material| {
        [RegionBinding {
            key: "a",
            material,
            point: &p,
            geometry: VolumeMoments::cuboid(Vec3::new(1., 1., 1.)).unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }]
    };
    let mount = bind(a);
    let first = bind(a);
    let second = bind(b);
    let hinge = JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap();
    resolve_articulated(
        &[
            MaterialLink {
                name: "mount",
                parent: None,
                parent_from_child_zero: Se3::identity(),
                joint: JointModel::FIXED,
                regions: &mount,
            },
            MaterialLink {
                name: "first",
                parent: Some(0),
                parent_from_child_zero: Se3::identity(),
                joint: hinge,
                regions: &first,
            },
            MaterialLink {
                name: "second",
                parent: Some(1),
                parent_from_child_zero: Se3::from_parts(So3::identity(), second_origin).unwrap(),
                joint: hinge,
                regions: &second,
            },
        ],
        3,
        || false,
    )
    .unwrap()
}

#[test]
fn expanding_joint_lever_uses_full_momentum_transfer_and_parameter_work() {
    let a = card(10., 1e-3, 500.);
    let b = card(10., 0., 500.);
    let ca = caloric(&a, 1e5, SolidLiquidPhase::Solid);
    let cb = caloric(&b, 1e5, SolidLiquidPhase::Solid);
    // The horizontal spinning case changes the full joint mass matrix. The
    // vertical resting case isolates the gravitational work of moving anchors.
    for (offset, velocity, gravity) in [
        (Vec3::new(0.5, 0., 0.), vec![3., -1.], Vec3::default()),
        (
            Vec3::new(0., 0.5, 0.),
            vec![0.; 2],
            Vec3::new(0., -9.81, 0.),
        ),
    ] {
        let reference = mechanism_at(&a, &b, offset);
        let thermal = HeatedArticulated::new(
            &reference,
            vec![
                LinkHeating {
                    link: "mount",
                    regions: vec![region(&a, &ca, "a", 0, 0.)],
                },
                LinkHeating {
                    link: "first",
                    regions: vec![region(&a, &ca, "a", 1, 40.)],
                },
                LinkHeating {
                    link: "second",
                    regions: vec![region(&b, &cb, "a", 2, 0.)],
                },
            ],
            vec![],
            budget(),
            config(),
            || Ok(()),
        )
        .unwrap()
        .with_joint_attachments(
            &[JointMaterialAttachment {
                link: "second".into(),
                parent_region: "a".into(),
                child_region: "a".into(),
            }],
            || Ok(()),
        )
        .unwrap();
        let make = |controls| {
            HeatedArticulatedMotion::new(
                thermal.clone(),
                Se3::identity(),
                vec![0.; 2],
                velocity.clone(),
                vec![0.; 2],
                gravity,
                controls,
                || Ok(()),
            )
            .unwrap()
        };
        let old = make(controls());
        let momentum = old.joint_momentum().unwrap();
        let (next, audit) = old.advance_to(2., config(), || Ok(())).unwrap();
        let stretch = 1.004;
        let point = next.configuration().model().links()[2]
            .parent_from_child_zero()
            .translation();
        close(point.x, offset.x * stretch);
        close(point.y, offset.y * stretch);
        for (actual, expected) in next.joint_momentum().unwrap().iter().zip(&momentum) {
            close(*actual, *expected);
        }
        if offset.x > 0. {
            // I_parent+m_child*L^2=(1/3+1/2)*s^2; absolute child spin is 2 rad/s.
            close(next.joint_velocity()[0], 3. / stretch.powi(2));
            close(next.joint_velocity()[1], 2. - 3. / stretch.powi(2));
            close(audit.parameter_work_j, 3.75 * (1. / stretch.powi(2) - 1.));
        } else {
            close(next.joint_velocity()[0], 0.);
            close(next.joint_velocity()[1], 0.);
            close(audit.parameter_work_j, 2. * 9.81 * 0.5 * (stretch - 1.));
        }
        close(audit.thermal.external_heat_j, 80.);
        close(audit.applied_joint_work_j, 0.);
        close(audit.combined_residual_j, 0.);
        let (middle, _) = old.advance_to(1., config(), || Ok(())).unwrap();
        let (restarted, _) = middle.advance_to(2., config(), || Ok(())).unwrap();
        assert_eq!(next.configuration(), restarted.configuration());
        assert_eq!(next.joint_position(), restarted.joint_position());
        assert_eq!(next.joint_velocity(), restarted.joint_velocity());
        let constrained = make(Controls {
            maximum_parameter_work_j: 0.,
            ..controls()
        });
        assert!(matches!(
            constrained.advance_to(2., config(), || Ok(())),
            Err(Error::Budget("absolute thermal-configuration work"))
        ));
        assert_eq!(constrained.thermal().time_s(), 0.);
        assert_eq!(constrained.joint_velocity(), velocity.as_slice());
    }
}
fn controls() -> Controls {
    Controls {
        maximum_step_s: 0.0625,
        relative_state_tolerance: 1e-10,
        absolute_position_tolerance: vec![1e-11; 2],
        absolute_velocity_tolerance: vec![1e-11; 2],
        maximum_rk_attempts_per_substep: 128,
        relative_momentum_tolerance: 1e-11,
        absolute_momentum_tolerance: vec![1e-12; 2],
        maximum_parameter_work_j: 10.,
        maximum_numerical_energy_change_j: 1e-7,
    }
}
fn heated<'a>(
    reference: &'a ResolvedArticulated,
    a: &'a MaterialCard,
    b: &'a MaterialCard,
    ca: &'a ResolvedHeatCapacityCurve,
    cb: &'a ResolvedHeatCapacityCurve,
    watts: [f64; 2],
    events: Vec<TimedEvent>,
    controls: Controls,
    effort: [f64; 2],
) -> HeatedArticulatedMotion<'a> {
    let thermal = HeatedArticulated::new(
        reference,
        vec![
            LinkHeating {
                link: "mount",
                regions: vec![region(a, ca, "a", 0, 0.)],
            },
            LinkHeating {
                link: "first",
                regions: vec![region(a, ca, "a", 1, watts[0])],
            },
            LinkHeating {
                link: "second",
                regions: vec![region(b, cb, "a", 2, watts[1])],
            },
        ],
        events,
        budget(),
        config(),
        || Ok(()),
    )
    .unwrap();
    HeatedArticulatedMotion::new(
        thermal,
        Se3::identity(),
        vec![0.; 2],
        vec![3., -1.],
        effort.to_vec(),
        Vec3::default(),
        controls,
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn heat_capacity_swap_changes_coupled_motion_at_conserved_joint_momentum() {
    let b = card(10., 0., 500.);
    let cb = caloric(&b, 1e5, SolidLiquidPhase::Solid);
    let mut first_rates = Vec::new();
    for cp in [10., 20.] {
        let a = card(cp, 1e-3, 500.);
        let ca = caloric(&a, 1e5, SolidLiquidPhase::Solid);
        let reference = mechanism(&a, &b);
        let old = heated(
            &reference,
            &a,
            &b,
            &ca,
            &cb,
            [40., 0.],
            vec![],
            controls(),
            [0.; 2],
        );
        // Both cubes have m=2, Izz=1/3. Coaxial serial joints give
        // M=[[I1+I2,I2],[I2,I2]], not two independent rotor inertias.
        let p = old.joint_momentum().unwrap();
        close(p[0], 5. / 3.);
        close(p[1], 2. / 3.);
        let (next, audit) = old.advance_to(2., config(), || Ok(())).unwrap();
        let stretch = 1. + 1e-3 * 80. / (2. * cp);
        let first = 3. / (stretch * stretch);
        close(next.joint_velocity()[0], first);
        close(next.joint_velocity()[1], 2. - first);
        for (&actual, &expected) in next.joint_momentum().unwrap().iter().zip(&p) {
            close(actual, expected);
        }
        assert!(next.joint_position()[0] > 0.);
        assert!(next.joint_position()[1] < 0.);
        close(audit.thermal.external_heat_j, 80.);
        close(audit.applied_joint_work_j, 0.);
        close(
            audit.parameter_work_j,
            1.5 * (1. / (stretch * stretch) - 1.),
        );
        close(audit.absolute_parameter_work_j, -audit.parameter_work_j);
        close(audit.combined_residual_j, 0.);
        assert!(audit.maximum_scaled_momentum_residual <= 1.);
        assert_eq!(audit.thermal.substeps, 32);
        for body in next.configuration().bodies() {
            close(body.inertia().mass(), 2.);
        }
        assert_eq!(old.joint_position(), &[0., 0.]);
        assert_eq!(old.joint_velocity(), &[3., -1.]);
        first_rates.push(first);
    }
    assert!(first_rates[0] < first_rates[1]);
}

#[test]
fn applied_joint_work_and_no_heat_follow_the_native_coupled_mass_law() {
    let a = card(10., 1e-3, 500.);
    let c = caloric(&a, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&a, &a);
    let old = heated(
        &reference,
        &a,
        &a,
        &c,
        &c,
        [0.; 2],
        vec![],
        controls(),
        [0.2, -0.1],
    );
    let (next, audit) = old.advance_to(0.5, config(), || Ok(())).unwrap();
    // M^-1=[[3,-3],[-3,6]], alpha=[0.9,-1.2].
    close(next.joint_velocity()[0], 3.45);
    close(next.joint_velocity()[1], -1.6);
    close(next.joint_position()[0], 1.5 + 0.5 * 0.9 * 0.25);
    close(next.joint_position()[1], -0.5 - 0.5 * 1.2 * 0.25);
    close(
        audit.applied_joint_work_j,
        0.2 * next.joint_position()[0] - 0.1 * next.joint_position()[1],
    );
    close(audit.parameter_work_j, 0.);
    close(audit.numerical_energy_change_j, 0.);
    close(audit.combined_residual_j, 0.);
}

#[derive(Clone)]
struct ExtraHeat {
    work: f64,
    fail: bool,
}
impl DepositionObserver for ExtraHeat {
    type Error = DepositionError;
    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        Ok(BTreeMap::from([(CellId(2), 10.)]))
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.work += 10. * dt;
        if self.fail {
            Err(DepositionError::InvalidInput("source refusal"))
        } else {
            Ok(())
        }
    }
}

#[test]
fn shared_heat_observer_and_late_refusals_publish_no_partial_joint_state() {
    let a = card(10., 1e-3, 500.);
    let c = caloric(&a, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&a, &a);
    let old = heated(
        &reference,
        &a,
        &a,
        &c,
        &c,
        [20., 0.],
        vec![],
        controls(),
        [0.; 2],
    );
    let source = ExtraHeat {
        work: 0.,
        fail: false,
    };
    let (_, accepted, audit) = old.advance_with(&source, 1., config(), || Ok(())).unwrap();
    close(accepted.work, 10.);
    close(audit.thermal.external_heat_j, 30.);
    let failing = ExtraHeat {
        fail: true,
        ..source.clone()
    };
    assert!(matches!(
        old.advance_with(&failing, 1., config(), || Ok(())),
        Err(ObserverCouplingError::First(_))
    ));
    let mut calls = 0;
    old.advance_to(1., config(), || {
        calls += 1;
        Ok(())
    })
    .unwrap();
    let mut count = 0;
    assert!(
        old.advance_to(1., config(), || {
            count += 1;
            if count == calls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        })
        .is_err()
    );
    let mut no_work = controls();
    no_work.maximum_parameter_work_j = 0.;
    let bounded = heated(
        &reference,
        &a,
        &a,
        &c,
        &c,
        [20., 0.],
        vec![],
        no_work,
        [0.; 2],
    );
    assert!(matches!(
        bounded.advance_to(1., config(), || Ok(())),
        Err(Error::Budget(_))
    ));
    assert_eq!(old.joint_position(), &[0., 0.]);
    assert_eq!(old.joint_velocity(), &[3., -1.]);
    assert_eq!(old.thermal().time_s(), 0.);
    assert_eq!(source.work, 0.);
}

#[test]
fn cooling_and_aligned_restart_keep_one_accepted_mechanism() {
    let a = card(10., 1e-3, 500.);
    let c = caloric(&a, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&a, &a);
    let old = heated(
        &reference,
        &a,
        &a,
        &c,
        &c,
        [40., 0.],
        vec![TimedEvent {
            time_s: 1.,
            event: Event::Boundary(CellId(1), boundary(-40.)),
        }],
        controls(),
        [0.; 2],
    );
    let (whole, full) = old.advance_to(2., config(), || Ok(())).unwrap();
    let (half, _) = old.advance_to(1., config(), || Ok(())).unwrap();
    let (split, _) = half.advance_to(2., config(), || Ok(())).unwrap();
    assert_eq!(whole.joint_position(), split.joint_position());
    assert_eq!(whole.joint_velocity(), split.joint_velocity());
    assert_eq!(whole.thermal().cells(), split.thermal().cells());
    close(whole.joint_velocity()[0], 3.);
    close(whole.joint_velocity()[1], -1.);
    assert!(full.absolute_parameter_work_j > full.parameter_work_j.abs());
    close(full.combined_residual_j, 0.);
}
