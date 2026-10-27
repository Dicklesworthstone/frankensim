//! Material-interface data -> native joint traction -> both finite bodies -> next friction.
#[path = "heating/common.rs"]
mod common;
use common::*;
use fs_evidence::ValidityDomain;
use fs_ga::{Se3, Vec3};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::phase::{ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_mbd::articulated::JointModel;
use fs_mbd_material::articulated::{MaterialLink, ResolvedArticulated, resolve_articulated};
use fs_mbd_material::thermal::heating::articulated_motion::{
    Controls, Error, HeatedArticulatedMotion,
    dry::{
        DryBrakeLimits, DryJointBrake, JointSlipGeometry, KINETIC_FRICTION_PROPERTY,
        STATIC_FRICTION_PROPERTY,
    },
};
use fs_mbd_material::thermal::heating::{HeatedArticulated, LinkHeating};
use fs_mbd_material::{RegionBinding, VolumeMoments};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_thermal::enthalpy::*;
use fs_tribo::{HeatPartition, InputAuthority};

fn interface(
    parent: &MaterialCard,
    child: &MaterialCard,
    mu: [f64; 2],
    upper: f64,
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for (name, values) in [
        (STATIC_FRICTION_PROPERTY, [1.0, 1.0]),
        (KINETIC_FRICTION_PROPERTY, mu),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, Dims::NONE),
                value: PropertyValue::Curve {
                    abscissa: "TA".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (upper, values[1])],
                    dims: Dims::NONE,
                },
                validity: ValidityDomain::unconstrained()
                    .with("TA", 300.0, upper)
                    .with("TB", 300.0, 500.0)
                    .with("pressure", 0.0, 1e6),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic brake fixture, not measured friction".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: child.id().clone(),
            texture_frame: "child surface".into(),
        },
        SurfaceSpec {
            material: parent.id().clone(),
            texture_frame: "parent surface".into(),
        },
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "synthetic fixture".into(),
            history: "fixed surface history".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("TA", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("TB", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("pressure", QuantitySpec::dimensional(Pressure::DIMS), 1e5)
        .unwrap()
}
fn mechanism(cards: &[&MaterialCard], joint: JointModel) -> ResolvedArticulated {
    let p = point(300.0, 1e5);
    let regions: Vec<_> = cards
        .iter()
        .map(|card| {
            vec![RegionBinding {
                key: "a",
                material: card,
                point: &p,
                geometry: VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap(),
                selection: MaterialPropertySelection::SingleClaimOnly,
            }]
        })
        .collect();
    let names = ["mount", "first", "second"];
    let links: Vec<_> = regions
        .iter()
        .enumerate()
        .map(|(i, regions)| MaterialLink {
            name: names[i],
            parent: if i == 0 { None } else { Some(i - 1) },
            parent_from_child_zero: Se3::identity(),
            joint: if i == 0 { JointModel::FIXED } else { joint },
            regions,
        })
        .collect();
    resolve_articulated(&links, links.len(), || false).unwrap()
}
fn moving<'a>(
    reference: &'a ResolvedArticulated,
    cards: &[&'a MaterialCard],
    heat: &[&'a ResolvedHeatCapacityCurve],
    speed: Vec<f64>,
    dt: f64,
) -> HeatedArticulatedMotion<'a> {
    let bindings = reference
        .model()
        .links()
        .iter()
        .enumerate()
        .map(|(i, link)| LinkHeating {
            link: link.name(),
            regions: vec![region(cards[i], heat[i], "a", i as u64, 0.0)],
        })
        .collect();
    let heated =
        HeatedArticulated::new(reference, bindings, vec![], budget(), config(), || Ok(())).unwrap();
    let n = speed.len();
    HeatedArticulatedMotion::new(
        heated,
        Se3::identity(),
        vec![0.0; n],
        speed,
        vec![0.0; n],
        Vec3::default(),
        Controls {
            maximum_step_s: dt,
            relative_state_tolerance: 1e-11,
            absolute_position_tolerance: vec![1e-12; n],
            absolute_velocity_tolerance: vec![1e-12; n],
            maximum_rk_attempts_per_substep: 256,
            relative_momentum_tolerance: 1e-11,
            absolute_momentum_tolerance: vec![1e-12; n],
            maximum_parameter_work_j: 0.0,
            maximum_numerical_energy_change_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
}
fn brake<'a>(
    card: &'a InterfaceSystemCard,
    parent: &'a MaterialCard,
    child: &'a MaterialCard,
    context: &'a QueryPoint,
    link: usize,
    geometry: JointSlipGeometry,
) -> DryJointBrake<'a> {
    DryJointBrake {
        link,
        cells: [CellId(link as u64), CellId(link as u64 - 1)],
        geometry,
        normal_force_n: 2.0,
        nominal_area_m2: 1e-5,
        pressure_axis: "pressure",
        point: ThermalInterfacePoint {
            interface: card,
            materials: [child, parent],
            context,
            temperature_axes: ["TA", "TB"],
            liquid_fraction_axes: [None, None],
            mean_temperature_axis: None,
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        authority: InputAuthority::SyntheticFixture,
        partition: HeatPartition::new(0.75, 0.25, 0.0).unwrap(),
        limits: DryBrakeLimits {
            regularization_m_s: 0.01,
            maximum_coefficient_change: 0.2,
            absolute_loss_tolerance_j: 1e-12,
        },
    }
}
fn hinge() -> JointModel {
    JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap()
}
fn rotor() -> JointSlipGeometry {
    JointSlipGeometry::Revolute {
        effective_radius_m: 0.2,
    }
}

#[test]
fn source_friction_swap_slows_real_motion_and_heats_both_ordered_bodies() {
    let parent = card(10.0, 0.0, 500.0);
    let child = card(5.0, 0.0, 500.0);
    let hp = caloric(&parent, 1e5, SolidLiquidPhase::Solid);
    let hc = caloric(&child, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&[&parent, &child], hinge());
    let point = context();
    for mu in [0.25, 0.5] {
        let source = interface(&parent, &child, [mu; 2], 500.0);
        for sign in [-1.0, 1.0] {
            let initial = moving(
                &reference,
                &[&parent, &child],
                &[&hp, &hc],
                vec![sign * 3.0],
                0.125,
            )
            .with_dry_joint_brakes(
                vec![brake(&source, &parent, &child, &point, 1, rotor())],
                || Ok(()),
            )
            .unwrap();
            let (next, audit) = initial.advance_to(1.0, config(), || Ok(())).unwrap();
            // Independent cube I = rho/6 = 1/3; |tau| = mu*N*r = 0.4 mu.
            let acceleration = 1.2 * mu;
            let speed = 3.0 - acceleration;
            let loss = (9.0 - speed * speed) / 6.0;
            close(next.joint_velocity()[0], sign * speed);
            close(next.joint_position()[0], sign * (3.0 - acceleration / 2.0));
            close(audit.friction_heat_j, loss);
            close(audit.viscous_heat_j, 0.0);
            close(audit.thermal.external_heat_j, loss);
            close(
                next.thermal().temperature_k(CellId(1)).unwrap(),
                300.0 + 0.75 * loss / 10.0,
            );
            close(
                next.thermal().temperature_k(CellId(0)).unwrap(),
                300.0 + 0.25 * loss / 20.0,
            );
            let state = next.dry_joint_brakes().next().unwrap();
            close(state.work.surface_a_heat_j(), 0.75 * loss);
            close(state.work.surface_b_heat_j(), 0.25 * loss);
            close(
                state.resolved.state_coordinate("TA").unwrap(),
                state.phases[0].temperature_k(),
            );
            close(
                state.resolved.state_coordinate("TB").unwrap(),
                state.phases[1].temperature_k(),
            );
            close(audit.combined_residual_j, 0.0);
            close(initial.thermal().time_s(), 0.0);
            close(
                initial
                    .dry_joint_brakes()
                    .next()
                    .unwrap()
                    .work
                    .dissipated_work_j(),
                0.0,
            );
        }
    }
}

#[test]
fn internal_brake_preserves_coupled_rotor_momentum_and_only_heats_its_bodies() {
    let solid = card(10.0, 0.0, 500.0);
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&[&solid; 3], hinge());
    let source = interface(&solid, &solid, [0.25; 2], 500.0);
    let point = context();
    let initial = moving(
        &reference,
        &[&solid; 3],
        &[&heat; 3],
        vec![3.0, -1.0],
        0.125,
    )
    .with_dry_joint_brakes(
        vec![brake(&source, &solid, &solid, &point, 2, rotor())],
        || Ok(()),
    )
    .unwrap();
    let (next, audit) = initial.advance_to(0.5, config(), || Ok(())).unwrap();
    // M = [[2/3,1/3],[1/3,1/3]], tau = [0,+0.1].
    // The first conjugate momentum must stay 5/3 despite both rates changing.
    close(next.joint_velocity()[0], 2.85);
    close(next.joint_velocity()[1], -0.7);
    close(
        (2.0 * next.joint_velocity()[0] + next.joint_velocity()[1]) / 3.0,
        5.0 / 3.0,
    );
    close(next.joint_position()[0], 1.4625);
    close(next.joint_position()[1], -0.425);
    let loss = 0.1 * 0.425;
    close(audit.friction_heat_j, loss);
    close(next.thermal().temperature_k(CellId(0)).unwrap(), 300.0);
    close(
        next.thermal().temperature_k(CellId(2)).unwrap(),
        300.0 + 0.75 * loss / 20.0,
    );
    close(
        next.thermal().temperature_k(CellId(1)).unwrap(),
        300.0 + 0.25 * loss / 20.0,
    );
    close(
        audit.mechanical_before.kinetic_j - audit.mechanical_after.kinetic_j,
        loss,
    );
}

#[test]
fn same_material_brake_handles_translation_and_native_low_speed_regularization() {
    let solid = card(10.0, 0.0, 500.0);
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let source = interface(&solid, &solid, [0.5; 2], 500.0);
    let point = context();
    let slider = mechanism(
        &[&solid; 2],
        JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
    );
    let initial = moving(&slider, &[&solid; 2], &[&heat; 2], vec![1.0], 0.125)
        .with_dry_joint_brakes(
            vec![brake(
                &source,
                &solid,
                &solid,
                &point,
                1,
                JointSlipGeometry::Prismatic,
            )],
            || Ok(()),
        )
        .unwrap();
    let (next, audit) = initial.advance_to(0.5, config(), || Ok(())).unwrap();
    // m = 2, F = mu*N = 1, v = 1-t/2, x = t-t²/4.
    close(next.joint_velocity()[0], 0.75);
    close(next.joint_position()[0], 0.4375);
    close(audit.friction_heat_j, 1.0 - 0.75 * 0.75);

    let reference = mechanism(&[&solid; 2], hinge());
    for speed in [-0.01, 0.0, 0.01] {
        let initial = moving(
            &reference,
            &[&solid; 2],
            &[&heat; 2],
            vec![speed],
            1.0 / 64.0,
        )
        .with_dry_joint_brakes(
            vec![brake(&source, &solid, &solid, &point, 1, rotor())],
            || Ok(()),
        )
        .unwrap();
        let (next, audit) = initial.advance_to(0.125, config(), || Ok(())).unwrap();
        // Below v_reg=.01: c=mu*N*r²/v_reg=4, c/I=12 s^-1.
        let decay = (-1.5_f64).exp();
        close(next.joint_velocity()[0], speed * decay);
        close(next.joint_position()[0], speed * (1.0 - decay) / 12.0);
        close(
            audit.friction_heat_j,
            speed * speed * (1.0 - decay * decay) / 6.0,
        );
        assert!(next.dry_joint_brakes().next().unwrap().dissipated_power_w >= 0.0);
    }
}

#[test]
fn self_heating_friction_fade_refines_to_independent_coupled_solution_and_restarts() {
    let parent = card(10.0, 0.0, 500.0);
    let child = card(0.1, 0.0, 500.0);
    let hp = caloric(&parent, 1e5, SolidLiquidPhase::Solid);
    let hc = caloric(&child, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&[&parent, &child], hinge());
    // mu(T_A) = 0.5 [1 - 0.05(T_A-300)].
    let source = interface(&parent, &child, [0.5, 0.25], 310.0);
    let point = context();
    // Eliminate heat: C_A=.2, I=1/3, child heat fraction=.75.
    // mu = .5(A+B*w²), B=.05*.75*I/(2*C_A)=.03125, A=1-9B.
    let b: f64 = 0.03125;
    let a: f64 = 1.0 - 9.0 * b;
    let k: f64 = 0.6; // mu0*N*r/I
    let target: f64 = 2.0;
    let angle = (3.0 * (b / a).sqrt()).atan() - k * (a * b).sqrt() * target;
    let exact = (a / b).sqrt() * angle.tan();
    let mut errors = Vec::new();
    for dt in [0.25, 0.125, 0.0625] {
        let initial = moving(&reference, &[&parent, &child], &[&hp, &hc], vec![3.0], dt)
            .with_dry_joint_brakes(
                vec![brake(&source, &parent, &child, &point, 1, rotor())],
                || Ok(()),
            )
            .unwrap();
        let (next, audit) = initial.advance_to(target, config(), || Ok(())).unwrap();
        errors.push((next.joint_velocity()[0] - exact).abs());
        let state = next.dry_joint_brakes().next().unwrap();
        let initial_state = initial.dry_joint_brakes().next().unwrap();
        let temperature = next.thermal().temperature_k(CellId(1)).unwrap();
        close(
            state.law.zero_speed_kinetic_coefficient().unwrap(),
            0.5 * (1.0 - 0.05 * (temperature - 300.0)),
        );
        for (new, old) in state
            .resolved
            .properties()
            .iter()
            .zip(initial_state.resolved.properties())
        {
            assert_eq!(new.answer().receipt.selected, old.answer().receipt.selected);
        }
        close(
            audit.mechanical_before.kinetic_j - audit.mechanical_after.kinetic_j,
            audit.friction_heat_j,
        );
        close(audit.combined_residual_j, 0.0);
        let (prefix, _) = initial.advance_to(1.0, config(), || Ok(())).unwrap();
        let (resumed, _) = prefix.advance_to(target, config(), || Ok(())).unwrap();
        assert_eq!(resumed.joint_position(), next.joint_position());
        assert_eq!(resumed.joint_velocity(), next.joint_velocity());
        assert_eq!(
            resumed.dry_joint_brakes().next(),
            next.dry_joint_brakes().next()
        );
    }
    assert!(
        errors[0] / errors[1] > 1.7 && errors[0] / errors[1] < 2.3,
        "{errors:?}"
    );
    assert!(
        errors[1] / errors[2] > 1.7 && errors[1] / errors[2] < 2.3,
        "{errors:?}"
    );
}

#[test]
fn source_refusal_wrong_body_mapping_and_late_cancellation_preserve_every_owner() {
    let solid = card(0.1, 0.0, 500.0);
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&[&solid; 2], hinge());
    let point = context();
    let source = interface(&solid, &solid, [0.5; 2], 500.0);
    let base = moving(&reference, &[&solid; 2], &[&heat; 2], vec![3.0], 0.125);
    // A value measured at one point is not an admitted constant friction law.
    let mut exact_claims = ClaimSet::new();
    for (_, original) in source.claims().claims_ordered() {
        let mut claim = original.clone();
        claim.value = PropertyValue::Scalar {
            value: if claim.key.name() == STATIC_FRICTION_PROPERTY {
                1.0
            } else {
                0.5
            },
            dims: Dims::NONE,
        };
        claim.interpolation = InterpolationPolicy::TabulatedOnly;
        claim.validity = ValidityDomain::unconstrained();
        exact_claims.insert_claim(claim).unwrap();
    }
    let exact = InterfaceSystemCard::assemble(
        source.surface_a().clone(),
        source.surface_b().clone(),
        source.context().clone(),
        exact_claims,
        vec![],
    )
    .unwrap();
    assert!(matches!(
        base.with_dry_joint_brakes(
            vec![brake(&exact, &solid, &solid, &point, 1, rotor())],
            || Ok(()),
        ),
        Err(Error::Invalid(_))
    ));
    let mut wrong = brake(&source, &solid, &solid, &point, 1, rotor());
    wrong.cells.swap(0, 1);
    assert!(base.with_dry_joint_brakes(vec![wrong], || Ok(())).is_err());
    let initial = base
        .with_dry_joint_brakes(
            vec![brake(&source, &solid, &solid, &point, 1, rotor())],
            || Ok(()),
        )
        .unwrap();
    let state = initial.dry_joint_brakes().next().unwrap().clone();
    let mut calls = 0usize;
    initial
        .advance_to(0.5, config(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    assert!(calls > 4);
    let mut attempted = 0usize;
    let cancelled = initial.advance_to(0.5, config(), || {
        attempted += 1;
        if attempted == calls - 1 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(cancelled, Err(Error::Heating(_))));
    assert_eq!(initial.dry_joint_brakes().next().unwrap(), &state);
    assert_eq!(initial.joint_velocity(), &[3.0]);
    assert_eq!(initial.joint_position(), &[0.0]);
    close(initial.thermal().time_s(), 0.0);
    close(initial.thermal().temperature_k(CellId(0)).unwrap(), 300.0);
    close(initial.thermal().temperature_k(CellId(1)).unwrap(), 300.0);

    let bounded = interface(&solid, &solid, [0.5; 2], 300.1);
    let narrow = base
        .with_dry_joint_brakes(
            vec![brake(&bounded, &solid, &solid, &point, 1, rotor())],
            || Ok(()),
        )
        .unwrap();
    let before = narrow.dry_joint_brakes().next().unwrap().clone();
    assert!(matches!(
        narrow.advance_to(0.5, config(), || Ok(())),
        Err(Error::Material(_))
    ));
    assert_eq!(narrow.dry_joint_brakes().next().unwrap(), &before);
    close(narrow.thermal().time_s(), 0.0);
    close(narrow.thermal().temperature_k(CellId(0)).unwrap(), 300.0);
    close(narrow.thermal().temperature_k(CellId(1)).unwrap(), 300.0);
}
