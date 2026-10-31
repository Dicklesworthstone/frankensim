//! G1/G3/G4: sourced speed/temperature friction drives actual motion and both bodies' heat.
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
    Controls, HeatedArticulatedMotion,
    dry::{
        DryBrakeLimits, DryJointBrake, JointSlipGeometry, KINETIC_FRICTION_PROPERTY,
        RateDependentDryJointBrake, STATIC_FRICTION_PROPERTY,
    },
};
use fs_mbd_material::thermal::heating::{HeatedArticulated, LinkHeating};
use fs_mbd_material::{RegionBinding, VolumeMoments};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature, Velocity};
use fs_thermal::enthalpy::*;
use fs_tribo::{HeatPartition, InputAuthority};

fn interface(
    cards: [&MaterialCard; 2],
    a: f64,
    b: f64,
    thermal_slope: f64,
    upper: f64,
    rate_surface: bool,
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for name in [STATIC_FRICTION_PROPERTY, KINETIC_FRICTION_PROPERTY] {
        let value = |t: f64, v: f64| {
            if name == STATIC_FRICTION_PROPERTY {
                2.0
            } else {
                a + b * v + thermal_slope * (t - 300.0)
            }
        };
        let response = if rate_surface {
            PropertyValue::Surface2D {
                axes: ["TA".into(), "speed".into()],
                axis_dims: [Temperature::DIMS, Velocity::DIMS],
                samples: [[300.0, 0.0], [upper, 0.0], [300.0, 2.0], [upper, 2.0]]
                    .into_iter()
                    .map(|p| (p, value(p[0], p[1])))
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: Dims::NONE,
            }
        } else {
            PropertyValue::Curve {
                abscissa: "TA".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, value(300.0, 0.0)), (upper, value(upper, 0.0))],
                dims: Dims::NONE,
            }
        };
        let validity = ValidityDomain::unconstrained()
            .with("TA", 300.0, upper)
            .with("TB", 300.0, 500.0)
            .with("pressure", 0.0, 1e6);
        let validity = if rate_surface {
            validity.with_quantity("speed", QuantitySpec::dimensional(Velocity::DIMS), 0.0, 2.0)
        } else {
            validity
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, Dims::NONE),
                value: response,
                validity,
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic speed/temperature friction fixture, not measured tribology"
                        .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "child".into(),
        },
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "parent".into(),
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
        .with_quantity("speed", QuantitySpec::dimensional(Velocity::DIMS), 0.0)
        .unwrap()
}

fn mechanism(cards: [&MaterialCard; 2]) -> ResolvedArticulated {
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
    resolve_articulated(
        &[
            MaterialLink {
                name: "mount",
                parent: None,
                parent_from_child_zero: Se3::identity(),
                joint: JointModel::FIXED,
                regions: &regions[0],
            },
            MaterialLink {
                name: "slider",
                parent: Some(0),
                parent_from_child_zero: Se3::identity(),
                joint: JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
                regions: &regions[1],
            },
        ],
        2,
        || false,
    )
    .unwrap()
}

fn moving<'a>(
    reference: &'a ResolvedArticulated,
    cards: [&'a MaterialCard; 2],
    heat: [&'a ResolvedHeatCapacityCurve; 2],
    speed: f64,
    dt: f64,
    drive: f64,
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
    HeatedArticulatedMotion::new(
        heated,
        Se3::identity(),
        vec![0.0],
        vec![speed],
        vec![drive],
        Vec3::default(),
        Controls {
            maximum_step_s: dt,
            relative_state_tolerance: 1e-11,
            absolute_position_tolerance: vec![1e-12],
            absolute_velocity_tolerance: vec![1e-12],
            maximum_rk_attempts_per_substep: 256,
            relative_momentum_tolerance: 1e-11,
            absolute_momentum_tolerance: vec![1e-12],
            maximum_parameter_work_j: 0.0,
            maximum_numerical_energy_change_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
}

fn binding<'a>(
    source: &'a InterfaceSystemCard,
    cards: [&'a MaterialCard; 2],
    point: &'a QueryPoint,
    maximum: f64,
) -> RateDependentDryJointBrake<'a> {
    RateDependentDryJointBrake {
        brake: DryJointBrake {
            link: 1,
            cells: [CellId(1), CellId(0)],
            geometry: JointSlipGeometry::Prismatic,
            normal_force_n: 2.0,
            nominal_area_m2: 1e-5,
            pressure_axis: "pressure",
            point: ThermalInterfacePoint {
                interface: source,
                materials: [cards[1], cards[0]],
                context: point,
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
        },
        slip_speed_axis: "speed",
        maximum_slip_speed_m_s: maximum,
    }
}

#[test]
fn material_speed_law_changes_signed_motion_and_conjugate_heat() {
    let parent = card(10.0, 0.0, 500.0);
    let child = card(5.0, 0.0, 500.0);
    let cards = [&parent, &child];
    let hp = caloric(&parent, 1e5, SolidLiquidPhase::Solid);
    let hc = caloric(&child, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(cards);
    let point = context();
    for b in [0.25, 0.5] {
        let a = 0.2;
        let source = interface(cards, a, b, 0.0, 500.0, true);
        for sign in [-1.0, 1.0] {
            let initial = moving(&reference, cards, [&hp, &hc], sign * 1.5, 0.125, 0.0)
                .with_rate_dependent_dry_joint_brakes(
                    vec![binding(&source, cards, &point, 2.0)],
                    || Ok(()),
                )
                .unwrap();
            let (next, audit) = initial.advance_to(1.0, config(), || Ok(())).unwrap();
            // m=N=2, so v' = -(a+bv); this independent exact solution
            // distinguishes the evolving rate law from a frozen initial mu.
            let decay = (-b).exp();
            let speed = (1.5 + a / b) * decay - a / b;
            let distance = (1.5 + a / b) * (1.0 - decay) / b - a / b;
            let loss = 1.5 * 1.5 - speed * speed;
            close(next.joint_velocity()[0], sign * speed);
            close(next.joint_position()[0], sign * distance);
            close(audit.friction_heat_j, loss);
            close(audit.thermal.external_heat_j, loss);
            close(audit.combined_residual_j, 0.0);
            close(
                next.thermal().temperature_k(CellId(1)).unwrap(),
                300.0 + 0.75 * loss / 10.0,
            );
            close(
                next.thermal().temperature_k(CellId(0)).unwrap(),
                300.0 + 0.25 * loss / 20.0,
            );
            let state = next.dry_joint_brakes().next().unwrap();
            let response = state.rate_response.as_ref().unwrap();
            close(response.resolved.state_coordinate("speed").unwrap(), speed);
            close(
                response.law.zero_speed_kinetic_coefficient().unwrap(),
                a + b * speed,
            );
            close(state.law.zero_speed_kinetic_coefficient().unwrap(), a);
            close(state.dissipated_power_w, 2.0 * (a + b * speed) * speed);
            close(
                state.work.surface_a_heat_j() + state.work.surface_b_heat_j(),
                loss,
            );
        }
    }
}

#[test]
fn constant_rate_data_matches_existing_brake_including_regularization() {
    let solid = card(10.0, 0.0, 500.0);
    let cards = [&solid; 2];
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(cards);
    let point = context();
    let rate_source = interface(cards, 0.4, 0.0, 0.0, 500.0, true);
    let old_source = interface(cards, 0.4, 0.0, 0.0, 500.0, false);
    for speed in [-1.5, -0.005, 0.0, 0.005, 1.5] {
        let base = moving(&reference, cards, [&heat; 2], speed, 0.03125, 0.0);
        let rate = base
            .with_rate_dependent_dry_joint_brakes(
                vec![binding(&rate_source, cards, &point, 2.0)],
                || Ok(()),
            )
            .unwrap();
        let old = base
            .with_dry_joint_brakes(vec![binding(&old_source, cards, &point, 2.0).brake], || {
                Ok(())
            })
            .unwrap();
        let (rate, ra) = rate.advance_to(0.125, config(), || Ok(())).unwrap();
        let (old, oa) = old.advance_to(0.125, config(), || Ok(())).unwrap();
        close(rate.joint_velocity()[0], old.joint_velocity()[0]);
        close(rate.joint_position()[0], old.joint_position()[0]);
        close(ra.friction_heat_j, oa.friction_heat_j);
        for cell in [CellId(0), CellId(1)] {
            close(
                rate.thermal().temperature_k(cell).unwrap(),
                old.thermal().temperature_k(cell).unwrap(),
            );
        }
    }
    // The original ramp uses mu(v_reg), including for a rate-dependent law.
    // Sampling mu(v) inside the ramp would produce a different nonlinear ODE.
    let varying = interface(cards, 0.4, 0.5, 0.0, 500.0, true);
    let initial = moving(&reference, cards, [&heat; 2], 0.005, 0.03125, 0.0)
        .with_rate_dependent_dry_joint_brakes(
            vec![binding(&varying, cards, &point, 2.0)],
            || Ok(()),
        )
        .unwrap();
    let (next, _) = initial.advance_to(0.125, config(), || Ok(())).unwrap();
    let speed = 0.005 * (-40.5_f64 * 0.125).exp();
    close(next.joint_velocity()[0], speed);
    let response = next
        .dry_joint_brakes()
        .next()
        .unwrap()
        .rate_response
        .as_ref()
        .unwrap();
    close(response.traction_speed_m_s, 0.01);
    close(
        response
            .traction_law
            .zero_speed_kinetic_coefficient()
            .unwrap(),
        0.405,
    );
    close(
        response.law.zero_speed_kinetic_coefficient().unwrap(),
        0.4 + 0.5 * speed,
    );
}

#[test]
fn temperature_and_speed_feedback_refines_to_independent_energy_reduction_and_restarts() {
    let parent = card(10.0, 0.0, 500.0);
    let child = card(0.1, 0.0, 500.0);
    let cards = [&parent, &child];
    let hp = caloric(&parent, 1e5, SolidLiquidPhase::Solid);
    let hc = caloric(&child, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(cards);
    let point = context();
    let source = interface(cards, 0.2, 0.5, 0.02, 310.0, true);
    // C_A=.2 and m=2 imply T_A-300=3.75(v0²-v²). Substitution
    // into v'=-mu gives the Riccati equation v'=.075v²-.5v-.36875.
    let k: f64 = 0.075;
    let discriminant = (0.25 + 4.0 * k * 0.36875).sqrt();
    let plus = (0.5 + discriminant) / (2.0 * k);
    let minus = (0.5 - discriminant) / (2.0 * k);
    let ratio = (1.5 - plus) / (1.5 - minus) * discriminant.exp();
    let exact = (ratio * minus - plus) / (ratio - 1.0);
    let mut errors = Vec::new();
    for dt in [0.25, 0.125, 0.0625] {
        let initial = moving(&reference, cards, [&hp, &hc], 1.5, dt, 0.0)
            .with_rate_dependent_dry_joint_brakes(
                vec![binding(&source, cards, &point, 2.0)],
                || Ok(()),
            )
            .unwrap();
        let (next, audit) = initial.advance_to(1.0, config(), || Ok(())).unwrap();
        let speed = next.joint_velocity()[0];
        errors.push((speed - exact).abs());
        let temperature = next.thermal().temperature_k(CellId(1)).unwrap();
        close(temperature, 300.0 + 3.75 * (2.25 - speed * speed));
        let state = next.dry_joint_brakes().next().unwrap();
        let response = state.rate_response.as_ref().unwrap();
        close(response.resolved.state_coordinate("speed").unwrap(), speed);
        close(
            response.resolved.state_coordinate("TA").unwrap(),
            temperature,
        );
        close(
            response.law.zero_speed_kinetic_coefficient().unwrap(),
            0.2 + 0.5 * speed + 0.02 * (temperature - 300.0),
        );
        for (new, old) in response.resolved.properties().iter().zip(
            initial
                .dry_joint_brakes()
                .next()
                .unwrap()
                .resolved
                .properties(),
        ) {
            assert_eq!(new.answer().receipt.selected, old.answer().receipt.selected);
        }
        close(audit.friction_heat_j, 2.25 - speed * speed);
        close(audit.combined_residual_j, 0.0);
        let (prefix, _) = initial.advance_to(0.5, config(), || Ok(())).unwrap();
        let (resumed, _) = prefix.advance_to(1.0, config(), || Ok(())).unwrap();
        assert_eq!(resumed.joint_velocity(), next.joint_velocity());
        assert_eq!(resumed.joint_position(), next.joint_position());
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
fn stage_speed_limit_temperature_support_and_late_cancel_preserve_original_state() {
    let solid = card(0.1, 0.0, 500.0);
    let cards = [&solid; 2];
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(cards);
    let point = context();
    let source = interface(cards, 0.2, 0.25, 0.0, 500.0, true);
    let base = moving(&reference, cards, [&heat; 2], 1.5, 0.125, 4.0);
    assert!(
        base.with_rate_dependent_dry_joint_brakes(
            vec![binding(&source, cards, &point, 1.0)],
            || Ok(())
        )
        .is_err()
    );
    let bounded = base
        .with_rate_dependent_dry_joint_brakes(vec![binding(&source, cards, &point, 1.6)], || Ok(()))
        .unwrap();
    let state = bounded.dry_joint_brakes().next().unwrap().clone();
    assert!(bounded.advance_to(0.5, config(), || Ok(())).is_err());
    assert_eq!(bounded.dry_joint_brakes().next().unwrap(), &state);
    assert_eq!(bounded.joint_velocity(), &[1.5]);
    close(bounded.thermal().time_s(), 0.0);
    let (short, _) = bounded.advance_to(0.03125, config(), || Ok(())).unwrap();
    assert!(short.joint_velocity()[0] > 1.5);

    let narrow_source = interface(cards, 0.2, 0.25, 0.0, 300.01, true);
    let initial = moving(&reference, cards, [&heat; 2], 1.5, 0.125, 0.0)
        .with_rate_dependent_dry_joint_brakes(vec![binding(&source, cards, &point, 2.0)], || Ok(()))
        .unwrap();
    let narrow = moving(&reference, cards, [&heat; 2], 1.5, 0.125, 0.0)
        .with_rate_dependent_dry_joint_brakes(
            vec![binding(&narrow_source, cards, &point, 2.0)],
            || Ok(()),
        )
        .unwrap();
    assert!(narrow.advance_to(0.25, config(), || Ok(())).is_err());
    close(narrow.thermal().time_s(), 0.0);
    close(narrow.thermal().temperature_k(CellId(1)).unwrap(), 300.0);
    let mut calls = 0;
    initial
        .advance_to(0.5, config(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    assert!(calls > 4);
    let original = initial.dry_joint_brakes().next().unwrap().clone();
    let mut attempted = 0;
    assert!(
        initial
            .advance_to(0.5, config(), || {
                attempted += 1;
                if attempted == calls - 1 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    assert_eq!(initial.dry_joint_brakes().next().unwrap(), &original);
    assert_eq!(initial.joint_position(), &[0.0]);
    assert_eq!(initial.joint_velocity(), &[1.5]);
    for cell in [CellId(0), CellId(1)] {
        close(initial.thermal().temperature_k(cell).unwrap(), 300.0);
    }
}
