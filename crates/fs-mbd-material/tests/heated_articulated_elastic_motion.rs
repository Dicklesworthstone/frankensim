//! Source elasticity -> native articulated motion -> storage and finite-body heat.
#[path = "heating/common.rs"]
mod common;
use common::{budget, caloric, point, region};
use fs_evidence::ValidityDomain;
use fs_ga::{Mat3, Se3, Vec3};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::phase::{ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section_temperature::{
    ConcentricSectionTemperatureLaw, thermoelastic::ConcentricThermoelasticLaw,
};
use fs_material::state_point::*;
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_mbd::articulated::JointModel;
use fs_mbd_material::articulated::{MaterialLink, ResolvedArticulated, resolve_articulated};
use fs_mbd_material::thermal::heating::{
    HeatedArticulated, LinkHeating,
    articulated_motion::{
        Controls, HeatedArticulatedMotion,
        dry::{
            DryBrakeLimits, DryJointBrake, JointSlipGeometry, KINETIC_FRICTION_PROPERTY,
            STATIC_FRICTION_PROPERTY,
        },
        elastic::AxialJointMember,
    },
};
use fs_mbd_material::{RegionBinding, VolumeMoments};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_thermal::enthalpy::*;
use fs_tribo::{HeatPartition, InputAuthority};

fn provenance() -> Provenance {
    Provenance {
        source: "synthetic elastic-member fixture, not measured data".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn claim(name: &str, dims: Dims, value: f64, slope: f64, upper: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: if slope == 0.0 {
            PropertyValue::Scalar { value, dims }
        } else {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, value), (upper, value + slope * (upper - 300.0))],
                dims,
            }
        },
        validity: ValidityDomain::unconstrained()
            .with("T", 300.0, upper)
            .with("pressure", 1e5, 2e5),
        interpolation: if slope == 0.0 {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        },
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: provenance(),
    }
}
fn card(young: f64, slope: f64, alpha: f64, upper: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, 2.0),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            10.0,
        ),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            alpha,
        ),
    ] {
        claims
            .insert_claim(claim(name, dims, value, 0.0, 500.0))
            .unwrap();
    }
    claims
        .insert_claim(claim(
            YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            young,
            slope,
            upper,
        ))
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic bonded member".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn member<'a>(cards: &[&'a MaterialCard]) -> AxialJointMember<'a> {
    let requirements = [
        (
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        ),
        (
            YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        ),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            ScalarAdmissibility::Finite,
        ),
    ]
    .map(|(name, dims, range)| ScalarPropertyRequirement::try_new(name, dims, range).unwrap());
    let layers = cards
        .iter()
        .enumerate()
        .map(|(i, card)| ConcentricMaterialLayer {
            // Each actual layer has one cubic metre at reference length one metre.
            outer_radius_m: ((i + 1) as f64 / core::f64::consts::PI).sqrt(),
            material: resolve_material_state_point(
                card,
                &point(300.0, 1e5),
                &requirements,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap(),
        })
        .collect();
    let section = ConcentricElasticSection::try_new(0.0, layers).unwrap();
    let temperature =
        ConcentricSectionTemperatureLaw::try_new(section, cards, &vec!["T"; cards.len()]).unwrap();
    AxialJointMember {
        link: 1,
        cells: (1..=cards.len()).map(|i| CellId(i as u64)).collect(),
        law: ConcentricThermoelasticLaw::try_new(temperature).unwrap(),
        reference_length_m: 1.0,
        reference_joint_position_m: 0.0,
        maximum_absolute_strain: 0.1,
    }
}
fn mechanism(cards: &[&MaterialCard]) -> ResolvedArticulated {
    let p = point(300.0, 1e5);
    let keys = ["a", "b"];
    let make = |i, card| RegionBinding {
        key: keys[i],
        material: card,
        point: &p,
        // Reference annuli of length one along the joint axis. Their mass is
        // lumped into this rigid link; no distributed axial inertia is claimed.
        geometry: {
            let radii_squared = (2 * i + 1) as f64 / core::f64::consts::PI;
            let transverse = (3.0 * radii_squared + 1.0) / 12.0;
            VolumeMoments::new(
                1.0,
                Vec3::default(),
                Mat3 {
                    m: [
                        0.5 * radii_squared,
                        0.0,
                        0.0,
                        0.0,
                        transverse,
                        0.0,
                        0.0,
                        0.0,
                        transverse,
                    ],
                },
            )
            .unwrap()
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    let mount = [make(0, cards[0])];
    let child: Vec<_> = cards[1..]
        .iter()
        .enumerate()
        .map(|(i, card)| make(i, *card))
        .collect();
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
                name: "member",
                parent: Some(0),
                parent_from_child_zero: Se3::identity(),
                joint: JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
                regions: &child,
            },
        ],
        3,
        || false,
    )
    .unwrap()
}
fn step() -> StepConfig {
    StepConfig {
        max_dt_s: 1.0,
        energy_tolerance_j: 1e-9,
    }
}
fn moving<'a>(
    reference: &'a ResolvedArticulated,
    cards: &[&'a MaterialCard],
    heat: &[&'a ResolvedHeatCapacityCurve],
    q: f64,
    v: f64,
    dt: f64,
    watts: f64,
) -> HeatedArticulatedMotion<'a> {
    let keys = ["a", "b"];
    let bindings = vec![
        LinkHeating {
            link: "mount",
            regions: vec![region(cards[0], heat[0], "a", 0, 0.0)],
        },
        LinkHeating {
            link: "member",
            regions: cards[1..]
                .iter()
                .zip(&heat[1..])
                .enumerate()
                .map(|(i, (card, heat))| region(*card, *heat, keys[i], (i + 1) as u64, watts))
                .collect(),
        },
    ];
    let heated =
        HeatedArticulated::new(reference, bindings, vec![], budget(), step(), || Ok(())).unwrap();
    HeatedArticulatedMotion::new(
        heated,
        Se3::identity(),
        vec![q],
        vec![v],
        vec![0.0],
        Vec3::default(),
        Controls {
            maximum_step_s: dt,
            relative_state_tolerance: 1e-11,
            absolute_position_tolerance: vec![1e-12],
            absolute_velocity_tolerance: vec![1e-12],
            maximum_rk_attempts_per_substep: 256,
            relative_momentum_tolerance: 1e-11,
            absolute_momentum_tolerance: vec![1e-12],
            maximum_parameter_work_j: 1.0,
            maximum_numerical_energy_change_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-9 * b.abs().max(1.0), "{a} != {b}");
}

#[test]
fn source_modulus_changes_native_oscillation_without_adding_mass_or_heat() {
    for young in [8.0_f64, 32.0] {
        let c = card(young, 0.0, 0.0, 500.0);
        let h = caloric(&c, 1e5, SolidLiquidPhase::Solid);
        let cards = [&c, &c];
        let reference = mechanism(&cards);
        let bare = moving(&reference, &cards, &[&h, &h], 0.02, 0.0, 0.125, 0.0);
        let old = bare
            .with_axial_joint_members(vec![member(&[&c])], || Ok(()))
            .unwrap();
        let (next, audit) = old.advance_to(0.5, step(), || Ok(())).unwrap();
        let omega = (young / 2.0).sqrt();
        near(next.joint_position()[0], 0.02 * (omega * 0.5).cos());
        near(
            next.joint_velocity()[0],
            -0.02 * omega * (omega * 0.5).sin(),
        );
        near(
            audit.mechanical_after.kinetic_j + audit.mechanical_after.elastic_j,
            0.5 * young * 0.02_f64.powi(2),
        );
        near(
            audit.mechanical_after.elastic_j,
            0.5 * young * next.joint_position()[0].powi(2),
        );
        near(audit.parameter_work_j, 0.0);
        near(audit.thermal.external_heat_j, 0.0);
        near(audit.combined_residual_j, 0.0);
        assert_eq!(next.thermal().mass_kg(), bare.thermal().mass_kg());
        let (middle, _) = old.advance_to(0.25, step(), || Ok(())).unwrap();
        let (restarted, _) = middle.advance_to(0.5, step(), || Ok(())).unwrap();
        assert_eq!(restarted.joint_position(), next.joint_position());
        assert_eq!(restarted.joint_velocity(), next.joint_velocity());
    }
}

#[test]
fn temperature_changes_real_restoring_force_and_accounts_parameter_work_separately() {
    let dt = 0.125;
    for (slope, alpha) in [(0.0, 0.0), (-4.0, 0.01)] {
        let c = card(8.0, slope, alpha, 301.0);
        let h = caloric(&c, 1e5, SolidLiquidPhase::Solid);
        let cards = [&c, &c];
        let reference = mechanism(&cards);
        let old = moving(&reference, &cards, &[&h, &h], 0.02, 0.0, dt, 20.0)
            .with_axial_joint_members(vec![member(&[&c])], || Ok(()))
            .unwrap();
        let (first, first_audit) = old.advance_to(dt, step(), || Ok(())).unwrap();
        let q1 = 0.02 * (2.0 * dt).cos();
        let v1 = -0.04 * (2.0 * dt).sin();
        let k1 = 8.0 + slope * dt;
        let eq1 = alpha * dt;
        near(first.joint_position()[0], q1);
        near(first.joint_velocity()[0], v1);
        near(
            first_audit.parameter_work_j,
            0.5 * k1 * (q1 - eq1).powi(2) - 4.0 * q1 * q1,
        );
        let (second, audit) = first.advance_to(2.0 * dt, step(), || Ok(())).unwrap();
        let w = (k1 / 2.0).sqrt();
        let q2 = eq1 + (q1 - eq1) * (w * dt).cos() + v1 / w * (w * dt).sin();
        let v2 = -(q1 - eq1) * w * (w * dt).sin() + v1 * (w * dt).cos();
        near(second.joint_position()[0], q2);
        near(second.joint_velocity()[0], v2);
        let k2 = 8.0 + slope * 2.0 * dt;
        let eq2 = alpha * 2.0 * dt;
        near(
            audit.parameter_work_j,
            0.5 * k2 * (q2 - eq2).powi(2) - 0.5 * k1 * (q2 - eq1).powi(2),
        );
        near(audit.thermal.external_heat_j, 20.0 * dt);
        near(
            second.thermal().energy_j() - first.thermal().energy_j(),
            20.0 * dt,
        );
        near(audit.combined_residual_j, 0.0);
        let endpoint = second.axial_joint_members().next().unwrap();
        near(endpoint.restoring_effort_n, -k2 * (q2 - eq2));
        near(endpoint.phases[0].temperature_k(), 300.0 + 2.0 * dt);
        assert_eq!(
            endpoint.section.section().layers()[0]
                .material
                .card_identity(),
            c.content_hash()
        );
    }
    // Two bonded layers retain unequal thermal stress even at zero resultant.
    let a = card(8.0, 0.0, 0.01, 500.0);
    let b = card(16.0, 0.0, 0.02, 500.0);
    let ca = caloric(&a, 1e5, SolidLiquidPhase::Solid);
    let cb = caloric(&b, 1e5, SolidLiquidPhase::Solid);
    let cards = [&a, &a, &b];
    let reference = mechanism(&cards);
    let old = moving(&reference, &cards, &[&ca, &ca, &cb], 0.0, 0.0, dt, 20.0)
        .with_axial_joint_members(vec![member(&[&a, &b])], || Ok(()))
        .unwrap();
    let (warm, _) = old.advance_to(dt, step(), || Ok(())).unwrap();
    let state = warm.axial_joint_members().next().unwrap();
    let eq = dt * (8.0 * 0.01 + 16.0 * 0.02) / 24.0;
    let (force, energy) = state
        .section
        .axial_force_and_energy_at_strain(eq, 0.1)
        .unwrap();
    near(force, 0.0);
    near(
        energy,
        4.0 * (eq - 0.01 * dt).powi(2) + 8.0 * (eq - 0.02 * dt).powi(2),
    );
    assert!(energy > 0.0);
    let (next, _) = warm.advance_to(2.0 * dt, step(), || Ok(())).unwrap();
    near(
        next.joint_position()[0],
        eq * (1.0 - (6.0_f64.sqrt() * dt).cos()),
    );
    near(
        next.joint_velocity()[0],
        eq * 6.0_f64.sqrt() * (6.0_f64.sqrt() * dt).sin(),
    );
}

fn interface(c: &MaterialCard) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for (name, value) in [
        (STATIC_FRICTION_PROPERTY, 0.3),
        (KINETIC_FRICTION_PROPERTY, 0.2),
    ] {
        let mut claim = claim(name, Dims::NONE, value, 0.0, 500.0);
        claim.validity = ValidityDomain::unconstrained()
            .with_quantity(
                "TA",
                QuantitySpec::dimensional(Temperature::DIMS),
                300.0,
                500.0,
            )
            .with_quantity(
                "TB",
                QuantitySpec::dimensional(Temperature::DIMS),
                300.0,
                500.0,
            )
            .with_quantity(
                "pressure",
                QuantitySpec::dimensional(Pressure::DIMS),
                2.0,
                2.0,
            );
        claims.insert_claim(claim).unwrap();
    }
    let surface = SurfaceSpec {
        material: c.id().clone(),
        texture_frame: "fixture".into(),
    };
    InterfaceSystemCard::assemble(
        surface.clone(),
        surface,
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "fixture".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}
#[test]
fn native_dry_brake_deposits_the_loss_of_elastic_plus_kinetic_energy_once() {
    let c = card(8.0, 0.0, 0.0, 500.0);
    let h = caloric(&c, 1e5, SolidLiquidPhase::Solid);
    let cards = [&c, &c];
    let reference = mechanism(&cards);
    let source = interface(&c);
    let point = QueryPoint::new()
        .with_quantity("TA", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("TB", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("pressure", QuantitySpec::dimensional(Pressure::DIMS), 2.0)
        .unwrap();
    let brake = DryJointBrake {
        link: 1,
        cells: [CellId(1), CellId(0)],
        geometry: JointSlipGeometry::Prismatic,
        normal_force_n: 2.0,
        nominal_area_m2: 1.0,
        pressure_axis: "pressure",
        point: ThermalInterfacePoint {
            interface: &source,
            materials: [&c, &c],
            context: &point,
            temperature_axes: ["TA", "TB"],
            liquid_fraction_axes: [None, None],
            mean_temperature_axis: None,
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        authority: InputAuthority::SyntheticFixture,
        partition: HeatPartition::new(0.75, 0.25, 0.0).unwrap(),
        limits: DryBrakeLimits {
            regularization_m_s: 1.0,
            maximum_coefficient_change: 0.1,
            absolute_loss_tolerance_j: 1e-12,
        },
    };
    let old = moving(&reference, &cards, &[&h, &h], 0.02, 0.0, 0.125, 0.0)
        .with_axial_joint_members(vec![member(&[&c])], || Ok(()))
        .unwrap()
        .with_dry_joint_brakes(vec![brake], || Ok(()))
        .unwrap();
    let (next, audit) = old.advance_to(0.5, step(), || Ok(())).unwrap();
    let beta = 0.1_f64;
    let w = (4.0 - beta * beta).sqrt();
    let q = 0.02 * (-beta * 0.5).exp() * ((w * 0.5).cos() + beta / w * (w * 0.5).sin());
    let v = -0.02 * (-beta * 0.5).exp() * 4.0 / w * (w * 0.5).sin();
    let loss = 0.0016 - v * v - 4.0 * q * q;
    near(next.joint_position()[0], q);
    near(next.joint_velocity()[0], v);
    near(audit.friction_heat_j, loss);
    near(audit.thermal.stored_energy_change_j, loss);
    near(
        next.thermal().cells()[&CellId(1)].specific_h_j_kg,
        0.75 * loss / 2.0,
    );
    near(
        next.thermal().cells()[&CellId(0)].specific_h_j_kg,
        0.25 * loss / 2.0,
    );
    near(audit.parameter_work_j, 0.0);
    near(audit.combined_residual_j, 0.0);
    assert!(loss > 0.0);
}

#[test]
fn actual_inventory_source_strain_and_cancellation_refusals_leave_clean_retry() {
    let c = card(8.0, 0.0, 0.0, 300.001);
    let h = caloric(&c, 1e5, SolidLiquidPhase::Solid);
    let cards = [&c, &c];
    let reference = mechanism(&cards);
    let bare = moving(&reference, &cards, &[&h, &h], 0.02, 0.0, 0.125, 20.0);
    let mut bad = member(&[&c]);
    bad.reference_length_m = 2.0;
    assert!(bare.with_axial_joint_members(vec![bad], || Ok(())).is_err());
    let good = member(&[&c]);
    assert!(
        bare.with_axial_joint_members(vec![good.clone(), good.clone()], || Ok(()))
            .is_err()
    );
    let old = bare
        .with_axial_joint_members(vec![good], || Ok(()))
        .unwrap();
    let saved = old.axial_joint_members().next().unwrap().clone();
    assert!(old.advance_to(0.125, step(), || Ok(())).is_err());
    let mut calls = 0;
    let (short, _) = old
        .advance_to(1.0 / 4096.0, step(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    let mut now = 0;
    assert!(
        old.advance_to(1.0 / 4096.0, step(), || {
            now += 1;
            if now == calls - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        })
        .is_err()
    );
    let (retry, _) = old.advance_to(1.0 / 4096.0, step(), || Ok(())).unwrap();
    assert_eq!(short.joint_position(), retry.joint_position());
    assert_eq!(short.thermal().cells(), retry.thermal().cells());
    assert_eq!(&saved, old.axial_joint_members().next().unwrap());
    assert_eq!(old.thermal().time_s(), 0.0);
    let mut narrow = member(&[&c]);
    narrow.maximum_absolute_strain = 0.021;
    let driven = moving(&reference, &cards, &[&h, &h], 0.02, 0.2, 0.125, 0.0)
        .with_axial_joint_members(vec![narrow], || Ok(()))
        .unwrap();
    assert!(driven.advance_to(0.125, step(), || Ok(())).is_err());
    assert!(driven.advance_to(1.0 / 4096.0, step(), || Ok(())).is_ok());
    assert_eq!(driven.joint_position(), &[0.02]);
}
