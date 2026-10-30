//! G1/G3/G4: original directional claims, native recoil and conduction-derived heat.
//! Material values are authored analytical fixtures, not experimental calibration.
use fs_contact::interface_binding::bind_dry_interface_system_card;
use fs_contact::interface_binding::thermal_friction::elliptic::SourceThermalEllipticCoulomb;
use fs_contact::interface_binding::thermal_friction::elliptic::motion::MotionQuerySolve;
use fs_contact::interface_binding::thermal_friction::{
    FrictionQueryAxis, FrictionQueryAxisKind, SourceThermalCoulomb, ThermalFrictionAxes,
};
use fs_couple::friction_heat::directional::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::*;
use fs_material::phase::{
    EquilibriumPhaseState, HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, SPECIFIC_HEAT_CAPACITY_DIMS,
    ScalarAdmissibility, ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS,
    resolve_material_state_point,
};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_material::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::directional::{RigidDirectionalContact, RigidDirectionalControls};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::{InputAuthority, TangentialLaw, TangentialSolveConfig};
use std::collections::BTreeMap;

#[path = "rigid_directional_heat/articulated.rs"]
mod articulated;
#[path = "rigid_directional_heat/driven.rs"]
mod driven;
#[path = "rigid_directional_heat/impacts.rs"]
mod impacts;
#[path = "rigid_directional_heat/sphere_plane.rs"]
mod sphere_plane;
#[path = "rigid_directional_heat/trajectory.rs"]
mod trajectory;

const CELLS: [CellId; 2] = [CellId(30), CellId(10)];
const CAPACITY: [f64; 2] = [0.02, 0.06];
const CONDUCTIVITY: [f64; 2] = [1.0, 3.0];
const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 2e-8 * a.abs().max(b.abs()).max(1.0),
        "{a:e} != {b:e}"
    );
}
fn provenance() -> Provenance {
    Provenance {
        source: "synthetic directional heat reference".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn scalar(x: f64) -> PropertyValue {
    PropertyValue::Scalar {
        value: x,
        dims: Dims::NONE,
    }
}
fn curve(axis: &str, dims: Dims, knots: [(f64, f64); 2]) -> PropertyValue {
    PropertyValue::Curve {
        abscissa: axis.into(),
        abscissa_dims: dims,
        knots: knots.to_vec(),
        dims: Dims::NONE,
    }
}
fn material_query(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn material(side: usize) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let cp = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let rho = PropertyKey::new("density", Density::DIMS);
    let mut claims = ClaimSet::new();
    for (key, value) in [
        (cp.clone(), 1.0),
        (rho.clone(), CAPACITY[side]),
        (
            PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            CONDUCTIVITY[side],
        ),
    ] {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value), (400.0, value)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    200.0,
                    400.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("synthetic-{side}"),
            phase: "solid".into(),
            process: "analytical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let heat = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &cp,
            density: &rho,
            temperature_axis: "T",
            lower: &material_query(200.0),
            upper: &material_query(400.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-8,
    )
    .unwrap();
    (card, heat)
}
fn query() -> QueryPoint {
    let mut q = QueryPoint::new();
    for axis in ["Ti", "Ta", "Tb"] {
        q = q
            .with_quantity(axis, absolute_temperature(), 300.0)
            .unwrap();
    }
    q.with_quantity("pressure", static_pressure(), 10.0)
        .unwrap()
        .with_quantity("speed", QuantitySpec::dimensional(SPEED), 1.0)
        .unwrap()
        .with_quantity("path", QuantitySpec::dimensional(LENGTH), 0.25)
        .unwrap()
}
fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn static_pressure() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}
fn interface(
    cards: &[MaterialCard; 2],
    kinetic: [PropertyValue; 2],
    temperature_limit: f64,
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for (name, value) in [
        ("static-x", scalar(0.9)),
        ("static-y", scalar(0.9)),
        ("kinetic-x", kinetic[0].clone()),
        ("kinetic-y", kinetic[1].clone()),
    ] {
        let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        };
        let mut validity = ValidityDomain::unconstrained();
        for axis in ["Ti", "Ta", "Tb"] {
            validity =
                validity.with_quantity(axis, absolute_temperature(), 200.0, temperature_limit);
        }
        validity = validity
            .with_quantity("pressure", static_pressure(), 0.0, 100.0)
            .with_quantity("speed", QuantitySpec::dimensional(SPEED), 0.0, 20.0)
            .with_quantity("path", QuantitySpec::dimensional(LENGTH), 0.0, 10.0);
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, Dims::NONE),
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "surface-A".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "surface-B".into(),
        },
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "synthetic".into(),
            history: "original".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn law(card: &InterfaceSystemCard, q: &QueryPoint, axis: [f64; 2]) -> SourceThermalEllipticCoulomb {
    let mut axes = ["Ta", "Tb"]
        .map(|name| FrictionQueryAxis {
            name: name.into(),
            kind: FrictionQueryAxisKind::Temperature,
        })
        .to_vec();
    axes.extend([
        FrictionQueryAxis {
            name: "speed".into(),
            kind: FrictionQueryAxisKind::SlidingSpeed,
        },
        FrictionQueryAxis {
            name: "path".into(),
            kind: FrictionQueryAxisKind::AccumulatedDistance,
        },
    ]);
    let sources = [0, 1].map(|i| {
        SourceThermalCoulomb::resolve_with_query_axes(
            bind_dry_interface_system_card(card, InputAuthority::SyntheticFixture).unwrap(),
            q.clone(),
            ScalarPropertyRequirement::try_new(
                ["static-x", "static-y"][i],
                Dims::NONE,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                ["kinetic-x", "kinetic-y"][i],
                Dims::NONE,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            ThermalFrictionAxes {
                temperature: "Ti".into(),
                nominal_pressure: Some("pressure".into()),
            },
            axes.clone(),
        )
        .unwrap()
    });
    SourceThermalEllipticCoulomb::from_principal_sources(sources, axis).unwrap()
}
fn binding<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    card: &'a InterfaceSystemCard,
    q: &'a QueryPoint,
) -> DirectionalHeatMaterial<'a> {
    DirectionalHeatMaterial {
        point: ThermalInterfacePoint {
            interface: card,
            materials: [&cards[0], &cards[1]],
            context: q,
            temperature_axes: ["Ta", "Tb"],
            liquid_fraction_axes: [None, None],
            mean_temperature_axis: None,
        },
        caloric: [curves[0].curve(), curves[1].curve()],
        interface_temperature_axis: Some("Ti"),
    }
}
fn thermal(
    curves: &[ResolvedHeatCapacityCurve; 2],
    temperatures: [f64; 2],
) -> DepositionSimulation<'_> {
    let mut events = Vec::new();
    for i in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELLS[i],
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: CAPACITY[i],
                specific_h_j_kg: temperatures[i] - 200.0,
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, i as f64 - 1.0],
                    end_m: [0.0, 0.0, i as f64],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: CELLS[0],
            b: CELLS[1],
            area_m2: 1.0,
            conductance_w_k: series_conductance_w_k(1.0, [0.5, 0.5], CONDUCTIVITY).unwrap(),
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn transport(
    cards: &[MaterialCard; 2],
    phases: [EquilibriumPhaseState; 2],
) -> Result<DirectionalHeatSource<[ResolvedMaterialStatePoint; 2]>, ThermalMaterialError> {
    let k = ScalarPropertyRequirement::try_new(
        "thermal_conductivity",
        THERMAL_CONDUCTIVITY_DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )
    .map_err(ThermalMaterialError::Material)?;
    let at = |i: usize| {
        let base = material_query(300.0);
        let q = ThermalMaterialPoint {
            material: &cards[i],
            context: &base,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        }
        .at(phases[i])?;
        resolve_material_state_point(
            &cards[i],
            &q,
            &[k.clone()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .map_err(ThermalMaterialError::Material)
    };
    let source = [at(0)?, at(1)?];
    Ok(DirectionalHeatSource {
        conductivity_w_m_k: source.each_ref().map(|s| s.properties()[0].value_si()),
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: source,
        },
    })
}
fn mechanics(
    full: bool,
) -> (
    [RigidBodyState; 3],
    [MassProperties; 3],
    RigidDirectionalContact,
) {
    let positions = if full {
        [Vec3::new(-0.2, -0.3, -0.7), Vec3::new(0.4, -0.1, 0.2)]
    } else {
        [Vec3::new(0.0, 0.0, -0.1), Vec3::new(0.0, 0.0, 0.1)]
    };
    let rotations = [
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.7).unwrap(),
        UnitQuaternion::from_axis_angle(Vec3::new(2.0, -1.0, 1.0), -0.4).unwrap(),
    ];
    let p = [0, 1].map(|i| {
        if full {
            MassProperties::from_principal_frame(
                2.0 + i as f64,
                Vec3::ZERO,
                Vec3::new(1.0, 1.5, 2.0),
                rotations[i],
            )
            .unwrap()
        } else {
            MassProperties::new(
                2.0 + i as f64,
                Vec3::ZERO,
                Vec3::new(2.0 + i as f64, 2.0 + i as f64, 2.0 + i as f64),
            )
            .unwrap()
        }
    });
    let bodies = [0, 1].map(|i| {
        RigidBodyState::new(
            Pose::new(positions[i], UnitQuaternion::IDENTITY).unwrap(),
            if i == 0 {
                Vec3::new(2.0, if full { 0.8 } else { 0.0 }, 0.0)
            } else {
                Vec3::ZERO
            },
            Vec3::ZERO,
        )
        .unwrap()
    });
    let spare =
        RigidBodyState::new(Pose::identity(), Vec3::new(3.0, 2.0, 1.0), Vec3::ZERO).unwrap();
    (
        [bodies[0], bodies[1], spare],
        [p[0], p[1], p[0]],
        RigidDirectionalContact {
            binding: RigidContactBinding {
                body_a: 0,
                body_b: 1,
                arm_a_body: positions[0].scale(-1.0),
                arm_b_body: positions[1].scale(-1.0),
            },
            frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
            normal_load_n: 10.0,
            contact_temperature_k: None,
            nominal_pressure_pa: None,
        },
    )
}
fn controls() -> DirectionalHeatSolve {
    DirectionalHeatSolve {
        mechanical: RigidDirectionalControls {
            maximum_bodies: 4,
            geometry_tolerance_m: 1e-10,
            velocity_tolerance_m_s: 1e-8,
            linear_momentum_tolerance_ns: 1e-9,
            angular_momentum_tolerance_nms: 1e-9,
            energy_tolerance_j: 1e-8,
            contact: TangentialSolveConfig {
                max_iterations: 128,
                maximum_force_n: 1e5,
                force_absolute_tolerance_n: 1e-11,
                force_relative_tolerance: 1e-13,
            },
        },
        motion: MotionQuerySolve {
            distance_bracket_m: [0.0, 0.02],
            absolute_distance_tolerance_m: 1e-12,
            relative_distance_tolerance: 1e-12,
            maximum_evaluations: 96,
        },
        maximum_cells: 2,
        maximum_trials: 300,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-6,
        temperature_scale_j_per_k: 0.02,
        heat: FrictionHeatTolerances {
            cell_j: 1e-8,
            thermal_j: 1e-8,
            interface_j: 1e-8,
            total_j: 1e-8,
        },
        combined_energy_tolerance_j: 2e-8,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-5,
            linear_restart: 3,
            ..NewtonKrylovConfig::default()
        },
    }
}

#[test]
fn principal_recoil_and_finite_heat_match_independent_closed_forms() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = interface(&cards, [scalar(0.2), scalar(0.4)], 400.0);
    let law = law(&card, &q, [1.0, 0.0]);
    let mut state = law.initial_state();
    let initial = state.clone();
    let mut thermal = thermal(&curves, [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, 7, CONDUCTIVITY).unwrap();
    let (mut bodies, props, c) = mechanics(false);
    let spare = bodies[2];
    let h = 0.01;
    let w = 1.0 / 2.0 + 1.0 / 3.0 + 0.01 / 2.0 + 0.01 / 3.0;
    let distance = h - 0.5 * h * h * w * 2.0;
    let work = 2.0 * distance;
    let report = apply_source_directional_heat(
        &mut bodies,
        &props,
        &mut state,
        &law,
        c,
        &mut thermal,
        &route,
        binding(&cards, &curves, &card, &q),
        h,
        controls(),
        &[0.0; 2],
        300.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    close(report.mechanical.mobility_per_kg[0][0], w);
    close(report.mechanical.contact.traction_n[0], 2.0);
    close(report.actual_distance_m, distance);
    close(report.heat.dissipated_j, work);
    close(report.mechanical.kinetic_energy_change_j, -work);
    close(report.combined_energy_residual_j, 0.0);
    close(report.heat.body_heat_j[0], 0.25 * work);
    close(report.heat.body_heat_j[1], 0.75 * work);
    close(report.heat.interface_temperature_k, 300.0 + work / 0.04);
    for id in CELLS {
        close(thermal.temperature_k(id).unwrap(), 300.0 + work / 0.08);
    }
    assert_eq!(bodies[2], spare);
    for (i, source) in state.friction_sources().iter().enumerate() {
        close(
            source.state_coordinate("Ti").unwrap(),
            report.heat.interface_temperature_k,
        );
        close(source.state_coordinate("pressure").unwrap(), 10.0);
        for p in source.properties() {
            assert_eq!(
                p.answer().receipt.selected,
                initial.friction_sources()[i]
                    .property(p.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            card.claims().verify_receipt(&p.answer().receipt).unwrap();
        }
    }
    for i in 0..2 {
        close(
            report.source[i].query().axes()["T"],
            report.phases[i].temperature_k(),
        );
    }
}

#[test]
fn solved_interface_temperature_weakens_the_original_friction_in_the_same_step() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = interface(
        &cards,
        [
            curve("Ti", Temperature::DIMS, [(299.0, 0.3), (301.0, 0.1)]),
            scalar(0.4),
        ],
        400.0,
    );
    let law = law(&card, &q, [1.0, 0.0]);
    let mut state = law.initial_state();
    let mut thermal = thermal(&curves, [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, 7, CONDUCTIVITY).unwrap();
    let (mut bodies, props, c) = mechanics(false);
    let h = 0.01;
    let w = 101.0 / 120.0;
    // Independent scalar elimination of the two finite cells and interface:
    // Ti-300 = D/0.04, F=2-(Ti-300), D=F*(h-h² W F/2).
    let mut lo = 0.0;
    let mut hi = 2.0;
    for _ in 0..80 {
        let f = 0.5 * (lo + hi);
        let d = f * (h - 0.5 * h * h * w * f);
        if f + d / 0.04 > 2.0 {
            hi = f;
        } else {
            lo = f;
        }
    }
    let f = 0.5 * (lo + hi);
    let report = apply_source_directional_heat(
        &mut bodies,
        &props,
        &mut state,
        &law,
        c,
        &mut thermal,
        &route,
        binding(&cards, &curves, &card, &q),
        h,
        controls(),
        &[0.0; 2],
        300.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    close(report.mechanical.contact.traction_n[0], f);
    assert!(f < 1.7);
    close(report.heat.interface_temperature_k, 302.0 - f);
    assert!(report.heat.interface_temperature_k > report.phases[0].temperature_k());
    close(report.combined_energy_residual_j, 0.0);
}

#[test]
fn full_tensor_motion_owned_queries_and_tangent_frame_rotation_agree() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = interface(
        &cards,
        [
            curve("speed", SPEED, [(0.0, 0.1), (20.0, 0.7)]),
            curve("path", LENGTH, [(0.0, 0.2), (10.0, 0.6)]),
        ],
        400.0,
    );
    let mut reference = None;
    for rotated in [false, true] {
        let (mut bodies, props, mut c) = mechanics(true);
        let axis = if rotated {
            c.frame = TangentFrame::new([0.0, 0.0, 1.0], [0.6, 0.8, 0.0]).unwrap();
            [1.0, 0.0]
        } else {
            [0.6, 0.8]
        };
        let law = law(&card, &q, axis);
        let mut state = law.initial_state();
        let mut heat = thermal(&curves, [300.0; 2]);
        let route = FrictionHeatRoute::from_link(&heat, 7, CONDUCTIVITY).unwrap();
        let h = 0.002;
        let report = apply_source_directional_heat(
            &mut bodies,
            &props,
            &mut state,
            &law,
            c,
            &mut heat,
            &route,
            binding(&cards, &curves, &card, &q),
            h,
            controls(),
            &[0.0; 2],
            300.0,
            |_, p| transport(&cards, p),
            &CancelGate::new(),
        )
        .unwrap();
        assert!(report.mechanical.mobility_per_kg[0][1].abs() > 1e-4);
        assert!(report.source_motion_residual_m <= 0.0);
        assert!(report.source_motion_residual_m.abs() < 2e-12);
        for source in state.friction_sources() {
            close(
                source.state_coordinate("speed").unwrap(),
                report.actual_distance_m / h,
            );
            close(
                source.state_coordinate("path").unwrap(),
                0.25 + report.actual_distance_m,
            );
            close(
                source.state_coordinate("Ta").unwrap(),
                heat.temperature_k(CELLS[0]).unwrap(),
            );
            close(
                source.state_coordinate("Tb").unwrap(),
                heat.temperature_k(CELLS[1]).unwrap(),
            );
        }
        let values = [
            report.mechanical.force_impulse_on_a_world_ns.x,
            report.mechanical.force_impulse_on_a_world_ns.y,
            bodies[0].linear_momentum_world().x,
            bodies[0].linear_momentum_world().y,
            report.heat.dissipated_j,
            report.heat.interface_temperature_k,
        ];
        if let Some(old) = reference {
            for (a, b) in values.into_iter().zip(old) {
                close(a, b);
            }
        } else {
            reference = Some(values);
        }
        close(report.combined_energy_residual_j, 0.0);
    }
}

#[test]
fn unequal_surface_temperatures_exchange_heat_without_friction_power() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = interface(&cards, [scalar(0.2), scalar(0.4)], 400.0);
    let law = law(&card, &q, [1.0, 0.0]);
    let mut state = law.initial_state();
    let mut heat = thermal(&curves, [280.0, 320.0]);
    let route = FrictionHeatRoute::from_link(&heat, 7, CONDUCTIVITY).unwrap();
    let (mut bodies, props, mut c) = mechanics(false);
    let original = bodies;
    c.normal_load_n = 0.0;
    let report = apply_source_directional_heat(
        &mut bodies,
        &props,
        &mut state,
        &law,
        c,
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        0.01,
        controls(),
        &[0.0; 2],
        310.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, original);
    close(report.heat.dissipated_j, 0.0);
    close(report.heat.body_heat_j[0], 0.3);
    close(report.heat.body_heat_j[1], -0.3);
    close(heat.temperature_k(CELLS[0]).unwrap(), 295.0);
    close(heat.temperature_k(CELLS[1]).unwrap(), 315.0);
    close(report.heat.interface_temperature_k, 310.0);
    close(report.combined_energy_residual_j, 0.0);
}

#[test]
fn failures_cancel_and_retry_preserve_all_original_owners() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = interface(&cards, [scalar(0.2), scalar(0.4)], 400.0);
    let law = law(&card, &q, [1.0, 0.0]);
    let initial = law.initial_state();
    let base = thermal(&curves, [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&base, 7, CONDUCTIVITY).unwrap();
    let (original, props, c) = mechanics(false);
    for mode in 0..5 {
        let mut bodies = original;
        let mut state = initial.clone();
        let mut heat = base.clone();
        let gate = CancelGate::new();
        let mut solve = controls();
        let mut geometry = c;
        if mode == 2 {
            solve.maximum_trials = 1;
        }
        if mode == 3 {
            geometry.binding.arm_a_body.x += 0.1;
        }
        let mut materials = binding(&cards, &curves, &card, &q);
        if mode == 4 {
            materials.caloric.swap(0, 1);
        }
        let result = apply_source_directional_heat(
            &mut bodies,
            &props,
            &mut state,
            &law,
            geometry,
            &mut heat,
            &route,
            materials,
            0.01,
            solve,
            &[0.0; 2],
            300.0,
            |_, p| {
                if mode == 0 {
                    return Err(ThermalMaterialError::Invalid("transport refused"));
                }
                if mode == 1 {
                    gate.request();
                }
                transport(&cards, p)
            },
            &gate,
        );
        assert!(result.is_err(), "failure mode {mode} unexpectedly accepted");
        assert_eq!(bodies, original);
        assert_eq!(state, initial);
        assert_eq!(heat.cells(), base.cells());
        assert_eq!(heat.time_s(), base.time_s());
    }
    let mut bodies = original;
    let mut state = initial.clone();
    let mut heat = base.clone();
    let report = apply_source_directional_heat(
        &mut bodies,
        &props,
        &mut state,
        &law,
        c,
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        0.01,
        controls(),
        &[0.0; 2],
        300.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    close(report.combined_energy_residual_j, 0.0);
    assert_ne!(state, initial);
    assert_eq!(bodies[2], original[2]);
}
