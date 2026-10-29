#![allow(dead_code)]
use fs_contact::capillary_bridge::BridgeLimits;
use fs_contact::capillary_bridge::asymmetric_contact_line::*;
use fs_contact::capillary_bridge::contact_line::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::liquid::surface::adhesion::{
    EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY, EquilibriumAngleOrigin,
};
use fs_material::liquid::surface::*;
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Density, Dims};

pub const VOLUME: f64 = 1e-10;
pub const WIDTH: f64 = 1e-3;
pub const SIGMA: f64 = 0.05;
pub const GAP: f64 = 3e-4;

fn claim(name: &str, value: f64, dims: Dims) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with("T", 280.0, 330.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: format!("manufactured asymmetric wetting: {name}"),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

pub fn sources(
    theta: [f64; 2],
    friction: [f64; 2],
    receding: [f64; 2],
    advancing: [f64; 2],
) -> (ResolvedLiquidSurfaceState, [ResolvedLiquidWettingState; 2]) {
    sources_with_origins(
        theta,
        friction,
        receding,
        advancing,
        [EquilibriumAngleOrigin::SourcedAngle; 2],
    )
}
pub fn sources_with_origins(
    theta: [f64; 2],
    friction: [f64; 2],
    receding: [f64; 2],
    advancing: [f64; 2],
    origins: [EquilibriumAngleOrigin; 2],
) -> (ResolvedLiquidSurfaceState, [ResolvedLiquidWettingState; 2]) {
    let id = |name: &str, phase: &str| MaterialStateId {
        chemistry: name.into(),
        phase: phase.into(),
        process: "manufactured unlike-surfaces fixture".into(),
        revision: 0,
    };
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density", 1000.0, Density::DIMS),
        (DYNAMIC_VISCOSITY_PROPERTY, 1e-3, DYNAMIC_VISCOSITY_DIMS),
        (SURFACE_TENSION_PROPERTY, SIGMA, SURFACE_TENSION_DIMS),
    ] {
        claims.insert_claim(claim(name, value, dims)).unwrap();
    }
    let card = MaterialCard::assemble(id("liquid", "liquid"), claims, vec![]).unwrap();
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    let liquid = resolve_liquid_surface_state(
        &card,
        &point,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let wetting = core::array::from_fn(|i| {
        let mut claims = ClaimSet::new();
        let equilibrium = match origins[i] {
            EquilibriumAngleOrigin::SourcedAngle => {
                (EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, theta[i], Dims::NONE)
            }
            EquilibriumAngleOrigin::YoungDupre => (
                EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY,
                SIGMA * (1.0 + theta[i].cos()),
                SURFACE_TENSION_DIMS,
            ),
        };
        for (name, value, dims) in [
            equilibrium,
            (
                CONTACT_LINE_FRICTION_PROPERTY,
                friction[i],
                DYNAMIC_VISCOSITY_DIMS,
            ),
            (RECEDING_CONTACT_ANGLE_PROPERTY, receding[i], Dims::NONE),
            (ADVANCING_CONTACT_ANGLE_PROPERTY, advancing[i], Dims::NONE),
        ] {
            claims.insert_claim(claim(name, value, dims)).unwrap();
        }
        let interface = InterfaceSystemCard::assemble(
            SurfaceSpec {
                material: id(&format!("plate {i}"), "solid"),
                texture_frame: format!("plate {i}"),
            },
            SurfaceSpec {
                material: card.id().clone(),
                texture_frame: "liquid".into(),
            },
            SystemContext {
                medium: "declared vapor".into(),
                third_body: None,
                environment: "isothermal asymmetric fixture".into(),
                history: format!("original surface {i}"),
            },
            claims,
            vec![],
        )
        .unwrap();
        resolve_contact_line_wetting(
            &interface,
            &liquid,
            origins[i],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
    });
    (liquid, wetting)
}
pub fn limits() -> BridgeLimits {
    BridgeLimits {
        minimum_gap_m: 1e-5,
        maximum_gap_m: 0.02,
        maximum_wetted_span_m: 0.02,
        maximum_abs_force_n: 0.1,
        maximum_abs_surface_energy_j: 1e-3,
        maximum_abs_step_work_j: 1e-3,
        maximum_transition_release_j: 1e-3,
        energy_absolute_tolerance_j: 3e-19,
        energy_relative_tolerance: 3e-11,
        gravity_m_s2: 0.0,
        maximum_bond_number: 0.01,
        maximum_capillary_number: 0.005,
        maximum_weber_number: 0.005,
    }
}
pub fn config() -> AsymmetricContactLineConfig {
    AsymmetricContactLineConfig {
        maximum_backtracks: 24,
        contact_line: ContactLineConfig {
            minimum_angle_radians: 0.2,
            maximum_angle_radians: 1.4,
            minimum_neck_width_m: 1e-7,
            maximum_gap_increment_m: 4e-5,
            maximum_span_increment_m: 8e-5,
            maximum_step_s: 0.2,
            maximum_contact_line_speed_m_s: 0.02,
            maximum_angle_departure_radians: 0.5,
            maximum_iterations: 32,
            maximum_geometry_iterations: 40,
            maximum_evaluations: 8192,
            geometry_absolute_tolerance_m: 1e-17,
            geometry_relative_tolerance: 5e-15,
            force_absolute_tolerance_n: 1e-13,
            force_relative_tolerance: 1e-10,
        },
    }
}
pub fn source(
    theta: [f64; 2],
    friction: [f64; 2],
    receding: [f64; 2],
    advancing: [f64; 2],
) -> SourceAsymmetricContactLineBridge {
    let (liquid, wetting) = sources(theta, friction, receding, advancing);
    SourceAsymmetricContactLineBridge::new(liquid, wetting, VOLUME, WIDTH, limits(), config())
        .unwrap()
}
pub fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

/// Independent root in circle radius/arc geometry for the ODE reference.
pub fn angles_at(h: f64, spans: [f64; 2]) -> [f64; 2] {
    let d = (spans[1] - spans[0]) / 2.0;
    let chord = h.hypot(d);
    let alpha = (d / h).atan();
    let mut lo = 0.05_f64;
    let mut hi = 1.45_f64;
    for _ in 0..64 {
        let angle = 0.5 * (lo + hi);
        let radius = chord / (2.0 * angle.sin());
        let segment = radius * radius * (2.0 * angle - (2.0 * angle).sin());
        let volume = WIDTH * (h * (spans[0] + spans[1]) / 2.0 - segment);
        if volume > VOLUME {
            lo = angle;
        } else {
            hi = angle;
        }
    }
    let phi = 0.5 * (lo + hi);
    [
        core::f64::consts::FRAC_PI_2 - phi + alpha,
        core::f64::consts::FRAC_PI_2 - phi - alpha,
    ]
}
