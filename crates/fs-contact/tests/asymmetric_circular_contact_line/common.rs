#![allow(dead_code)]
//! Reuse original manufactured source claims, never their planar mechanics.
#[path = "../asymmetric_contact_line/common.rs"]
mod materials;
#[path = "../axisymmetric_contact_line/common.rs"]
mod shared;
use fs_contact::capillary_bridge::axisymmetric::asymmetric_contact_line::*;
use fs_contact::capillary_bridge::axisymmetric::*;
use fs_material::liquid::surface::{ResolvedLiquidSurfaceState, ResolvedLiquidWettingState};

pub const VOLUME: f64 = materials::VOLUME;
pub const SIGMA: f64 = materials::SIGMA;
pub const INITIAL_Q: f64 = 2.0;
pub fn config() -> AsymmetricCircularContactLineConfig {
    shared::config()
}
pub fn numerical() -> AxisymmetricNumerics {
    shared::numerical()
}
pub fn original(
    theta: [f64; 2],
    friction: [f64; 2],
    receding: [f64; 2],
    advancing: [f64; 2],
) -> SourceAxisymmetricCapillaryBridge {
    let (liquid, wetting) = materials::sources(theta, friction, receding, advancing);
    from_parts(liquid, wetting).unwrap()
}
fn from_parts(
    liquid: ResolvedLiquidSurfaceState,
    wetting: [ResolvedLiquidWettingState; 2],
) -> Result<SourceAxisymmetricCapillaryBridge, AxisymmetricError> {
    SourceAxisymmetricCapillaryBridge::new(
        liquid,
        wetting,
        VOLUME,
        AxisymmetricBranch {
            minimum_suction_parameter: 0.0,
            maximum_suction_parameter: 32.0,
            minimum_stability_margin: 0.01,
        },
        AxisymmetricLimits {
            minimum_gap_m: 1e-6,
            maximum_gap_m: 0.01,
            maximum_wetted_radius_m: 0.01,
            maximum_suction_pa: 1e7,
            maximum_abs_force_n: 1.0,
            maximum_abs_surface_energy_j: 0.01,
            maximum_abs_step_work_j: 0.001,
            gravity_m_s2: 0.0,
            maximum_bond_number: 0.01,
        },
    )
}

pub fn different_environment() -> Result<SourceAxisymmetricCapillaryBridge, AxisymmetricError> {
    use fs_contact::capillary_bridge::contact_line::{
        ADVANCING_CONTACT_ANGLE_PROPERTY, CONTACT_LINE_FRICTION_PROPERTY,
        RECEDING_CONTACT_ANGLE_PROPERTY, resolve_contact_line_wetting,
    };
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterfaceSystemCard, InterpolationPolicy, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, UncertaintyModel,
    };
    use fs_material::liquid::DYNAMIC_VISCOSITY_DIMS;
    use fs_material::liquid::surface::EQUILIBRIUM_CONTACT_ANGLE_PROPERTY;
    use fs_material::liquid::surface::adhesion::EquilibriumAngleOrigin;
    use fs_material::state_point::MaterialPropertySelection;
    use fs_qty::Dims;
    let (liquid, mut wetting) =
        materials::sources([1.05, 1.15], [2.0, 4.0], [1.0, 1.1], [1.1, 1.2]);
    let mut claims = ClaimSet::new();
    for (name, dims) in [
        (EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, Dims::NONE),
        (ADVANCING_CONTACT_ANGLE_PROPERTY, Dims::NONE),
        (RECEDING_CONTACT_ANGLE_PROPERTY, Dims::NONE),
        (CONTACT_LINE_FRICTION_PROPERTY, DYNAMIC_VISCOSITY_DIMS),
    ] {
        let value = wetting[1].resolved().property(name).unwrap().value_si();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 330.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured unequal environment".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let mut context = wetting[1].context().clone();
    context.environment = "different manufactured vapor environment".into();
    let card = InterfaceSystemCard::assemble(
        wetting[1].surfaces()[0].clone(),
        wetting[1].surfaces()[1].clone(),
        context,
        claims,
        vec![],
    )
    .unwrap();
    wetting[1] = resolve_contact_line_wetting(
        &card,
        &liquid,
        EquilibriumAngleOrigin::SourcedAngle,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    from_parts(liquid, wetting)
}
pub fn source(
    theta: [f64; 2],
    friction: [f64; 2],
    receding: [f64; 2],
    advancing: [f64; 2],
) -> SourceAsymmetricCircularContactLineBridge {
    SourceAsymmetricCircularContactLineBridge::new(
        original(theta, friction, receding, advancing),
        config(),
    )
    .unwrap()
}

/// Manufactured temperature curves retain one claim per original property.
/// These verify source causality and continuation, not material calibration.
pub fn temperature_cards(
    density_factor: f64,
) -> (fs_matdb::MaterialCard, [fs_matdb::InterfaceSystemCard; 2]) {
    use fs_contact::capillary_bridge::contact_line::{
        ADVANCING_CONTACT_ANGLE_PROPERTY, CONTACT_LINE_FRICTION_PROPERTY,
        RECEDING_CONTACT_ANGLE_PROPERTY,
    };
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, PropertyClaim,
        PropertyKey, PropertyValue, Provenance, UncertaintyModel,
    };
    use fs_material::liquid::surface::{
        EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, SURFACE_TENSION_DIMS, SURFACE_TENSION_PROPERTY,
    };
    use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
    use fs_qty::{Density, Dims, Temperature};
    let (liquid, wetting) = materials::sources([1.05, 1.15], [2.0, 4.0], [1.0, 1.1], [1.1, 1.2]);
    let curve = |name: &str, dims: Dims, cold: f64, hot: f64| PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, cold), (320.0, hot)],
            dims,
        },
        validity: ValidityDomain::unconstrained().with("T", 300.0, 320.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: format!("manufactured circular wetting thermal continuation: {name}"),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    };
    let mut claims = ClaimSet::new();
    for (name, dims, cold, hot) in [
        ("density", Density::DIMS, 1000.0, 1000.0 * density_factor),
        (
            DYNAMIC_VISCOSITY_PROPERTY,
            DYNAMIC_VISCOSITY_DIMS,
            1e-3,
            8e-4,
        ),
        (
            SURFACE_TENSION_PROPERTY,
            SURFACE_TENSION_DIMS,
            SIGMA,
            0.8 * SIGMA,
        ),
    ] {
        claims.insert_claim(curve(name, dims, cold, hot)).unwrap();
    }
    let material =
        MaterialCard::assemble(liquid.material().material().clone(), claims, vec![]).unwrap();
    let interfaces = core::array::from_fn(|i| {
        let mut claims = ClaimSet::new();
        for (name, dims, cold, hot) in [
            (
                EQUILIBRIUM_CONTACT_ANGLE_PROPERTY,
                Dims::NONE,
                [1.05, 1.15][i],
                [0.98, 1.12][i],
            ),
            (
                RECEDING_CONTACT_ANGLE_PROPERTY,
                Dims::NONE,
                [1.0, 1.1][i],
                [0.93, 1.07][i],
            ),
            (
                ADVANCING_CONTACT_ANGLE_PROPERTY,
                Dims::NONE,
                [1.1, 1.2][i],
                [1.03, 1.17][i],
            ),
            (
                CONTACT_LINE_FRICTION_PROPERTY,
                DYNAMIC_VISCOSITY_DIMS,
                [2.0, 4.0][i],
                [1.0, 8.0][i],
            ),
        ] {
            claims.insert_claim(curve(name, dims, cold, hot)).unwrap();
        }
        InterfaceSystemCard::assemble(
            wetting[i].surfaces()[0].clone(),
            wetting[i].surfaces()[1].clone(),
            wetting[i].context().clone(),
            claims,
            vec![],
        )
        .unwrap()
    });
    (material, interfaces)
}

pub fn temperature_source(
    material: &fs_matdb::MaterialCard,
    interfaces: &[fs_matdb::InterfaceSystemCard; 2],
    temperature: f64,
) -> SourceAsymmetricCircularContactLineBridge {
    use fs_contact::capillary_bridge::contact_line::resolve_contact_line_wetting;
    use fs_material::liquid::surface::adhesion::EquilibriumAngleOrigin;
    use fs_material::liquid::surface::{LiquidCapillarity, resolve_liquid_surface_state};
    use fs_material::liquid_temperature::LiquidViscosityBasis;
    use fs_material::state_point::MaterialPropertySelection;
    let point = fs_matdb::QueryPoint::new().with("T", temperature).unwrap();
    let liquid = resolve_liquid_surface_state(
        material,
        &point,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let wetting = interfaces.each_ref().map(|interface| {
        resolve_contact_line_wetting(
            interface,
            &liquid,
            EquilibriumAngleOrigin::SourcedAngle,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
    });
    SourceAsymmetricCircularContactLineBridge::new(from_parts(liquid, wetting).unwrap(), config())
        .unwrap()
}

pub fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
