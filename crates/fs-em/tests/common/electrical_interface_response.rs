//! Manufactured R''(T,|J_n|), shared by the native planar/tetrahedral checks.
use fs_em::conduction::response::contact::ContactResponseControls;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::{
    conductor::ResolvedConductorState,
    electrical_interface::{
        ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
        resolve_electrical_interface_state,
        response::{ElectricalInterfaceCoordinate, ScalarElectricalInterfaceResponse},
    },
    state_point::MaterialPropertySelection,
};
use fs_qty::Temperature;

pub fn controls() -> ContactResponseControls {
    ContactResponseControls {
        absolute_resistance_tolerance_ohm_m2: 1e-11,
        relative_resistance_tolerance: 1e-11,
        voltage_tolerance_v: 1e-10,
        current_density_tolerance_a_m2: 1e-10,
    }
}

pub fn law(
    first: &ResolvedConductorState,
    second: &ResolvedConductorState,
    beta: f64,
) -> ScalarElectricalInterfaceResponse {
    build(first, second, beta, true)
}

pub fn fixed_law(
    first: &ResolvedConductorState,
    second: &ResolvedConductorState,
) -> ScalarElectricalInterfaceResponse {
    build(first, second, 0.0, false)
}

fn build(
    first: &ResolvedConductorState,
    second: &ResolvedConductorState,
    beta: f64,
    electrical: bool,
) -> ScalarElectricalInterfaceResponse {
    let coordinate = ElectricalInterfaceCoordinate::CurrentDensityMagnitude;
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
                ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            ),
            value: if electrical {
                PropertyValue::Surface2D {
                    axes: ["T".into(), "normal_current".into()],
                    axis_dims: [Temperature::DIMS, coordinate.quantity().dims()],
                    samples: [[300., 0.], [400., 0.], [300., 10.], [400., 10.]]
                        .into_iter()
                        .map(|p| (p, 0.5 + 0.002 * (p[0] - 300.) + beta * p[1]))
                        .collect(),
                    triangles: vec![[0, 1, 2], [1, 3, 2]],
                    dims: ELECTRICAL_CONTACT_RESISTANCE_DIMS,
                }
            } else {
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, 0.5), (400.0, 0.7)],
                    dims: ELECTRICAL_CONTACT_RESISTANCE_DIMS,
                }
            },
            validity: if electrical {
                ValidityDomain::unconstrained()
                    .with("T", 300., 400.)
                    .with_quantity("normal_current", coordinate.quantity(), 0., 10.)
            } else {
                ValidityDomain::unconstrained().with("T", 300., 400.)
            },
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured interface law; not empirical material validation".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: first.material().material().clone(),
            texture_frame: "first".into(),
        },
        SurfaceSpec {
            material: second.material().material().clone(),
            texture_frame: "second".into(),
        },
        SystemContext {
            medium: "manufactured film".into(),
            third_body: None,
            environment: "explicit".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new().with("T", 300.).unwrap();
    let query = if electrical {
        query
            .with_quantity("normal_current", coordinate.quantity(), 0.)
            .unwrap()
    } else {
        query
    };
    let source = resolve_electrical_interface_state(
        &card,
        &query,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let axes = if electrical {
        vec![("normal_current", coordinate)]
    } else {
        vec![]
    };
    ScalarElectricalInterfaceResponse::new(&card, source, &axes).unwrap()
}

// Two unit-length/unit-resistivity bulk regions, plus R''=rc0+beta*|J|.
pub fn density(voltage: f64, temperature: f64, beta: f64) -> f64 {
    let linear = 2.5 + 0.002 * (temperature - 300.);
    2. * voltage.abs() / (linear + (linear * linear + 4. * beta * voltage.abs()).sqrt())
}

/// Manufactured R''=|Delta V|, giving |J_n|=1 at every admitted local point.
/// This is a bounded numerical model, not an empirical contact measurement.
pub fn voltage_jump_law(
    first: &ResolvedConductorState,
    second: &ResolvedConductorState,
) -> ScalarElectricalInterfaceResponse {
    let coordinate = ElectricalInterfaceCoordinate::VoltageJumpMagnitude;
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
                ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Curve {
                abscissa: "jump".into(),
                abscissa_dims: coordinate.quantity().dims(),
                knots: vec![(0.5, 0.5), (3.0, 3.0)],
                dims: ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "jump",
                coordinate.quantity(),
                0.5,
                3.0,
            ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured nonuniform interface; not empirical validation".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: first.material().material().clone(),
            texture_frame: "first".into(),
        },
        SurfaceSpec {
            material: second.material().material().clone(),
            texture_frame: "second".into(),
        },
        SystemContext {
            medium: "manufactured film".into(),
            third_body: None,
            environment: "explicit".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new()
        .with_quantity("jump", coordinate.quantity(), 1.5)
        .unwrap();
    let source = resolve_electrical_interface_state(
        &card,
        &query,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ScalarElectricalInterfaceResponse::new(&card, source, &[("jump", coordinate)]).unwrap()
}
