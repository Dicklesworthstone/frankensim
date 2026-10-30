#![allow(dead_code)]
use fs_evidence::ValidityDomain;
use fs_ga::{Se3, So3, Twist, Vec3};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fluid_transport::{FluidTransportState, ThermalFluidTransportState};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::*;
use fs_material::thermal::transport::{PreparedThermalTransport, ThermalTransportSourceQuery};
use fs_mbd::{Gravity, RigidBodyIntegrator};
use fs_mbd_material::creeping_flow::CreepingFlowLimits;
use fs_mbd_material::principal::transfer::ConfigurationTransferControls;
use fs_mbd_material::thermal::heating::motion::creeping_flow::HeatedCreepingControls;
use fs_mbd_material::thermal::heating::motion::creeping_flow::heat_exchange::{
    SphereConductionLimits, SphereHeatExchange,
};
use fs_mbd_material::thermal::heating::motion::{HeatedMotionControls, HeatedRigidBody};
use fs_mbd_material::thermal::heating::{HeatedBody, ThermalRegion};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec, Temperature};
use fs_thermal::enthalpy::*;

pub const RADIUS: f64 = 0.1;
pub const DENSITY: f64 = 1000.0;
pub const CP: f64 = 100.0;
pub fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
pub fn point(temperature: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_quantity(), temperature)
        .unwrap()
        .with("pressure", pressure)
        .unwrap()
}
pub fn card(cp: f64, expansion: f64, conductivity: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, DENSITY),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            cp,
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            conductivity,
        ),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            expansion,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, value), (360.0, value)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature_quantity(), 280.0, 360.0)
                    .with("pressure", 1e5, 2e5),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic conductive sphere fixture, not measurements".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
pub fn caloric(material: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        material,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &point(280.0, 1e5),
            upper: &point(360.0, 1e5),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-10,
    )
    .unwrap()
}
pub fn transport(material: &MaterialCard, pressure: f64, upper: f64) -> PreparedThermalTransport {
    PreparedThermalTransport::try_new(
        material,
        ThermalTransportSourceQuery {
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            specific_heat: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            conductivity: &PropertyKey::new(
                THERMAL_CONDUCTIVITY_PROPERTY,
                THERMAL_CONDUCTIVITY_DIMS,
            ),
            temperature_axis: "T",
            lower: &point(280.0, pressure),
            upper: &point(upper, pressure),
            selections: [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 3],
        },
    )
    .unwrap()
}
pub fn reference(material: &MaterialCard, temperature: f64, radius: f64) -> ResolvedBody {
    resolve_body(
        &[RegionBinding {
            key: "sphere",
            material,
            point: &point(temperature, 1e5),
            geometry: VolumeMoments::sphere(radius).unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap()
}
pub fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-8,
    }
}
pub fn body<'a>(
    reference: &'a ResolvedBody,
    material: &'a MaterialCard,
    chart: &'a ResolvedHeatCapacityCurve,
    cp: f64,
    temperature: f64,
    radius: f64,
    dt: f64,
    gravity: Gravity,
    spin: f64,
) -> HeatedRigidBody<'a> {
    body_with_events(
        reference,
        material,
        chart,
        cp,
        temperature,
        radius,
        dt,
        gravity,
        spin,
        vec![],
    )
}

pub fn body_with_events<'a>(
    reference: &'a ResolvedBody,
    material: &'a MaterialCard,
    chart: &'a ResolvedHeatCapacityCurve,
    cp: f64,
    temperature: f64,
    radius: f64,
    dt: f64,
    gravity: Gravity,
    spin: f64,
    events: Vec<TimedEvent>,
) -> HeatedRigidBody<'a> {
    let region = ThermalRegion {
        key: "sphere",
        card: material,
        caloric: chart,
        cell_id: CellId(0),
        material_id: MaterialId(1),
        segment_id: 1,
        body_id: 1,
        geometry: CellGeometry {
            start_m: [0.0, 0.0, -radius],
            end_m: [0.0, 0.0, radius],
            volume_m3: reference.regions()[0].geometry().volume_m3(),
            surface_area_m2: 4.0 * std::f64::consts::PI * radius * radius,
        },
        initial_specific_enthalpy_j_kg: cp * (temperature - 280.0),
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
        expansion_selection: MaterialPropertySelection::SingleClaimOnly,
        fixed_point_m: Vec3::default(),
    };
    let thermal = HeatedBody::new(
        reference,
        vec![region],
        events,
        DepositionBudget {
            max_cells: 1,
            max_links: 1,
            max_events: 8,
            max_substeps: 10000,
        },
        config(dt),
        || Ok(()),
    )
    .unwrap();
    HeatedRigidBody::new(
        thermal,
        Se3::from_parts(So3::identity(), Vec3::default()).unwrap(),
        Twist::new(Vec3::new(0.0, 0.0, spin), Vec3::default()),
        RigidBodyIntegrator::new(gravity),
        HeatedMotionControls {
            maximum_step_s: dt,
            relative_tensor_tolerance: 1e-12,
            transfer: ConfigurationTransferControls {
                relative_momentum_tolerance: 1e-12,
                maximum_parameter_work_j: 10.0,
            },
            maximum_numerical_energy_change_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
}
pub fn fluid(temperature: f64, conductivity: f64, beta: Option<f64>) -> ThermalFluidTransportState {
    ThermalFluidTransportState::try_new(
        FluidTransportState::try_new(0.001, 1.0, conductivity, 0.1).unwrap(),
        temperature,
        beta,
    )
    .unwrap()
}
pub fn flow_limits() -> CreepingFlowLimits {
    CreepingFlowLimits {
        maximum_reynolds: 0.1,
        maximum_diffusion_ratio: 0.1,
    }
}
pub fn flow_controls() -> HeatedCreepingControls {
    HeatedCreepingControls {
        maximum_relaxation_fraction: 0.5,
        maximum_step_energy_defect_j: 1e-11,
    }
}
pub fn exchange(solid: &PreparedThermalTransport) -> SphereHeatExchange<'_> {
    SphereHeatExchange {
        solid,
        limits: SphereConductionLimits {
            maximum_peclet: 0.1,
            maximum_rayleigh: 0.1,
            maximum_biot: 0.1,
            maximum_diffusion_ratio: 0.1,
            maximum_relative_temperature_contrast: 0.1,
        },
    }
}
pub fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}
