//! Synthetic source-bound coupon; these numbers are not measured material data.
use std::sync::Arc;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::source::{
    SnLifeMeasure, SnSourceRequest, SnStressMeasure, resolve_log_log_sn,
};
use fs_material::fatigue::{MeanStressCorrection, RainflowLimits, SnCurve};
use fs_material::state_point::{
    IsotropicElasticStatePoint, MaterialPropertySelection, resolve_isotropic_elastic_state_point,
};
use fs_plate::fatigue::{PlateFatigueConfig, PlateFatigueGauge, PlateFatigueSection};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};

pub const E: f64 = 200e9;
pub const NU: f64 = 0.3;
pub const H: f64 = 0.001;
pub const LIFE_STRESS: f64 = 1e10;

pub struct Fixture {
    pub card: MaterialCard,
    pub point: QueryPoint,
    pub material: Arc<IsotropicElasticStatePoint>,
    pub curve: Arc<SnCurve>,
}
pub fn fixture(young: f64, humidity: f64) -> Fixture {
    let domain = ValidityDomain::unconstrained()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            290.0,
            310.0,
        )
        .with_quantity("humidity", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0);
    let point = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("humidity", QuantitySpec::dimensional(Dims::NONE), humidity)
        .unwrap();
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, 7800.0),
        ("young_modulus", Pressure::DIMS, young),
        ("poisson_ratio", Dims::NONE, NU),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: domain.clone(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic plate fatigue fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let request = SnSourceRequest {
        stress_measure: SnStressMeasure::Amplitude,
        life_measure: SnLifeMeasure::Cycles,
        temperature_axis: "T".into(),
        stress_ratio_axis: "R".into(),
        surface_and_geometry: "synthetic smooth nominal-stress coupon".into(),
        environment: "synthetic fixed environment".into(),
        load_basis: "fully reversed full cycles".into(),
        maximum_points: 2,
    };
    let life = request.life_measure.axis_name();
    claims
        .insert_claim(PropertyClaim {
            key: request.stress_measure.property_key(),
            value: PropertyValue::Curve {
                abscissa: life.into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(100.0, 1e8), (1e14, 1e-4)],
                dims: Pressure::DIMS,
            },
            validity: domain
                .with_quantity("R", QuantitySpec::dimensional(Dims::NONE), -1.0, -1.0)
                .with_quantity(life, QuantitySpec::dimensional(Dims::NONE), 100.0, 1e14),
            interpolation: InterpolationPolicy::TabulatedOnly,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic N = 1e10 / stress amplitude".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic plate alloy".into(),
            phase: "solid".into(),
            process: "synthetic fixed process".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let material = Arc::new(
        resolve_isotropic_elastic_state_point(
            &card,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    );
    let fit = point
        .clone()
        .with_quantity("R", QuantitySpec::dimensional(Dims::NONE), -1.0)
        .unwrap()
        .with_quantity(life, QuantitySpec::dimensional(Dims::NONE), 100.0)
        .unwrap();
    let curve = Arc::new(
        resolve_log_log_sn(
            &card,
            &fit,
            &request,
            ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        )
        .unwrap(),
    );
    Fixture {
        card,
        point,
        material,
        curve,
    }
}
pub fn config(samples: u64) -> PlateFatigueConfig {
    PlateFatigueConfig {
        history: RainflowLimits {
            maximum_samples: samples,
            maximum_residue: 256,
        },
        maximum_elements: 1000,
        maximum_full_dofs: 2000,
        maximum_gauges: 100,
        maximum_total_residue: 25600,
        maximum_slope: 0.05,
        maximum_absolute_strain: 0.01,
    }
}
pub fn section(f: &Fixture) -> PlateFatigueSection {
    PlateFatigueSection {
        material: f.material.clone(),
        thickness_m: H,
    }
}
pub fn gauge(f: &Fixture, element: usize, depth: f64, normal: [f64; 2]) -> PlateFatigueGauge {
    PlateFatigueGauge {
        element,
        surface_point: 0,
        depth_m: depth,
        normal_in_mesh: normal,
        temperature_axis: "T".into(),
        curve: f.curve.clone(),
        mean_stress: MeanStressCorrection::Uncorrected,
    }
}
pub fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-15),
        "{actual:.15e} != {expected:.15e}"
    );
}
