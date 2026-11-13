//! Axial/torsional native rod cycles observed on fixed material planes.
//! All material and S-N values are authored demonstrations, not measured life.
//! Run: cargo run -p fs-solid --example rod_fatigue -- --demo 20

use std::sync::Arc;

use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fatigue::{
    FatigueScope, MeanStressCorrection, NormalStressDirection, RainflowLimits, SnCurve, SnPoint,
};
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::{Density, Pressure, QuantitySpec, Temperature};
use fs_solid::rod::TipLoad;
use fs_solid::rod::material_field::SpatialMaterialRod;
use fs_solid::rod::material_field::fatigue::{
    RodFatigueConfig, RodFatigueGauge, RodFatiguePoint, SpatialRodFatiguePath,
};

fn section(radius_m: f64) -> Result<ConcentricShearSection, Box<dyn std::error::Error>> {
    let temperature = QuantitySpec::dimensional(Temperature::DIMS);
    let domain = ValidityDomain::unconstrained().with_quantity("T", temperature, 300.0, 300.0);
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [
        ("density", 7800.0, Density::DIMS),
        ("young_modulus", 200e9, Pressure::DIMS),
        ("shear_modulus", 200e9 / 2.6, Pressure::DIMS),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: domain.clone(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "authored synthetic rod demonstration".into(),
                license: "demonstration-only".into(),
                artifact: None,
            },
        })?;
        requirements.push(ScalarPropertyRequirement::try_new(
            name,
            dims,
            ScalarAdmissibility::StrictlyPositive,
        )?);
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic elastic metal".into(),
            phase: "solid".into(),
            process: "authored demonstration".into(),
            revision: 0,
        },
        claims,
        vec![],
    )?;
    let query = QueryPoint::new().with_quantity("T", temperature, 300.0)?;
    let material = resolve_material_state_point(
        &card,
        &query,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )?;
    Ok(ConcentricElasticSection::try_new(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: radius_m,
            material,
        }],
    )?
    .with_shear_modulus("shear_modulus", 0.8)?)
}

fn run(cycles: usize) -> Result<(), Box<dyn std::error::Error>> {
    let gate = CancelGate::new();
    let radius_m: f64 = 0.01;
    let source = section(radius_m)?;
    let rod = SpatialMaterialRod::try_nonuniform_straight(vec![0.15, 0.35], vec![source; 2])?;
    let curve = Arc::new(SnCurve::new(
        FatigueScope {
            source_reference: "authored synthetic table; not measured".into(),
            // An S-N scope names the material "chemistry; phase".
            material: "synthetic elastic metal; solid".into(),
            process: "authored demonstration".into(),
            surface_and_geometry: "smooth synthetic circular shaft".into(),
            environment: "synthetic fixed atmosphere".into(),
            load_basis: "signed normal stress amplitude on a fixed material plane".into(),
            temperature_k: 300.0,
        },
        vec![
            SnPoint {
                amplitude_pa: 1.0,
                cycles_to_failure: 1e12,
            },
            SnPoint {
                amplitude_pa: 300e6,
                cycles_to_failure: 1.0,
            },
        ],
        None,
    )?);
    let point = RodFatiguePoint {
        segment: 0,
        layer: 0,
        coordinates_yz_m: [radius_m, 0.0],
    };
    let gauges = [[1.0, 0.0, 0.0], [1.0, 0.0, 1.0]]
        .into_iter()
        .map(|normal| {
            Ok(RodFatigueGauge {
                point,
                normal_in_material_frame: NormalStressDirection::new(normal)?,
                temperature_axis: "T".into(),
                curve: Arc::clone(&curve),
                mean_stress: MeanStressCorrection::Uncorrected,
            })
        })
        .collect::<Result<Vec<_>, fs_material::fatigue::FatigueError>>()?;
    let mut path = SpatialRodFatiguePath::new(
        rod,
        gauges,
        vec![],
        RodFatigueConfig {
            history: RainflowLimits {
                maximum_samples: (2 * cycles + 2) as u64,
                maximum_residue: 64,
            },
            maximum_gauges: 2,
            maximum_total_residue: 128,
            maximum_total_findley_samples: 0,
            maximum_absolute_axial_strain: 0.005,
            maximum_engineering_shear: 0.005,
            maximum_unassessed_transverse_force_n: 0.0,
        },
        &gate,
    )?;
    let area = core::f64::consts::PI * radius_m.powi(2);
    let polar = 0.5 * core::f64::consts::PI * radius_m.powi(4);
    let load = |sign: f64| TipLoad {
        force: [sign * 40e6 * area, 0.0, 0.0],
        moment: [sign * 15e6 * polar / radius_m, 0.0, 0.0],
    };
    path.advance(&[load(1.0)], 1e-4, &gate)?;
    for _ in 0..cycles {
        // A block boundary does not close or discard unresolved rainflow cycles.
        path.advance(&[load(-1.0), load(1.0)], 1e-4, &gate)?;
    }
    let sample = path.observation().normal_samples()[0];
    let committed = path.summaries();
    let preview = path.preview_fatigue(&gate)?;
    assert_eq!(path.summaries(), committed);
    let completed = path.finish(&gate)?;
    assert_eq!(preview, completed);
    println!("synthetic native shaft; no service-life qualification");
    println!("cycles={cycles} accepted_samples={}", path.sample_count());
    println!(
        "local_axial_stress_pa={:.8e} local_torsional_shear_pa={:.8e}",
        sample.stress_pa[0], sample.stress_pa[5]
    );
    for (name, result) in ["axial plane", "oblique plane"].into_iter().zip(completed) {
        println!(
            "{name}: damage_if_stopped_now={:.8e} full_cycles={} half_cycles={}",
            result.miner_damage, result.full_cycles, result.half_cycles
        );
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.first().map(String::as_str) != Some("--demo") || args.len() > 2 {
        return Err(
            "usage: rod_fatigue --demo [cycles: 1..5000]; synthetic calibration only".into(),
        );
    }
    let cycles = args
        .get(1)
        .map(|value| value.parse::<usize>())
        .transpose()?
        .unwrap_or(20);
    if !(1..=5000).contains(&cycles) {
        return Err("cycles must be in 1..=5000".into());
    }
    run(cycles)
}
