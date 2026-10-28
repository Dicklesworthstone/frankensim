//! Run the SAME frame geometry and solvers with different source material data.
//!
//! cargo run -p fs-solid --example material_frame -- 1000 2e9 4000 2e9
//!
//! Arguments are pairs of density [kg/m3] and Young's modulus [Pa]. These are
//! explicit synthetic data, not measured claims for a named metal. Replace the
//! input cards with sourced cards for a real specimen; no solver needs changing.

use std::error::Error;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::{Density, Pressure};
use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::{Frame2d, Geometry, Rayleigh};

fn response(density: f64, young: f64) -> Result<[f64; 4], Box<dyn Error>> {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [
        ("density", density, Density::DIMS),
        ("young_modulus", young, Pressure::DIMS),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 300.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "explicit synthetic frame material".into(),
                license: "CC0-1.0".into(),
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
            chemistry: "synthetic material".into(),
            phase: "solid".into(),
            process: "illustrative isotropic elastic member".into(),
            revision: 0,
        },
        claims,
        vec![],
    )?;
    let point = QueryPoint::new().with("T", 300.0)?;
    let material = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )?;
    let section = ConcentricElasticSection::try_new(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: 0.02,
            material,
        }],
    )?;
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    let base = frame.add_node(0.0, 0.0);
    frame.fix(base, [true; 3]);
    let mut previous = base;
    for i in 1..=4 {
        let node = frame.add_node(i as f64 / 4.0, 0.0);
        frame.add_material_element(previous, node, 5, section.clone())?;
        previous = node;
    }
    let first_omega = frame.natural_frequencies()?[0];
    let acceleration: Vec<_> = (0..frame.ndof())
        .map(|d| if d % 3 == 1 { -9.81 } else { 0.0 })
        .collect();
    let gravity = frame.inertial_forces(&acceleration)?;
    let weight: f64 = -gravity.iter().skip(1).step_by(3).sum::<f64>();
    let mut tip_load = vec![0.0; frame.ndof()];
    tip_load[3 * previous + 1] = 0.01;
    let mut static_frame = frame.clone();
    static_frame.static_load(&tip_load, 1)?;
    let history = frame.newmark_with_load(
        &[0.0; 201],
        1e-5,
        (0.0, 0.0),
        Rayleigh { a0: 0.0, a1: 0.0 },
        &tip_load,
    )?;
    Ok([
        weight,
        static_frame.displacements()[3 * previous + 1],
        first_omega / std::f64::consts::TAU,
        history.u.last().ok_or("empty time history")?[3 * previous + 1],
    ])
}

fn main() -> Result<(), Box<dyn Error>> {
    let values: Vec<f64> = std::env::args()
        .skip(1)
        .map(|s| s.parse::<f64>())
        .collect::<Result<_, _>>()?;
    let values = if values.is_empty() {
        vec![1000.0, 2e9, 4000.0, 2e9, 1000.0, 8e9]
    } else {
        values
    };
    if values.len() % 2 != 0 || values.len() > 32 {
        return Err(
            "usage: material_frame [density_kg_m3 young_modulus_pa]... (at most 16 specimens)"
                .into(),
        );
    }
    println!(
        "Synthetic circular cantilevers: L=1m, radius=0.02m, tip force=0.01N; unchanged solver and geometry."
    );
    println!("density_kg_m3,young_pa,weight_n,static_tip_m,first_mode_hz,tip_at_0_002s_m");
    for pair in values.chunks_exact(2) {
        let [weight, displacement, frequency, transient] = response(pair[0], pair[1])?;
        println!(
            "{:.9e},{:.9e},{weight:.9e},{displacement:.9e},{frequency:.9e},{transient:.9e}",
            pair[0], pair[1]
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn runnable_material_swap_changes_weight_compliance_frequency_and_motion() {
        let base = response(1000.0, 2e9).unwrap();
        let dense = response(4000.0, 2e9).unwrap();
        assert!((dense[0] / base[0] - 4.0).abs() < 1e-12);
        assert_eq!(dense[1], base[1]);
        assert!((dense[2] / base[2] - 0.5).abs() < 1e-9);
        assert!((dense[3] - base[3]).abs() > 0.1 * base[3].abs());
        assert!(response(0.0, 2e9).is_err());
    }
}
