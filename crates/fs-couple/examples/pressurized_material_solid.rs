//! Pressure-driven finite-deformation motion from a normalized material card.
//!
//! The specimen is a one-metre reference tetrahedron with three free axial
//! coordinates. The nine remaining coordinates are fixed. This small affine
//! coupon demonstrates the native pressure/mass/material path, not a resolved
//! pressure vessel, an incompressibility treatment, or measured validation.
//!
//! ```text
//! cargo run -p fs-couple --example pressurized_material_solid -- --demo 0.8 20 0.01
//! cargo run -p fs-couple --example pressurized_material_solid -- \
//!   --pack MATERIAL.matpack STATE MODEL_INDEX TEMPERATURE_K PRESSURE_PA STEPS DT_S
//! ```
//!
//! Positive pressure is inward on the complete current boundary. Temperature
//! is held; pressure is an external load and does not silently replace the
//! material's source-query pressure coordinate. The card must support the
//! supplied fixed temperature without requiring any further missing axis.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::finite_solid::{FiniteSolidBody, FiniteSolidControls, FiniteSolidMaterial};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::hyper3::pressure::ClosedTetPressureSurface;
use fs_solid::hyper3::{HyperTetProblem, HyperTetSettings, HyperelasticLawField};
use fs_solid::linear3::TetAssemblyBudget;
use std::path::Path;
use std::sync::Arc;

#[path = "../src/bin/music_render/material_packs.rs"]
mod material_packs;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];

fn synthetic_card() -> Result<MaterialCard> {
    let source =
        fs_blake3::hash_bytes(b"authored compressible pressure coupon; not measured material");
    let provenance = Provenance {
        source: "synthetic compressible Ogden pressure coupon, not measured data".into(),
        license: "CC0-1.0".into(),
        artifact: Some(source),
    };
    let domain = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        290.0,
        310.0,
    );
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new("density", Density::DIMS),
        value: PropertyValue::Scalar {
            value: 6.0,
            dims: Density::DIMS,
        },
        validity: domain.clone(),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: provenance.clone(),
    })?;
    let model = ConstitutiveModelCard {
        law: LawId("ogden-decoupled".into()),
        law_version: 1,
        parameters: [
            ("stretch_limit", 3.0, Dims::NONE),
            ("volumetric_kind", 0.0, Dims::NONE),
            ("term_count", 1.0, Dims::NONE),
            ("mu_0", 15.0, Pressure::DIMS),
            ("alpha_0", 2.0, Dims::NONE),
            ("kappa", 40.0, Pressure::DIMS),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: domain,
        sources: vec![source],
        provenance,
    };
    Ok(MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "authored-pressure-coupon".into(),
            revision: 0,
        },
        claims,
        vec![model],
    )?)
}

#[derive(Debug, Clone, PartialEq)]
struct Sample {
    time_s: f64,
    stretch: [f64; 3],
    volume_m3: f64,
    kinetic_j: f64,
    strain_j: f64,
    pressure_work_j: f64,
    integration_defect_j: f64,
}

fn simulate(
    material: &FiniteSolidMaterial,
    pressure_pa: f64,
    steps: usize,
    dt_s: f64,
) -> Result<(f64, Vec<Sample>)> {
    if !pressure_pa.is_finite()
        || !(1..=4096).contains(&steps)
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !(dt_s * steps as f64).is_finite()
    {
        return Err(
            "finite pressure, 1..=4096 steps and a finite positive time interval required".into(),
        );
    }
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 58,
                kernel_id: 58,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let fixed: Vec<_> = (0..12)
            .filter(|dof| ![3, 7, 11].contains(dof))
            .map(|dof| (dof, 0.0))
            .collect();
        let problem = HyperTetProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            material: HyperelasticLawField::Uniform(material),
            prescribed_m: &fixed,
            nodal_forces_n: &[],
            body_force_n_m3: [0.0; 3],
            budget: TetAssemblyBudget::standard(),
            settings: HyperTetSettings::default(),
        };
        let surface = Arc::new(ClosedTetPressureSurface::new(
            &NODES,
            &TETS,
            problem.budget,
            4,
            &cx,
        )?);
        let mut body = FiniteSolidBody::new_with_pressure(
            problem,
            surface.clone(),
            pressure_pa,
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            &cx,
        )?;
        // A declared integration-error gate, not heat, damping or a continuum
        // accuracy certificate. A refusal asks the caller to reduce the step.
        let mut controls = FiniteSolidControls {
            energy_tolerance_j: 1.0e-5,
            maximum_steps: steps,
            ..FiniteSolidControls::default()
        };
        controls.solve.max_newton_iterations = 32;
        let mut samples = Vec::with_capacity(steps);
        for _ in 0..steps {
            let (next, report) =
                body.step_candidate_with_pressure(dt_s, &[], [0.0; 3], pressure_pa, controls, &cx)?;
            let u = next.displacement_m();
            let full: Vec<f64> = u.iter().flatten().copied().collect();
            samples.push(Sample {
                time_s: next.native_state().t,
                stretch: [1.0 + u[1][0], 1.0 + u[2][1], 1.0 + u[3][2]],
                volume_m3: surface.volume_m3(&full, &cx)?,
                kinetic_j: report.kinetic_after_j,
                strain_j: report.strain_after_j,
                pressure_work_j: report.pressure_work_j,
                integration_defect_j: report.integration_defect_j,
            });
            body = next;
        }
        Ok((body.reference_mass_kg(), samples))
    })
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (card, model_index, temperature, pressure, steps, dt) = match args.as_slice() {
        [mode, pressure, steps, dt] if mode == "--demo" =>
            (synthetic_card()?, 0, 300.0, pressure.parse()?, steps.parse()?, dt.parse()?),
        [mode, path, state, model, temperature, pressure, steps, dt] if mode == "--pack" => {
            let path = Path::new(path);
            let parent = path.parent().unwrap_or_else(|| Path::new("."));
            let name = path.file_name().and_then(|name| name.to_str())
                .ok_or("material pack needs a UTF-8 file name")?;
            let pack = material_packs::section_pack_loader(parent)(name)
                .map_err(std::io::Error::other)?;
            if pack.card().id().process != *state {
                return Err("material pack does not have the explicitly requested process state".into());
            }
            (pack.card().clone(), model.parse()?, temperature.parse()?,
                pressure.parse()?, steps.parse()?, dt.parse()?)
        }
        _ => return Err("usage: pressurized_material_solid --demo PRESSURE_PA STEPS DT_S | --pack MATERIAL.matpack STATE MODEL_INDEX TEMPERATURE_K PRESSURE_PA STEPS DT_S".into()),
    };
    let point = QueryPoint::new().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        temperature,
    )?;
    let material = FiniteSolidMaterial::resolve(
        &card,
        model_index,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )?;
    let (mass, samples) = simulate(&material, pressure, steps, dt)?;
    println!(
        "# source_card={}; model_index={model_index}; held_temperature_K={temperature:e}; pressure_Pa={pressure:e}; reference_mass_kg={mass:e}; fixed affine tetrahedral coupon; no measured or mesh-convergence claim",
        card.content_hash().to_hex()
    );
    println!(
        "time_s,stretch_x,stretch_y,stretch_z,volume_m3,kinetic_J,strain_J,step_pressure_work_J,signed_integration_defect_J"
    );
    for s in samples {
        println!(
            "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            s.time_s,
            s.stretch[0],
            s.stretch[1],
            s.stretch[2],
            s.volume_m3,
            s.kinetic_j,
            s.strain_j,
            s.pressure_work_j,
            s.integration_defect_j
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pressure_drives_native_motion_and_replays_from_the_same_source() {
        let card = synthetic_card().unwrap();
        let point = QueryPoint::new()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
            .unwrap();
        let material = FiniteSolidMaterial::resolve(
            &card,
            0,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let (mass, inward) = simulate(&material, 0.8, 3, 0.01).unwrap();
        let (same_mass, outward) = simulate(&material, -0.8, 3, 0.01).unwrap();
        assert_eq!(mass, same_mass);
        assert!(inward[0].volume_m3 < 1.0 / 6.0);
        assert!(outward[0].volume_m3 > 1.0 / 6.0);
        assert!(inward[0].pressure_work_j > 0.0);
        assert!(outward[0].pressure_work_j > 0.0);
        assert_eq!(simulate(&material, 0.8, 3, 0.01).unwrap().1, inward);
    }
}
