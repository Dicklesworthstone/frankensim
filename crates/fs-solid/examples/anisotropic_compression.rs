//! Compress and release the same coupon in four fixed material orientations.
//!
//! cargo run -p fs-solid --example anisotropic_compression -- 0.004
//!
//! The optional argument is a positive engineering compression magnitude,
//! bounded by 0.02. Output is CSV. All material data are explicit synthetic
//! source claims; they do not calibrate a named alloy. Change their values in
//! `source_card` to compare authored material data through the SAME solver.
//! The source path retains plastic r-values r0=1.8, r45=1.1, r90=2.2 and
//! independent transverse shear-yield ratios R23=0.8 and R31=1.3. The
//! reference yield stress is rolling-direction yield. This symmetric Hill
//! model applies the fitted surface to both tensile and compressive loading.
//!
//! The 100 x 20 x 10 mm coupon has frictionless axial rollers and free lateral
//! faces. P1 tetrahedra solve lateral expansion and oblique plastic shear;
//! removing the end drive solves the zero-force permanent deformation
//! using the retained history.
//! Material orientation changes, while geometry, loads and solver stay fixed.
//! Hill alpha is the work-conjugate accumulation of this anisotropic surface.
//! It is not generally the axial or von Mises equivalent plastic strain.
//! This infinitesimal-strain example does not model finite forging, contact,
//! texture evolution, buckling, necking, damage, or a forming-limit envelope.

#![deny(unsafe_code)]

use std::error::Error;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::plastic::hill::source::planar::{
    HILL_PLANAR_R_PROPERTIES, HILL_TRANSVERSE_SHEAR_RATIO_PROPERTIES,
    resolve_hill48_from_planar_r_values,
};
use fs_material::plastic::hill::source::{
    HILL_HARDENING_MODULUS_PROPERTY, HILL_REFERENCE_YIELD_STRESS_PROPERTY,
};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Dims, Pressure};
use fs_solid::{IncrementSettings, SmallStrainLawField, SmallStrainTetProblem, TetAssemblyBudget};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [0.1, 0.0, 0.0],
    [0.0, 0.02, 0.0],
    [0.1, 0.02, 0.0],
    [0.0, 0.0, 0.01],
    [0.1, 0.0, 0.01],
    [0.0, 0.02, 0.01],
    [0.1, 0.02, 0.01],
];
const TETS: [[usize; 4]; 6] = [
    [0, 1, 3, 7],
    [0, 1, 7, 5],
    [0, 2, 7, 3],
    [0, 2, 6, 7],
    [0, 4, 5, 7],
    [0, 4, 7, 6],
];

fn source_card() -> Result<MaterialCard> {
    let mut properties = vec![
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 70e9),
        (POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.3),
        (HILL_REFERENCE_YIELD_STRESS_PROPERTY, Pressure::DIMS, 100e6),
        (HILL_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 1e9),
    ];
    for (name, value) in HILL_PLANAR_R_PROPERTIES
        .into_iter()
        .zip([1.8, 1.1, 2.2])
        .chain(
            HILL_TRANSVERSE_SHEAR_RATIO_PROPERTIES
                .into_iter()
                .zip([0.8, 1.3]),
        )
    {
        properties.push((name, Dims([0; 6]), value));
    }
    let mut claims = ClaimSet::new();
    for (name, dims, value) in properties {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "explicit synthetic anisotropic compression example".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })?;
    }
    Ok(MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic solid".into(),
            phase: "solid".into(),
            process: "illustrative fixed-axis Hill48 plasticity".into(),
            revision: 0,
        },
        claims,
        vec![],
    )?)
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 41,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn supports(driven: bool) -> Vec<(usize, f64)> {
    let mut fixed = Vec::new();
    for (node, x) in NODES.iter().enumerate() {
        if driven || x[0] == 0.0 {
            fixed.push((3 * node, if driven { x[0] } else { 0.0 }));
        }
        if *x == [0.0, 0.0, 0.0] {
            fixed.extend([(3 * node + 1, 0.0), (3 * node + 2, 0.0)]);
        }
        if *x == [0.0, 0.02, 0.0] {
            fixed.push((3 * node + 2, 0.0));
        }
    }
    fixed
}

fn run(card: &MaterialCard, label: &str, rotation: [[f64; 3]; 3], compression: f64) -> Result<()> {
    let material = resolve_hill48_from_planar_r_values(
        card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        0.05,
        rotation,
    )?;
    let fixed = supports(true);
    let path = [0.0625, 0.125, 0.25, 0.5, 1.0].map(|factor| -compression * factor);
    let problem = SmallStrainTetProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        laws: SmallStrainLawField::Uniform(&material),
        prescribed_m: &fixed,
        nodal_forces_n: &[],
        load_path: &path,
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    };
    let compressed = with_cx(|cx| problem.solve(cx))?;
    let record = compressed
        .increments
        .last()
        .ok_or("empty compression path")?;
    let reaction: f64 = fixed
        .iter()
        .zip(&record.reactions_n)
        .filter(|((dof, _), _)| dof % 3 == 0 && NODES[dof / 3][0] == 0.1)
        .map(|(_, force)| force)
        .sum();
    let released = supports(false);
    let unloading = SmallStrainTetProblem {
        prescribed_m: &released,
        load_path: &[0.0],
        ..problem
    };
    let unloaded = with_cx(|cx| {
        unloading.solve_from(
            compressed.displacement_m.as_flattened(),
            &compressed.element_states,
            cx,
        )
    })?;
    let release_stress = unloaded
        .element_stress
        .iter()
        .flatten()
        .fold(0.0_f64, |largest, stress| largest.max(stress.abs()));
    let strain = compressed.element_strain[0];
    let plastic = &compressed.element_states[0];
    let u = &unloaded.displacement_m;
    // Measure width/thickness at one section, independently of the lateral
    // offset between the end faces caused by plastic shear at 45 degrees.
    let shortening = -(u[7][0] - u[6][0]);
    let width_change = u[7][1] - u[5][1];
    let thickness_change = u[7][2] - u[3][2];
    let end_shear_y = u[5][1] - u[4][1];
    println!(
        "{label},{:.9e},{reaction:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{shortening:.9e},{width_change:.9e},{thickness_change:.9e},{end_shear_y:.9e},{release_stress:.9e}",
        compressed.element_stress[0][0], strain[1], strain[2], strain[3], plastic.alpha,
    );
    Ok(())
}

fn main() -> Result<()> {
    let mut args = std::env::args().skip(1);
    let compression = args
        .next()
        .map_or(Ok(0.004), |value| value.parse::<f64>())?;
    if args.next().is_some() || !compression.is_finite() || compression <= 0.0 || compression > 0.02
    {
        return Err(
            "usage: anisotropic_compression [compression_magnitude], with 0 < magnitude <= 0.02"
                .into(),
        );
    }
    let card = source_card()?;
    println!(
        "orientation,peak_axial_stress_pa,peak_reaction_n,width_strain,thickness_strain,tensor_xy_strain,hill_alpha,permanent_shortening_m,permanent_width_change_m,permanent_thickness_change_m,permanent_end_shear_y_m,released_max_abs_stress_pa"
    );
    let c = std::f64::consts::FRAC_1_SQRT_2;
    for (label, rotation) in [
        (
            "rolling",
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        ),
        (
            "cross",
            [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        ),
        (
            "thickness",
            [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        ),
        (
            "in_plane_45_deg",
            [[c, -c, 0.0], [c, c, 0.0], [0.0, 0.0, 1.0]],
        ),
    ] {
        run(&card, label, rotation, compression)?;
    }
    Ok(())
}
