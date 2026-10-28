//! Source-pack-driven, small-strain mechanical coupon.
//!
//! Compile a committed manifest with `xtask matdb-pack`, then run:
//! ```text
//! cargo run -p fs-solid --example material_coupon -- --template phosphor-bronze-c52100 > coupon.tsv
//! cargo run -p fs-solid --example material_coupon -- material.fsmatpk coupon.tsv
//! ```
//! Templates also exist for `aluminum-5005` and `loblolly-pine`. They contain
//! scenario choices and source-coordinate mappings; every material value comes
//! from the supplied normalized pack. Edit the explicit property mappings to
//! use another supported pack. All numerical inputs are coherent SI, including
//! moisture fractions (0.12 means 12 percent), not source-table display units.
//!
//! The TSV starts with `frankensim.material-coupon.v1`. Required records are
//! `pack ID`, `mode NAME`, `geometry-m L Y Z`, `cells N`, `traction-pa S`,
//! and `max-strain B`, separated by tabs. Repeat `property ROLE EXACT_NAME`
//! and `axis EXACT_NAME VALUE_SI` as needed. Modes are `young-poisson`,
//! `young-shear`, `bulk-shear`, and `orthotropic-reference-ratios`.
//! The first three require density and their named pair. The last explicitly
//! selects the existing complete FPL L/R/T reference-ratio schema and requires
//! `orientation` followed by nine row-major entries of the proper rotation
//! whose columns are L/R/T in coupon xyz coordinates. No plane-only fallback.
//!
//! Loading is uniform axial traction on x=L, with x-normal rollers on x=0,
//! a yz pin at the origin, and a y roller at (0,0,Z). Other faces are free.
//! The mesh has N axial cells, one cell across each transverse direction, and
//! six conforming linear tetrahedra per cell. Aligned homogeneous uniaxial
//! response is affine; general rotated anisotropy is a coarse FE calculation.
//! `load-displacement-product-j` reports the final load pairing F dot u,
//! not work integrated over a loading history. For zero initial strain in
//! this linear model it equals twice the stored elastic energy.
//! The selected linear isotropy/orthotropy, zero initial strain and strain
//! ceiling are explicit scenario assumptions, not new measured source claims.
//! This is a fixed-state static comparison, not fracture, creep, moisture
//! evolution, a structural design allowable, or a mesh-error certificate.

#![deny(unsafe_code)]

use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::io::Read;
use std::path::Path;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    MaterialCard, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack, QueryPoint,
};
use fs_material::state_point::{
    IsotropicElasticParameters, MaterialPropertySelection, OrthotropicRequirementSet,
    PoissonReciprocityRule, ReferenceRatioOrthotropicStatePoint, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, is_fixed_external_coordinate,
    resolve_isotropic_elastic_state_point_from_parameters,
    resolve_orthotropic_elastic_state_point_from_reference_ratios,
};
use fs_qty::QuantitySpec;
use fs_solid::linear3::eigenstrain::{TetEigenstrainField, TetEigenstrainState};
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
    TetStaticSolveConfig,
};
use fs_solid::loads3::{SurfaceLoad3, TetLoadGeometry3, VolumeLoad3};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const UNIT_ROTATION: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const WOOD_PROPERTIES: [&str; 13] = [
    "density",
    "young_modulus_longitudinal",
    "er_over_el",
    "et_over_el",
    "glr_over_el",
    "grt_over_el",
    "glt_over_el",
    "nu_lr",
    "nu_lt",
    "nu_rt",
    "nu_rl",
    "nu_tl",
    "nu_tr",
];

const ALUMINUM_INPUT: &str = "frankensim.material-coupon.v1
# Geometry, traction, linear isotropy and strain bound are scenario assumptions.
pack\taluminum-5005-austral-wright-reference
mode\tyoung-shear
property\tdensity\tdensity
property\tyoung\tyoung-modulus
property\tshear\tshear-modulus
axis\ttemperature\t293.15
axis\tengineering-reference-model\t1
axis\tsource-temperatures-matched\t0
geometry-m\t0.1\t0.02\t0.02
cells\t4
traction-pa\t1000000
max-strain\t0.001
";
const BRONZE_INPUT: &str = "frankensim.material-coupon.v1
# E/G test temperature, temper and pressure are unknown in the source.
# The explicit isotropic model uses those constants; it does not qualify a temper.
pack\tphosphor-bronze-c52100-cda-reference
mode\tyoung-shear
property\tdensity\tdensity
property\tyoung\tyoung-modulus
property\tshear\tshear-modulus
axis\ttemperature\t293.15
axis\tsource-test-temperature-known\t0
axis\tsource-temper-known\t0
axis\tsource-pressure-known\t0
geometry-m\t0.1\t0.02\t0.02
cells\t4
traction-pa\t1000000
max-strain\t0.001
";
const WOOD_INPUT: &str = "frankensim.material-coupon.v1
# Clear-specimen species averages at 12% moisture, not graded lumber.
# The source does not specify test temperature; no temperature is inserted.
pack\tpine-loblolly-fpl-gtr282
mode\torthotropic-reference-ratios
axis\tmoisture-content\t0.12
axis\tclear-specimen\t1
axis\tsource-test-temperature-known\t0
geometry-m\t0.1\t0.02\t0.02
cells\t4
traction-pa\t1000000
max-strain\t0.001
# Columns are L/R/T in coupon xyz. To put R along x use: 0 1 0 -1 0 0 0 0 1.
orientation\t1\t0\t0\t0\t1\t0\t0\t0\t1
";

fn invalid(message: impl Into<String>) -> Box<dyn Error> {
    std::io::Error::new(std::io::ErrorKind::InvalidInput, message.into()).into()
}
fn number(text: &str) -> Result<f64> {
    let value = text.parse::<f64>()?;
    if !value.is_finite() {
        return Err(invalid(
            "all scenario coordinates and controls must be finite",
        ));
    }
    Ok(value)
}
fn once<T>(slot: &mut Option<T>, value: T, name: &str) -> Result<()> {
    if slot.replace(value).is_some() {
        return Err(invalid(format!("duplicate {name} record")));
    }
    Ok(())
}
fn required<T>(slot: Option<T>, name: &str) -> Result<T> {
    slot.ok_or_else(|| invalid(format!("missing {name} record")))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Mode {
    YoungPoisson,
    YoungShear,
    BulkShear,
    Wood,
}
impl Mode {
    fn parse(name: &str) -> Result<Self> {
        match name {
            "young-poisson" => Ok(Self::YoungPoisson),
            "young-shear" => Ok(Self::YoungShear),
            "bulk-shear" => Ok(Self::BulkShear),
            "orthotropic-reference-ratios" => Ok(Self::Wood),
            _ => Err(invalid(
                "unknown mode; select an explicit supported elastic parameterization",
            )),
        }
    }
    fn roles(self) -> &'static [&'static str] {
        match self {
            Self::YoungPoisson => &["density", "young", "poisson"],
            Self::YoungShear => &["density", "young", "shear"],
            Self::BulkShear => &["density", "bulk", "shear"],
            Self::Wood => &[],
        }
    }
}

#[derive(Debug)]
struct Inputs {
    pack: String,
    mode: Mode,
    properties: BTreeMap<String, String>,
    axes: BTreeMap<String, f64>,
    size_m: [f64; 3],
    cells: usize,
    traction_pa: f64,
    max_strain: f64,
    orientation: Option<[[f64; 3]; 3]>,
}
impl Inputs {
    fn parse(text: &str) -> Result<Self> {
        if text.len() > 65_536 {
            return Err(invalid("coupon input exceeds 64 KiB"));
        }
        let mut rows = text
            .lines()
            .map(str::trim)
            .filter(|s| !s.is_empty() && !s.starts_with('#'));
        if rows.next() != Some("frankensim.material-coupon.v1") {
            return Err(invalid("expected frankensim.material-coupon.v1 header"));
        }
        let (mut pack, mut mode, mut size, mut cells) = (None, None, None, None);
        let (mut traction, mut strain, mut orientation) = (None, None, None);
        let (mut properties, mut axes) = (BTreeMap::new(), BTreeMap::new());
        for line in rows {
            let fields: Vec<_> = line.split('\t').map(str::trim).collect();
            if fields.iter().any(|field| field.is_empty()) {
                return Err(invalid("empty TSV field"));
            }
            match fields.as_slice() {
                ["pack", value] => once(&mut pack, (*value).to_owned(), "pack")?,
                ["mode", value] => once(&mut mode, Mode::parse(value)?, "mode")?,
                ["property", role, name] => {
                    if properties
                        .insert((*role).to_owned(), (*name).to_owned())
                        .is_some()
                    {
                        return Err(invalid(format!("duplicate property role {role}")));
                    }
                }
                ["axis", name, value] => {
                    if axes.insert((*name).to_owned(), number(value)?).is_some() {
                        return Err(invalid(format!("duplicate source axis {name}")));
                    }
                }
                ["geometry-m", x, y, z] => once(
                    &mut size,
                    [number(x)?, number(y)?, number(z)?],
                    "geometry-m",
                )?,
                ["cells", value] => once(&mut cells, value.parse::<usize>()?, "cells")?,
                ["traction-pa", value] => once(&mut traction, number(value)?, "traction-pa")?,
                ["max-strain", value] => once(&mut strain, number(value)?, "max-strain")?,
                ["orientation", a, b, c, d, e, f, g, h, i] => once(
                    &mut orientation,
                    [
                        [number(a)?, number(b)?, number(c)?],
                        [number(d)?, number(e)?, number(f)?],
                        [number(g)?, number(h)?, number(i)?],
                    ],
                    "orientation",
                )?,
                _ => {
                    return Err(invalid(format!(
                        "unknown record or wrong field count: {line}"
                    )));
                }
            }
            if axes.len() > 64 || properties.len() > 10 {
                return Err(invalid("source-coordinate or property-role limit exceeded"));
            }
        }
        let result = Self {
            pack: required(pack, "pack")?,
            mode: required(mode, "mode")?,
            properties,
            axes,
            size_m: required(size, "geometry-m")?,
            cells: required(cells, "cells")?,
            traction_pa: required(traction, "traction-pa")?,
            max_strain: required(strain, "max-strain")?,
            orientation,
        };
        if result.size_m.iter().any(|x| *x <= 0.0)
            || !(1..=64).contains(&result.cells)
            || result.max_strain <= 0.0
            || result.traction_pa == 0.0
        {
            return Err(invalid(
                "positive dimensions/strain bound, 1..64 cells and nonzero traction required",
            ));
        }
        if result.properties.len() != result.mode.roles().len()
            || result
                .mode
                .roles()
                .iter()
                .any(|role| !result.properties.contains_key(*role))
        {
            return Err(invalid(format!(
                "mode requires exactly these property roles: {:?}",
                result.mode.roles()
            )));
        }
        let distinct: BTreeSet<_> = result.properties.values().collect();
        if distinct.len() != result.properties.len() {
            return Err(invalid(
                "independent parameter roles require distinct source properties",
            ));
        }
        if (result.mode == Mode::Wood) != result.orientation.is_some() {
            return Err(invalid(
                "orientation is required only for orthotropic-reference-ratios",
            ));
        }
        Ok(result)
    }

    fn source_names(&self) -> Vec<&str> {
        if self.mode == Mode::Wood {
            WOOD_PROPERTIES.to_vec()
        } else {
            self.properties.values().map(String::as_str).collect()
        }
    }
}

/// Axis names and values are caller inputs. Their exact quantity descriptors
/// come from the selected source properties, never from spelling heuristics.
fn source_point(card: &MaterialCard, inputs: &Inputs) -> Result<QueryPoint> {
    let mut types: BTreeMap<String, Option<QuantitySpec>> = BTreeMap::new();
    for name in inputs.source_names() {
        for (_, claim) in card.claims().claims_for(name) {
            for axis in claim.validity.bounds().keys() {
                let quantity = claim.validity.axis_quantities().get(axis).copied();
                if let Some(previous) = types.insert(axis.clone(), quantity) {
                    if previous != quantity {
                        return Err(invalid(format!(
                            "source claims disagree on the exact axis type for {axis}"
                        )));
                    }
                }
            }
        }
    }
    let mut point = QueryPoint::new();
    for (name, value) in &inputs.axes {
        let quantity = types.get(name).ok_or_else(|| {
            invalid(format!(
                "axis {name} is absent from the selected source properties"
            ))
        })?;
        point = match quantity {
            Some(quantity) => point.with_quantity(name, *quantity, *value)?,
            None => point.with(name, *value)?,
        };
    }
    // Canonical material resolution below checks every required coordinate
    // against each selected claim's complete validity domain.
    Ok(point)
}

struct BoundMaterial {
    material: TetElasticMaterial,
    source: ResolvedMaterialStatePoint,
    elastic_identity: ContentHash,
    orientation_identity: Option<ContentHash>,
}
fn requirement(
    card: &MaterialCard,
    name: &str,
    admissibility: ScalarAdmissibility,
) -> Result<ScalarPropertyRequirement> {
    let claims = card.claims().claims_for(name);
    let (_, first) = claims
        .first()
        .ok_or_else(|| invalid(format!("missing source property {name}")))?;
    if claims.iter().any(|(_, claim)| claim.key != first.key) {
        return Err(invalid(format!(
            "property {name} has ambiguous source contexts"
        )));
    }
    Ok(ScalarPropertyRequirement::try_with_key(
        &first.key,
        admissibility,
    )?)
}
fn bind_material(card: &MaterialCard, inputs: &Inputs) -> Result<BoundMaterial> {
    let point = source_point(card, inputs)?;
    let result = if inputs.mode == Mode::Wood {
        let state = resolve_orthotropic_elastic_state_point_from_reference_ratios(
            card,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
            inputs.max_strain,
            OrthotropicRequirementSet::Full3d,
            PoissonReciprocityRule::MajorRatios,
        )?;
        let ReferenceRatioOrthotropicStatePoint::Full3d(state) = state else {
            return Err(invalid("a complete 3-D reference-ratio law is required"));
        };
        let orientation = inputs
            .orientation
            .ok_or_else(|| invalid("missing L/R/T orientation"))?;
        let mut hash =
            DomainHasher::new("org.frankensim.examples.material-coupon.lrt-orientation.v1");
        for x in orientation.iter().flatten() {
            hash.update(&x.to_bits().to_le_bytes());
        }
        let orientation_identity = hash.finalize();
        let elastic_identity = state
            .derivation()
            .ok_or_else(|| invalid("reference-ratio derivation was not retained"))?
            .identity;
        // Pass the derivation identity, including the selected reciprocity rule,
        // into the physical operator instead of discarding it at the source seam.
        let material = TetElasticMaterial::try_new_oriented_orthotropic(
            state.density_kg_m3(),
            state.law(),
            orientation,
            *elastic_identity.as_bytes(),
            *orientation_identity.as_bytes(),
        )?;
        BoundMaterial {
            material,
            source: state.resolved().clone(),
            elastic_identity,
            orientation_identity: Some(orientation_identity),
        }
    } else {
        let positive = |role: &str| -> Result<ScalarPropertyRequirement> {
            requirement(
                card,
                &inputs.properties[role],
                ScalarAdmissibility::StrictlyPositive,
            )
        };
        let density = positive("density")?;
        let parameters = match inputs.mode {
            Mode::YoungPoisson => IsotropicElasticParameters::YoungPoisson {
                young_modulus: positive("young")?,
                poisson_ratio: requirement(
                    card,
                    &inputs.properties["poisson"],
                    ScalarAdmissibility::Finite,
                )?,
            },
            Mode::YoungShear => IsotropicElasticParameters::YoungShear {
                young_modulus: positive("young")?,
                shear_modulus: positive("shear")?,
            },
            Mode::BulkShear => IsotropicElasticParameters::BulkShear {
                bulk_modulus: positive("bulk")?,
                shear_modulus: positive("shear")?,
            },
            Mode::Wood => return Err(invalid("unreachable isotropic/wood mode mismatch")),
        };
        let state = resolve_isotropic_elastic_state_point_from_parameters(
            card,
            &point,
            &density,
            &parameters,
            MaterialPropertySelection::SingleClaimOnly,
        )?;
        BoundMaterial {
            material: TetElasticMaterial::from_resolved_elastic_state(&state),
            source: state.resolved().clone(),
            elastic_identity: state.identity(),
            orientation_identity: None,
        }
    };
    // Freezing a strain/rate/history response table at one sampled value does
    // not turn it into the constant elastic law used by this static operator.
    for property in result.source.properties() {
        for (axis, _) in property.source_axes() {
            if !is_fixed_external_coordinate(
                axis,
                point.axis_quantities().get(axis).copied(),
                point.axes().get(axis).copied(),
            ) {
                return Err(invalid(format!(
                    "{} varies on {axis}; this coupon requires a fixed external-state elastic law",
                    property.requirement().name()
                )));
            }
        }
    }
    Ok(result)
}

#[derive(Debug)]
struct CouponResult {
    material_identity: ContentHash,
    solution_identity: ContentHash,
    mass_kg: f64,
    mean_tip_m: [f64; 3],
    applied_n: [f64; 3],
    reaction_n: [f64; 3],
    energy_j: f64,
    external_work_j: f64,
    max_strain: f64,
    free_residual_n: f64,
    relative_residual: f64,
    iterations: usize,
}
fn coupon(material: &TetElasticMaterial, inputs: &Inputs, cx: &Cx<'_>) -> Result<CouponResult> {
    let [length, width, depth] = inputs.size_m;
    let mut nodes = Vec::with_capacity(4 * (inputs.cells + 1));
    let mut tets = Vec::with_capacity(6 * inputs.cells);
    for i in 0..=inputs.cells {
        let x = length * (i as f64 / inputs.cells as f64);
        nodes.extend([
            [x, 0.0, 0.0],
            [x, 0.0, depth],
            [x, width, 0.0],
            [x, width, depth],
        ]);
        if i < inputs.cells {
            let offset = 4 * i;
            for tet in [
                [0, 4, 6, 7],
                [0, 6, 2, 7],
                [0, 2, 3, 7],
                [0, 3, 1, 7],
                [0, 1, 5, 7],
                [0, 5, 4, 7],
            ] {
                tets.push(tet.map(|node| node + offset));
            }
        }
    }
    let fixed = [0, 1, 2, 3, 4, 6, 9];
    let budget = TetAssemblyBudget {
        maximum_nodes: 260,
        maximum_tetrahedra: 384,
        maximum_free_dofs: 780,
        minimum_scaled_jacobian: 1.0e-12,
    };
    let problem = TetLinearElasticProblem {
        nodes_m: &nodes,
        tetrahedra: &tets,
        materials: TetMaterialField::Uniform(material),
        fixed_dofs: &fixed,
        budget,
    };
    let geometry = TetLoadGeometry3::new(&nodes, &tets, budget, 4 * tets.len(), cx)?;
    let tip_faces: Vec<_> = geometry
        .faces()
        .iter()
        .enumerate()
        .filter(|(_, face)| face.nodes.iter().all(|&node| node >= 4 * inputs.cells))
        .map(|(i, _)| i)
        .collect();
    let load = geometry.assemble(
        &[SurfaceLoad3::Traction {
            faces: &tip_faces,
            pascals: [inputs.traction_pa, 0.0, 0.0],
        }],
        VolumeLoad3::None,
        &[],
        cx,
    )?;
    let frame =
        DomainHasher::new("org.frankensim.examples.material-coupon.xyz-frame.v1").finalize();
    let zero_source =
        DomainHasher::new("org.frankensim.examples.material-coupon.zero-initial-strain.v1")
            .finalize();
    let initial = TetEigenstrainState::try_from_principal_strains(
        material,
        [0.0; 3],
        UNIT_ROTATION,
        frame,
        frame,
        zero_source,
        inputs.max_strain,
    )?;
    let field = TetEigenstrainField::Uniform(&initial);
    let solution = problem.solve_eigenstrain_displacement_with_nodal_forces(
        field,
        load.nodal_forces_n(),
        TetStaticSolveConfig {
            relative_residual_tolerance: 1.0e-10,
            maximum_iterations: 10_000,
        },
        cx,
    )?;
    let stress = problem.recover_eigenstrain_stress_with_nodal_forces(
        &solution,
        field,
        load.nodal_forces_n(),
        frame,
        inputs.max_strain,
        cx,
    )?;
    let assembly = problem.assemble(cx)?;
    let displacement = solution.displacement_m();
    let mut imbalance = stress.nodal_internal_force_n().to_vec();
    let mut external_work_j = 0.0;
    for &(dof, force) in load.nodal_forces_n() {
        imbalance[dof / 3][dof % 3] -= force;
        external_work_j += force * displacement[dof / 3][dof % 3];
    }
    let mut reaction_n = [0.0; 3];
    for &dof in &fixed {
        reaction_n[dof % 3] += imbalance[dof / 3][dof % 3];
        imbalance[dof / 3][dof % 3] = 0.0;
    }
    let free_residual_n = imbalance.iter().flatten().fold(0.0_f64, |n, f| n.hypot(*f));
    let mut mean_tip_m = [0.0; 3];
    let area: f64 = tip_faces.iter().map(|&i| geometry.faces()[i].area_m2).sum();
    for &i in &tip_faces {
        let face = &geometry.faces()[i];
        for &node in &face.nodes {
            for c in 0..3 {
                mean_tip_m[c] += (face.area_m2 / area / 3.0) * displacement[node][c];
            }
        }
    }
    let result = CouponResult {
        material_identity: ContentHash(initial.elastic_material_state_identity()),
        solution_identity: solution.identity(),
        mass_kg: assembly.total_mass_kg,
        mean_tip_m,
        applied_n: load.resultant([0.0; 3], cx)?.0,
        reaction_n,
        energy_j: stress.elastic_energy_j(),
        external_work_j,
        max_strain: stress
            .elements()
            .iter()
            .flat_map(|s| s.elastic_strain_mandel())
            .map(|s| s.abs())
            .fold(0.0_f64, f64::max),
        free_residual_n,
        relative_residual: solution.true_relative_residual(),
        iterations: solution.iterations(),
    };
    if [
        result.mass_kg,
        result.energy_j,
        result.external_work_j,
        result.max_strain,
        result.free_residual_n,
        result.relative_residual,
    ]
    .into_iter()
    .chain(result.mean_tip_m)
    .chain(result.applied_n)
    .chain(result.reaction_n)
    .any(|v| !v.is_finite())
    {
        return Err(invalid("coupon output arithmetic is not representable"));
    }
    Ok(result)
}

fn read_bounded(path: &Path, limit: usize) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)?
        .take((limit + 1) as u64)
        .read_to_end(&mut bytes)?;
    if bytes.len() > limit {
        return Err(invalid("input file exceeds the example's resource limit"));
    }
    Ok(bytes)
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() == 2 && args[0] == "--template" {
        let template = match args[1].to_str() {
            Some("aluminum-5005") => ALUMINUM_INPUT,
            Some("phosphor-bronze-c52100") => BRONZE_INPUT,
            Some("loblolly-pine") => WOOD_INPUT,
            _ => {
                return Err(invalid(
                    "template must be aluminum-5005, phosphor-bronze-c52100, or loblolly-pine",
                ));
            }
        };
        print!("{template}");
        return Ok(());
    }
    if args.len() != 2 {
        return Err(invalid(
            "usage: material_coupon <normalized-pack.fsmatpk> <coupon.tsv> | --template <aluminum-5005|phosphor-bronze-c52100|loblolly-pine>",
        ));
    }
    let input_bytes = read_bounded(Path::new(&args[1]), 65_536)?;
    let inputs = Inputs::parse(std::str::from_utf8(&input_bytes)?)?;
    let normalized =
        NormalizedPack::from_bytes(&read_bounded(Path::new(&args[0]), 16 * 1024 * 1024)?)?;
    if normalized.pack_id() != inputs.pack {
        return Err(invalid(format!(
            "requested pack {}, received {}",
            inputs.pack,
            normalized.pack_id()
        )));
    }
    let pack_identity = normalized.content_hash();
    let pack = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: inputs.pack.clone(),
            phase: "solid".into(),
            process: "caller-selected small-strain coupon; exact source conditions retained".into(),
            revision: 0,
        },
        normalized,
    )?;
    let bound = bind_material(pack.card(), &inputs)?;
    let gate = CancelGate::new_clock_free();
    let result = ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 813,
                kernel_id: 813,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        coupon(&bound.material, &inputs, &cx)
    })?;
    // Publish only after source admission, solve and strain-bound recovery pass.
    println!("frankensim.material-coupon-result.v1");
    println!("pack\t{}", inputs.pack);
    println!("pack-identity\t{pack_identity}");
    println!("card-identity\t{}", bound.source.card_identity());
    println!("source-state-identity\t{}", bound.source.identity());
    println!("elastic-derivation-identity\t{}", bound.elastic_identity);
    if let Some(orientation) = bound.orientation_identity {
        println!("orientation-identity\t{orientation}");
    }
    println!("material-identity\t{}", result.material_identity);
    println!("solution-identity\t{}", result.solution_identity);
    for (name, value) in bound.source.query_point() {
        println!(
            "source-axis\t{name}\t{value:.17e}\t{:?}",
            bound.source.axis_quantities().get(name)
        );
    }
    for property in bound.source.properties() {
        println!(
            "source-property\t{}\t{:.17e}\t{}",
            property.requirement().name(),
            property.value_si(),
            property.answer().receipt.selected.0
        );
    }
    println!("mass-kg\t{:.17e}", result.mass_kg);
    println!(
        "mean-loaded-face-displacement-m\t{:.17e}\t{:.17e}\t{:.17e}",
        result.mean_tip_m[0], result.mean_tip_m[1], result.mean_tip_m[2]
    );
    println!(
        "applied-force-n\t{:.17e}\t{:.17e}\t{:.17e}",
        result.applied_n[0], result.applied_n[1], result.applied_n[2]
    );
    println!(
        "support-reaction-n\t{:.17e}\t{:.17e}\t{:.17e}",
        result.reaction_n[0], result.reaction_n[1], result.reaction_n[2]
    );
    println!("elastic-energy-j\t{:.17e}", result.energy_j);
    println!(
        "load-displacement-product-j\t{:.17e}",
        result.external_work_j
    );
    println!("maximum-absolute-mandel-strain\t{:.17e}", result.max_strain);
    println!("free-force-residual-n\t{:.17e}", result.free_residual_n);
    println!("true-relative-residual\t{:.17e}", result.relative_residual);
    println!("iterations\t{}", result.iterations);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_native_coupon_matches_uniaxial_mass_force_displacement_and_work() {
        let gate = CancelGate::new_clock_free();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 814,
                    kernel_id: 814,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            // Explicit numerical fixture; no catalog or source-data claim.
            let material = TetElasticMaterial::try_new(2700.0, 120.0e9, 0.25, [41; 32]).unwrap();
            for cells in [1, 4] {
                for sign in [-1.0, 1.0] {
                    let mut input = Inputs::parse(ALUMINUM_INPUT).unwrap();
                    input.cells = cells;
                    input.traction_pa = sign * 1.0e6;
                    let result = coupon(&material, &input, &cx).unwrap();
                    let expected_tip = input.traction_pa * 0.1 / 120.0e9;
                    let expected_force = input.traction_pa * 0.02 * 0.02;
                    assert!((result.mean_tip_m[0] / expected_tip - 1.0).abs() < 2.0e-9);
                    assert!((result.mass_kg / 0.108 - 1.0).abs() < 2.0e-13);
                    assert!((result.applied_n[0] / expected_force - 1.0).abs() < 2.0e-13);
                    assert!((result.reaction_n[0] / -expected_force - 1.0).abs() < 2.0e-9);
                    let expected_energy = 0.5 * expected_force * expected_tip;
                    assert!((result.energy_j / expected_energy - 1.0).abs() < 2.0e-9);
                    assert!(
                        (result.external_work_j / (2.0 * expected_energy) - 1.0).abs() < 2.0e-9
                    );
                    assert!(result.free_residual_n < expected_force.abs() * 2.0e-9);
                    input.max_strain = 1.0e-7;
                    assert!(coupon(&material, &input, &cx).is_err());
                }
            }
        });
    }

    #[test]
    fn g0_source_templates_require_explicit_complete_model_inputs() {
        for template in [ALUMINUM_INPUT, BRONZE_INPUT, WOOD_INPUT] {
            Inputs::parse(template).unwrap();
        }
        assert!(
            Inputs::parse(&ALUMINUM_INPUT.replace("property\tshear\tshear-modulus\n", "")).is_err()
        );
        assert!(
            Inputs::parse(&WOOD_INPUT.replace("orientation\t1\t0\t0\t0\t1\t0\t0\t0\t1\n", ""))
                .is_err()
        );
        assert!(Inputs::parse(&(ALUMINUM_INPUT.to_owned() + "axis\ttemperature\t294\n")).is_err());
        assert!(Inputs::parse(&ALUMINUM_INPUT.replace("cells\t4", "cells\t65")).is_err());
        assert!(
            Inputs::parse(&ALUMINUM_INPUT.replace("traction-pa\t1000000", "traction-pa\tNaN"))
                .is_err()
        );
    }
}
