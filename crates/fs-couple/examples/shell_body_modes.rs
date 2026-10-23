//! Closed instrument body modes from a shell assembly table: separately sourced shell components (plates, walls, rib strips,
//! posts, bridge) tied by the rigid work-conjugate mortar, one assembled eigensolve, then the enclosed air (collective
//! compression spring plus the neck plug) on those same modes. Components of one plate may be several uniform-thickness
//! zones; a component's group is its name up to the first '-' (`top-z3` is in group `top`).
//!
//! The geometry is data (`frankensim.shell-assembly.v1`, written by `shell_body_modes/gen_violin_shell_assembly.py`); every
//! wood constant comes from a normalized matdb pack named by the component. Nothing is fitted: no frequency, damping or gain
//! is edited here. `drive` records name the node nearest a point on a part and the global direction of a unit force there.
//!
//!   shell_body_modes ASSEMBLY.tsv PACK_DIR [--modes N] [--fmin HZ] [--fmax HZ] [--omega RAD_S] [--drilling-alpha A] [--shapes PATH]
//!                    [--node-shapes PATH]
//!
//! `PACK_DIR/<pack>.fsmatpk` is read for each pack name in the table. `--shapes` writes the displacement of the top and back
//! nodes along their directors; `--node-shapes` writes the three displacement components of every node of every part.
use std::{collections::BTreeMap, error::Error, fs, io::Write, path::Path};

use fs_couple::cavity_modal::{CavityNeckLoss, HelmholtzModalCavity};
use fs_couple::shear_shell::assembly::{
    MaterialShellAssembly, MaterialShellComponent, ShellAssemblyForcePort,
    ShellAssemblyReductionOptions,
};
use fs_couple::shear_shell::pressure::ShellPressureSurfaces;
use fs_couple::shear_shell::{
    ShellProjectionBudget, ShellRegionMaterial, TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES,
    compile_material_shear_shell,
};
use fs_couple::thin_plate::{PlateMaterialModel, PlateRegionMaterial, PlateThicknessConstraint};
use fs_exec::CancelGate;
use fs_matdb::{MaterialCard, MaterialStateId, NormalizedPack, QueryPoint};
use fs_material::gas::{GasSpec, GasState};
use fs_material::plane_stress_viscosity::PLANE_STRESS_VISCOSITY_PROPERTIES as VISCOSITY;
use fs_material::state_point::{
    MaterialPropertySelection as Pick, ResolvedMaterialStatePoint, ScalarAdmissibility as Admit,
    ScalarPropertyRequirement as Requirement, resolve_material_state_point,
};
use fs_plate::PlateRegion;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::junction::{
    EdgeTie, RotationCoupling, TieBudget, TieKind, TieOffset, TieSide,
};
use fs_plate::shell::mitc::{MitcBudget, NodalDirectors, ShellConstraint, ShellMass};
use fs_plate::shell::pressure::{OrientedShellFacet, ShellFace};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use fs_scenario::{HelmholtzCavity, NeckEnd};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

const TEMPERATURE_K: f64 = 293.15;
const SHEAR_CORRECTION: f64 = 5.0 / 6.0;
/// Default drilling penalty of the shell elements; `--drilling-alpha` changes it for a declared sensitivity (the junction stiffness of a fully tied T-junction depends on it).
const DRILLING_ALPHA: f64 = 1e-3;
const MIN_DIRECTOR_DOT: f64 = 0.9;

struct Part {
    name: String,
    pack: String,
    thickness_m: f64,
    axes: Vec<[f64; 3]>,
    nodes: Vec<[f64; 3]>,
    directors: Vec<[f64; 3]>,
    tris: Vec<[usize; 3]>,
    region_of_tri: Vec<usize>,
    exterior: bool,
    cavity: bool,
}

struct Tie {
    master: String,
    slave: String,
    side: i8,
    master_nodes: Vec<usize>,
    slave_nodes: Vec<usize>,
}

struct Drive {
    label: String,
    part: String,
    position: [f64; 3],
    direction: [f64; 3],
}

struct Table {
    meta: Vec<(String, String)>,
    parts: Vec<Part>,
    ties: Vec<Tie>,
    drives: Vec<Drive>,
    cavity_volume_m3: f64,
    neck_radius_m: f64,
    neck_length_m: f64,
}

fn floats(fields: &[&str]) -> Result<Vec<f64>> {
    fields.iter().map(|f| Ok(f.parse::<f64>()?)).collect()
}

fn indices(field: &str) -> Result<Vec<usize>> {
    field
        .split_ascii_whitespace()
        .map(|f| Ok(f.parse::<usize>()?))
        .collect()
}

fn unit(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
}

fn parse(text: &str) -> Result<Table> {
    let mut lines = text.lines();
    if lines.next() != Some("frankensim.shell-assembly.v1") {
        return Err("not a frankensim.shell-assembly.v1 table".into());
    }
    let mut table = Table {
        meta: Vec::new(),
        parts: Vec::new(),
        ties: Vec::new(),
        drives: Vec::new(),
        cavity_volume_m3: 0.0,
        neck_radius_m: 0.0,
        neck_length_m: 0.0,
    };
    let mut lines = lines.peekable();
    while let Some(line) = lines.next() {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "meta" => table.meta.push((f[1].into(), f[2].into())),
            "component" => {
                let n_regions: usize = f[4].parse()?;
                let n_nodes: usize = f[5].parse()?;
                let n_tris: usize = f[6].parse()?;
                let mut part = Part {
                    name: f[1].into(),
                    pack: f[2].into(),
                    thickness_m: f[3].parse()?,
                    axes: Vec::new(),
                    nodes: Vec::new(),
                    directors: Vec::new(),
                    tris: Vec::new(),
                    region_of_tri: Vec::new(),
                    exterior: f[7] == "1",
                    cavity: f[8] == "1",
                };
                for _ in 0..n_regions {
                    let r = floats(
                        &lines
                            .next()
                            .ok_or("missing region")?
                            .split('\t')
                            .skip(1)
                            .collect::<Vec<_>>(),
                    )?;
                    part.axes.push([r[0], r[1], r[2]]);
                }
                for _ in 0..n_nodes {
                    let l = lines.next().ok_or("missing node")?;
                    let r = floats(&l.split('\t').skip(1).collect::<Vec<_>>())?;
                    part.nodes.push([r[0], r[1], r[2]]);
                    part.directors.push(unit([r[3], r[4], r[5]]));
                }
                for _ in 0..n_tris {
                    let l = lines.next().ok_or("missing triangle")?;
                    let r: Vec<usize> = l
                        .split('\t')
                        .skip(1)
                        .map(|x| x.parse())
                        .collect::<std::result::Result<_, _>>()?;
                    part.tris.push([r[0], r[1], r[2]]);
                    part.region_of_tri.push(r[3]);
                }
                table.parts.push(part);
            }
            "tie" => {
                let master = lines.next().ok_or("missing master")?;
                let slave = lines.next().ok_or("missing slave")?;
                table.ties.push(Tie {
                    master: f[1].into(),
                    slave: f[2].into(),
                    side: f[3].parse()?,
                    master_nodes: indices(master.split('\t').nth(1).ok_or("empty master")?)?,
                    slave_nodes: indices(slave.split('\t').nth(1).ok_or("empty slave")?)?,
                });
            }
            "drive" => {
                let v = floats(&f[3..9])?;
                table.drives.push(Drive {
                    label: f[1].into(),
                    part: f[2].into(),
                    position: [v[0], v[1], v[2]],
                    direction: [v[3], v[4], v[5]],
                });
            }
            "cavity" => {
                table.cavity_volume_m3 = f[1].parse()?;
                table.neck_radius_m = f[2].parse()?;
                table.neck_length_m = f[3].parse()?;
            }
            other => return Err(format!("unknown record {other}").into()),
        }
    }
    Ok(table)
}

fn load_pack(dir: &Path, name: &str) -> Result<(NormalizedPack, MaterialCard)> {
    let bytes = fs::read(dir.join(format!("{name}.fsmatpk")))?;
    let pack = NormalizedPack::from_bytes(&bytes)?;
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("source-pack:{}", pack.pack_id()),
            phase: "declared solid model".into(),
            process: "frozen source-pack realization".into(),
            revision: 0,
        },
        pack.claims().clone(),
        vec![],
    )?;
    Ok((pack, card))
}

fn requirements() -> Result<Vec<Requirement>> {
    let mut r = vec![Requirement::try_new(
        "density",
        Density::DIMS,
        Admit::StrictlyPositive,
    )?];
    for (key, dims) in [
        ("young_modulus_1", Pressure::DIMS),
        ("young_modulus_2", Pressure::DIMS),
        ("poisson_ratio_12", Dims::NONE),
        ("shear_modulus_12", Pressure::DIMS),
    ] {
        r.push(Requirement::try_new(key, dims, Admit::Finite)?);
    }
    for key in VISCOSITY {
        r.push(Requirement::try_new(
            key,
            DynViscosity::DIMS,
            Admit::Finite,
        )?);
    }
    for key in ["shear_modulus_13", "shear_modulus_23"] {
        r.push(Requirement::try_new(
            key,
            Pressure::DIMS,
            Admit::StrictlyPositive,
        )?);
    }
    for key in TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES {
        r.push(Requirement::try_new(
            key,
            DynViscosity::DIMS,
            Admit::Finite,
        )?);
    }
    Ok(r)
}

fn component(
    part: &Part,
    material: &ResolvedMaterialStatePoint,
    drilling_alpha: f64,
    gate: &CancelGate,
) -> Result<MaterialShellComponent> {
    let mesh = ShellMesh::new(part.nodes.clone(), part.tris.clone())?;
    let directors = NodalDirectors::new(&mesh, part.directors.clone(), MIN_DIRECTOR_DOT)?;
    let budget = MitcBudget {
        max_nodes: part.nodes.len(),
        max_triangles: part.tris.len(),
    };
    let regions = part
        .axes
        .iter()
        .enumerate()
        .map(|(k, axis)| ShellRegionMaterial {
            region: PlateRegionMaterial {
                region: PlateRegion {
                    name: format!("{}-r{k}", part.name),
                    triangle_indices: (0..part.tris.len())
                        .filter(|&t| part.region_of_tri[t] == k)
                        .collect(),
                },
                material: material.clone(),
                model: PlateMaterialModel::Orthotropic12,
                material_angle_rad: 0.0,
                thickness_constraint: PlateThicknessConstraint::FixedThickness(part.thickness_m),
            },
            shear_correction: SHEAR_CORRECTION,
            axis1_global: *axis,
        })
        .collect();
    let facets = |on: bool| -> Vec<OrientedShellFacet> {
        if on {
            (0..part.tris.len())
                .map(|triangle| OrientedShellFacet {
                    triangle,
                    side: ShellFace::Positive,
                })
                .collect()
        } else {
            Vec::new()
        }
    };
    let surfaces = ShellPressureSurfaces {
        exterior: facets(part.exterior),
        cavity: facets(part.cavity),
    };
    let specimen = compile_material_shear_shell(mesh, regions, budget, gate)?
        .with_pressure_surfaces(surfaces, gate)?;
    Ok(MaterialShellComponent {
        name: part.name.clone(),
        specimen,
        directors,
        constraints: Vec::<ShellConstraint>::new(),
        mass: ShellMass::Consistent,
        drilling_alpha,
        assembly_budget: budget,
    })
}

/// Eigen-decomposition of a small dense symmetric matrix by cyclic Jacobi rotations.
fn jacobi(mut a: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut v = vec![vec![0.0; n]; n];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for _sweep in 0..100 {
        let mut off = 0.0;
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    off += a[i][j] * a[i][j];
                }
            }
        }
        let diag: f64 = (0..n).map(|i| a[i][i] * a[i][i]).sum();
        if off <= 1e-28 * diag.max(1e-300) {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for row in v.iter_mut() {
                    let (vkp, vkq) = (row[p], row[q]);
                    row[p] = c * vkp - s * vkq;
                    row[q] = s * vkp + c * vkq;
                }
            }
        }
    }
    ((0..n).map(|i| a[i][i]).collect(), v)
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        return Err("usage: shell_body_modes ASSEMBLY.tsv PACK_DIR [--modes N] [--fmin HZ] [--fmax HZ] [--omega RAD_S] [--drilling-alpha A] [--shapes PATH] [--node-shapes PATH]".into());
    }
    let (mut n_modes, mut fmin, mut fmax, mut omega) = (
        10usize,
        15.0f64,
        600.0f64,
        2.0 * std::f64::consts::PI * 200.0,
    );
    let mut shapes: Option<String> = None;
    let mut node_shapes: Option<String> = None;
    let mut drilling_alpha = DRILLING_ALPHA;
    let mut it = args[3..].iter();
    while let Some(flag) = it.next() {
        let value = it.next().ok_or("flag without a value")?;
        match flag.as_str() {
            "--modes" => n_modes = value.parse()?,
            "--fmin" => fmin = value.parse()?,
            "--fmax" => fmax = value.parse()?,
            "--omega" => omega = value.parse()?,
            "--drilling-alpha" => drilling_alpha = value.parse()?,
            "--shapes" => shapes = Some(value.clone()),
            "--node-shapes" => node_shapes = Some(value.clone()),
            other => return Err(format!("unknown flag {other}").into()),
        }
    }
    let table = parse(&fs::read_to_string(&args[1])?)?;
    let dir = Path::new(&args[2]);
    for (k, v) in &table.meta {
        println!("meta {k} = {v}");
    }
    if drilling_alpha != DRILLING_ALPHA {
        println!(
            "drilling_alpha = {drilling_alpha:e} (declared sensitivity; the default is {DRILLING_ALPHA:e})"
        );
    }
    let gate = CancelGate::new_clock_free();
    let requested = requirements()?;
    let query = QueryPoint::new()
        .with("T", TEMPERATURE_K)?
        .with("omega", omega)?;
    let mut cards: BTreeMap<String, (NormalizedPack, MaterialCard)> = BTreeMap::new();
    let mut parts = Vec::with_capacity(table.parts.len());
    let mut densities = Vec::with_capacity(table.parts.len());
    for part in &table.parts {
        if !cards.contains_key(&part.pack) {
            cards.insert(part.pack.clone(), load_pack(dir, &part.pack)?);
        }
        let (pack, card) = &cards[&part.pack];
        let material =
            resolve_material_state_point(card, &query, &requested, Pick::SingleClaimOnly)?;
        println!(
            "part {:<28} pack {} ({}) t={:.3} mm nodes={} triangles={}",
            part.name,
            part.pack,
            pack.content_hash(),
            part.thickness_m * 1e3,
            part.nodes.len(),
            part.tris.len()
        );
        densities.push(
            material
                .property("density")
                .ok_or("pack without a unique density")?
                .value_si(),
        );
        parts.push(component(part, &material, drilling_alpha, &gate)?);
    }
    let position = |name: &str| -> Result<usize> {
        table
            .parts
            .iter()
            .position(|p| p.name == name)
            .ok_or_else(|| format!("tie names unknown part {name}").into())
    };
    let mut edge_ties = Vec::with_capacity(table.ties.len());
    for tie in &table.ties {
        edge_ties.push(EdgeTie {
            a: TieSide {
                component: position(&tie.master)?,
                nodes: tie.master_nodes.clone(),
            },
            b: TieSide {
                component: position(&tie.slave)?,
                nodes: tie.slave_nodes.clone(),
            },
            offset: if tie.side == 0 {
                TieOffset::Midsurface
            } else {
                TieOffset::Skin { side: tie.side }
            },
            rotations: RotationCoupling::Full,
            kind: TieKind::Rigid,
        });
    }
    let total_nodes: usize = table.parts.iter().map(|p| p.nodes.len()).sum();
    let total_tris: usize = table.parts.iter().map(|p| p.tris.len()).sum();
    let slave_nodes: usize = table.ties.iter().map(|t| t.slave_nodes.len()).sum();
    let assembly = MaterialShellAssembly {
        components: parts,
        edge_ties,
        face_ties: Vec::new(),
        tie_budget: TieBudget {
            max_segments: 4 * slave_nodes + 64,
            max_slave_dofs: 6 * slave_nodes + 64,
            gap_tolerance_m: 2e-4,
            angle_tolerance_rad: std::f64::consts::PI,
        },
    };
    // Unit forces named by the table's `drive` records: the node of the named part nearest the stated point, along the stated direction.
    if table.drives.is_empty() {
        return Err("the table has no drive records".into());
    }
    let mut drives: Vec<(String, usize, usize, [f64; 3])> = Vec::new();
    for d in &table.drives {
        let part = position(&d.part)?;
        let node = (0..table.parts[part].nodes.len())
            .min_by(|&a, &b| {
                let dist = |n: usize| {
                    let p = table.parts[part].nodes[n];
                    (0..3).map(|k| (p[k] - d.position[k]).powi(2)).sum::<f64>()
                };
                dist(a).total_cmp(&dist(b))
            })
            .ok_or("drive part has no nodes")?;
        drives.push((d.label.clone(), part, node, unit(d.direction)));
    }
    let drive_labels = drives
        .iter()
        .map(|d| d.0.clone())
        .collect::<Vec<_>>()
        .join(" | ");
    let ports: Vec<ShellAssemblyForcePort> = drives
        .iter()
        .map(
            |(_, drive_part, drive_node, direction)| ShellAssemblyForcePort {
                components: table
                    .parts
                    .iter()
                    .enumerate()
                    .map(|(i, p)| {
                        let mut row = vec![0.0; 6 * p.nodes.len()];
                        if i == *drive_part {
                            row[6 * drive_node..6 * drive_node + 3].copy_from_slice(direction);
                        }
                        row
                    })
                    .collect(),
            },
        )
        .collect();
    let two_pi = 2.0 * std::f64::consts::PI;
    let options = ShellAssemblyReductionOptions {
        eigenvalue_window: ((two_pi * fmin).powi(2), (two_pi * fmax).powi(2)),
        n_modes,
        projection_budget: ShellProjectionBudget {
            max_modes: n_modes.max(8),
            max_projection_terms: 100_000_000_000,
            reference_resistance_n_s_m: 1.0,
        },
        max_total_nodes: total_nodes,
        max_total_triangles: total_tris,
    };
    println!(
        "assembly: {} parts, {total_nodes} nodes, {total_tris} triangles, {} ties",
        table.parts.len(),
        table.ties.len()
    );
    let started = std::time::Instant::now();
    let image = assembly
        .reduce(&options, &ports, &gate)
        .map_err(|e| format!("shell assembly refused: {e:?}"))?;
    println!(
        "reduce: {:.1} s; tie receipt {:?}",
        started.elapsed().as_secs_f64(),
        image.tie_receipt()
    );
    let omegas: Vec<f64> = image
        .eigenvalue_intervals()
        .iter()
        .map(|&(lo, hi)| (0.5 * (lo + hi)).sqrt())
        .collect();
    let n = omegas.len();
    let s_cavity = image.cavity_volume_shapes_m2_per_sqrt_kg();
    let s_exterior = image.volume_shapes_m2_per_sqrt_kg();
    // Modal Kelvin-Voigt damping matrix [1/s] of the shell on these modes. The packs set each viscosity to eta x stiffness / omega_ref
    // (omega_ref = the --omega argument, 2 pi 200 rad/s by default), so a constant loss factor eta_eff = 2 zeta is recovered per mode by
    // C_hyst(w) = (omega_ref / w) C_KV: Q = w^2 / (omega_ref c_eff) with c_eff = v^T C v. Solid loss only: no radiation, no sound-hole loss.
    let damping = image.damping().matrix_per_s();
    if damping.len() != n * n {
        return Err("modal damping matrix has an unexpected size".into());
    }
    let q_hysteretic = |v: &[f64], w: f64| -> f64 {
        let mut c_eff = 0.0;
        for i in 0..n {
            for j in 0..n {
                c_eff += v[i] * v[j] * damping[i * n + j];
            }
        }
        if c_eff > 0.0 {
            w * w / (omega * c_eff)
        } else {
            f64::INFINITY
        }
    };
    println!("in-vacuo shell modes (no air):");
    println!(
        "  mode   f_Hz   rel_interval   cavity_shape_m2/sqrt(kg)   exterior_shape_m2/sqrt(kg)   Q_solid"
    );
    for i in 0..n {
        let (lo, hi) = image.eigenvalue_intervals()[i];
        let mut unit = vec![0.0; n];
        unit[i] = 1.0;
        println!(
            "  {:>4} {:>8.2} {:>12.2e} {:>26.5e} {:>26.5e} {:>9.1}",
            i + 1,
            omegas[i] / two_pi,
            (hi - lo) / (0.5 * (hi + lo)),
            s_cavity[i],
            s_exterior[i],
            q_hysteretic(&unit, omegas[i])
        );
    }
    // Enclosed air on the same modes: U = sum(w_i^2 q_i^2)/2 + (sum s_j q_j + a z)^2 k_c / 2, z = sqrt(M_neck) x.
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), TEMPERATURE_K, 101_325.0)?;
    let cavity = HelmholtzModalCavity::new(
        HelmholtzCavity {
            volume_m3: table.cavity_volume_m3,
            neck_radius_m: table.neck_radius_m,
            neck_length_m: table.neck_length_m,
            inner_end: NeckEnd::Flanged,
            outer_end: NeckEnd::Flanged,
        },
        gas,
        CavityNeckLoss::Resistance(0.0),
    )?;
    let area = cavity.neck_area_m2();
    let k_c = cavity.compression_stiffness_n_m() / (area * area);
    let a_neck = cavity.neck_volume_shape_m2_per_sqrt_kg();
    println!(
        "cavity: V = {:.4} L, neck radius {:.3} mm, length {:.3} mm, isolated Helmholtz {:.2} Hz (rigid walls)",
        table.cavity_volume_m3 * 1e3,
        table.neck_radius_m * 1e3,
        table.neck_length_m * 1e3,
        cavity.isolated_omega_rad_s() / two_pi
    );
    // Mean displacement along the nodal directors (outward from the cavity) of the nodes of one group of components, for each shell mode:
    // the phase relation of the plates.
    let group_of = |name: &str| name.split('-').next().unwrap_or("").to_string();
    let group_w = |group: &str, mode: usize| -> f64 {
        let (mut sum, mut count) = (0.0, 0usize);
        for (p, part) in table.parts.iter().enumerate() {
            if group_of(&part.name) != group {
                continue;
            }
            let z = &image.physical_modes()[p][mode];
            for node in 0..part.nodes.len() {
                let d = part.directors[node];
                sum += d[0] * z[6 * node] + d[1] * z[6 * node + 1] + d[2] * z[6 * node + 2];
                count += 1;
            }
        }
        if count == 0 { 0.0 } else { sum / count as f64 }
    };
    let top_w: Vec<f64> = (0..n).map(|m| group_w("top", m)).collect();
    let back_w: Vec<f64> = (0..n).map(|m| group_w("back", m)).collect();
    let mut open_order: Vec<usize> = Vec::new();
    let mut open_vectors: Vec<Vec<f64>> = Vec::new();
    let mut open_values: Vec<f64> = Vec::new();
    for sealed in [false, true] {
        let neck = if sealed { None } else { Some(a_neck) };
        let dim = n + usize::from(neck.is_some());
        let mut w = s_cavity.to_vec();
        if let Some(a) = neck {
            w.push(a);
        }
        let mut k = vec![vec![0.0; dim]; dim];
        for i in 0..dim {
            for j in 0..dim {
                k[i][j] = k_c * w[i] * w[j];
            }
            if i < n {
                k[i][i] += omegas[i] * omegas[i];
            }
        }
        let (values, vectors) = jacobi(k);
        let mut order: Vec<usize> = (0..dim).collect();
        order.sort_by(|&a, &b| values[a].total_cmp(&values[b]));
        if sealed {
            println!("sealed sound hole (collective air spring on the shell modes, no plug):");
        } else {
            println!("open sound hole (shell modes + collective air spring + sound-hole plug):");
        }
        println!(
            "  mode   f_Hz   neck_share   net_cavity_dV   top_mean_w   back_mean_w   Q_solid   Q_rad(baffled)   Q_total   m_eff[kg] at {}",
            drive_labels
        );
        for (rank, &idx) in order.iter().enumerate() {
            let f = values[idx].max(0.0).sqrt() / two_pi;
            let v: Vec<f64> = (0..dim).map(|r| vectors[r][idx]).collect();
            let neck_share = if neck.is_some() { v[n] * v[n] } else { 0.0 };
            let dv: f64 = (0..dim).map(|r| w[r] * v[r]).sum();
            let tw: f64 = (0..n).map(|r| top_w[r] * v[r]).sum();
            let bw: f64 = (0..n).map(|r| back_w[r] * v[r]).sum();
            let m_eff: Vec<f64> = image
                .force_ports()
                .iter()
                .map(|row| {
                    let phi: f64 = (0..n).map(|r| row[r] * v[r]).sum();
                    1.0 / (phi * phi)
                })
                .collect();
            let q = q_hysteretic(&v[..n], f * two_pi);
            // Compact baffled monopole: P = R |U|^2 / 2 with R = rho w^2 / (2 pi c) gives Q_rad = 2 pi c / (rho w u^2) for a mass-normalised
            // mode of net outward volume u; a free-space monopole radiates half the power, so its Q is twice this.
            let u_net: f64 = (0..n).map(|r| s_exterior[r] * v[r]).sum::<f64>()
                + if neck.is_some() { a_neck * v[n] } else { 0.0 };
            let w_mode = f * two_pi;
            let q_rad = if u_net != 0.0 {
                two_pi * gas.sound_speed / (gas.density * w_mode * u_net * u_net)
            } else {
                f64::INFINITY
            };
            let q_tot = 1.0 / (1.0 / q + 1.0 / q_rad);
            let m_eff_text = m_eff
                .iter()
                .map(|m| format!("{m:.4}"))
                .collect::<Vec<_>>()
                .join(" | ");
            println!(
                "  {:>4} {:>8.2} {:>10.3} {:>13.3e} {:>12.3e} {:>13.3e} {:>9.1} {:>14.1} {:>9.1}   {}",
                rank + 1,
                f,
                neck_share,
                dv,
                tw,
                bw,
                q,
                q_rad,
                q_tot,
                m_eff_text
            );
        }
        if !sealed {
            open_order = order;
            open_vectors = vectors;
            open_values = values;
        }
    }
    let (order, vectors) = (open_order, open_vectors);
    // Share of the shell kinetic energy that each group of parts carries in every open-sound-hole mode, with lumped nodal masses
    // (a third of each triangle's area x thickness x density): a mode carried by one small part, such as the bridge, shows here.
    let lumped: Vec<Vec<f64>> = table
        .parts
        .iter()
        .zip(&densities)
        .map(|(part, &rho)| {
            let mut m = vec![0.0; part.nodes.len()];
            for tri in &part.tris {
                let [a, b, c] = tri.map(|i| part.nodes[i]);
                let (e1, e2) = (
                    [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
                    [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
                );
                let cross = [
                    e1[1] * e2[2] - e1[2] * e2[1],
                    e1[2] * e2[0] - e1[0] * e2[2],
                    e1[0] * e2[1] - e1[1] * e2[0],
                ];
                let area = 0.5 * cross.iter().map(|v| v * v).sum::<f64>().sqrt();
                for &i in tri {
                    m[i] += area * part.thickness_m * rho / 3.0;
                }
            }
            m
        })
        .collect();
    let mut groups: Vec<String> = table.parts.iter().map(|p| group_of(&p.name)).collect();
    groups.sort();
    groups.dedup();
    println!(
        "kinetic energy share of the shell parts by group, open sound hole (lumped nodal mass):"
    );
    println!(
        "  mode   f_Hz{}",
        groups
            .iter()
            .map(|g| format!(" {g:>10}"))
            .collect::<String>()
    );
    for (rank, &idx) in order.iter().enumerate() {
        let mut energy = vec![0.0; groups.len()];
        for (p, part) in table.parts.iter().enumerate() {
            let g = groups
                .iter()
                .position(|g| *g == group_of(&part.name))
                .ok_or("part group missing")?;
            for node in 0..part.nodes.len() {
                let u2: f64 = (0..3)
                    .map(|k| {
                        let u: f64 = (0..n)
                            .map(|r| image.physical_modes()[p][r][6 * node + k] * vectors[r][idx])
                            .sum();
                        u * u
                    })
                    .sum();
                energy[g] += lumped[p][node] * u2;
            }
        }
        let total: f64 = energy.iter().sum();
        println!(
            "  {:>4} {:>8.2}{}",
            rank + 1,
            open_values[idx].max(0.0).sqrt() / two_pi,
            energy
                .iter()
                .map(|e| format!(" {:>10.3}", e / total))
                .collect::<String>()
        );
    }
    if let Some(path) = shapes {
        // Displacement along the nodal director of the top and back plate nodes for the first 36 air-coupled modes.
        let mut out = fs::File::create(&path)?;
        let wanted: Vec<usize> = order.iter().copied().take(36).collect();
        write!(out, "part\tnode\tx_m\ty_m\tz_m")?;
        for rank in 0..wanted.len() {
            write!(out, "\tmode{}", rank + 1)?;
        }
        writeln!(out)?;
        for (p, part) in table.parts.iter().enumerate() {
            let group = group_of(&part.name);
            if group != "top" && group != "back" {
                continue;
            }
            for (node, pos) in part.nodes.iter().enumerate() {
                write!(
                    out,
                    "{}\t{node}\t{:.6}\t{:.6}\t{:.6}",
                    part.name, pos[0], pos[1], pos[2]
                )?;
                let d = part.directors[node];
                for &idx in &wanted {
                    let value: f64 = (0..n)
                        .map(|r| {
                            let z = &image.physical_modes()[p][r];
                            (d[0] * z[6 * node] + d[1] * z[6 * node + 1] + d[2] * z[6 * node + 2])
                                * vectors[r][idx]
                        })
                        .sum();
                    write!(out, "\t{value:.6e}")?;
                }
                writeln!(out)?;
            }
        }
        println!("wrote director displacement of the top and back nodes to {path}");
    }
    if let Some(path) = node_shapes {
        // The three global displacement components of every node of every part [m/sqrt(kg): the shell modes are mass-normalized]
        // for the first 36 air-coupled modes, in the order of the open-sound-hole list.
        let mut out = fs::File::create(&path)?;
        let wanted: Vec<usize> = order.iter().copied().take(36).collect();
        write!(out, "part\tnode\tx_m\ty_m\tz_m")?;
        for rank in 0..wanted.len() {
            write!(out, "\tmode{0}_ux\tmode{0}_uy\tmode{0}_uz", rank + 1)?;
        }
        writeln!(out)?;
        for (p, part) in table.parts.iter().enumerate() {
            for (node, pos) in part.nodes.iter().enumerate() {
                write!(
                    out,
                    "{}\t{node}\t{:.6}\t{:.6}\t{:.6}",
                    part.name, pos[0], pos[1], pos[2]
                )?;
                for &idx in &wanted {
                    for k in 0..3 {
                        let value: f64 = (0..n)
                            .map(|r| image.physical_modes()[p][r][6 * node + k] * vectors[r][idx])
                            .sum();
                        write!(out, "\t{value:.6e}")?;
                    }
                }
                writeln!(out)?;
            }
        }
        println!("wrote the displacement of every node to {path}");
    }
    Ok(())
}
