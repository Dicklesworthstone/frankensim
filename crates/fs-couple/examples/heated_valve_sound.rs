//! A finite material heat history changes a 3-D valve and its actual sound.
//!
//! Native input is a normalized material pack; --demo uses explicitly authored
//! integration data. This is a pair of frozen-state acoustic experiments, not
//! continuous audio during heating or a validated cane/instrument prediction.
use std::io::{Read, Write};

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_couple::material_assembly::{
    AcousticMaterialBindings, valve::compile_pressure_loaded_reed_assembly,
};
use fs_couple::pcm_wav::encode_pcm16_wav_interleaved;
use fs_couple::solid_enthalpy::{
    SolidEnthalpy, TetCellMaterial, ThermalExpansionModel, ThermalStiffnessModel,
    valve::HeatedValveReduction,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, NormalizedMaterialCardPack,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    resolve_isotropic_elastic_state_point,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Density, Dims, Pressure, Temperature};
use fs_scenario::{
    AcousticAssembly, AmbientGas, BeatingReed, CylinderSegment, Listener, ViscothermalDuct,
    WaveguideEnd,
};
use fs_solid::linear3::{TetAssemblyBudget, TetStaticSolveConfig};
use fs_solid::loads3::{SurfaceLoad3, TetLoadGeometry3, VolumeLoad3};
use fs_solid::reduce::{OrificeSpec, SingleModeValveAdmission};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, TimedEvent,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const RATE: u32 = 48_000;
const FRAMES: usize = 4096;
const FRAME: ContentHash = ContentHash([71; 32]);
const LENGTH: f64 = 0.010;
const WIDTH: f64 = 0.012;
const THICKNESS: f64 = 0.004;

fn point(t: f64) -> Result<QueryPoint> {
    Ok(QueryPoint::new().with("T", t)?)
}

struct Mesh {
    nodes: Vec<[f64; 3]>,
    tets: Vec<[usize; 4]>,
    fixed: Vec<usize>,
    lip: Vec<usize>,
}
fn mesh() -> Mesh {
    let node = |i: usize, j: usize, k: usize| (i * 3 + j) * 3 + k;
    let mut nodes = Vec::new();
    for i in 0..=2 {
        for j in 0..=2 {
            for k in 0..=2 {
                nodes.push([
                    LENGTH * i as f64 / 2.0,
                    WIDTH * j as f64 / 2.0,
                    THICKNESS * k as f64 / 2.0,
                ]);
            }
        }
    }
    let mut tets = Vec::new();
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                let v = [
                    node(i, j, k),
                    node(i + 1, j, k),
                    node(i, j + 1, k),
                    node(i + 1, j + 1, k),
                    node(i, j, k + 1),
                    node(i + 1, j, k + 1),
                    node(i, j + 1, k + 1),
                    node(i + 1, j + 1, k + 1),
                ];
                for local in [
                    [0, 1, 3, 7],
                    [0, 3, 2, 7],
                    [0, 2, 6, 7],
                    [0, 6, 4, 7],
                    [0, 4, 5, 7],
                    [0, 5, 1, 7],
                ] {
                    tets.push(local.map(|id| v[id]));
                }
            }
        }
    }
    let mut fixed = Vec::new();
    let mut lip = Vec::new();
    for (id, p) in nodes.iter().enumerate() {
        if p[0] == 0.0 {
            fixed.extend([3 * id, 3 * id + 1, 3 * id + 2]);
        }
        if p[2] == 0.0 && p[0] >= LENGTH / 2.0 {
            lip.push(id);
        }
    }
    Mesh {
        nodes,
        tets,
        fixed,
        lip,
    }
}

// Geometric cell surface, including internal faces. Every face is insulated in
// this uniform-heating experiment; no convection/contact coefficient uses it.
fn tetrahedron_area(nodes: &[[f64; 3]], tet: [usize; 4]) -> f64 {
    let mut area = 0.0;
    for face in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
        let [a, b, c] = face.map(|i| nodes[tet[i]]);
        let u: [f64; 3] = core::array::from_fn(|i| b[i] - a[i]);
        let v: [f64; 3] = core::array::from_fn(|i| c[i] - a[i]);
        let cross = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        area += 0.5 * cross[0].hypot(cross[1]).hypot(cross[2]);
    }
    area
}

fn template() -> AcousticAssembly {
    AcousticAssembly {
        ambient: AmbientGas {
            temperature_k: 300.0,
            pressure_pa: 101_325.0,
            relative_humidity: 0.0,
        },
        string: None,
        pluck: None,
        bow: None,
        blow: None,
        duct: Some(ViscothermalDuct {
            segments: vec![CylinderSegment::cylinder(0.007, 0.5)],
            tone_holes: vec![],
            termination: WaveguideEnd::FlangedOpen,
            wall: None,
        }),
        // Only blowing pressure and attack survive source compilation. Every
        // mechanical coefficient below is replaced by the admitted solid valve.
        reed: Some(BeatingReed {
            rest_opening_m: 1.0,
            width_m: 1.0,
            closing_pressure_pa: 1.0,
            blowing_pressure_pa: 100.0,
            attack_s: 0.008,
            mass_kg: 1.0,
            stiffness_n_m: 1.0,
            damping_ratio: 1.0,
        }),
        soundboard: None,
        body_modes: vec![],
        plate: None,
        cavity: None,
        obstacles: vec![],
        contact_texture: None,
        listener: Listener { distance_m: 1.0 },
        sample_rate_hz: RATE,
        duration_s: FRAMES as f64 / f64::from(RATE),
    }
}

struct Comparison {
    cold: Vec<f64>,
    hot: Vec<f64>,
    cold_frequency: f64,
    hot_frequency: f64,
    cold_gap: f64,
    hot_gap: f64,
    cold_mass: f64,
    hot_mass: f64,
    temperature: f64,
    heat_j: f64,
    stored_heat_j: f64,
    residual_j: f64,
    displacement_m: f64,
    cold_identity: ContentHash,
    hot_identity: ContentHash,
}

fn synthesize(
    card: &MaterialCard,
    power_w: f64,
    duration_s: f64,
    gate: &CancelGate,
) -> Result<Comparison> {
    if !power_w.is_finite()
        || power_w < 0.0
        || !duration_s.is_finite()
        || duration_s <= 0.0
        || duration_s > 3600.0
    {
        return Err(
            "finite nonnegative heat power and duration in (0,3600] seconds required".into(),
        );
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| -> Result<_> {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 0x485641,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        cx.checkpoint().map_err(|_| "heated valve cancelled")?;
        let context = point(300.0)?;
        let elastic = resolve_isotropic_elastic_state_point(
            card,
            &context,
            MaterialPropertySelection::SingleClaimOnly,
        )?;
        let caloric = ResolvedHeatCapacityCurve::try_new(
            card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new(
                    "specific_heat_capacity",
                    SPECIFIC_HEAT_CAPACITY_DIMS,
                ),
                density: &PropertyKey::new("density", Density::DIMS),
                temperature_axis: "T",
                lower: &context,
                upper: &point(400.0)?,
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            0.0,
            1e-6,
        )?;
        let mesh = mesh();
        let geometry = TetLoadGeometry3::new(
            &mesh.nodes,
            &mesh.tets,
            TetAssemblyBudget::standard(),
            4 * mesh.tets.len(),
            &cx,
        )?;
        let total_volume: f64 = geometry.element_volumes_m3().iter().sum();
        let density = elastic.density_kg_m3();
        let events = geometry
            .element_volumes_m3()
            .iter()
            .enumerate()
            .map(|(id, &volume)| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(id as u64),
                    segment_id: id as u64,
                    body_id: 1,
                    material_id: MaterialId(0),
                    mass_kg: density * volume,
                    specific_h_j_kg: 0.0,
                    geometry: CellGeometry {
                        start_m: mesh.nodes[mesh.tets[id][0]],
                        end_m: mesh.nodes[mesh.tets[id][3]],
                        volume_m3: volume,
                        surface_area_m2: tetrahedron_area(&mesh.nodes, mesh.tets[id]),
                    },
                    // Uniform volumetric heating on one homogeneous insulated
                    // specimen gives equal T in every cell and zero conduction.
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: power_w * volume / total_volume,
                    },
                }),
            })
            .collect();
        let config = StepConfig {
            max_dt_s: duration_s / 128.0,
            energy_tolerance_j: 1e-8,
        };
        let thermal = DepositionSimulation::new(
            events,
            &[(MaterialId(0), caloric.curve() as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: mesh.tets.len(),
                max_links: 0,
                max_events: mesh.tets.len(),
                max_substeps: 1024,
            },
        )?
        .advance_to(0.0, config, || {
            cx.checkpoint()
                .map_err(|_| fs_thermal::enthalpy::DepositionError::Cancelled)
        })?
        .0;
        let source = SolidEnthalpy::try_new(
            &thermal,
            &mesh.nodes,
            &mesh.tets,
            (0..mesh.tets.len())
                .map(|id| TetCellMaterial {
                    cell: CellId(id as u64),
                    curve: caloric.curve(),
                    point: ThermalMaterialPoint {
                        material: card,
                        context: &context,
                        temperature_axis: "T",
                        liquid_fraction_axis: None,
                    },
                    mechanical_selection: MaterialPropertySelection::SingleClaimOnly,
                    expansion_selection: MaterialPropertySelection::SingleClaimOnly,
                    stiffness_model: ThermalStiffnessModel::Isotropic,
                    expansion_model: ThermalExpansionModel::Isotropic,
                })
                .collect(),
            TetAssemblyBudget::standard(),
            0.02,
            FRAME,
            &cx,
        )?;
        // Unit inward pressure on the reference underside opens the aperture.
        // The original geometric owner computes the signed nodal areas.
        let wetted: Vec<_> = source
            .geometry()
            .faces()
            .iter()
            .enumerate()
            .filter(|(_, face)| {
                face.nodes
                    .iter()
                    .all(|&id| mesh.nodes[id][2] == 0.0 && mesh.nodes[id][0] >= LENGTH / 2.0)
            })
            .map(|(id, _)| id)
            .collect();
        let loading = source.geometry().assemble(
            &[SurfaceLoad3::Pressure {
                faces: &wetted,
                pascals: 1.0,
            }],
            VolumeLoad3::None,
            &[],
            &cx,
        )?;
        let mut pressure = vec![[0.0; 3]; mesh.nodes.len()];
        for &(dof, area) in loading.nodal_forces_n() {
            pressure[dof / 3][dof % 3] = area;
        }
        let reduction = HeatedValveReduction {
            fixed_dofs: &mesh.fixed,
            orifice: OrificeSpec {
                face_nodes: &mesh.lip,
                opening_axis: [0.0, 0.0, 1.0],
                width_axis: [0.0, 1.0, 0.0],
                opposing_plane_offset_m: -0.0004,
            },
            equilibrium_nodal_forces_n: &[],
            nodal_load_per_pa_m2: &pressure,
            window_hz: (1.0, 50_000.0),
            retained_compliance_floor: 0.5,
            loss_factor: 0.1,
            source_id: "heated-valve-sound/coarse-reference-slab/v1",
            static_solve: TetStaticSolveConfig::default(),
            maximum_nodal_displacement_m: 0.0002,
        };
        let cold = source
            .resolve(&thermal, &cx)?
            .reduce_pressure_valve(&reduction, &cx)?;
        let heated =
            source.advance_pressure_valve_to(&thermal, duration_s, config, &reduction, &cx)?;
        let admission = SingleModeValveAdmission {
            minimum_opening_compliance_fraction: 0.9,
            minimum_pressure_compliance_fraction: 0.9,
            maximum_static_opening_relative_error: 0.1,
        };
        let template = template();
        let cold = compile_pressure_loaded_reed_assembly(
            &template,
            cold,
            admission,
            &AcousticMaterialBindings::default(),
        )?;
        let hot = compile_pressure_loaded_reed_assembly(
            &template,
            heated.valve,
            admission,
            &AcousticMaterialBindings::default(),
        )?;
        let cold_card = cold.valve().valve().card();
        let hot_card = hot.valve().valve().card();
        let displacement_m = hot
            .valve()
            .equilibrium()
            .displacement_m()
            .iter()
            .map(|u| u[0].hypot(u[1]).hypot(u[2]))
            .fold(0.0, f64::max);
        cx.checkpoint().map_err(|_| "heated valve cancelled")?;
        let cold_sound = cold.realize()?;
        cx.checkpoint().map_err(|_| "heated valve cancelled")?;
        let hot_sound = hot.realize()?;
        cx.checkpoint().map_err(|_| "heated valve cancelled")?;
        if cold_sound.gas != hot_sound.gas
            || cold_sound.pressure_pa.len() != FRAMES
            || hot_sound.pressure_pa.len() != FRAMES
            || cold_sound
                .pressure_pa
                .iter()
                .chain(&hot_sound.pressure_pa)
                .any(|p| !p.is_finite())
        {
            return Err(
                "finite pressure samples on one unchanged medium and observation clock required"
                    .into(),
            );
        }
        Ok(Comparison {
            cold: cold_sound.pressure_pa,
            hot: hot_sound.pressure_pa,
            cold_frequency: cold_card.modes[0].frequency_hz,
            hot_frequency: hot_card.modes[0].frequency_hz,
            cold_gap: cold_card.rest_gap_m,
            hot_gap: hot_card.rest_gap_m,
            cold_mass: cold_card.total_mass_kg,
            hot_mass: hot_card.total_mass_kg,
            temperature: heated.thermal.temperature_k(CellId(0))?,
            heat_j: heated.thermal_audit.external_heat_j,
            stored_heat_j: heated.thermal_audit.stored_energy_change_j,
            residual_j: heated.thermal_audit.residual_j,
            displacement_m,
            cold_identity: cold.valve().identity(),
            hot_identity: hot.valve().identity(),
        })
    })
}

fn demo(heat_capacity_scale: f64) -> Result<MaterialCard> {
    if !heat_capacity_scale.is_finite() || heat_capacity_scale <= 0.0 {
        return Err("demo heat-capacity scale must be finite and positive".into());
    }
    let mut claims = ClaimSet::new();
    for (name, dims, a, b) in [
        ("density", Density::DIMS, 1050.0, 1050.0),
        ("young_modulus", Pressure::DIMS, 1e6, 0.6e6),
        ("poisson_ratio", Dims::NONE, 0.3, 0.3),
        (
            "linear_thermal_expansion_coefficient",
            INVERSE_TEMPERATURE_DIMS,
            1e-5,
            3e-5,
        ),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            1000.0 * heat_capacity_scale,
            1000.0 * heat_capacity_scale,
        ),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, a), (400.0, b)],
                dims,
            },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured heated-valve integration data, not measurements".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })?;
    }
    Ok(MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "manufactured integration data".into(),
            revision: 0,
        },
        claims,
        vec![],
    )?)
}

fn run(args: &[String]) -> Result<()> {
    if args.len() != 6 {
        return Err("usage: heated_valve_sound --demo CP_SCALE POWER_W HEAT_SECONDS FULL_SCALE_PA OUTPUT_PREFIX | --pack MATERIAL.pack POWER_W HEAT_SECONDS FULL_SCALE_PA OUTPUT_PREFIX".into());
    }
    let card = match args[0].as_str() {
        "--demo" => demo(args[1].parse()?)?,
        "--pack" => {
            let mut bytes = Vec::new();
            std::fs::File::open(&args[1])?
                .take(16 * 1024 * 1024 + 1)
                .read_to_end(&mut bytes)?;
            if bytes.len() > 16 * 1024 * 1024 {
                return Err("material pack exceeds 16 MiB".into());
            }
            NormalizedMaterialCardPack::from_bytes(&bytes)?
                .card()
                .clone()
        }
        _ => return Err("select --demo or --pack explicitly".into()),
    };
    let scale: f64 = args[4].parse()?;
    if !scale.is_finite() || scale <= 0.0 {
        return Err("positive finite PCM full scale in pascals required".into());
    }
    let result = synthesize(
        &card,
        args[2].parse()?,
        args[3].parse()?,
        &CancelGate::new(),
    )?;
    let interleaved: Vec<_> = result
        .cold
        .iter()
        .zip(&result.hot)
        .flat_map(|(&a, &b)| [a, b])
        .collect();
    let (wav, clips) = encode_pcm16_wav_interleaved(&interleaved, RATE, 2, scale)?;
    let csv_path = format!("{}.csv", args[5]);
    let wav_path = format!("{}.wav", args[5]);
    if std::path::Path::new(&csv_path).exists() || std::path::Path::new(&wav_path).exists() {
        return Err("output CSV and WAV paths must both be new".into());
    }
    // All physical admission and both complete sound experiments precede I/O.
    let mut csv = std::io::BufWriter::new(
        std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(csv_path)?,
    );
    writeln!(csv, "time_s,cold_pressure_pa,heated_pressure_pa")?;
    for (i, (&a, &b)) in result.cold.iter().zip(&result.hot).enumerate() {
        writeln!(csv, "{},{a},{b}", i as f64 / f64::from(RATE))?;
    }
    csv.flush()?;
    std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(wav_path)?
        .write_all(&wav)?;
    eprintln!(
        "material_card={}; cold_valve={}; heated_valve={}",
        card.content_hash().to_hex(),
        result.cold_identity.to_hex(),
        result.hot_identity.to_hex()
    );
    eprintln!(
        "T_K={}; heat_J={}; stored_heat_J={}; residual_J={}; maximum_displacement_m={}",
        result.temperature,
        result.heat_j,
        result.stored_heat_j,
        result.residual_j,
        result.displacement_m
    );
    eprintln!(
        "cold/heated: mass_kg={}/{}; first_mode_Hz={}/{}; gap_m={}/{}; stereo left=cold,right=heated; full_scale_Pa={scale}; clips={clips}",
        result.cold_mass,
        result.hot_mass,
        result.cold_frequency,
        result.hot_frequency,
        result.cold_gap,
        result.hot_gap
    );
    Ok(())
}
fn main() -> Result<()> {
    run(&std::env::args().skip(1).collect::<Vec<_>>())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn changing_only_heat_capacity_changes_heating_mechanics_and_actual_pressure() {
        let a = synthesize(&demo(1.0).unwrap(), 1.0, 10.0, &CancelGate::new()).unwrap();
        let b = synthesize(&demo(2.0).unwrap(), 1.0, 10.0, &CancelGate::new()).unwrap();
        let mass = 1050.0 * LENGTH * WIDTH * THICKNESS;
        assert!((a.temperature - (300.0 + 10.0 / (mass * 1000.0))).abs() < 1e-8);
        assert!((b.temperature - (300.0 + 10.0 / (mass * 2000.0))).abs() < 1e-8);
        for result in [&a, &b] {
            assert!((result.cold_mass - mass).abs() < 1e-12);
            assert_eq!(result.cold_mass, result.hot_mass);
            assert!((result.heat_j - 10.0).abs() < 1e-8);
            assert!((result.stored_heat_j - 10.0).abs() < 1e-8);
            assert!(result.residual_j.abs() < 1e-8);
            assert!(result.displacement_m > 0.0);
            assert!(result.hot_frequency < result.cold_frequency);
            assert!(result.hot.iter().any(|p| *p != 0.0));
            assert_ne!(result.cold, result.hot);
        }
        assert_eq!(a.cold, b.cold);
        assert_ne!(a.hot, b.hot);
        let replay = synthesize(&demo(1.0).unwrap(), 1.0, 10.0, &CancelGate::new()).unwrap();
        assert_eq!(a.hot, replay.hot);
        assert_eq!(a.hot_identity, replay.hot_identity);
    }
    #[test]
    fn unsupported_heat_or_cancelled_run_returns_no_sound_comparison() {
        let card = demo(1.0).unwrap();
        assert!(synthesize(&card, 100.0, 10.0, &CancelGate::new()).is_err());
        assert!(synthesize(&card, 1.0, 0.0, &CancelGate::new()).is_err());
        let gate = CancelGate::new();
        gate.request();
        assert!(synthesize(&card, 1.0, 10.0, &gate).is_err());
    }
}
