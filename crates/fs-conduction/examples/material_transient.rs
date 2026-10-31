//! Spatial transient heating with real Cp(T), rho(T), k(T) source curves.
//!
//! Compile catalog manifests with xtask matdb-pack, then run, for example:
//! cargo run -p fs-conduction --example material_transient -- \
//!   nasa316 steel.fsmp 0.01 8 120 120
//! cargo run -p fs-conduction --example material_transient -- \
//!   aa2024 aluminum.fsmp 0.01 8 120 120 nasa316 steel.fsmp 0.01 0.0001
//! cargo run -p fs-conduction --example material_transient -- \
//!   source crates/fs-conduction/examples/material_transient/nasa316.source 0.01 8 120 120
//!
//! Arguments: profile pack thickness_m cells_per_layer duration_s steps,
//! optionally followed by profile2 pack2 thickness2_m contact_R_m2_K_per_W.
//! `source descriptor` replaces either profile/pack pair with explicit material
//! data, including arbitrary property names, typed source axes, and a scalar or
//! three-principal conductivity selection. Run `--help` for descriptor syntax.
//! Profiles select explicit source keys/conditions, not representative properties.
//! nasapeek selects the separately declared NASA PEEK linearized thermal model,
//! including its constant-density assumption and unknown density test temperature.
//! The stated experiment has a 20x20 mm section, initial/reference 300 K,
//! left/right reservoirs at 360/300 K and 1000 W/(m2 K) end-film coefficients.
//! Sides are adiabatic. Contact R is a caller-declared scenario assumption.
//! Geometry and reference mass stay fixed; this does not model expansion,
//! convection inside a fluid, phase change, or first-heating microstructure.
//! Output includes every accepted endpoint; a failed endpoint is never emitted.

use std::collections::BTreeMap;
use std::error::Error;
use std::io::Read;
use std::path::Path;

use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyBackwardEuler;
use fs_conduction::transient::enthalpy::material::{
    MaterialEnthalpySource, MaterialEnthalpyTransport,
};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    ConductionMesh, ConductionProblem, ElementMaterials, InterfaceFacePair, InterfaceResistance,
    InterfaceSurface, MaterialId, MaterialTable, ScalarField, ThermalBc, ThermalBoundaryBuilder,
    ThermalInterfaces,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    NormalizedMaterialCardPack, NormalizedPack, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use fs_qty::{
    Density, Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;

#[path = "material_transient/source.rs"]
mod source;

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const LOW_K: f64 = 300.0;
const HIGH_K: f64 = 360.0;
const WIDTH_M: f64 = 0.02;
const END_HTC: f64 = 1000.0;

fn positive(text: &str, name: &str) -> Result<f64> {
    let value: f64 = text.parse()?;
    if !value.is_finite() || value <= 0.0 {
        return Err(format!("{name} must be finite and positive").into());
    }
    Ok(value)
}

fn profile(name: &str) -> Result<(&'static str, &'static [(&'static str, f64)])> {
    match name {
        "nasa316" => Ok((
            "stainless-316-nasa-tp216435",
            &[
                ("source-generic-316-identity", 1.0),
                ("source-uns-designation-known", 0.0),
                ("source-heat-treatment-known", 0.0),
                ("source-product-form-known", 0.0),
                ("source-pressure-known", 0.0),
            ],
        )),
        "aa2024" => Ok((
            "aluminum-2024-t3-first-heating-km2024",
            &[("first-heating", 1.0), ("source-pressure-known", 0.0)],
        )),
        "nasapeek" => Ok((
            "peek-nasa-thermic-plate-linearized-model",
            &[
                ("source-pressure-atmospheric", 1.0),
                ("source-constant-density-model", 1.0),
                ("source-test-temperature-known", 0.0),
            ],
        )),
        _ => Err("profile must be nasa316, aa2024 or nasapeek".into()),
    }
}

fn query(temperature_k: f64, flags: &[(&str, f64)]) -> Result<QueryPoint> {
    let mut point = QueryPoint::new().with_quantity(
        "temperature",
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        )),
        temperature_k,
    )?;
    for &(name, value) in flags {
        point = point.with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(point)
}

fn load(
    name: &str,
    path: &Path,
) -> Result<(NormalizedMaterialCardPack, MaterialEnthalpyTransport)> {
    if name == "source" {
        return source::load(path);
    }
    let (slug, flags) = profile(name)?;
    let pack = read_pack(path)?;
    if pack.pack_id() != slug {
        return Err(format!(
            "profile {name} requires {slug}, received {}",
            pack.pack_id()
        )
        .into());
    }
    let pack = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: slug.into(),
            phase: "solid".into(),
            process:
                "source-selected stationary thermal specimen; original source conditions retained"
                    .into(),
            revision: 0,
        },
        pack,
    )?;
    let lower = query(LOW_K, flags)?;
    let upper = query(HIGH_K, flags)?;
    let source = MaterialEnthalpyTransport::from_isotropic_card(
        pack.card(),
        MaterialEnthalpySource {
            caloric: HeatCapacitySourceQuery {
                density: &PropertyKey::new("density", Density::DIMS),
                heat_capacity: &PropertyKey::new(
                    "specific-heat-capacity",
                    Dims([2, 0, -2, -1, 0, 0]),
                ),
                temperature_axis: "temperature",
                lower: &lower,
                upper: &upper,
                selection: SelectionPolicy::SingleClaimOnly,
            },
            conductivity: &PropertyKey::with_quantity(
                "thermal-conductivity",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::ThermalConductivity,
                    ValueForm::Static,
                )),
            ),
            reference: &lower,
            phase: SolidLiquidPhase::Solid,
            lower_specific_enthalpy_j_kg: 0.0,
            maximum_interpolation_error_j_kg: 0.01,
        },
    )?;
    Ok((pack, source))
}

fn read_pack(path: &Path) -> Result<NormalizedPack> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)?
        .take(16 * 1024 * 1024 + 1)
        .read_to_end(&mut bytes)?;
    if bytes.len() > 16 * 1024 * 1024 {
        return Err("normalized source pack exceeds 16 MiB".into());
    }
    Ok(NormalizedPack::from_bytes(&bytes)?)
}

fn contact(a: &MaterialCard, b: &MaterialCard, resistance: f64) -> Result<InterfaceResistance> {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
            AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS),
        value: PropertyValue::Scalar {
            value: resistance, dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
        },
        validity: ValidityDomain::unconstrained(),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: Provenance {
            source: "explicit caller-declared contact resistance for this numerical scenario; not a measurement".into(),
            license: "CC0-1.0".into(), artifact: None,
        },
    })?;
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: a.id().clone(),
            texture_frame: "normal-plus-x".into(),
        },
        SurfaceSpec {
            material: b.id().clone(),
            texture_frame: "normal-minus-x".into(),
        },
        SystemContext {
            medium: "caller-declared stationary thermal contact".into(),
            third_body: None,
            environment: "reference thermal scenario".into(),
            history: "fixed resistance; no adhesion or pressure law".into(),
        },
        claims,
        Vec::new(),
    )?;
    Ok(InterfaceResistance::from_card(
        "declared-layer-contact",
        &card,
        &QueryPoint::new(),
        SelectionPolicy::SingleClaimOnly,
    )?)
}

fn mesh(thicknesses: &[f64], cells: usize) -> Result<(ConductionMesh, Vec<usize>, Vec<usize>)> {
    let mut positions = Vec::new();
    let mut tets = Vec::new();
    let mut element_layers = Vec::new();
    let mut vertex_layers = Vec::new();
    let mut left = 0.0;
    for (layer, &thickness) in thicknesses.iter().enumerate() {
        let right = left + thickness;
        let (complex, points) = box_grid([cells, 1, 1], [thickness, WIDTH_M, WIDTH_M]);
        let offset = u32::try_from(positions.len())?;
        element_layers.extend(std::iter::repeat_n(layer, complex.tets.len()));
        vertex_layers.extend(std::iter::repeat_n(layer, points.len()));
        tets.extend(complex.tets.into_iter().map(|tet| tet.map(|v| v + offset)));
        positions.extend(points.into_iter().map(|[x, y, z]| {
            // Exact duplicate trace coordinates for the existing contact owner.
            let x = if x == 0.0 {
                left
            } else if (x - thickness).abs() <= 1e-12 * thickness {
                right
            } else {
                left + x
            };
            [x, y, z]
        }));
        left = right;
    }
    Ok((
        ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions)?,
        element_layers,
        vertex_layers,
    ))
}

fn print_profile(time_s: f64, mesh: &ConductionMesh, layers: &[usize], temperatures: &[f64]) {
    let mut sections = BTreeMap::<(usize, u64), (f64, usize)>::new();
    for ((position, &layer), &t) in mesh.positions().iter().zip(layers).zip(temperatures) {
        let entry = sections.entry((layer, position[0].to_bits())).or_default();
        entry.0 += t;
        entry.1 += 1;
    }
    for ((layer, x), (sum, count)) in sections {
        println!(
            "profile\t{time_s:.12e}\t{layer}\t{:.12e}\t{:.12e}",
            f64::from_bits(x),
            sum / count as f64
        );
    }
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() == 1 && args[0] == "--help" {
        println!("{}", source::USAGE);
        return Ok(());
    }
    if args.len() != 6 && args.len() != 10 {
        return Err(source::USAGE.into());
    }
    let cells: usize = args[3].parse()?;
    let steps: usize = args[5].parse()?;
    let duration = positive(&args[4], "duration")?;
    if !(1..=64).contains(&cells) || !(1..=2000).contains(&steps) || duration > 86400.0 {
        return Err(
            "bounded run requires 1..64 cells/layer, 1..2000 steps, duration <= 86400s".into(),
        );
    }
    let mut thicknesses = vec![positive(&args[2], "thickness")?];
    let mut loaded = vec![load(&args[0], Path::new(&args[1]))?];
    let contact_resistance = if args.len() == 10 {
        thicknesses.push(positive(&args[8], "second thickness")?);
        loaded.push(load(&args[6], Path::new(&args[7]))?);
        Some(positive(&args[9], "contact resistance")?)
    } else {
        None
    };
    if thicknesses.iter().any(|&x| !(1e-5..=1.0).contains(&x)) {
        return Err("layer thickness must be between 1e-5 and 1 metre".into());
    }
    let (mesh, element_layers, vertex_layers) = mesh(&thicknesses, cells)?;
    let length: f64 = thicknesses.iter().sum();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot-reservoir",
            |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::robin(END_HTC, HIGH_K)?,
        )?
        .region(
            "cold-reservoir",
            |face| on_box_face(face.centroid[0], length),
            ThermalBc::robin(END_HTC, LOW_K)?,
        )?
        .adiabatic_remainder()
        .finish()?;
    let interfaces = if let Some(resistance) = contact_resistance {
        let pairs = ThermalInterfaces::coincident_face_pairs(&mesh)?
            .into_iter()
            .map(|pair| {
                if mesh.boundary()[pair.side_a].outward_normal[0] > 0.0 {
                    pair
                } else {
                    InterfaceFacePair {
                        side_a: pair.side_b,
                        side_b: pair.side_a,
                    }
                }
            })
            .collect();
        Some(ThermalInterfaces::new(
            &mesh,
            &boundary,
            vec![InterfaceSurface::new(
                "declared-layer-contact",
                pairs,
                contact(loaded[0].0.card(), loaded[1].0.card(), resistance)?,
            )?],
        )?)
    } else {
        None
    };
    let conductivity = ElementMaterials::new(
        MaterialTable::new(
            loaded
                .iter()
                .enumerate()
                .map(|(i, (_, source))| (MaterialId(i as u32), source.conductivity().clone())),
        )?,
        element_layers
            .iter()
            .map(|&i| MaterialId(i as u32))
            .collect(),
    )?;
    let reference_materials: Vec<_> = loaded
        .iter()
        .map(|(_, source)| source.reference_material())
        .collect();
    let heat = ScalarField::Uniform(0.0);
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| -> Result<()> {
        let cx = Cx::new(&gate, arena, StreamKey {
            seed: 3162024, kernel_id: 3162024, tile: 0, iteration: 0,
        }, Budget::INFINITE, ExecMode::Deterministic);
        let stepper = HeterogeneousEnthalpyBackwardEuler::new(&cx, &mesh, &reference_materials,
            &element_layers, EnthalpyBudget {
                max_vertices: mesh.vertex_count(), max_elements: mesh.element_count(),
            })?;
        let config = EnthalpyStepConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-9, relative_tolerance: 1e-9,
                linear_restart: 32, max_linear_cycles: 16,
                forcing_minimum: 1e-12, forcing_maximum: 1e-3,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 32, energy_tolerance_j: 1e-6,
        };
        let mut h: Vec<_> = vertex_layers.iter().map(|&i|
            loaded[i].1.caloric().curve().knots()[0].specific_enthalpy_j_kg).collect();
        let initial = h.clone();
        println!("frankensim.material-transient.v1");
        println!("scenario\tinitial_K={LOW_K}\treference_K={LOW_K}\tleft_ambient_K={HIGH_K}\tright_ambient_K={LOW_K}\tend_htc_W_m2_K={END_HTC}\twidth_m={WIDTH_M}\tcontact_R_m2_K_W={contact_resistance:?}");
        println!("mesh\tvertices={}\ttetrahedra={}\tcells_per_layer={cells}",
            mesh.vertex_count(), mesh.element_count());
        for (i, (_, source)) in loaded.iter().enumerate() {
            println!("material\t{i}\tcard={}\tpaired={}\trho0_kg_m3={:.12e}\tthickness_m={:.12e}\tCp_rho_receipts={}\tk_receipts={}",
                source.card_identity(), source.identity(), source.reference_density_kg_m3(),
                thicknesses[i], source.caloric().receipts().len(), source.conductivity().receipts().len());
            for (axis, value) in source.reference_query().axes() {
                println!("source-axis\t{i}\t{axis}\t{value:.12e}\t{:?}",
                    source.reference_query().axis_quantities().get(axis));
            }
        }
        println!("step-columns\ttime_s\tstorage_change_J\texternal_input_J\tenergy_residual_J\tmin_K\tmax_K");
        println!("profile-columns\ttime_s\tlayer\tx_m\tmean_temperature_K");
        print_profile(0.0, &mesh, &vertex_layers, &vec![LOW_K; h.len()]);
        let mut external_heat = 0.0;
        for index in 0..steps {
            let t0 = duration * index as f64 / steps as f64;
            let t1 = duration * (index + 1) as f64 / steps as f64;
            let dt = t1 - t0;
            let next = stepper.advance(&cx, ConductionProblem {
                mesh: &mesh, boundary: &boundary, material: loaded[0].1.conductivity(),
                element_materials: Some(&conductivity), source: &heat,
            }, interfaces.as_ref(), &h, dt, config)?;
            let input = dt * (next.source_w + next.dirichlet_in_w
                - next.neumann_out_w - next.robin_out_w);
            external_heat += input;
            let low = next.temperature.iter().copied().fold(f64::INFINITY, f64::min);
            let high = next.temperature.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            println!("step\t{t1:.12e}\t{:.12e}\t{input:.12e}\t{:.12e}\t{low:.12e}\t{high:.12e}",
                next.stored_energy_change_j, next.energy_residual_j);
            for flux in &next.robin_fluxes {
                println!("boundary-out-W\t{t1:.12e}\t{}\t{:.12e}", flux.region, flux.heat_rate_w);
            }
            for flux in &next.contact_fluxes {
                println!("contact-A-to-B-W\t{t1:.12e}\t{}\t{:.12e}\tjump_K={:.12e}",
                    flux.interface, flux.heat_rate_a_to_b_w, flux.mean_jump_k);
            }
            print_profile(t1, &mesh, &vertex_layers, &next.temperature);
            h = next.specific_enthalpy_j_kg;
        }
        let stored: f64 = stepper.reference_nodal_masses_kg().iter().zip(h.iter().zip(initial))
            .map(|(mass, (now, old))| mass * (now - old)).sum();
        println!("total\treference_mass_kg={:.12e}\tstored_change_J={stored:.12e}\texternal_input_J={external_heat:.12e}\tenergy_residual_J={:.12e}",
            stepper.reference_nodal_masses_kg().iter().sum::<f64>(), stored - external_heat);
        Ok(())
    })
}
