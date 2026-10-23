//! Source packs + an authored stationary gas field -> physical DG pressure traces.
//! See gas-field-pulse.md for source preparation, file grammar and interpretation.
#[path = "../src/bin/music_render/json.rs"]
mod json;
#[path = "../src/bin/music_render/thermochemical.rs"]
mod thermochemical;
use json::json_string;

use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, DgWallResolution};
use fs_duct::{Duct, Segment};
use fs_matdb::QueryPoint;
use fs_material::gas_field::{AxialGasField, PrimitiveGasState};
use fs_material::thermochemical_gas::packs::from_packs;
use fs_thermochem::{Composition, CompositionBasis};
use std::{
    error::Error,
    fs::{File, OpenOptions},
    io::{Read, Write},
    path::Path,
    str::{FromStr, SplitAsciiWhitespace},
};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const MAX_INPUT_BYTES: usize = 64 * 1024;
const MAX_OUTPUT_BYTES: usize = 4 * 1024 * 1024;

fn word<'a>(words: &mut SplitAsciiWhitespace<'a>, label: &str) -> Result<&'a str> {
    words
        .next()
        .ok_or_else(|| format!("missing {label}").into())
}
fn tag(words: &mut SplitAsciiWhitespace<'_>, expected: &str) -> Result<()> {
    if word(words, expected)? != expected {
        return Err(format!("expected {expected}").into());
    }
    Ok(())
}
fn number<T: FromStr>(words: &mut SplitAsciiWhitespace<'_>, label: &str) -> Result<T> {
    word(words, label)?
        .parse()
        .map_err(|_| format!("invalid {label}").into())
}

fn main() -> Result<()> {
    let (args, options) =
        thermochemical::Options::extract_os(&std::env::args_os().collect::<Vec<_>>())?;
    if args.len() != 3 {
        return Err("usage: gas_field_pulse INPUT.pulse NEW_OUTPUT.csv --gas-component SPECIES.pack NASA9.pack TRANSPORT.pack ... --gas-water none|ID --gas-phase frozen|liquid-water".into());
    }
    if options.record.is_some() {
        return Err("use the pulse file's stations, not --gas-record".into());
    }
    let output_path = Path::new(&args[2]);
    if output_path.try_exists()? {
        return Err("output already exists; choose a new CSV path".into());
    }
    let mut bytes = Vec::new();
    File::open(&args[1])?
        .take((MAX_INPUT_BYTES + 1) as u64)
        .read_to_end(&mut bytes)?;
    if bytes.len() > MAX_INPUT_BYTES {
        return Err("pulse input exceeds 64 KiB".into());
    }
    let text = std::str::from_utf8(&bytes)?;
    let mut words = text.split_ascii_whitespace();
    tag(&mut words, "gas_field_pulse_v1")?;
    let shape = word(&mut words, "cylinder or cone")?;
    if !matches!(shape, "cylinder" | "cone") {
        return Err("pulse geometry must be cylinder or cone".into());
    }
    let radius: f64 = number(&mut words, "inlet radius [m]")?;
    let outlet_radius: f64 = if shape == "cone" {
        number(&mut words, "outlet radius [m]")?
    } else {
        radius
    };
    let length: f64 = number(&mut words, "length [m]")?;
    tag(&mut words, "pressure")?;
    let pressure: f64 = number(&mut words, "static pressure [Pa]")?;
    tag(&mut words, "mesh")?;
    let elements: usize = number(&mut words, "element count")?;
    let order: usize = number(&mut words, "polynomial order")?;
    tag(&mut words, "time")?;
    let dt: f64 = number(&mut words, "step [s]")?;
    let steps: usize = number(&mut words, "step count")?;
    tag(&mut words, "pulse")?;
    let amplitude: f64 = number(&mut words, "pressure amplitude [Pa]")?;
    let center: f64 = number(&mut words, "pulse center [m]")?;
    let half_width: f64 = number(&mut words, "pulse half-width [m]")?;
    if [
        radius,
        outlet_radius,
        length,
        pressure,
        dt,
        amplitude,
        half_width,
    ]
    .iter()
    .any(|v| !v.is_finite() || *v <= 0.0)
        || !center.is_finite()
        || center <= half_width
        || center + half_width >= length
    {
        return Err("geometry, pressure, time and pulse scales must be positive and finite; pulse support must lie strictly inside the duct".into());
    }
    if !(1..=64).contains(&elements)
        || !(1..=8).contains(&order)
        || !(1..=10_000).contains(&steps)
        || 2 * elements * (order + 1) > 256
        || !(dt * steps as f64).is_finite()
    {
        return Err("example budget: 1..=64 elements, order 1..=8, at most 256 DOFs, 1..=10000 steps and finite duration".into());
    }
    tag(&mut words, "species")?;
    let count: usize = number(&mut words, "species count")?;
    if !(1..=16).contains(&count) {
        return Err("example species count must be in 1..=16".into());
    }
    let ids = (0..count)
        .map(|_| word(&mut words, "canonical species ID").map(str::to_owned))
        .collect::<Result<Vec<_>>>()?;
    let mut canonical: Vec<_> = (0..count).collect();
    canonical.sort_by(|&a, &b| ids[a].cmp(&ids[b]));
    if canonical
        .windows(2)
        .any(|pair| ids[pair[0]] == ids[pair[1]])
    {
        return Err("species IDs must be distinct".into());
    }
    tag(&mut words, "stations")?;
    let stations: usize = number(&mut words, "station count")?;
    if !(2..=64).contains(&stations) {
        return Err("example station count must be in 2..=64".into());
    }
    let mut positions = Vec::with_capacity(stations);
    let mut states = Vec::with_capacity(stations);
    for _ in 0..stations {
        tag(&mut words, "station")?;
        positions.push(number::<f64>(&mut words, "station coordinate [m]")?);
        let temperature_k = number(&mut words, "temperature [K]")?;
        let fractions = (0..count)
            .map(|_| number::<f64>(&mut words, "mole fraction"))
            .collect::<Result<Vec<_>>>()?;
        states.push(PrimitiveGasState {
            temperature_k,
            composition: Composition::new(
                CompositionBasis::MoleFraction,
                canonical.iter().map(|&i| fractions[i]).collect(),
            )?,
        });
    }
    let viscothermal = match words.next() {
        None => None,
        Some(kind @ ("viscothermal" | "spatial_viscothermal")) => {
            let poles = number::<usize>(&mut words, "viscothermal pole count")?;
            let band = [
                number::<f64>(&mut words, "viscothermal minimum frequency [Hz]")?,
                number::<f64>(&mut words, "viscothermal maximum frequency [Hz]")?,
            ];
            let error = number::<f64>(&mut words, "viscothermal maximum relative operator error")?;
            let quadrature = if kind == "spatial_viscothermal" {
                Some(number::<usize>(
                    &mut words,
                    "wall Gauss points per element",
                )?)
            } else {
                None
            };
            Some((poles, band, error, quadrature))
        }
        Some(_) => {
            return Err(
                "expected optional viscothermal or spatial_viscothermal record after stations"
                    .into(),
            );
        }
    };
    if words.next().is_some() {
        return Err("pulse input has trailing fields".into());
    }

    let loaded = options.load()?;
    let bindings = loaded.bindings();
    let registry = loaded
        .registry(&bindings)
        .ok_or("explicit gas source packs are required")?;
    if registry.species.len() != count {
        return Err("load exactly one source-pack triple for each declared species".into());
    }
    let sources = canonical
        .iter()
        .map(|&i| {
            registry
                .species
                .iter()
                .find(|source| source.species.association().species().as_str() == ids[i])
                .cloned()
                .ok_or_else(|| format!("missing source packs for {}", ids[i]))
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    // Positive dyadic seed retains every source, including zero-inlet species.
    // It is never evaluated at a physical coordinate; stations own composition.
    let fraction = 1.0 / count.next_power_of_two() as f64;
    let mut seed = vec![fraction; count];
    seed[count - 1] = 1.0 - fraction * (count - 1) as f64;
    let model = from_packs(
        &sources,
        Composition::new(CompositionBasis::MoleFraction, seed)?,
        registry.water_species,
    )?;
    let field = AxialGasField::new(
        positions,
        pressure,
        states,
        model,
        QueryPoint::new(),
        registry.phase,
    )?;
    let duct = Duct {
        segments: vec![if shape == "cone" {
            Segment::Cone {
                inlet_radius: radius,
                outlet_radius,
                length,
            }
        } else {
            Segment::Cylinder { radius, length }
        }],
    };
    let dg = DgDuct::new_with_gas_field(
        &duct,
        &field,
        DgMesh { elements, order },
        DgBoundary::Closed,
        DgBoundary::Matched,
    )?;
    let dg = if let Some((poles, band, error, quadrature)) = viscothermal {
        if let Some(quadrature_points) = quadrature {
            dg.with_spatial_viscothermal(DgWallResolution {
                poles_per_field: poles,
                quadrature_points,
                band_rad_s: band.map(|frequency| core::f64::consts::TAU * frequency),
                maximum_relative_error: error,
            })?
        } else {
            dg.with_viscothermal(poles, band, error)?
        }
    } else {
        dg
    };
    if dg.total_dofs() > 256 {
        return Err(
            "example's 256-DOF budget includes all acoustic and viscothermal memory states".into(),
        );
    }
    let propagator = dg.propagator(dt, 1e-13)?;
    let mut state = dg.project_initial_state(
        |x| {
            let r = (x - center) / half_width;
            if r.abs() < 1.0 {
                let shape = (1.0 - r * r) * (1.0 - r * r);
                amplitude * shape * shape
            } else {
                0.0
            }
        },
        |_| 0.0,
    );
    if !(state.initial_energy_j().is_finite() && state.initial_energy_j() > 0.0) {
        return Err(
            "initial pulse is non-finite or unresolved; refine the mesh or widen it".into(),
        );
    }

    // Finish numerical admission and the bounded render before creating output.
    let mut csv = Vec::with_capacity((steps + 1) * 140);
    writeln!(
        csv,
        "time_s,inlet_pressure_pa,outlet_pressure_pa,stored_energy_j,{}",
        if dg.has_viscothermal_memory() {
            "inferred_total_loss_j"
        } else {
            "inferred_matched_loss_j"
        }
    )?;
    let mut inferred_loss = 0.0;
    for step in 0..=steps {
        if step > 0 {
            let row = propagator.step(&mut state, &[])?;
            inferred_loss = row.ledger.cumulative_dissipated_j;
        }
        propagator.sync_nodal(&mut state)?;
        // LGL endpoint traces at x=0 and x=L, without spatial interpolation.
        let p = state.pressure();
        let row = [
            step as f64 * dt,
            p[0],
            p[p.len() - 1],
            state.stored_energy_j(),
            inferred_loss,
        ];
        if row.iter().any(|v| !v.is_finite()) {
            return Err("non-finite physical output; no CSV was created".into());
        }
        writeln!(
            csv,
            "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            row[0], row[1], row[2], row[3], row[4]
        )?;
        if csv.len() > MAX_OUTPUT_BYTES {
            return Err("CSV exceeds the example's 4 MiB output budget".into());
        }
    }
    let wall_power_w = dg.wall_dissipated_power_w(&state)?;
    let mut output = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output_path)?;
    output.write_all(&csv)?;
    if let Some(model) = dg.viscothermal_model() {
        eprintln!(
            "Bessel wall poles={}; sampled complex error={:e}, real loss error={:e}; wall memory is included in stored_energy_j",
            model.series_terms().len(),
            model.max_complex_error(),
            model.max_real_loss_error()
        );
    }
    if let Some(admission) = dg.spatial_viscothermal_admission() {
        eprintln!(
            "spatial Bessel wall poles={}, Gauss points/element={}; sampled complex error={:e}, real loss error={:e}; physical wall energy={:e} J, instantaneous wall loss={:e} W",
            admission.resolution.poles_per_field,
            admission.resolution.quadrature_points,
            admission.max_complex_error,
            admission.max_real_loss_error,
            state.wall_stored_energy_j(),
            wall_power_w
        );
    }
    eprintln!(
        "wrote {} rows; boundary and any wall loss is the inferred E(0)-E(t) deficit, not an independent energy-balance check",
        steps + 1
    );
    Ok(())
}
