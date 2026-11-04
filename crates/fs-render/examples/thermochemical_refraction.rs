//! Original thermochemical source packs -> continuous gas optics -> native rays.
//!
//! Run with `--features material-state,volumes --example thermochemical_refraction`
//! followed by FIELD OPTICS SPAN ARC_LENGTH_M WATER_ID_OR_none and one or more
//! SPECIES.pack NASA9.pack TRANSPORT.material-card.pack triples. Each species
//! requires its own original normalized sources; no optical species presets or
//! constant-heat-capacity substitute are supplied by this example.
//!
//! FIELD uses the original primitive-column reader (SI pressure/length/temperature):
//! ```text
//! pressure P_PA stations COUNT mole SPECIES_COUNT ID...
//! station X_M T_K MOLE_FRACTION...
//! ...
//! ```
//! Stations start at zero. The phase policy is explicitly FrozenGas. OPTICS is:
//! ```text
//! molecular COUNT
//! species ID N_MINUS_1 N_REF_PER_M3 LAMBDA_REF_NM DEPOLARIZATION B_UM2 C_UM4 LO_NM HI_NM
//! ...
//! ```
//! Values must come from the caller's declared optical sources. The optical file
//! hash identifies those input bytes; it does not attest their measured accuracy.
//! All field species must have optical rows and be present at the selected span
//! midpoint used to bind the original source catalogue. Zero endpoint fractions
//! are supported. No trace species is inserted to make an attachment succeed.
//!
//! The axial field is embedded along world x. Each ray starts at the selected
//! span midpoint and travels initially along +z. The original analytic field
//! supplies temperature, partial densities, refractive gradients/Hessians and
//! molecular scattering. Five supported wavelengths each reuse one prepared
//! spectral law across native adaptive ray/bundle stages. CSV reports physical
//! optical/group path lengths and the local angular footprint on the z-normal
//! plane THROUGH THE ACTUAL ENDPOINT, not a fabricated target-plane intersection.
//! The footprint is a local geometric ray tube, not a finite-pixel caustic.
//!
//! A path leaving the one smooth span refuses. Station knots need explicit
//! interface treatment, and source/step refinement remains the caller's task.
//! This example reports local elastic scattering; it does not sample collisions
//! or interpret scattering as true absorption, thermal emission or lost heat.

use fs_blake3::{ContentHash, DomainHasher};
use fs_geom::{Point3, Vec3};
use fs_matdb::{
    NormalizedMaterialCardPack, NormalizedModelPack, NormalizedSpeciesPack, QueryPoint,
};
use fs_material::{
    gas_field::{AxialGasField, input::PrimitiveGasColumns},
    gas_scattering::{
        MolecularRayleigh,
        composition::{MolecularCauchyDispersion, MolecularRefractionModel},
        thermochemical::{SpeciesMolecularRayleigh, ThermochemicalRayleighState},
    },
    gas_species::input::{ThermochemicalPhase, ThermochemicalSources},
    state_point::MaterialPropertySelection,
    thermochemical_gas::packs::SpeciesPackRef,
};
use fs_qty::chemistry::SpeciesId;
use fs_render::{
    charts::Ray,
    gradient_index::{
        GradientIndexOptions,
        bundle::{BundleOptions, ProjectionOptions, propagate_bundle},
    },
    material_state::gas_scattering::AxialGasCompositionField,
};
use std::{
    error::Error,
    ffi::OsString,
    fmt::Write as _,
    io::{Read, Write as _},
    path::Path,
    str::SplitAsciiWhitespace,
};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const MAX_PACK_BYTES: usize = 16 * 1024 * 1024;
const MAX_TOTAL_BYTES: usize = 64 * 1024 * 1024;
const USAGE: &str = "thermochemical_refraction FIELD OPTICS SPAN ARC_LENGTH_M WATER_ID_OR_none \
    SPECIES.pack NASA9.pack TRANSPORT.material-card.pack [SPECIES.pack NASA9.pack TRANSPORT.material-card.pack ...]";

struct Source {
    species: NormalizedSpeciesPack,
    nasa9: NormalizedModelPack,
    transport: NormalizedMaterialCardPack,
}

fn read(path: &Path, total: &mut usize) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)?
        .take((MAX_PACK_BYTES + 1) as u64)
        .read_to_end(&mut bytes)?;
    if bytes.is_empty() || bytes.len() > MAX_PACK_BYTES {
        return Err(format!("{} must contain 1 byte..16 MiB", path.display()).into());
    }
    *total = total
        .checked_add(bytes.len())
        .ok_or("input byte count overflow")?;
    if *total > MAX_TOTAL_BYTES {
        return Err("inputs exceed the aggregate 64 MiB budget".into());
    }
    Ok(bytes)
}

fn word<'a>(words: &mut SplitAsciiWhitespace<'a>) -> Result<&'a str> {
    words.next().ok_or_else(|| "missing input value".into())
}

fn expect(words: &mut SplitAsciiWhitespace<'_>, expected: &str) -> Result<()> {
    let actual = word(words)?;
    if actual != expected {
        return Err(format!("expected {expected}, found {actual}").into());
    }
    Ok(())
}

fn scalar(words: &mut SplitAsciiWhitespace<'_>) -> Result<f64> {
    let value: f64 = word(words)?.parse()?;
    if !value.is_finite() {
        return Err("input scalars must be finite".into());
    }
    Ok(value)
}

fn field(text: &str, sources: &ThermochemicalSources<'_>) -> Result<AxialGasField> {
    let mut words = text.split_ascii_whitespace();
    expect(&mut words, "pressure")?;
    let pressure = scalar(&mut words)?;
    expect(&mut words, "stations")?;
    let count: usize = word(&mut words)?.parse()?;
    if !(2..=64).contains(&count) {
        return Err("this example admits 2..=64 original stations".into());
    }
    let columns = PrimitiveGasColumns::read(&mut words, Some(sources), count)?;
    let mut positions = Vec::with_capacity(count);
    let mut states = Vec::with_capacity(count);
    for _ in 0..count {
        expect(&mut words, "station")?;
        positions.push(scalar(&mut words)?);
        states.push(columns.read_state(&mut words)?);
    }
    if words.next().is_some() {
        return Err("unexpected suffix in field input".into());
    }
    let (model, phase) = columns.into_parts();
    Ok(AxialGasField::new(
        positions,
        pressure,
        states,
        model,
        QueryPoint::new(),
        phase,
    )?)
}

fn optical_rows(
    text: &str,
    identity: ContentHash,
) -> Result<Vec<(SpeciesMolecularRayleigh, MolecularCauchyDispersion)>> {
    let mut words = text.split_ascii_whitespace();
    expect(&mut words, "molecular")?;
    let count: usize = word(&mut words)?.parse()?;
    if !(1..=128).contains(&count) {
        return Err("optical species count must be in 1..=128".into());
    }
    let mut rows = Vec::with_capacity(count);
    for _ in 0..count {
        expect(&mut words, "species")?;
        let species = SpeciesId::new(word(&mut words)?)?;
        let optical = MolecularRayleigh::try_new(
            scalar(&mut words)?,
            scalar(&mut words)?,
            scalar(&mut words)?,
            scalar(&mut words)?,
            identity,
        )?;
        let dispersion = MolecularCauchyDispersion::try_new(
            optical,
            scalar(&mut words)?,
            scalar(&mut words)?,
            [scalar(&mut words)?, scalar(&mut words)?],
            identity,
        )?;
        rows.push((SpeciesMolecularRayleigh { species, optical }, dispersion));
    }
    if words.next().is_some() {
        return Err("unexpected suffix in optical input".into());
    }
    // Original attachments are canonical. Keep the spectral rows attached to
    // their actual species while ordering; never sort optical values alone.
    rows.sort_by(|(a, _), (b, _)| a.species.cmp(&b.species));
    Ok(rows)
}

fn utf8(argument: &OsString) -> Result<&str> {
    argument
        .to_str()
        .ok_or_else(|| "scalar/identifier argument is not UTF-8".into())
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() < 8 || (args.len() - 5) % 3 != 0 || (args.len() - 5) / 3 > 128 {
        return Err(USAGE.into());
    }
    let span: usize = utf8(&args[2])?.parse()?;
    let distance: f64 = utf8(&args[3])?.parse()?;
    if !distance.is_finite() || distance <= 0.0 {
        return Err("ARC_LENGTH_M must be finite and positive".into());
    }
    let water = match utf8(&args[4])? {
        "none" => None,
        id => Some(SpeciesId::new(id)?),
    };
    let mut total = 0;
    let field_text = String::from_utf8(read(Path::new(&args[0]), &mut total)?)?;
    let optics_text = String::from_utf8(read(Path::new(&args[1]), &mut total)?)?;
    let mut sources = Vec::with_capacity((args.len() - 5) / 3);
    for paths in args[5..].chunks_exact(3) {
        sources.push(Source {
            species: NormalizedSpeciesPack::from_bytes(&read(Path::new(&paths[0]), &mut total)?)?,
            nasa9: NormalizedModelPack::from_bytes(&read(Path::new(&paths[1]), &mut total)?)?,
            transport: NormalizedMaterialCardPack::from_bytes(&read(
                Path::new(&paths[2]),
                &mut total,
            )?)?,
        });
    }
    let bindings: Vec<_> = sources
        .iter()
        .map(|source| SpeciesPackRef {
            species: &source.species,
            nasa9: Some(&source.nasa9),
            transport: Some((
                &source.transport,
                MaterialPropertySelection::SingleClaimOnly,
            )),
        })
        .collect();
    let registry = ThermochemicalSources {
        species: &bindings,
        water_species: water.as_ref(),
        phase: ThermochemicalPhase::FrozenGas,
    };
    let gas = field(&field_text, &registry)?;
    let end = span.checked_add(1).ok_or("span index overflow")?;
    let lo = *gas
        .stations_m()
        .get(span)
        .ok_or("span index outside field")?;
    let hi = *gas
        .stations_m()
        .get(end)
        .ok_or("span index outside field")?;
    let midpoint = lo + (hi - lo) * 0.5;
    let mut optical_hash = DomainHasher::new("org.frankensim.example.molecular-optical-input.v1");
    optical_hash.update(optics_text.as_bytes());
    let rows = optical_rows(&optics_text, optical_hash.finalize())?;
    let references: Vec<_> = rows.iter().map(|(row, _)| row.clone()).collect();
    let reference_state = gas.state_at(midpoint)?;
    if reference_state.thermodynamics().receipt().species().len()
        != gas.model().thermodynamics().components().len()
    {
        return Err("all catalogue species must be present at the selected span midpoint".into());
    }
    let attached = ThermochemicalRayleighState::try_new(reference_state, &references)?;
    let dispersion: Vec<_> = rows.iter().map(|(_, law)| *law).collect();
    let model =
        MolecularRefractionModel::from_thermochemical(&attached)?.with_dispersion(&dispersion)?;
    let mut snapshot = DomainHasher::new("org.frankensim.example.thermochemical-refraction.v1");
    snapshot.update(&(field_text.len() as u64).to_le_bytes());
    snapshot.update(field_text.as_bytes());
    snapshot.update(model.content_identity().as_bytes());
    snapshot.update(&(span as u64).to_le_bytes());
    let field = AxialGasCompositionField::try_new(
        &gas,
        &model,
        span,
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        snapshot.finalize(),
    )?;
    let refractive = field.bind_refraction()?;
    let scattering =
        refractive.bind_scattering_field(field.component_number_density_upper_bounds())?;
    let ray = Ray {
        origin: Point3::new(midpoint, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    };
    let step = distance.min(hi - lo) / 128.0;
    let options = BundleOptions {
        integration: GradientIndexOptions {
            initial_step_m: step,
            min_step_m: 1e-10_f64.min(step / 1024.0),
            max_step_m: step,
            ..GradientIndexOptions::default()
        },
        ..BundleOptions::default()
    };
    let [lambda_lo, lambda_hi] = model.wavelength_range_nm();
    let mut csv = String::from(
        "wavelength_nm,arc_length_m,x_m,y_m,z_m,dir_x,dir_y,dir_z,phase_path_m,group_path_m,group_delay_s,n,group_index,temperature_k,density_kg_m3,scattering_per_m,depolarization,scattering_majorant_per_m,area_m2_per_sr,projection_condition,accepted_steps,rejected_steps,field_evaluations\n",
    );
    for i in 0..5 {
        let wavelength = lambda_lo + (lambda_hi - lambda_lo) * (i as f64 / 4.0);
        let prepared = refractive.at_wavelength(wavelength)?;
        let bundle = propagate_bundle(
            &prepared,
            ray,
            Vec3::new(1.0, 0.0, 0.0),
            wavelength,
            distance,
            options,
        )?;
        let path = bundle.propagation();
        let projected =
            bundle.project_at_endpoint(Vec3::new(0.0, 0.0, 1.0), ProjectionOptions::default())?;
        let position = path.ray.origin;
        let local = field.sample_state(position)?;
        let scatter = scattering.sample_scattering_optics(position, wavelength)?;
        let majorant = scattering.scattering_majorant_optics(wavelength)?;
        let group_path = path
            .group_optical_path_length_m
            .ok_or("missing group path")?;
        let scalars = [
            wavelength,
            path.distance_m,
            position.x,
            position.y,
            position.z,
            path.ray.dir.x,
            path.ray.dir.y,
            path.ray.dir.z,
            path.optical_path_length_m,
            group_path,
            group_path / 299_792_458.0,
            path.final_sample.phase_index(),
            path.final_sample
                .group_index()
                .ok_or("missing group index")?,
            local.state().temperature,
            local.state().density,
            scatter.scattering_per_m,
            scatter.depolarization_ratio,
            majorant.scattering_per_m,
            projected.area_per_solid_angle_m2(),
            projected.condition_number(),
        ];
        for value in scalars {
            write!(&mut csv, "{value:.16e},")?;
        }
        writeln!(
            &mut csv,
            "{},{},{}",
            path.accepted_steps, path.rejected_steps, path.field_evaluations
        )?;
    }
    // Emit the complete finite result only after every original-source and
    // native propagation query succeeds. Diagnostics stay off the CSV stream.
    eprintln!(
        "source_field={} optical_model={} frozen_span={span} [{lo},{hi}] m",
        field.content_identity().to_hex(),
        model.content_identity().to_hex()
    );
    std::io::stdout().lock().write_all(csv.as_bytes())?;
    Ok(())
}
