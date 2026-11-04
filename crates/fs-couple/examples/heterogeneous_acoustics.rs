//! Source packs and an authored primitive gas gradient -> full 3D harmonic pressure.
//!
//! cargo run -p fs-couple --example heterogeneous_acoustics -- FIELD HZ GRID
//!   --gas-component SPECIES.pack NASA9.pack TRANSPORT.pack ...
//!   --gas-water none --gas-phase frozen
//!
//! Supply those arguments on one command line. FIELD uses the existing
//! primitive-column grammar, for example:
//! pressure 100000 stations 2 mole 2 A B
//! station 0 300 1 0
//! station 1 600 0 1
//!
//! Species IDs refer only to the explicitly supplied source packs. Output is
//! complex peak pressure [Pa], particle velocity [m/s] and active/reactive
//! intensity [W/m^2] at original element centroids, in exp(+i omega t).
//! Each observation uses the same original local material as the wave equation.
//! The rectangular cavity has length L and width/height L/5. A phased velocity on x=0
//! excites transverse structure; x=L has a local passive impedance, while side
//! walls are rigid. This outlet is not an exact multimode radiation condition.
//! Refine GRID and frequency/wavelength resolution before interpreting results.

#[path = "../src/bin/music_render/json.rs"]
mod json;
#[path = "../src/bin/music_render/thermochemical.rs"]
mod thermochemical;
use json::json_string;

use fs_couple::cavity_modal::interior::heterogeneous::{AxialGasEmbedding, AxialGasInterior};
use fs_feec::acoustics::helmholtz::heterogeneous::{BoundaryAdmittance, HarmonicLimits};
use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_matdb::QueryPoint;
use fs_material::gas_field::{AxialGasField, input::PrimitiveGasColumns};
use fs_math::{c64::C64, det};
use std::{error::Error, fs::File, io::Read, str::SplitAsciiWhitespace};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

fn word<'a>(words: &mut SplitAsciiWhitespace<'a>) -> Result<&'a str> {
    words.next().ok_or_else(|| "missing field input".into())
}

fn tag(words: &mut SplitAsciiWhitespace<'_>, expected: &str) -> Result<()> {
    if word(words)? != expected {
        return Err(format!("expected {expected}").into());
    }
    Ok(())
}

fn main() -> Result<()> {
    let (args, options) =
        thermochemical::Options::extract_os(&std::env::args_os().collect::<Vec<_>>())?;
    if args.len() != 4 || options.record.is_some() {
        return Err("usage: heterogeneous_acoustics FIELD HZ GRID --gas-component SPECIES NASA9 TRANSPORT ... --gas-water none|ID --gas-phase frozen|liquid-water".into());
    }
    let frequency: f64 = args[2].to_str().ok_or("frequency must be UTF-8")?.parse()?;
    let grid: usize = args[3].to_str().ok_or("GRID must be UTF-8")?.parse()?;
    if !(1..=6).contains(&grid) || !frequency.is_finite() || frequency <= 0.0 {
        return Err("example requires GRID in 1..=6 and positive finite HZ".into());
    }
    let loaded = options.load()?;
    let bindings = loaded.bindings();
    let registry = loaded
        .registry(&bindings)
        .ok_or("explicit thermochemical packs are required")?;
    let mut text = String::new();
    File::open(&args[1])?
        .take(65_537)
        .read_to_string(&mut text)?;
    if text.len() > 65_536 {
        return Err("field input exceeds 64 KiB".into());
    }
    let mut words = text.split_ascii_whitespace();
    tag(&mut words, "pressure")?;
    let pressure = word(&mut words)?.parse()?;
    tag(&mut words, "stations")?;
    let count: usize = word(&mut words)?.parse()?;
    if !(2..=64).contains(&count) {
        return Err("example needs 2..=64 stations".into());
    }
    let columns = PrimitiveGasColumns::read(&mut words, Some(&registry), count)?;
    let mut positions = Vec::with_capacity(count);
    let mut states = Vec::with_capacity(count);
    for _ in 0..count {
        tag(&mut words, "station")?;
        positions.push(word(&mut words)?.parse()?);
        states.push(columns.read_state(&mut words)?);
    }
    if words.next().is_some() {
        return Err("trailing fields after final gas station".into());
    }
    let (model, phase) = columns.into_parts();
    let field = AxialGasField::new(positions, pressure, states, model, QueryPoint::new(), phase)?;
    let length = field.length_m();
    let (complex, mut positions) = fs_feec::fixtures::kuhn_cube(grid);
    for point in &mut positions {
        point[0] *= length;
        point[1] *= 0.2 * length;
        point[2] *= 0.2 * length;
    }
    let inlet = field.state_at(0.0)?;
    let inlet = inlet.state();
    let geometry = HelmholtzFe::new(
        &complex,
        &positions,
        1,
        Medium {
            density_kg_m3: inlet.density,
            sound_speed_m_s: inlet.sound_speed,
        },
    )?;
    let acoustic = AxialGasInterior::new(
        &geometry,
        &field,
        AxialGasEmbedding::new([0.0; 3], [1.0, 0.0, 0.0])?,
        5,
        &mut || false,
    )?;
    let faces = geometry.boundary_faces();
    let outlet_faces = faces
        .iter()
        .copied()
        .filter(|face| face.iter().all(|&v| positions[v as usize][0] == length))
        .collect();
    let outlet = field.state_at(length)?;
    let outlet = outlet.state();
    let boundary = BoundaryAdmittance {
        faces: outlet_faces,
        admittance_m_pa_s: C64::from_re(1.0 / (outlet.density * outlet.sound_speed)),
    };
    let prepared = acoustic.operator().prepare(
        frequency,
        &[boundary],
        HarmonicLimits {
            max_dofs: 512,
            max_dense_bytes: 64 * 1024 * 1024,
            residual_tolerance: 1e-10,
        },
        &mut || false,
    )?;
    let transverse_wavenumber = 0.4 * core::f64::consts::TAU * frequency / inlet.sound_speed;
    let velocities = faces
        .iter()
        .copied()
        .filter(|face| face.iter().all(|&v| positions[v as usize][0] == 0.0))
        .map(|face| {
            (
                face,
                face.map(|vertex| {
                    let angle = -transverse_wavenumber * positions[vertex as usize][1];
                    C64::new(det::cos(angle), det::sin(angle)).scale(0.001)
                }),
            )
        })
        .collect::<Vec<_>>();
    let result = prepared.solve(&velocities, &mut || false)?;
    eprintln!(
        "{} pressure DOFs, f={} Hz, input={} W, boundary={} W, balance={} W, backward_error={}",
        geometry.ndof(),
        frequency,
        result.input_power_w,
        result.boundary_dissipation_w,
        result.power_balance_residual_w(),
        result.backward_error
    );
    println!(concat!(
        "tet,x_m,y_m,z_m,pressure_re_pa,pressure_im_pa,",
        "velocity_x_re_m_s,velocity_x_im_m_s,velocity_y_re_m_s,velocity_y_im_m_s,",
        "velocity_z_re_m_s,velocity_z_im_m_s,",
        "active_intensity_x_w_m2,active_intensity_y_w_m2,active_intensity_z_w_m2,",
        "reactive_intensity_x_w_m2,reactive_intensity_y_w_m2,reactive_intensity_z_w_m2"
    ));
    for element in 0..geometry.tetrahedra().len() {
        let probe = acoustic.probe(element, [0.25; 4])?;
        let sample = probe.harmonic(&result.pressure_pa, prepared.frequency_hz())?;
        let center = probe.position_m();
        print!(
            "{element},{:e},{:e},{:e},{:e},{:e}",
            center[0], center[1], center[2], sample.pressure_pa.re, sample.pressure_pa.im
        );
        for velocity in sample.particle_velocity_m_s {
            print!(",{:e},{:e}", velocity.re, velocity.im);
        }
        for intensity in sample.active_intensity_w_m2.iter().chain(&sample.reactive_intensity_w_m2) {
            print!(",{intensity:e}");
        }
        println!();
    }
    Ok(())
}
