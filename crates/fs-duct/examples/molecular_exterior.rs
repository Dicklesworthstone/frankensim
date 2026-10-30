//! Solve one molecular, thermoviscous duct against a complete exterior matrix.
//!
//! `cargo run -p fs-duct --example molecular_exterior -- \
//!   '<gas declaration>' '<relaxation declaration>' '<source identity>' matrix.txt`
//!
//! The first two arguments use the same material grammar as `molecular_duct`.
//! The matrix file has one frequency per nonempty line: frequency [Hz], then
//! nine complex impedances as real/imaginary pairs in row-major order [Pa s/m³].
//! `#` starts a comment. Supply at least two strictly increasing frequencies.
//! Matrix order is `hole_a`, `hole_b`, `bell`; all volume flows point OUTWARD.
//! These are physical mouth-plane loads, so they include exterior radiation
//! only, without another chimney or inner end correction. For an existing
//! `fs_bem::helmholtz::aperture::ApertureImpedance` artifact, export each
//! `acoustic(row)` result, never its force/velocity `impedance` directly.
//!
//! The actual geometry is a 7.5 mm-radius, 450 mm-long rigid bore. Open holes
//! at 180 and 300 mm each have radius 2.5 mm and physical chimney height 3 mm.
//! The exterior data must describe these openings and the terminal mouth in
//! the same order, at these mouth planes, for the declared exterior medium.
//! The example cannot establish that binding from a matrix alone; its explicit
//! source identity is a caller declaration, not a geometry certificate.
//!
//! The cold wall preparation uses 32 positive poles per field and a 1% sampled
//! equilibrium wall-error ceiling over the supplied band. The source is a
//! 1e-6 m³/s peak inlet flow. CSV includes every solved mouth pressure/flow and
//! net radiation/internal power. Exterior mutual terms enter the feedback
//! solve. They are not an observation filter. This is harmonic plane-wave
//! engineering analysis, not a nonlinear performance or measured validation.

use fs_duct::compose::ComposeOptions;
use fs_duct::molecular::{MolecularDuct, MolecularWallResolution};
use fs_duct::network::CoupledTermination;
use fs_duct::{ConeSlices, Duct, HoleState, Segment};
use fs_material::gas_species::input::read_gas_record;
use fs_material::relaxation::input::{medium_from_record, read_relaxation};
use fs_math::c64::C64;
use std::io::{BufRead, BufReader};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 4 {
        return Err("usage: molecular_exterior '<gas declaration>' '<relaxation declaration>' '<source identity>' matrix.txt".into());
    }
    let mut fields = args[0].split_ascii_whitespace();
    let gas = read_gas_record(&mut fields, None, None)?;
    if fields.next().is_some() {
        return Err("extra gas fields".into());
    }
    let mut fields = args[1].split_ascii_whitespace();
    let law = read_relaxation(&mut fields)?;
    if fields.next().is_some() {
        return Err("extra relaxation fields".into());
    }
    let state = law.evaluate(&medium_from_record(&gas)?)?;
    let mut rows: Vec<(f64, Vec<C64>)> = Vec::new();
    for line in BufReader::new(std::fs::File::open(&args[3])?).lines() {
        let line = line?;
        let line = line.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        let values: Vec<f64> = line
            .split_ascii_whitespace()
            .map(str::parse)
            .collect::<Result<_, _>>()?;
        if values.len() != 19
            || values.iter().any(|x| !x.is_finite())
            || values[0] <= 0.0
            || rows.last().is_some_and(|row| values[0] <= row.0)
        {
            return Err(
                "each row needs increasing positive Hz and nine finite complex impedances".into(),
            );
        }
        rows.push((
            values[0],
            values[1..]
                .chunks_exact(2)
                .map(|z| C64::new(z[0], z[1]))
                .collect(),
        ));
    }
    if rows.len() < 2 {
        return Err("supply at least two frequency rows to declare the wall band".into());
    }
    let radius = 0.0075;
    let hole = || Segment::ToneHole {
        hole_radius: 0.0025,
        chimney_height: 0.003,
        bore_radius: radius,
        state: HoleState::Open,
    };
    let model = MolecularDuct::new(
        &Duct {
            segments: vec![
                Segment::Cylinder {
                    radius,
                    length: 0.18,
                },
                hole(),
                Segment::Cylinder {
                    radius,
                    length: 0.12,
                },
                hole(),
                Segment::Cylinder {
                    radius,
                    length: 0.15,
                },
            ],
        },
        state,
        Some(MolecularWallResolution {
            poles_per_field: 32,
            band_rad_s: [
                rows[0].0 * core::f64::consts::TAU,
                rows[rows.len() - 1].0 * core::f64::consts::TAU,
            ],
            maximum_relative_error: 0.01,
        }),
        ConeSlices::Fixed(1),
        9,
    )?;
    let labels: Vec<String> = ["hole_a", "hole_b", "bell"].map(str::to_owned).into();
    let amplitude: f64 = 1e-6;
    println!(
        "frequency_hz,zin_re_pa_s_m3,zin_im_pa_s_m3,input_power_w,exterior_power_w,internal_loss_w,hole_a_p_re_pa,hole_a_p_im_pa,hole_a_u_re_m3_s,hole_a_u_im_m3_s,hole_b_p_re_pa,hole_b_p_im_pa,hole_b_u_re_m3_s,hole_b_u_im_m3_s,bell_p_re_pa,bell_p_im_pa,bell_u_re_m3_s,bell_u_im_m3_s"
    );
    for (hz, z) in rows {
        let omega = core::f64::consts::TAU * hz;
        let ports = model.prepare_exterior_ports(omega, &labels, ComposeOptions::default())?;
        let response = ports.response(&CoupledTermination {
            omega,
            labels: labels.clone(),
            z,
            source_id: args[2].clone(),
        })?;
        let scale = amplitude * amplitude;
        print!(
            "{hz:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            response.z_in.re,
            response.z_in.im,
            0.5 * response.z_in.re * scale,
            response.exterior_power_per_squared_inlet_flow * scale,
            response.dissipated_power_per_squared_inlet_flow * scale
        );
        for (p, u) in response
            .exterior_pressures_over_u_in
            .iter()
            .zip(&response.exterior_flows_over_u_in)
        {
            print!(
                ",{:.12e},{:.12e},{:.12e},{:.12e}",
                p.re * amplitude,
                p.im * amplitude,
                u.re * amplitude,
                u.im * amplitude
            );
        }
        println!();
    }
    Ok(())
}
