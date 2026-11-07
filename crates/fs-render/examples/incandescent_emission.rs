//! Absolute visible emission from an explicitly authored gray thermal source.
//!
//! Run with `cargo run -p fs-render --example incandescent_emission --features photometry`.
//! This demonstration declares its area and wavelength-independent emissivity;
//! it does not calibrate tungsten or infer filament electrical/thermal dynamics.

use fs_blake3::ContentHash;
use fs_render::thermal_surface::{DiffuseThermalSurface, GrayDiffuseThermalSurface};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let emitting_area_m2 = 1.0e-4;
    let emissivity = 0.8;
    // One declared body identity, with temperature retained by each source view.
    let source_identity = ContentHash([0x49; 32]);
    println!("# Authored gray/diffuse model; no tungsten calibration.");
    println!("# Emitting area: {emitting_area_m2} m^2; emissivity: {emissivity}.");
    println!("# Visible band: 380-780 nm; gross emission from one emitting face.");
    println!(
        "temperature_k,bolometric_emitted_power_w,visible_emitted_power_w,luminous_flux_lm,luminance_cd_m2,x,y"
    );
    for temperature_k in [1000.0, 2000.0, 2856.0, 3500.0] {
        let surface =
            GrayDiffuseThermalSurface::try_new(temperature_k, emissivity, source_identity)?;
        let visible = DiffuseThermalSurface::from(surface).visible_observables()?;
        let [x, y] = visible
            .chromaticity_xy()
            .ok_or("declared incandescent source has no represented visible response")?;
        println!(
            "{temperature_k:.0},{:.9e},{:.9e},{:.9e},{:.9e},{x:.9},{y:.9}",
            emitting_area_m2 * surface.bolometric_exitance_w_m2(),
            emitting_area_m2 * visible.visible_exitance_w_m2(),
            emitting_area_m2 * visible.luminous_exitance_lm_m2(),
            visible.luminance_cd_m2(),
        );
    }
    Ok(())
}
