//! Headless SI photometry from ordinary Rust scene and fixture data.
//!
//! Run: cargo run -p fs-render --features photometry --example photometric_stage
//! Add `-- --csv` to emit the direct-illuminance calculation grid as CSV.
//! These declared fixture profiles and gray haze are examples, not measurements.

use std::error::Error;

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::camera::{Aperture, CameraProjection, PhysicalCamera};
use fs_render::charts::TriMesh;
use fs_render::instances::RigidTransform;
use fs_render::photometry::image::LuminanceJob;
use fs_render::photometry::medium::{HomogeneousMedium, PhaseFunction, ScatterSettings};
use fs_render::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
use fs_render::photometry::{CalculationGrid, IlluminanceJob, PhotometricScene, TransportOptions};
use fs_render::spectral::LiftedSpectrum;
use fs_render::tracer::{Material, Primitive, Shape};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let csv = match args.as_slice() {
        [] => false,
        [flag] if flag == "--csv" => true,
        _ => return Err("usage: photometric_stage [--csv]".into()),
    };
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x5052_4556_495a,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        run(&cx, csv)
    })
}

fn run(cx: &Cx<'_>, csv: bool) -> Result<(), Box<dyn Error>> {
    // Right-handed world, Z up, metres. Local source +Z is its optical axis.
    // This proper rotation turns +Z downward; changing these transforms and
    // constructing a new scene supplies another pan/tilt/fixture snapshot.
    let profile = AngularDistribution::beam_field(
        24.0_f64.to_radians(),
        40.0_f64.to_radians(),
        60.0_f64.to_radians(),
    )?;
    let mut sources = Vec::new();
    for (x_m, dimmer) in [(-1.5, 1.0), (1.5, 0.65)] {
        sources.push(PhotometricSource::try_new(
            RigidTransform::try_new([1.0, 0.0, 0.0, 0.0], [x_m, 0.0, 5.0])?,
            12_000.0,
            dimmer,
            profile.clone(),
            SourceSpectrum::equal_energy(),
        )?);
    }
    // Original spectral-tracer geometry and material, with a constant 0.5
    // spectral reflectance. The observer does not create a parallel mesh model.
    let surfaces = [Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [-5.0, -4.0, 0.0],
                [5.0, -4.0, 0.0],
                [5.0, 4.0, 0.0],
                [-5.0, 4.0, 0.0],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )),
        material: Material::Lambertian {
            reflectance: LiftedSpectrum { c: [0.0; 3] },
        },
        emission: None,
    }];
    let medium = HomogeneousMedium::try_new(
        0.005, // absorption, m^-1
        0.025, // scattering, m^-1
        PhaseFunction::HenyeyGreenstein { g: 0.7 },
    )?;
    let scene = PhotometricScene::try_new(
        &surfaces,
        &sources,
        medium,
        TransportOptions {
            max_outputs: 4096,
            ..TransportOptions::default()
        },
        cx,
    )?;
    let grid = CalculationGrid::try_new(
        Point3::new(-3.0, -2.0, 0.0),
        Vec3::new(0.5, 0.0, 0.0),
        Vec3::new(0.0, 0.5, 0.0),
        13,
        9,
    )?;
    let mut illuminance = IlluminanceJob::for_grid(&scene, grid)?;
    while !illuminance.is_complete() {
        illuminance.advance(cx, 32)?; // Completed prefixes are usable immediately.
    }
    if csv {
        println!("x_m,y_m,z_m,direct_lux");
        for (index, estimate) in illuminance.results().iter().enumerate() {
            let point = grid.probe(index)?.position_m();
            println!("{},{},{},{:.9}", point.x, point.y, point.z, estimate.lux);
        }
        return Ok(());
    }
    let peak_lux = illuminance
        .results()
        .iter()
        .map(|v| v.lux)
        .fold(0.0, f64::max);
    let mean_lux = illuminance.results().iter().map(|v| v.lux).sum::<f64>() / grid.len() as f64;
    println!("2 fixtures, 19,800 emitted lm; 13 x 9 stage grid at 0.5 m spacing");
    println!("Direct illuminance: mean {mean_lux:.3} lx, peak {peak_lux:.3} lx");

    let camera = PhysicalCamera::try_look_at(
        Point3::new(7.0, -10.0, 5.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 0.0, 1.0),
        CameraProjection::try_vertical_fov(50.0_f64.to_radians())?,
        10.0,
        Aperture::try_circular(0.0)?, // Pinhole; no exposure or tone mapping.
    )?;
    let mut image = LuminanceJob::new(
        &scene,
        &camera,
        32,
        24,
        20.0,
        ScatterSettings {
            initial_panels: 16,
            max_panels: 512,
            absolute_tolerance_cd_m2: 1.0e-4,
            relative_tolerance: 1.0e-3,
        },
    )?;
    while image.results().len() < 32 * 24 {
        image.advance(cx, 64)?;
    }
    let peak = image
        .results()
        .iter()
        .map(|v| v.luminance_cd_m2)
        .fold(0.0, f64::max);
    let haze = image
        .results()
        .iter()
        .map(|v| v.single_scatter.luminance_cd_m2)
        .fold(0.0, f64::max);
    let unconverged = image
        .results()
        .iter()
        .filter(|v| !v.single_scatter.tolerance_met)
        .count();
    println!(
        "32 x 24 luminance image: peak {peak:.6} cd/m^2; peak single-scattered haze {haze:.6} cd/m^2"
    );
    println!(
        "Haze refinement tolerance unmet at {unconverged} pixels; estimates are not certified bounds."
    );
    println!(
        "Total error: unknown (fixture/observer/geometry error, interreflection, and multiple scattering)."
    );
    Ok(())
}
