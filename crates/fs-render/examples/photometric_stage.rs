//! Headless SI photometry from ordinary Rust scene and fixture data.
//!
//! Run: cargo run -p fs-render --features photometry --example photometric_stage
//! Add `-- --csv` to emit the direct-illuminance calculation grid as CSV.
//! Add `-- --paths` for reflected and multiply scattered spectral fixture paths
//! through a borrowed heterogeneous physical medium.
//! Add `-- --total-csv` for direct plus reflected/scattered surface lux through
//! that same physical medium, with estimated IID sampling error per grid point.
//! These declared fixture profiles and gray haze are examples, not measurements.

use std::error::Error;

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::camera::{Aperture, CameraProjection, PhysicalCamera};
use fs_render::charts::TriMesh;
use fs_render::instances::RigidTransform;
use fs_render::photometry::image::LuminanceJob;
use fs_render::photometry::medium::{HomogeneousMedium, PhaseFunction, ScatterSettings};
use fs_render::photometry::path::{
    IlluminanceSampling, PathIlluminanceJob, PathLuminanceJob, PathMedium,
};
use fs_render::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
use fs_render::photometry::{CalculationGrid, IlluminanceJob, PhotometricScene, TransportOptions};
use fs_render::spectral::LiftedSpectrum;
use fs_render::tracer::{
    Camera, DirectStrategy, Material, Primitive, Sampler, Scene, Settings, Shape,
};
use fs_render::volumes::VolumeGrid;
use fs_render::volumes::transport::scattering::{
    MassScatteringSpectrum, PhaseFunction as GridPhase, ScatteringGrid,
};
use fs_render::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (csv, paths, total_csv) = match args.as_slice() {
        [] => (false, false, false),
        [flag] if flag == "--csv" => (true, false, false),
        [flag] if flag == "--paths" => (false, true, false),
        [flag] if flag == "--total-csv" => (false, false, true),
        _ => return Err("usage: photometric_stage [--csv | --paths | --total-csv]".into()),
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
        run(&cx, csv, paths, total_csv)
    })
}

fn run(cx: &Cx<'_>, csv: bool, paths: bool, total_csv: bool) -> Result<(), Box<dyn Error>> {
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
    if !total_csv {
        println!("2 fixtures, 19,800 emitted lm; 13 x 9 stage grid at 0.5 m spacing");
        println!("Direct illuminance: mean {mean_lux:.3} lx, peak {peak_lux:.3} lx");
    }

    let camera = PhysicalCamera::try_look_at(
        Point3::new(7.0, -10.0, 5.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 0.0, 1.0),
        CameraProjection::try_vertical_fov(50.0_f64.to_radians())?,
        10.0,
        Aperture::try_circular(0.0)?, // Pinhole; no exposure or tone mapping.
    )?;
    if paths || total_csv {
        return run_paths(cx, Vec::from(surfaces), &sources, &camera, grid, total_csv);
    }
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

fn run_paths(
    cx: &Cx<'_>,
    mut primitives: Vec<Primitive>,
    sources: &[PhotometricSource],
    camera: &PhysicalCamera,
    grid: CalculationGrid,
    total_csv: bool,
) -> Result<(), Box<dyn Error>> {
    // A back wall receives reflected floor light through the same material law.
    primitives.push(Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [-5.0, 3.8, 0.0],
                [-5.0, 3.8, 4.0],
                [5.0, 3.8, 4.0],
                [5.0, 3.8, 0.0],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )),
        material: Material::Lambertian {
            reflectance: LiftedSpectrum { c: [0.0; 3] },
        },
        emission: None,
    });
    let scene = Scene {
        primitives,
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: camera.eye(),
            forward: camera.forward(),
            up: camera.up(),
            half_tan: camera.projection().vertical_half_tan(),
        },
    };
    // The original physical field is borrowed. rho*kappa has units m^-1;
    // these illustrative coefficients are not measured haze-machine data.
    let density: Vec<f64> = (0..32)
        .map(|i| if i % 4 < 2 { 0.25 } else { 1.0 })
        .collect();
    let absorption = AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([4, 4, 2], &density, [-5.0, -4.0, 0.0], [2.5, 2.0, 3.0]),
        MassAbsorptionSpectrum::table(&[(380.0, 0.01), (780.0, 0.002)])?,
        None,
        AbsorbingGridBudget {
            max_grid_cells: 32,
            max_steps_per_segment: 32,
        },
    )?;
    let medium = ScatteringGrid::try_new(
        cx,
        absorption,
        MassScatteringSpectrum::constant(0.025)?,
        GridPhase::henyey_greenstein(0.7)?,
    )?;
    if total_csv {
        let mut job = PathIlluminanceJob::for_grid(
            &scene,
            sources,
            PathMedium::Scattering(&medium),
            grid,
            TransportOptions {
                max_outputs: 4096,
                ..TransportOptions::default()
            },
            IlluminanceSampling {
                samples_per_probe: 1024,
                max_depth: 6,
                sampler: Sampler::Iid,
                seed: 0x5052_4556_495a,
            },
            cx,
        )?;
        while !job.is_complete() {
            job.advance(cx, 768)?;
        }
        println!(
            "x_m,y_m,z_m,direct_lux,continuous_lux,total_lux,estimated_standard_error_lux,samples,total_error_bound_lux"
        );
        for (index, estimate) in job.estimates().enumerate() {
            let estimate = estimate.ok_or("completed probe has no estimate")?;
            let point = grid.probe(index)?.position_m();
            let stderr = estimate
                .estimated_standard_error_lux
                .map(|value| format!("{value:.9}"))
                .unwrap_or_default();
            println!(
                "{},{},{},{:.9},{:.9},{:.9},{},{},unknown",
                point.x,
                point.y,
                point.z,
                estimate.direct_lux,
                estimate.continuous_lux,
                estimate.illuminance_lux,
                stderr,
                estimate.samples,
            );
        }
        return Ok(());
    }
    let settings = Settings {
        width: 16,
        height: 12,
        spp: 128,
        max_depth: 6,
        sampler: Sampler::Iid,
        strategy: DirectStrategy::Mis,
        seed: 0x5052_4556_495a,
    };
    let mut image = PathLuminanceJob::new(
        &scene,
        sources,
        PathMedium::Scattering(&medium),
        TransportOptions {
            max_outputs: 4096,
            ..TransportOptions::default()
        },
        settings,
        cx,
    )?;
    while !image.is_complete() {
        image.advance(cx, 768)?;
    }
    let peak = image
        .estimates()
        .flatten()
        .map(|p| p.luminance_cd_m2)
        .fold(0.0, f64::max);
    let stderr = image
        .estimates()
        .flatten()
        .filter_map(|p| p.estimated_standard_error_cd_m2)
        .fold(0.0, f64::max);
    println!("16 x 12 spectral path image, 128 samples/pixel, maximum 6 surface/volume vertices");
    println!("Reflections and multiple scattering: peak {peak:.6} cd/m^2");
    println!(
        "Largest estimated IID sampling standard error: {stderr:.6} cd/m^2 (not an error bound)"
    );
    println!(
        "Total error unknown: calibration, observer/discretization, and finite-depth remainder."
    );
    Ok(())
}
