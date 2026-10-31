//! Headless illumination reference with a budget for omitted reflection orders.
//!
//! Run: cargo run -p fs-render --features photometry --example photometric_error_budget
//! Stdout is per-probe CSV; stderr names the independent analytic reference.
//! Add `-- --direct` for direct-lux arithmetic enclosures at the ten CIE-derived
//! diffuse point-source positions, or `-- --direct-haze` for an authored gray
//! extinction extension. `--direct-shadow` and `--direct-shadow-haze` add a
//! native opaque triangular screen and a receiving floor, including two probes
//! on the screen's shadow edge. Their unresolved visibility widens the bound.
//! These use the original source, geometry, material and medium owners.
//!
//! An inward-facing 2 m cube emits uniform Le=8 cd/m² and reflects rho=0.5.
//! Its ideal vacuum equilibrium is L=Le/(1-rho), E=pi*L for every interior
//! probe orientation. This follows by summing the independent geometric series
//! Le*(1+rho+rho²+...), not from a second rendering implementation.
//!
//! The 0.01 lx budget covers only the admitted model's omitted reflection
//! orders. Fixed-count IID sampling error is estimated separately; arithmetic,
//! spectral, geometry, physical-model and calibration errors remain unknown.
//! The optional direct modes instead enclose the retained numerical model,
//! including supported mesh visibility. Source construction, physical geometry
//! fidelity, calibration and total physical error remain unknown. Haze and
//! shadows are authored extensions, not additional CIE validation cases.
//! No mode establishes full CIE 171 conformance or measured-fixture accuracy.

use std::error::Error;
use std::f64::consts::{FRAC_PI_2, PI};
use std::io::{self, Write};

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::charts::TriMesh;
use fs_render::instances::RigidTransform;
use fs_render::photometry::area::{
    EmissionSidedness, PhotometricAreaEmitter, PhotometricEmitters,
};
use fs_render::photometry::medium::{HomogeneousMedium, PhaseFunction};
use fs_render::photometry::path::{
    IlluminanceSampling, PathErrorBudget, PathIlluminanceJob, PathMedium,
};
use fs_render::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
use fs_render::photometry::source::validation::{EnclosedIlluminance, EnclosedIlluminanceJob};
use fs_render::photometry::{IlluminanceProbe, PhotometricScene, TransportOptions};
use fs_render::spectral::SpectralReflectance;
use fs_render::tracer::{Camera, Material, Primitive, Sampler, Scene, Shape};

const SOURCE_LUMINANCE: f64 = 8.0;
const REFLECTANCE: f64 = 0.5;
const TAIL_BUDGET_LUX: f64 = 0.01;
const DEPTH_CAP: u32 = 64;
const SAMPLES: u32 = 1024;
const SEED: u64 = 173;

fn scene() -> Result<Scene, Box<dyn Error>> {
    // Every triangle faces the interior; six 4 m² faces give exactly 24 m².
    let mesh = TriMesh::new(
        vec![
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ],
        vec![
            [0, 1, 2], [0, 2, 3],
            [4, 6, 5], [4, 7, 6],
            [0, 5, 1], [0, 4, 5],
            [3, 2, 6], [3, 6, 7],
            [0, 3, 7], [0, 7, 4],
            [1, 6, 2], [1, 5, 6],
        ],
    );
    Ok(Scene {
        primitives: vec![Primitive {
            shape: Shape::Mesh(mesh),
            material: Material::SpectralDiffuse {
                reflectance: SpectralReflectance::constant(REFLECTANCE)?,
            },
            emission: None,
        }],
        lights: Vec::new(),
        environment: None,
        // Probe jobs do not use the camera.
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 0.0),
            forward: Vec3::new(0.0, 0.0, 1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.0,
        },
    })
}

fn calculate(cx: &Cx<'_>) -> Result<(), Box<dyn Error>> {
    let scene = scene()?;
    let surfaces = [PhotometricAreaEmitter::try_new(
        0,
        24.0 * PI * SOURCE_LUMINANCE,
        1.0,
        SourceSpectrum::equal_energy(),
        EmissionSidedness::Front,
    )?];
    let emitters = PhotometricEmitters {
        points: &[],
        surfaces: &surfaces,
    };
    let probes = [
        IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0))?,
        IlluminanceProbe::try_new(Point3::new(0.2, -0.3, 0.1), Vec3::new(1.0, 0.0, 0.0))?,
    ];
    let mut sampling = IlluminanceSampling {
        samples_per_probe: SAMPLES,
        max_depth: 1,
        sampler: Sampler::Iid,
        seed: SEED,
    };
    let initial = PathIlluminanceJob::new_with_emitters(
        &scene,
        emitters,
        PathMedium::Vacuum,
        &probes,
        TransportOptions::default(),
        sampling,
        cx,
    )?;
    let envelope = initial.diffuse_vacuum_remainder().ok_or_else(|| {
        io::Error::other("scene has no admitted diffuse-vacuum remainder")
    })?;
    let plan = envelope
        .depth_for_illuminance_budget(TAIL_BUDGET_LUX, DEPTH_CAP)
        .ok_or_else(|| io::Error::other("depth cap cannot establish the model-tail budget"))?;
    sampling.max_depth = plan.max_depth();
    let mut job = PathIlluminanceJob::new_with_emitters(
        &scene,
        emitters,
        PathMedium::Vacuum,
        &probes,
        TransportOptions::default(),
        sampling,
        cx,
    )?;
    let actual = job
        .diffuse_vacuum_remainder()
        .ok_or_else(|| io::Error::other("reconstructed job lost its model envelope"))?;
    if actual.illuminance_bound_lux() > TAIL_BUDGET_LUX {
        return Err(io::Error::other("reconstructed job exceeds the model-tail budget").into());
    }
    // Progressive partitions retain the original sample addresses and sums.
    while !job.is_complete() {
        job.advance(cx, 127)?;
    }
    eprintln!(
        "analytic_all_order_lux={} depth={} omitted_model_bound_lux={} samples={} seed={SEED}",
        PI * SOURCE_LUMINANCE / (1.0 - REFLECTANCE),
        actual.max_depth(),
        actual.illuminance_bound_lux(),
        SAMPLES,
    );
    eprintln!("model-tail budget only; sampling is estimated; total physical error is unknown");
    let stdout = io::stdout();
    let mut writer = stdout.lock();
    PathErrorBudget::write_csv_header(&mut writer)?;
    for index in 0..job.probe_count() {
        job.error_budget_lux(index)
            .ok_or_else(|| io::Error::other("completed probe has no estimate"))?
            .write_csv_row(&mut writer)?;
    }
    writer.flush()?;
    Ok(())
}

fn shadow_geometry() -> Result<Vec<Primitive>, Box<dyn Error>> {
    let material = Material::SpectralDiffuse {
        reflectance: SpectralReflectance::constant(0.3)?,
    };
    Ok(vec![
        Primitive {
            shape: Shape::Mesh(TriMesh::new(
                vec![[-0.25, -0.25, 1.5], [0.75, -0.25, 1.5], [-0.25, 0.75, 1.5]],
                vec![[0, 1, 2]],
            )),
            material,
            emission: None,
        },
        Primitive {
            // The floor passes through the probes; receiver endpoint exclusion
            // must prevent self-blocking without discarding the separate screen.
            shape: Shape::Mesh(TriMesh::new(
                vec![[-2.0, -2.0, 0.0], [2.0, -2.0, 0.0], [2.0, 2.0, 0.0], [-2.0, 2.0, 0.0]],
                vec![[0, 1, 2], [0, 2, 3]],
            )),
            material,
            emission: None,
        },
    ])
}

fn direct_reference(cx: &Cx<'_>, haze: bool, shadow: bool) -> Result<(), Box<dyn Error>> {
    // Factual setup and provenance are retained in
    // tests/photometric_path/cie171_analytical.rs. No rounded table is copied.
    let positions: [[f64; 2]; 10] = [
        [0.0, 0.0], [0.5, 0.0], [1.0, 0.0], [1.5, 0.0], [0.5, 0.5],
        [1.0, 0.5], [1.5, 0.5], [1.0, 1.0], [1.5, 1.0], [1.5, 1.5],
    ];
    let sources = [PhotometricSource::try_new(
        RigidTransform::try_new([1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 3.0])?,
        1000.0 * PI,
        1.0,
        AngularDistribution::cosine_power(1.0, FRAC_PI_2)?,
        SourceSpectrum::equal_energy(),
    )?];
    let absorption = if haze { 0.07 } else { 0.0 };
    let scattering = if haze { 0.13 } else { 0.0 };
    let medium = HomogeneousMedium::try_new(absorption, scattering, PhaseFunction::Isotropic)?;
    let geometry = if shadow { shadow_geometry()? } else { Vec::new() };
    let scene = PhotometricScene::try_new(
        &geometry, &sources, medium, TransportOptions::default(), cx,
    )?;
    let probes: Vec<_> = positions.iter().map(|[x, y]| {
        IlluminanceProbe::try_new(Point3::new(*x, *y, 0.0), Vec3::new(0.0, 0.0, 1.0))
    }).collect::<Result<_, _>>()?;
    let mut job = EnclosedIlluminanceJob::new(&scene, &probes, cx)?;
    while !job.is_complete() {
        job.advance(cx, 3)?;
    }
    eprintln!("direct unscattered lux; source construction, geometry fidelity, calibration and total physical error unknown");
    eprintln!("unblocked reference: 9000/(9+x*x+y*y)^2 * exp(-(absorption+scattering)*sqrt(9+x*x+y*y))");
    eprintln!("gray absorption_per_m={absorption} scattering_per_m={scattering} shadow={shadow}; haze and shadows are authored extensions, not CIE cases");
    if shadow {
        eprintln!("screen intersection is (x/2,y/2,1.5); blocked when x>=-0.5, y>=-0.5, x+y<=1; the two edge probes may retain both visibility branches");
    }
    let stdout = io::stdout();
    let mut writer = stdout.lock();
    EnclosedIlluminance::write_csv_header(&mut writer)?;
    for (index, (report, [x, y])) in job.results().iter().zip(positions).enumerate() {
        let distance_squared = 9.0 + x * x + y * y;
        // These independent formula/plane tests are diagnostics, NOT the
        // source of the production interval or a measured-fixture reference.
        let blocked = shadow && x >= -0.5 && y >= -0.5 && x + y <= 1.0;
        let reference = if blocked {
            0.0
        } else {
            9000.0 / (distance_squared * distance_squared)
                * (-(absorption + scattering) * distance_squared.sqrt()).exp()
        };
        eprintln!("probe={index} x_m={x} y_m={y} analytic_blocked={blocked} analytic_reference_f64_lux={reference}");
        report.write_csv_row(index, &mut writer)?;
    }
    writer.flush()?;
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let mode = match args.as_slice() {
        [] => 0,
        [arg] if arg == "--direct" => 1,
        [arg] if arg == "--direct-haze" => 2,
        [arg] if arg == "--direct-shadow" => 3,
        [arg] if arg == "--direct-shadow-haze" => 4,
        _ => return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "usage: photometric_error_budget [--direct | --direct-haze | --direct-shadow | --direct-shadow-haze]",
        ).into()),
    };
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: SEED,
                kernel_id: 0x5048_4552,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        match mode {
            1 => direct_reference(&cx, false, false),
            2 => direct_reference(&cx, true, false),
            3 => direct_reference(&cx, false, true),
            4 => direct_reference(&cx, true, true),
            _ => calculate(&cx),
        }
    })
}
