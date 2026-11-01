//! Deterministic uniform-panel reference; no GPU, random sampling or scene tracing.
//!
//! cargo run -p fs-render --features photometry --example photometric_area_reference
//!
//! A 1 m x 1 m downward Lambertian panel carries 1000 declared lumens. The
//! receiving grid is in vacuum with no occluders or reflected contribution.
//! Stdout is scoped numerical-bound CSV; stderr gives an independently derived
//! axial rectangle integral. Inputs are authored, not measured calibration.

#[cfg(feature = "photometry")]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    use fs_geom::{Point3, Vec3};
    use fs_render::charts::TriMesh;
    use fs_render::photometry::area::{EmissionSidedness, PhotometricAreaEmitter};
    use fs_render::photometry::area::reference::{
        AreaReferenceBudget, AreaReferenceIlluminance, UnoccludedAreaJob,
        UnoccludedAreaReference,
    };
    use fs_render::photometry::source::SourceSpectrum;
    use fs_render::photometry::{CalculationGrid, IlluminanceProbe};
    use fs_render::spectral::SpectralReflectance;
    use fs_render::tracer::{Material, Primitive, Shape};
    use std::io::{self, Write};

    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| -> Result<(), Box<dyn std::error::Error>> {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 173, kernel_id: 0x4152_4546, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let primitives = [Primitive {
            shape: Shape::Mesh(TriMesh::new(
                vec![[-0.5,-0.5,1.0],[0.5,-0.5,1.0],[0.5,0.5,1.0],[-0.5,0.5,1.0]],
                vec![[0,2,1],[0,3,2]],
            )),
            material: Material::SpectralDiffuse {
                reflectance: SpectralReflectance::constant(0.0)?,
            },
            emission: None,
        }];
        let sources = [PhotometricAreaEmitter::try_new(
            0, 1000.0, 1.0, SourceSpectrum::equal_energy(), EmissionSidedness::Front,
        )?];
        let reference = UnoccludedAreaReference::new(
            &primitives, &sources, AreaReferenceBudget::default(), &cx,
        )?;
        let grid = CalculationGrid::try_new(
            Point3::new(-1.0,-1.0,0.0), Vec3::new(0.5,0.0,0.0),
            Vec3::new(0.0,0.5,0.0), 5, 5,
        )?;
        let mut job = UnoccludedAreaJob::for_grid(&reference, grid, &cx)?;
        while !job.is_complete() { job.advance(&cx, 3)?; }
        // Direct area integral, independently from the spherical edge formula.
        // This f64 diagnostic assumes unit photopic spectrum normalization;
        // production enclosures retain the actual normalized spectral bins.
        let q = 0.5 / 1.25_f64.sqrt();
        let analytic = 1000.0 / std::f64::consts::PI * 4.0 * q * q.atan();
        eprintln!("scope=unoccluded-vacuum-direct-area; central_analytic_f64_lux={analytic}");
        eprintln!("authored input; fixture calibration and total physical error unknown");
        let stdout = io::stdout();
        let mut writer = stdout.lock();
        AreaReferenceIlluminance::write_csv_header(&mut writer)?;
        for (index, value) in job.results().iter().enumerate() {
            value.write_csv_row(index, &mut writer)?;
        }
        // A vertical receiver sees only half the panel: horizon clipping is
        // necessary even though the entire source is unobstructed.
        let vertical = IlluminanceProbe::try_new(
            Point3::new(0.0,0.0,0.0), Vec3::new(1.0,0.0,0.0),
        )?;
        reference.illuminance(&cx, &vertical)?.write_csv_row(25, &mut writer)?;
        writer.flush()?;
        Ok(())
    })
}

#[cfg(not(feature = "photometry"))]
fn main() {
    eprintln!("enable --features photometry to run this conditional area reference");
}
