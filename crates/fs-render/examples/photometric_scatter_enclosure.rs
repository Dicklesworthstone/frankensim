//! Headless numerical bounds on homogeneous single-scattered haze glow.
//!
//! cargo run -p fs-render --features photometry --example photometric_scatter_enclosure
//! Add `-- --shadow`, `-- --coarse`, or `-- --shadow --coarse`.
//! Authored analytical inputs, NOT measured fixture/medium calibration. Stdout
//! is component-scoped CSV. Stderr holds independent integral diagnostics.
//! Surface radiance and higher scattering orders remain outside these bounds.

#[cfg(feature = "photometry")]
mod demo {
    use std::{error::Error, io::{self, Write}};
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    use fs_geom::Point3;
    use fs_render::charts::TriMesh;
    use fs_render::instances::RigidTransform;
    use fs_render::photometry::{PhotometricScene, TransportOptions};
    use fs_render::photometry::medium::{HomogeneousMedium, PhaseFunction};
    use fs_render::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
    use fs_render::photometry::source::validation::scattering::{
        EnclosedSingleScatter, HomogeneousScatterReference, ScatterEnclosureBudget,
        ScatterEnclosureJob, ScatterSegment,
    };
    use fs_render::spectral::SpectralReflectance;
    use fs_render::tracer::{Material, Primitive, Shape};

    fn calculate(cx: &Cx<'_>, shadow: bool, coarse: bool) -> Result<(), Box<dyn Error>> {
        // A nominal 100 cd isotropic fixture, two metres to the side of the camera.
        let sources = [PhotometricSource::try_new(
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [2.0, 0.0, 0.0])?,
            400.0 * std::f64::consts::PI,
            1.0,
            AngularDistribution::isotropic(),
            SourceSpectrum::equal_energy(),
        )?];
        let geometry = if shadow {
            // At x=1 the source connection to camera-ray position z=s has z=s/2.
            // Hence this strip blocks .4 <= s <= 1.2, independently of a ray tracer.
            vec![Primitive {
                shape: Shape::Mesh(TriMesh::new(
                    vec![[1.0, -1.0, 0.2], [1.0, 1.0, 0.2],
                         [1.0, 1.0, 0.6], [1.0, -1.0, 0.6]],
                    vec![[0, 1, 2], [0, 2, 3]],
                )),
                material: Material::SpectralDiffuse {
                    reflectance: SpectralReflectance::constant(0.0)?,
                },
                emission: None,
            }]
        } else { Vec::new() };
        let medium = HomogeneousMedium::try_new(
            0.1, 0.2, PhaseFunction::HenyeyGreenstein { g: 0.7 },
        )?;
        let scene = PhotometricScene::try_new(
            &geometry, &sources, medium, TransportOptions::default(), cx,
        )?;
        let reference = HomogeneousScatterReference::new(&scene, cx)?;
        let segments = [
            ScatterSegment::try_new(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 1.0))?,
            ScatterSegment::try_new(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 2.0))?,
            ScatterSegment::try_new(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 3.0))?,
        ];
        let budget = ScatterEnclosureBudget {
            max_evaluations: if coarse { 1 } else { 16_384 },
            ..ScatterEnclosureBudget::default()
        };
        let mut job = ScatterEnclosureJob::new(&reference, &segments, budget, cx)?;
        while !job.is_complete() { job.advance(cx, 1)?; }
        // Independently evaluated at 80 decimal digits in physical s coordinates:
        // r=sqrt(4+s²), mu=-s/r, integral 20*p_HG(mu)*exp(-.3*(r+s))/r² ds.
        // Shadow mode integrates only outside the analytically projected strip.
        // Rounded diagnostics are NOT inputs to the enclosure or accuracy proof.
        let diagnostic = if shadow {
            [0.020_026_788_327_591_796, 0.028_145_024_504_812_26, 0.031_789_885_537_136_57]
        } else {
            [0.037_406_117_541_195_62, 0.049_078_971_532_979_35, 0.052_723_832_565_303_65]
        };
        eprintln!("authored gray single-scattering model; absorption=0.1/m, scattering=0.2/m, HG g=0.7");
        eprintln!("shadow={shadow}; coarse={coarse}; multiple scattering, surface light, calibration and total physical error unknown");
        let stdout = io::stdout();
        let mut writer = stdout.lock();
        EnclosedSingleScatter::write_csv_header(&mut writer)?;
        for (index, result) in job.results().iter().enumerate() {
            eprintln!("segment_end_m={} independent_integral_f64_cd_m2={}", index + 1, diagnostic[index]);
            result.write_csv_row(index, &mut writer)?;
        }
        writer.flush()?;
        Ok(())
    }

    pub fn run() -> Result<(), Box<dyn Error>> {
        let (mut shadow, mut coarse) = (false, false);
        for arg in std::env::args().skip(1) {
            match arg.as_str() {
                "--shadow" if !shadow => shadow = true,
                "--coarse" if !coarse => coarse = true,
                _ => return Err(io::Error::new(io::ErrorKind::InvalidInput,
                    "usage: photometric_scatter_enclosure [--shadow] [--coarse]").into()),
            }
        }
        let gate = CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(&gate, arena,
                StreamKey { seed: 173, kernel_id: 0x5343_4558, tile: 0, iteration: 0 },
                Budget::INFINITE, ExecMode::Deterministic);
            calculate(&cx, shadow, coarse)
        })
    }
}

#[cfg(feature = "photometry")]
fn main() -> Result<(), Box<dyn std::error::Error>> { demo::run() }

#[cfg(not(feature = "photometry"))]
fn main() {
    eprintln!("photometric_scatter_enclosure requires --features photometry");
    std::process::exit(2);
}
