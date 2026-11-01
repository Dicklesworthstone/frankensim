//! Infer a scalar fixture gain and check additional lux data without refitting.
//!
//! cargo run -p fs-render --features photometry --example photometric_flux_calibration
//!
//! ALL readings here are authored assumptions, not measured fixture data.
//! One 900 cd isotropic source gives 100 and 25 nominal lux at 3 and 6 m.
//! Synthetic observations describe a gain of 0.8 and a 2 lx background. The
//! independent necessary gain intersection is [0.795, 0.805]. A fourth-metre
//! observation is checked separately; a conflicting observation is also shown.
//! Stdout uses the existing forward flux CSV, with unknown total physical error.

#[cfg(feature = "photometry")]
mod enabled {
    use std::error::Error;
    use std::f64::consts::PI;
    use std::io::{self, Write};
    use fs_blake3::DomainHasher;
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    use fs_geom::{Point3, Vec3};
    use fs_ivl::Interval;
    use fs_render::instances::RigidTransform;
    use fs_render::photometry::medium::HomogeneousMedium;
    use fs_render::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
    use fs_render::photometry::source::validation::flux::{FluxBoundedIlluminance, FluxEnclosureJob};
    use fs_render::photometry::source::validation::flux::calibration::{FluxCalibrationJob, IlluminanceMeasurement};
    use fs_render::photometry::{IlluminanceProbe, PhotometricScene, TransportOptions};

    fn reading(z: f64, lo: f64, hi: f64) -> Result<IlluminanceMeasurement, Box<dyn Error>> {
        let mut evidence = DomainHasher::new("org.frankensim.example.authored-lux-observation.v1");
        for value in [z, lo, hi, 2.0] { evidence.update(&value.to_bits().to_le_bytes()); }
        Ok(IlluminanceMeasurement::try_new(
            IlluminanceProbe::try_new(Point3::new(0.0, 0.0, z), Vec3::new(0.0, 0.0, -1.0))?,
            Interval::new(lo, hi), Interval::point(2.0), evidence.finalize(),
        )?)
    }

    fn calculate(cx: &Cx<'_>) -> Result<(), Box<dyn Error>> {
        let sources = [PhotometricSource::try_new(RigidTransform::identity(), 3600.0 * PI,
            1.0, AngularDistribution::isotropic(), SourceSpectrum::equal_energy())?];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), TransportOptions::default(), cx)?;
        let data = [reading(3.0, 81.0, 83.0)?, reading(6.0, 21.875, 22.125)?];
        let scope = DomainHasher::new("org.frankensim.example.authored-common-flux-scope.v1").finalize();
        let mut calibration = FluxCalibrationJob::new(&scene, &data, scope, cx)?;
        while !calibration.is_complete() { calibration.advance(cx, 1)?; }
        let result = calibration.result().ok_or_else(|| io::Error::other("incomplete calibration"))?;
        let gain = result.flux_scale_bound().ok_or_else(|| io::Error::other("no finite compatible gain"))?;
        eprintln!("AUTHORED observations; not measured calibration; all physical totals remain unknown");
        eprintln!("analytical gain range=[0.795,0.805]; inferred={:?}; status={:?}", gain.scale(), result.status());
        for observation in [reading(4.0, 46.9, 47.1)?, reading(4.0, 52.0, 53.0)?] {
            let check = calibration.check_measurement(cx, &observation)?;
            eprintln!("additional observed_lux={:?} predicted_with_background_lux={:?} compatibility={:?}",
                check.measured_lux, check.predicted_measured_lux, check.compatibility);
        }
        let probes = [data[0].probe(), reading(4.0, 46.9, 47.1)?.probe(), data[1].probe()];
        let bounds = [Some(gain)];
        let mut forward = FluxEnclosureJob::new(&scene, &probes, &bounds, cx)?;
        forward.advance(cx, usize::MAX)?;
        let stdout = io::stdout();
        let mut writer = stdout.lock();
        FluxBoundedIlluminance::write_csv_header(&mut writer)?;
        for (index, report) in forward.results().iter().enumerate() { report.write_csv_row(index, &mut writer)?; }
        writer.flush()?;
        Ok(())
    }

    pub fn run() -> Result<(), Box<dyn Error>> {
        let gate = CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| calculate(&Cx::new(&gate, arena,
            StreamKey { seed: 173, kernel_id: 0x4643_414c, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic)))
    }
}

#[cfg(feature = "photometry")]
fn main() -> Result<(), Box<dyn std::error::Error>> { enabled::run() }

#[cfg(not(feature = "photometry"))]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    Err(std::io::Error::other("photometric_flux_calibration requires --features photometry").into())
}
