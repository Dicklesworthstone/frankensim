//! Independently calibrate two overlapping emitters from controlled lux readings.
//!
//! cargo run -p fs-render --features photometry --example photometric_emitter_calibration
//!
//! ALL readings and drive levels are authored assumptions, not measured data.
//! Reference contributions are 100 and 50 lx at the central receiver. Outputs
//! [1,1/2] and [1/2,1] distinguish gains 0.8 and 1.2 with a 2 lx background.
//! Independent inversion gives g0 in [.798,.802], g1 in [1.196,1.204]. The
//! all-on observation is held out, never fitted. Marginal gain ranges lose joint
//! correlations and may widen its prediction. No total physical accuracy claim.
//! Stdout is original forward-flux CSV at full reference output, WITHOUT the
//! meter observations' 2 lx background. Stderr gives inverse/holdout diagnostics.

#[cfg(feature = "photometry")]
mod enabled {
    use std::error::Error;
    use std::io::{self, Write};
    use fs_blake3::DomainHasher;
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    use fs_geom::{Point3, Vec3};
    use fs_ivl::Interval;
    use fs_render::instances::RigidTransform;
    use fs_render::photometry::{IlluminanceProbe, PhotometricScene, TransportOptions};
    use fs_render::photometry::medium::HomogeneousMedium;
    use fs_render::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
    use fs_render::photometry::source::validation::flux::{FluxBoundedIlluminance, FluxEnclosureJob};
    use fs_render::photometry::source::validation::flux::calibration::IlluminanceMeasurement;
    use fs_render::photometry::source::validation::flux::calibration::emitters::{
        EmitterCalibrationBudget, EmitterCalibrationJob, EmitterCalibrationObservation,
    };

    fn probe(x: f64) -> Result<IlluminanceProbe, Box<dyn Error>> {
        Ok(IlluminanceProbe::try_new(Point3::new(x, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0))?)
    }
    fn reading(lo: f64, hi: f64, output: [f64; 2]) -> Result<IlluminanceMeasurement, Box<dyn Error>> {
        let mut hash = DomainHasher::new("org.frankensim.example.authored-emitter-reading.v1");
        for value in [lo, hi, 2.0, output[0], output[1]] { hash.update(&value.to_bits().to_le_bytes()); }
        Ok(IlluminanceMeasurement::try_new(probe(0.0)?, Interval::new(lo, hi),
            Interval::point(2.0), hash.finalize())?)
    }

    fn calculate(cx: &Cx<'_>) -> Result<(), Box<dyn Error>> {
        let source = |height, intensity| -> Result<PhotometricSource, Box<dyn Error>> {
            Ok(PhotometricSource::try_new(
                RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, height])?,
                4.0 * std::f64::consts::PI * intensity, 1.0,
                AngularDistribution::isotropic(), SourceSpectrum::equal_energy())?)
        };
        let sources = [source(3.0, 900.0)?, source(6.0, 1800.0)?];
        let scene = PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(),
            TransportOptions::default(), cx)?;
        let data = [
            EmitterCalibrationObservation { measurement: reading(111.9, 112.1, [1.0, 0.5])?, relative_output: &[1.0, 0.5] },
            EmitterCalibrationObservation { measurement: reading(101.9, 102.1, [0.5, 1.0])?, relative_output: &[0.5, 1.0] },
        ];
        let scope = DomainHasher::new("org.frankensim.example.authored-independent-emitter-scope.v1").finalize();
        let mut job = EmitterCalibrationJob::new(&scene, &data, scope,
            EmitterCalibrationBudget { max_sweeps: 64, ..Default::default() }, cx)?;
        // Mix preparation/refinement chunks; no callback or sample depends on them.
        while !job.is_complete() { job.advance(cx, 3)?; }
        let fit = job.result().ok_or_else(|| io::Error::other("no completed calibration"))?;
        let gains = fit.scale_enclosures().ok_or_else(|| io::Error::other("inconsistent observations"))?;
        eprintln!("authored assumptions only; not measured calibration; termination={:?} sweeps={}", fit.termination, fit.completed_sweeps);
        eprintln!("original emitter gains={gains:?}; calibration_identity={}", fit.evidence_identity);
        eprintln!("independent reference: g0=(u-v/2)/75 in [.798,.802]; g1=(2v-u)/75 in [1.196,1.204]");
        for (lo, hi) in [(141.9, 142.1), (160.0, 161.0)] {
            let observation = EmitterCalibrationObservation {
                measurement: reading(lo, hi, [1.0, 1.0])?, relative_output: &[1.0, 1.0] };
            let check = job.check_observation(cx, observation)?;
            eprintln!("held_out={:?} predicted_with_background={:?} compatibility={:?}",
                check.measured_lux, check.predicted_measured_lux, check.compatibility);
        }
        let bounds = [fit.flux_scale_bound(0), fit.flux_scale_bound(1)];
        if bounds.iter().any(Option::is_none) {
            return Err(io::Error::other("an emitter remains unbounded").into());
        }
        let probes = [probe(0.0)?, probe(1.0)?, probe(2.0)?];
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
            StreamKey { seed: 173, kernel_id: 0x454d_4341, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic)))
    }
}

#[cfg(feature = "photometry")]
fn main() -> Result<(), Box<dyn std::error::Error>> { enabled::run() }

#[cfg(not(feature = "photometry"))]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    Err(std::io::Error::other("photometric_emitter_calibration requires --features photometry").into())
}
