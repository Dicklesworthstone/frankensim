//! Headless bounds on haze glow in native spectrally selective physical cells.
//!
//! cargo run -p fs-render --features photometry --example photometric_grid_scatter
//! Options: --shadow, --coarse, --volume-tail, and one of --blue / --red.
//! Authored numerical inputs, not measured fixture/medium calibration. The
//! single-scattering result excludes surface light, thermal light and higher
//! scattering orders. Stdout is CSV; independent diagnostics go to stderr.

#[cfg(feature = "photometry")]
mod example {
    use std::{error::Error, io::{self, Write}, sync::Arc};
    use core::f64::consts::PI;
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    use fs_geom::Point3;
    use fs_render::charts::TriMesh;
    use fs_render::instances::RigidTransform;
    use fs_render::photometry::{PhotometricScene, TransportOptions, transport::PhotometricMedium};
    use fs_render::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum,
        SPECTRAL_BINS, LUMINOUS_EFFICACY_LM_W, photopic_weight_nm};
    use fs_render::photometry::source::validation::scattering::{SingleScatterReference,
        ScatterEnclosureJob, ScatterEnclosureBudget, ScatterSegment, EnclosedSingleScatter};
    use fs_render::spectral::SpectralReflectance;
    use fs_render::tracer::{Material, Primitive, Shape};
    use fs_render::volumes::VolumeGrid;
    use fs_render::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget,
        CellPalette, CellPaletteBudget, MassAbsorptionSpectrum};
    use fs_render::volumes::transport::scattering::{ScatteringGrid, ScatteringOptics,
        MassScatteringSpectrum, PhaseFunction};

    fn calculate(cx: &Cx<'_>, shadow: bool, coarse: bool, blue: f64, volume_tail: bool) -> Result<(), Box<dyn Error>> {
        let density = [1.0, 2.0, 0.5];
        let cap = CellPaletteBudget { max_models: 2, max_cells: 3 };
        let absorption = CellPalette::try_new(cx, vec![
            MassAbsorptionSpectrum::constant(0.05)?,
            MassAbsorptionSpectrum::table(&[(380.0, 0.01), (780.0, 0.2)])?,
        ], Arc::from([0, 1, 0]), cap)?;
        let absorption = AbsorbingGrid::try_new_indexed(cx,
            VolumeGrid::new([1, 1, 3], &density, [-2.0, -0.5, 0.0], [4.0, 1.0, 1.0]),
            absorption, None, AbsorbingGridBudget { max_grid_cells: 3, max_steps_per_segment: 64 })?;
        let scattering = CellPalette::try_new(cx, vec![
            ScatteringOptics { scattering: MassScatteringSpectrum::constant(0.1)?, phase: PhaseFunction::Isotropic },
            ScatteringOptics { scattering: MassScatteringSpectrum::rayleigh(400.0, 0.2)?,
                phase: PhaseFunction::spectral_table(&[(380.0, PhaseFunction::Isotropic), (780.0, PhaseFunction::Rayleigh)])? },
        ], Arc::from([0, 1, 0]), cap)?;
        let medium = ScatteringGrid::try_new_indexed(cx, absorption, scattering)?;
        let mut bins = [0.0; SPECTRAL_BINS];
        bins[4] = blue / (LUMINOUS_EFFICACY_LM_W * photopic_weight_nm(4));
        bins[44] = (1.0 - blue) / (LUMINOUS_EFFICACY_LM_W * photopic_weight_nm(44));
        let sources = [PhotometricSource::try_new(
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [if shadow { 1.0 } else { 0.0 }, 0.0, 4.0])?,
            400.0 * PI, 1.0, AngularDistribution::isotropic(), SourceSpectrum::try_from_bins(bins)?)?];
        let geometry = if shadow {
            vec![Primitive { shape: Shape::Mesh(TriMesh::new(
                vec![[0.25, -0.5, 2.0], [0.375, -0.5, 2.0], [0.375, 0.5, 2.0], [0.25, 0.5, 2.0]],
                vec![[0, 1, 2], [0, 2, 3]])),
                material: Material::SpectralDiffuse { reflectance: SpectralReflectance::constant(0.0)? }, emission: None }]
        } else { Vec::new() };
        let scene = PhotometricScene::try_with_medium(&geometry, &sources,
            PhotometricMedium::Scattering(&medium), TransportOptions::default(), cx)?;
        let reference = SingleScatterReference::new(&scene, cx)?;
        let reference = if volume_tail { reference.with_volume_remainder(3, cx)? } else { reference };
        let segments = [ScatterSegment::try_new(Point3::new(0.0, 0.0, -1.0), Point3::new(0.0, 0.0, 3.5))?];
        let limits = ScatterEnclosureBudget { absolute_error_cd_m2: 1e-3,
            max_evaluations: if coarse { 1 } else { 16_384 }, max_depth: 40 };
        let mut job = ScatterEnclosureJob::new(&reference, &segments, limits, cx)?;
        while !job.is_complete() { job.advance(cx, 1)?; }
        // Independent physical-z integrals at 90 digits. The unshadowed source
        // is collinear: total two-leg optical depth is constant per wavelength.
        // For the off-axis source, project the screen to z in [4/5,4/3] and
        // integrate separately on native layer and shadow boundaries. These
        // rounded diagnostics are not the proof used by the production bound.
        let (blue_reference, red_reference) = if shadow {
            (0.258_804_792_582_718_86, 0.146_501_140_091_266_15)
        } else {
            (0.413_502_705_572_333, 0.232_716_744_046_849_13)
        };
        let analytic = blue * blue_reference + (1.0 - blue) * red_reference;
        eprintln!("authored native-grid single scattering; density and optical/phase maps are shared with physical transport");
        eprintln!("blue_fraction={blue} shadow={shadow} independent_reference_cd_m2={analytic}");
        eprintln!("calibration, discretization, surface light and total physical error remain unknown");
        match reference.volume_remainder() {
            Some(tail) => eprintln!("volume-only orders >=2: upper_cd_m2={} global_first_order_cd_m2={} contraction={} source_distance_lower_m={:?} grid_diameter_upper_m={} proof_version={}",
                tail.upper_cd_m2(), tail.first_order_upper_cd_m2(), tail.max_contraction(),
                tail.min_source_distance_m(), tail.grid_diameter_m(), tail.semantics_version()),
            None => eprintln!("higher-order volume bound unknown or not requested"),
        }
        let stdout = io::stdout();
        let mut writer = stdout.lock();
        EnclosedSingleScatter::write_csv_header(&mut writer)?;
        for (index, result) in job.results().iter().enumerate() { result.write_csv_row(index, &mut writer)?; }
        writer.flush()?;
        Ok(())
    }

    pub fn run() -> Result<(), Box<dyn Error>> {
        let (mut shadow, mut coarse, mut color, mut volume_tail) = (false, false, None, false);
        for arg in std::env::args().skip(1) {
            match arg.as_str() {
                "--shadow" if !shadow => shadow = true,
                "--coarse" if !coarse => coarse = true,
                "--volume-tail" if !volume_tail => volume_tail = true,
                "--blue" if color.is_none() => color = Some(1.0),
                "--red" if color.is_none() => color = Some(0.0),
                _ => return Err(io::Error::new(io::ErrorKind::InvalidInput,
                    "usage: photometric_grid_scatter [--shadow] [--coarse] [--volume-tail] [--blue | --red]").into()),
            }
        }
        let gate = CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(&gate, arena, StreamKey { seed: 173, kernel_id: 0x4753_4558, tile: 0, iteration: 0 },
                Budget::INFINITE, ExecMode::Deterministic);
            calculate(&cx, shadow, coarse, color.unwrap_or(0.75), volume_tail)
        })
    }
}

#[cfg(feature = "photometry")]
fn main() -> Result<(), Box<dyn std::error::Error>> { example::run() }
#[cfg(not(feature = "photometry"))]
fn main() {
    eprintln!("photometric_grid_scatter requires --features photometry");
    std::process::exit(2);
}
