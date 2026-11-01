//! Finite-panel illumination through spatially and spectrally varying native haze.
//!
//! cargo run -p fs-render --features photometry --example photometric_area_grid
//! Optional flags: --shadow, --coarse, and one of --blue or --red.
//! Default source: 75/25 photopic power at 402.5/602.5 nm. All inputs are authored
//! assumptions, not measured fixture/medium calibration. The CSV bounds DIRECT
//! unscattered light only; thermal, reflected and in-scattered light are omitted.

#[cfg(feature = "photometry")]
mod demo {
    use std::{error::Error, io::{self, Write}, sync::Arc};
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    use fs_geom::{Point3, Vec3};
    use fs_render::charts::TriMesh;
    use fs_render::photometry::{IlluminanceProbe, observer};
    use fs_render::photometry::area::{EmissionSidedness, PhotometricAreaEmitter};
    use fs_render::photometry::area::reference::AreaReferenceBudget;
    use fs_render::photometry::area::reference::occlusion::{
        AreaOcclusionBudget, OccludedAreaIlluminance, OccludedAreaJob, OccludedAreaReference,
    };
    use fs_render::photometry::source::{SourceSpectrum, SPECTRAL_BINS,
        BIN_WIDTH_NM, LUMINOUS_EFFICACY_LM_W};
    use fs_render::photometry::transport::PhotometricMedium;
    use fs_render::spectral::SpectralReflectance;
    use fs_render::tracer::{Material, Primitive, Shape};
    use fs_render::volumes::VolumeGrid;
    use fs_render::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget,
        CellPalette, CellPaletteBudget, MassAbsorptionSpectrum};
    use fs_render::volumes::transport::enclosure::GridEnclosureBudget;
    use fs_render::volumes::transport::scattering::{ScatteringGrid, ScatteringOptics,
        MassScatteringSpectrum, PhaseFunction};

    fn panel(x0: f64, x1: f64, y0: f64, y1: f64, z: f64, material: Material) -> Primitive {
        Primitive { shape: Shape::Mesh(TriMesh::new(
            vec![[x0,y0,z], [x1,y0,z], [x1,y1,z], [x0,y1,z]],
            vec![[0,2,1], [0,3,2]],
        )), material, emission: None }
    }

    fn calculate(cx: &Cx<'_>, shadow: bool, coarse: bool, blue_fraction: f64) -> Result<(), Box<dyn Error>> {
        let material = Material::SpectralDiffuse { reflectance: SpectralReflectance::constant(0.0)? };
        let mut primitives = vec![panel(-1.0,1.0,-1.0,1.0,2.0,material)];
        if shadow { primitives.push(panel(-0.25,0.25,-2.0,2.0,1.0,material)); }
        primitives.push(panel(-4.0,4.0,-4.0,4.0,0.0,material));
        let mut bins = [0.0; SPECTRAL_BINS];
        bins[4] = blue_fraction / observer::xyz_at_bin(4)[1];
        bins[44] = (1.0-blue_fraction) / observer::xyz_at_bin(44)[1];
        let sources = [PhotometricAreaEmitter::try_new(0,1000.0,1.0,
            SourceSpectrum::try_from_bins(bins)?,EmissionSidedness::Front)?];

        // Four actual cells: x-major left/right, then lower/upper z layers.
        // Absorption and scattering have independently indexed native palettes.
        let density = [1.0,4.0,2.0,8.0];
        let limits = CellPaletteBudget { max_models:2,max_cells:4 };
        let absorption = CellPalette::try_new(cx,vec![
            MassAbsorptionSpectrum::constant(0.125)?,
            MassAbsorptionSpectrum::table(&[(380.0,0.0),(780.0,1.0)])?,
        ],Arc::from([0,1,1,0]),limits)?;
        let absorption = AbsorbingGrid::try_new_indexed(cx,
            VolumeGrid::new([2,1,2],&density,[-2.0,-2.0,0.5],[2.0,4.0,0.5]),
            absorption,None,AbsorbingGridBudget {max_grid_cells:4,max_steps_per_segment:32})?;
        let scattering = CellPalette::try_new(cx,vec![
            ScatteringOptics {scattering:MassScatteringSpectrum::rayleigh(400.0,0.125)?,phase:PhaseFunction::Rayleigh},
            ScatteringOptics {scattering:MassScatteringSpectrum::constant(0.25)?,phase:PhaseFunction::Isotropic},
        ],Arc::from([1,0,1,0]),limits)?;
        let medium = ScatteringGrid::try_new_indexed(cx,absorption,scattering)?;
        let reference = OccludedAreaReference::new_with_medium(&primitives,&sources,
            PhotometricMedium::Scattering(&medium),AreaReferenceBudget::default(),
            GridEnclosureBudget::default(),cx)?;
        let probes = [IlluminanceProbe::try_new(Point3::new(0.0,0.0,0.0),Vec3::new(0.0,0.0,1.0))?];
        let budget = AreaOcclusionBudget {absolute_error_lux:2.0,
            max_patch_evaluations:if coarse {2} else {4096},..AreaOcclusionBudget::default()};
        let mut job = OccludedAreaJob::new(&reference,&probes,budget,cx)?;
        while !job.is_complete() { job.advance(cx,1)?; }

        // Independent 70-digit 2D surface quadrature, checked with two orders.
        // These diagnostic constants never participate in the production bound.
        // The spectral weights retain the admitted normalization, not an
        // assumed exact one-lumen round trip through binary64 arithmetic.
        let references = if shadow { [10.819_576_768_452_719,6.415_476_278_584_347] }
            else { [24.759_454_271_862_572,14.911_501_983_981_905] };
        let analytic: f64 = [4,44].into_iter().zip(references).map(|(bin,value)| value
            *sources[0].spectrum().bins_w_per_lumen_nm()[bin]*observer::xyz_at_bin(bin)[1]
            *BIN_WIDTH_NM*LUMINOUS_EFFICACY_LM_W).sum();
        eprintln!("authored spectral-grid direct-emission case; calibration and total physical error unknown");
        eprintln!("independent_surface_integral_f64_lux={analytic}; shadow={shadow}; blue_photopic_fraction={blue_fraction}");
        eprintln!("medium_identity={}; field density, absorption and out-scattering are retained; in-scattering is omitted",medium.content_identity());
        let stdout = io::stdout(); let mut writer = stdout.lock();
        OccludedAreaIlluminance::write_csv_header(&mut writer)?;
        for (index,report) in job.results().iter().enumerate() { report.write_csv_row(index,&mut writer)?; }
        writer.flush()?;
        Ok(())
    }

    pub fn run() -> Result<(), Box<dyn Error>> {
        let mut shadow = false; let mut coarse = false; let mut color = None;
        for arg in std::env::args().skip(1) {
            match arg.as_str() {
                "--shadow" if !shadow => shadow = true,
                "--coarse" if !coarse => coarse = true,
                "--blue" if color.is_none() => color = Some(1.0),
                "--red" if color.is_none() => color = Some(0.0),
                _ => return Err(io::Error::new(io::ErrorKind::InvalidInput,
                    "usage: photometric_area_grid [--shadow] [--coarse] [--blue | --red]").into()),
            }
        }
        let gate = CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(&gate,arena,StreamKey {seed:173,kernel_id:0x41524752,tile:0,iteration:0},
                Budget::INFINITE,ExecMode::Deterministic);
            calculate(&cx,shadow,coarse,color.unwrap_or(0.75))
        })
    }
}

#[cfg(feature = "photometry")]
fn main() -> Result<(), Box<dyn std::error::Error>> { demo::run() }
#[cfg(not(feature = "photometry"))]
fn main() {
    eprintln!("photometric_area_grid requires --features photometry");
    std::process::exit(2);
}
