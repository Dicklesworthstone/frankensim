//! Headless bounds on point-fixture light reflected by an original spectral mesh.
//!
//! cargo run -p fs-render --features photometry --example photometric_surface_reference
//! Options: --grid, --shadow, and one of --blue / --red.
//! Authored model inputs, not measured fixture/material/medium calibration.
//! Stdout contains component-scoped cd/m² CSV; stderr contains independent
//! formula diagnostics. This is not a full image, first-hit or pixel-area bound.

#[cfg(feature = "photometry")]
mod example {
    use std::{error::Error, io::{self, Write}, sync::Arc};
    use core::f64::consts::PI;
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    use fs_geom::Point3;
    use fs_render::charts::TriMesh;
    use fs_render::instances::RigidTransform;
    use fs_render::photometry::{PhotometricScene, TransportOptions};
    use fs_render::photometry::medium::{HomogeneousMedium, PhaseFunction};
    use fs_render::photometry::transport::PhotometricMedium;
    use fs_render::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum,
        SPECTRAL_BINS, LUMINOUS_EFFICACY_LM_W, photopic_weight_nm, wavelength_nm};
    use fs_render::photometry::source::validation::surface::{DiffuseSurfaceReference,
        DiffuseSurfaceSample, DiffuseSurfaceJob, EnclosedDiffuseLuminance};
    use fs_render::spectral::SpectralReflectance;
    use fs_render::tracer::{Material, Primitive, Shape};
    use fs_render::volumes::VolumeGrid;
    use fs_render::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget,
        MassAbsorptionSpectrum, CellPalette, CellPaletteBudget};
    use fs_render::volumes::transport::scattering::{ScatteringGrid, MassScatteringSpectrum, ScatteringOptics};

    fn screen(x0: f64, x1: f64) -> Result<Primitive, Box<dyn Error>> {
        Ok(Primitive { shape: Shape::Mesh(TriMesh::new(vec![
            [x0,-1.0,1.0], [x1,-1.0,1.0], [x1,2.0,1.0], [x0,2.0,1.0]],
            vec![[0,1,2],[0,2,3]])),
            material: Material::SpectralDiffuse { reflectance: SpectralReflectance::constant(0.0)? },
            emission: None })
    }

    fn calculate(cx: &Cx<'_>, grid: bool, shadow: bool, blue: f64) -> Result<(), Box<dyn Error>> {
        let mut primitives = vec![Primitive {
            shape: Shape::Mesh(TriMesh::new(
                vec![[-2.0,-2.0,0.0],[2.0,-2.0,0.0],[0.0,2.0,0.0]],vec![[0,1,2]])),
            material: Material::SpectralDiffuse {
                reflectance: SpectralReflectance::table(&[(380.0,0.125),(780.0,0.875)])?,
            }, emission: None,
        }];
        if shadow {
            // At x=1/8 the first strip blocks only the observer leg. At x=1/4
            // the second strip blocks only the incident leg. Neither is a
            // rounded shadow estimate used to construct the production bound.
            primitives.push(screen(0.0625,0.09375)?);
            primitives.push(screen(0.1875,0.21875)?);
        }
        let density = [1.0,2.0,4.0];
        let cap = CellPaletteBudget { max_models: 2, max_cells: 3 };
        let a = CellPalette::try_new(cx, vec![MassAbsorptionSpectrum::constant(0.125)?,
            MassAbsorptionSpectrum::table(&[(380.0,0.0),(780.0,0.25)])?], Arc::from([0,1,0]), cap)?;
        let absorption = AbsorbingGrid::try_new_indexed(cx,
            VolumeGrid::new([1,1,3], &density, [-0.5,-0.5,1.0], [1.0;3]), a, None,
            AbsorbingGridBudget { max_grid_cells: 3, max_steps_per_segment: 64 })?;
        let s = CellPalette::try_new(cx, vec![
            ScatteringOptics { scattering: MassScatteringSpectrum::rayleigh(400.0,0.25)?, phase: PhaseFunction::Rayleigh },
            ScatteringOptics { scattering: MassScatteringSpectrum::constant(0.125)?, phase: PhaseFunction::Isotropic },
        ], Arc::from([1,0,1]), cap)?;
        let physical_grid = ScatteringGrid::try_new_indexed(cx, absorption, s)?;
        let medium = if grid { PhotometricMedium::Scattering(&physical_grid) }
            else { PhotometricMedium::Homogeneous(HomogeneousMedium::vacuum()) };
        // Fractions are nominal photopic fractions. SourceSpectrum retains its
        // actual binary64 normalization, and the new reference encloses that sum.
        let mut bins = [0.0; SPECTRAL_BINS];
        bins[4] = blue / (LUMINOUS_EFFICACY_LM_W * photopic_weight_nm(4));
        bins[44] = (1.0-blue) / (LUMINOUS_EFFICACY_LM_W * photopic_weight_nm(44));
        let sources = [PhotometricSource::try_new(
            RigidTransform::try_new([0.0,0.0,0.0,1.0],[0.0,0.0,5.0])?,
            4000.0*PI, 1.0, AngularDistribution::isotropic(), SourceSpectrum::try_from_bins(bins)?)?];
        let scene = PhotometricScene::try_with_medium(&primitives, &sources, medium, TransportOptions::default(), cx)?;
        let reference = DiffuseSurfaceReference::new(&scene, 0, cx)?;
        let positions = [-0.25,0.0,0.125,0.25];
        let samples: Vec<_> = positions.into_iter().map(|x| DiffuseSurfaceSample::try_new(
            0, [0.25-x/4.0,0.25+x/4.0,0.5], Point3::new(0.0,0.0,2.5),
        )).collect::<Result<_,_>>()?;
        let mut job = DiffuseSurfaceJob::new(&reference, &samples, cx)?;
        while !job.is_complete() { job.advance(cx, 2)?; }
        eprintln!("authored direct diffuse surface reference; no first-hit/pixel-area or total physical-accuracy claim");
        eprintln!("spectral_grid={grid} native_shadows={shadow}; omitted interreflection/in-scattering/emission and calibration remain unknown");
        let stdout = io::stdout(); let mut out = stdout.lock();
        EnclosedDiffuseLuminance::write_csv_header(&mut out)?;
        for (index, (&x, result)) in positions.iter().zip(job.results()).enumerate() {
            let light_distance = (25.0_f64+x*x).sqrt();
            let observer_distance = (6.25_f64+x*x).sqrt();
            let mut physical_reference = 0.0;
            for (bin, fraction) in [(4,blue),(44,1.0-blue)] {
                let lambda = wavelength_nm(bin);
                let rho = 0.125+0.75*(lambda-380.0)/400.0;
                let middle = (lambda-380.0)/1600.0+0.25*(400.0/lambda).powi(4);
                let tau = if grid { (1.25+2.0*middle)*light_distance/5.0
                    +(0.25+middle)*observer_distance/2.5 } else { 0.0 };
                physical_reference += fraction*rho*(-tau).exp();
            }
            physical_reference *= 5000.0/(PI*light_distance.powi(3));
            // Independent projection of the two strips onto the floor samples.
            let hit = |coordinate: f64| (0.0625..=0.09375).contains(&coordinate)
                || (0.1875..=0.21875).contains(&coordinate);
            let light_blocked = shadow && hit(x*4.0/5.0);
            let observer_blocked = shadow && hit(x*3.0/5.0);
            if light_blocked || observer_blocked { physical_reference = 0.0; }
            eprintln!("sample={index} x_m={x} light_blocked={light_blocked} observer_blocked={observer_blocked} analytic_f64_cd_m2={physical_reference}");
            result.write_csv_row(index,&mut out)?;
        }
        out.flush()?;
        Ok(())
    }

    pub fn run() -> Result<(), Box<dyn Error>> {
        let (mut grid,mut shadow,mut color) = (false,false,None);
        for arg in std::env::args().skip(1) {
            match arg.as_str() {
                "--grid" => grid=true,
                "--shadow" => shadow=true,
                "--blue" if color.is_none() => color=Some(1.0),
                "--red" if color.is_none() => color=Some(0.0),
                _ => return Err(io::Error::new(io::ErrorKind::InvalidInput,
                    "usage: photometric_surface_reference [--grid] [--shadow] [--blue | --red]").into()),
            }
        }
        let gate=CancelGate::new_clock_free();
        let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| calculate(&Cx::new(&gate,arena,
            StreamKey{seed:173,kernel_id:0x5355_5246,tile:0,iteration:0},
            Budget::INFINITE,ExecMode::Deterministic),grid,shadow,color.unwrap_or(0.5)))
    }
}
#[cfg(feature = "photometry")]
fn main() -> Result<(), Box<dyn std::error::Error>> { example::run() }
#[cfg(not(feature = "photometry"))]
fn main() { eprintln!("photometric_surface_reference requires --features photometry"); std::process::exit(2); }
