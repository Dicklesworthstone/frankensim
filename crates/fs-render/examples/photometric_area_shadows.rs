//! Headless finite-panel shadows with an explicit numerical lux budget.
//!
//! cargo run -p fs-render --features photometry --example photometric_area_shadows
//! Add `-- --coarse` to stop at the two initial source triangles. The coarse
//! output still encloses the result but explicitly reports unmet accuracy.
//! Add `-- --haze` for native homogeneous absorption/out-scattering, or combine
//! `-- --haze --coarse` to see unresolved shadow AND attenuation bounds.
//!
//! A 2 m by 2 m downward panel at z=2 m emits 1000 lm. An opaque strip at z=1 m
//! blocks part of it; the receiving floor lies at z=0. The emitter, blockers,
//! material, source spectrum and probes are the original native owners.
//! Stdout is scoped direct-area bound CSV. Stderr prints independent f64
//! references, not measured data or a substitute for the production enclosure.
//! Haze uses sigma_a=sigma_s=0.125 m^-1. In-scattered light, reflections,
//! fixture calibration and total physical error remain unknown.

#[cfg(feature = "photometry")]
mod demo {
    use std::error::Error;
    use std::f64::consts::PI;
    use std::io::{self, Write};

    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    use fs_geom::{Point3, Vec3};
    use fs_render::charts::TriMesh;
    use fs_render::photometry::IlluminanceProbe;
    use fs_render::photometry::medium::{HomogeneousMedium, PhaseFunction};
    use fs_render::photometry::area::{EmissionSidedness, PhotometricAreaEmitter};
    use fs_render::photometry::area::reference::AreaReferenceBudget;
    use fs_render::photometry::area::reference::occlusion::{
        AreaOcclusionBudget, OccludedAreaIlluminance, OccludedAreaJob, OccludedAreaReference,
    };
    use fs_render::photometry::source::SourceSpectrum;
    use fs_render::spectral::SpectralReflectance;
    use fs_render::tracer::{Material, Primitive, Shape};

    fn panel(x0: f64, x1: f64, y0: f64, y1: f64, z: f64, material: Material) -> Primitive {
        Primitive {
            shape: Shape::Mesh(TriMesh::new(
                vec![[x0, y0, z], [x1, y0, z], [x1, y1, z], [x0, y1, z]],
                vec![[0, 2, 1], [0, 3, 2]],
            )),
            material,
            emission: None,
        }
    }

    // Independent double antiderivative of h^2 / (x^2 + y^2 + h^2)^2.
    fn rectangle(x0: f64, x1: f64, y0: f64, y1: f64, h: f64) -> f64 {
        let f = |x: f64, y: f64| {
            let ax = (x * x + h * h).sqrt();
            let by = (y * y + h * h).sqrt();
            0.5 * (x / ax * (y / ax).atan() + y / by * (x / by).atan())
        };
        f(x1, y1) - f(x0, y1) - f(x1, y0) + f(x0, y0)
    }

    fn calculate(cx: &Cx<'_>, coarse: bool, haze: bool) -> Result<(), Box<dyn Error>> {
        let material = Material::SpectralDiffuse {
            reflectance: SpectralReflectance::constant(0.0)?,
        };
        let primitives = [
            panel(-1.0, 1.0, -1.0, 1.0, 2.0, material),
            panel(-0.25, 0.25, -2.0, 2.0, 1.0, material),
            panel(-4.0, 4.0, -4.0, 4.0, 0.0, material),
        ];
        let sources = [PhotometricAreaEmitter::try_new(
            0, 1000.0, 1.0, SourceSpectrum::equal_energy(), EmissionSidedness::Front,
        )?];
        let xs = [-1.0, -0.5, 0.0, 0.5, 1.0];
        let probes: Vec<_> = xs.iter().map(|&x| {
            IlluminanceProbe::try_new(Point3::new(x, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0))
        }).collect::<Result<_, _>>()?;
        let medium = if haze {
            HomogeneousMedium::try_new(0.125, 0.125, PhaseFunction::HenyeyGreenstein { g: 0.7 })?
        } else {
            HomogeneousMedium::vacuum()
        };
        let reference = OccludedAreaReference::new(
            &primitives, &sources, AreaReferenceBudget::default(), cx,
        )?.with_homogeneous_medium(&medium);
        let budget = AreaOcclusionBudget {
            absolute_error_lux: 1.0,
            max_patch_evaluations: if coarse { 2 } else { 16_384 },
            ..AreaOcclusionBudget::default()
        };
        let mut job = OccludedAreaJob::new(&reference, &probes, budget, cx)?;
        while !job.is_complete() {
            job.advance(cx, 2)?;
        }
        eprintln!("authored opaque-mesh direct-emission case; in-scattered/reflected light and total physical error unknown");
        eprintln!("absorption_per_m={} scattering_per_m={}; attenuation covers the full source-to-receiver path",
            medium.absorption_per_m(), medium.scattering_per_m());
        eprintln!("requested_midpoint_error_lux={} max_patch_evaluations={}; inspect actual termination per probe",
            budget.absolute_error_lux, budget.max_patch_evaluations);
        eprintln!("analytic shadow projection at source plane: x_source in [-0.5-x_probe,0.5-x_probe]; y is fully covered");
        let stdout = io::stdout();
        let mut writer = stdout.lock();
        OccludedAreaIlluminance::write_csv_header(&mut writer)?;
        for (index, (report, x)) in job.results().iter().zip(xs).enumerate() {
            let all = rectangle(-1.0 - x, 1.0 - x, -1.0, 1.0, 2.0);
            let left = (-0.5 - x).max(-1.0);
            let right = (0.5 - x).min(1.0);
            let blocked = if left < right {
                rectangle(left - x, right - x, -1.0, 1.0, 2.0)
            } else {
                0.0
            };
            // Nominal 1000/(4*pi) cd/m2; the production enclosure separately
            // includes the retained source spectrum's actual observer integral.
            // Haze references come from independent 70-digit 2D surface
            // integration of L*h²/r^4 * exp(-r/4), split at the projected strip.
            // Orders 32 and 48 agreed within 8e-38 lx; they are rounded here.
            // This is not a centre-ray attenuation factor on the vacuum value.
            const HAZE_REFERENCE_LUX: [f64; 5] = [
                22.595_096_099_976_296, 19.623_377_792_975_71,
                15.435_136_245_708_192, 19.623_377_792_975_71,
                22.595_096_099_976_296,
            ];
            let independent = if haze { HAZE_REFERENCE_LUX[index] }
                else { 1000.0 / (4.0 * PI) * (all - blocked) };
            eprintln!("probe={index} x_m={x} independent_reference_f64_lux={independent} termination={:?}", report.termination);
            report.write_csv_row(index, &mut writer)?;
        }
        writer.flush()?;
        Ok(())
    }

    pub fn run() -> Result<(), Box<dyn Error>> {
        let args: Vec<_> = std::env::args().skip(1).collect();
        let (mut coarse, mut haze) = (false, false);
        for flag in &args {
            match flag.as_str() {
                "--coarse" if !coarse => coarse = true,
                "--haze" if !haze => haze = true,
                _ => return Err(io::Error::new(io::ErrorKind::InvalidInput,
                    "usage: photometric_area_shadows [--haze] [--coarse]").into()),
            }
        }
        let gate = CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate, arena,
                StreamKey { seed: 173, kernel_id: 0x4152_5348, tile: 0, iteration: 0 },
                Budget::INFINITE, ExecMode::Deterministic,
            );
            calculate(&cx, coarse, haze)
        })
    }
}

#[cfg(feature = "photometry")]
fn main() -> Result<(), Box<dyn std::error::Error>> { demo::run() }

#[cfg(not(feature = "photometry"))]
fn main() {
    eprintln!("photometric_area_shadows requires --features photometry");
}
