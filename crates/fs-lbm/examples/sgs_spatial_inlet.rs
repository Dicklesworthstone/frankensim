//! Advected, solenoidal Fourier inlet through the native D3Q19 SGS boundary.
//! Two explicit modes exercise spatial/temporal inlet coupling; they do not
//! represent a measured turbulent spectrum or a developed-flow validation.
//! Transverse periods equal the domain periods, avoiding an artificial seam.
//! All quantities are lattice units (one time step = one lattice time unit).

use fs_exec::CancelGate;
use fs_lbm::d3q19::{BoundaryGrid3, BoundarySpec3, Face3, FaceBoundary3, InletFlux3};
use fs_lbm::sgs3::{SgsConfig, SgsModel, SgsUse};
use fs_math::det;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let steps = std::env::args()
        .nth(1)
        .map_or(Ok(200_usize), |s| s.parse())?;
    if steps > 100_000 {
        return Err("steps must be <=100000".into());
    }
    let periodic = FaceBoundary3::Periodic;
    let mut grid = BoundaryGrid3::new(
        32,
        16,
        16,
        0.56,
        [0.0; 3],
        BoundarySpec3::new([
            FaceBoundary3::Velocity {
                velocity: [0.04, 0.0, 0.0],
            },
            FaceBoundary3::Pressure { density: 1.0 },
            periodic,
            periodic,
            periodic,
            periodic,
        ]),
    );
    let config = SgsConfig {
        model: SgsModel::Wale { cw: 0.5 },
        filter_width_lu: 1.0,
        maximum_eddy_to_molecular: 50.0,
        usage: SgsUse::Exploration,
    };
    let mut workspace = grid.sgs_workspace()?;
    let gate = CancelGate::new_clock_free();
    println!(
        "step,prescribed_mean_u,normal_correction,max_eddy_viscosity,molecular_shear_dissipation,sgs_shear_dissipation,total_mass"
    );
    for step in 1..=steps {
        // Both polarizations satisfy k dot e = 0. The same frozen field
        // translates at the mean speed, instead of redrawing node noise.
        let inlet =
            grid.set_velocity_profile(Face3::XMin, InletFlux3::MeanNormal(-0.04), &gate, |p| {
                let x = -0.04 * step as f64;
                let y = p[1] as f64;
                let z = p[2] as f64;
                let k = core::f64::consts::TAU / 16.0;
                let a = 0.004 * det::sin(k * (x + y + z));
                let b = 0.003 * det::sin(k * (2.0 * x - y + z));
                Ok([0.04 + a + b, -a, -2.0 * b])
            })?;
        let report = grid.step_with_sgs(&gate, config, &mut workspace, None)?.sgs;
        println!(
            "{step},{},{},{},{},{},{}",
            inlet.prescribed_mean_velocity_lu[0],
            inlet.normal_velocity_correction_lu,
            report.maximum_eddy_viscosity_lu,
            report.resolved_dissipation_lu,
            report.sgs_dissipation_lu,
            grid.total_mass()
        );
    }
    Ok(())
}
