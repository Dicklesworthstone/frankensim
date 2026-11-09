//! Signed vibration estimates for a restrained bar with a held free strain.
//! Run `cargo run -p fs-solid --example cutfem_prestress_modes`.
//!
//! This synthetic fixture uses E=30 MPa, nu=0, rho=1000 kg/m^3 and a
//! 1 x 0.05 x 0.037 m bar cut from a 1 x 0.05 x 0.05 m background grid.
//! Both end faces are fixed. Uniform axial manufacturing extension is supplied
//! explicitly as free strain; no temperature-to-expansion law is invented.
//! With nu=0 the resulting axial compression has no lateral traction, so zero
//! displacement is the exact continuum equilibrium. The native free-row force
//! balance is checked independently before preparing geometric stiffness.
//!
//! Material, geometry and discretization are illustrative assumptions, not a
//! measured source card or a converged beam benchmark. The reported signed
//! squared-frequency estimates use the original consistent physical inertia.
//! Negative values indicate instability of this frozen reference-geometry
//! approximation; there is no completeness, critical-load or nonlinear claim.

use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::{CutElasticity3, ElasticMaterialField3};
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_exec::CancelGate;
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solid::cutfront::dynamics3::CutElasticityRef3;
use fs_solid::cutfront::prestress3::{
    CutPrestressLimits3, CutPrestressModeOptions3, CutPrestressModes3,
};
use fs_solver::LinearOp;

const LENGTH: f64 = 1.0;
const WIDTH: f64 = 0.05;
const HEIGHT: f64 = 0.037;
const YOUNG: f64 = 30e6;

struct Top;
impl CutSdf3 for Top {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - HEIGHT
    }

    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(HEIGHT, HEIGHT)
    }

    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let value = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(value, value)
    }
}

fn proceed() -> ControlFlow<()> {
    ControlFlow::Continue(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let law = IsotropicElastic::new(YOUNG, 0.0, 0.02)?;
    // Explicit synthetic identity: this is not a material-database citation.
    let material = ElasticAcousticState::from_isotropic(&law, 1000.0, ContentHash([73; 32]))?;
    let mut callback = |_| proceed();
    let mut quadrature = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 1,
            ..Default::default()
        },
        &mut callback,
    )?;
    let operator = CutElasticity3::build_with_material_field(
        HexCell::try_new([0.0; 3], [LENGTH, WIDTH, WIDTH])?,
        [12, 2, 2],
        &Top,
        ElasticMaterialField3::Uniform(&material),
        &|p| p[0] == 0.0 || p[0] == LENGTH,
        Default::default(),
        &mut quadrature,
    )?;
    let mass = operator.prepare_material_mass(100_000, proceed)?;
    let occupancy = vec![1.0; operator.cells()];
    let displacement = vec![0.0; operator.n()];
    let gate = CancelGate::new_clock_free();
    let limits = CutPrestressLimits3 {
        max_native_dofs: 2_000,
        max_cells: 100,
        max_mass_entries: 600_000,
        max_geometric_entries: 100_000,
        max_workspace_components: 200_000,
        max_operator_applications: 100_000,
        ..Default::default()
    };
    let query = CutPrestressModeOptions3 {
        modes: 3,
        max_ticks: 2_048,
        // Absolute tolerances are expressed in this SI pencil's s^-2 units.
        spectral_tolerance: 1e-4,
        residual_relative_tolerance: 1e-7,
        residual_absolute_tolerance: 1e-4,
        mass_orthogonality_tolerance: 1e-8,
    };

    println!("# Synthetic restrained CutFEM bar; fixed reference geometry");
    println!("# Signed Ritz estimates; spectral completeness is not certified");
    println!(
        "free_strain,axial_stress_pa,free_force_residual_n,mode,lambda_s_minus_2,whitened_residual_s_minus_2"
    );
    for extension in [0.0, 0.003, 0.007] {
        let free = operator.prepare_eigenstrain(
            &|_, _| Ok([extension, 0.0, 0.0, 0.0, 0.0, 0.0]),
            100_000,
            proceed,
        )?;
        let force = operator.eigenstrain_load(&free, proceed)?;
        let residual = force.iter().fold(0.0_f64, |value, f| value.max(f.abs()));
        let force_scale = YOUNG * extension.abs() * WIDTH * HEIGHT;
        if !residual.is_finite() || residual > 1e-10 * force_scale.max(1.0) {
            return Err("zero state failed the native free-row equilibrium check".into());
        }
        let prestress = operator.prepare_eigenstrain_prestress(
            &free,
            &displacement,
            0.05,
            100_000,
            100_000,
            proceed,
        )?;
        let modes = CutPrestressModes3::new(
            CutElasticityRef3::Cartesian(&operator),
            &prestress,
            &mass,
            &occupancy,
            limits,
            &gate,
        )?
        .lowest_modes(query)?;
        for (i, mode) in modes.modes.iter().enumerate() {
            println!(
                "{extension:.6},{:.9e},{residual:.9e},{},{:.9e},{:.9e}",
                -YOUNG * extension,
                i + 1,
                mode.eigenvalue,
                mode.whitened_residual_norm,
            );
        }
    }
    Ok(())
}
