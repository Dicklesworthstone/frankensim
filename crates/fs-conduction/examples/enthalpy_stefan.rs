//! A computed melting front in a tetrahedral slab, checked against the
//! one-phase Stefan similarity solution. All material numbers are explicitly
//! declared numerical-reference inputs, not experimentally validated data.
//!
//! `cargo run -p fs-conduction --example enthalpy_stefan -- 40 120`
//!
//! The optional arguments are axial cells and time steps. Geometry is fixed;
//! this example tests heat transport and phase state, not liquid motion.
//! With eta=x/(2 sqrt(alpha t)) and s=2 lambda sqrt(alpha t), the liquid
//! temperature is Tm+A(erf(lambda)-erf(eta)); the solid remains at Tm.
//! A=(L/cp) lambda sqrt(pi) exp(lambda^2) satisfies the Stefan jump balance.
//! Only this initial field and q_in(t)=k A/sqrt(pi alpha t) are inputs to the
//! production solve. The front position is never supplied after initialization.

use fs_blake3::ContentHash;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyStepConfig,
};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, ScalarField, ThermalBc,
    ThermalBoundaryBuilder,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_solver::NewtonKrylovConfig;
use std::error::Error;

const LENGTH: f64 = 1.0;
const WIDTH: f64 = 0.02;
const RHO: f64 = 1.0;
const CP: f64 = 1.0;
const CONDUCTIVITY: f64 = 1.0;
const LATENT: f64 = 1.0;
const MELTING: f64 = 300.0;
const LAMBDA: f64 = 0.5;
const START: f64 = 0.04;
const END: f64 = 0.16;

fn alpha() -> f64 {
    CONDUCTIVITY / (RHO * CP)
}
fn amplitude() -> f64 {
    (LATENT / CP) * LAMBDA * std::f64::consts::PI.sqrt() * (LAMBDA * LAMBDA).exp()
}
fn front(time: f64) -> f64 {
    2.0 * LAMBDA * (alpha() * time).sqrt()
}
fn exact_temperature(x: f64, time: f64) -> f64 {
    if x >= front(time) {
        return MELTING;
    }
    MELTING
        + amplitude()
            * (fs_math::det::erf(LAMBDA) - fs_math::det::erf(x / (2.0 * (alpha() * time).sqrt())))
}

fn chart() -> Result<EquilibriumEnthalpyPhaseCurve, Box<dyn Error>> {
    let knots = [
        (-5.0, MELTING - 5.0, 0.0),
        (0.0, MELTING, 0.0),
        (LATENT, MELTING, 1.0),
        (LATENT + 5.0, MELTING + 5.0, 1.0),
    ]
    .into_iter()
    .map(
        |(h, temperature_k, liquid_mass_fraction)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k,
            liquid_mass_fraction,
            bulk_density_kg_m3: RHO,
        },
    )
    .collect();
    // Explicit identity for this synthetic reference chart; no sourced claim.
    Ok(EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash::from_slice(&[19; 32]).ok_or("invalid reference identity")?,
        knots,
    )?)
}

/// Quantities measured from the computed spatial field and independent solution.
#[derive(Debug)]
pub struct StefanResult {
    /// Number of actual thermal degrees of freedom.
    pub vertices: usize,
    /// Number of tetrahedral finite elements.
    pub elements: usize,
    /// Location of the cross-section-averaged 50% liquid-fraction contour [m].
    pub front_m: f64,
    /// Similarity front at the final time [m].
    pub reference_front_m: f64,
    /// Liquid mass divided by reference density and slab cross section [m].
    pub equivalent_molten_length_m: f64,
    /// Reference-mass-weighted nodal RMS temperature discrepancy [K].
    pub temperature_rms_error_k: f64,
    /// Change of total nodal enthalpy over the complete run [J].
    pub stored_energy_change_j: f64,
    /// Independent exact time integral of the supplied boundary heat [J].
    pub applied_heat_j: f64,
    /// Stored change minus independently integrated heat [J].
    pub energy_residual_j: f64,
    /// Largest accepted per-step energy mismatch [J].
    pub max_step_energy_residual_j: f64,
    /// Actual total inner iterations over accepted Newton steps.
    pub krylov_iterations: usize,
}

/// Run one bounded refinement case. No front or temperature is clamped during
/// evolution. Insulated side and far-end boundaries close the energy account.
pub fn run(cx: &Cx<'_>, cells: usize, steps: usize) -> Result<StefanResult, Box<dyn Error>> {
    if !(4..=256).contains(&cells) || !(1..=4096).contains(&steps) {
        return Err("expected 4..=256 cells and 1..=4096 steps".into());
    }
    let (complex, positions) = box_grid([cells, 1, 1], [LENGTH, WIDTH, WIDTH]);
    let mesh = ConductionMesh::new(complex, positions)?;
    let curve = chart()?;
    let material = ConductivityModel::isotropic_declared(CONDUCTIVITY)?;
    let source = ScalarField::uniform("no volumetric heating", 0.0)?;
    let integrator = EnthalpyBackwardEuler::uniform(
        cx,
        &mesh,
        &curve,
        RHO,
        EnthalpyBudget {
            max_vertices: mesh.vertex_count(),
            max_elements: mesh.element_count(),
        },
    )?;
    let initial: Vec<f64> = mesh
        .positions()
        .iter()
        .map(|position| {
            let x = position[0];
            if (x - front(START)).abs() <= 1e-14 {
                // Midpoint value of the initial enthalpy jump. Its point value does
                // not affect the continuum solution; the FEM field resolves it.
                0.5 * LATENT
            } else if x < front(START) {
                LATENT + CP * (exact_temperature(x, START) - MELTING)
            } else {
                0.0
            }
        })
        .collect();
    let mut enthalpy = initial.clone();
    let mut temperature = Vec::new();
    let mut liquid = Vec::new();
    let mut krylov_iterations = 0;
    let mut max_step_energy_residual_j = 0.0_f64;
    let flux_coefficient = CONDUCTIVITY * amplitude() / (std::f64::consts::PI * alpha()).sqrt();
    let config = EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-13,
            relative_tolerance: 1e-11,
            linear_restart: 24,
            max_linear_cycles: 16,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-3,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 32,
        energy_tolerance_j: 2e-11,
    };
    for step in 0..steps {
        let t0 = (END - START).mul_add(step as f64 / steps as f64, START);
        let t1 = (END - START).mul_add((step + 1) as f64 / steps as f64, START);
        // The exact mean of q0/sqrt(t) over the step removes boundary-time
        // quadrature error from the independent whole-run energy comparison.
        let mean_inward_flux = 2.0 * flux_coefficient / (t1.sqrt() + t0.sqrt());
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "heated end",
                |face| on_box_face(face.centroid[0], 0.0),
                ThermalBc::neumann(-mean_inward_flux)?,
            )?
            .adiabatic_remainder()
            .finish()?;
        let next = integrator.advance(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                element_materials: None,
                source: &source,
            },
            None,
            &enthalpy,
            t1 - t0,
            config,
        )?;
        krylov_iterations += next
            .newton
            .history
            .iter()
            .map(|row| row.linear_iterations)
            .sum::<usize>();
        max_step_energy_residual_j = max_step_energy_residual_j.max(next.energy_residual_j.abs());
        enthalpy = next.specific_enthalpy_j_kg;
        temperature = next.temperature;
        liquid = next.liquid_mass_fraction;
    }
    let masses = integrator.reference_nodal_masses_kg();
    let stored_energy_change_j: f64 = masses
        .iter()
        .zip(&enthalpy)
        .zip(&initial)
        .map(|((mass, h), old)| mass * (h - old))
        .sum();
    let applied_heat_j = WIDTH * WIDTH * 2.0 * flux_coefficient * (END.sqrt() - START.sqrt());
    let equivalent_molten_length_m = masses
        .iter()
        .zip(&liquid)
        .map(|(mass, fraction)| mass * fraction)
        .sum::<f64>()
        / (RHO * WIDTH * WIDTH);
    let total_mass: f64 = masses.iter().sum();
    let temperature_rms_error_k = (masses
        .iter()
        .zip(&temperature)
        .zip(mesh.positions())
        .map(|((mass, t), x)| mass * (t - exact_temperature(x[0], END)).powi(2))
        .sum::<f64>()
        / total_mass)
        .sqrt();
    let averages: Vec<f64> = liquid
        .as_chunks::<4>()
        .0
        .iter()
        .map(|plane| plane.iter().sum::<f64>() / 4.0)
        .collect();
    let crossing = averages
        .windows(2)
        .position(|pair| pair[0] >= 0.5 && pair[1] < 0.5)
        .ok_or("computed liquid field does not contain a 50% front")?;
    let fraction = (averages[crossing] - 0.5) / (averages[crossing] - averages[crossing + 1]);
    let front_m = (crossing as f64 + fraction) * LENGTH / cells as f64;
    Ok(StefanResult {
        vertices: mesh.vertex_count(),
        elements: mesh.element_count(),
        front_m,
        reference_front_m: front(END),
        equivalent_molten_length_m,
        temperature_rms_error_k,
        stored_energy_change_j,
        applied_heat_j,
        energy_residual_j: stored_energy_change_j - applied_heat_j,
        max_step_energy_residual_j,
        krylov_iterations,
    })
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut arguments = std::env::args().skip(1);
    let cells = arguments
        .next()
        .map(|value| value.parse())
        .transpose()?
        .unwrap_or(40);
    let steps = arguments
        .next()
        .map(|value| value.parse())
        .transpose()?
        .unwrap_or(120);
    if arguments.next().is_some() {
        return Err("usage: enthalpy_stefan [cells] [steps]".into());
    }
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 19, kernel_id: 19, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let result = run(&cx, cells, steps)?;
        println!("Stefan melting: cells={cells}, steps={steps}, vertices={}, tetrahedra={}",
            result.vertices, result.elements);
        println!("computed_front_m={:.9}; reference_front_m={:.9}; molten_length_m={:.9}",
            result.front_m, result.reference_front_m, result.equivalent_molten_length_m);
        println!("temperature_rms_error_K={:.9e}; stored_change_J={:.9e}; applied_heat_J={:.9e}",
            result.temperature_rms_error_k, result.stored_energy_change_j, result.applied_heat_j);
        println!("whole_run_energy_residual_J={:.9e}; max_step_energy_residual_J={:.9e}; krylov_iterations={}",
            result.energy_residual_j, result.max_step_energy_residual_j, result.krylov_iterations);
        Ok(())
    })
}
