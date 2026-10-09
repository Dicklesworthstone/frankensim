//! Recover two heater-pulse amplitudes from spatial temperature histories.
//!
//! `cargo run -p fs-ascent --example enthalpy_calibration`
//! `cargo run -p fs-ascent --example enthalpy_calibration -- --ambient-radiation`
//!
//! The same fixed tetrahedral slab resolves sensible and latent heat in every
//! trial. The source is spatially localized; the phase front is computed.
//! A discrete enthalpy adjoint carries observation seeds through all accepted
//! steps. Existing SQP owns parameter bounds, line search and the KKT stop.
//! Synthetic observations demonstrate inverse computation, not experimental
//! validation or general identifiability. Chart, mass, mesh and time grid are
//! fixed. This small example retains at most 24 endpoint linearizations.
//! The explicit radiation mode includes the same ambient gray surface law in
//! observation generation, every candidate trajectory, and the physical adjoint.

use fs_ascent::sqp::SqpSample;
use fs_ascent::{SqpRunReport, SqpState};
use fs_blake3::ContentHash;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyStepConfig,
};
use fs_conduction::{
    AmbientRadiationConfig, AmbientRadiationPatch, ConductionMesh, ConductionProblem,
    ConductivityModel, EMISSIVITY_DIMS, LinearConfig, SURFACE_EMISSIVITY_PROPERTY, ScalarField,
    SurfaceEmissivity, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_solver::NewtonKrylovConfig;
use std::error::Error;

/// Declared synthetic heater densities [W/m3], before the spatial profile.
pub const TRUTH: [f64; 2] = [6.3, 4.7];
/// Admitted starting guess in the same physical coordinates.
pub const START: [f64; 2] = [4.5, 6.5];
const BOUNDS: [[f64; 2]; 2] = [[2.0, 9.0], [2.0, 9.0]];
const CELLS: usize = 6;
const STEPS: usize = 24;
const DT: f64 = 0.05;
const INITIAL_H: f64 = -0.17;
const SURFACE: &str = "radiating-heater-face";
const EMISSIVITY: f64 = 0.04;
const AMBIENT_K: f64 = 299.0;
const CONVECTION_W_M2_K: f64 = 0.01;

/// One cross-section-mean temperature at an accepted endpoint.
#[derive(Debug, Clone)]
pub struct Observation {
    pub endpoint: usize,
    pub plane: usize,
    pub temperature_k: f64,
}

/// Fixed physical experiment; observations may be supplied independently.
pub struct EnthalpyCalibration {
    mesh: ConductionMesh,
    curve: EquilibriumEnthalpyPhaseCurve,
    material: ConductivityModel,
    boundary: ThermalBoundary,
    source_profile: Vec<f64>,
    radiation: Option<AmbientRadiationPatch>,
    pub observations: Vec<Observation>,
}

/// One complete physical trajectory evaluation.
#[derive(Debug)]
pub struct Evaluation {
    pub loss_k2: f64,
    pub gradient: [f64; 2],
    pub final_enthalpy_j_kg: Vec<f64>,
    pub final_temperature_k: Vec<f64>,
    pub final_liquid_fraction: Vec<f64>,
    pub max_step_energy_residual_j: f64,
}

fn config() -> EnthalpyStepConfig {
    EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            linear_restart: 24,
            max_linear_cycles: 16,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-3,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 32,
        energy_tolerance_j: 1e-9,
    }
}

fn sensor(temperature: &[f64], plane: usize) -> f64 {
    temperature[4 * plane..4 * plane + 4].iter().sum::<f64>() / 4.0
}

fn radiation_config() -> AmbientRadiationConfig {
    AmbientRadiationConfig {
        max_iterations: 48,
        relaxation: 1.0,
        ..AmbientRadiationConfig::default()
    }
}

/// Explicit numerical surface declaration; this is not measured material data.
fn radiation_patch() -> Result<AmbientRadiationPatch, Box<dyn Error>> {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Scalar {
            value: EMISSIVITY,
            dims: EMISSIVITY_DIMS,
        },
        validity: ValidityDomain::unconstrained().with("T", 295.0, 310.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: Provenance {
            source: "declared synthetic radiative heater calibration".into(),
            license: "internal-example-use".into(),
            artifact: None,
        },
    })?;
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic gray surface".into(),
            phase: "declared fixed surface law".into(),
            process: "numerical calibration example".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )?;
    let emissivity =
        SurfaceEmissivity::from_card(SURFACE, &card, 300.0, SelectionPolicy::SingleClaimOnly)?;
    Ok(AmbientRadiationPatch::new(SURFACE, emissivity, AMBIENT_K)?)
}

impl EnthalpyCalibration {
    /// Build the fixed experiment and generate declared same-model observations.
    pub fn synthetic(cx: &Cx<'_>) -> Result<Self, Box<dyn Error>> {
        Self::build_synthetic(cx, None)
    }

    /// Add known ambient radiation and convection at the heated slab face.
    /// Both pulse amplitudes remain unknown; surface/environment data are fixed.
    pub fn synthetic_with_ambient_radiation(cx: &Cx<'_>) -> Result<Self, Box<dyn Error>> {
        Self::build_synthetic(cx, Some(radiation_patch()?))
    }

    fn build_synthetic(
        cx: &Cx<'_>,
        radiation: Option<AmbientRadiationPatch>,
    ) -> Result<Self, Box<dyn Error>> {
        let (complex, positions) = box_grid([CELLS, 1, 1], [1.0, 0.2, 0.2]);
        let mesh = ConductionMesh::new(complex, positions)?;
        let boundary = if radiation.is_some() {
            ThermalBoundaryBuilder::new(&mesh)
                .region(
                    SURFACE,
                    |face| on_box_face(face.centroid[0], 0.0),
                    ThermalBc::robin(CONVECTION_W_M2_K, AMBIENT_K)?,
                )?
                .adiabatic_remainder()
                .finish()?
        } else {
            ThermalBoundaryBuilder::new(&mesh)
                .adiabatic_remainder()
                .finish()?
        };
        let knots = [
            (-5.0, 295.0, 0.0),
            (0.0, 300.0, 0.0),
            (1.0, 300.0, 1.0),
            (11.0, 310.0, 1.0),
        ]
        .into_iter()
        .map(
            |(h, temperature_k, liquid_mass_fraction)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k,
                liquid_mass_fraction,
                bulk_density_kg_m3: 1.0,
            },
        )
        .collect();
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(
            ContentHash::from_slice(&[29; 32]).ok_or("invalid synthetic chart identity")?,
            knots,
        )?;
        let source_profile = mesh
            .positions()
            .iter()
            .map(|x| (1.0 - x[0] / 0.4).max(0.0))
            .collect();
        let mut experiment = Self {
            mesh,
            curve,
            boundary,
            source_profile,
            radiation,
            material: ConductivityModel::isotropic_declared(0.02)?,
            observations: Vec::new(),
        };
        let temperatures = experiment.forward_temperatures(cx, &TRUTH)?;
        for endpoint in [6, 12, 18, 24] {
            for plane in [0, 1, 2] {
                experiment.observations.push(Observation {
                    endpoint,
                    plane,
                    temperature_k: sensor(&temperatures[endpoint - 1], plane),
                });
            }
        }
        Ok(experiment)
    }

    fn integrator<'a>(
        &'a self,
        cx: &Cx<'_>,
    ) -> Result<EnthalpyBackwardEuler<'a, 'a>, Box<dyn Error>> {
        Ok(EnthalpyBackwardEuler::uniform(
            cx,
            &self.mesh,
            &self.curve,
            1.0,
            EnthalpyBudget {
                max_vertices: 28,
                max_elements: 36,
            },
        )?)
    }

    fn source(&self, point: &[f64], step: usize) -> Result<ScalarField, Box<dyn Error>> {
        let amplitude = point[usize::from(step >= STEPS / 2)];
        Ok(ScalarField::nodal(
            "localized heater density",
            self.mesh.vertex_count(),
            self.source_profile
                .iter()
                .map(|shape| amplitude * shape)
                .collect(),
        )?)
    }

    fn problem<'a>(&'a self, source: &'a ScalarField) -> ConductionProblem<'a> {
        ConductionProblem {
            mesh: &self.mesh,
            boundary: &self.boundary,
            material: &self.material,
            source,
            element_materials: None,
        }
    }

    fn validate(&self, point: &[f64]) -> Result<(), Box<dyn Error>> {
        if point.len() != 2 || point.iter().any(|v| !v.is_finite()) {
            return Err("two finite heater-density parameters required".into());
        }
        if self.observations.len() > 72
            || self.observations.iter().any(|o| {
                o.endpoint == 0
                    || o.endpoint > STEPS
                    || o.plane > CELLS
                    || !o.temperature_k.is_finite()
                    || o.temperature_k <= 0.0
            })
        {
            return Err("observations must fit the fixed mesh, clock and 72-sample cap".into());
        }
        Ok(())
    }

    /// Forward-only route used for data generation and finite-difference checks.
    pub fn forward_temperatures(
        &self,
        cx: &Cx<'_>,
        point: &[f64],
    ) -> Result<Vec<Vec<f64>>, Box<dyn Error>> {
        self.validate(point)?;
        let integrator = self.integrator(cx)?;
        let mut h = vec![INITIAL_H; self.mesh.vertex_count()];
        let mut temperatures = Vec::with_capacity(STEPS);
        for step in 0..STEPS {
            let source = self.source(point, step)?;
            let next = if let Some(patch) = &self.radiation {
                integrator
                    .advance_with_ambient_radiation(
                        cx,
                        self.problem(&source),
                        None,
                        &h,
                        DT,
                        config(),
                        std::slice::from_ref(patch),
                        radiation_config(),
                    )?
                    .conduction
            } else {
                integrator.advance(cx, self.problem(&source), None, &h, DT, config())?
            };
            h = next.specific_enthalpy_j_kg;
            temperatures.push(next.temperature);
        }
        Ok(temperatures)
    }

    /// Mean squared temperature misfit using only the forward production solver.
    pub fn forward_loss(&self, cx: &Cx<'_>, point: &[f64]) -> Result<f64, Box<dyn Error>> {
        if self.observations.is_empty() {
            return Err("temperature observations required".into());
        }
        let history = self.forward_temperatures(cx, point)?;
        Ok(self
            .observations
            .iter()
            .map(|o| 0.5 * (sensor(&history[o.endpoint - 1], o.plane) - o.temperature_k).powi(2))
            .sum::<f64>()
            / f64::from(u32::try_from(self.observations.len())?))
    }

    /// Accumulate observation and history derivatives through one reverse sweep.
    pub fn evaluate(&self, cx: &Cx<'_>, point: &[f64]) -> Result<Evaluation, Box<dyn Error>> {
        self.validate(point)?;
        if self.observations.is_empty() {
            return Err("temperature observations required".into());
        }
        let integrator = self.integrator(cx)?;
        let n = self.mesh.vertex_count();
        let mut h = vec![INITIAL_H; n];
        let mut tape = Vec::with_capacity(STEPS);
        let mut max_step_energy_residual_j = 0.0_f64;
        for step in 0..STEPS {
            let source = self.source(point, step)?;
            let next = if let Some(patch) = &self.radiation {
                integrator
                    .linearize_step_with_ambient_radiation(
                        cx,
                        self.problem(&source),
                        None,
                        &h,
                        DT,
                        config(),
                        std::slice::from_ref(patch),
                        radiation_config(),
                        2 * n,
                    )?
                    .into_transport()
            } else {
                integrator.linearize_step(cx, self.problem(&source), None, &h, DT, config())?
            };
            max_step_energy_residual_j =
                max_step_energy_residual_j.max(next.primal().energy_residual_j.abs());
            h.clone_from(&next.primal().specific_enthalpy_j_kg);
            tape.push(next);
        }
        let mut history_bar = vec![0.0; n];
        let mut gradient = [0.0; 2];
        let mut loss = 0.0;
        let weight = 1.0 / f64::from(u32::try_from(self.observations.len())?);
        for (step, linearization) in tape.iter().enumerate().rev() {
            let mut temperature_bar = vec![0.0; n];
            for observation in self.observations.iter().filter(|o| o.endpoint == step + 1) {
                let error = sensor(&linearization.primal().temperature, observation.plane)
                    - observation.temperature_k;
                loss += 0.5 * weight * error * error;
                for slot in &mut temperature_bar[4 * observation.plane..4 * observation.plane + 4] {
                    *slot += weight * error / 4.0;
                }
            }
            let temperature_h_bar = linearization.temperature_pullback(cx, &temperature_bar)?;
            for (total, seed) in history_bar.iter_mut().zip(temperature_h_bar) {
                *total += seed;
            }
            let adjoint = linearization.pullback(
                cx,
                &history_bar,
                LinearConfig {
                    tolerance: 1e-11,
                    max_iterations: 192,
                    restart: 24,
                },
            )?;
            gradient[usize::from(step >= STEPS / 2)] += adjoint
                .source_density
                .iter()
                .zip(&self.source_profile)
                .map(|(bar, shape)| bar * shape)
                .sum::<f64>();
            history_bar = adjoint.previous_specific_enthalpy;
        }
        let final_step = tape.last().ok_or("empty thermal trajectory")?.primal();
        Ok(Evaluation {
            loss_k2: loss,
            gradient,
            final_enthalpy_j_kg: h,
            final_temperature_k: final_step.temperature.clone(),
            final_liquid_fraction: final_step.liquid_mass_fraction.clone(),
            max_step_energy_residual_j,
        })
    }

    /// SQP inequalities enforce both pulse bounds in physical W/m3 coordinates.
    pub fn sample(&self, cx: &Cx<'_>, point: &[f64]) -> Result<Option<SqpSample>, Box<dyn Error>> {
        self.validate(point)?;
        if point
            .iter()
            .zip(BOUNDS)
            .any(|(x, b)| *x < b[0] || *x > b[1])
        {
            return Ok(None);
        }
        let evaluation = self.evaluate(cx, point)?;
        Ok(Some(SqpSample {
            f: evaluation.loss_k2,
            gradient: evaluation.gradient.to_vec(),
            ce: Vec::new(),
            je: Vec::new(),
            ci: vec![
                BOUNDS[0][0] - point[0],
                point[0] - BOUNDS[0][1],
                BOUNDS[1][0] - point[1],
                point[1] - BOUNDS[1][1],
            ],
            ji: vec![-1.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 1.0],
        }))
    }
}

/// Run the existing optimizer with cumulative evaluation and accepted-step caps.
pub fn fit(
    cx: &Cx<'_>,
    experiment: &EnthalpyCalibration,
) -> Result<(SqpState, SqpRunReport), Box<dyn Error>> {
    let mut evaluate = |point: &[f64]| experiment.sample(cx, point);
    let mut state =
        SqpState::try_new(&START, 6, &mut evaluate, Some(cx)).map_err(|error| error.to_string())?;
    let report = state
        .try_run(&mut evaluate, 1e-8, 50, 160, Some(cx))
        .map_err(|error| error.to_string())?;
    Ok((state, report))
}

fn main() -> Result<(), Box<dyn Error>> {
    let arguments: Vec<_> = std::env::args().skip(1).collect();
    let radiation = match arguments.as_slice() {
        [] => false,
        [flag] if flag == "--ambient-radiation" => true,
        _ => return Err("usage: enthalpy_calibration [--ambient-radiation]".into()),
    };
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 29,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let experiment = if radiation {
            EnthalpyCalibration::synthetic_with_ambient_radiation(&cx)?
        } else {
            EnthalpyCalibration::synthetic(&cx)?
        };
        let initial = experiment.forward_loss(&cx, &START)?;
        let (state, report) = fit(&cx, &experiment)?;
        let final_state = experiment.evaluate(&cx, state.point())?;
        println!(
            "spatial_enthalpy_calibration: vertices=28, steps={STEPS}, observations={}",
            experiment.observations.len()
        );
        if radiation {
            println!(
                "declared_ambient_radiation: emissivity={EMISSIVITY}; reservoir_K={AMBIENT_K}; convection_W_m2_K={CONVECTION_W_M2_K}; fixed_surface={SURFACE}"
            );
        }
        println!(
            "stop={:?}; iterations={}; evaluations={}",
            report.stop,
            state.iterations(),
            state.evaluations()
        );
        println!(
            "initial_loss_K2={initial:.9e}; final_loss_K2={:.9e}",
            state.sample().f
        );
        println!(
            "heater_density_W_m3={:?}; synthetic_reference_W_m3={TRUTH:?}",
            state.point()
        );
        println!(
            "final_liquid_fraction={:?}",
            final_state.final_liquid_fraction
        );
        println!(
            "max_step_energy_residual_J={:.9e}",
            final_state.max_step_energy_residual_j
        );
        if !report.solution.converged {
            return Err("calibration stopped before local KKT convergence".into());
        }
        Ok(())
    })
}
