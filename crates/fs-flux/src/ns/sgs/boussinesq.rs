//! Coupled temperature and planar SGS momentum on the original mixed owner.
//!
//! Each Picard candidate resolves temperature from the SAME old thermal state
//! using its own velocity and eddy diffusivity. That P1 temperature supplies
//! -beta*(T-Tref)*g to the original cell-indexed momentum quadrature. Final
//! native residuals/reactions reassemble both fields at the accepted velocity.
//! This is implicit partitioned BDF1, not a monolithic Newton or energy proof.

use super::{PlanarSgs, PlanarSgsField, PlanarSgsScalarError};
use crate::bdm::{edge_gauss_pub, eval_basis, tri_quad};
use crate::ns::FluxParams;
use crate::ns::boundary::{
    NavierStokesBoundary, NavierStokesControls, NavierStokesStep, StokesBoundaryError,
};
use crate::ns::transport::scalar::{
    ScalarBoundary, ScalarSolveControls, ScalarTransportControls, ScalarTransportError,
    ScalarTransportStep, assembly::ScalarBasis,
};
use fs_math::dd::Dd;
use std::cell::RefCell;

/// Authored linear density response about a fixed positive reference state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarBoussinesqConfig {
    pub reference_temperature_k: f64,
    /// Signed volumetric thermal expansion coefficient [1/K]. Zero disables
    /// buoyancy; negative values represent an explicitly authored anomaly.
    pub thermal_expansion_per_k: f64,
    /// Physical signed gravity [m/s²], e.g. [0,-9.81] with upward y.
    pub gravity_m_s2: [f64; 2],
    /// Maximum |beta*(T-Tref)|, admitted in [0,0.1]. The caller's chosen
    /// small-departure limit is not an approximation-error certificate.
    pub maximum_relative_density_departure: f64,
}

/// Frozen thermal data for one requested physical interval. Temperature is in
/// kelvin, diffusivity in m²/s and source in K/s (Q/(rho Cp)). Density and Cp
/// are fixed. All callbacks must be pure endpoint data throughout the solve.
#[derive(Clone, Copy)]
pub struct PlanarSgsThermalInput<'a> {
    pub previous_temperature_k: &'a [f64],
    pub molecular_diffusivity_m2_s: &'a [f64],
    pub turbulent_prandtl: &'a [f64],
    pub boundaries: &'a [ScalarBoundary<'a>],
    pub source_k_s: &'a dyn Fn(usize, [f64; 2]) -> f64,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct PlanarSgsThermalControls {
    pub flow: NavierStokesControls,
    pub transport: ScalarTransportControls,
    pub thermal_solve: ScalarSolveControls,
}

#[derive(Clone, Debug, PartialEq)]
pub enum PlanarSgsBoussinesqError {
    Invalid(&'static str),
    Temperature {
        dof: usize,
        temperature_k: f64,
        relative_density_departure: f64,
    },
    BoundaryTemperature {
        edge: usize,
        temperature_k: f64,
    },
    SubgridScalar(PlanarSgsScalarError),
    Thermal(ScalarTransportError),
    Flow(StokesBoundaryError),
    Cancelled,
}
impl core::fmt::Display for PlanarSgsBoussinesqError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar SGS Boussinesq: {self:?}")
    }
}
impl std::error::Error for PlanarSgsBoussinesqError {}
type Result<T> = core::result::Result<T, PlanarSgsBoussinesqError>;

/// Joint accepted endpoint. The original momentum and scalar residuals retain
/// their separate tolerances and budgets; neither history is modified.
pub struct PlanarSgsBoussinesqStep {
    pub flow: NavierStokesStep,
    pub temperature: ScalarTransportStep,
    pub sgs: PlanarSgsField,
    pub turbulent_diffusivity_m2_s: Vec<f64>,
    pub maximum_relative_density_departure: f64,
    /// Integral of accepted buoyancy acceleration dot BDM velocity [m⁴/s³].
    /// Multiply by constant rho and depth for watts. Not a heat source,
    /// gravitational-potential change or a whole-step mechanical balance.
    pub buoyancy_power_m4_s3: f64,
}

struct Candidate {
    temperature: ScalarTransportStep,
    sgs: PlanarSgsField,
    turbulent_diffusivity: Vec<f64>,
    departure: f64,
}

impl PlanarBoussinesqConfig {
    fn validate(self) -> Result<()> {
        if !self.reference_temperature_k.is_finite()
            || self.reference_temperature_k <= 0.0
            || !self.thermal_expansion_per_k.is_finite()
            || self.gravity_m_s2.iter().any(|g| !g.is_finite())
            || !self.maximum_relative_density_departure.is_finite()
            || !(0.0..=0.1).contains(&self.maximum_relative_density_departure)
        {
            return Err(PlanarSgsBoussinesqError::Invalid(
                "positive reference kelvin, finite beta/gravity and departure limit in [0,0.1] required",
            ));
        }
        Ok(())
    }

    fn departure(self, temperature: f64) -> f64 {
        self.thermal_expansion_per_k * (temperature - self.reference_temperature_k)
    }

    fn admits(self, temperature: f64) -> bool {
        temperature.is_finite()
            && temperature > 0.0
            && self.departure(temperature).is_finite()
            && self.departure(temperature).abs() <= self.maximum_relative_density_departure
    }

    fn admit_field(self, values: &[f64]) -> Result<f64> {
        let mut maximum = 0.0_f64;
        for (dof, &temperature) in values.iter().enumerate() {
            let departure = self.departure(temperature).abs();
            if !self.admits(temperature) {
                return Err(PlanarSgsBoussinesqError::Temperature {
                    dof,
                    temperature_k: temperature,
                    relative_density_departure: departure,
                });
            }
            maximum = maximum.max(departure);
        }
        Ok(maximum)
    }
}

impl PlanarSgs<'_, '_> {
    /// Couple implicit planar SGS flow to actual P1-DG temperature evolution.
    ///
    /// At each ORIGINAL mixed Picard candidate, the original sparse thermal
    /// BDF1 solve starts from immutable old temperature, with alpha=alpha_m+nu_t/Pr_t.
    /// Its affine temperature supplies buoyancy at each original momentum
    /// quadrature point. The same candidate's molecular+SGS stress, convection,
    /// scalar transport and buoyancy are reassembled for final acceptance.
    /// No independent outer coupling loop or lagged buoyancy force is used.
    ///
    /// `force` is additional nonbuoyant acceleration. Pressure and open
    /// traction use p_reduced/rho = p_physical/rho - g.(x-x_ref); do not add
    /// reference-density gravity again. This method supplies no hydrostatic
    /// datum. Gravity has its physical sign, so hot fluid rises opposite g
    /// for positive beta. No buoyancy work, molecular/SGS dissipation or
    /// numerical loss is silently deposited as heat.
    ///
    /// Molecular viscosity, thermal diffusivity, density and Cp stay frozen;
    /// only the authored linear density anomaly enters momentum. The separate
    /// source-backed laminar enthalpy owner is unchanged. Temperature BDF1 is
    /// not bound-preserving: inadmissible kelvin/departure or sparse-solve loss
    /// refuses the WHOLE joint interval. Picard convergence is conditional on
    /// coupling strength and timestep. Original flow/thermal work caps remain.
    /// Cancellation is checked at thermal and native assembly boundaries;
    /// the existing momentum linear solve is not made preemptible here.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf1_boussinesq_step(
        &self,
        params: FluxParams,
        config: PlanarBoussinesqConfig,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        thermal: PlanarSgsThermalInput<'_>,
        timestep_s: f64,
        controls: PlanarSgsThermalControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PlanarSgsBoussinesqStep> {
        if cancelled() {
            return Err(PlanarSgsBoussinesqError::Cancelled);
        }
        config.validate()?;
        let mesh = self.system.mesh;
        if thermal.previous_temperature_k.len() != 3 * mesh.tris.len() {
            return Err(PlanarSgsBoussinesqError::Invalid(
                "three original P1 temperatures per triangle required",
            ));
        }
        config.admit_field(thermal.previous_temperature_k)?;
        let mut bases = Vec::with_capacity(mesh.tris.len());
        for (cell, tri) in mesh.tris.iter().enumerate() {
            if cancelled() {
                return Err(PlanarSgsBoussinesqError::Cancelled);
            }
            bases.push(ScalarBasis::new(
                tri.map(|v| mesh.verts[v]),
                mesh.areas[cell],
            ));
        }
        // Admit the actual affine boundary trace, including its extrema. A
        // flux callback supplies K*m/s, not temperature; only its donor is read.
        for &boundary in thermal.boundaries {
            let (edge, value) = match boundary {
                ScalarBoundary::Value { edge, value } => (edge, Some(value)),
                ScalarBoundary::Flux { edge, inflow, .. } => (edge, inflow),
            };
            let Some(value) = value else { continue };
            let geometry = mesh
                .edges
                .get(edge)
                .ok_or(PlanarSgsBoussinesqError::Invalid(
                    "thermal boundary edge out of range",
                ))?;
            if cancelled() {
                return Err(PlanarSgsBoussinesqError::Cancelled);
            }
            let samples =
                edge_gauss_pub(mesh.verts[geometry.verts.0], mesh.verts[geometry.verts.1])
                    .map(|(p, _)| value(p));
            for fraction in [
                0.0,
                1.0,
                0.5 - 0.5 * 3.0_f64.sqrt(),
                0.5 + 0.5 * 3.0_f64.sqrt(),
            ] {
                let temperature = (1.0 - fraction) * samples[0] + fraction * samples[1];
                if !config.admits(temperature) {
                    return Err(PlanarSgsBoussinesqError::BoundaryTemperature {
                        edge,
                        temperature_k: temperature,
                    });
                }
            }
        }
        let interpolate = |values: &[f64], cell: usize, point| {
            let phi = bases[cell].value(point);
            (0..3).map(|i| phi[i] * values[3 * cell + i]).sum::<f64>()
        };
        let prepared = RefCell::<Option<Candidate>>::new(None);
        let failure = RefCell::<Option<PlanarSgsBoussinesqError>>::new(None);
        let cancel = RefCell::new(cancelled);
        let prepare = |velocity: &[f64]| {
            let result = (|| -> Result<Candidate> {
                let coupled = self
                    .prepare_scalar_transport(
                        &velocity[..self.system.n_u],
                        thermal.molecular_diffusivity_m2_s,
                        thermal.turbulent_prandtl,
                        thermal.boundaries,
                        controls.transport,
                        || (*cancel.borrow_mut())(),
                    )
                    .map_err(|error| match error {
                        PlanarSgsScalarError::Scalar(ScalarTransportError::Cancelled) => {
                            PlanarSgsBoussinesqError::Cancelled
                        }
                        error => PlanarSgsBoussinesqError::SubgridScalar(error),
                    })?;
                let temperature = coupled
                    .transport
                    .bdf1_step(
                        thermal.previous_temperature_k,
                        timestep_s,
                        thermal.source_k_s,
                        controls.thermal_solve,
                        || (*cancel.borrow_mut())(),
                    )
                    .map_err(|error| match error {
                        ScalarTransportError::Cancelled => PlanarSgsBoussinesqError::Cancelled,
                        error => PlanarSgsBoussinesqError::Thermal(error),
                    })?;
                let departure = config.admit_field(&temperature.values)?;
                Ok(Candidate {
                    temperature,
                    sgs: coupled.sgs,
                    turbulent_diffusivity: coupled.turbulent_diffusivity_m2_s,
                    departure,
                })
            })();
            match result {
                Ok(candidate) => {
                    let viscosities = candidate
                        .sgs
                        .cells
                        .iter()
                        .map(|cell| cell.effective_viscosity_m2_s)
                        .collect();
                    *prepared.borrow_mut() = Some(candidate);
                    Ok(viscosities)
                }
                Err(error) => {
                    *failure.borrow_mut() = Some(error);
                    Err(StokesBoundaryError::InvalidInput(
                        "coupled SGS temperature refused",
                    ))
                }
            }
        };
        let load = |cell, point| {
            let prepared = prepared.borrow();
            let values = &prepared
                .as_ref()
                .expect("prepared cell force")
                .temperature
                .values;
            let departure = config.departure(interpolate(values, cell, point));
            let external = force(point);
            core::array::from_fn(|a| external[a] - departure * config.gravity_m_s2[a])
        };
        let flow = self
            .system
            .bdf1_navier_stokes_mixed_prepared_cell_forces_checked(
                params,
                &prepare,
                &load,
                boundaries,
                previous_velocity,
                timestep_s,
                controls.flow,
                &|| {
                    if failure.borrow().is_some() {
                        return false;
                    }
                    if (*cancel.borrow_mut())() {
                        *failure.borrow_mut() = Some(PlanarSgsBoussinesqError::Cancelled);
                        false
                    } else {
                        true
                    }
                },
            );
        if let Some(error) = failure.into_inner() {
            return Err(error);
        }
        let flow = flow.map_err(PlanarSgsBoussinesqError::Flow)?;
        // The original successful endpoint/reaction assembly is the LAST
        // preparation. Return its thermal solve, never a previous Picard field.
        let candidate = prepared
            .into_inner()
            .ok_or(PlanarSgsBoussinesqError::Invalid(
                "missing accepted thermal endpoint",
            ))?;
        let velocity = &flow.solution.field.x[..self.system.n_u];
        let mut power = Dd::from_f64(0.0);
        for (cell, tri) in mesh.tris.iter().enumerate() {
            if (*cancel.borrow_mut())() {
                return Err(PlanarSgsBoussinesqError::Cancelled);
            }
            for (point, weight) in tri_quad(tri.map(|v| mesh.verts[v]), mesh.areas[cell]) {
                let departure =
                    config.departure(interpolate(&candidate.temperature.values, cell, point));
                let mut u = [0.0; 2];
                for local in 0..6 {
                    let dof = 2 * mesh.tri_edges[cell][local / 2].0 + local % 2;
                    let basis = eval_basis(&self.system.bases[cell], local, point);
                    for a in 0..2 {
                        u[a] += velocity[dof] * basis[a];
                    }
                }
                power = power
                    + Dd::from_f64(
                        -weight
                            * departure
                            * (u[0] * config.gravity_m_s2[0] + u[1] * config.gravity_m_s2[1]),
                    );
            }
        }
        let power = power.to_f64();
        if !power.is_finite() {
            return Err(PlanarSgsBoussinesqError::Invalid(
                "unrepresentable buoyancy power",
            ));
        }
        if (*cancel.borrow_mut())() {
            return Err(PlanarSgsBoussinesqError::Cancelled);
        }
        Ok(PlanarSgsBoussinesqStep {
            flow,
            temperature: candidate.temperature,
            sgs: candidate.sgs,
            turbulent_diffusivity_m2_s: candidate.turbulent_diffusivity,
            maximum_relative_density_departure: candidate.departure,
            buoyancy_power_m4_s3: power,
        })
    }
}
