//! Total endpoint derivatives of the spatial enthalpy/ambient-radiation balance.
//!
//! For patch surface mass B, area-mean weights w, reservoir A and secant s,
//! radiation contributes s(w^T T) B(T-A). Its temperature Jacobian adds the
//! rank-one term s'(w^T T) B(T-A) w^T to the frozen-secant spatial tangent.
//! The enthalpy tangent scales the right factor by dT/dh and the left by dt.
//! Neither the radiative iterations nor the Krylov iterations are differentiated.

use fs_exec::Cx;

use super::{boundary_at_temperature, initial_residual_tolerance_j};
use crate::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution, StepContext,
    adjoint::{EnthalpyAdjointError, EnthalpyStepGradient, EnthalpyStepLinearization},
    poll,
};
use crate::{
    AmbientRadiationConfig, AmbientRadiationPatch, ConductionProblem, LinearConfig,
    ThermalInterfaces, assemble::ASSEMBLY_TILE,
};

/// Physical derivatives of one endpoint objective, with fixed mesh, phase
/// charts, masses, conductivity laws, contact and original convection data.
#[derive(Debug, Clone, PartialEq)]
pub struct EnthalpyRadiationStepGradient {
    /// History/source gradients and checked multiplier of the COMPLETE
    /// radiative enthalpy residual, not the frozen-secant approximation.
    pub transport: EnthalpyStepGradient,
    /// Derivative per second of dt, with endpoint sources/boundaries held fixed.
    /// A time-dependent schedule must add its own endpoint-data derivatives.
    pub time_step_s: f64,
    /// Radiating regions, in the original patch declaration order.
    pub radiation_regions: Vec<String>,
    /// Derivative per kelvin of each black-reservoir temperature. Includes
    /// the secant-coefficient derivative and the explicit reservoir reference.
    pub reservoir_temperatures: Vec<f64>,
    /// Derivative per unit absolute hemispherical emissivity. Claim identity,
    /// validity and query temperature are fixed; at epsilon=1 only feasible
    /// parameter directions are meaningful.
    pub emissivities: Vec<f64>,
}

#[derive(Debug)]
struct PatchDerivative {
    patch: AmbientRadiationPatch,
    faces: Vec<([usize; 3], f64)>,
    secant: f64,
    reservoir_partial: f64,
}

/// A rechecked physical endpoint with a sparse-plus-low-rank enthalpy tangent.
/// The two retained nodal factors per patch are capped explicitly. Face records
/// are disjoint subsets of the admitted boundary; no dense nodal matrix exists.
#[derive(Debug)]
pub struct EnthalpyRadiationStepLinearization<'m> {
    transport: EnthalpyStepLinearization<'m>,
    patches: Vec<PatchDerivative>,
    transport_w: Vec<f64>,
    dt_s: f64,
}

impl<'m> EnthalpyBackwardEuler<'m, '_> {
    /// Solve the physical step and bind its total implicit derivatives.
    /// `max_feedback_entries` caps `2 * vertices * patches` retained scalar
    /// entries, not process RSS. At most 64 distinct patches are admitted.
    /// A tighter derivative request may refuse an insufficient primal solve;
    /// the returned field is never silently replaced during derivative binding.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_step_with_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        step_config: EnthalpyStepConfig,
        patches: &[AmbientRadiationPatch],
        radiation_config: AmbientRadiationConfig,
        max_feedback_entries: usize,
    ) -> Result<EnthalpyRadiationStepLinearization<'m>, EnthalpyAdjointError> {
        poll(cx, 0)?;
        admit_feedback(self.masses.len(), patches.len(), max_feedback_entries)?;
        let accepted = self.advance_with_ambient_radiation(
            cx,
            problem,
            interfaces,
            old_h,
            dt_s,
            step_config,
            patches,
            radiation_config,
        )?;
        self.linearize_accepted_with_ambient_radiation(
            cx,
            problem,
            interfaces,
            old_h,
            dt_s,
            step_config,
            patches,
            accepted.conduction,
            max_feedback_entries,
        )
    }

    /// Bind an existing endpoint without another forward solve. The supplied
    /// specific enthalpy is authoritative input; temperatures, phases, residual
    /// flags and energy telemetry are recomputed rather than trusted.
    ///
    /// `problem.boundary` contains ORIGINAL convection. Radiation is rebuilt
    /// at the actual chart-resolved endpoint, and the full residual must meet
    /// the original physical initial-residual target. Chart/material slope
    /// kinks and emissivity-validity endpoints refuse classical derivatives.
    /// Plateaus strictly inside their enthalpy interval remain differentiable.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn linearize_accepted_with_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        step_config: EnthalpyStepConfig,
        patches: &[AmbientRadiationPatch],
        accepted: EnthalpyStepSolution,
        max_feedback_entries: usize,
    ) -> Result<EnthalpyRadiationStepLinearization<'m>, EnthalpyAdjointError> {
        self.admit_step(cx, problem, old_h, dt_s, step_config)?;
        let n = self.masses.len();
        admit_feedback(n, patches.len(), max_feedback_entries)?;
        if accepted.specific_enthalpy_j_kg.len() != n {
            return Err(EnthalpyAdjointError::InvalidInput(
                "radiative endpoint length mismatch",
            ));
        }
        let target = initial_residual_tolerance_j(
            self,
            cx,
            problem,
            interfaces,
            old_h,
            dt_s,
            step_config,
            patches,
        )?;
        let context = StepContext {
            storage: self,
            cx,
            problem,
            interfaces,
            old: old_h,
            dt: dt_s,
        };
        let temperature = context.temperatures(&accepted.specific_enthalpy_j_kg)?;
        let boundary = boundary_at_temperature(cx, problem, patches, &temperature)?;
        let physical_problem = ConductionProblem {
            boundary: &boundary,
            ..problem
        };
        let mut transport = self.linearize_accepted_with_target(
            cx,
            physical_problem,
            interfaces,
            old_h,
            dt_s,
            step_config,
            accepted,
            Some(target),
        )?;
        let mut retained = Vec::with_capacity(patches.len());
        for patch in patches {
            poll(cx, retained.len())?;
            let region = problem
                .boundary
                .region_names()
                .iter()
                .position(|name| name == patch.region())
                .ok_or(EnthalpyAdjointError::InvalidInput(
                    "unknown radiative derivative region",
                ))?;
            let mut faces = Vec::new();
            let mut area = 0.0;
            let mut integral = 0.0;
            for (slot, face) in problem.mesh.boundary().iter().enumerate() {
                if slot % ASSEMBLY_TILE == 0 {
                    poll(cx, slot)?;
                }
                if problem.boundary.region_for(slot) != Some(region) {
                    continue;
                }
                let vertices = face.vertices.map(|v| v as usize);
                area = finite(area + face.area)?;
                let mean = vertices.iter().map(|&v| temperature[v] / 3.0).sum::<f64>();
                integral = finite(integral + face.area * mean)?;
                faces.push((vertices, face.area));
            }
            if area <= 0.0 {
                return Err(EnthalpyAdjointError::InvalidInput(
                    "radiative derivative has no trace area",
                ));
            }
            let mean = finite(integral / area)?;
            let secant = patch.secant_coefficient_w_m2_k(mean)?;
            let [wall_partial, reservoir_partial] = patch.secant_partials_w_m2_k2(mean)?;
            let mut left = vec![0.0; n];
            let mut mean_weights = vec![0.0; n];
            for (index, (vertices, face_area)) in faces.iter().enumerate() {
                if index % ASSEMBLY_TILE == 0 {
                    poll(cx, index)?;
                }
                for (a, &row) in vertices.iter().enumerate() {
                    mean_weights[row] = finite(mean_weights[row] + face_area / area / 3.0)?;
                    for (b, &column) in vertices.iter().enumerate() {
                        let mass = (face_area / 12.0) * if a == b { 2.0 } else { 1.0 };
                        left[row] = finite(
                            left[row]
                                + dt_s
                                    * wall_partial
                                    * mass
                                    * (temperature[column] - patch.ambient_temperature_k()),
                        )?;
                    }
                }
            }
            // The constitutive derivative scales COLUMNS, including exact zero
            // temperature sensitivity in latent plateau interiors.
            let right = transport.temperature_pullback(cx, &mean_weights)?;
            transport.add_radiation_feedback(cx, left, right)?;
            retained.push(PatchDerivative {
                patch: patch.clone(),
                faces,
                secant,
                reservoir_partial,
            });
        }
        let physical = StepContext {
            problem: physical_problem,
            ..context
        };
        let system = physical.assemble(&temperature)?;
        let mut transport_w = Vec::with_capacity(n);
        for row in 0..n {
            if row % ASSEMBLY_TILE == 0 {
                poll(cx, row)?;
            }
            let (columns, entries) = system.operator.row(row);
            let mut flux = 0.0;
            for (&column, &entry) in columns.iter().zip(entries) {
                flux = finite(entry.mul_add(temperature[column], flux))?;
            }
            transport_w.push(finite(flux - system.load[row])?);
        }
        // Recheck physical heat separately from the numerically combined Robin
        // reference. This is not inferred from storage or from a small residual.
        let original = context.assemble(&temperature)?;
        let dofs = crate::assemble::DofMap::new(problem.boundary, n)?;
        let (energy, _) = crate::solve::energy_balance(
            problem.mesh,
            problem.boundary,
            problem.source,
            &original,
            &dofs,
            &temperature,
        );
        let mut radiation_w = 0.0;
        for entry in &retained {
            let mut area = 0.0;
            let mut integral = 0.0;
            for (index, (vertices, face_area)) in entry.faces.iter().enumerate() {
                if index % ASSEMBLY_TILE == 0 {
                    poll(cx, index)?;
                }
                area = finite(area + face_area)?;
                let mean = vertices.iter().map(|&v| temperature[v] / 3.0).sum::<f64>();
                integral = finite(integral + face_area * mean)?;
            }
            radiation_w =
                finite(radiation_w + area * entry.patch.heat_flux_w_m2(integral / area)?)?;
        }
        let net =
            finite(energy.source_w - energy.neumann_out_w - energy.robin_out_w - radiation_w)?;
        let defect = finite(transport.primal().stored_energy_change_j - dt_s * net)?;
        if defect.abs() > step_config.energy_tolerance_j {
            return Err(EnthalpyError::EnergyBalance {
                residual_j: defect,
                tolerance_j: step_config.energy_tolerance_j,
            }
            .into());
        }
        // Coupled air sees only the original convective ports. The physical
        // tangent above still contains all radiative state feedback.
        transport.retain_convection_boundary(problem.boundary);
        poll(cx, n)?;
        Ok(EnthalpyRadiationStepLinearization {
            transport,
            patches: retained,
            transport_w,
            dt_s,
        })
    }
}

impl<'m> EnthalpyRadiationStepLinearization<'m> {
    /// Rechecked endpoint with the physical combined convection/radiation rows.
    #[must_use]
    pub fn primal(&self) -> &EnthalpyStepSolution {
        self.transport.primal()
    }

    /// Existing enthalpy interface with ALL radiative state feedback retained.
    #[must_use]
    pub const fn transport(&self) -> &EnthalpyStepLinearization<'m> {
        &self.transport
    }

    /// Retain the complete state/history/source derivatives in an existing
    /// enthalpy tape. Only the extra dt/reservoir/emissivity metadata is dropped.
    #[must_use]
    pub fn into_transport(self) -> EnthalpyStepLinearization<'m> {
        self.transport
    }

    /// Convert temperature observation seeds through the exact nodal dT/dh.
    pub fn temperature_pullback(
        &self,
        cx: &Cx<'_>,
        seed: &[f64],
    ) -> Result<Vec<f64>, EnthalpyAdjointError> {
        self.transport.temperature_pullback(cx, seed)
    }

    /// Apply the complete implicit endpoint Jacobian in residual-joule units.
    pub fn apply_jacobian(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
    ) -> Result<Vec<f64>, EnthalpyAdjointError> {
        self.transport.apply_jacobian(cx, direction)
    }

    /// Apply its genuine transpose; a nonuniform patch is not treated as SPD.
    pub fn apply_jacobian_transpose(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
    ) -> Result<Vec<f64>, EnthalpyAdjointError> {
        self.transport.apply_jacobian_transpose(cx, direction)
    }

    /// Solve the complete transpose once, then contract all supported physical
    /// controls. The inner-column cap, actual transpose residual and cancellation
    /// gates are owned by the existing enthalpy pullback. No coefficient or
    /// reservoir reference is frozen when differentiating radiative feedback.
    pub fn pullback(
        &self,
        cx: &Cx<'_>,
        seed: &[f64],
        config: LinearConfig,
    ) -> Result<EnthalpyRadiationStepGradient, EnthalpyAdjointError> {
        let transport = self.transport.pullback(cx, seed, config)?;
        let mut time_step_s = 0.0;
        for (index, (&lambda, &flux)) in transport.adjoint.iter().zip(&self.transport_w).enumerate()
        {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            time_step_s = finite((-lambda).mul_add(flux, time_step_s))?;
        }
        let mut radiation_regions = Vec::with_capacity(self.patches.len());
        let mut reservoir_temperatures = Vec::with_capacity(self.patches.len());
        let mut emissivities = Vec::with_capacity(self.patches.len());
        for entry in &self.patches {
            let ambient = entry.patch.ambient_temperature_k();
            let mut flux_pair = 0.0;
            let mut load_pair = 0.0;
            for (index, (vertices, area)) in entry.faces.iter().enumerate() {
                if index % ASSEMBLY_TILE == 0 {
                    poll(cx, index)?;
                }
                for (a, &row) in vertices.iter().enumerate() {
                    let lambda = transport.adjoint[row];
                    load_pair = finite(lambda.mul_add(area / 3.0, load_pair))?;
                    for (b, &column) in vertices.iter().enumerate() {
                        let mass = (area / 12.0) * if a == b { 2.0 } else { 1.0 };
                        flux_pair = finite(
                            (lambda * mass)
                                .mul_add(self.primal().temperature[column] - ambient, flux_pair),
                        )?;
                    }
                }
            }
            radiation_regions.push(entry.patch.region().to_string());
            reservoir_temperatures.push(finite(
                self.dt_s * (entry.secant * load_pair - entry.reservoir_partial * flux_pair),
            )?);
            emissivities.push(finite(
                -self.dt_s * (entry.secant / entry.patch.emissivity().value()) * flux_pair,
            )?);
        }
        poll(cx, transport.iterations)?;
        Ok(EnthalpyRadiationStepGradient {
            transport,
            time_step_s,
            radiation_regions,
            reservoir_temperatures,
            emissivities,
        })
    }
}

fn admit_feedback(n: usize, patches: usize, limit: usize) -> Result<(), EnthalpyAdjointError> {
    let entries = n.checked_mul(2).and_then(|v| v.checked_mul(patches));
    if patches == 0
        || patches > 64
        || entries.is_none_or(|v| v > limit)
        || entries
            .and_then(|v| v.checked_mul(std::mem::size_of::<f64>()))
            .is_none_or(|bytes| isize::try_from(bytes).is_err())
    {
        return Err(EnthalpyAdjointError::InvalidInput(
            "radiative enthalpy adjoint exceeds its two-vector entry cap or 1..=64 patches",
        ));
    }
    Ok(())
}

fn finite(value: f64) -> Result<f64, EnthalpyAdjointError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(EnthalpyAdjointError::NonFiniteArithmetic)
    }
}
