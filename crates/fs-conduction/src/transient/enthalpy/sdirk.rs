//! Alexander's two-stage, stiffly accurate SDIRK2 on the original h residual.
//!
//! With gamma=1-1/sqrt(2), M(h1-h0)=gamma*dt*Q(h1) and
//! M(h2-h0)=dt*((1-gamma)*Q(h1)+gamma*Q(h2)). Both stages retain the
//! ORIGINAL h0. A dual-joule carry enters the existing residual without ever
//! evaluating an extrapolated fictitious enthalpy as a material state.
//! Stage heat quadrature is distinct from endpoint power and from internal
//! stage carry. This method is second order for smooth admitted histories;
//! phase/source corners and forcing jumps need event resolution/refinement.
//! L-stability does not establish unconditional positivity or an entropy theorem.

use super::*;
use super::heterogeneous::EnthalpyRadiation;
use crate::interface::temperature::{TemperatureContactEvaluation, TemperatureContacts};
use crate::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig;
use fs_math::dd::Dd;

/// Diagonal coefficient and first physical stage time, matching the existing
/// SDIRK2 gas-chemistry tableau without a dependency on that peer domain.
pub const SDIRK2_GAMMA: f64 = 1.0 - core::f64::consts::FRAC_1_SQRT_2;

/// Physical laws sampled at one tableau time, t+gamma*dt or t+dt.
/// Each prepared radiation owner must retain this exact mesh and boundary.
/// Fixed interfaces and source-dependent contacts are mutually exclusive.
#[derive(Clone, Copy)]
pub struct EnthalpySdirk2Input<'a> {
    pub problem: ConductionProblem<'a>,
    pub interfaces: Option<&'a ThermalInterfaces>,
    /// None or an original pointwise/material ambient/enclosure owner.
    /// Mean-temperature outer radiation iterations are explicitly unsupported.
    pub radiation: EnthalpyRadiation<'a>,
    pub contacts: Option<(&'a TemperatureContacts<'a, 'a>, ContactEnthalpyNewtonConfig)>,
}

/// An actual implicit stage, with original endpoint sources and fluxes.
/// Its power fields are instantaneous physical watts, never weighted averages.
/// The second stage is the accepted endpoint; the first is not an accepted
/// subinterval and must not be published or charged as one.
#[derive(Debug)]
pub struct EnthalpySdirk2Stage {
    pub specific_enthalpy_j_kg: Vec<f64>,
    pub temperature: Vec<f64>,
    pub liquid_mass_fraction: Vec<f64>,
    pub source_w: f64,
    pub neumann_out_w: f64,
    pub robin_out_w: f64,
    pub radiation_out_w: f64,
    /// Physical prescribed-bath reaction, including the original storage and
    /// prior stage carry. That carry is not an additional external source.
    pub dirichlet_in_w: f64,
    pub robin_fluxes: Vec<RobinFlux>,
    pub contact_fluxes: Vec<InterfaceFlux>,
    pub pointwise_radiation: Option<PointwiseRadiationReport>,
    pub pointwise_enclosure_radiation: Option<PointwiseEnclosureReport>,
    pub physical_contacts: Option<TemperatureContactEvaluation>,
    pub newton: NewtonReport,
    /// Free-node residual norm of THIS stage equation, in joules.
    pub stage_residual_norm_j: f64,
    /// Global stage residual after the internal carry is subtracted, in joules.
    pub stage_energy_residual_j: f64,
    // Authentic assembled net input, including internal contact/conduction
    // and prescribed-node reaction; retained for final tableau verification.
    net_nodal_input_w: Vec<f64>,
}

/// One accepted interval with both authentic stages and integrated heat.
/// No stage carry is included in these external physical joule accounts.
#[derive(Debug)]
pub struct EnthalpySdirk2Step {
    pub duration_s: f64,
    pub first_stage: EnthalpySdirk2Stage,
    pub second_stage: EnthalpySdirk2Stage,
    pub source_j: f64,
    pub neumann_out_j: f64,
    pub robin_out_j: f64,
    pub radiation_out_j: f64,
    pub dirichlet_in_j: f64,
    /// Independently reconstructed original-to-final reference-mass inventory.
    pub stored_energy_change_j: f64,
    pub energy_residual_j: f64,
    /// Original-to-final balance with both physical stage rates. Constrained
    /// rows include their actual stage reservoir reactions and are zero to
    /// arithmetic accuracy. No residual is inferred from global cancellation.
    pub nodal_residual_j: Vec<f64>,
}

impl EnthalpyBackwardEuler<'_, '_> {
    /// Advance with two nonlinear stages on the existing storage, spatial
    /// assembly, exact tangents and bounded Newton/FGMRES solvers.
    ///
    /// The caller supplies physical data at the two tableau times. Source k,
    /// caloric state, material contact and pointwise reflection are reevaluated
    /// at each original stage h. Both nodewise and global joule gates apply;
    /// a failed second stage returns no partially accepted interval.
    /// Prescribed temperatures must have unique h on their original charts.
    /// A second-order adjoint and unconditional phase monotonicity are not
    /// supplied by the first-order endpoint adjoint or by L-stability.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sdirk2(
        &self, cx: &Cx<'_>, inputs: [EnthalpySdirk2Input<'_>; 2],
        old_h: &[f64], dt_s: f64, config: EnthalpyStepConfig,
    ) -> Result<EnthalpySdirk2Step, EnthalpyError> {
        validate_work(self.masses.len(), dt_s, config)?;
        config.newton.validate().map_err(EnthalpyError::Newton)?;
        config.max_newton_iterations.checked_mul(2)
            .and_then(|n| n.checked_mul(config.newton.linear_restart))
            .and_then(|n| n.checked_mul(config.newton.max_linear_cycles))
            .ok_or(EnthalpyError::InvalidInput("SDIRK2 work budget overflow"))?;
        self.masses.len().checked_mul(256).filter(|&bytes| isize::try_from(bytes).is_ok())
            .ok_or(EnthalpyError::InvalidInput("SDIRK2 stage storage overflow"))?;
        let diagonal_dt = finite(SDIRK2_GAMMA * dt_s)?;
        let first_weight = finite((1.0 - SDIRK2_GAMMA) * dt_s)?;
        if diagonal_dt <= 0.0 || first_weight <= 0.0 {
            return Err(EnthalpyError::InvalidInput("positive SDIRK2 weights must be representable"));
        }
        // Admit both samplings before the first costly Newton stage.
        for input in inputs {
            self.admit_step(cx, input.problem, old_h, diagonal_dt, config)?;
            self.sdirk_radiation(cx, input)?;
            if let Some((_, contact)) = input.contacts {
                if input.interfaces.is_some() || contact.step != config {
                    return Err(EnthalpyError::InvalidInput(
                        "SDIRK2 requires exclusive contact ownership and unchanged step controls"));
                }
                contact.validate()?;
            }
        }
        let first_stage = self.sdirk_stage(cx, inputs[0], old_h, diagonal_dt, config, None, None)?;
        let mut carry = Vec::with_capacity(old_h.len());
        for (vertex, &input) in first_stage.net_nodal_input_w.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 { poll(cx, vertex)?; }
            let joules = Dd::from_f64(first_weight) * Dd::from_f64(input);
            checked(joules)?;
            carry.push(joules);
        }
        let second_stage = self.sdirk_stage(cx, inputs[1], old_h, diagonal_dt, config,
            Some(&carry), Some(&first_stage.specific_enthalpy_j_kg))?;
        let integral = |a: f64, b: f64| checked(Dd::from_f64(first_weight) * Dd::from_f64(a)
            + Dd::from_f64(diagonal_dt) * Dd::from_f64(b));
        let source_j = integral(first_stage.source_w, second_stage.source_w)?;
        let neumann_out_j = integral(first_stage.neumann_out_w, second_stage.neumann_out_w)?;
        let robin_out_j = integral(first_stage.robin_out_w, second_stage.robin_out_w)?;
        let radiation_out_j = integral(first_stage.radiation_out_w, second_stage.radiation_out_w)?;
        let dirichlet_in_j = integral(first_stage.dirichlet_in_w, second_stage.dirichlet_in_w)?;
        let mut stored = Dd::ZERO;
        let mut nodal_residual_j = Vec::with_capacity(old_h.len());
        for (i, (&old, &new)) in old_h.iter().zip(&second_stage.specific_enthalpy_j_kg).enumerate() {
            if i % ASSEMBLY_TILE == 0 { poll(cx, i)?; }
            let delta = Dd::from_f64(self.masses[i]) * (Dd::from_f64(new) - Dd::from_f64(old));
            stored = stored + delta;
            nodal_residual_j.push(checked(delta - carry[i]
                - Dd::from_f64(diagonal_dt) * Dd::from_f64(second_stage.net_nodal_input_w[i]))?);
        }
        let residual = checked(stored - Dd::from_f64(source_j) - Dd::from_f64(dirichlet_in_j)
            + Dd::from_f64(neumann_out_j) + Dd::from_f64(robin_out_j) + Dd::from_f64(radiation_out_j))?;
        check_residual(&nodal_residual_j, config.energy_tolerance_j)?;
        check_energy(residual, config.energy_tolerance_j)?;
        poll(cx, old_h.len())?;
        Ok(EnthalpySdirk2Step { duration_s: dt_s, first_stage, second_stage, source_j,
            neumann_out_j, robin_out_j, radiation_out_j, dirichlet_in_j,
            stored_energy_change_j: checked(stored)?, energy_residual_j: residual, nodal_residual_j })
    }

    fn sdirk_radiation<'a>(&self, cx: &Cx<'_>, input: EnthalpySdirk2Input<'a>)
        -> Result<Option<PointwiseBoundary<'a>>, EnthalpyError>
    {
        let p = input.problem;
        let (matches, native) = match input.radiation {
            EnthalpyRadiation::None => return Ok(None),
            EnthalpyRadiation::PointwiseAmbient(owner) =>
                (owner.matches(p.mesh, p.boundary), PointwiseBoundary::Ambient(owner)),
            EnthalpyRadiation::PointwiseEnclosure(owner) =>
                (owner.matches(p.mesh, p.boundary), PointwiseBoundary::Enclosure(owner)),
            EnthalpyRadiation::MaterialAmbient(owner) =>
                (owner.matches(p.mesh, p.boundary)
                    && owner.matches_curves(cx, |v| self.curve_for_vertex(v))?,
                    PointwiseBoundary::MaterialAmbient(owner)),
            EnthalpyRadiation::MaterialEnclosure(owner) =>
                (owner.matches(p.mesh, p.boundary)
                    && owner.matches_curves(cx, |v| self.curve_for_vertex(v))?,
                    PointwiseBoundary::MaterialEnclosure(owner)),
            EnthalpyRadiation::Ambient { .. } | EnthalpyRadiation::Enclosure { .. } =>
                return Err(EnthalpyError::InvalidInput(
                    "SDIRK2 uses the coupled pointwise radiation owners, not mean-temperature outer iterations")),
        };
        if !matches { return Err(EnthalpyError::InvalidInput(
            "SDIRK2 radiation must retain exact stage mesh, boundary and caloric objects")); }
        Ok(Some(native))
    }

    #[allow(clippy::too_many_arguments)]
    fn sdirk_stage(&self, cx: &Cx<'_>, input: EnthalpySdirk2Input<'_>, old_h: &[f64],
        dt: f64, config: EnthalpyStepConfig, carry: Option<&[Dd]>, guess: Option<&[f64]>)
        -> Result<EnthalpySdirk2Stage, EnthalpyError>
    {
        let radiation = self.sdirk_radiation(cx, input)?;
        let base = StepContext { storage: self, cx, problem: input.problem,
            interfaces: input.interfaces, old: old_h, dt, pointwise_radiation: radiation };
        let (h, newton, contacts) = if let Some((owner, contact_config)) = input.contacts {
            let (h, report, physical) = self.solve_with_temperature_contacts(cx, input.problem,
                owner, old_h, dt, contact_config, radiation, carry, guess)?;
            (h, report, Some(physical))
        } else {
            let (h, report) = base.solve_newton_from(config, carry, guess)?;
            (h, report, None)
        };
        let context = StepContext { interfaces: contacts.as_ref()
            .map(TemperatureContactEvaluation::frozen_interfaces).or(input.interfaces), ..base };
        let mut stage = finish_stage(&context, h, newton, carry, config.energy_tolerance_j)?;
        if let Some((_, contact_config)) = input.contacts {
            if stage.stage_residual_norm_j > contact_config.residual_tolerance_j {
                return Err(EnthalpyError::EndpointResidual { residual_j: stage.stage_residual_norm_j,
                    tolerance_j: contact_config.residual_tolerance_j });
            }
        }
        stage.physical_contacts = contacts;
        Ok(stage)
    }
}

fn checked(value: Dd) -> Result<f64, EnthalpyError> {
    if value.hi.is_finite() && value.lo.is_finite() { finite(value.to_f64()) }
    else { Err(EnthalpyError::NonFiniteArithmetic) }
}
fn check_energy(residual: f64, tolerance_j: f64) -> Result<(), EnthalpyError> {
    if residual.abs() > tolerance_j {
        Err(EnthalpyError::EnergyBalance { residual_j: residual, tolerance_j })
    } else { Ok(()) }
}
fn check_residual(values: &[f64], tolerance_j: f64) -> Result<f64, EnthalpyError> {
    let norm = finite(fs_solver::norm2(values))?;
    if norm > tolerance_j { Err(EnthalpyError::EndpointResidual { residual_j: norm, tolerance_j }) }
    else { Ok(norm) }
}

fn finish_stage(context: &StepContext<'_, '_, '_, '_>, h: Vec<f64>, newton: NewtonReport,
    carry: Option<&[Dd]>, tolerance_j: f64) -> Result<EnthalpySdirk2Stage, EnthalpyError>
{
    context.check_prescribed_values(&h)?;
    let temperature = context.temperatures(&h)?;
    let scales = context.phase_scales(&h)?;
    let system = context.assemble(&temperature, scales.as_deref())?;
    let radiation = context.pointwise_radiation
        .map(|owner| owner.load_and_report(context.cx, &temperature, &h)).transpose()?;
    let mut liquid_mass_fraction = Vec::with_capacity(h.len());
    let mut net_nodal_input_w = Vec::with_capacity(h.len());
    let mut nodal_residual = Vec::with_capacity(h.len());
    let mut storage = Dd::ZERO;
    let mut carried = Dd::ZERO;
    let mut reaction = Dd::ZERO;
    for (i, &value) in h.iter().enumerate() {
        if i % ASSEMBLY_TILE == 0 { poll(context.cx, i)?; }
        let delta = Dd::from_f64(context.storage.masses[i])
            * (Dd::from_f64(value) - Dd::from_f64(context.old[i]));
        let prior = carry.map_or(Dd::ZERO, |carry| carry[i]);
        storage = storage + delta; carried = carried + prior;
        liquid_mass_fraction.push(context.storage.curve_for_vertex(i).state_at_specific_enthalpy(value)
            .map_err(|source| EnthalpyError::Phase { vertex: i, source })?.liquid_mass_fraction());
        let (columns, entries) = system.operator.row(i);
        let outward: f64 = columns.iter().zip(entries).map(|(&j, &v)| v * temperature[j]).sum();
        let outward = match &radiation { Some((load, _)) => finite(outward + load[i])?, None => finite(outward)? };
        let known = finite(system.load[i] - outward)?;
        let residual = delta - prior - Dd::from_f64(context.dt) * Dd::from_f64(known);
        if context.is_prescribed(i) {
            reaction = reaction + residual / Dd::from_f64(context.dt);
            net_nodal_input_w.push(checked((delta - prior) / Dd::from_f64(context.dt))?);
            nodal_residual.push(0.0);
        } else {
            net_nodal_input_w.push(known);
            nodal_residual.push(checked(residual)?);
        }
    }
    let stage_residual_norm_j = check_residual(&nodal_residual, tolerance_j)?;
    let dofs = DofMap::new_allow_fully_prescribed(context.problem.boundary, h.len())?;
    let (energy, robin_fluxes) = energy_balance(context.storage.mesh, context.problem.boundary,
        context.problem.source, &system, &dofs, &temperature);
    let radiation_out_w = radiation.as_ref().map_or(0.0, |(_, report)| report.outward_heat_w());
    let net = Dd::from_f64(energy.source_w) + reaction - Dd::from_f64(energy.neumann_out_w)
        - Dd::from_f64(energy.robin_out_w) - Dd::from_f64(radiation_out_w);
    let stage_energy_residual_j = checked(storage - carried - Dd::from_f64(context.dt) * net)?;
    check_energy(stage_energy_residual_j, tolerance_j)?;
    let contact_fluxes = match context.interfaces {
        Some(interfaces) => interfaces.fluxes(&temperature)?, None => Vec::new(),
    };
    let (pointwise_radiation, pointwise_enclosure_radiation) = radiation
        .map_or((None, None), |(_, report)| report.into_reports());
    poll(context.cx, h.len())?;
    Ok(EnthalpySdirk2Stage { specific_enthalpy_j_kg: h, temperature, liquid_mass_fraction,
        source_w: finite(energy.source_w)?, neumann_out_w: finite(energy.neumann_out_w)?,
        robin_out_w: finite(energy.robin_out_w)?, radiation_out_w: finite(radiation_out_w)?,
        dirichlet_in_w: checked(reaction)?, robin_fluxes, contact_fluxes, pointwise_radiation,
        pointwise_enclosure_radiation, physical_contacts: None, newton, stage_residual_norm_j,
        stage_energy_residual_j, net_nodal_input_w })
}
