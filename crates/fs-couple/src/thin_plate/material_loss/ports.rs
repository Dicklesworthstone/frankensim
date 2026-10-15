//! Spatial force/displacement ports on ONE source-derived plate eigenbasis.
//! No per-footprint eigensolves: signs and degenerate eigenspaces must be shared.

use super::*;
use crate::plate_viscosity::{PlateViscosity, PlateViscosityLimits, ModalPlateViscosity, ViscousPlateImage};
use fs_exec::CancelGate;

impl ResolvedPlateChart {
    /// Reduce the body once and project additional nodal force footprints on
    /// that SAME basis, with the selected regional material loss preserved.
    ///
    /// The first result contains ordinary radiators whose reference drive is
    /// `options.unit_force_weights`. The second contains one mass-normalized
    /// attachment vector per additional footprint, in input order. For vector
    /// `b`, displacement is `b^T q` and modal force is `b F`; this transpose is
    /// essential for two-way work and cross-port reciprocity. Use the vectors
    /// only with these returned modes, never a separately solved body.
    ///
    /// Every footprint obeys the existing nodal unit-total-force convention.
    /// Signed shape coefficients are retained, including acoustically silent
    /// modes. Equal and opposite forces at distinct ports can excite rocking
    /// motion even though their net force vanishes. All ports share one signed
    /// volume/radiation projection, not one independently tuned output filter.
    ///
    /// `max_projection_terms` bounds additional footprint construction and
    /// projection by p*(3*triangles + nodes*(1+requested_modes)); the existing
    /// geometry/modal budgets still govern the single underlying eigensolve.
    /// This is a frozen linear DKT image, not a rigid bridge, joint, shell,
    /// torque connector, or full nonproportional damping model.
    ///
    /// # Errors
    /// Invalid footprints, projection budget/overflow, source-loss/domain or
    /// ordinary geometry/modal refusal. No partial set of ports is returned.
    pub fn certified_radiators_with_force_ports(
        &self,
        options: &PlateChartRadiation,
        additional_force_weights: &[Vec<f64>],
        max_projection_terms: usize,
    ) -> Result<(Vec<CompactBody>, Vec<Vec<f64>>), AcousticRealizeError> {
        reduce(self, options, additional_force_weights, max_projection_terms)
    }
}

pub(super) fn reduce(
    specimen: &ResolvedPlateChart,
    options: &PlateChartRadiation,
    additional: &[Vec<f64>],
    max_projection_terms: usize,
) -> Result<(Vec<CompactBody>, Vec<Vec<f64>>), AcousticRealizeError> {
    reduce_with_tensor(specimen, options, additional, max_projection_terms, None)
        .map(|(radiators, ports, _)| (radiators, ports))
}

impl ResolvedPlateChart {
    /// Project the complete source-derived viscous tensor on the SAME modes and
    /// physical ports. Radiators retain zero scalar loss: the returned passive
    /// connections own ALL material damping, including cross-mode reactions.
    /// No selected isotropic or authored damping can be counted a second time.
    /// The original eigensolver is synchronous; gates bracket that bounded call.
    pub fn certified_viscous_radiators(&self, options: &PlateChartRadiation,
        viscosity: &PlateViscosity, limits: PlateViscosityLimits, gate: &CancelGate)
        -> Result<ViscousPlateImage, AcousticRealizeError>
    {
        self.certified_viscous_radiators_with_force_ports(options, &[], 0, viscosity, limits, gate)
            .map(|(image, _)| image)
    }

    /// Keep full tensor damping and multiple physical force footprints on one
    /// eigenbasis. Additional force columns use the original unit-total-force
    /// and projection-budget contract; no new modal solve or normalization.
    pub fn certified_viscous_radiators_with_force_ports(&self, options: &PlateChartRadiation,
        additional_force_weights: &[Vec<f64>], max_force_projection_terms: usize,
        viscosity: &PlateViscosity, limits: PlateViscosityLimits, gate: &CancelGate)
        -> Result<(ViscousPlateImage, Vec<Vec<f64>>), AcousticRealizeError>
    {
        let (radiators, ports, damping) = reduce_with_tensor(self, options,
            additional_force_weights, max_force_projection_terms, Some((viscosity, limits, gate)))?;
        Ok((ViscousPlateImage { source: viscosity.clone(), radiators,
            damping: damping.expect("selected tensor projection returns its full matrix") }, ports))
    }
}

fn reduce_with_tensor(
    specimen: &ResolvedPlateChart, options: &PlateChartRadiation,
    additional: &[Vec<f64>], max_projection_terms: usize,
    tensor: Option<(&PlateViscosity, PlateViscosityLimits, &CancelGate)>,
) -> Result<(Vec<CompactBody>, Vec<Vec<f64>>, Option<ModalPlateViscosity>), AcousticRealizeError> {
    if let Some((law, limits, gate)) = tensor {
        law.preflight(specimen, options.n_modes, limits, gate)?;
        if options.damping_ratio != 0.0 { return Err(refuse("full material viscosity requires zero authored modal damping")); }
    }
    let window = options.eigenvalue_window;
    if options.n_modes == 0 || !options.damping_ratio.is_finite() || options.damping_ratio < 0.0
        || !window.0.is_finite() || !window.1.is_finite() || window.0 < 0.0 || window.1 <= window.0 {
        return Err(refuse("plate needs positive modes, nonnegative finite damping and an ordered finite window"));
    }
    if specimen.bending_viscosity().is_some() && options.damping_ratio != 0.0 {
        return Err(refuse("source-bound regional bending loss requires zero authored damping"));
    }
    let chart = specimen.chart();
    if !additional.is_empty() {
        let terms = options.n_modes.checked_add(1)
            .and_then(|n| chart.mesh.nodes.len().checked_mul(n))
            .and_then(|n| chart.mesh.tris.len().checked_mul(3).and_then(|t| n.checked_add(t)))
            .and_then(|n| n.checked_mul(additional.len()))
            .ok_or_else(|| refuse("plate force-port projection work overflows"))?;
        if terms > max_projection_terms {
            return Err(refuse("plate force ports exceed their projection budget"));
        }
    }
    let reference = PlateProjection::with_force_weights(&chart.mesh, &options.unit_force_weights)?;
    let projections = additional.iter().map(|weights|
        PlateProjection::with_force_weights(&chart.mesh, weights)
    ).collect::<Result<Vec<_>, _>>()?;
    let model = chart.assemble(&[], &options.assembly).map_err(map_plate)?;
    let loss = specimen.bending_viscosity().map(|laws|
        RegionalBendingLoss::new(specimen, laws, &model)
    ).transpose()?;
    let report = modes(&model, window, &fs_modal::SliceOptions::default()).map_err(map_plate)?;
    let keep = options.n_modes.min(report.modes.len());
    let mut radiators = Vec::with_capacity(keep);
    let mut ports = vec![Vec::with_capacity(keep); projections.len()];
    let mut bending_zetas = Vec::with_capacity(keep);
    let mut normalized = Vec::new();
    for pair in report.modes.iter().take(keep) {
        let omega = pair.lambda.max(0.0).sqrt();
        if !(omega > 0.0) { continue; }
        // The existing owner computes the ACTUAL phi^T M phi, signed area and
        // reference drive. Every additional column uses precisely that mass.
        let radiator = reference.mode(&model, &pair.phi, omega, options.damping_ratio)?;
        let root_mass = fs_math::det::sqrt(radiator.mass_kg);
        if let Some((_, _, gate)) = tensor {
            if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
            normalized.push(pair.phi.iter().map(|p| p / root_mass).collect::<Vec<_>>());
        }
        for (projection, port) in projections.iter().zip(&mut ports) {
            let mut drive = 0.0;
            for (node, &weight) in projection.unit_force.iter().enumerate() {
                if let Some(dof) = model.dof_map[3 * node] { drive += weight * pair.phi[dof]; }
            }
            let coefficient = drive / root_mass;
            if !coefficient.is_finite() { return Err(refuse("plate force-port normalization is not finite")); }
            port.push(coefficient);
        }
        bending_zetas.push(match &loss {
            Some(loss) => loss.modal_zeta(&model, &pair.phi, omega)?,
            None => 0.0,
        });
        radiators.push(radiator);
    }
    if radiators.is_empty() { return Err(refuse("plate modal window returned no usable modes")); }
    // Preserve the original authored Rayleigh convention, with material loss
    // added separately. This does not give authored damping source authority.
    if radiators.len() >= 2 {
        let w0 = radiators[0].omega;
        if let Ok(rayleigh) = RayleighDamping::from_two_points(
            w0, options.damping_ratio, 4.0 * w0, options.damping_ratio,
        ) {
            for body in &mut radiators { body.zeta = rayleigh.zeta_at(body.omega); }
        }
    }
    for (body, bending_zeta) in radiators.iter_mut().zip(bending_zetas) { body.zeta += bending_zeta; }
    let damping = tensor.map(|(law, limits, gate)| {
        let omega: Vec<_> = radiators.iter().map(|r|r.omega).collect();
        law.project(specimen, &model, &normalized, &omega, limits, gate)
    }).transpose()?;
    Ok((radiators, ports, damping))
}

#[cfg(test)]
mod tests;

mod pressure;

#[cfg(test)]
mod tensor_tests;
