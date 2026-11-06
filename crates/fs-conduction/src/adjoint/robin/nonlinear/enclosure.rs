//! Bind the original thermal tangent only after the radiation residual owner
//! has checked the unchanged field with explicit physical nodal loads.

use super::*;
use crate::interface::temperature::tangent::ContactTemperatureLinearization;
use crate::radiation::{RadiationSurface, RadiosityReport, equilibrium};

impl super::super::RobinResponse {
    /// Internal steady radiation producer. No shifted Robin reference enters
    /// admission or the base tangent; the caller composes the complete physical
    /// flux feedback separately after this unchanged-field watt audit.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn for_steady_enclosure(
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        contact: Option<&ContactTemperatureLinearization>,
        linear: LinearConfig,
        temperature: &[f64],
        regions: &[&str],
        surfaces: &[RadiationSurface],
        radiosity: RadiosityReport,
        residual_tolerance_w: f64,
        energy_tolerance_w: f64,
    ) -> Result<(Self, equilibrium::Endpoint), ConductionError> {
        poll(cx, 0)?;
        super::super::admit_linear(linear)?;
        if linear.restart == 0
            || !(residual_tolerance_w.is_finite() && residual_tolerance_w > 0.0)
            || !(energy_tolerance_w.is_finite() && energy_tolerance_w > 0.0)
        {
            return Err(invalid(
                "steady enclosure derivatives require a positive FGMRES restart and finite positive watt tolerances",
            ));
        }
        super::super::vector(cx, temperature, problem.mesh.vertex_count())?;
        let dofs = DofMap::new(problem.boundary, problem.mesh.vertex_count())?;
        if dofs.fixed().is_empty() && !problem.boundary.has_robin() {
            return Err(ConductionError::SingularPureNeumann);
        }
        for &vertex in dofs.fixed() {
            poll(cx, vertex)?;
            if temperature[vertex].to_bits() != dofs.prescribed()[vertex].to_bits() {
                return Err(invalid(
                    "steady enclosure endpoint differs from its original prescribed temperatures",
                ));
            }
        }
        let interfaces = if let Some(contact) = contact {
            if temperature
                .iter()
                .map(|v| v.to_bits())
                .ne(contact.temperature().iter().map(|v| v.to_bits()))
            {
                return Err(invalid(
                    "contact tangent differs from the enclosure endpoint field",
                ));
            }
            Some(contact.physical_contacts().frozen_interfaces())
        } else {
            interfaces
        };
        let endpoint = equilibrium::endpoint(
            cx,
            problem,
            interfaces,
            temperature,
            residual_tolerance_w,
            radiosity.clone(),
            &radiosity,
            &|_| Ok(problem.boundary.clone()),
            Some(surfaces),
        )?;
        for value in [
            endpoint.energy.source_w,
            endpoint.energy.neumann_out_w,
            endpoint.energy.robin_out_w,
            endpoint.energy.dirichlet_in_w,
            endpoint.energy.closure_w,
            endpoint.energy.scale_w,
        ] {
            checked(value)?;
        }
        if endpoint.report.residual_w > residual_tolerance_w
            || endpoint.energy.closure_w.abs() > energy_tolerance_w
        {
            return Err(invalid(
                "steady enclosure endpoint fails its physical free-node residual or heat-balance target",
            ));
        }
        // Only the physical endpoint owner above can admit this tangent.
        // Original convection supplies consistent face mass; K(T) and the
        // frozen contact block come from the production material assembly.
        let ports = super::super::bind_ports(cx, problem, regions)?;
        let (matrix, nonlinear) = prepare(cx, problem, interfaces, temperature, &dofs)?;
        let mut response = Self {
            temperature: temperature.to_vec(),
            robin_fluxes: endpoint.robin_fluxes.clone(),
            matrix,
            dofs,
            ports,
            linear,
            nonlinear: Some(nonlinear),
        };
        if let Some(contact) = contact {
            response = response.with_contact_temperature_feedback(cx, contact)?;
        }
        poll(cx, 0)?;
        Ok((response, endpoint))
    }
}
