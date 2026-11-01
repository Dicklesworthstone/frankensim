//! Attach the exact in-crate source-law correction to a checked endpoint.
use super::*;
use crate::interface::temperature::tangent::ContactTemperatureLinearization;

impl super::super::RobinResponse {


    /// Bind a complete steady contact response only after rechecking the
    /// supplied physical field, prescribed values and watt balance.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub(crate) fn for_steady_contact(
        cx: &Cx<'_>, problem: ConductionProblem<'_>,
        contact: &ContactTemperatureLinearization, linear: LinearConfig,
        residual_tolerance_w: f64, energy_tolerance_w: f64, regions: &[&str],
    ) -> Result<(Self, f64, crate::EnergyBalance, Vec<(usize, f64)>), ConductionError> {
        poll(cx, 0)?;
        super::super::admit_linear(linear)?;
        if linear.restart == 0
            || !(residual_tolerance_w.is_finite() && residual_tolerance_w > 0.0)
            || !(energy_tolerance_w.is_finite() && energy_tolerance_w > 0.0)
        {
            return Err(invalid("steady contact derivatives require a positive FGMRES restart and finite positive watt tolerances"));
        }
        let temperature = contact.temperature();
        super::super::vector(cx, temperature, problem.mesh.vertex_count())?;
        let dofs = DofMap::new(problem.boundary, problem.mesh.vertex_count())?;
        if dofs.fixed().is_empty() && !problem.boundary.has_robin() {
            return Err(ConductionError::SingularPureNeumann);
        }
        for &vertex in dofs.fixed() {
            poll(cx, vertex)?;
            if temperature[vertex].to_bits() != dofs.prescribed()[vertex].to_bits() {
                return Err(invalid("steady contact endpoint differs from its prescribed temperature values"));
            }
        }
        let interfaces = contact.physical_contacts().frozen_interfaces();
        let system = crate::assemble::assemble_operator_scaled_with_interfaces(
            cx, problem.mesh, problem.boundary, problem.material, problem.source,
            temperature, None, Some(interfaces), problem.element_materials,
        )?;
        let residual = crate::assemble::full_residual(&system, temperature);
        for (i, &value) in residual.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            checked(value)?;
        }
        let residual_norm_w = checked(fs_solver::norm2(&dofs.gather(&residual)))?;
        let (energy, robin_fluxes) = crate::solve::energy_balance(
            problem.mesh, problem.boundary, problem.source, &system, &dofs, temperature,
        );
        for value in [energy.source_w, energy.neumann_out_w, energy.robin_out_w,
            energy.dirichlet_in_w, energy.closure_w, energy.scale_w] { checked(value)?; }
        if residual_norm_w > residual_tolerance_w || energy.closure_w.abs() > energy_tolerance_w {
            return Err(invalid("steady contact endpoint fails its physical free-node residual or heat-balance target"));
        }
        let reactions = dofs.fixed().iter().map(|&vertex| (vertex, residual[vertex])).collect();
        let ports = super::super::bind_ports(cx, problem, regions)?;
        // Original material assembly supplies K + K'(T)T and frozen contact.
        // The source-law owner adds only dK_contact/dT * T, exactly once.
        let (matrix, nonlinear) = prepare(cx, problem, Some(interfaces), temperature, &dofs)?;
        let response = Self { temperature:temperature.to_vec(), robin_fluxes,
            matrix, dofs, ports, linear, nonlinear:Some(nonlinear) }
            .with_contact_temperature_feedback(cx, contact)?;
        poll(cx, 0)?;
        Ok((response, residual_norm_w, energy, reactions))
    }
    pub(crate) fn with_contact_temperature_feedback(
        mut self, cx: &Cx<'_>, contact: &ContactTemperatureLinearization,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        if self.temperature.iter().map(|v| v.to_bits())
            .ne(contact.temperature().iter().map(|v| v.to_bits()))
            || self.linear.restart == 0
        {
            return Err(invalid("contact feedback must match the checked endpoint and admit a positive FGMRES restart"));
        }
        let (extra, _) = reduce_matrix_and_lift(contact.feedback(), &self.dofs);
        let n = self.matrix.nrows();
        if extra.nrows() != n { return Err(invalid("contact feedback has a different free-node layout")); }
        let mut coo = fs_sparse::Coo::new(n, n);
        for row in 0..n {
            poll(cx, row)?;
            for matrix in [&self.matrix, &extra] {
                let (columns, values) = matrix.row(row);
                for (&column, &value) in columns.iter().zip(values) {
                    coo.push(row, column, checked(value)?);
                }
            }
        }
        let matrix = coo.assemble();
        for row in 0..n {
            poll(cx, row)?;
            for &value in matrix.row(row).1 { checked(value)?; }
        }
        // R'(T) makes contact nonsymmetric even with a constant bulk tensor.
        // Preserve the old material smoothness gate; do not fall back to PCG.
        let smooth = self.has_smooth_material_tangent();
        self.nonlinear = Some(TangentSystem { transpose: ops::transpose(&matrix), smooth });
        self.matrix = matrix;
        poll(cx, n)?;
        Ok(self)
    }
}
