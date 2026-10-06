//! Lift derivatives for prescribed temperatures on an unchanged thermal field.
use super::{ConductionError, ConductionProblem, Cx, ThermalInterfaces,
    assemble_jacobian_with_optional_interfaces, checked, invalid, material_is_smooth, poll};
use super::super::{RobinResponse, bind_ports, vector};

impl RobinResponse {
    /// Contract a complete steady nodal-temperature dual with prescribed values.
    ///
    /// Returns `w_D - J_FD^T lambda_F` at prescribed vertices, zero elsewhere.
    /// `nodal_weights` is the FULL temperature objective, including fixed nodes;
    /// the full load-adjoint vector must instead be zero on fixed nodes. Omitting
    /// the first term incorrectly gives zero for a goal at a prescribed vertex.
    /// The production full material/contact Jacobian supplies the lift, including
    /// the nonsymmetric K'(T) contribution. No additional linear solve occurs.
    ///
    /// Bind `problem.boundary` to the physical Robin values AT `temperature`.
    /// Selected regions additionally have local h/reference mean slopes and an
    /// m-by-m row-major ADDITIONAL reference-feedback matrix. Mean weights include
    /// fixed vertices: discarding those weights would drop boundary feedback from
    /// the lift. Radiation callers must include its coefficient and weighted
    /// reference slopes; air callers must include upstream reference feedback.
    /// Empty region/slope/matrix slices admit fixed Robin or Dirichlet-only laws.
    /// The caller owns the physical law, full primal/dual residual checks, and any
    /// objective dependence other than nodal temperature. This contraction is
    /// not a solver, proof of external dual binding, or derivative certificate.
    ///
    /// `max_entries` bounds each field, element/face traversal estimates, and
    /// the admitted assembled CSR entries. The existing contact/assembly producer
    /// owns its construction workspace; this is NOT a peak allocator/RSS bound.
    /// At most 64 feedback regions are admitted. No dense boundary matrix or
    /// explicit transpose is constructed. No field or prescribed value is edited.
    ///
    /// # Errors
    /// Invalid sizes/work limits, nonfinite inputs/arithmetic, inconsistent fixed
    /// values/dual, material kinks/domains, invalid interfaces, or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn prescribed_temperature_pullback(
        cx: &Cx<'_>, problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
        temperature: &[f64], nodal_weights: &[f64], nodal_load_adjoint: &[f64],
        regions: &[&str], h_slopes: &[f64], reference_slopes: &[f64],
        reference_feedback: &[f64], max_entries: usize,
    ) -> Result<Vec<f64>, ConductionError> {
        poll(cx, 0)?;
        let n = problem.mesh.vertex_count();
        let m = regions.len();
        let work = problem.mesh.element_count().checked_mul(32)
            .and_then(|v| problem.mesh.boundary().len().checked_mul(27).and_then(|b| v.checked_add(b)));
        if n > max_entries || m > 64 || m.checked_mul(m).is_none_or(|v| v > max_entries)
            || work.is_none_or(|v| v > max_entries) {
            return Err(invalid("prescribed-temperature pullback exceeds its operator/trace work allowance"));
        }
        vector(cx, temperature, n)?;
        vector(cx, nodal_weights, n)?;
        vector(cx, nodal_load_adjoint, n)?;
        vector(cx, h_slopes, m)?;
        vector(cx, reference_slopes, m)?;
        vector(cx, reference_feedback, m*m)?;
        if let Some(materials) = problem.element_materials { materials.validate_for(problem.mesh)?; }
        let mut fixed = vec![false; n];
        for &(v, value) in problem.boundary.dirichlet() {
            poll(cx, v)?;
            if v >= n || temperature[v] != value || nodal_load_adjoint[v] != 0.0 {
                return Err(invalid("prescribed-temperature pullback requires unchanged fixed values and a zero fixed-node dual"));
            }
            fixed[v] = true;
        }
        if !material_is_smooth(cx, problem, temperature)? {
            return Err(invalid("prescribed-temperature pullback cannot select a derivative at a material kink or endpoint"));
        }
        let ports = bind_ports(cx, problem, regions)?;
        let matrix = assemble_jacobian_with_optional_interfaces(cx, problem.mesh,
            problem.boundary, problem.material, temperature, interfaces, problem.element_materials)?;
        let mut entries = 0_usize;
        for v in 0..n {
            if v % 512 == 0 { poll(cx, v)?; }
            entries = entries.checked_add(matrix.row(v).0.len())
                .ok_or_else(|| invalid("prescribed-temperature operator count overflow"))?;
            if entries > max_entries { return Err(invalid("prescribed-temperature CSR exceeds its entry allowance")); }
        }
        let mut result = Vec::new();
        result.try_reserve_exact(n).map_err(|_| invalid("prescribed-temperature output allocation refused"))?;
        result.extend(nodal_weights.iter().copied());
        for (v, &lambda) in nodal_load_adjoint.iter().enumerate() {
            if v % 512 == 0 { poll(cx, v)?; }
            let (columns, values) = matrix.row(v);
            for (&w, &coefficient) in columns.iter().zip(values) {
                if fixed[w] { result[w] = checked(result[w] - checked(lambda*checked(coefficient)?)?)?; }
            }
        }
        let mut mean_bars = vec![0.0; m];
        let mut reference_bars = vec![0.0; m];
        for (i, port) in ports.iter().enumerate() {
            let mut coefficient_bar = 0.0;
            for (vertices, area) in &port.faces {
                poll(cx, i)?;
                for (a, &v) in vertices.iter().enumerate() {
                    reference_bars[i] = checked(reference_bars[i]
                        + checked(nodal_load_adjoint[v]*port.htc_w_m2_k*(area/3.0))?)?;
                    for (b, &w) in vertices.iter().enumerate() {
                        let mass = (area/12.0)*if a == b {2.0} else {1.0};
                        coefficient_bar = checked(coefficient_bar + checked(nodal_load_adjoint[v]
                            *mass*checked(temperature[w]-port.reference_k)?)?)?;
                    }
                }
            }
            mean_bars[i] = checked(checked(h_slopes[i]*coefficient_bar)?
                - checked(reference_slopes[i]*reference_bars[i])?)?;
        }
        for i in 0..m {
            poll(cx, i)?;
            for j in 0..m {
                mean_bars[j] = checked(mean_bars[j] - checked(reference_bars[i]*reference_feedback[i*m+j])?)?;
            }
        }
        for (port, bar) in ports.iter().zip(mean_bars) {
            for (vertices, area) in &port.faces {
                poll(cx, 0)?;
                for &v in vertices {
                    if fixed[v] { result[v] = checked(result[v] - checked(bar*(area/port.area_m2/3.0))?)?; }
                }
            }
        }
        for (v, value) in result.iter_mut().enumerate() {
            if v % 512 == 0 { poll(cx, v)?; }
            if !fixed[v] { *value = 0.0; }
        }
        poll(cx, n)?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
