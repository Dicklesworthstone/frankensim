//! Analytic resistance controls for the CLI's uniform matching-P1 contacts.
//!
//! For p=ln(R''), the production contact block obeys dK_contact/dp=-K_contact.
//! The existing coupled nodal-load pullback is lambda, including air feedback
//! and the solid K'(T) Jacobian. Hence dQoI/dp=lambda^T K_contact T.
//! This is an exact derivative of the declared discrete operator, not a
//! finite-difference estimate, a frozen-air derivative or a mesh error bound.
//! The CLI's objectives depend on temperatures, not directly on contact flux.
use super::*;

/// Cached correspondence for faces ALREADY admitted by ThermalInterfaces::new.
/// Names, ownership and coincidence policy remain the production binder's.
#[derive(Debug)]
pub(super) struct Trace {
    a: [usize; 3],
    b: [usize; 3],
    area: f64,
}

impl Trace {
    pub(super) fn bind(mesh: &ConductionMesh, pairs: &[InterfaceFacePair]) -> Result<Vec<Self>> {
        // Match the production bound-face reduction order, not input ordering.
        let mut pairs = pairs.to_vec();
        pairs.sort_by_key(|p| (p.side_a.min(p.side_b), p.side_a.max(p.side_b)));
        pairs.iter().map(|pair| {
            let first = mesh.boundary().get(pair.side_a).ok_or_else(|| bad("missing admitted contact face"))?;
            let second = mesh.boundary().get(pair.side_b).ok_or_else(|| bad("missing admitted contact face"))?;
            let a = first.vertices.map(|v| v as usize);
            let mut b = [0; 3];
            for (i, &vertex) in a.iter().enumerate() {
                let key = mesh.positions()[vertex].map(f64::to_bits);
                b[i] = second.vertices.iter().find_map(|&other| {
                    (mesh.positions()[other as usize].map(f64::to_bits) == key).then_some(other as usize)
                }).ok_or_else(|| bad("admitted contact vertex correspondence disappeared"))?;
            }
            Ok(Self { a, b, area: first.area })
        }).collect()
    }

    fn contraction(&self, temperature: &[f64], adjoint: &[f64], resistance: f64) -> Result<f64> {
        let mut t = [0.0; 3];
        let mut lambda = [0.0; 3];
        for i in 0..3 {
            t[i] = checked(temperature[self.a[i]] - temperature[self.b[i]])?;
            lambda[i] = checked(adjoint[self.a[i]] - adjoint[self.b[i]])?;
        }
        // Integral of two P1 traces: M_ii=A/6, M_ij=A/12. Product
        // of mean jumps would omit the nonuniform within-face contribution.
        let conductance = checked(1.0 / resistance)?;
        let mut result = 0.0;
        for (i, &weight) in lambda.iter().enumerate() {
            for (j, &jump) in t.iter().enumerate() {
                let mass = if i == j { self.area / 6.0 } else { self.area / 12.0 };
                let contribution = checked(checked(conductance * mass)? * weight * jump)?;
                result = checked(result + contribution)?;
            }
        }
        Ok(result)
    }
}

impl Contacts {
    /// The supplied gradient must be the coupled load pullback at THIS field.
    /// Production calls this only from Evaluation's retained primal/adjoint.
    pub(crate) fn log_resistance_gradients(&self, temperature: &[f64], adjoint: &[f64]) -> Result<Vec<f64>> {
        if temperature.len() != self.vertex_count || adjoint.len() != self.vertex_count {
            return Err(bad("contact sensitivity requires the complete primal and coupled nodal-load adjoint"));
        }
        for &value in temperature.iter().chain(adjoint) { checked(value)?; }
        self.declarations.iter().map(|row| {
            row.traces.iter().try_fold(0.0, |sum, trace| {
                checked(sum + trace.contraction(temperature, adjoint, row.resistance)?)
            })
        }).collect()
    }

    pub(crate) fn sensitivity_json(&self, temperature: &[f64], gradient: Option<&CoupledGradient>) -> Result<String> {
        let Some(gradient) = gradient else { return Ok("null".into()); };
        let values = self.log_resistance_gradients(temperature, &gradient.nodal_load)?;
        let rows = self.declarations.iter().zip(values).map(|(row, value)| {
            Ok(format!("{{\"contact\":{},\"resistance_m2_k_w\":{},\"dobjective_dlog_resistance_k\":{},\"dobjective_dresistance_w_m2\":{}}}",
                quote(&row.name), num(row.resistance)?, num(value)?, num(checked(value / row.resistance)?)?))
        }).collect::<Result<Vec<_>>>()?.join(",");
        Ok(format!("{{\"method\":\"coupled-adjoint-contact-bilinear-form\",\"rows\":[{rows}],\"scope\":\"steady temperature objective; one scalar area-specific resistance per declared matching contact; full coupled air feedback and solid material tangent; geometry, materials, loads, hydraulic flow and convection laws fixed; no transient or contact-flux objective derivative\"}}"))
    }
}

fn checked(value: f64) -> Result<f64> {
    if value.is_finite() { Ok(value) } else { Err(producer("nonfinite contact sensitivity arithmetic")) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn consistent_mass_retains_nonuniform_contact_modes() {
        let trace = Trace { a: [0,1,2], b: [3,4,5], area: 6.0 };
        let t = [301.0,299.0,300.0,300.0,300.0,300.0];
        let lambda = [1.0,-1.0,0.0,0.0,0.0,0.0];
        // Both mean jumps vanish, yet lambda^T M jump / R = 1/2.
        assert_eq!(trace.contraction(&t, &lambda, 2.0).unwrap(), 0.5);
        let reversed = Trace { a: trace.b, b: trace.a, area: trace.area };
        assert_eq!(trace.contraction(&t, &lambda, 2.0).unwrap(),
            reversed.contraction(&t, &lambda, 2.0).unwrap());
    }

    #[test]
    fn constant_jump_matches_integrated_heat_and_rejects_nonfinite_values() {
        let trace = Trace { a: [0,1,2], b: [3,4,5], area: 3.0 };
        let t = [304.0,304.0,304.0,300.0,300.0,300.0];
        let lambda = [2.0,2.0,2.0,1.0,1.0,1.0];
        assert_eq!(trace.contraction(&t, &lambda, 2.0).unwrap(), 6.0);
        let mut poisoned = t; poisoned[0] = f64::NAN;
        assert!(trace.contraction(&poisoned, &lambda, 2.0).is_err());
    }
}
