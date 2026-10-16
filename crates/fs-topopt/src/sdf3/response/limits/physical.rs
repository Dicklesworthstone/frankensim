//! Shared independent equilibria for differentiable designs and binary endpoints.
//! Caller owns admission, material installation, preparation and rollback.
use super::*;
use fs_sparse::precond::Precond;

pub(super) struct LimitPrimals3 {
    pub displacements: Vec<Vec<f64>>,
    pub measurements: Vec<Vec<f64>>,
    pub residuals: Vec<Vec<f64>>,
    pub all: Vec<f64>,
}
impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    pub(super) fn response_limit_primals(&self, cases: &[ResponseLimitCase3<'_>],
        prepared: &impl Precond, control: &mut SolveControl<'_>) -> Result<LimitPrimals3, ResponseError3> {
        let op = self.operator.adaptive();
        let mut displacements = Vec::with_capacity(cases.len());
        let mut measurements = Vec::with_capacity(cases.len());
        let mut residuals = Vec::with_capacity(cases.len());
        let mut all = Vec::new();
        for case in cases {
            control.checkpoint("response-limits-case")?;
            let mut rhs: Vec<f64> = case.force.iter().enumerate()
                .map(|(i, &f)| if op.fixed()[i/3] { 0.0 } else { f }).collect();
            if let Some(g) = case.prescribed {
                let lifting = op.prescribed_displacement_load(g, || poll(control))?;
                for (r, b) in rhs.iter_mut().zip(lifting) { *r = finite(*r+b)?; }
            }
            let u = checked_solve_preconditioned(&self.operator, prepared, &rhs, 1e-11,
                "response-limits-primal", control)?;
            let mut values = Vec::with_capacity(case.limits.len());
            let mut row = Vec::with_capacity(case.limits.len());
            for limit in case.limits {
                control.checkpoint("response-limits-measure")?;
                let value = finite(match limit.response {
                    BoundedResponse3::Displacement(q) => dot(q, &u),
                    BoundedResponse3::Reaction(mode) =>
                        op.embedded_reaction(&u, case.prescribed, mode, || poll(control))?.value,
                })?;
                let r = finite(limit.sense.sign()*finite(value-limit.bound)?/limit.scale)?;
                values.push(value); row.push(r); all.push(r);
            }
            displacements.push(u); measurements.push(values); residuals.push(row);
        }
        Ok(LimitPrimals3 { displacements, measurements, residuals, all })
    }
}
