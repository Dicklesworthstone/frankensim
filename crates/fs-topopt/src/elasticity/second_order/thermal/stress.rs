//! Source-bound physical stress constraints on a retained thermal equilibrium.
//!
//! The dimensionless objective is the UNNORMALIZED moment
//! J = sum_case,cell (von_mises(sigma_eff)/allowable)^p. Every case participates,
//! including a case whose compliance weight is zero. Thus J <= 1 implies each
//! sampled effective von Mises stress is below the specified allowance in exact
//! arithmetic; neither small cells nor low load weights hide a hotspot.
//! Sampling/count dependence is explicit. This is not a continuum bound,
//! a material failure law, or solid-phase stress in a perforated microstructure.
//!
//! The moment, rather than its p-th root, has continuous second derivatives at
//! zero deviatoric stress for integer p >= 2. No artificial stress floor is used.
//! For s_e=m_e C_e(B_e u-e_th), one adjoint satisfies K lambda=J_u, and
//! g_e=J_s . C_e(B_e u-e_th) + lambda . (q_e-K_e u).
//! Exact Hessian actions differentiate both this adjoint and the native thermal
//! equilibrium. No finite differences or global Jacobian/Hessian are assembled.

use super::{ComplianceError, ThermalComplianceGoal, finite, gather, material_at};
use crate::SolveControl;

/// Explicit stress allowance and integer moment power. A larger power approaches
/// a sampled maximum after taking the p-th root, at greater numerical contrast.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalStressOptions {
    /// Positive physical von Mises allowance [Pa], shared by the sampled cells.
    /// The caller supplies this criterion; it is not inferred from stiffness.
    pub allowable_pa: f64,
    /// Integer power in 2..=64. The unnormalized moment is constrained to <= 1.
    pub power: u32,
}
impl ThermalStressOptions {
    pub(crate) fn admit(self) -> Result<(), ComplianceError> {
        if !(self.allowable_pa.is_finite() && self.allowable_pa > 0.0)
            || !(2..=64).contains(&self.power)
        {
            return Err(ComplianceError::InvalidInput(
                "thermal stress requires a positive finite Pa allowance and power in 2..=64",
            ));
        }
        Ok(())
    }
}

// Owned numerical data for the existing exact-point optimizer cache. Only this
// crate can detach it from the borrowing facade, and must keep the same owner.
pub(crate) struct StressLinearization {
    pub(crate) options: ThermalStressOptions,
    pub(crate) value: f64,
    pub(crate) maximum_ratio: f64,
    pub(crate) gradient: Vec<f64>,
    adjoints: Vec<Vec<f64>>,
    residuals: Vec<f64>,
}

/// A stress objective tied by its borrow to the original immutable thermal
/// design point. One adjoint per independent case is retained and reused.
pub struct ThermalStressGoal<'g, 'a> {
    source: &'g ThermalComplianceGoal<'a>,
    linearization: StressLinearization,
}
impl<'g, 'a> ThermalStressGoal<'g, 'a> {
    pub fn value(&self) -> f64 {
        self.linearization.value
    }
    pub fn gradient(&self) -> &[f64] {
        &self.linearization.gradient
    }
    /// Largest actual sampled von Mises stress divided by the specified allowance.
    pub fn maximum_ratio(&self) -> f64 {
        self.linearization.maximum_ratio
    }
    pub fn options(&self) -> ThermalStressOptions {
        self.linearization.options
    }
    pub fn adjoint_relative_residuals(&self) -> &[f64] {
        &self.linearization.residuals
    }
    pub fn source(&self) -> &'g ThermalComplianceGoal<'a> {
        self.source
    }

    /// Two native solves per case (primal tangent and adjoint tangent), with no
    /// new equilibrium or adjoint preparation. A zero direction costs no solve.
    pub fn hessian_vector_product(
        &self,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        self.source.stress_hessian(&self.linearization, direction, control)
    }
}

impl<'a> ThermalComplianceGoal<'a> {
    /// Differentiate the physical stress at this original thermal equilibrium.
    /// Compliance case weights are deliberately excluded from a strength limit.
    /// Prescribed temperature, source state, mesh and supports stay fixed.
    pub fn stress_goal<'g>(
        &'g self,
        options: ThermalStressOptions,
        control: &mut SolveControl<'_>,
    ) -> Result<ThermalStressGoal<'g, 'a>, ComplianceError> {
        let linearization = self.prepare_stress(options, control)?;
        Ok(ThermalStressGoal { source: self, linearization })
    }

    pub(crate) fn prepare_stress(
        &self,
        options: ThermalStressOptions,
        control: &mut SolveControl<'_>,
    ) -> Result<StressLinearization, ComplianceError> {
        control.checkpoint("thermal-stress-admission")?;
        options.admit()?;
        let el = self.base.elasticity;
        let mut point = StressLinearization {
            options,
            value: 0.0,
            maximum_ratio: 0.0,
            gradient: vec![0.0; el.cells()],
            adjoints: Vec::with_capacity(self.cases.len()),
            residuals: Vec::with_capacity(self.cases.len()),
        };
        for index in 0..self.cases.len() {
            let mut rhs = vec![0.0; el.n];
            for cell in 0..el.cells() {
                if cell % 256 == 0 {
                    control.checkpoint("thermal-stress-adjoint-load")?;
                }
                let s = &self.cases[index].response.elements[cell].effective_stress_mandel_pa;
                let (value, ratio, derivative) = moment(s, options)?;
                point.value = finite(point.value + value, "thermal stress moment")?;
                point.maximum_ratio = point.maximum_ratio.max(ratio);
                let reference = self.reference_stress(index, cell);
                point.gradient[cell] = finite(
                    point.gradient[cell] + dot(&derivative, &reference),
                    "thermal stress explicit gradient",
                )?;
                self.scatter_stress(cell, &derivative, self.base.moduli[cell], &mut rhs)?;
            }
            let (adjoint, residual) = self.base.solve(&rhs, control, "thermal-stress-adjoint")?;
            for (cell, tet) in el.tets.iter().enumerate() {
                if cell % 256 == 0 {
                    control.checkpoint("thermal-stress-gradient")?;
                }
                let local = gather(el, tet, &adjoint);
                let balance = self.local_design_load(index, cell);
                point.gradient[cell] = finite(
                    point.gradient[cell] + dot(&local, &balance),
                    "thermal stress total gradient",
                )?;
            }
            point.adjoints.push(adjoint);
            point.residuals.push(residual);
        }
        control.checkpoint("thermal-stress-publish")?;
        Ok(point)
    }

    pub(crate) fn stress_hessian(
        &self,
        point: &StressLinearization,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        self.admit_direction(direction, control)?;
        let el = self.base.elasticity;
        let mut output = vec![0.0; el.cells()];
        if direction.iter().all(|&v| v == 0.0) {
            return Ok(output);
        }
        for (index, adjoint) in point.adjoints.iter().enumerate() {
            let du = self.case_tangent(index, direction, control)?;
            let mut rhs = vec![0.0; el.n];
            for (cell, tet) in el.tets.iter().enumerate() {
                if cell % 256 == 0 {
                    control.checkpoint("thermal-stress-adjoint-tangent-load")?;
                }
                let local_du = gather(el, tet, &du);
                let elastic_tangent = self.stress_action(cell, &local_du);
                let reference = self.reference_stress(index, cell);
                let ds = core::array::from_fn(|i| {
                    direction[cell].mul_add(
                        reference[i],
                        self.base.moduli[cell] * elastic_tangent[i],
                    )
                });
                let s = &self.cases[index].response.elements[cell].effective_stress_mandel_pa;
                let (_, _, derivative) = moment(s, point.options)?;
                let derivative_tangent = moment_hessian(s, &ds, point.options)?;
                let force_derivative = core::array::from_fn(|i| {
                    direction[cell].mul_add(
                        derivative[i],
                        self.base.moduli[cell] * derivative_tangent[i],
                    )
                });
                self.scatter_stress(cell, &force_derivative, 1.0, &mut rhs)?;
                let local_adjoint = gather(el, tet, adjoint);
                let k_lambda = element_action(&el.ke[cell], &local_adjoint);
                for (local, &value) in k_lambda.iter().enumerate() {
                    let dof = 3 * tet[local / 3] as usize + local % 3;
                    if el.free[dof] {
                        rhs[dof] = finite(
                            rhs[dof] - direction[cell] * value,
                            "thermal stress adjoint tangent RHS",
                        )?;
                    }
                }
                let k_du = element_action(&el.ke[cell], &local_du);
                output[cell] = finite(
                    output[cell]
                        + dot(&derivative_tangent, &reference)
                        + dot(&derivative, &elastic_tangent)
                        - dot(&local_adjoint, &k_du),
                    "thermal stress Hessian explicit action",
                )?;
            }
            let (dlambda, _) = self.base.solve(&rhs, control, "thermal-stress-adjoint-tangent")?;
            for (cell, tet) in el.tets.iter().enumerate() {
                if cell % 256 == 0 {
                    control.checkpoint("thermal-stress-hessian")?;
                }
                output[cell] = finite(
                    output[cell] + dot(
                        &gather(el, tet, &dlambda),
                        &self.local_design_load(index, cell),
                    ),
                    "thermal stress Hessian action",
                )?;
            }
        }
        control.checkpoint("thermal-stress-hessian-publish")?;
        Ok(output)
    }

    fn reference_stress(&self, index: usize, cell: usize) -> [f64; 6] {
        let c = material_at(self.material_field(), cell).stiffness_mandel_pa();
        let strain = &self.cases[index].response.elements[cell].elastic_strain_mandel;
        core::array::from_fn(|i| dot(&c[i], strain))
    }

    fn local_design_load(&self, index: usize, cell: usize) -> [f64; 12] {
        let el = self.base.elasticity;
        let local = gather(el, &el.tets[cell], &self.base.displacements[index]);
        let ku = element_action(&el.ke[cell], &local);
        core::array::from_fn(|i| self.cases[index].element_forces_n[cell][i] - ku[i])
    }

    // Existing CB rows return physical shear stress; scale the OUTPUT rows to
    // Mandel. Its transpose then pairs with a Mandel stress derivative exactly.
    fn stress_action(&self, cell: usize, local: &[f64; 12]) -> [f64; 6] {
        core::array::from_fn(|row| {
            self.base.elasticity.cb[cell][12 * row..12 * (row + 1)]
                .iter().zip(local).map(|(c, u)| c * u).sum::<f64>()
                * mandel_factor(row)
        })
    }

    fn scatter_stress(
        &self,
        cell: usize,
        derivative: &[f64; 6],
        weight: f64,
        rhs: &mut [f64],
    ) -> Result<(), ComplianceError> {
        let el = self.base.elasticity;
        for local in 0..12 {
            let dof = 3 * el.tets[cell][local / 3] as usize + local % 3;
            if el.free[dof] {
                let mut force = 0.0;
                for (row, &d) in derivative.iter().enumerate() {
                    force += d * mandel_factor(row) * el.cb[cell][12 * row + local];
                }
                rhs[dof] = finite(
                    weight.mul_add(force, rhs[dof]),
                    "thermal stress adjoint load",
                )?;
            }
        }
        Ok(())
    }
}

fn mandel_factor(row: usize) -> f64 {
    if row < 3 { 1.0 } else { core::f64::consts::SQRT_2 }
}
fn dot<const N: usize>(a: &[f64; N], b: &[f64; N]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn element_action(k: &[f64; 144], x: &[f64; 12]) -> [f64; 12] {
    core::array::from_fn(|i| {
        k[12 * i..12 * (i + 1)].iter().zip(x).map(|(k, x)| k * x).sum()
    })
}

// Q = (3/2) times the Mandel deviatoric projector, vm(s)^2 = s^T Q s.
fn q_action(s: &[f64; 6]) -> [f64; 6] {
    [
        s[0] - 0.5 * s[1] - 0.5 * s[2],
        s[1] - 0.5 * s[0] - 0.5 * s[2],
        s[2] - 0.5 * s[0] - 0.5 * s[1],
        1.5 * s[3], 1.5 * s[4], 1.5 * s[5],
    ]
}

// Return dimensionless von Mises ratio and its derivative with respect to s/a.
// Scaling only evaluates the algebraic norm; it does not freeze its derivative.
fn ratio(s: &[f64; 6], allowable: f64) -> Result<(f64, [f64; 6]), ComplianceError> {
    let mut scaled = [0.0; 6];
    for i in 0..6 {
        scaled[i] = finite(s[i] / allowable, "normalized physical stress")?;
    }
    let scale = scaled.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return Ok((0.0, [0.0; 6]));
    }
    for value in &mut scaled {
        *value /= scale;
    }
    let square = 0.5 * (
        (scaled[0] - scaled[1]).powi(2)
            + (scaled[1] - scaled[2]).powi(2)
            + (scaled[2] - scaled[0]).powi(2)
    ) + 1.5 * (scaled[3].powi(2) + scaled[4].powi(2) + scaled[5].powi(2));
    let root = fs_math::det::sqrt(square);
    if root == 0.0 {
        return Ok((0.0, [0.0; 6]));
    }
    let q = q_action(&scaled);
    Ok((
        finite(scale * root, "physical stress ratio")?,
        core::array::from_fn(|i| q[i] / root),
    ))
}

fn moment(
    s: &[f64; 6],
    options: ThermalStressOptions,
) -> Result<(f64, f64, [f64; 6]), ComplianceError> {
    let (r, unit) = ratio(s, options.allowable_pa)?;
    let power = f64::from(options.power);
    let value = finite(fs_math::det::pow(r, power), "physical stress moment")?;
    let slope = finite(
        power * fs_math::det::pow(r, power - 1.0) / options.allowable_pa,
        "physical stress moment slope",
    )?;
    let mut gradient = [0.0; 6];
    for i in 0..6 {
        gradient[i] = finite(slope * unit[i], "physical stress derivative")?;
    }
    Ok((value, r, gradient))
}

fn moment_hessian(
    s: &[f64; 6],
    ds: &[f64; 6],
    options: ThermalStressOptions,
) -> Result<[f64; 6], ComplianceError> {
    let (r, unit) = ratio(s, options.allowable_pa)?;
    let power = f64::from(options.power);
    let mut dr = [0.0; 6];
    for i in 0..6 {
        dr[i] = finite(ds[i] / options.allowable_pa, "physical stress tangent")?;
    }
    let scale = if options.power == 2 {
        2.0
    } else {
        power * fs_math::det::pow(r, power - 2.0)
    };
    let q = q_action(&dr);
    let along = dot(&unit, &dr);
    let mut output = [0.0; 6];
    for i in 0..6 {
        output[i] = finite(
            scale * (q[i] + (power - 2.0) * along * unit[i]) / options.allowable_pa,
            "physical stress curvature",
        )?;
    }
    Ok(output)
}
