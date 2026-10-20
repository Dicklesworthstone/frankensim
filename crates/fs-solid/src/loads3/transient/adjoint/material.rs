//! Material-design derivatives of accepted P1 dynamic equilibrium equations.
//!
//! Each element has independent multipliers sK and sRho, evaluated at one:
//! C_e(sK)=sK*C_e and rho_e(sRho)=sRho*rho_e. This retains anisotropy, Poisson
//! ratio and material axes; it is not a derivative of each tensor entry. Geometry,
//! restraints, authored support kinematics, free thermal strain and the external
//! NODAL force remain fixed. Full prescribed-support inertia and damping are
//! included. For mass-dependent external forcing, add its explicit force
//! derivative through the returned force cotangent.
use super::*;

/// Bounds for full-node scratch and element-gradient output, checked before
/// allocation. Scratch retains four full xyz fields (three motion, one seed).
#[derive(Debug, Clone, Copy)]
pub struct MaterialGradientBudget3 {
    pub maximum_full_dofs: usize,
    pub maximum_elements: usize,
}
impl Default for MaterialGradientBudget3 {
    fn default() -> Self { Self { maximum_full_dofs: 150_000, maximum_elements: 100_000 } }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MaterialScaleGradient3 {
    /// [dJ/dsK, dJ/dsRho] per original tetrahedron, at scales one.
    /// A common material multiplier uses the SUM of its element entries.
    pub element_scales: Vec<[f64; 2]>,
    /// [dJ/deta, dJ/dzeta] for C_damping=eta*M+zeta*K. The units are
    /// objective*s and objective/s respectively, unlike the scale gradients.
    pub rayleigh: [f64; 2],
}
impl MaterialScaleGradient3 {
    /// Chain independent per-element scale directions and global damping
    /// directions without constructing a parameter Jacobian. For example,
    /// d/dE with fixed nu uses dsK/dE=1/E at this material point.
    pub fn contract(&self, element_directions: &[[f64; 2]], rayleigh_direction: [f64; 2],
        cx: &Cx<'_>) -> Result<f64, TransientAdjointError3> {
        check(cx)?;
        if element_directions.len() != self.element_scales.len() {
            return Err(TransientAdjointError3::Invalid("one stiffness/density direction per element required"));
        }
        let mut result = 0.0;
        for (g, d) in self.element_scales.iter().zip(element_directions) {
            check(cx)?;
            for i in 0..2 { result = number(number(d[i])?.mul_add(g[i], result))?; }
        }
        for i in 0..2 { result = number(number(rayleigh_direction[i])?.mul_add(self.rayleigh[i], result))?; }
        check(cx)?; Ok(result)
    }
}

impl TetTransientRun3<'_, '_, '_, '_> {
    /// Contract a residual multiplier at this EXACT accepted endpoint. Supply
    /// `TransientPullback3::force` from the corresponding interval, initial mass
    /// solve or finite jump. This adds MATERIAL terms to that state/load
    /// pullback; it does not replace its previous-state or external-force terms.
    ///
    /// One bounded element sweep reuses canonical K_e, M_e and Mandel B. No
    /// new equilibrium solve, finite difference or global assembly is performed.
    /// Frozen free thermal strain is retained but its equivalent force changes
    /// with stiffness; both Rayleigh terms change with their respective matrices.
    /// Support-enabled runs use the actual full q/v/a, including the initial
    /// prescribed acceleration and all consistent-mass cross terms. Supply the
    /// force cotangent from the matching prescribed step or initialization.
    /// Direct parameter dependence of an objective must be added by its owner.
    pub fn material_scale_pullback(&self, force_cotangent: &[f64], budget: MaterialGradientBudget3,
        cx: &Cx<'_>) -> Result<MaterialScaleGradient3, TransientAdjointError3> {
        if self.prescribed.is_some() {
            self.admit_material_gradient(force_cotangent, budget, cx)?;
            let motion = self.nodal_motion(cx)?;
            return self.material_scale_pullback_nodal_at(&motion, force_cotangent, budget, cx);
        }
        self.material_scale_pullback_at(self.state(), force_cotangent, budget, cx)
    }

    pub(in crate::loads3::transient::adjoint) fn material_scale_pullback_at(
        &self, state: &DynamicState, force_cotangent: &[f64], budget: MaterialGradientBudget3,
        cx: &Cx<'_>,
    ) -> Result<MaterialScaleGradient3, TransientAdjointError3> {
        check(cx)?;
        if self.prescribed.is_some() {
            return Err(TransientAdjointError3::Invalid("prescribed material pullbacks require full nodal kinematics"));
        }
        let n = self.free_dofs().len();
        self.admit_material_gradient(force_cotangent, budget, cx)?;
        for values in [&state.displacement, &state.velocity, &state.acceleration] {
            validate_values(values, n, cx)?;
        }
        let nodes = self.geometry.nodes.len();
        let mut motion: [Vec<[f64; 3]>; 3] = std::array::from_fn(|_| vec![[0.0; 3]; nodes]);
        for (i, &d) in self.free_dofs().iter().enumerate() {
            if i % 256 == 0 { check(cx)?; }
            motion[0][d/3][d%3] = state.displacement[i];
            motion[1][d/3][d%3] = state.velocity[i];
            motion[2][d/3][d%3] = state.acceleration[i];
        }
        self.material_scale_pullback_nodal_at(&motion, force_cotangent, budget, cx)
    }
}
