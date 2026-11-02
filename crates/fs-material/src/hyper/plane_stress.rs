//! Finite-strain plane stress from the original three-dimensional energy.
//!
//! For in-plane `a = [F00,F01,F10,F11]`, relax the entire transverse column
//! `t = [F02,F12,F22]` in `F = [a0,a1,t0; a2,a3,t1; 0,0,t2]` until
//! `P[:,2] = dW/dt = 0`. Since `F[2,:] = [0,0,t2]`, these are exactly the
//! three traction-free conditions `sigma[:,2] = 0`. Relaxing only F22 would
//! leave shear tractions for generally oriented anisotropic materials.
//!
//! The condensed tangent is `Aaa - Aat Att^-1 Ata`, from differentiating
//! the transverse equilibrium, not from discarding rows of the 3-D tangent.
//! Quadrature-local Newton and consistent condensation follow the procedure
//! of Ahmadi et al. (2025), doi:10.1007/s00707-025-04310-z. Here the same
//! variational elimination includes both transverse shear components.
//! This is a local stable stationary branch, not a global minimization or
//! a shell, bending, wrinkling, incompressibility or uniqueness certificate.

use super::Hyperelastic;
use crate::MaterialError;
use fs_ad::{Dual, Dual64};

const IN_PLANE: [usize; 4] = [0, 1, 3, 4];
const TRANSVERSE: [usize; 3] = [2, 5, 8];

/// Numerical budgets for the local transverse-equilibrium solve.
#[derive(Clone, Copy, Debug)]
pub struct PlaneStressSettings {
    /// Traction norm divided by the initial characteristic transverse stiffness
    /// `1/trace(Att^-1)` times `max(1, ||F_in_plane||)`. The stiffness stays
    /// bounded by the softest transverse mode as the bulk modulus grows.
    pub relative_tolerance: f64,
    /// Maximum accepted Newton updates, excluding the initial evaluation.
    pub max_iterations: usize,
    /// Maximum trial step lengths per Newton update (1, 1/2, 1/4, ...).
    pub max_backtracks: usize,
}

impl Default for PlaneStressSettings {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-11,
            max_iterations: 32,
            max_backtracks: 24,
        }
    }
}

/// Converged energy and stress without an in-plane material Hessian.
#[derive(Clone, Copy, Debug)]
pub struct PlaneStressState {
    /// Full row-major three-dimensional gradient, including relaxed thickness.
    pub deformation_gradient: [f64; 9],
    /// Stored energy per original reference volume (Pa = J/m³).
    pub energy_density_pa: f64,
    /// First Piola stress in order `[P00,P01,P10,P11]`, in Pa.
    pub piola: [f64; 4],
    /// Final transverse traction norm divided by the fixed initial scale.
    pub relative_transverse_residual: f64,
    /// Number of accepted Newton updates.
    pub iterations: usize,
}

/// Converged plane-stress response and its consistent in-plane tangent.
#[derive(Clone, Copy, Debug)]
pub struct PlaneStressResponse {
    /// Full row-major three-dimensional gradient, including relaxed thickness.
    pub deformation_gradient: [f64; 9],
    /// Stored energy per original reference volume (Pa = J/m³).
    pub energy_density_pa: f64,
    /// First Piola stress in order `[P00,P01,P10,P11]`, in Pa.
    pub piola: [f64; 4],
    /// `dP_in_plane/dF_in_plane`, with the transverse equilibrium eliminated.
    pub tangent: [[f64; 4]; 4],
    /// Final transverse traction norm divided by the fixed initial scale.
    pub relative_transverse_residual: f64,
    /// Number of accepted Newton updates.
    pub iterations: usize,
}

fn state_error(what: &str) -> MaterialError {
    MaterialError::State {
        what: format!("hyperelastic plane stress: {what}"),
    }
}

// Factor the 3x3 eliminated block after scalar normalization. No inverse,
// regularization, stiffness floor, or pivot clipping is introduced.
struct PositiveBlock {
    lower: [[f64; 3]; 3],
    scale: f64,
}
impl PositiveBlock {
    fn new(a: [[f64; 3]; 3]) -> Result<Self, MaterialError> {
        let scale = a.iter().flatten().map(|v| v.abs()).fold(0.0, f64::max);
        if !(scale.is_finite() && scale > 0.0) || a.iter().flatten().any(|v| !v.is_finite()) {
            return Err(state_error("nonfinite or zero transverse stiffness"));
        }
        let mut lower = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..=i {
                if (a[i][j] / scale - a[j][i] / scale).abs() > 128.0 * f64::EPSILON {
                    return Err(state_error("nonsymmetric transverse stiffness"));
                }
                let value =
                    a[i][j] / scale - (0..j).map(|k| lower[i][k] * lower[j][k]).sum::<f64>();
                if i == j {
                    if !(value.is_finite() && value > 0.0) {
                        return Err(state_error("transverse stationary branch is not stable"));
                    }
                    lower[i][j] = fs_math::det::sqrt(value);
                } else {
                    lower[i][j] = value / lower[j][j];
                }
            }
        }
        Ok(Self { lower, scale })
    }

    fn solve(&self, rhs: [f64; 3]) -> Result<[f64; 3], MaterialError> {
        let mut x = [0.0; 3];
        for i in 0..3 {
            x[i] = (rhs[i] / self.scale - (0..i).map(|j| self.lower[i][j] * x[j]).sum::<f64>())
                / self.lower[i][i];
        }
        for i in (0..3).rev() {
            x[i] = (x[i] - (i + 1..3).map(|j| self.lower[j][i] * x[j]).sum::<f64>())
                / self.lower[i][i];
        }
        if x.iter().any(|v| !v.is_finite()) {
            return Err(state_error("unrepresentable transverse solve"));
        }
        Ok(x)
    }

    fn traction_scale(&self) -> Result<f64, MaterialError> {
        let mut compliance = 0.0;
        for i in 0..3 {
            let mut rhs = [0.0; 3];
            rhs[i] = self.scale;
            compliance += self.solve(rhs)?[i];
        }
        // Inverse trace of compliance is <= the smallest stiffness eigenvalue.
        // Unlike max |Att|, it cannot hide traction behind a large bulk penalty.
        let scale = self.scale / compliance;
        if !(scale.is_finite() && scale > 0.0) {
            return Err(state_error("unrepresentable transverse compliance scale"));
        }
        Ok(scale)
    }
}

fn relative_residual(p: &[f64; 9], scale: f64) -> f64 {
    let r = TRANSVERSE.map(|i| p[i] / scale);
    let largest = r.iter().map(|v| v.abs()).fold(0.0, f64::max);
    if largest == 0.0 || !largest.is_finite() {
        return largest;
    }
    largest * fs_math::det::sqrt(r.iter().map(|v| (v / largest).powi(2)).sum())
}

impl Hyperelastic {
    // Three-direction AD on the ORIGINAL stored energy. A local Newton step
    // does not need nine sweeps with nine gradient lanes. Domain admission is
    // performed by the accompanying native `piola` evaluation.
    fn transverse_tangent(&self, f: &[f64; 9]) -> [[f64; 3]; 3] {
        let mut h = [[0.0; 3]; 3];
        for col in 0..3 {
            let mut vars = f.map(|v| Dual::<Dual64<3>, 1>::constant(Dual64::constant(v)));
            for (row, &index) in TRANSVERSE.iter().enumerate() {
                vars[index] = Dual {
                    re: Dual64::variable(f[index], row),
                    eps: [Dual64::constant(if row == col { 1.0 } else { 0.0 })],
                };
            }
            let w = self.energy(&vars);
            for row in 0..3 {
                h[row][col] = w.eps[0].eps[row];
            }
        }
        h
    }

    /// Relax all three transverse tractions using bounded, residual-decreasing
    /// Newton steps, without forming the full 9×9 constitutive Hessian.
    ///
    /// The initial transverse column is `[0,0,1]`. Every trial must retain
    /// positive orientation and the original material's principal-stretch
    /// domain. An unstable transverse block or exhausted solve budget refuses.
    /// Iteration and backtracking are bounded; spatial callers poll cancellation
    /// between quadrature points. No history or warm-start state is mutated.
    pub fn plane_stress_state(
        &self,
        in_plane: &[f64; 4],
        settings: PlaneStressSettings,
    ) -> Result<PlaneStressState, MaterialError> {
        if !(settings.relative_tolerance.is_finite()
            && settings.relative_tolerance > 0.0
            && settings.relative_tolerance < 1.0)
            || settings.max_iterations == 0
            || settings.max_backtracks == 0
        {
            return Err(MaterialError::Parameters {
                what: "plane stress needs 0 < relative_tolerance < 1 and positive iteration/backtracking budgets".into(),
            });
        }
        let determinant = in_plane[0] * in_plane[3] - in_plane[1] * in_plane[2];
        if in_plane.iter().any(|v| !v.is_finite())
            || !(determinant.is_finite() && determinant > 0.0)
        {
            return Err(state_error(
                "finite orientation-preserving in-plane gradient required",
            ));
        }
        let mut f = [
            in_plane[0],
            in_plane[1],
            0.0,
            in_plane[2],
            in_plane[3],
            0.0,
            0.0,
            0.0,
            1.0,
        ];
        let mut p = self.piola(&f)?;
        let mut block = PositiveBlock::new(self.transverse_tangent(&f))?;
        let gradient_norm = in_plane[0]
            .hypot(in_plane[1])
            .hypot(in_plane[2])
            .hypot(in_plane[3]);
        let scale = block.traction_scale()? * gradient_norm.max(1.0);
        if !(scale.is_finite() && scale > 0.0) {
            return Err(state_error("unrepresentable traction normalization"));
        }
        for iteration in 0..=settings.max_iterations {
            let residual = relative_residual(&p, scale);
            if residual <= settings.relative_tolerance {
                let energy = self.energy(&f);
                if !energy.is_finite() {
                    return Err(state_error("nonfinite stored energy"));
                }
                return Ok(PlaneStressState {
                    deformation_gradient: f,
                    energy_density_pa: energy,
                    piola: IN_PLANE.map(|i| p[i]),
                    relative_transverse_residual: residual,
                    iterations: iteration,
                });
            }
            if iteration == settings.max_iterations || !residual.is_finite() {
                break;
            }
            let direction = block.solve(TRANSVERSE.map(|i| -p[i]))?;
            let mut alpha = 1.0;
            let mut accepted = None;
            for _ in 0..settings.max_backtracks {
                let mut trial = f;
                for i in 0..3 {
                    trial[TRANSVERSE[i]] += alpha * direction[i];
                }
                if trial[8] > 0.0 {
                    match self.piola(&trial) {
                        Ok(stress) => {
                            if relative_residual(&stress, scale) <= (1.0 - 1e-4 * alpha) * residual
                            {
                                // Every accepted iterate remains on a stable
                                // transverse branch, not merely a traction root.
                                if let Ok(next) =
                                    PositiveBlock::new(self.transverse_tangent(&trial))
                                {
                                    accepted = Some((trial, stress, next));
                                    break;
                                }
                            }
                        }
                        Err(MaterialError::State { .. }) => (),
                        Err(error) => return Err(error),
                    }
                }
                alpha *= 0.5;
            }
            let Some((next_f, next_p, next_block)) = accepted else {
                return Err(state_error(
                    "transverse line search exhausted its budget or material domain",
                ));
            };
            f = next_f;
            p = next_p;
            block = next_block;
        }
        Err(state_error(
            "transverse Newton solve exhausted its iteration budget",
        ))
    }

    /// Traction-free energy/stress and the exact Schur-complement tangent of
    /// the locally converged transverse equilibrium.
    pub fn plane_stress(
        &self,
        in_plane: &[f64; 4],
        settings: PlaneStressSettings,
    ) -> Result<PlaneStressResponse, MaterialError> {
        let state = self.plane_stress_state(in_plane, settings)?;
        let a = self.tangent(&state.deformation_gradient)?;
        let block = PositiveBlock::new(TRANSVERSE.map(|i| TRANSVERSE.map(|j| a[i][j])))?;
        let mut tangent = [[0.0; 4]; 4];
        for (j, &col) in IN_PLANE.iter().enumerate() {
            let extension = block.solve(TRANSVERSE.map(|k| a[k][col]))?;
            for (i, &row) in IN_PLANE.iter().enumerate() {
                tangent[i][j] = a[row][col]
                    - (0..3)
                        .map(|k| a[row][TRANSVERSE[k]] * extension[k])
                        .sum::<f64>();
            }
        }
        if tangent.iter().flatten().any(|v| !v.is_finite()) {
            return Err(state_error("unrepresentable condensed tangent"));
        }
        Ok(PlaneStressResponse {
            deformation_gradient: state.deformation_gradient,
            energy_density_pa: state.energy_density_pa,
            piola: state.piola,
            tangent,
            relative_transverse_residual: state.relative_transverse_residual,
            iterations: state.iterations,
        })
    }
}
